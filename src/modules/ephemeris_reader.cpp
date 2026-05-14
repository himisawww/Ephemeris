#include"ephemeris_reader.h"
#include<algorithm>
#include"ephemeris_generator.h"
#include"configs.h"
#include"utils/logger.h"

ephemeris_reader::chapter::chapter(msystem &ms,const std::string &_chname)
    :izippack(_chname),chname(_chname),_interp_size(0),_bselect(-1){
    std::string failure;
    do{
        if(!izippack::operator bool())
            break;
        auto mf_files=load_central_directory();
        MFILE mf_index,mf_offset;
        bool index_loaded=false,chpt_loaded=false,offset_loaded=false;
        for(izipfile zf:mf_files){
            const std::string &zfn=zf.name();
            if(zfn==Configs::SaveNameIndex){
                if(index_loaded)failure+="    Duplicated index;\n";
                else{
                    index_loaded=true;
                    zf.fetch();
                    zf.dumpfile(mf_index);
                }
            }
            else if(zfn==Configs::SaveNameCheckpoint){
                if(chpt_loaded)failure+="    Duplicated checkpoint;\n";
                else{
                    chpt_loaded=true;
                    MFILE mf_chpt;
                    zf.fetch();
                    zf.dumpfile(mf_chpt);
                    msystem msc;
                    if(!msc.load_checkpoint(&mf_chpt))
                        failure+="    Invalid checkpoint;\n";
                    else if(!msc.is_same(ms))
                        failure+="    Incompatible celestial system between ephemerides;\n";
                }
            }
            else if(zfn==Configs::SaveNameBarycentricOffsetIndex){
                if(offset_loaded)failure+="    Duplicated barycentric offset index;\n";
                else{
                    offset_loaded=true;
                    zf.fetch();
                    zf.dumpfile(mf_offset);
                }
            }
            if(failure.size()||index_loaded&&chpt_loaded&&offset_loaded)break;
        }
        if(failure.size()||!(index_loaded&&chpt_loaded)){
            if(failure.empty())failure+="    Missing index/checkpoint;\n";
            break;
        }

        htl::map<std::string,izipfile> allfiles;
        for(const auto &zf:mf_files){
            if(!allfiles.try_emplace(get_file_name(zf.name()),zf).second){
                failure+="    Duplicated file;\n";
                break;
            }
        }
        if(failure.size())break;

        int dir=0;
        const size_t mn=ms.size();
        ephm_index.resize(mn);
        htl::vector<ephemeris_entry> bindices;
        while(failure.empty()){
            ephemeris_entry index;
            if(1!=fread(&index,sizeof(index),1,&mf_index))
                break;
            bool fwd=index.t_start<index.t_end;
            bool bak=index.t_start>index.t_end;
            int _dir=fwd?1:-1;
            if(fwd==bak||dir&&dir!=_dir){
                failure+="    Incompatible direction;\n";
                break;
            }
            int_t cur_min=fwd?index.t_start:index.t_end;
            int_t cur_max=fwd?index.t_end:index.t_start;
            if(dir){
                t_start=std::min(t_start,cur_min);
                t_end=std::max(t_end,cur_max);
            }
            else{
                dir=_dir;
                t_start=cur_min;
                t_end=cur_max;
            }
            if(index.fid==0){
                index.fid=blists.size();
                auto it=blist_index.try_emplace(dir*index.t_end,index);
                if(!it.second||!blists.emplace_back().load_barycen_structure(&mf_index,index.sid))
                    failure+="    Error loading system structure;\n";
                else if(blists.back().compatible_size()!=mn)
                    failure+="    Invalid structure;\n";
                else{
                    auto &bindex=bindices.emplace_back(index);
                    bindex.t_start*=dir;
                    bindex.t_end*=dir;
                }
            }
            else{
                auto itorb=allfiles.find(index.entry_name(false,false));
                auto itrot=allfiles.find(index.entry_name(true,false));
                if(itorb==allfiles.end()||itrot==allfiles.end()){
                    failure+="    Missing file;\n";
                    break;
                }
                size_t mid=ms.get_mid(index.sid);
                if(mid>=mn){
                    failure+="    Invalid sid;\n";
                    break;
                }
                index.fid=ephm_files.size();
                if(!ephm_index[mid].try_emplace(dir*index.t_end,index).second){
                    failure+="    Duplicated entry;\n";
                    break;
                }
                ephm_files.push_back(itorb->second);
                ephm_files.push_back(itrot->second);
            }
        }
        if(failure.size()||!dir){
            if(failure.empty())failure+="    Missing direction;\n";
            break;
        }
        
        if(offset_loaded){
            MFILE *mf_cache=&mf_offset;

            int_t bremains=blists.size();
            htl::set<int_t> vuse;
            for(const bsystem &blist:blists){
                size_t csize,bn=blist.size();
                if(fread(&csize,sizeof(csize),1,mf_cache)!=1)break;
                if(csize>bn)break;
                auto &offset_map=offset_maps.emplace_back();
                htl::set<int_t> vs;
                for(int_t i=0;i<csize;++i){
                    int_t k,v;
                    if(fread(&k,sizeof(int_t),1,mf_cache)!=1||k>=bn
                     ||fread(&v,sizeof(int_t),1,mf_cache)!=1||v==0)
                        break;
                    const barycen &b=blist[k];
                    if(b.children.size()<2)break;
                    offset_map[k].fid=v;
                    vs.insert(v);
                }
                if(offset_map.size()!=csize||vs.size()!=csize)
                    break;
                vuse.insert(vs.begin(),vs.end());
                --bremains;
            }

            bool success=false;
            size_t old_size=ephm_files.size();
            htl::map<int_t,ephemeris_entry> offset_index;
            do{
                if(bremains)break;
                ephemeris_entry index;
                int_t fremain;
                while((fremain=fread(&index,sizeof(index),1,mf_cache))==1){
                    if(index.sid!=0||vuse.erase(index.fid)!=1)
                        break;

                    auto itoffset=allfiles.find(index.offset_name(false));
                    if(itoffset==allfiles.end())
                        break;
                    auto &oindex=offset_index.try_emplace(index.fid,index).first->second;
                    oindex.fid=ephm_files.size();
                    ephm_files.push_back(itoffset->second);
                }
                success=!fremain&&vuse.empty();
                if(success)for(int_t i=0,nbsys=blists.size();i<nbsys;++i){
                    const auto &bindex=bindices[i];
                    for(auto &p:offset_maps[i]){
                        p.second=offset_index.at(p.second.fid);
                        if(!(dir*p.second.t_start<=bindex.t_start&&bindex.t_end<=dir*p.second.t_end))
                            success=false;
                    }
                }
            } while(0);
            if(!success){
                offset_maps.clear();
                ephm_files.erase(ephm_files.begin()+old_size,ephm_files.end());
                LogWarning("Warning: Invalid barycentric offset cache. Ignored.\n");
            }
        }

        htl::set<izipfile> dedup_file(ephm_files.begin(),ephm_files.end());
        dedup_file.insert(izippack::end());
        int_t n_files=ephm_files.size();
        if(dedup_file.size()!=n_files+1){
            failure+="    Mismatch ephemeris;\n";
            break;
        }

        if(dir<0)
            std::swap(t_start,t_end);

        for(size_t i=0;i<=mn;++i){
            const auto &range=i==mn?blist_index:ephm_index[i];
            int_t t_min,t_max,n_count=0;
            for(const auto &p:range){
                const auto &eid=p.second;
                if(n_count==0)
                    t_min=eid.t_start;
                else if(t_max!=eid.t_start)
                    break;
                t_max=eid.t_end;
                ++n_count;
            }
            if(!n_count||n_count!=range.size()||t_min!=t_start||t_max!=t_end){
                failure+="    Invalid coverage;\n";
                break;
            }
        }
        if(failure.size())break;

        ephm_interps.resize(n_files);
        return;
    } while(0);

    if(failure.size())
        LogError("Error loading %s:\n%s",_chname.c_str(),failure.c_str());
    close();
}

ephemeris_reader::ephemeris_reader(const char *ephemeris_path){
    cur_chid=-1;
    set_memory_limit();
    update_physics_parallel_option=0;
    update_physics=false;
    update_bsystem=false;
    update_orbits=false;
    {
        using Configs::MAX_LINESIZE;
        //load 0.zip
        //Note here we require the 1st and 4th file is checkpoint & initial file
        //see msystem::load(const char *,const char *);
        //that means user is prohibited to extract/modify/rezip the output .zips
        izippack ez(strprintf("%s.0.zip",ephemeris_path));
        int_t i_file=0,n_names=0;
        std::string failure;
        for(const izipfile &zf:ez){
            ++i_file;
            if(i_file!=1&&i_file!=4)continue;
            MFILE mf;
            zf.dumpfile(mf);
            if(i_file==1){
                if(!ms.load_checkpoint(&mf)){
                    failure+="    Invalid checkpoint;\n";
                    break;
                }
                for(const auto &mi:ms)
                    minfos.emplace_back(mi);
                continue;
            }
            char sname[MAX_LINESIZE],sid[MAX_LINESIZE];
            while(failure.empty()){
                std::string chbuf=fgetstr(&mf,true);
                size_t lsize=chbuf.size();
                if(lsize==0)break;
                if(lsize>=MAX_LINESIZE){
                    failure+="    Line too long;\n";
                    break;
                }
                if(2!=sscanf(chbuf.c_str(),"%[^\t]%s",sname,sid)){
                    failure+="    Invalid line in initial file;\n";
                    break;
                }
                int_t mid=ms.get_mid(sid);
                if(mid<0)
                    failure+="    Invalid sid;\n";
                else if(minfos[mid].name.size())
                    failure+="    Duplicate sid;\n";
                else if((minfos[mid].name=sname).empty())
                    failure+="    Empty name;\n";
                else
                    ++n_names;
            }
            break;
        }
        if(failure.empty()){
            if(ms.empty())
                failure+=ez?"    Failed to load system;\n":"    Failed to open zip file;\n";
            else if(n_names!=ms.size())
                failure+="    Failed to load names of celestials;\n";
        }
        if(failure.size()){
            LogError("Error when loading %s.0.zip:\n%s%s",
                ephemeris_path,failure.c_str(),ez?"\n"
                "  NOTE: For this program to work properly, the output zip packages\n"
                "        shall not be tampered with, i.e. do not modify or replace them\n"
                "        with an extracted and repacked version.\n\n":"");
            ms.clear();
            minfos.clear();
            return;
        }
    }

    int_t t_middle=ms.ephemeris_time();
    int_t t_fwd=t_middle,t_bak=t_middle;
    int_t n_fwd=0,n_bak=0;
    for(int dir=1;dir>=-1;dir-=2){
        bool fwd=dir>0;
        const char *fwdbak=fwd?"fwd":"bak";
        size_t cur_index=0;
        int_t &t_last=fwd?t_fwd:t_bak;
        int_t &n_chpt=fwd?n_fwd:n_bak;
        std::string failure;
        do{
            std::string chname=strprintf("%s.%llu.%s.zip",ephemeris_path,++cur_index,fwdbak);
            chapter curchpt(ms,chname);
            if(!curchpt)
                break;
            if(dir>0?curchpt.t_start>=curchpt.t_end:curchpt.t_start<=curchpt.t_end){
                failure+="    Incompatible direction;\n";
                break;
            }
            if(curchpt.t_start!=t_last){
                failure+=cur_index==1?"    Mismatch epoch;\n":"    Uncontinuous ephemeris;\n";
                break;
            }
            t_last=curchpt.t_end;
            n_chpt+=1;
            chapters.emplace_back(std::move(curchpt));
        } while(1);
        if(failure.size())
            LogWarning("Warning: %s ephemerides of %s is interrupted at %llu.%s.zip:\n%s",
                fwd?"Forward":"Backward",ephemeris_path,cur_index,fwdbak,failure.c_str());
        htl::vector<chapter> revchpts;
        for(auto it=chapters.rbegin();it!=chapters.rend();++it)
            revchpts.emplace_back(std::move(*it));
        revchpts.swap(chapters);
    }
}

int_t ephemeris_reader::interpolator_size() const{
    int_t ret=0;
    for(int_t ich:active_chapters)
        ret+=chapters[ich].interpolator_size();
    return ret;
}

bool ephemeris_reader::checkout(real t_eph){
    int_t chids=0,chide=chapters.size();
    if(!(chids<=cur_chid&&cur_chid<chide))
        cur_chid=chide/2;
    do{
        if(cur_chid<chids||chide<=cur_chid)
            return false;
        const auto &cur_chapter=chapters[cur_chid];
        bool left=!(t_eph>=cur_chapter.t_min());
        bool right=!(t_eph<=cur_chapter.t_max());
        if(left==right){
            if(left)return false;
            break;
        }
        cur_chid+=left?-1:1;
    } while(1);

    lru();
    return chapters[cur_chid].checkout(*this,t_eph);
}

//see bsystem::compose
int_t ephemeris_reader::chapter::compose_active(int_t bid){
    if(bid<0)return 0;
    auto amask=active_map[bid];
    if(!(amask&(REQUIRE_RVC|REQUIRE_RVSYSC)))return 0;
    bsystem &blist=blists[_bselect];
    barycen &b=blist[bid];
    int_t nret=1;

    mpvec cr,cv;
    if(amask&LOADED_OFFSET){
        cr=b.r_sys;
        cv=b.v_sys;
    }

    if(b.pid<0){
        b.r_sys=b.r;
        b.v_sys=b.v;
    }
    else{
        barycen &p=blist[b.pid];
        if(bid==p.gid){
            b.r_sys=b.r+blist[p.hid].r_sys;
            b.v_sys=b.v+blist[p.hid].v_sys;
        }
        else if(bid==p.hid){
            barycen &g=blist[p.gid];
            real gdm=GM_map[p.gid].GM_sys(_ft_eph)/GM_map[b.pid].GM(_ft_eph);
            b.r_sys=p.r-gdm*g.r;
            b.v_sys=p.v-gdm*g.v;
        }
        else{
            b.r_sys=b.r+p.r;
            b.v_sys=b.v+p.v;
        }
    }

    if(!(amask&REQUIRE_RVC))
        return nret;

    if(b.children.empty()){
        b.r=b.r_sys;
        b.v=b.v_sys;
    }
    else{
        if(!(amask&LOADED_OFFSET)){
            mpvec cracc(0),cvacc(0);
            for(const auto cid:b.children){
                const barycen &c=blist[cid];
                real csys=GM_map[cid].GM_sys(_ft_eph);
                cracc+=c.r*csys;
                cvacc+=c.v*csys;
            }
            real bsys=GM_map[bid].GM_sys(_ft_eph);
            cr=cracc/bsys;
            cv=cvacc/bsys;
        }
        b.r=b.r_sys-cr;
        b.v=b.v_sys-cv;
    }

    if(b.gid>=0){
        nret+=compose_active(b.hid);
        nret+=compose_active(b.gid);
    }

    for(const auto cid:b.children)
        nret+=compose_active(cid);

    return nret;
}

ephemeris_reader::active_info::active_info(active_type _type,int_t _bid,ephemeris_entry _e)
    :type(_type),bid(_bid),fid(_e.fid),t_start(_e.t_start),t_end(_e.t_end){
    fid+=_type==LOADED_ROTATION;
}

bool ephemeris_reader::chapter::checkout(ephemeris_reader &ereader,real t_eph){
    msystem &ms=ereader.ms;
    int dir=t_start<t_end?1:-1;
    int_t t_key=-int_t(dir>0?-t_eph:+t_eph);//ceil(dir*t_eph)
    if(t_key<dir*t_start||dir*t_end<t_key)
        return false;
    const int_t bsysid=blist_index.lower_bound(t_key)->second.fid;
    bsystem &blist=blists[bsysid];
    const int_t bn=blist.size(),mn=ms.size();
    if(bsysid!=_bselect){
        _bselect=-1;
        _brootid=blist.root_id();
        if(_brootid<0)return false;

        bool full_checkout=ereader.update_physics;
        auto *pcmap=offset_maps.empty()?nullptr:&offset_maps[bsysid];
        //active_type
        htl::vector<uint8_t>(bn,0).swap(active_map);
        active_composed.clear();
        for(int_t mid=0;mid<mn;++mid){
            uint8_t mconfig=ereader[mid].config;
            if(!full_checkout&&!mconfig)
                continue;
            active_composed.push_back(mid);
            int_t bid=mid;
            do{
                /*
                target: get rvc (r,v,composed) of selected blist[bid=mid]
                bid.rvc depend on:
                    bid.rv_sysc, depend on:
                        if bid.pid<0:
                            bid.rvd,
                        else if bid==pid.gid:
                            bid.rvd,
                            [hid.rv_sysc], done.
                        else if bid==pid.hid:
                            [pid.rvc], done.
                            gid.rvd,
                            gid.GM_sys,
                            pid.GM = gid.GM_sys + hid.GM_sys,
                        else (bid==pid.cid)
                            [pid.rvc], done.
                            bid.rvd,
                  if bid.children:
                    crvd(cid.rvd and cid.GM_sys and bid.GM_sys), or, bid.offset,
                */
                active_map[bid]|=REQUIRE_RVC;
                const barycen &b=blist[bid];
                if(b.pid<0)
                    active_map[bid]|=REQUIRE_RVD;
                else{
                    const barycen &p=blist[b.pid];
                    if(bid==p.gid){
                        active_map[bid]|=REQUIRE_RVD;
                        active_map[p.hid]|=REQUIRE_RVSYSC;
                    }
                    else if(bid==p.hid)
                        active_map[p.gid]|=REQUIRE_RVD;
                    else
                        active_map[bid]|=REQUIRE_RVD;
                }
                if(b.children.size()){
                    //always load offset if applicable,
                    //avoid parent position influenced by children selection.
                    if(pcmap&&pcmap->contains(bid))
                        active_map[bid]|=REQUIRE_OFFSET;
                    else for(int_t c:b.children)
                        active_map[c]|=REQUIRE_RVD;
                }
                bid=b.pid;
            } while(bid>=0);
        }

        active_files.clear();
        //load offset
        for(int_t bid=0;bid<bn;++bid){
            auto &reqmask=active_map[bid];
            if(active_map[bid]&REQUIRE_OFFSET){
                reqmask^=REQUIRE_OFFSET|LOADED_OFFSET;
                active_files.emplace_back(LOADED_OFFSET,bid,pcmap->at(bid));
            }
        }
        //load rvd & rotation
        for(int_t mid=0;mid<mn;++mid){
            auto &reqmask=active_map[blist[mid].tid];
            bool load_rvd=reqmask&REQUIRE_RVD;
            bool load_rot=full_checkout||ereader[mid].config&ROTATION;
            if(!load_rot&&!load_rvd)continue;
            const auto &eindex=ephm_index[mid].lower_bound(t_key)->second;
            if(load_rvd){
                reqmask^=REQUIRE_RVD|LOADED_RVD;
                active_files.emplace_back(LOADED_RVD,mid,eindex);
            }
            if(load_rot)
                active_files.emplace_back(LOADED_ROTATION,mid,eindex);
        }
        std::sort(active_files.begin(),active_files.end(),[this](const active_info &lhs,const active_info &rhs){
            if((lhs.type==LOADED_ROTATION)!=(rhs.type==LOADED_ROTATION))
                return rhs.type==LOADED_ROTATION;
            return ephm_files[lhs.fid]<ephm_files[rhs.fid];
        });

        int_t full_size=0,n_files=active_files.size();
        htl::vector<std::pair<int_t,int_t>> ephm_sizes;
        for(int_t i=0;i<n_files;++i)
            full_size+=ephm_sizes.emplace_back(ephm_files[active_files[i].fid].size(),i).first;
        std::sort(ephm_sizes.begin(),ephm_sizes.end());

        htl::vector<bool> ephm_expand(n_files,false);
        int_t _cache_bytes=0,used_budget=0,memory_budget=std::max(int_t(0),ereader.memory_limit/2);
        for(int_t i=0,partial_size=0;i<n_files;++i){
            auto &es=ephm_sizes[i];
            int_t cur_size=es.first;
            int_t cur_expect=partial_size+(n_files-i)*cur_size;
            if(cur_expect>memory_budget)
                break;
            used_budget=cur_expect;
            _cache_bytes=cur_size;
            partial_size+=cur_size;
        }
        memory_budget-=used_budget;
        for(int_t i=0,partial_size=0;i<n_files;++i){
            auto &es=ephm_sizes[i];
            int_t cur_size=es.first;
            int_t cur_expect=(partial_size+cur_size)*(ephemeris_compressor::max_bspline_degree+1);
            if(cur_expect>memory_budget)
                break;
            ephm_expand[es.second]=true;
            partial_size+=cur_size;
        }

        for(auto &einterp:ephm_interps)
            einterp.clear();
        for(int_t i=0;i<n_files;++i){
            const auto &ainfo=active_files[i];
            auto &einterp=ephm_interps[ainfo.fid];
            auto &efile=ephm_files[ainfo.fid];
            if(efile.fetch())
                einterp=ephemeris_interpolator(izippack::get_file(),ainfo.t_end-ainfo.t_start,
                    efile.offset(),efile.size(),_cache_bytes);
            if(!einterp)
                return false;
            if(ephm_expand[i])
                einterp.expand();
        }

        GM_map.resize(bn);
        for(mass &mi:ms)std::swap(mi.GM0,mi.dGM);
        ms.update(0,&blist);
        for(int_t i=0;i<bn;++i){
            GM_map[i].dGM=blist[i].GM;
            GM_map[i].dGM_sys=blist[i].GM_sys;
        }
        for(mass &mi:ms)std::swap(mi.GM0,mi.dGM);
        ms.update(0,&blist);
        for(int_t i=0;i<bn;++i){
            GM_map[i].GM0=blist[i].GM;
            GM_map[i].GM0_sys=blist[i].GM_sys;
        }

        _bselect=bsysid;
    }
    if(_bselect<0)
        return false;
    _interp_size=0;
    _ft_eph=fast_real(t_eph);
    ms.t_eph=t_eph;
    const bool full_update=ereader.update_bsystem||ereader.update_physics;
    if(full_update)
        ms.update(_ft_eph,&blist);
    for(const auto &ainfo:active_files){
        int_t i=ainfo.bid;
        fast_real t_offset(t_eph-real(ainfo.t_start));
        auto &einterp=ephm_interps[ainfo.fid];
        if(ainfo.type==LOADED_OFFSET){
            barycen &bi=blist[i];
            orbital_state_t orb;
            einterp(t_offset,&orb);
            bi.r_sys=orb.r;
            bi.v_sys=orb.v;
        }
        else if(ainfo.type==LOADED_RVD){
            mass &m=ms[i];
            barycen &bt=blist[blist[i].tid];
            orbital_state_t orb;
            if(!ereader.update_orbits)
                einterp(t_offset,&orb);
            else{
                massinfo &mi=ereader.minfos[i];
                mi.keplerian_GM=einterp(t_offset,&orb,mi.parameters);
                mi.state_vectors=orb;
            }
            bt.r=orb.r;
            bt.v=orb.v;
        }
        else{
            mass &m=ms[i];
            rotational_state_t rot;
            if(einterp.requires_orbital_state()){
                barycen &bt=blist[blist[i].tid];
                einterp.set_orbital_state(bt.r,bt.v);
            }
            einterp(t_offset,&rot);
            m.s.x=rot.x;
            m.s.z=rot.z;
            m.s.y=m.s.z*m.s.x;
            m.w=rot.w;
            m.GL=NAN;
        }
        _interp_size+=einterp.memory_size();
    }
    compose_active(_brootid);
    for(int_t i:active_composed){
        barycen &b=blist[i];
        mass &m=ms[i];
        m.r=b.r;
        m.v=b.v;
        if(!full_update)
            m.update(_ft_eph);
    }
    if(ereader.update_physics)
        ms.accel(ereader.update_physics_parallel_option);
    if(ereader.update_bsystem){
        ms.blist=blist;
        ms.t_update=ms.t_eph;
    }
    return true;
}

void ephemeris_reader::unload(){
    for(auto &ch:chapters)
        ch.unload();
    active_chapters.clear();
}

void ephemeris_reader::chapter::unload(){
    for(auto &einterp:ephm_interps)
        einterp.clear();
    _bselect=-1;
    _interp_size=0;
}

void ephemeris_reader::lru(){
    auto it=std::find(active_chapters.begin(),active_chapters.end(),cur_chid);
    if(it!=active_chapters.end())
        active_chapters.erase(it);
    active_chapters.push_back(cur_chid);

    int_t n_active=active_chapters.size();
    if(n_active>2){
        int_t unload_target=interpolator_size()-chapters[active_chapters.back()].interpolator_size()+memory_limit/2;
        unload_target-=memory_limit;
        int_t n_unload=0;
        for(int_t i=0;i+2<n_active;++i){
            if(unload_target<=0)
                break;
            auto &oldch=chapters[active_chapters[i]];
            unload_target-=oldch.interpolator_size();
            oldch.unload();
            ++n_unload;
        }
        active_chapters.erase(active_chapters.begin(),active_chapters.begin()+n_unload);
    }
}

void ephemeris_reader::reset_selection(){
    for(int_t ich:active_chapters)
        chapters[ich]._bselect=-1;
}

size_t ephemeris_reader::deselect_all(select_type _s){
    size_t retval=0;
    for(massinfo &mi:minfos){
        retval+=(_s&mi.config)+1>>1;
        mi.config&=~_s;
    }
    if(retval)
        reset_selection();
    return retval;
}

size_t ephemeris_reader::select_all(select_type _s){
    size_t retval=0;
    for(massinfo &mi:minfos){
        if((_s&ORBIT)>(mi.config&ORBIT)){
            mi.config|=ORBIT;
            ++retval;
        }
        if((_s&ROTATION)>(mi.config&ROTATION)){
            mi.config|=ROTATION;
            ++retval;
        }
    }
    if(retval)
        reset_selection();
    return retval;
}

ephemeris_reader::select_type ephemeris_reader::select(const massinfo &minfo,select_type _s){
    uint8_t &config=minfos[&minfo-minfos.data()].config;
    uint8_t oldconfig=config;
    if(_s&ORBIT)
        config|=ORBIT;
    if(_s&ROTATION)
        config|=ROTATION;
    if(config!=oldconfig)
        reset_selection();
    return select_type(oldconfig);
}

ephemeris_reader::select_type ephemeris_reader::deselect(const massinfo &minfo,select_type _s){
    uint8_t &config=minfos[&minfo-minfos.data()].config;
    uint8_t oldconfig=config;
    config&=~_s;
    if(config!=oldconfig)
        reset_selection();
    return select_type(oldconfig);
}

bool ephemeris_reader::set_update_physics(bool new_setting){
    bool retval=update_physics;
    update_physics=new_setting;
    if(retval!=new_setting)
        reset_selection();
    return retval;
}

