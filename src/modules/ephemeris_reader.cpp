#include"ephemeris_reader.h"
#include<algorithm>
#include"ephemeris_generator.h"
#include"configs.h"
#include"utils/logger.h"

ephemeris_reader::chapter::chapter(msystem &ms,const std::string &_,int_t memory_budget):izippack(_){
    std::string failure;
    do{
        _interp_size=0;
        if(!izippack::operator bool())
            break;
        auto mf_files=load_central_directory();
        MFILE mf_index;
        bool index_loaded=false,chpt_loaded=false;
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
            if(failure.size()||index_loaded&&chpt_loaded)break;
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
        ephm_expand.resize(n_files);
        memory_budget=std::max(int_t(0),memory_budget);
        htl::vector<std::pair<int_t,int_t>> ephm_sizes;
        int_t full_size=0;
        for(int_t i=0;i<n_files;++i)
            full_size+=ephm_sizes.emplace_back(ephm_files[i].size(),i).first;
        std::sort(ephm_sizes.begin(),ephm_sizes.end());
        _cache_bytes=0;
        int_t used_budget=0;
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

        return;
    } while(0);

    if(failure.size())
        LogError("Error loading %s:\n%s",_.c_str(),failure.c_str());
    close();
}

ephemeris_reader::ephemeris_reader(const char *ephemeris_path,int_t memory_budget){
    cur_chid=-1;
    mem_budget=memory_budget;
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
                minfos.resize(ms.size());
                continue;
            }
            char sname[MAX_LINESIZE],sid[MAX_LINESIZE];
            while(failure.empty()){
                std::string chbuf=readline(&mf);
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
            chapter curchpt(ms,chname,memory_budget/2);
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

bool ephemeris_reader::chapter::checkout(ephemeris_reader &ereader,real t_eph){
    msystem &ms=ereader.ms;
    fast_real t_range=t_end-t_start;
    int dir=t_start<t_end?1:-1;
    int_t t_key=-int_t(dir>0?-t_eph:+t_eph);//ceil(dir*t_eph)
    if(t_key<dir*t_start||dir*t_end<t_key)
        return false;
    _interp_size=0;
    bool failed=false;
    bsystem &blist=blists[blist_index.lower_bound(t_key)->second.fid];
    ms.update(fast_real(t_eph),&blist);
    int_t bn=blist.size();
    for(int_t i=0;i<bn;++i){
        barycen &b=blist[i];
        if(b.hid<0){
            const ephemeris_entry &index=ephm_index[b.mid].lower_bound(t_key)->second;
            int_t fid=index.fid;
            fast_real t_offset(t_eph-real(index.t_start));
            orbital_state_t orb;
            rotational_state_t rot;
            for(int_t k=0;k<2;++k){
                auto &einterp=ephm_interps[fid+k];
                if(!einterp){
                    auto &efile=ephm_files[fid+k];
                    if(efile.fetch())
                        einterp=ephemeris_interpolator(get_file(),index.t_end-index.t_start,
                            efile.offset(),efile.size(),_cache_bytes);
                    if(!einterp){
                        efile=izippack::end();
                        failed=true;
                        continue;
                    }
                    if(ephm_expand[fid+k])
                        einterp.expand();
                }
                _interp_size+=einterp.memory_size();
                if(k==0){
                    if(!ereader.update_orbits)
                        einterp(t_offset,&orb);
                    else{
                        massinfo &mi=ereader.minfos[b.mid];
                        mi.keplerian_GM=einterp(t_offset,&orb,mi.parameters);
                        mi.state_vectors=orb;
                    }
                }
                else{
                    einterp.set_orbital_state(orb.r,orb.v);
                    einterp(t_offset,&rot);
                }
            }
            mass &m=ms[b.mid];
            barycen &bt=blist[b.tid];
            bt.r=orb.r;
            bt.v=orb.v;
            m.s.x=rot.x;
            m.s.z=rot.z;
            m.s.y=m.s.z*m.s.x;
            m.w=rot.w;
            m.GL=NAN;
        }
    }
    if(failed)
        ms.t_eph=NAN;
    else{
        blist.compose();
        for(int_t i=0;i<bn;++i){
            barycen &b=blist[i];
            if(b.hid<0){
                mass &m=ms[b.mid];
                m.r=b.r;
                m.v=b.v;
            }
        }
        ms.t_eph=t_eph;
        if(ereader.update_physics)
            ms.accel(ereader.update_physics_parallel_option);
    }
    if(ereader.update_bsystem){
        ms.blist=blist;
        ms.t_update=ms.t_eph;
    }
    return !failed;
}

void ephemeris_reader::unload(){
    for(auto &ch:chapters)
        ch.unload();
    active_chapters.clear();
}

void ephemeris_reader::chapter::unload(){
    for(auto &einterp:ephm_interps)
        einterp.clear();
    _interp_size=0;
}

void ephemeris_reader::lru(){
    auto it=std::find(active_chapters.begin(),active_chapters.end(),cur_chid);
    if(it!=active_chapters.end())
        active_chapters.erase(it);
    active_chapters.push_back(cur_chid);

    int_t n_active=active_chapters.size();
    if(n_active>2){
        int_t unload_target=interpolator_size()-chapters[active_chapters.back()].interpolator_size()+mem_budget/2;
        unload_target-=mem_budget;
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
