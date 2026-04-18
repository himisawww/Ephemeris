#include<algorithm>
#include"ephemeris_generator.h"
#include"configs.h"
#include"utils/zipio.h"
#include"utils/logger.h"
#include"utils/calctime.h"
#include"utils/threadpool.h"
#include"modules/ephemeris_compressor.h"

int ephemeris_collector::convert_format(const char *path,int_t fix_interval,htl::vector<const char*> *psid_subset){
    bool is_broken=false;
    std::string sop=path;
    for(int dir=1;dir>=-1;dir-=2){
        const char *fwdbak=dir>0?"fwd":"bak";

        std::string zckpt;
        size_t cur_index=0;
        MFILE *fout=mopen(sop+"."+fwdbak,MFILE_STATE::WRITE_FILE);
        if(!fout){
            LogError("Error: Cannot open %s for write!\n",(sop+"."+fwdbak).c_str());
            break;
        }
        do{
            ++cur_index;
            zckpt=strprintf("%s.%llu.%s.zip",sop.c_str(),cur_index,fwdbak);
            if(!file_exist(zckpt))break;

            msystem ms;
            izippack zp(zckpt);
            MFILE mf_index;
            int version=-1;
            htl::vector<MFILE> mf_mlist;
            for(const izipfile &zf:zp){
                const std::string &zfn=zf.name();
                if(zfn==Configs::SaveNameTimestamps){
                    zf.dumpfile(mf_index);
                    version=0;
                    continue;
                }
                if(zfn==Configs::SaveNameIndex){
                    zf.dumpfile(mf_index);
                    version=1;
                    continue;
                }
                if(zfn==Configs::SaveNameCheckpoint){
                    if(version<1)continue;
                    MFILE mf_ckpt;
                    zf.dumpfile(mf_ckpt);
                    ms.load_checkpoint(&mf_ckpt);
                    continue;
                }
                std::string fext=get_file_extension(zfn);
                if(fext=="json"||fext=="txt")
                    continue;
                zf.dumpfile(mf_mlist.emplace_back());
            }

            if(mf_mlist.empty()){
                is_broken=true;
                break;
            }

            if(version==0){//index is timestamps
                vec v5[5];
                int_t it_eph;
                int_t n_data=0;
                while(1==fread(&it_eph,sizeof(int_t),1,&mf_index)){
                    double mst_eph=it_eph;
                    bool output=n_data||dir==1&&cur_index==1;
                    if(output)fwrite(&mst_eph,sizeof(double),1,fout);
                    for(auto &mf:mf_mlist){
                        fread(&v5,sizeof(vec),5,&mf);
                        if(output)fwrite(&v5,sizeof(vec),5,fout);
                    }
                    ++n_data;
                }
            }
            else if(version==1){
                if(!fix_interval){
                    is_broken=true;
                    break;
                }
                int_t mn=ms.size();
                if(!mn){
                    is_broken=true;
                    break;
                }

                //{fname, fid}
                htl::map<std::string,int_t> fidmap;
                //{fid, ephemeris_entry of fid.dat}
                htl::map<int_t,ephemeris_entry> indices;
                //{t_end, blist over [t_start,t_end]}
                htl::map<int_t,bsystem> bss;
                //[mid]={t_end, fid over [t_start,t_end]}
                htl::vector<htl::map<int_t,int_t>> fids(mn);
                //{fid, fid.dat}
                htl::map<int_t,ephemeris_interpolator> ephm_files;
                do{
                    ephemeris_entry index;
                    if(1!=fread(&index,sizeof(index),1,&mf_index))
                        break;
                    if(index.fid==0){
                        auto &blist=bss[dir*index.t_end];
                        blist.load_barycen_structure(&mf_index,index.sid);
                    }
                    else{
                        int_t fid=indices.size();
                        if(!fidmap.insert({index.entry_name(false,false),fid}).second
                         ||!fidmap.insert({index.entry_name(true,false),fid+1}).second)
                            is_broken=true;
                        indices[fid]=index;
                        index.fid=-index.fid;
                        indices[fid+1]=index;
                    }
                } while(1);

                if(is_broken)
                    break;

                typedef vec _orb_t[2];
                typedef vec _rot_t[3];
                typedef struct{ _orb_t _orb;_rot_t _rot; } _data_t;
                int_t t_start=LLONG_MAX,t_end=LLONG_MIN,t_interval=fix_interval;
                for(auto &mf:mf_mlist){
                    auto itfid=fidmap.find(get_file_name(mf.get_name()));
                    if(itfid==fidmap.end())
                        continue;
                    int_t fid=itfid->second;
                    auto it=indices.find(fid);
                    if(it==indices.end()){
                        is_broken=true;
                        break;
                    }
                    ephemeris_entry &index=it->second;

                    t_start=std::min(t_start,dir>0?index.t_start:index.t_end);
                    t_end=std::max(t_end,dir>0?index.t_end:index.t_start);
                    int_t mid=ms.get_mid(index.sid);
                    if(index.fid>0)
                        fids[mid][dir*index.t_end]=fid;
                    ephm_files.insert({fid,ephemeris_interpolator(&mf,index.t_end-index.t_start)});
                    mf.reset();
                    index.fid=fid;
                }

                if(is_broken||t_interval<=0||(t_end-t_start)%t_interval||(t_end-t_start)/t_interval<1){
                    is_broken=true;
                    continue;
                }

                if(dir<0){
                    t_interval=-t_interval;
                    std::swap(t_start,t_end);
                }

                auto it_barycen=bss.end();
                bsystem curblist;
                htl::vector<int_t> tids(mn),bids(mn);
                htl::vector<_data_t> ephm_data(mn);
                t_end+=t_interval;
                for(int_t it_eph=t_start;it_eph!=t_end;it_eph+=t_interval){
                    double mst_eph=it_eph;
                    bool output=it_eph!=t_start||dir==1&&cur_index==1;
                    if(!output)
                        continue;
                    auto it=bss.lower_bound(dir*it_eph);
                    if(it==bss.end()){
                        is_broken=true;
                        break;
                    }
                    if(it!=it_barycen){
                        curblist=it->second;
                        it_barycen=it;
                        int_t bn=curblist.size();
                        for(int_t i=0;i<bn;++i){
                            auto &b=curblist[i];
                            if(b.hid<0){
                                tids[b.mid]=b.tid;
                                bids[b.mid]=i;
                            }
                        }
                    }
                    fwrite(&mst_eph,sizeof(double),1,fout);
                    ms.update(mst_eph,&curblist);
                    for(int_t i=0;i<mn;++i){
                        auto it_data=fids[i].lower_bound(dir*it_eph);
                        if(it_data==fids[i].end()){
                            is_broken=true;
                            continue;
                        }
                        ephemeris_entry &oindex=indices[it_data->second];
                        ephemeris_entry &rindex=indices[it_data->second+1];
                        ephemeris_interpolator &fodata=ephm_files.at(oindex.fid);
                        ephemeris_interpolator &frdata=ephm_files.at(rindex.fid);
                        if(!fodata||!frdata){
                            is_broken=true;
                            continue;
                        }
                        _data_t &v5=ephm_data[i];
                        fodata(it_eph-oindex.t_start,&v5._orb);
                        frdata.set_orbital_state(v5._orb[0],v5._orb[1]);
                        frdata(it_eph-oindex.t_start,&v5._rot);
                        barycen &b=curblist[tids[i]];
                        b.r=v5._orb[0];
                        b.v=v5._orb[1];
                    }
                    curblist.compose();
                    for(int_t i=0;i<mn;++i){
                        barycen &b=curblist[bids[i]];
                        _data_t &v5=ephm_data[i];
                        v5._orb[0]=vec(b.r);
                        v5._orb[1]=vec(b.v);
                    }
                    if(!psid_subset)
                        fwrite(ephm_data.data(),sizeof(_data_t),mn,fout);
                    else for(auto ssid:*psid_subset){
                        int_t mid=ms.get_mid(ssid);
                        if(mid<0){
                            is_broken=true;
                            break;
                        }
                        fwrite(&ephm_data[mid],sizeof(_data_t),1,fout);
                    }
                }
            }

        } while(!is_broken);
        fclose(fout);
    }
    if(is_broken)
        LogError("Error: Ephemeris %s is broken!\n",path);
    return 0;
}

int_t ephemeris_compressor::compress_work::priority() const{
    int_t sumsize=morb->size()+(msuborb?msuborb->size():0);
    if(pindex->sid)sumsize+=mrot->size()+(msubrot?msubrot->size():0);
    return sumsize;
}

void ephemeris_compressor::compress_work::run(){
    const auto &index=*pindex;
    const bool is_mass=index.sid;
    // for debug
    htl::vector<orbital_state_t> sorb,ssuborb;
    htl::vector<rotational_state_t> srot,ssubrot;
    sorb.insert(sorb.begin(),
        (orbital_state_t*)morb->data(),
        (orbital_state_t*)(morb->data()+morb->size()));
    if(msuborb)
        ssuborb.insert(ssuborb.begin(),
            (orbital_state_t*)msuborb->data(),
            (orbital_state_t*)(msuborb->data()+msuborb->size()));
    if(is_mass){
        srot.insert(srot.begin(),
            (rotational_state_t*)mrot->data(),
            (rotational_state_t*)(mrot->data()+mrot->size()));
        if(msubrot)
            ssubrot.insert(ssubrot.begin(),
                (rotational_state_t*)msubrot->data(),
                (rotational_state_t*)(msubrot->data()+msubrot->size()));
    }

    //orbital,rotational
    double trange=double(index.t_end-index.t_start);
    for(int k=0;k<1+is_mass;++k){
        MFILE *&mbase=k==0?morb:mrot;
        MFILE *&msub=k==0?msuborb:msubrot;

        newsize[k]=oldsize[k]=mbase->size();
        header_base *&pheader=pheaders[k];
        pheader=nullptr;
        int_t target_clevel=1;

        int_t clevel=
            k==0?compress_orbital_data(*mbase,trange,is_mass)
             :compress_rotational_data(*mbase,trange,morb);
        bool use_substep=false;
        if(clevel){
            newsize[k]=mbase->size();
            pheader=(header_base*)mbase->data();
            double ecrit=pheader->relative_error;
            ecrit=std::log10(ecrit/epsilon_relative_error);
            if(ecrit>0)
                target_clevel=(int_t)std::ceil(ecrit*ecrit);
        }
        if(msub&&clevel<target_clevel){
            int_t subclevel=
                k==0?compress_orbital_data(*msub,trange,is_mass)
                 :compress_rotational_data(*msub,trange,morb);
            if(subclevel){
                auto *psubheader=(header_base*)msub->data();
                use_substep=!(pheader&&pheader->relative_error<=psubheader->relative_error);
                if(use_substep){
                    newsize[k]=msub->size();
                    pheader=psubheader;
                    clevel=subclevel;
                }
            }
        }
        clevels[k]=use_substep?-clevel:clevel;
        if(!clevel)
            continue;
        if(use_substep){
            std::swap(mbase,msub);
            mbase->set_name(msub->get_name());
        }
        if(msub){
            msub->reset();
            msub->close();
        }
    }

    //debug
    ephemeris_interpolator iorb(morb,trange);
    ephemeris_interpolator irot(is_mass?mrot:nullptr,trange);
    max_r=max_v=max_xz=max_w=0;
    end_r=end_v=end_xz=end_w=0;
    max_r_relative=0;
    auto *state_error=iorb.data_format()==STATE_VECTORS?absolute_state_error:relative_state_error;
    for(size_t i=0;i<sorb.size();++i){
        orbital_state_t os;
        rotational_state_t rs;
        double t=double(i)/(sorb.size()-1)*trange;
        iorb(t,&os);
        checked_maximize(max_r_relative,state_error(&sorb[i].r,&os.r));
        checked_maximize(max_r,(sorb[i].r-os.r).norm());
        checked_maximize(max_v,(sorb[i].v-os.v).norm());
        if(is_mass){
            irot.set_orbital_state(os.r,os.v);
            irot(t,&rs);
            checked_maximize(max_w,(srot[i].w-rs.w).norm());
            checked_maximize(max_xz,(srot[i].x-rs.x).norm());
            checked_maximize(max_xz,(srot[i].z-rs.z).norm());
        }
        if(i==0||i+1==sorb.size()){
            checked_maximize(end_r,(sorb[i].r-os.r).norm());
            checked_maximize(end_v,(sorb[i].v-os.v).norm());
            if(is_mass){
                checked_maximize(end_w,(srot[i].w-rs.w).norm());
                checked_maximize(end_xz,(srot[i].x-rs.x).norm());
                checked_maximize(end_xz,(srot[i].z-rs.z).norm());
            }
        }
    }
    for(size_t i=0;i<ssuborb.size();++i){
        orbital_state_t os;
        rotational_state_t rs;
        double t=double(i)/(ssuborb.size()-1)*trange;
        iorb(t,&os);
        checked_maximize(max_r_relative,state_error(&ssuborb[i].r,&os.r));
        checked_maximize(max_r,(ssuborb[i].r-os.r).norm());
        checked_maximize(max_v,(ssuborb[i].v-os.v).norm());
        if(is_mass){
            irot.set_orbital_state(os.r,os.v);
            irot(t,&rs);
            checked_maximize(max_w,(ssubrot[i].w-rs.w).norm());
            checked_maximize(max_xz,(ssubrot[i].x-rs.x).norm());
            checked_maximize(max_xz,(ssubrot[i].z-rs.z).norm());
        }
        if(i==0||i+1==ssuborb.size()){
            checked_maximize(end_r,(ssuborb[i].r-os.r).norm());
            checked_maximize(end_v,(ssuborb[i].v-os.v).norm());
            if(is_mass){
                checked_maximize(end_w,(ssubrot[i].w-rs.w).norm());
                checked_maximize(end_xz,(ssubrot[i].x-rs.x).norm());
                checked_maximize(end_xz,(ssubrot[i].z-rs.z).norm());
            }
        }
    }
}

int_t ephemeris_compressor::compress(htl::vector<MFILE> &ephemeris_data){
    MFILE *mf_readme=nullptr,*mf_cache=nullptr;

    htl::vector<bsystem> blists;
    htl::vector<ephemeris_entry> indices;
    htl::map<std::string,MFILE*> indexmap;
    for(MFILE &mf:ephemeris_data){
        std::string namestr=get_file_name(mf.get_name());
        if(namestr==Configs::SaveNameReadme){
            mf_readme=&mf;
            continue;
        }
        if(namestr==Configs::SaveNameBarycentricOffsetIndex){
            mf_cache=&mf;
            continue;
        }
        if(namestr!=Configs::SaveNameIndex){
            indexmap[namestr]=&mf;
            continue;
        }
        mf.publish();
        ephemeris_entry index;
        while(fread(&index,sizeof(index),1,&mf)==1){
            if(index.fid==0){
                bsystem &blist=blists.emplace_back();
                if(!blist.load_barycen_structure(&mf,index.sid)){
                    LogError("\nInvalid barycenter list.\n");
                    return -1;
                }
            }
            else if(index.sid==0||index.sid>mass::max_sid){
                LogError("\nInvalid sid <%llu>.\n",index.sid);
                return -1;
            }
            else
                indices.push_back(index);
        }
    }
    htl::map<int_t,htl::set<std::pair<int_t,int_t>>> cmmap;
    if(mf_cache){
        mf_cache->publish();
        int_t bremains=blists.size();
        htl::set<int_t> vuse;
        for(const bsystem &blist:blists){
            size_t csize,bn=blist.size();
            if(fread(&csize,sizeof(csize),1,mf_cache)!=1)break;
            if(csize>bn)break;
            htl::set<int_t> ks,vs;
            for(int_t i=0;i<csize;++i){
                int_t k,v;
                if(fread(&k,sizeof(int_t),1,mf_cache)!=1||k>=bn
                 ||fread(&v,sizeof(int_t),1,mf_cache)!=1||v==0)
                    break;
                const barycen &b=blist[k];
                if(b.hid<0){
                    if(b.mid<0)
                        break;
                    cmmap[v].emplace(b.mid,-1);
                }
                else if(b.hid<bn&&b.gid<bn){
                    int_t gmid=blist[b.gid].mid;
                    if(gmid<0)
                        break;
                    cmmap[v].emplace(b.mid,gmid);
                }
                else break;
                ks.insert(k);
                vs.insert(v);
            }
            if(ks.size()!=csize||vs.size()!=csize)
                break;
            vuse.insert(vs.begin(),vs.end());
            --bremains;
        }
        bool success=false;
        size_t n_indices=indices.size();
        do{
            if(bremains)break;
            ephemeris_entry index;
            int_t fremain;
            while((fremain=fread(&index,sizeof(index),1,mf_cache))==1){
                if(index.sid!=0||vuse.erase(index.fid)!=1)
                    break;
                indices.push_back(index);
            }
            success=!fremain&&vuse.empty();
        } while(0);
        if(!success){
            indices.resize(n_indices);
            mf_cache=nullptr;
            LogWarning("Warning: Invalid barycentric offset cache. Ignored.\n");
        }
    }

    htl::vector<compress_work> tasks;
    for(const ephemeris_entry &index:indices)if(index.sid>0){
        MFILE *mrot=indexmap[index.entry_name(true,false)];
        MFILE *msubrot=indexmap[index.entry_name(true,true)];
        MFILE *morb=indexmap[index.entry_name(false,false)];
        MFILE *msuborb=indexmap[index.entry_name(false,true)];
        if(!morb||!mrot){
            LogError("\nError: Missing data for <%s>.\n",&index.sid);
            return -1;
        }

        morb->publish();
        mrot->publish();
        if(msuborb)msuborb->publish();
        if(msubrot)msubrot->publish();

        auto &work=tasks.emplace_back();
        work.pindex=&index;
        work.morb=morb;
        work.mrot=mrot;
        work.msuborb=msuborb;
        work.msubrot=msubrot;
    }
    else{
        MFILE *moffset=indexmap[index.offset_name(false)];
        MFILE *msuboffset=indexmap[index.offset_name(true)];
        if(!moffset){
            LogError("\nError: Missing offset file <%llu>.\n",index.fid);
            return -1;
        }

        moffset->publish();
        if(msuboffset)msuboffset->publish();
        auto &work=tasks.emplace_back();
        work.pindex=&index;
        work.morb=moffset;
        work.msuborb=msuboffset;
    }

    const size_t n_tasks=tasks.size();
    htl::map<uint64_t,htl::vector<size_t>> taskids_map;
    htl::vector<uint64_t> key_orders;
    for(size_t it=n_tasks;it>0;){
        const auto &w=tasks[--it];
        uint64_t sid=w.pindex->sid;
        if(sid==0)
            taskids_map[sid].push_back(it);
        else{
            size_t oldsize=taskids_map.size();
            taskids_map[sid].push_back(it);
            if(taskids_map.size()>oldsize)
                key_orders.push_back(sid);
        }
    }
    const int_t mn=key_orders.size();
    key_orders.push_back(0);
    std::reverse(key_orders.begin(),key_orders.end());
    for(bsystem &blist:blists){
        int_t cbn=blist.compatible_size();
        if(cbn!=mn){
            LogError("\nError: Uncompatible barycenter/mass system size: <%lld/%lld>.\n",cbn,mn);
            return -1;
        }
    }

    ThreadPool *pthread_pool=ThreadPool::get_thread_pool();
    if(pthread_pool){
        htl::vector<std::pair<int_t,void*>> priorities;
        for(auto &w:tasks)
            priorities.push_back({-w.priority(),&w});
        std::sort(priorities.begin(),priorities.end());
        ThreadPool::TaskGroup task_group;
        for(auto &p:priorities)
            pthread_pool->add_task(do_compress_work,p.second,&task_group);
        auto callback=[&](){
            LogInfo(" Compressing ephemerides. [%llu/%llu]\r",n_tasks-task_group.load(),n_tasks);
        };
        pthread_pool->wait_for_all(&task_group,callback,1.5);
    }
    else{
        for(size_t i=0;i<n_tasks;++i){
            LogInfo(" Compressing ephemerides. [%llu/%llu]\r",i+1,n_tasks);
            tasks[i].run();
        }
    }

    LogInfo("\n");

    int_t error_count=0;
    for(const auto &w:tasks){
        error_count+=!w.pheaders[0];
        if(w.pindex->sid)
            error_count+=!w.pheaders[1];
    }

    do{
        if(!mf_readme){
            LogWarning("Warning: Missing %s in ephemerides pack.\n",Configs::SaveNameReadme);
            break;
        }
        mf_readme->publish();
        const char *pdata=(const char *)mf_readme->data();
        const std::string readmestr(pdata,pdata+mf_readme->size());
        const char search_pattern[]="\n"
            "  Object List (index & sid):  \n";
        const char *plocate=strstr(readmestr.c_str(),search_pattern);
        if(!plocate){
            LogWarning("Warning: Unrecognized %s.\n",Configs::SaveNameReadme);
            break;
        }
        mf_readme->reset();
        mf_readme->set_name(Configs::SaveNameReadme);
        fwrite(readmestr.c_str(),1,plocate-readmestr.c_str(),mf_readme);
        fprintf(mf_readme,"\n"
            "Compressed Format:\n"
            "       data_file : method(degree, segments)\n"
            "                   [(substep*)compress_level, size/original_size, ratio @ original_sample_count in [t_start, t_end]]\n"
            "                 : relative_error [endpoints:max_state_error, endpoints:max_rate_error]:\n");
        int_t mi=0;
        const char *fexts[3]={Configs::SaveBarycentricOffsetDataExtension,Configs::SaveOrbitalDataExtension,Configs::SaveRotationalDataExtension};
        for(uint64_t sid:key_orders){
            const auto &v=taskids_map.at(sid);
            const bool is_mass=sid>0;
            if(!is_mass){
                fprintf(mf_readme,"\n"
                    "Barycentric Offset List:\n"
                    "<sid[-companion], ...> :\n");
            }
            else{
                fprintf(mf_readme,"%s[%7s]%7lld :\n",
                    mi?"":"\n"
                    "Object List:\n"
                    "[    sid]  index :\n",(char*)&sid,mi);
                ++mi;
            }
            for(auto it=v.rbegin();it!=v.rend();++it){
                const auto &w=tasks[*it];
                const auto &index=*w.pindex;
                for(int k=0;k<1+is_mass;++k){
                    auto *pheader=w.pheaders[k];
                    const char *fext=fexts[k+is_mass];
                    if(!is_mass){
                        std::string ssids;
                        for(const auto &mip:cmmap[index.fid]){
                            if(!ssids.empty())ssids+=", ";
                            ssids+=(char*)&key_orders[mip.first+1];
                            if(mip.second<0)continue;
                            ssids+='-';
                            ssids+=(char*)&key_orders[mip.second+1];
                        }
                        fprintf(mf_readme,"<%-15s :\n",(ssids+='>').c_str());
                    }
                    if(!pheader){
                        fprintf(mf_readme,"%12lld%s : %18s(FAILED)\n",
                            index.fid,fext,
                            format_name(ephemeris_format::NONE));
                        LogError("Error: Failed to compress ephemeris file <%s>\n",get_file_name((k?w.mrot:w.morb)->get_name()).c_str());
                    }
                    else{
                        bool use_substep=w.clevels[k]<0;
                        int_t clevel=std::abs(w.clevels[k]);
                        size_t newsize=w.newsize[k],oldsize=w.oldsize[k];
                        // verbose info
                        int_t samplesize=(2+k)*sizeof(vec);
                        fprintf(mf_readme,"%12lld%s : %18s(%d, %lld)\n",
                            index.fid,fext,
                            format_name(ephemeris_format(~pheader->uformat)),
                            (int)pheader->degree,pheader->n);
                        fprintf(mf_readme,
                            "                   [%c%3lld, %8llu/%8llu, %6.2f%% @ %6llu in [%lld, %lld]]\n",
                            use_substep[" *"],clevel,
                            newsize,oldsize,100.*newsize/oldsize,oldsize/samplesize,
                            index.t_start,index.t_end);
                        double err_rel=pheader->relative_error;
                        //this maximize shall not change err_rel by much,
                        // otherwise suspect erroneous subsystem link, recording inconsistent mainstep/substep data.
                        checked_maximize(err_rel,k==0?w.max_r_relative:w.max_xz/2);
                        fprintf(mf_readme,
                            k==0?"                 : %.3e [%.3e:%.3e m  , %.3e:%.3e m/s  ]\n"
                                :"                 : %.3e [%.3e:%.3e rad, %.3e:%.3e rad/s]\n",
                            err_rel,
                            k==0?w.end_r:w.end_xz,
                            k==0?w.max_r:w.max_xz,
                            k==0?w.end_v:w.end_w,
                            k==0?w.max_v:w.max_w);
                        if(!(err_rel<relative_error_warning_threshold))
                            LogWarning("Warning: Relative fit error (%.3e) too large for <%s>\n",
                                err_rel,get_file_name((k?w.mrot:w.morb)->get_name()).c_str());
                    }
                }
            }
        }
    } while(0);
    return error_count;
}
