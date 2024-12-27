#include<algorithm>
#include"ephemeris_generator.h"
#include"configs.h"
#include"utils/zipio.h"
#include"utils/logger.h"
#include"utils/calctime.h"
#include"utils/threadpool.h"
#include"modules/ephemeris_compressor.h"

int ephemeris_collector::convert_format(const char *path,int_t fix_interval){
    bool is_broken=false;
    std::string sop=path;
    for(int dir=1;dir>=-1;dir-=2){
        const char *fwdbak=dir>0?"fwd":"bak";

        std::string zckpt;
        size_t cur_index=0;
        MFILE *fout=mopen(sop+"."+fwdbak,MFILE_STATE::WRITE_FILE);
        do{
            ++cur_index;
            zckpt=strprintf("%s.%llu.%s.zip",sop.c_str(),cur_index,fwdbak);
            if(!file_exist(zckpt))break;

            msystem ms;
            izippack zp(zckpt);
            MFILE mf_index;
            int version=-1;
            std::vector<MFILE> mf_mlist;
            int_t mi=-Configs::ExportHeaderCount;
            for(const auto &zf:zp){
                std::string zfn=zf.name();
                if(mi>=0){
                    mf_mlist.resize(mi+1);
                    zf.dumpfile(mf_mlist[mi]);
                }
                else if(zfn==Configs::SaveNameTimestamps){
                    zf.dumpfile(mf_index);
                    version=0;
                }
                else if(zfn==Configs::SaveNameIndex){
                    zf.dumpfile(mf_index);
                    version=1;
                }
                else if(zfn==Configs::SaveNameCheckpoint){
                    MFILE mf_ckpt;
                    zf.dumpfile(mf_ckpt);
                    ms.load_checkpoint(&mf_ckpt);
                }
                ++mi;
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
                ephemeris_collector ephc(ms);
                int_t mn=ms.size();
                if(!mn){
                    is_broken=true;
                    break;
                }

                //{fname, fid}
                std::map<std::string,int_t> fidmap;
                //{fid, ephemeris_entry of fid.dat}
                std::map<int_t,ephemeris_entry> indices;
                //{t_end, blist over [t_start,t_end]}
                std::map<int_t,std::vector<barycen>> bss;
                //[mid]={t_end, fid over [t_start,t_end]}
                std::vector<std::map<int_t,int_t>> fids(mn);
                //{fid, fid.dat}
                std::map<int_t,ephemeris_interpolator> ephm_files;
                do{
                    ephemeris_entry index;
                    if(1!=fread(&index,sizeof(index),1,&mf_index))
                        break;
                    if(index.fid==0){
                        auto &blist=bss[dir*index.t_end];
                        blist.resize(index.sid);
                        msystem::load_barycen_structure(&mf_index,blist);
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
                    auto itfid=fidmap.find(mf.get_name());
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
                std::vector<int_t> tids(mn),bids(mn);
                std::vector<_data_t> ephm_data(mn);
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
                        ephc.blist=it->second;
                        barycen::update_barycens(ms,ephc.blist);
                        it_barycen=it;
                        int_t bn=ephc.blist.size();
                        for(int_t i=0;i<bn;++i){
                            auto &b=ephc.blist[i];
                            if(b.hid<0){
                                tids[b.mid]=b.tid;
                                bids[b.mid]=i;
                            }
                        }
                    }
                    fwrite(&mst_eph,sizeof(double),1,fout);
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
                        barycen &b=ephc.blist[tids[i]];
                        b.r=v5._orb[0];
                        b.v=v5._orb[1];
                    }
                    barycen::compose(ephc.blist);
                    for(int_t i=0;i<mn;++i){
                        barycen &b=ephc.blist[bids[i]];
                        _data_t &v5=ephm_data[i];
                        v5._orb[0]=vec(b.r);
                        v5._orb[1]=vec(b.v);
                    }
                    fwrite(ephm_data.data(),sizeof(_data_t),mn,fout);
                }
            }

        } while(1);
        fclose(fout);
    }
    if(is_broken)
        LogError("Error: Ephemeris %s is broken!\n",path);
    return 0;
}

int_t ephemeris_compressor::compress_work::priority() const{
    int_t sumsize=morb->size()+mrot->size();
    if(msuborb)sumsize+=msuborb->size();
    if(msubrot)sumsize+=msubrot->size();
    return sumsize;
}

void ephemeris_compressor::compress_work::run(){
    const auto &index=*pindex;

    // for debug
    std::vector<orbital_state_t> sorb,ssuborb;
    std::vector<rotational_state_t> srot,ssubrot;
    sorb.insert(sorb.begin(),
        (orbital_state_t*)morb->data(),
        (orbital_state_t*)(morb->data()+morb->size()));
    srot.insert(srot.begin(),
        (rotational_state_t*)mrot->data(),
        (rotational_state_t*)(mrot->data()+mrot->size()));
    if(msuborb)
        ssuborb.insert(ssuborb.begin(),
            (orbital_state_t*)msuborb->data(),
            (orbital_state_t*)(msuborb->data()+msuborb->size()));
    if(msubrot)
        ssubrot.insert(ssubrot.begin(),
            (rotational_state_t*)msubrot->data(),
            (rotational_state_t*)(msubrot->data()+msubrot->size()));

    //orbital,rotational
    double trange=double(index.t_end-index.t_start);
    for(int k=0;k<2;++k){
        MFILE *&mbase=k==0?morb:mrot;
        MFILE *&msub=k==0?msuborb:msubrot;

        newsize[k]=oldsize[k]=mbase->size();
        header_base *&pheader=pheaders[k];
        pheader=nullptr;
        int_t target_clevel=1;

        int_t clevel=
            k==0?compress_orbital_data(*mbase,trange)
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
                k==0?compress_orbital_data(*msub,trange)
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
            mbase->set_name(index.entry_name(k,false));
        }
        if(msub){
            msub->reset();
            msub->close();
        }
    }

    //debug
    ephemeris_interpolator iorb(morb,trange);
    ephemeris_interpolator irot(mrot,trange);
    max_r=max_v=max_xz=max_w=0;
    end_r=end_v=end_xz=end_w=0;
    for(size_t i=0;i<sorb.size();++i){
        orbital_state_t os;
        rotational_state_t rs;
        double t=double(i)/(sorb.size()-1)*trange;
        iorb(t,&os);
        irot.set_orbital_state(os.r,os.v);
        irot(t,&rs);
        checked_maximize(max_r,(sorb[i].r-os.r).norm());
        checked_maximize(max_v,(sorb[i].v-os.v).norm());
        checked_maximize(max_w,(srot[i].w-rs.w).norm());
        checked_maximize(max_xz,(srot[i].x-rs.x).norm());
        checked_maximize(max_xz,(srot[i].z-rs.z).norm());
        if(i==0||i+1==sorb.size()){
            checked_maximize(end_r,(sorb[i].r-os.r).norm());
            checked_maximize(end_v,(sorb[i].v-os.v).norm());
            checked_maximize(end_w,(srot[i].w-rs.w).norm());
            checked_maximize(end_xz,(srot[i].x-rs.x).norm());
            checked_maximize(end_xz,(srot[i].z-rs.z).norm());
        }
    }
    for(size_t i=0;i<ssuborb.size();++i){
        orbital_state_t os;
        rotational_state_t rs;
        double t=double(i)/(ssuborb.size()-1)*trange;
        iorb(t,&os);
        irot.set_orbital_state(os.r,os.v);
        irot(t,&rs);
        checked_maximize(max_r,(ssuborb[i].r-os.r).norm());
        checked_maximize(max_v,(ssuborb[i].v-os.v).norm());
        checked_maximize(max_w,(ssubrot[i].w-rs.w).norm());
        checked_maximize(max_xz,(ssubrot[i].x-rs.x).norm());
        checked_maximize(max_xz,(ssubrot[i].z-rs.z).norm());
        if(i==0||i+1==ssuborb.size()){
            checked_maximize(end_r,(ssuborb[i].r-os.r).norm());
            checked_maximize(end_v,(ssuborb[i].v-os.v).norm());
            checked_maximize(end_w,(ssubrot[i].w-rs.w).norm());
            checked_maximize(end_xz,(ssubrot[i].x-rs.x).norm());
            checked_maximize(end_xz,(ssubrot[i].z-rs.z).norm());
        }
    }
}

int_t ephemeris_compressor::compress(std::vector<MFILE> &ephemeris_data){
    struct entry_info{
        int_t entry_id;
        MFILE *pmfile;
        bool rotational,substep;

        entry_info():entry_id(-1),pmfile(nullptr){}
        void set_info(int_t _entry_id,bool _rotational,bool _substep){
            entry_id=_entry_id;
            rotational=_rotational;
            substep=_substep;
        }
    };

    MFILE *mf_readme=nullptr;

    std::vector<ephemeris_entry> indices;
    std::map<std::string,entry_info> indexmap;
    for(MFILE &mf:ephemeris_data){
        const std::string &namestr=mf.get_name();
        if(namestr==Configs::SaveNameReadme){
            mf_readme=&mf;
            continue;
        }
        if(namestr!=Configs::SaveNameIndex){
            indexmap[namestr].pmfile=&mf;
            continue;
        }
        mf.publish();
        ephemeris_entry index;
        while(fread(&index,sizeof(index),1,&mf)==1){
            if(index.fid==0){
                if(index.sid>2*ephemeris_data.size()){
                    LogError("\nInvalid barycenter list size.\n");
                    return -1;
                }
                std::vector<barycen> blist(index.sid);
                msystem::load_barycen_structure(&mf,blist);
            }
            else{
                int_t entry_id=indices.size();
                indexmap[index.entry_name(true,false)].set_info(entry_id,true,false);
                indexmap[index.entry_name(false,false)].set_info(entry_id,false,false);
                indexmap[index.entry_name(true,true)].set_info(entry_id,true,true);
                indexmap[index.entry_name(false,true)].set_info(entry_id,false,true);
                indices.push_back(index);
            }
        }
    }

    std::vector<compress_work> tasks;
    for(const ephemeris_entry &index:indices){
        MFILE *mrot=indexmap[index.entry_name(true,false)].pmfile;
        MFILE *msubrot=indexmap[index.entry_name(true,true)].pmfile;
        MFILE *morb=indexmap[index.entry_name(false,false)].pmfile;
        MFILE *msuborb=indexmap[index.entry_name(false,true)].pmfile;
        if(!morb||!mrot){
            LogError("\nError: Missing data for <%s>\n",&index.sid);
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

    const size_t n_tasks=tasks.size();
    ThreadPool *pthread_pool=ThreadPool::get_thread_pool();
    if(pthread_pool){
        std::vector<std::pair<int_t,void*>> priorities;
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

    std::map<uint64_t,std::vector<size_t>> compress_info_map;
    std::vector<uint64_t> key_orders;
    for(size_t it=n_tasks;it>0;){
        const auto &w=tasks[--it];
        error_count+=!w.pheaders[0];
        error_count+=!w.pheaders[1];
        uint64_t sid=w.pindex->sid;
        size_t oldsize=compress_info_map.size();
        compress_info_map[sid].push_back(it);
        if(compress_info_map.size()>oldsize)
            key_orders.push_back(sid);
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
            "Object List:\n"
            "[    sid]  index :\n"
            "       data_file : method(degree, segments)\n"
            "                   [(substep*)compress_level, size/original_size, ratio @ original_sample_count in [t_start, t_end]]\n"
            "                 : relative_error [endpoints:max_state_error, endpoints:max_rate_error]:\n");
        int_t mi=0;
        for(auto sit=key_orders.rbegin();sit!=key_orders.rend();++sit){
            uint64_t sid=*sit;
            const auto &v=compress_info_map.at(sid);
            fprintf(mf_readme,"[%7s]%7lld :\n",(char*)&sid,mi);
            ++mi;
            for(auto it=v.rbegin();it!=v.rend();++it){
                const auto &w=tasks[*it];
                const auto &index=*w.pindex;
                for(int k=0;k<2;++k){
                    auto *pheader=w.pheaders[k];
                    if(!pheader){
                        fprintf(mf_readme,"%12lld%s : %18s(FAILED)\n",
                            index.fid,k==0?Configs::SaveOrbitalDataExtension:Configs::SaveRotationalDataExtension,
                            format_name(ephemeris_format::NONE));
                        LogError("Error: Failed to compress ephemeris file <%s>\n",index.entry_name(k==1,false).c_str());
                    }
                    else{
                        bool use_substep=w.clevels[k]<0;
                        int_t clevel=std::abs(w.clevels[k]);
                        size_t newsize=w.newsize[k],oldsize=w.oldsize[k];
                        // verbose info
                        int_t samplesize=(2+k)*sizeof(vec);
                        fprintf(mf_readme,"%12lld%s : %18s(%d, %lld)\n",
                            index.fid,k==0?Configs::SaveOrbitalDataExtension:Configs::SaveRotationalDataExtension,
                            format_name(ephemeris_format(~pheader->uformat)),
                            (int)pheader->degree,pheader->n);
                        fprintf(mf_readme,
                            "                   [%c%3lld, %8llu/%8llu, %6.2f%% @ %6llu in [%lld, %lld]]\n",
                            use_substep[" *"],clevel,
                            newsize,oldsize,100.*newsize/oldsize,oldsize/samplesize,
                            index.t_start,index.t_end);
                        fprintf(mf_readme,
                            k==0?"                 : %.3e [%.3e:%.3e m  , %.3e:%.3e m/s  ]\n"
                                :"                 : %.3e [%.3e:%.3e rad, %.3e:%.3e rad/s]\n",
                            pheader->relative_error,
                            k==0?w.end_r:w.end_xz,
                            k==0?w.max_r:w.max_xz,
                            k==0?w.end_v:w.end_w,
                            k==0?w.max_v:w.max_w);
                        if(!(pheader->relative_error<relative_error_warning_threshold))
                            LogWarning("Warning: Relative fit error (%.3e) too large for <%s>\n",
                                pheader->relative_error,index.entry_name(k==1,false).c_str());
                    }
                }
            }
        }
    } while(0);
    return error_count;
}
