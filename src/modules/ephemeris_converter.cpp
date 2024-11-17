#include"ephemeris_generator.h"
#include"configs.h"
#include"utils/zipio.h"
#include"utils/logger.h"
#include"utils/calctime.h"
#include"modules/ephemeris_compressor.h"

int ephemeris_collector::convert_format(const char *path){
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
                continue;
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
                ephemeris_collector ephc(ms);
                int_t mn=ms.size();
                if(!mn){
                    is_broken=true;
                    continue;
                }

                //{fname, fid}
                std::map<std::string,int_t> fidmap;
                //{fid, ephemeris_entry of fid.dat}
                std::map<int_t,ephemeris_entry> indices;
                //{t_end, blist over [t_start,t_end]}
                std::map<int_t,std::vector<barycen>> bss;
                //[mid]={t_end, fid over [t_start,t_end]}
                std::vector<std::map<int_t,int_t>> orb_fids(mn);
                std::vector<int_t> rot_fids(mn);
                //{fid, fid.dat}
                std::map<int_t,MFILE*> ephm_files;
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
                        int_t fid=fidmap.size();
                        std::string fname=index.fid>0
                            ?strprintf("%s.%lld%s",&index.sid,index.fid,Configs::SaveOrbitalDataExtension)
                            :strprintf("%s%s",&index.sid,Configs::SaveRotationalDataExtension);
                        if(!fidmap.insert({fname,fid}).second)
                            is_broken=true;
                        indices[fid]=index;
                    }
                } while(1);

                typedef vec _orb_t[2];
                typedef vec _rot_t[3];
                typedef struct{ _orb_t _orb;_rot_t _rot; } _data_t;
                int_t t_start,t_end,t_interval=0;
                for(auto &mf:mf_mlist){
                    auto itfid=fidmap.find(mf.get_name());
                    if(itfid==fidmap.end())
                        continue;
                    int_t fid=itfid->second;
                    auto it=indices.find(fid);
                    if(it==indices.end()){
                        is_broken=true;
                        continue;
                    }
                    ephemeris_entry &index=it->second;
                    bool is_orb=index.fid>0;
                    mf.seek(0,SEEK_END);
                    int_t fsize=mf.tell()/(is_orb?sizeof(_orb_t):sizeof(_rot_t));
                    mf.seek(0,SEEK_SET);
                    if(fsize<=1){
                        is_broken=true;
                        continue;
                    }
                    int_t finterval=(index.t_end-index.t_start)/(fsize-1);
                    if(t_interval==0){
                        t_interval=finterval;
                        t_start=index.t_start;
                    }
                    else{
                        if(t_interval*(fsize-1)!=index.t_end-index.t_start){
                            is_broken=true;
                            continue;
                        }
                        t_end=index.t_end;
                    }
                    int_t mid=ms.get_mid(index.sid);
                    if(is_orb)
                        orb_fids[mid][dir*index.t_end]=fid;
                    else
                        rot_fids[mid]=fid;
                    ephm_files[fid]=&mf;
                    index.fid=fid;
                }

                if(t_interval==0||(t_end-t_start)%t_interval||(t_end-t_start)/t_interval<1){
                    is_broken=true;
                    continue;
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
                        auto it_data=orb_fids[i].lower_bound(dir*it_eph);
                        if(it_data==orb_fids[i].end()){
                            is_broken=true;
                            continue;
                        }
                        ephemeris_entry &oindex=indices[it_data->second];
                        ephemeris_entry &rindex=indices[rot_fids[i]];
                        MFILE *fodata=ephm_files[oindex.fid];
                        MFILE *frdata=ephm_files[rindex.fid];
                        _data_t &v5=ephm_data[i];
                        fodata->seek((it_eph-oindex.t_start)/t_interval*sizeof(_orb_t),SEEK_SET);
                        frdata->seek((it_eph-rindex.t_start)/t_interval*sizeof(_rot_t),SEEK_SET);
                        if(1!=fread(v5._orb,sizeof(_orb_t),1,fodata)
                         ||1!=fread(v5._rot,sizeof(_rot_t),1,frdata)){
                            is_broken=true;
                            continue;
                        }
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

int_t ephemeris_compressor::compress(std::vector<MFILE> &ephemeris_data){
    int_t error_count=0;
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

    std::vector<ephemeris_entry> indices;
    std::map<std::string,entry_info> indexmap;
    for(MFILE &mf:ephemeris_data){
        const std::string &namestr=mf.get_name();
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

            int_t oldsize=mbase->size();
            int_t newsize=oldsize;
            header_base *pheader=nullptr;
            int_t target_clevel=1;
            double tstart=CalcTime();
            int_t clevel=
                k==0?compress_orbital_data(*mbase,trange)
                 :compress_rotational_data(*mbase,trange,morb);
            bool use_substep=false;
            if(clevel){
                newsize=mbase->size();
                pheader=(header_base*)mbase->data();
                double ecrit=pheader->relative_error;
                ecrit=std::log10(ecrit/epsilon_relative_error);
                if(ecrit>0)
                    target_clevel=(int_t)std::ceil(ecrit*ecrit);
            }
            if(msub&&clevel<target_clevel){
                clevel=
                    k==0?compress_orbital_data(*msub,trange)
                     :compress_rotational_data(*msub,trange,morb);
                if(clevel){
                    auto *psubheader=(header_base*)msub->data();
                    use_substep=!(pheader&&pheader->relative_error<=psubheader->relative_error);
                    if(use_substep){
                        newsize=msub->size();
                        pheader=psubheader;
                    }
                }
            }
            if(!clevel){
                ++error_count;
                LogWarning("\nWarning: Cannot compress %s data for <%s>, ignored.\n",
                    k==0?"orbital":"rotational",&index.sid);
                continue;
            }
            if(use_substep){
                std::swap(mbase,msub);
                mbase->set_name(index.entry_name(k,false));
            }
            if(msub){
                msub->reset();
                msub->close();
            }
            // verbose info
            int_t samplesize=(2+k)*sizeof(vec);
            LogInfo("%8s %c:%d[%c%3lld, %8llu/%8llu, %6.2f%% @ %6llu] : %.17e in %fs\n",
                &index.sid,k["or"],(int)ephemeris_format(~pheader->uformat),use_substep[" *"],clevel,
                newsize,oldsize,100.*newsize/oldsize,oldsize/samplesize,
                pheader->relative_error,CalcTime()-tstart
            );
        }

        //debug
        ephemeris_interpolator iorb(morb,trange);
        ephemeris_interpolator irot(mrot,trange);
        double max_r=0,max_v=0,max_xz=0,max_w=0;
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
        }
        printf(" max_errs:[%.6e, %.6e, %.6e, %.6e]\n",max_r,max_v,max_xz,max_w);
        if(ssuborb.empty())continue;
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
        }
        printf("*max_errs:[%.6e, %.6e, %.6e, %.6e]\n",max_r,max_v,max_xz,max_w);
    }
    return error_count;
}
