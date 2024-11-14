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

int_t ephemeris_collector::compress(std::vector<MFILE> &ephemeris_data){
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
                    LogError("Invalid barycenter list size.\n");
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

    //compress orbital data
    for(MFILE &mf:ephemeris_data){
        auto &einfo=indexmap[mf.get_name()];
        if(einfo.entry_id<0||!einfo.pmfile||einfo.rotational||einfo.substep)
            continue;
        ephemeris_entry &index=indices[einfo.entry_id];
        MFILE *psubstep=indexmap[index.entry_name(false,true)].pmfile;
        mf.publish();
        if(psubstep)psubstep->publish();

        do{
            int_t oldsize=mf.size();
            int_t samplesize=2*sizeof(vec);
            double tstart=CalcTime();
            int_t clevel=ephemeris_compressor::compress_orbital_data(mf,index.t_end-index.t_start);
            printf("%8s[%3lld, ",&index.sid,clevel);
            if(!clevel)printf("failed] in %fs\n",CalcTime()-tstart);
            else{
                auto *pheader=(ephemeris_compressor::header_base*)mf.data();
                printf("%8llu/%8llu, %.2f%% @ %6llu] : %.17e in %fs\n",
                    mf.size(),oldsize,100.*mf.size()/oldsize,(oldsize)/samplesize,
                    pheader->relative_error,CalcTime()-tstart
                );
            }
        } while(0);
    }

    //compress rotational data


    //remove substep data

    return error_count;
}
