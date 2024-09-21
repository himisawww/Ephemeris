#include"ephemeris_generator.h"
#include<queue>
#include"configs.h"
#include"utils/zipio.h"
#include"utils/calctime.h"
#include"utils/logger.h"

std::mutex ephemeris_generator::io_mutex;

int ephemeris_generator::make_ephemeris(int dir){
    using namespace Configs;
    constexpr int FAILED_FLAG=0xff;

    msystem ms;
    std::string sop=op;
    std::string ickpt=sop+".0.zip";
    const char *fwdbak=dir>0?"fwd":"bak";
    size_t cur_index=1;
    bool success=false;

    //loading process
    io_mutex.lock();
  do{
    if(fix_dir==FAILED_FLAG)
        break;

    if(file_exist(ickpt)){
        // load file from ckpt
        std::string zckpt;

        do{
            zckpt=strprintf("%s.%llu.%s.zip",sop.c_str(),cur_index,fwdbak);
            if(!file_exist(zckpt))break;
            ickpt.swap(zckpt);
            ++cur_index;
        } while(1);
        //ickpt exists
        izippack zp(ickpt);
        for(const auto &zf:zp){
            if(zf.name()==SaveNameCheckpoint){
                MFILE mf;
                zf.dumpfile(mf);
                ms.load_checkpoint(&mf);
                break;
            }
        }
        LogInfo("Loaded %lld bodies from checkpoint %s\n",ms.mlist.size(),ickpt.c_str());
    }
    else if(ip){
        ms.load(ip,ickpt.c_str());
    }

    if(!ms.mlist.size()){
        LogCritical(
            "Failed to Load System.\n"
            "The program will exit.\n");
        fix_dir=FAILED_FLAG;
        break;
    }

    success=true;
  } while(0);
    io_mutex.unlock();

    if(!success)
        return -1;

    bool use_cpu     =ms.integrator==int_t(msystem::     CPU_RK12);
    bool use_gpu     =ms.integrator==int_t(msystem::     GPU_RK12);
    bool use_combined=ms.integrator==int_t(msystem::COMBINED_RK12);
    
    fast_real dt=use_combined?ms.combined_delta_t:ms.delta_t;
    int_t jsize=std::round(ms.data_cadence/dt);
    if(jsize<1)jsize=1;
    dt=ms.data_cadence/jsize;

    int_t n_combine;
    if(use_combined){
        n_combine=std::round(dt/ms.delta_t);
        if(n_combine<1)n_combine=1;
        dt/=n_combine;
    }

    if(dt!=std::round(dt)){
        LogError(
            "Non-integer Delta_t(%fs):\n    This may cause round-off errors on timestamps and is disallowed.\n",
            dt);
        return -3;
    }
    
    dt*=dir;

    int_t isize=t_years*Constants::year/ms.data_cadence;
    int_t iunit=ms.max_ephm_length/ms.data_cadence;

    const int_t min_iunit=4096;
    if(isize<1){
        LogWarning("Nothing to do.\n");
        return 0;
    }
    if(iunit<min_iunit){
        LogWarning(
            "Too less Data Points per File(%lld):\n    Changed to %lld.\n",
            iunit,min_iunit);
        iunit=min_iunit;
    }

    int_t time_idx=0;

    msystem::thread_local_pool_alloc();
    do{
        if(time_idx+iunit>isize)iunit=isize-time_idx;

        //prepare mem files to save
        int_t mn=ms.mlist.size();
        int_t max_dp_perfile=1+iunit;
        std::vector<MFILE> zms(ExportHeaderCount);

        int_t t_eph_start=ms.ephemeris_time();

        bool barycen_updated=false;
        ephemeris_collector ephc(ms);
        //ephemeris integration
        for(int_t i=0;i<=iunit;++i){
            if(i>0){
                if(use_cpu     )ms.integrate         (dt,          jsize,0);
           else if(use_gpu     )ms.integrate         (dt,          jsize,1);
           else if(use_combined)ms.combined_integrate(dt,n_combine,jsize,1);
                barycen_updated=ms.analyse();
            }

            ephc.record();

            if(i==iunit)
                ephc.extract(zms,true);
            else if(barycen_updated)
                ephc.extract(zms,false);

            int_t mst_eph=ms.ephemeris_time();
            
            io_mutex.lock();
            if(barycen_updated)LogInfo(
                "\nInfo: At Ephemeris Time: %lld s\n"
                "   System orbital structure is updated.\n",
                mst_eph);
            if(fix_dir||dir>0){
                static double s=CalcTime();
                static double oldt=-INFINITY;
                static int_t skip_count=0,skip_size=1;
                if(++skip_count>=skip_size||i==iunit){
                    double yr=mst_eph/Constants::year;
                    double t=CalcTime();
                    LogInfo(" Integrating. t_eph: %.6fyr, time: %.6fs\r",yr,t-s);
                    int_t new_skip_size=std::round(skip_size/(t-oldt));
                    if(new_skip_size<=0)new_skip_size=1;
                    if(new_skip_size>2*skip_size)new_skip_size=2*skip_size;
                    skip_size=new_skip_size;
                    skip_count=0;
                    oldt=t;
                    if(time_idx+i==isize&&!fix_dir)fix_dir=-dir;
                }
            }
            io_mutex.unlock();
        }

        int_t t_eph_end=ms.ephemeris_time();

        MFILE &mf_index=zms[0];
        MFILE &mf_ckpt=zms[1];
        MFILE &mf_readme=zms[2];
        MFILE &mf_struct=zms[3];

        //prepare checkpoint/readme/structure files
        ms.save_checkpoint(&mf_ckpt);

        fprintf(&mf_readme,
            "Calculated & Generated by Ephemeris Integrator %s\n"
            "   Github: https://github.com/himisawww/Ephemeris \n"
            "   Author: %s\n\n",
            VersionString,AuthorName);
        fprintf(&mf_readme,
            "Number of Objects:  %lld\n"
            "   Time Range (s): [%lld, %lld]\n"
            "   Time Step  (s):  %lld\n\n",
            mn,t_eph_start,t_eph_end,(t_eph_end-t_eph_start)/iunit);
        fprintf(&mf_readme,
            "Checkpoint Format: version %lld\n"
            "      Data Format: version %lld\n\n"
            "  Object List (index & sid):  \n",
            CheckPointVersion,DataFormatVersion);

        ms.print_structure(&mf_struct);
        
        mf_index.set_name(SaveNameIndex);
        mf_ckpt.set_name(SaveNameCheckpoint);
        mf_readme.set_name(SaveNameReadme);
        mf_struct.set_name(SaveNameStructure);

        for(int_t mi=0;mi<mn;++mi){
            fprintf(&mf_readme,
                "%12lld : %s\n",mi,(char*)&ms[mi].sid
            );
        }

        //save .zip
        std::string zckpt;
        zckpt=strprintf("%s.%llu.%s.zip",sop.c_str(),cur_index,fwdbak);
        ++cur_index;

        io_mutex.lock();
        {
            ozippack zp(zckpt);
            zp.swap(zms);
            LogInfo("\nSaving ephemeris & checkpoint %s\n",zckpt.c_str());
        }
        io_mutex.unlock();
        
        time_idx+=iunit;
    }while(time_idx<isize);
    msystem::thread_local_pool_free();

    return 0;
}

ephemeris_collector::ephemeris_collector(msystem &_ms):ms(_ms){
    t_start=ms.ephemeris_time();
    rebind();
}

int_t ephemeris_collector::decompose(int_t bid){
    if(bid<0){
        int_t bn=blist.size();
        int_t rootid=-1;
        int_t nroot=0;
        for(int_t i=0;i<bn;++i){
            if(blist[i].pid<0){
                rootid=i;
                ++nroot;
            }
        }
        return nroot==1?decompose(rootid):-1;
    }

    barycen &b=blist[bid];
    int_t nret=0;
    for(const auto cid:b.children)
        nret+=decompose(cid);

    if(b.gid>=0){
        nret+=decompose(b.gid);
        nret+=decompose(b.hid);
    }

    if(b.pid>=0){
        barycen &p=blist[b.pid];
        if(bid==p.hid){
            b.r=NAN;
            b.v=NAN;
        }
        else if(bid==p.gid){
            b.r=b.r_sys-blist[p.hid].r_sys;
            b.v=b.v_sys-blist[p.hid].v_sys;
        }
        else{
            b.r=b.r_sys-p.r;
            b.v=b.v_sys-p.v;
        }
    }
    else{
        b.r=b.r_sys;
        b.v=b.v_sys;
    }

    b.r_sys=NAN;
    b.v_sys=NAN;

    return nret+1;
}

int_t ephemeris_collector::compose(int_t bid){
    if(bid<0){
        int_t bn=blist.size();
        int_t rootid=-1;
        int_t nroot=0;
        for(int_t i=0;i<bn;++i){
            if(blist[i].pid<0){
                rootid=i;
                ++nroot;
            }
        }
        return nroot==1?compose(rootid):-1;
    }

    barycen &b=blist[bid];
    int_t nret=1;

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
            real gdm=g.GM_sys/p.GM;
            b.r_sys=p.r-gdm*g.r;
            b.v_sys=p.v-gdm*g.v;
        }
        else{
            b.r_sys=b.r+p.r;
            b.v_sys=b.v+p.v;
        }
    }

    mpvec cracc(0),cvacc(0);
    for(const auto cid:b.children){
        barycen &c=blist[cid];
        cracc+=c.r*c.GM_sys;
        cvacc+=c.v*c.GM_sys;
    }
    b.r=b.r_sys-cracc/b.GM_sys;
    b.v=b.v_sys-cvacc/b.GM_sys;

    if(b.gid>=0){
        nret+=compose(b.hid);
        nret+=compose(b.gid);
    }

    for(const auto cid:b.children)
        nret+=compose(cid);

    return nret;
}

void ephemeris_collector::rebind(){
    blist=ms.get_barycens();
    int_t bn=blist.size();
    if(orbital_data.empty()){
        int_t mn=ms.size();
        int_t t_eph=ms.ephemeris_time();
        orbital_data.resize(mn);
        for(orbital_datapack_t &d:orbital_data)
            d.t_start=t_eph;
        rotational_data.resize(mn);
        for(rotational_datapack_t &d:rotational_data)
            d.t_start=t_eph;
    }

    std::vector<std::vector<int_t>> barycen_mids(bn);
    for(const barycen &b:blist)if(b.hid<0){
        orbital_data[b.mid].tid=b.tid;
        int_t bp=b.mid;
        do{
            barycen_mids[bp].push_back(b.mid);
            bp=blist[bp].pid;
        } while(bp>=0);
    }

    for(orbital_datapack_t &d:orbital_data){
        int_t pid=blist[d.tid].pid;
        if(pid<0){
            d.parent_barycen_id=-1;
            continue;
        }
        std::set<int_t> pids,tids;
        barycen &p=blist[pid];
        if(p.hid<0)pids.insert(p.mid);
        else{
            auto &h=barycen_mids[p.hid];
            auto &g=barycen_mids[p.gid];
            pids.insert(h.begin(),h.end());
            pids.insert(g.begin(),g.end());
        }
        auto &b=barycen_mids[d.tid];
        tids.insert(b.begin(),b.end());
        d.parent_barycen_id=barycen_ids.insert({{pids,tids},(int_t)barycen_ids.size()}).first->second;
    }
}

void ephemeris_collector::record(){
    update_barycens();
    decompose();

    int_t mn=ms.size();
    int_t t_eph=ms.ephemeris_time();
    for(int_t i=0;i<mn;++i){
        orbital_datapack_t &od=orbital_data[i];
        const mass &m=ms[i];
        const barycen &b=blist[od.tid];
        od.t_end=t_eph;
        vec r=b.r,v=b.v,w=m.w,x=m.s.x,z=m.s.z;
        fwrite(&r,sizeof(vec),1,&od.data);
        fwrite(&v,sizeof(vec),1,&od.data);

        rotational_datapack_t &rd=rotational_data[i];
        rd.t_end=t_eph;
        fwrite(&w,sizeof(vec),1,&rd.data);
        fwrite(&x,sizeof(vec),1,&rd.data);
        fwrite(&z,sizeof(vec),1,&rd.data);
    }
}

void ephemeris_collector::extract(std::vector<MFILE> &ephm_files,bool force){
    int_t mn=ms.size();
    std::vector<int_t> ex_entry;
    const auto *pold_blist=&blist;
    std::vector<barycen> _old_blist_slot;
    if(!force){
        std::vector<int_t> polds;
        polds.reserve(mn);
        for(int_t i=0;i<mn;++i)
            polds.push_back(orbital_data[i].parent_barycen_id);
        _old_blist_slot.swap(blist);
        rebind();
        //not necessary: when(!force=analyse()) (blist=ms.blist) is updated by analyse()
        //but make sure
        update_barycens();
        decompose();
        pold_blist=&_old_blist_slot;
        for(int_t i=0;i<mn;++i)
            if(orbital_data[i].parent_barycen_id!=polds[i])
                ex_entry.push_back(i);
    }
    else{
        ex_entry.reserve(mn);
        for(int_t i=0;i<mn;++i)
            ex_entry.push_back(i);
    }
    int_t t_eph=ms.ephemeris_time();
    //avoid realloc && MFILE *invalidate
    ephm_files.reserve(ephm_files.size()+(1+force)*ex_entry.size());
    MFILE *findex=&ephm_files[0];
    for(auto i:ex_entry){
        orbital_datapack_t &d=orbital_data[i];
        const mass &m=ms[i];
        const barycen &b=blist[d.tid];

        index_entry_t idat;
        idat.fid=ephm_files.size()-Configs::ExportHeaderCount+1;
        idat.sid=m.sid;
        idat.t_start=d.t_start;
        idat.t_end=d.t_end;
        fwrite(&idat,sizeof(idat),1,findex);

        MFILE *mfp=&d.data;
        mfp->set_name(strprintf("%lld%s",idat.fid,Configs::SaveDataExtension));
        ephm_files.push_back(std::move(*mfp));
        if(force)
            continue;
        mfp->reset();
        d.t_start=t_eph;
        d.t_end=t_eph;
        vec r=b.r,v=b.v;
        fwrite(&r,sizeof(vec),1,mfp);
        fwrite(&v,sizeof(vec),1,mfp);
    }

    //record blist for [t_start,t_end]
    index_entry_t idat;
    idat.fid=0;
    idat.sid=blist.size();
    idat.t_start=t_start;
    idat.t_end=t_eph;
    fwrite(&idat,sizeof(idat),1,findex);

    msystem::save_barycen_structure(findex,*pold_blist);

    if(force)for(auto i:ex_entry){
        rotational_datapack_t &d=rotational_data[i];
        const mass &m=ms[i];
        vec w=m.w,x=m.s.x,z=m.s.z;

        index_entry_t idat;
        idat.fid=-(i+1);
        idat.sid=m.sid;
        idat.t_start=d.t_start;
        idat.t_end=d.t_end;
        fwrite(&idat,sizeof(idat),1,findex);

        MFILE *mfp=&d.data;
        mfp->set_name(strprintf("%lld%s",idat.fid,Configs::SaveDataExtension));
        ephm_files.push_back(std::move(*mfp));
    }

    t_start=t_eph;
}
