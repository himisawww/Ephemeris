#include"ephemeris_generator.h"
#include"configs.h"
#include"utils/zipio.h"
#include"utils/calctime.h"
#include"utils/logger.h"
#include"utils/threadpool.h"
#include"modules/ephemeris_compressor.h"

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
    else
        ms.load_internal(ickpt.c_str());

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
    int_t jsize=(int_t)std::round(ms.data_cadence/dt);
    if(jsize<1)jsize=1;
    dt=ms.data_cadence/jsize;

    int_t n_combine;
    if(use_combined){
        n_combine=(int_t)std::round(dt/ms.delta_t);
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

    int_t isize=int_t(t_years*Constants::year/ms.data_cadence);
    int_t iunit=int_t(ms.max_ephm_length/ms.data_cadence);

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

    ThreadPool::thread_local_pool_alloc();
    static bool newline=true;
    auto check_newline=[](){
        if(!newline){
            LogInfo("\n");
            newline=true;
        }
        return true;
    };
    do{
        if(time_idx+iunit>isize)iunit=isize-time_idx;

        //prepare mem files to save
        int_t mn=ms.mlist.size();
        int_t max_dp_perfile=1+iunit;
        htl::vector<MFILE> zms(ExportHeaderCount);

        int_t t_eph_start=ms.ephemeris_time();

        bool resync=false;
        ephemeris_collector ephc(ms);
        //ephemeris integration
        for(int_t i=0;i<=iunit;++i){
            if(i>0){
                if(use_cpu     )ms.integrate         (dt,          jsize,0);
           else if(use_gpu     )ms.integrate         (dt,          jsize,1);
           else if(use_combined)ms.combined_integrate(dt,n_combine,jsize,1,&ephc);
                resync=!ephc.synchronized();
            }

            ephc.record();

            std::string exerr;
            if(i==iunit)
                exerr=ephc.extract(zms,true);
            else if(resync)
                exerr=ephc.extract(zms,false);

            int_t mst_eph=ms.ephemeris_time();
            
            io_mutex.lock();
            if(resync&&check_newline())LogInfo(
                "Info: At Ephemeris Time: %lld s\n"
                "   System orbital structure is updated.\n",
                mst_eph);
            if(exerr.size()&&check_newline())LogWarning("%s",exerr.c_str());
            if(fix_dir||dir>0){
                static double s=CalcTime();
                static double oldt=-INFINITY;
                static int_t skip_count=0,skip_size=1;
                if(++skip_count>=skip_size||i==iunit){
                    double yr=mst_eph/Constants::year;
                    double t=CalcTime();
                    LogInfo(" Integrating. t_eph: %.6fyr, time: %.6fs\r",yr,t-s);
                    newline=false;
                    int_t new_skip_size=(int_t)std::round(skip_size/(t-oldt));
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
        MFILE &mf_datdir=zms[3];

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
        
        mf_index.set_name(SaveNameIndex);
        mf_ckpt.set_name(SaveNameCheckpoint);
        mf_readme.set_name(SaveNameReadme);
        mf_datdir.set_name(SaveNameDirectory);

        for(int_t mi=0;mi<mn;++mi){
            fprintf(&mf_readme,
                "%12lld : %s\n",mi,ms[mi].get_ssid()
            );
        }

        //save .zip
        std::string zckpt;
        zckpt=strprintf("%s.%llu.%s.zip",sop.c_str(),cur_index,fwdbak);
        ++cur_index;

        //sort .zip
        std::vector<std::pair<uint64_t,int_t>> index_sort;
        for(int_t i=0,n=zms.size();i<n;++i){
            const auto &zname=zms[i].get_name();
            uint64_t k=!zname.find(SaveNameDirectory);
            if(k){
                auto fext='.'+get_file_extension(zname);
                if(!fext.find(SaveBarycentricOffsetDataExtension))
                    k=2;
                else if(!fext.find(SaveRotationalDataExtension)
                      ||!fext.find(SaveOrbitalDataExtension))
                    k=3;
            }
            index_sort.emplace_back(k,i);
        }
        std::sort(index_sort.begin(),index_sort.end());
        htl::vector<MFILE> szms;
        szms.reserve(zms.size());
        for(const auto &p:index_sort)
            szms.emplace_back(std::move(zms[p.second]));
        szms.swap(zms);

        io_mutex.lock();
        {
            check_newline();
            ephemeris_compressor::compress(zms);
            ozippack zp(zckpt);
            zp.swap_content(zms);
            LogInfo("Saving ephemeris & checkpoint %s\n",zckpt.c_str());
        }
        io_mutex.unlock();
        
        time_idx+=iunit;
    }while(time_idx<isize);
    ThreadPool::thread_local_pool_free();

    return 0;
}

ephemeris_collector::ephemeris_collector(msystem &_ms):ms(_ms){
    t_start=ms.ephemeris_time();
    t_link=NAN;
    rebind();
}

bool ephemeris_collector::synchronized(){
    return ms.analyse()==t_bind;
}

void ephemeris_collector::rebind(){
    t_bind=ms.analyse();
    blist=ms.get_barycens();
    int_t bn=blist.size();
    int_t t_eph=ms.ephemeris_time();
    if(data.empty()){
        data.resize(ms.size());
        for(datapack_t &d:data)
            d.t_start=t_eph;
    }

    htl::vector<htl::vector<int_t>> barycen_mids(bn);
    for(const barycen &b:blist)if(b.hid<0){
        data[b.mid].tid=b.tid;
        int_t bp=b.mid;
        do{
            barycen_mids[bp].push_back(b.mid);
            bp=blist[bp].pid;
        } while(bp>=0);
    }

    for(datapack_t &d:data){
        int_t pid=blist[d.tid].pid;
        if(pid<0){
            d.parent_barycen_id=-1;
            continue;
        }
        htl::set<int_t> pids,tids;
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
        d.parent_barycen_id=barycen_ids.try_emplace({pids,tids},barycen_ids.size()).first->second;
    }

    auto &bcache_map=bcache_maps.emplace_back();
    int_t ttag=t_eph+1;
    for(auto &c:cache)c.t_end=t_eph;
    for(int_t bi=0;bi<bn;++bi){
        const barycen &b=blist[bi];
        if(b.children.size()<=1)
            continue;
        htl::set<int_t> bids,chids;
        if(b.hid<0)bids.insert(b.mid);
        else{
            auto &h=barycen_mids[b.hid];
            auto &g=barycen_mids[b.gid];
            bids.insert(h.begin(),h.end());
            bids.insert(g.begin(),g.end());
        }
        for(int_t cid:b.children){
            auto &c=barycen_mids[cid];
            chids.insert(c.begin(),c.end());
        }
        auto ires=cache_ids.try_emplace({bids,chids},cache.size());
        int_t cache_id=ires.first->second;
        bcache_map.try_emplace(bi,cache_id);
        if(ires.second)
            cache.emplace_back().t_start=t_eph;
        cache[cache_id].t_end=ttag;
    }
    for(auto it=cache_ids.begin();it!=cache_ids.end();){
        const auto &c=cache[it->second];
        if(c.t_end==ttag)
            ++it;
        else
            it=cache_ids.erase(it);
    }
}

void ephemeris_collector::record(){
    blist.update_barycens(ms);
    blist.decompose();

    int_t mn=ms.size();
    int_t t_eph=ms.ephemeris_time();
    for(int_t i=0;i<mn;++i){
        datapack_t &d=data[i];
        const mass &m=ms[i];
        const barycen &b=blist[d.tid];
        d.t_end=t_eph;
        vec r=b.r,v=b.v,w=m.w,x=m.s.x,z=m.s.z;
        fwrite(&r,sizeof(vec),1,&d.orbital_data);
        fwrite(&v,sizeof(vec),1,&d.orbital_data);
        fwrite(&w,sizeof(vec),1,&d.rotational_data);
        fwrite(&x,sizeof(vec),1,&d.rotational_data);
        fwrite(&z,sizeof(vec),1,&d.rotational_data);
    }

    for(const auto &pbc:bcache_maps.back()){
        cachepack_t &c=cache[pbc.second];
        const barycen &b=blist[pbc.first];
        mpvec cr,cv;
        blist.get_children_offset(b,cr,cv);
        c.t_end=t_eph;
        vec r=cr,v=cv;
        fwrite(&r,sizeof(vec),1,&c.offset_data);
        fwrite(&v,sizeof(vec),1,&c.offset_data);
    }
}

void msystem::record_substeps(fast_real dt,bool initialize){
    if(tidal_childlist.empty()||!p_collector)return;
    auto &mc=*p_collector;
    if(mc.t_substep!=dt){
        p_collector=nullptr;
        return;
    }

    auto &subsys=mc.sublists.at(mlist[tidal_parent].sid);
    bsystem &blist=subsys.sublist;

    blist.update_barycens(*this);
    blist.decompose();
    
    for(auto &bp:subsys.subcache_map){
        mpvec cr,cv;
        const auto &b=blist[bp.first];
        blist.get_children_offset(b,cr,cv);
        MFILE *mcache=&mc.cache[bp.second].offset_subdata;
        if(!(initialize&&mcache->size())){
            vec r=cr,v=cv;
            fwrite(&r,sizeof(vec),1,mcache);
            fwrite(&v,sizeof(vec),1,mcache);
        }
    }

    int_t bn=blist.size();
    for(int_t i=0;i<bn;++i){
        const barycen &bi=blist[i];
        if(bi.hid>=0)
            continue;
        const mass &m=mlist[bi.mid];
        const barycen &b=blist[bi.tid];
        if(b.pid<0)
            continue;
        int_t sid=m.sid;
        auto &sd=initialize?mc.subdata[sid]:mc.subdata.at(sid);
        MFILE *morb=&sd.orbital_data;
        if(!(initialize&&morb->size())){
            vec r=b.r,v=b.v;
            fwrite(&r,sizeof(vec),1,morb);
            fwrite(&v,sizeof(vec),1,morb);
        }
        MFILE *mrot=&sd.rotational_data;
        if(!(initialize&&mrot->size())){
            vec w=m.w,x=m.s.x,z=m.s.z;
            fwrite(&w,sizeof(vec),1,mrot);
            fwrite(&x,sizeof(vec),1,mrot);
            fwrite(&z,sizeof(vec),1,mrot);
        }
    }
}

std::string ephemeris_entry::entry_name(bool rotational,bool substep) const{
    if(fid<=0)
        return std::string();
    return strprintf("%s.%lld%s%s",&sid,fid,
        rotational?Configs::SaveRotationalDataExtension:Configs::SaveOrbitalDataExtension,
        substep?Configs::SaveSubstepDataExtension:"");
}
std::string ephemeris_entry::offset_name(bool substep) const{
    return strprintf("%lld%s%s",fid,
        Configs::SaveBarycentricOffsetDataExtension,
        substep?Configs::SaveSubstepDataExtension:"");
}

std::string ephemeris_collector::extract(htl::vector<MFILE> &ephm_files,bool force){
    std::string err_msg;
    auto push_back_substep=[this,&ephm_files,&err_msg](ephemeris_entry &idata,MFILE &mf,size_t channel_size){
        size_t msize=mf.size();
        if(msize%channel_size==0&&(int_t(msize/channel_size)-1)*t_substep==idata.t_end-idata.t_start)
            ephm_files.push_back(std::move(mf));
        else
            err_msg+=strprintf("Unexpected size %llu when extracting <%s>.\n",msize,mf.get_name().c_str());
        mf.reset();
    };

    int_t mn=ms.size();
    htl::vector<int_t> ex_entry;
    const auto *pold_blist=&blist;
    bsystem _old_blist_slot;
    if(!force){
        htl::vector<int_t> polds;
        polds.reserve(mn);
        for(int_t i=0;i<mn;++i)
            polds.push_back(data[i].parent_barycen_id);
        _old_blist_slot.swap(blist);
        rebind();
        //not necessary: when(!force=analyse()) (blist=ms.blist) is updated by analyse()
        //but make sure
        blist.update_barycens(ms);
        blist.decompose();
        pold_blist=&_old_blist_slot;
        for(int_t i=0;i<mn;++i)
            if(data[i].parent_barycen_id!=polds[i])
                ex_entry.push_back(i);
    }
    else{
        ex_entry.reserve(mn);
        for(int_t i=0;i<mn;++i)
            ex_entry.push_back(i);
        //also extract offset_cache here
        MFILE mf_cindex;
        for(const auto &v:bcache_maps){
            size_t csize=v.size();
            fwrite(&csize,sizeof(csize),1,&mf_cindex);
            for(const auto &p:v){
                int_t bi=p.first,cfid=p.second+1;
                fwrite(&bi,sizeof(int_t),1,&mf_cindex);
                fwrite(&cfid,sizeof(int_t),1,&mf_cindex);
            }
        }

        int_t cache_index=0;
        for(cachepack_t &c:cache){
            ephemeris_entry icache;
            icache.fid=++cache_index;
            icache.sid=0;
            icache.t_start=c.t_start;
            icache.t_end=c.t_end;
            fwrite(&icache,sizeof(icache),1,&mf_cindex);

            c.offset_data.set_name(Configs::SaveNameDirectory+icache.offset_name(false));
            ephm_files.push_back(std::move(c.offset_data));

            if(c.offset_subdata.size()){
                c.offset_subdata.set_name(Configs::SaveNameDirectory+icache.offset_name(true));
                push_back_substep(icache,c.offset_subdata,2*sizeof(vec));
            }
        }
        mf_cindex.set_name(Configs::SaveNameBarycentricOffsetIndex);
        ephm_files.push_back(std::move(mf_cindex));
    }
    int_t t_eph=ms.ephemeris_time();
    for(auto i:ex_entry){
        datapack_t &d=data[i];
        const mass &m=ms[i];
        const barycen &b=blist[d.tid];

        ephemeris_entry idat;
        idat.fid=++file_index[m.sid];
        idat.sid=m.sid;
        idat.t_start=d.t_start;
        idat.t_end=d.t_end;
        fwrite(&idat,sizeof(idat),1,&ephm_files[0]);

        d.orbital_data.set_name(Configs::SaveNameDirectory+idat.entry_name(false,false));
        ephm_files.push_back(std::move(d.orbital_data));

        d.rotational_data.set_name(Configs::SaveNameDirectory+idat.entry_name(true,false));
        ephm_files.push_back(std::move(d.rotational_data));

        //save substep data, if exist
        auto it=subdata.find(m.sid);
        if(it!=subdata.end()){
            MFILE &morb=it->second.orbital_data;
            if(morb.size()){
                morb.set_name(Configs::SaveNameDirectory+idat.entry_name(false,true));
                push_back_substep(idat,morb,2*sizeof(vec));
            }
            MFILE &mrot=it->second.rotational_data;
            if(mrot.size()){
                mrot.set_name(Configs::SaveNameDirectory+idat.entry_name(true,true));
                push_back_substep(idat,mrot,3*sizeof(vec));
            }
        }

        if(force)
            continue;
        d.t_start=t_eph;
        d.t_end=t_eph;
        d.orbital_data.reset();
        d.rotational_data.reset();
        vec r=b.r,v=b.v,w=m.w,x=m.s.x,z=m.s.z;
        fwrite(&r,sizeof(vec),1,&d.orbital_data);
        fwrite(&v,sizeof(vec),1,&d.orbital_data);
        fwrite(&w,sizeof(vec),1,&d.rotational_data);
        fwrite(&x,sizeof(vec),1,&d.rotational_data);
        fwrite(&z,sizeof(vec),1,&d.rotational_data);
    }

    MFILE *findex=&ephm_files[0];
    //record blist for [t_start,t_end]
    ephemeris_entry idat;
    idat.fid=0;
    idat.sid=pold_blist->size();
    idat.t_start=t_start;
    idat.t_end=t_eph;
    fwrite(&idat,sizeof(idat),1,findex);
    pold_blist->save_barycen_structure(findex);

    MFILE &mf_struct=ephm_files.emplace_back();
    pold_blist->print_structure(&mf_struct,ms);
    mf_struct.set_name(strprintf("structure[%lld,%lld].json",idat.t_start,idat.t_end));

    t_start=t_eph;
    return err_msg;
}
