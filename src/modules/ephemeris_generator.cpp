#include"ephemeris_generator.h"
#include"utils/zipio.h"
#include"configs.h"
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

    do{
        if(time_idx+iunit>isize)iunit=isize-time_idx;

        //prepare mem files to save
        int_t mn=ms.mlist.size();
        int_t max_dp_perfile=1+iunit;
        std::vector<MFILE> zms(ExportHeaderCount+mn);
        MFILE &mf_ckpt=zms[0];
        MFILE &mf_readme=zms[1];
        MFILE &mf_struct=zms[2];
        MFILE &mf_time=zms[3];
        mf_time.reserve(max_dp_perfile*sizeof(int_t));
        MFILE *mf_mlist=&zms[ExportHeaderCount];
        for(int_t mi=0;mi<mn;++mi){
            MFILE &mf=mf_mlist[mi];
            mf.reserve(max_dp_perfile*(5*sizeof(vec)));
        }

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

            ephc.update_barycens();
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

            fwrite(&mst_eph,sizeof(int_t),1,&mf_time);

            for(int_t mi=0;mi<mn;++mi){
                const mass &m=ms.mlist[mi];
                MFILE *mfp=&mf_mlist[mi];
                vec r=m.r,v=m.v,w=m.w,x=m.s.x,z=m.s.z;
                fwrite(&r,sizeof(vec),1,mfp);
                fwrite(&v,sizeof(vec),1,mfp);
                fwrite(&w,sizeof(vec),1,mfp);
                fwrite(&x,sizeof(vec),1,mfp);
                fwrite(&z,sizeof(vec),1,mfp);
            }
        }

        int_t t_eph_end=ms.ephemeris_time();

        //prepare checkpoint/readme/structure files
        ms.save_checkpoint(&mf_ckpt);

        fprintf(&mf_readme,
            "Calculated & Generated by Ephemeris Integrator %s\n"
            "   Github: https://github.com/himisawww/Ephemeris \n"
            "   Author: %s\n\n"
            "Number of Objects:  %lld\n"
            "   Time Range (s): [%lld, %lld]\n"
            "   Time Step  (s):  %lld\n\n"
            "      Time Format:   < t_eph(s): int64 >\n"
            "      Data Format:   < r(m), v(m/s), w(rad/s), x_axis(direction), z_axis(direction): vec3{double x,y,z;} >\n\n"
            "  Object List (index & sid):  \n"
            ,
            VersionString,AuthorName, mn, t_eph_start,t_eph_end, (t_eph_end-t_eph_start)/iunit
            );

        ms.print_structure(&mf_struct);
        
        mf_ckpt.set_name(SaveNameCheckpoint);
        mf_readme.set_name(SaveNameReadme);
        mf_struct.set_name(SaveNameStructure);
        mf_time.set_name(SaveNameTimestamps);
        for(int_t mi=0;mi<mn;++mi){
            MFILE &mf=mf_mlist[mi];
            std::string sid((char*)&ms.mlist[mi].sid);
            fprintf(&mf_readme,
                "%12lld : %s\n",mi,sid.c_str()
            );
            mf.set_name(sid+SaveDataExtension);
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
    return 0;
}

ephemeris_collector::ephemeris_collector(msystem &_ms):ms(_ms),blist(_ms.get_barycens()){

}



void ephemeris_collector::extract(std::vector<MFILE> &ephm_files,bool force){

}
