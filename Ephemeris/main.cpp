#include"ephemeris.h"
#include<iostream>
#include"utils.h"
#include<thread>
#include"ziptree.h"
#include"memfile.h"
#include<mutex>

const char *readme="readme.txt";
const char *checkpoint="checkpoint.dat";
const char *timestamps="timestamps.dat";
const char *version="v0.1.0";
const char author[]={104, 105, 109, 196, 171, 197, 155, 196, 129, 0};

std::mutex io_mutex,calc_mutex;

const char *ip;
const char *op;

//  1 year = 8766 h
double t_years;
//  1: only do forward integration
// -1: only do backward integration
//  0: do both, default
int fix_dir=0;

bool file_exist(const std::string &path){
    FILE *fin=fopen(path.c_str(),"rb");
    if(!fin)return false;

    fclose(fin);
    return true;
}

bool sanity(double dt){
    return dt>0&&dt!=INFINITY;
}

void dumpfile(mem_file &mf,const zipfile &zf){
    mf.wdata.resize(zf.filesize);
    zf.dumpfile(mf.wdata.data());
    mf.idata=mf.wdata.data();
    mf.isize=zf.filesize;
    mf.offset=0;
}

int de_worker(int dir){
    msystem ms;
    std::string sop=op;
    std::string ickpt=sop+".0.zip";
    const char *fwdbak=dir>0?"fwd":"bak";
    size_t cur_index=1;
    bool is_main_thread=(dir>0||fix_dir<0);

    //loading process
    io_mutex.lock();

    if(file_exist(ickpt)){
        // load file from ckpt
        std::string zckpt;

        do{
            zckpt.resize(sop.size()+30);
            zckpt.resize(
                sprintf(zckpt.data(),"%s.%llu.%s.zip",sop.data(),cur_index,fwdbak)
            );
            if(!file_exist(zckpt))break;
            ickpt.swap(zckpt);
            ++cur_index;
        } while(1);
        //ickpt exists
        zippack zp(strtowcs(ickpt));
        for(const auto &zf:zp){
            if(wcstostr(zf.filename)==checkpoint){
                mem_file mf;
                dumpfile(mf,zf);
                ms.load_checkpoint(&mf);
                break;
            }
        }
    }
    else if(ip){
        std::string sip=ip;
        size_t spos=1+sip.find_last_of("/\\");
        ms.load_dir(sip.substr(0,spos).c_str(),sip.substr(spos).c_str());
        if(ms.mlist.size()){
            zippack zp(strtowcs(ickpt),true);
            zp.zipmems.resize(1);
            zipmem &zm=zp.zipmems[0];
            mem_file mf;
            ms.save_checkpoint(&mf);
            zm.data.swap(mf.wdata);
            zm.filename=strtowcs(checkpoint);
        }
    }

    io_mutex.unlock();
    
    if(!ms.mlist.size()){
        fprintf(stderr,"Failed to Load System.\n");
        return -1;
    }

    //check config
    if(!sanity(ms.delta_t)||!sanity(ms.data_cadence)||!sanity(ms.max_ephm_length)||!sanity(ms.combined_delta_t)){
        fprintf(stderr,"Delta_t/Cadence/Max_Ephemeris_Length/Combined_Delta_t_Max should be finity positive real numbers.\n");
        return -2;
    }

    bool use_cpu     =ms.integrator==int_t(msystem::     CPU_RK12);
    bool use_gpu     =ms.integrator==int_t(msystem::     GPU_RK12);
    bool use_combined=ms.integrator==int_t(msystem::COMBINED_RK12);
    bool no_parallel =fix_dir||!use_cpu;
    
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
        fprintf(stderr,
            "Non-integer Delta_t(%fs):\n    This may cause round-off errors on timestamps and is disallowed.\n",
            dt);
        return -3;
    }
    
    dt*=dir;

    const double year_seconds=8766*3600;

    int_t isize=t_years*year_seconds/ms.data_cadence;
    int_t iunit=ms.max_ephm_length/ms.data_cadence;

    const int_t min_iunit=4096;
    if(isize<1){
        fprintf(stderr,"Nothing to do.\n");
        return 0;
    }
    if(iunit<min_iunit){
        fprintf(stderr,
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
        mem_file mf_time;
        mf_time.wdata.reserve(max_dp_perfile*sizeof(int_t));
        std::vector<mem_file> mf_mlist(mn);
        for(auto &mf:mf_mlist){
            mf.wdata.reserve(max_dp_perfile*(5*sizeof(vec)));
        }

        int_t t_eph_start=int_t(ms.t_eph.hi)+int_t(ms.t_eph.lo);

        //ephemeris integration
        if(no_parallel)calc_mutex.lock();
        for(int_t i=0;i<=iunit;++i){
            if(i>0){
                if(use_cpu     )ms.integrate         (dt,          jsize,0);
           else if(use_gpu     )ms.integrate         (dt,          jsize,1);
           else if(use_combined)ms.combined_integrate(dt,n_combine,jsize,1);
            }

            int_t mst_eph=int_t(ms.t_eph.hi)+int_t(ms.t_eph.lo);
            if(no_parallel||dir==1){
                static double s=CalcTime();
                double yr=mst_eph/year_seconds;
                double t=CalcTime()-s;
                printf(" Integrating. t_eph: %.6fyr, time: %.6fs\r",yr,t);
            }
            fwrite(&mst_eph,sizeof(int_t),1,&mf_time);

            for(int_t mi=0;mi<mn;++mi){
                const mass &m=ms.mlist[mi];
                mem_file *mfp=&mf_mlist[mi];
                vec r=m.r,v=m.v,w=m.w,x=m.s.x,z=m.s.z;
                fwrite(&r,sizeof(vec),1,mfp);
                fwrite(&v,sizeof(vec),1,mfp);
                fwrite(&w,sizeof(vec),1,mfp);
                fwrite(&x,sizeof(vec),1,mfp);
                fwrite(&z,sizeof(vec),1,mfp);
            }
        }
        if(no_parallel)calc_mutex.unlock();

        int_t t_eph_end=int_t(ms.t_eph.hi)+int_t(ms.t_eph.lo);

        //prepare checkpoint/readme files
        mem_file mf_ckpt,mf_readme;
        ms.save_checkpoint(&mf_ckpt);

        char buffer[512];
        int buffer_chars=sprintf(buffer,
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
            version,author, mn, t_eph_start,t_eph_end, (t_eph_end-t_eph_start)/iunit
            );
        fwrite(buffer,buffer_chars,1,&mf_readme);

        std::vector<zipmem> zms(3+mn);
        for(int_t mi=0;mi<mn;++mi){
            std::string sid((char*)&ms.mlist[mi].sid);
            buffer_chars=sprintf(buffer,
                "%12lld : %s\n",mi,sid.c_str()
            );
            fwrite(buffer,buffer_chars,1,&mf_readme);
            zms[3+mi].filename=strtowcs(sid+".dat");
            zms[3+mi].data.swap(mf_mlist[mi].wdata);
        }

        zms[0].filename=strtowcs(checkpoint);
        zms[0].data.swap(mf_ckpt.wdata);
        zms[1].filename=strtowcs(readme);
        zms[1].data.swap(mf_readme.wdata);
        zms[2].filename=strtowcs(timestamps);
        zms[2].data.swap(mf_time.wdata);

        //save .zip
        std::string zckpt;
        zckpt.resize(sop.size()+30);
        zckpt.resize(
            sprintf(zckpt.data(),"%s.%llu.%s.zip",sop.data(),cur_index,fwdbak)
        );
        ++cur_index;

        io_mutex.lock();
        {
            zippack zp(strtowcs(zckpt),true);
            zp.zipmems.swap(zms);
        }
        io_mutex.unlock();
        
        time_idx+=iunit;
    }while(time_idx<isize);
    return 0;
}

int main_fun(int argc,const char **argv){
    do{
        if(argc<3||argc>4)break;

        const char *t_str=argv[argc-1];

        if(*t_str=='+')fix_dir=1;
        else if(*t_str=='-')fix_dir=-1;
        else fix_dir=0;

        t_str+=std::abs(fix_dir);

        t_years=-1;
        if(1!=sscanf(t_str,"%lf",&t_years)||t_years<0)break;
        
        ip=argc>3?argv[argc-3]:nullptr;
        op=argv[argc-2];

#ifdef NDEBUG
        std::thread th_future;
        std::thread th_past;

        if(fix_dir>=0)th_future=std::thread(de_worker,1);
        if(fix_dir<=0)th_past=std::thread(de_worker,-1);

        if(fix_dir>=0)th_future.join();
        if(fix_dir<=0)th_past.join();
#else
        if(fix_dir>=0)de_worker(1);
        if(fix_dir<=0)de_worker(-1);
#endif
        return 0;
    } while(0);

    printf("%s%s\n%s",
        "Ephemeris Integrator ",version,
        "Github: https://github.com/himisawww/Ephemeris \n\n"
        "command line usage:\n\n"
        "   exe_name [[ip]] [op] [[[dir]]t] \n\n"
        "   ip: full path to configuration file (initial values & parameters\n"
        "           should be under same directory as configuration file);\n\n"
        "   op: full path to output ephemerides and checkpoints;\n"
        "       when [op] contains checkpoints of previous run,\n"
        "       [ip] can be omitted (and will be ignored if presence);\n\n"
        "  dir: direction of integration, can be +/-, optional;\n"
        "       +: forward, -: backward, default: both;\n\n"
        "    t: integrate the system for [t]-years;\n\n"
        "examples:\n\n"
        "   // first run, integrate SolarSystem 20 years forward and backward (1980~2020):\n"
        "   exe_name  .\\SolarSystem\\SolarSystem_Config.txt  .\\results\\dat  20\n\n"
        "   // resume previous run, integrate 20 years backward (1960~1980):\n"
        "   exe_name  .\\results\\dat  -20\n\n"
        "press Enter to exit..."
    );
    getchar();
    return 0;
}

//convert zips to old data pack
int convert_format(const char *path){
    std::string sop=path;
    for(int dir=1;dir>=-1;dir-=2){
        const char *fwdbak=dir>0?"fwd":"bak";

        std::string zckpt;
        size_t cur_index=1;
        FILE *fout=fopen((sop+"."+fwdbak).c_str(),"wb");
        do{
            zckpt.resize(sop.size()+30);
            zckpt.resize(
                sprintf(zckpt.data(),"%s.%llu.%s.zip",sop.data(),cur_index,fwdbak)
            );
            if(!file_exist(zckpt))break;

            zippack zp(strtowcs(zckpt));
            mem_file mf_time;
            std::vector<mem_file> mf_mlist;
            int_t mi=-3;
            for(const auto &zf:zp){
                std::string zfn=wcstostr(zf.filename);
                if(mi>=0){
                    mf_mlist.resize(mi+1);
                    dumpfile(mf_mlist[mi],zf);
                }
                else if(zfn==timestamps){
                    dumpfile(mf_time,zf);
                }
                ++mi;
            }

            vec v5[5];
            int_t it_eph;
            int_t n_data=0;
            while(1==fread(&it_eph,sizeof(int_t),1,&mf_time)){
                double mst_eph=it_eph;
                bool output=n_data||dir==1&&cur_index==1;
                if(output)fwrite(&mst_eph,sizeof(double),1,fout);
                for(auto &mf:mf_mlist){
                    fread(&v5,sizeof(vec),5,&mf);
                    if(output)fwrite(&v5,sizeof(vec),5,fout);
                }
                ++n_data;
            }

            ++cur_index;
        } while(1);
        fclose(fout);
    }
    return 0;
}

int main(int argc,const char **argv){
    return main_fun(argc,argv);
    
    /*
    const char *m_argv[]={
        argv[0],
        "E:\\Prog\\Ephemeris\\SolarSystem\\SolarSystem_Config.txt",
        "F:\\Temp\\ephm\\Ephemeris\\combine\\result",
        "0.01"
    };
    const int m_argc=sizeof(m_argv)/sizeof(char *);
    return main_fun(m_argc,m_argv);
    */
    
    //convert_format("R:\\testcg\\result");
    //return 0;
}