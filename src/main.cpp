#include"physics/CelestialSystem.h"
#include<iostream>
#include<thread>
#include<mutex>
#include"utils/zipio.h"
#include"utils/wcs_convert.h"
#include"utils/calctime.h"
#include"tests/tests.h"
#include"configs.h"

namespace Configs{

const int ExportHeaderCount=4;
const char *SaveNameCheckpoint="checkpoint.dat";
const char *SaveNameReadme="readme.txt";
const char *SaveNameStructure="structure.json";
const char *SaveNameTimestamps="timestamps.dat";
const char *SaveDataExtension=".dat";
const char *SaveNameInitialDirectory="system_initial";
const char *VersionString="v0.2.0.beta";
const char AuthorName[]={104, 105, 109, 196, 171, 197, 155, 196, 129, 0};

}

std::mutex io_mutex;

const char *ip;
const char *op;

double t_years;
//  1: only do forward integration
// -1: only do backward integration
//  0: do both, default
int fix_dir=0;


void print_string(MFILE *mf,int_t level,const std::string &str){
    if(mf->tell()==0||mf->data()[mf->tell()-1]=='\n'){
        //indent
        std::string indent(level,' ');
        fwrite(indent.data(),indent.size(),1,mf);
    }
    fwrite(str.data(),str.size(),1,mf);
}
void print_structure(const msystem &ms,MFILE *mf,int_t root=-1,int_t level=0){
    if(root<0){
        root=ms.blist.size();
        do{
            if(--root<0)return;
        } while(!(ms.blist[root].pid<0));
    }

    const barycen& br=ms.blist[root];
    int_t cn=br.children.size();
    std::string barycen_name;

    if(br.hid<0){
        barycen_name=(const char *)&ms.mlist[br.mid].sid;
        if(cn==0){
            print_string(mf,level,"\""+barycen_name+"\"");
            return;
        }
        print_string(mf,level,"{\n");
        print_string(mf,level,"       \"ID\":\""+barycen_name+"\"");
    }
    else{
        barycen_name="Barycenter[";
        barycen_name+=(const char *)&ms.mlist[ms.blist[br.hid].mid].sid;
        if(ms.blist[br.hid].hid>=0)barycen_name+=" System";
        barycen_name+=", ";
        barycen_name+=(const char *)&ms.mlist[ms.blist[br.gid].mid].sid;
        if(ms.blist[br.gid].hid>=0)barycen_name+=" System";
        barycen_name+="]";
        print_string(mf,level,"{\n");
        print_string(mf,level,"       \"ID\":\""+barycen_name+"\",\n");
        print_string(mf,level,"  \"Primary\":");
        print_structure(ms,mf,br.hid,level+12);
        print_string(mf,level,",\n");
        print_string(mf,level,"\"Secondary\":");
        print_structure(ms,mf,br.gid,level+12);
    }

    if(cn){
        print_string(mf,level,",\n");
        print_string(mf,level," \"Children\":[\n");
        for(int_t i=0;i<cn;){
            print_structure(ms,mf,br.children[i],level+16);
            print_string(mf,level,(++i==cn)+",\n");
        }

        print_string(mf,level,"            ]\n");
    }
    else{
        print_string(mf,level,"\n");
    }
    print_string(mf,level,"}");
}

int de_worker(int dir){
    using namespace Configs;

    msystem ms;
    std::string sop=op;
    std::string ickpt=sop+".0.zip";
    const char *fwdbak=dir>0?"fwd":"bak";
    size_t cur_index=1;

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
        izippack zp(ickpt);
        for(const auto &zf:zp){
            if(zf.name()==SaveNameCheckpoint){
                MFILE mf;
                zf.dumpfile(mf);
                ms.load_checkpoint(&mf);
                break;
            }
        }
        printf("Loaded %lld bodies from checkpoint %s\n",ms.mlist.size(),ickpt.c_str());
    }
    else if(ip){
        ms.load(ip,ickpt.c_str());
    }

    if(!ms.mlist.size()){
        fprintf(stderr,
            "Failed to Load System.\n"
            "The program will exit.\n");
        exit(-1);
    }

    io_mutex.unlock();

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
        fprintf(stderr,
            "Non-integer Delta_t(%fs):\n    This may cause round-off errors on timestamps and is disallowed.\n",
            dt);
        return -3;
    }
    
    dt*=dir;

    int_t isize=t_years*Constants::year/ms.data_cadence;
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

        int_t t_eph_start=int_t(ms.t_eph.hi)+int_t(ms.t_eph.lo);

        //ephemeris integration
        for(int_t i=0;i<=iunit;++i){
            if(i>0){
                if(use_cpu     )ms.integrate         (dt,          jsize,0);
           else if(use_gpu     )ms.integrate         (dt,          jsize,1);
           else if(use_combined)ms.combined_integrate(dt,n_combine,jsize,1);
            }

            int_t mst_eph=int_t(ms.t_eph.hi)+int_t(ms.t_eph.lo);
            
            io_mutex.lock();
            if(fix_dir||dir>0){
                static double s=CalcTime();
                static double oldt=-INFINITY;
                static int_t skip_count=0,skip_size=1;
                if(++skip_count>=skip_size||i==iunit){
                    double yr=mst_eph/Constants::year;
                    double t=CalcTime();
                    printf(" Integrating. t_eph: %.6fyr, time: %.6fs\r",yr,t-s);
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

        int_t t_eph_end=int_t(ms.t_eph.hi)+int_t(ms.t_eph.lo);

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

        print_structure(ms,&mf_struct);
        
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
        zckpt.resize(sop.size()+30);
        zckpt.resize(
            sprintf(zckpt.data(),"%s.%llu.%s.zip",sop.c_str(),cur_index,fwdbak)
        );
        ++cur_index;

        io_mutex.lock();
        {
            ozippack zp(zckpt);
            zp.swap(zms);
            printf("\nSaving ephemeris & checkpoint %s\n",zckpt.c_str());
        }
        io_mutex.unlock();
        
        time_idx+=iunit;
    }while(time_idx<isize);
    return 0;
}

int main_fun(int argc,const char **argv){

    printf("%s%s\n%s",
        "Ephemeris Integrator ",Configs::VersionString,
        "Github: https://github.com/himisawww/Ephemeris \n\n");

    do{
        if(argc==2&&strcmp(argv[1],"RUN_TEST")==0)return test_fun();

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

        bool fwd=fix_dir>=0,bak=fix_dir<=0;

#ifdef NDEBUG
        std::thread th_future;
        std::thread th_past;

        if(fwd)th_future=std::thread(de_worker,1);
        if(bak)th_past=std::thread(de_worker,-1);

        if(fwd)th_future.join();
        if(bak)th_past.join();
#else
        if(fwd)de_worker(1);
        if(bak)de_worker(-1);
#endif
        return 0;
    } while(0);

    printf("%s",
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
        MFILE *fout=mopen(sop+"."+fwdbak,MFILE_STATE::WRITE_FILE);
        do{
            zckpt.resize(sop.size()+30);
            zckpt.resize(
                sprintf(zckpt.data(),"%s.%llu.%s.zip",sop.data(),cur_index,fwdbak)
            );
            if(!file_exist(zckpt))break;

            izippack zp(zckpt);
            MFILE mf_time;
            std::vector<MFILE> mf_mlist;
            int_t mi=-Configs::ExportHeaderCount;
            for(const auto &zf:zp){
                std::string zfn=zf.name();
                if(mi>=0){
                    mf_mlist.resize(mi+1);
                    zf.dumpfile(mf_mlist[mi]);
                }
                else if(zfn==Configs::SaveNameTimestamps){
                    zf.dumpfile(mf_time);
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

    const char *m_argv[]={
        argv[0],
        //"F:\\Temp\\ephm\\Ephemeris\\SolarSystem\\SolarSystem_Config.txt",
        //"f:\\Temp\\ephm\\Ephemeris\\TestNew\\test5",
        //"F:\\Temp\\ephm\\MoonsFit\\grad\\Test401",
        //"0.05"
        "RUN_TEST"
    };
    const int m_argc=sizeof(m_argv)/sizeof(char *);
    return main_fun(m_argc,m_argv);
    
    //convert_format("R:\\testcg\\result");
    //return 0;
    return 0;
}