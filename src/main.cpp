#include"physics/mass.h"
#include<iostream>
#include<thread>
#include"modules/ephemeris_generator.h"
#include"utils/zipio.h"
#include"tests/tests.h"
#include"configs.h"
#include"utils/logger.h"

int de_worker(ephemeris_generator *egen,int dir){
    return egen->make_ephemeris(dir);
}

int main_fun(int argc,const char **argv){

    LogAnnouncement("%s%s\n%s",
        "Ephemeris Integrator ",Configs::VersionString,
        "Github: https://github.com/himisawww/Ephemeris \n\n");

    do{
        if(argc==2&&strcmp(argv[1],"RUN_TEST")==0)return test_all();

        if(argc<3||argc>4)break;

        const char *t_str=argv[argc-1];

        ephemeris_generator egen;
        auto &t_years=egen.t_years;
        auto &fix_dir=egen.fix_dir;

        if(*t_str=='+')fix_dir=1;
        else if(*t_str=='-')fix_dir=-1;
        else fix_dir=0;

        t_str+=std::abs(fix_dir);

        t_years=-1;
        if(1!=sscanf(t_str,"%lf",&t_years)||t_years<0)break;
        
        egen.ip=argc>3?argv[argc-3]:nullptr;
        egen.op=argv[argc-2];

        bool fwd=fix_dir>=0,bak=fix_dir<=0;
#ifdef NDEBUG
        std::thread th_future;
        std::thread th_past;

        if(fwd)th_future=std::thread(de_worker,&egen,1);
        if(bak)th_past=std::thread(de_worker,&egen,-1);

        if(fwd)th_future.join();
        if(bak)th_past.join();
#else
        if(fwd)de_worker(&egen,1);
        if(bak)de_worker(&egen,-1);
#endif
        return 0;
    } while(0);

    LogAnnouncement("%s",
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

/*/convert zips to old data pack
int convert_format(const char *path){
    std::string sop=path;
    for(int dir=1;dir>=-1;dir-=2){
        const char *fwdbak=dir>0?"fwd":"bak";

        std::string zckpt;
        size_t cur_index=1;
        MFILE *fout=mopen(sop+"."+fwdbak,MFILE_STATE::WRITE_FILE);
        do{
            zckpt=strprintf("%s.%llu.%s.zip",sop.c_str(),cur_index,fwdbak);
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
*/
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