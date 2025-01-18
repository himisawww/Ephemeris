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

    LogAnnouncement("%s%s (cpu version)\n%s",
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


int main(int argc,const char **argv){
#if 0
    ephemeris_collector::convert_format("F:\\Temp\\ephm\\Ephemeris\\v0.3\\testold",3600);
    return 0;
#endif
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