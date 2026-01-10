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
        
        if(argc==4&&strcmp(argv[1],"CONVERT_FORMAT")==0&&fix_dir>=0&&t_years>0)
            return ephemeris_collector::convert_format(argv[2],t_years);

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
        "   ip: full path to configuration file, optional;\n"
        "       default: use built-in solar system initial at J2000;\n\n"
        "   op: full path to output ephemerides and checkpoints;\n"
        "       when [op] contains previous checkpoints, [ip] takes no effect;\n\n"
        "  dir: direction of integration, can be +/-, optional;\n"
        "       +: forward, -: backward, default: both;\n"
        "    t: integrate for [t]-years;\n\n"
        "examples:\n\n"
        "   // first run, integrate default SolarSystem 20 years forward and backward (1980~2020):\n"
        "   exe_name  .\\results\\dat  20\n"
        "   // resume previous run, integrate 20 years backward (1960~1980):\n"
        "   exe_name  .\\results\\dat  -20\n\n"
        "   // load the default initial, do not integrate, save to output:\n"
        "   exe_name  .\\internal.output  0\n"
        "   // unzip the output, edit system_initial, then use it as custom initial:\n"
        "   exe_name  .\\system_initial\\Edited_Config.txt  .\\custom\\dat  20\n\n"
        "press Enter to exit, or input [t/T] to run tests:"
    );
    if(int i=getchar();i=='t'||i=='T')
        return test_all();
    return 0;
}

int main_for_ksp(const char *eph_path,const char *principia_config,const char *export_path);
int main(int argc,const char **argv){
    return main_for_ksp(
        "F:\\Temp\\ephm\\Ephemeris\\Ephemeris\\SolarSystem",
        "E:\\SteamLibrary\\steamapps\\KSP_RSS\\GameData\\Principia\\real_solar_system\\gravity_model.cfg",
        "R:\\"
        );
#if 0
    std::vector<const char*> subset{
        "10", "199", "299", "301", "399", "401", "402", "499", "501", "502",
        "503", "504", "505", "514", "515", "516", "599", "601", "602", "603",
        "604", "605", "606", "607", "608", "610", "611", "612", "613", "614",
        "615", "616", "617", "618", "632", "633", "634", "635", "649", "653",
        "699", "701", "702", "703", "704", "705", "706", "707", "708", "709",
        "710", "711", "712", "713", "714", "715", "725", "726", "727", "799",
        "801", "803", "804", "805", "806", "807", "808", "814", "899", "901",
        "999", "N01", "N01S1", "N02", "N02S1", "N02S2"
    };
    ephemeris_collector::convert_format("F:\\Temp\\ephm\\Ephemeris\\EphemerisCompressed\\SolarSystem",3600,&subset);
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