#include"physics/mass.h"
#include<iostream>
#include<thread>
#include"modules/ephemeris_generator.h"
#include"utils/zipio.h"
#include"tests/tests.h"
#include"configs.h"
#include"utils/logger.h"
#include"modules/ephemeris_reader.h"
#include"utils/calctime.h"

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

int main(int argc,const char **argv){
    struct check_version{
        bool pass;
        check_version():pass(true){
            const char *vstr=Configs::VersionString;
            size_t vsize=strlen(vstr);
            if(vsize==0||'0'>vstr[vsize-1]||vstr[vsize-1]>'9'){
                LogAnnouncement("Warning: This executable is compiled from development branch of code.\n");
                if(!file_exist("./_NOTES/DEVELOP")){
                    LogCritical("\n         User shall either find a release version, or compile an executable using main branch.\n\n");
                    pass=false;
                }
            }
        }
        ~check_version(){
            LogAnnouncement("The program is about to exit. Press Enter to continue...");
            getchar();
        }
    } chkv;
    if(!chkv.pass)return 0;

#if 0
    double s=CalcTime();
    ephemeris_reader ereader("f:\\temp\\ephm\\ephemeris\\Ephemeris\\SolarSystem");
    if(!ereader)
        return -1;
    int_t errs=ereader.make_cache();
    if(errs)printf("cache error: %lld\n",errs);
    printf("Loaded %llu objects in [%lld, %lld]\n",ereader.size(),ereader.t_min(),ereader.t_max());
    MFILE *fout=mopen("r:\\test.bin",MFILE_STATE::WRITE_FILE);
    if(!fout)
        return -3;
    const auto &earth=ereader["399"];
    const auto &moon=ereader["301"];
    earth.select(ereader.ORBIT);
    moon.select(ephemeris_reader::ORBIT);
    double rmin=INFINITY,rmax=-INFINITY,vavg=0,vcount=0;
    for(double t=0;t<Constants::year*40;t+=3600){
        if(!ereader.checkout(t))
            return -2;
        vec r(moon->r-earth->r),v(moon->v-earth->v);
        fwrite(&r,sizeof(vec),1,fout);
        fwrite(&v,sizeof(vec),1,fout);
        vcount+=1;
        vavg+=v.norm();
        double rn=r.norm();
        checked_minimize(rmin,rn);
        checked_maximize(rmax,rn);
    }
    printf("%fs\n",CalcTime()-s);
    printf("[%f, %f] km @ %f m/s\n",rmin/1000,rmax/1000,vavg/vcount);
    fclose(fout);
    /* not cached:
    Loaded 504 objects in [-65008656000, 65008656000]
    78.708161s,71.835861s
    [356445.428445, 406706.956274] km @ 1022.351862 m/s
    */
    return 0;


#if 0
    htl::vector<const char*> subset{
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

#else
#endif
    const char *m_argv[]={
        argv[0],
        //"F:\\Temp\\ephm\\Ephemeris\\SolarSystem\\SolarSystem_Config.txt",
        "f:\\Temp\\ephm\\Ephemeris\\Ephemeris\\SolarSystem",
        //"F:\\Temp\\ephm\\MoonsFit\\grad\\Test401",
        "60"
        //"RUN_TEST"
    };
    const int m_argc=sizeof(m_argv)/sizeof(char *);
    return main_fun(m_argc,m_argv);
    
    //convert_format("R:\\testcg\\result");
    //return 0;
    return 0;
}