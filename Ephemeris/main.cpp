#include"ephemeris.h"
#include<iostream>
#include"utils.h"
#include<thread>

//  0: use CPU RungeKutta
//  1: use GPU RungeKutta
//  2: use CPU/GPU Combined RungeKutta
int USE_METHOD;


//  1 year = 8766 h
int_t t_years;
// dt = 10s * dt_ratio
// should not greater than 8640
fast_real dt_ratio;
//  print interval in hours
int_t print_interval;
//  1: only do forward integration
// -1: only do backward integration
//  0: do both, default
int fix_dir=0;
//fast_real solar_mass_loss_rate=-3.771963849337985E-1*(1+0.33812787124240545*0.01);
//fast_real solar_luminosity=3.828E+26;
//fast_real lunar_recpt=1.457;
//fast_real earth_recpt=1.224;
//fast_real earth_J2rate=-1.55772E-18;

const char *ip;
const char *op;

const char *extra_params="E:\\Temp\\ephm\\MoonsFit\\SolarSystem_ExtraParameters.txt";
const char *gp_dir="E:\\Temp\\ephm\\MoonsFit\\Geopotentials";
const char *ring_dir="E:\\Temp\\ephm\\MoonsFit\\Rings";


bool USE_MULTITHREAD;
//128M
#define MAX_RAMBUFFER_SIZE (1ull<<27)
//64k
#define REDIRECT_BUFFER_SIZE (1ull<<16)
void redirect(FILE *fdst,FILE *fsrc,size_t target){
    fseek(fsrc,0,SEEK_SET);
    size_t bufsize=REDIRECT_BUFFER_SIZE;
    char binbuff[REDIRECT_BUFFER_SIZE];
    size_t readsize;
    while(target){
        if(bufsize>target)bufsize=target;
        fread(binbuff,1,bufsize,fsrc);
        fwrite(binbuff,1,bufsize,fdst);
        target-=bufsize;
    }
    fseek(fsrc,0,SEEK_SET);
}

int de_worker(int dir){
    msystem ms;
    char fbuf[257];
    sprintf(fbuf,"%s%s",ip,dir==1?".fwd":".bak");
    bool from_ckpt=ms.load_checkpoint(fbuf);
    if(!from_ckpt&&!ms.load(ip,extra_params,gp_dir,ring_dir))exit(-1);
    int_t jsize=round(print_interval*360/dt_ratio);
    if(jsize<1)jsize=1;
    fast_real dt=fast_real(print_interval)*(3600*dir)/jsize;
    sprintf(fbuf,"%s%s",op,dir==1?".fwd":".bak");
    FILE *fout=fopen(fbuf,from_ckpt?"ab":"wb");

    if(!fout){
        fprintf(stderr,"Cannot open output file.");
        exit(-1);
    }
    bool useRambuffer=false;
    const char *pfbase=op+strlen(op);
    while(pfbase>op&&!(pfbase[-1]=='\\'||pfbase[-1]=='/'))
        --pfbase;
    sprintf(fbuf,"del R:\\%s%s",pfbase,dir==1?".fwd":".bak");
    FILE *frout=fopen(fbuf+4,"wb+");
    if(frout){
        useRambuffer=true;
        FILE *ftemp=fout;
        fout=frout;
        frout=ftemp;
    }

    int_t isize=t_years*8766/print_interval;
    for(int_t i=(dir==1&&!from_ckpt?0:1);i<=isize;++i){
        if(i>0){
            if(USE_METHOD==0)ms.integrate(dt,jsize,0);
       else if(USE_METHOD==1)ms.integrate(dt,jsize,1);
       else if(USE_METHOD==2){
                const fast_real mdt=300.;
                int_t nc=std::round(dt/mdt);
                ms.combined_integrate(dt/nc,nc,jsize,1);
            }
        }

            double mst_eph=(double)ms.t_eph;
            if(!USE_MULTITHREAD||dir==1){
                static double s=CalcTime();
                double yr=mst_eph/(8766*3600);
                double t=CalcTime()-s;
                printf(" %.6lfyr %.6lfs\r",yr,t);
            }
            fwrite(&mst_eph,sizeof(double),1,fout);

            for(const auto &m:ms.mlist){
                vec r=m.r,v=m.v,w=m.w,x=m.s.x,z=m.s.z;
                fwrite(&r,sizeof(vec),1,fout);
                fwrite(&v,sizeof(vec),1,fout);
                fwrite(&w,sizeof(vec),1,fout);
                fwrite(&x,sizeof(vec),1,fout);
                fwrite(&z,sizeof(vec),1,fout);
            }

        if(i%24==0){
            fflush(fout);
            if(useRambuffer){
                size_t frsize=_ftelli64(fout);
                if(frsize>MAX_RAMBUFFER_SIZE){
                    redirect(frout,fout,frsize);
                }
            }
        }
    }
    
    if(useRambuffer){
        size_t frsize=_ftelli64(fout);
        redirect(frout,fout,frsize);
        fclose(fout);
        fout=frout;
        system(fbuf);
    }

    fclose(fout);

    sprintf(fbuf,"%s%s",op,dir==1?".ckpt.fwd":".ckpt.bak");
    ms.save_checkpoint(fbuf);
    return 0;
}

int main(int argc,char **argv){
#if 0
    ip="E:\\Temp\\ephm\\MoonsFit\\grad\\Test2_CGC_300s.txt";
#ifdef NDEBUG
    //Combined Performance Test
    op="E:\\Temp\\ephm\\MoonsFit\\grad\\SolarSystem_MoonsBase_CombinedPerformance_2";
#else
    op="E:\\Temp\\ephm\\MoonsFit\\grad\\Debug";
#endif
    USE_METHOD=0;
    t_years=2;
    dt_ratio=30;
    print_interval=24;

#else  


    if(argc<3||argc>8){
        return -1;
    }

    ip=argv[1];
    op=argv[2];
    if(argc>3)USE_METHOD=atoi(argv[3]);
    if(argc>4)t_years=atoi(argv[4]);
    if(argc>5)dt_ratio=atof(argv[5]);
    if(argc>6)print_interval=atoi(argv[6]);
    if(argc>7)fix_dir=atoi(argv[7]);
#endif

#ifdef NDEBUG
    USE_MULTITHREAD=!fix_dir&&(USE_METHOD==0);
#else
    USE_MULTITHREAD=0;
#endif
    double s=CalcTime();
    if(USE_MULTITHREAD){
        std::thread th_future(de_worker,1);
        std::thread th_past(de_worker,-1);

        th_future.join();
        th_past.join();
    }
    else{
        if(fix_dir>=0)de_worker(1);
        if(fix_dir<=0)de_worker(-1);
    }
    printf("\n%lfs\n",CalcTime()-s);
    return 0;
}

double convert(int i){
    if(i==0)return 0;
    return (i>0?1:-1)*exp(abs(i)*0.00858291+17.7275);
}
int main_testring(){
    /*
    msystem ms;
    if(!ms.load(
        "F:\\Temp\\ephm\\MoonsFit\\RingTest\\SolarSystem_MoonsBase.txt",
        "F:\\Temp\\ephm\\MoonsFit\\RingTest\\SolarSystem_ExtraParameters.txt",
        "F:\\Temp\\ephm\\MoonsFit\\RingTest\\Geopotentials",
        "F:\\Temp\\ephm\\MoonsFit\\RingTest\\Rings")
        )exit(-1);
    */

    //test conservation
    /*
    FILE *fout=fopen("e:\\temp\\conserv_ring.bin","wb");

    for(int_t i=0;i<100000;++i){
        ms.integrate(300,1);

        fast_mpvec moment=0,momentring=0,amoment=0,amomentring=0;

        for(const auto &m:ms.mlist){
            fast_mpvec dm=m.GM*fast_mpvec(m.v);
            moment+=dm;
            momentring+=dm;

            fast_mpvec mGL(m.GL);
            fast_mpvec da=m.GM*fast_mpvec(m.r*m.v);
            fast_mpvec das=mGL*m.GM;
            amoment+=da+das;
            amomentring+=da+das;

            if(m.ringmodel){
                const ring &mr=*m.ringmodel;
                momentring+=dm*(mr.GM_ratio/(1-mr.GM_ratio));
                amomentring+=da*(mr.GM_ratio/(1-mr.GM_ratio));
                amomentring+=(mGL/mGL.norm())*(mr.GL*m.GM);
            
                //changing direction of disk will also produce an angular momentum
                //however this cannot be modeled ...
              //fast_mpvec angular_accel;
              //angular_accel=m.dtorque;
              ////direction changing rate of disk pole, perpendicular to disk pole
              //angular_accel-=angular_accel%mGL/(mGL%mGL)*mGL;
              //angular_accel/=mGL.norm();
              ////disk inertia along this axis
              //fast_real din=(2*mr.A-mr.J2)/3*m.GM*m.R2;
              //amomentring+=angular_accel*din;
            }
        }

        fwrite(&moment,sizeof(fast_mpvec),1,fout);
        fwrite(&amoment,sizeof(fast_mpvec),1,fout);
        fwrite(&momentring,sizeof(fast_mpvec),1,fout);
        fwrite(&amomentring,sizeof(fast_mpvec),1,fout);
        if(i%10000==0)printf(" %d\r",i);
    }
    printf("\n");

    fclose(fout);*/
    
    //test scaling
    /*
    mass &s=ms.mlist[ms.get_mid("699")];
    for(int i=0;i<10;++i){
        fast_mpvec ma=s.daccel;
        ring &mr=*s.ringmodel;
        printf("%.16le %.16le %.16le ",ma.x,ma.y,ma.z);

        double gm=0;
        for(int_t k=0;k<mr.N;++k){
            gm+=mr.c_table[k].Gs*(pi*mr.c_table[k].R*mr.c_table[k].R);
        }
        printf("%.16le\n",gm/(mr.GM_ratio/(1-mr.GM_ratio)));

        s.scale(0.1);
        ms.accel();
    }
    */

    //test force
    /*
    msystem ms;
    if(!ms.load(
        "F:\\Temp\\ephm\\MoonsFit\\SolarSystem_MoonsBase.txt",
        "F:\\Temp\\ephm\\MoonsFit\\SolarSystem_ExtraParameters.txt",
        "F:\\Temp\\ephm\\MoonsFit\\Geopotentials",
        "F:\\Temp\\ephm\\MoonsFit\\Rings")
        )exit(-1);

    ring &mr=*ms.mlist[ms.get_mid("699")].ringmodel;
    FILE *fout=fopen("e:\\temp\\ringtest.bin","wb");
    for(int i=-128;i<=128;++i)
    for(int j=-128;j<=128;++j)
    for(int k=-128;k<=128;++k){
        fast_mpvec r=randomdirection();
        r*=convert(i);
        fast_mpvec f=mr.sum(r);
        fwrite(&r,sizeof(fast_mpvec),1,fout);
        fwrite(&f,sizeof(fast_mpvec),1,fout);
    }
    fclose(fout);
    */
    return 0;
}

int main_test70678(){//test 70678
    msystem ms;
    if(!ms.load("F:\\Temp\\ephm\\MoonsFit\\70678\\70678.txt"))exit(-1);

    mass &u=ms.mlist[ms.get_mid("799")];
    mass &u6=ms.mlist[ms.get_mid("706")];
    mass &u7=ms.mlist[ms.get_mid("707")];
    mass &u8=ms.mlist[ms.get_mid("708")];

    double maxrn=0,minrn=INFINITY;
    double maxang=0;
    double maxangz=0;
    for(int_t i=0;i<100000;++i){
        ms.integrate(300,1);

        mass &us=u6;

        vec r=us.r-u.r,v=us.v-u.v,j=r*v;
        vec x=us.s.x,z=us.s.z;
        double rn=r.norm();
        if(rn<minrn)minrn=rn;
        if(rn>maxrn)maxrn=rn;
        double e=(maxrn-minrn)/(maxrn+minrn);
        r/=-rn;
        j/=j.norm();
        double ang=atan2((r*x).norm(),r%x)/degree;
        double angz=atan2((z*j).norm(),z%j)/degree;
        if(ang>maxang)maxang=ang;
        if(angz>maxangz)maxangz=angz;
        printf("%lfd, %lf, %lf, %lf\n",
            (double)ms.t_eph/86400.,
            maxang,maxangz,
            e);
    }
    return 0;
}