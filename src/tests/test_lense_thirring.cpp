#include"tests/tests.h"
#include"utils/logger.h"

#define TEST_LTCON_MA       3
#define TEST_LTCON_MB       1
#define TEST_LTCON_WA       vec( 0.00, 0.05, 0.10)
#define TEST_LTCON_WB       vec( 0.05,-0.10, 0.20)
#define TEST_LTCON_P        100
#define TEST_LTCON_E        0.1
#define TEST_LTCON_ORBITS   16
#define TEST_LTCON_EPSILON_MOMENTUM             4e-15
#define TEST_LTCON_EPSILON_ANGULAR_MOMENTUM     4e-17
#define TEST_LTCON_TRANSFER_ANGULAR_MOMENTUM    1e-5

#define TEST_MERCURY_ORBITS       10
#define TEST_MERCURY_GR_TARGET    42.9799
#define TEST_MERCURY_GR_EPSILON    0.0020
#define TEST_MERCURY_LT_TARGET   (-0.0020)
#define TEST_MERCURY_LT_EPSILON    0.0001

#define TEST_GPB_HEIGHT       656255.0
#define TEST_GPB_VELOCITY_FACTOR   1.000223
#define TEST_GPB_W              1000.0
#define TEST_GPB_ORBITS           60
#define TEST_GPB_PERIOD_TARGET  5862.6
#define TEST_GPB_PERIOD_EPSILON    1.0
#define TEST_GPB_LT_LS_TARGET   6600.0
#define TEST_GPB_LT_LS_EPSILON    15.0
#define TEST_GPB_LT_SS_TARGET     40.0
#define TEST_GPB_LT_SS_EPSILON     5.0
#define TEST_GPB_MISALIGN          1e-10

static int test_lt_conservation(){
    msystem ms;
    ms.push_back("A",mass().initialize(TEST_LTCON_MA,1.0));
    ms.push_back("B",mass().initialize(TEST_LTCON_MB,1.0));
    mass &mA=ms["A"];
    mass &mB=ms["B"];
    mA.w=Constants::c*TEST_LTCON_WA;
    mB.w=Constants::c*TEST_LTCON_WB;
    fast_real mu=mA.GM+mB.GM;
    mA.r.x=-TEST_LTCON_P*mB.GM/mu;
    mB.r.x=TEST_LTCON_P*mA.GM/mu;
    fast_real vmu=std::sqrt(mu/TEST_LTCON_P*(1+TEST_LTCON_E));
    mA.v.y=-vmu*mB.GM/mu;
    mB.v.y=vmu*mA.GM/mu;
    ms.accel();

    struct vec_range{
        fast_mpvec min,max;
        fast_real ref;
        vec_range():min(INFINITY),max(-INFINITY),ref(0){}

        void update(const fast_mpvec &l,fast_real _ref){
            checked_minimize(min.x,l.x);
            checked_minimize(min.y,l.y);
            checked_minimize(min.z,l.z);
            checked_maximize(max.x,l.x);
            checked_maximize(max.y,l.y);
            checked_maximize(max.z,l.z);
            checked_maximize(ref,_ref);
        }
    };
    vec_range lrange,lsrange,prange;
    mpvec s0=mA.GM*mA.GL+mB.GM*mB.GL;
    fast_real T=2*mu/100-vmu*vmu;
    T=Constants::pi_mul2*mu/(T*sqrt(T));

    for(int_t i=0;i<8*TEST_LTCON_ORBITS;++i){
        ms.integrate(T/128,16,0);

        mpvec lA=mA.GM*(mA.r*mA.v),sA=mA.GM*mA.GL;
        mpvec lB=mB.GM*(mB.r*mB.v),sB=mB.GM*mB.GL;

        lrange.update(lA+lB,fast_mpvec(lA).norm());
        lsrange.update(lA+lB+sA+sB-s0,fast_mpvec(sA).norm());
        if((i+1)%8==0){
            mpvec pA=mA.GM*mA.v,pB=mB.GM*mB.v;
            prange.update(pA+pB,fast_mpvec(pA).norm());
        }
    }

    fast_real lvar=(lrange.max-lrange.min).norm()/lrange.ref;
    fast_real lsvar=(lsrange.max-lsrange.min).norm()/lsrange.ref;
    fast_real pvar=(prange.max.norm()+prange.min.norm())/prange.ref;

    if(!(pvar<TEST_LTCON_EPSILON_MOMENTUM))
        return 6;
    if(!(lsvar<TEST_LTCON_EPSILON_ANGULAR_MOMENTUM))
        return 7;
    if(!(lvar>TEST_LTCON_TRANSFER_ANGULAR_MOMENTUM))
        return 8;
    return 0;
}

static int test_GPB(fast_real &wg,fast_real &wlt){
    msystem ms;
    ms.push_back("Earth",mass().initialize(3.986004355070227E14,6371000,0,3.306557720310232E-1,mat().from_harmonics(1.085060817270971E-3,0,0,0,0)));
    ms.push_back("GPB",mass().initialize(Constants::G,0.1));
    mass &Earth=ms["Earth"];
    mass &GPB=ms["GPB"];
    Earth.w.z=Constants::pi_mul2/(23.93447209*3600);
    GPB.r.x=Earth.R+TEST_GPB_HEIGHT;
    GPB.v.z=TEST_GPB_VELOCITY_FACTOR*std::sqrt(Earth.GM/fast_real(GPB.r.x));
    GPB.w.x=TEST_GPB_W;
    ms.accel();

    fast_real lastta=0,misalign=0;
    for(int_t i=0;i<8*TEST_GPB_ORBITS;++i){
        ms.integrate(TEST_GPB_PERIOD_TARGET/128,16,0);
        fast_real ta=std::atan2((fast_real)(GPB.r.z-Earth.r.z),(fast_real)(GPB.r.x-Earth.r.x));
        ta=lastta+angle_reduce(ta-lastta);
        lastta=ta;
        fast_mpvec sx(GPB.s.x);
        checked_maximize(misalign,std::atan2((sx*GPB.w).norm(),sx%GPB.w));
    }
    fast_real period=Constants::pi_mul2*ms.ephemeris_time()/lastta;
    fast_mpvec l=GPB.GL;
    // milli-arcsec/year
    fast_real wf=Constants::year*3600000/(Constants::degree*ms.ephemeris_time());
    wg =wf*std::asin(l.z/l.x);
    wlt=wf*std::asin(l.y/l.x);

    if(!(std::abs(period-TEST_GPB_PERIOD_TARGET)<TEST_GPB_PERIOD_EPSILON))
        return 1;
    if(!(std::abs(    wg-TEST_GPB_LT_LS_TARGET )<TEST_GPB_LT_LS_EPSILON ))
        return 2;
    if(!(std::abs(   wlt-TEST_GPB_LT_SS_TARGET )<TEST_GPB_LT_SS_EPSILON ))
        return 3;
    if(!(misalign<TEST_GPB_MISALIGN))
        return 9;
    return 0;
}

static int test_mercury(fast_real &wmp,fast_real &wmplt){
    const msystem &mssrc=get_test_msystem();
    const mass &ssrc=mssrc["10"];
    const mass &msrc=mssrc["199"];
    const fast_mpvec jm=(msrc.r-ssrc.r)*(msrc.v-ssrc.v);
    const fast_mpmat jmat(jm.asc_node(),0,jm.unit());
    for(int k=0;k<2;++k){
        msystem ms;
        ms.push_back("10",mass().initialize(ssrc.GM,ssrc.R,0,ssrc.inertia));
        ms.push_back("199",mass().initialize(msrc.GM,msrc.R,0,msrc.inertia));
        mass &Sun=ms["10"];
        mass &Mercury=ms["199"];
        Sun.r=ssrc.r;
        Sun.v=ssrc.v;
        Sun.w=k==0?fast_mpvec(0):ssrc.w;
        Sun.s=ssrc.s;
        Mercury.r=msrc.r;
        Mercury.v=msrc.v;
        Mercury.w=msrc.w;
        Mercury.s=msrc.s;
        ms.accel();

        auto llr=[&](){
            fast_mpvec r=Mercury.r-Sun.r,v=Mercury.v-Sun.v;
            return jmat.lon(v*(r*v)/(Mercury.GM+Sun.GM)-r.unit());
        };
        fast_real llr0=llr();
        ms.integrate(36193,210*TEST_MERCURY_ORBITS,0);//mercury orbit period = 36193*210 s
        fast_real llr1=llr();
        // arcsec/century
        (k==0?wmp:wmplt)=360000*Constants::year*(llr1-llr0)/(Constants::degree*ms.ephemeris_time());
    }
    wmplt-=wmp;
    if(!(std::abs(wmp  -TEST_MERCURY_GR_TARGET)<TEST_MERCURY_GR_EPSILON))
        return 4;
    if(!(std::abs(wmplt-TEST_MERCURY_LT_TARGET)<TEST_MERCURY_LT_EPSILON))
        return 5;
    return 0;
}

int test_lense_thirring(){
    int res_conserv=test_lt_conservation();
    if(res_conserv)return res_conserv;

    fast_real wg,wlt;
    int res_GPB=test_GPB(wg,wlt);
    if(res_GPB)return res_GPB;

    fast_real wmp,wmplt;
    int res_mercury=test_mercury(wmp,wmplt);
    if(res_mercury)return res_mercury;

    LogInfo("\n      Passed([%.2f,%.2f], [%.4f,%.4f]]), ",wg,wlt,wmp,wmplt);
    return 0;
}
