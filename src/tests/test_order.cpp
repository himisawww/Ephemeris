#include"tests/tests.h"
#include"modules/logger.h"

#define TEST_DELTA_T_CHOICE     1200,1080,960,880,800,720,640,600
#define TEST_DELTA_T_LCM        950400
#define TEST_PLANET             "499"
#define TEST_MOON               "401"

static double order_solve(double dratio,double dt0,double dt1,double dt2){
    const double b0=dt2/dt1;
    const double b1=dt0/dt1;
    //solve dratio=(1-b0^order)/(b1^order-1) for order
    const double logb0=std::log(b0);
    const double logb1=std::log(b1);
    double order=std::log(dratio)/(logb0-logb1)*2;
    if(order<3)order=3;
    double best_order=NAN;
    double last_diff=INFINITY;
    int n_iter=0;
    do{
        double b1order=std::pow(b1,order);
        double b0order=std::pow(b0,order);
        double diff=dratio*(b1order-1)-(1-b0order);
        if(!(std::abs(diff)<last_diff))break;
        best_order=order;
        last_diff=std::abs(diff);
        order-=diff/(b0order*logb0+b1order*dratio*logb1);
        ++n_iter;
    } while(1);
    if(n_iter<2)//means Newton method is diverging
        return NAN;
    return best_order;
}

int test_order(){
    std::vector<const char *> test_ssids={
        "10",TEST_PLANET,TEST_MOON
    };
    msystem mssrc=get_test_subsystem(test_ssids);
    if(mssrc.size()!=test_ssids.size())return 1;
    
    const auto moon_id=mssrc.get_mid(TEST_MOON);
    //mssrc.scale(moon_id,10);
    //mssrc.scale_geopotential(moon_id,0.1);
    //mssrc.accel();

    const int_t test_dt[]={TEST_DELTA_T_CHOICE};
    const int_t test_size=sizeof(test_dt)/sizeof(test_dt[0]);

    for(int_t i=0;i<test_size;++i){
        if(test_dt[i]<=0||TEST_DELTA_T_LCM%test_dt[i]||i&&test_dt[i]>=test_dt[i-1])
            return 2;
    }

    typedef vec moon_state[5];

    int_t order_weight=0;
    double order_translational=0,order_rotational=0;
    {
        ScopedLogLevelSettings loglvset(global_logger,LogLevel::ERROR);

        moon_state states[test_size];
        for(int_t i=0;i<test_size;++i){
            msystem ms=mssrc;
            ms.integrate(test_dt[i],TEST_DELTA_T_LCM/test_dt[i]);
            const mass &planet=ms[TEST_PLANET],&moon=ms[TEST_MOON];
            states[i][0]=vec(moon.r-planet.r);
            states[i][1]=vec(moon.v-planet.v);
            states[i][2]=vec(moon.w);
            states[i][3]=vec(moon.s.x);
            states[i][4]=vec(moon.s.z);
        }

        for(int_t i=1;i+1<test_size;++i){
            const moon_state &s0=states[i-1],&s1=states[i],&s2=states[i+1];
            double orders[5];
            int_t iweight=i*i;
            order_weight+=iweight;
            for(int_t j=0;j<5;++j){
                double dratio=(s1[j]-s2[j]).norm()/(s0[j]-s1[j]).norm();
                if(dratio<0||dratio>1)
                    return 3;
                orders[j]=order_solve(dratio,test_dt[i-1],test_dt[i],test_dt[i+1]);
                if(!(orders[j]>0))
                    return 4;
                if(!(orders[j]>(j<2?11.375:2.5))){
                    LogError(
                        "\nIntegrator Order @%c[%d] %f < %d Too Low",
                        j["rvwxz"],i,orders[j],j<2?12:3);
                    return 5;
                }
                (j<2?order_translational:order_rotational)+=orders[j]*iweight;
            }
        }
    }

    order_translational/=2*order_weight;
    order_rotational/=3*order_weight;

    LogInfo("\n      Passed([%f, %f]/[12, 3]), ",order_translational,order_rotational);
    return 0;
}
