#include<cstdio>
#include<algorithm>
#include"math/Keplerians.h"
#include"math/random.h"
#include"utils/calctime.h"
#include"utils/logger.h"

static double max_relative_error=0;
// max allowed relative error
#define TEST_EPSILON 1e-13
#define TEST_N 71

double test_rv_reproduce(const vec &r,const vec &v){
    vec nr,nv;
    ephem_orb k(0,r,v);
    k.rv(0,nr,nv);
    double rn=r.norm();
    double vref=checked_max(1/std::sqrt(rn),v.norm());
    double rerr=checked_max((nr-r).norm()/rn,(nv-v).norm()/vref);
    checked_maximize(max_relative_error,rerr);
    return rerr;
}

int test_kepler(){
    {
        double vn=1e50;
        do{
            double s2=std::sqrt(2);
            double cn=1e-1;
            const double cdiv=std::pow(1e-18/cn,double(1)/TEST_N);
            for(int_t i=0;i<TEST_N;++i){
                cn*=cdiv;
                vec r=randomdirection();
                vec v=randomdirection();
                vec dv=v*cn;
                //test random/up/down-ward throw & free falling
                test_rv_reproduce(r,v*vn);
                test_rv_reproduce(r,r*vn);
                test_rv_reproduce(r,-r*vn);
                test_rv_reproduce(r,(r+dv)*vn);
                test_rv_reproduce(r,-(r+dv)*vn);
                //test parabolic
                test_rv_reproduce(r,dv+s2*randomdirection());
                test_rv_reproduce(r,dv+s2*r);
                test_rv_reproduce(r,dv-s2*r);
                vec vp=r.asc_node();
                vp+=rotation_matrix(r,randomreal()*Constants::pi_mul2)%vp;
                vp+=dv;
                test_rv_reproduce(r,vp+r);
                //test appe
                test_rv_reproduce(r,vp*vn);
                //test circular
                test_rv_reproduce(r,vp);
            }
            vn/=2;
        } while(vn);
    }

    if(!(max_relative_error<TEST_EPSILON)){
        LogError(
            "\nMax Rel. Error %.16le Too Large",
            max_relative_error);
        return 1;
    }

    LogInfo("\n      Passed(%f), ",max_relative_error/TEST_EPSILON);
    return 0;
}
