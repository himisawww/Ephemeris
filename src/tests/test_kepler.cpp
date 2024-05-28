#include"math/Keplerians.h"
#include"math/random.h"
#include"utils/calctime.h"
#include<cstdio>
#include<algorithm>

static double max_relative_error=0;
// max allowed relative error
#define TEST_EPSILON 1e-13
#define TEST_N 17

double test_rv_reproduce(const vec &r,const vec &v){
    vec nr,nv;
    ephem_orb k(0,r,v);
    k.rv(0,nr,nv);
    double rn=r.norm();
    double vref=std::max(1/std::sqrt(rn),v.norm());
    double rerr=std::max((nr-r).norm()/rn,(nv-v).norm()/vref);
    max_relative_error=std::max(max_relative_error,rerr);
    return rerr;
}

int test_kepler(){
    double start_time=CalcTime();
    {
        double vn=1e50;
        do{
            fast_real s2=std::sqrt(2);
            for(int_t i=0;i<TEST_N;++i){
                //test random/up/down-ward throw & free falling
                vec r=randomdirection();
                test_rv_reproduce(r,randomdirection()*vn);
                test_rv_reproduce(r,r*vn);
                test_rv_reproduce(r,-r*vn);
                //test parabolic;
                vec v=s2*randomdirection();
                test_rv_reproduce(r,v+randomdirection()*vn);
                v=r+r.perpunit();
                v+=rotation_matrix(r,randomreal()*(2*pi))%v;
                test_rv_reproduce(r,v+randomdirection()*vn);
            }
            vn/=2;
        } while(vn);
    }

    if(!(max_relative_error<TEST_EPSILON)){
        fprintf(stderr,
            "\nMax Rel. Error %.16le Too Large",
            max_relative_error);
        return 1;
    }

    printf("Passed(%f), Done in %fs",max_relative_error/TEST_EPSILON,CalcTime()-start_time);
    return 0;
}
