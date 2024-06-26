#include"math/Keplerians.h"
#include"math/random.h"
#include"utils/calctime.h"
#include<cstdio>
#include<algorithm>

static double max_relative_error=0;
// max allowed relative error
#define TEST_EPSILON 1e-13
#define TEST_N 71

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
                vec vp=r.perpunit();
                vp+=rotation_matrix(r,randomreal()*(2*pi))%vp;
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
        fprintf(stderr,
            "\nMax Rel. Error %.16le Too Large",
            max_relative_error);
        return 1;
    }

    printf("Passed(%f), Done in %fs",max_relative_error/TEST_EPSILON,CalcTime()-start_time);
    return 0;
}
