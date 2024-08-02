#include"tests/tests.h"
#include"physics/ring.h"
#include"physics/geopotential.h"
#include"utils/logger.h"

#define TEST_EPSILON_CONSERVATION       1e-12
#define TEST_EPSILON_RING_CORRECTION    1e-8
#define TEST_INTEGRATION_STEP           6000

int test_conservation(){
    std::vector<const char *> test_ssids={
        "699","601","603","605","606"
    };
    msystem ms=get_test_subsystem(test_ssids);
    if(ms.size()!=test_ssids.size())return 1;

    //push adjust
    int_t saturn_id=ms.get_mid("699");
    ms.scale_geopotential(saturn_id,3);
    ms.scale_ring(saturn_id,1e7);
    ms.scale(saturn_id,1/(1+ms[saturn_id].ringmodel->GM_ratio));
    ms.accel();

    union momentum{
    public:
        struct{
            fast_mpvec translational_noring;
            fast_mpvec angular_noring;
            fast_mpvec translational;
            fast_mpvec angular;
        };
        fast_real data[12];
        momentum(){}
        momentum(fast_real p){
            for(int_t i=0;i<12;++i)
                data[i]=p;
        }
        void maximize(const momentum &other){
            for(int_t i=0;i<12;++i)
                checked_maximize(data[i],other.data[i]);
        }
        void minimize(const momentum &other){
            for(int_t i=0;i<12;++i)
                checked_minimize(data[i],other.data[i]);
        }
    };

    momentum max_momentum(-INFINITY),min_momentum(INFINITY);

    for(int_t i=0;i<TEST_INTEGRATION_STEP;++i){
        ms.integrate(900,1);

        momentum p(0);

        for(const auto &m:ms){
            fast_mpvec dm=m.GM*fast_mpvec(m.v);
            p.translational_noring+=dm;
            p.translational+=dm;

            fast_mpvec mGL(m.GL);
            fast_mpvec da=m.GM*fast_mpvec(m.r*m.v);
            fast_mpvec das=mGL*m.GM;
            p.angular_noring+=da+das;
            p.angular+=da+das;

            if(m.ringmodel){
                const ring &mr=*m.ringmodel;
                p.translational+=dm*(mr.GM_ratio/(1-mr.GM_ratio));
                p.angular+=da*(mr.GM_ratio/(1-mr.GM_ratio));
                p.angular+=(mGL/mGL.norm())*(mr.GL_R2*m.GM*m.R2);

                //changing direction of disk will also produce an angular momentum
                //however this cannot be modeled ...
#if 0
                fast_mpvec angular_accel;
                angular_accel=m.dtorque;
                //direction changing rate of disk pole, perpendicular to disk pole
                angular_accel-=angular_accel%mGL/(mGL%mGL)*mGL;
                angular_accel/=mGL.norm();
                //disk inertia along this axis
                fast_real din=(2*mr.A-mr.J2)/3*m.GM*m.R2;
                p.angular+=angular_accel*din;
#endif
            }
        }

        max_momentum.maximize(p);
        min_momentum.minimize(p);
    }

    momentum difference;
    for(int_t i=0;i<12;++i)
        difference.data[i]=max_momentum.data[i]-min_momentum.data[i];
    fast_real max_conserv_err=checked_max(
        difference.translational.norm()/max_momentum.translational.norm(),
        difference.angular.norm()/max_momentum.angular.norm());
    if(!(max_conserv_err<TEST_EPSILON_CONSERVATION)){
        LogError(
            "\nMax Relative Difference of Momentum %.16le Too Large",
            max_conserv_err);
        return 2;
    }

    fast_real max_conserv_err_ring_factor=max_conserv_err/checked_max(
        difference.translational_noring.norm()/max_momentum.translational.norm(),
        difference.angular_noring.norm()/max_momentum.angular.norm());
    if(!(max_conserv_err_ring_factor<TEST_EPSILON_RING_CORRECTION)){
        LogError(
            "\nMax Ring Correction Factor %.16le Too Large",
            max_conserv_err_ring_factor);
        return 3;
    }

    LogInfo("\n      Passed(%f), ",max_conserv_err/TEST_EPSILON_CONSERVATION);
    return 0;
}
