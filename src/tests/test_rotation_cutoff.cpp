#include"tests/tests.h"
#include"utils/logger.h"
#include"utils/memio.h"

#if 1
#define TEST_CUTOFF_804_YEARS   1
#define TEST_CUTOFF_804_DT      600
#else   //longer, but more robust
#define TEST_CUTOFF_804_YEARS   200
#define TEST_CUTOFF_804_DT      300
#endif
#define TEST_CUTOFF_804_LON     0.0007
#define TEST_CUTOFF_804_LAT     0.000014

#define TEST_CUTOFF_HAUMEA_DT       600
#define TEST_CUTOFF_HAUMEA_STEP     100000
#define TEST_CUTOFF_HAUMEA_MISALIGN 1e-5

int test_rotation_cutoff(){
    fast_real max_misalign=0,max_lon=0,max_lat=0;
    {
        msystem ms=get_test_subsystem({"N02"});
        if(ms.size()!=1)
            return 1;
        const mass &Haumea=ms["N02"];
        for(int_t i=0;i<TEST_CUTOFF_HAUMEA_STEP;++i){
            ms.integrate(TEST_CUTOFF_HAUMEA_DT,1);
            fast_real misalign=std::atan2(fast_mpvec(Haumea.s.z*Haumea.GL).norm(),fast_real(Haumea.s.z%Haumea.GL));
            checked_maximize(max_misalign,misalign);
        }
        if(!(max_misalign<TEST_CUTOFF_HAUMEA_MISALIGN))
            return 2;
    }
    {
        msystem ms=get_test_subsystem({"899","804"});
        if(ms.size()!=2)
            return 3;
        const mass &Neptune=ms["899"];
        mass &Thalassa=ms["804"];
        Logger::ScopedSettings loglv(LogLevel::ERROR);
        for(int_t i=0;i<8766*TEST_CUTOFF_804_YEARS;++i){
            ms.integrate(TEST_CUTOFF_804_DT,3600/TEST_CUTOFF_804_DT);
            fast_mpvec r(Thalassa.s.tolocal(Neptune.r-Thalassa.r));
            r.normalize();
            fast_real ry=std::abs(r.y),rz=std::abs(r.z);
            checked_maximize(max_lon,ry);
            checked_maximize(max_lat,rz);
        }
        if(!(max_lon<TEST_CUTOFF_804_LON))
            return 4;
        if(!(max_lat<TEST_CUTOFF_804_LAT))
            return 5;
    }

    LogInfo("\n      Passed(%f, [%f, %f]), ",
        max_misalign/TEST_CUTOFF_HAUMEA_MISALIGN,
        max_lon/TEST_CUTOFF_804_LON,
        max_lat/TEST_CUTOFF_804_LAT);
    return 0;
}
