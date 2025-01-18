#include"tests/tests.h"
#include<thread>
#include"utils/logger.h"
#include"utils/calctime.h"
#include"utils/threadpool.h"
#include"utils/memio.h"
#include"utils/crc32.h"

#define TEST_DELTA_T        300
#define TEST_COMBINED_T     28800
#define TEST_TOTAL_T        86400

typedef fast_real msystem_err_t[5];

static constexpr msystem_err_t
    ref_cgpu_err     ={ 1e-6,   1e-8,   1e-16,  1e-13,  1e-13   },
    ref_combined_err ={ 1e-4,   1e-8,   1e-11,  1e-07,  1e-07   },
    ref_cutoff_err   ={ 1e-3,   1e-7,   1e-12,  1e-07,  1e-07   };

static constexpr fast_real ref_maxrec_err=1e-15;

static fast_real msystem_compare(const msystem_err_t &ref_err,const msystem &ms1,const msystem &ms2,bool comp_maxrecs=false){
    msystem_err_t merr={0};
    const size_t msize=ms1.size();
    if(msize!=ms2.size()||ms1.ephemeris_time()!=ms2.ephemeris_time())
        return NAN;
    fast_real ret=0;
    for(size_t i=0;i<msize;++i){
        const mass &m1=ms1[i];
        const mass &m2=ms2[i];
        checked_maximize(merr[0],fast_mpvec(m1.r-m2.r).norm());
        checked_maximize(merr[1],fast_mpvec(m1.v-m2.v).norm());
        checked_maximize(merr[2],fast_mpvec(m1.w-m2.w).norm());
        checked_maximize(merr[3],fast_mpvec(m1.s.x-m2.s.x).norm());
        checked_maximize(merr[4],fast_mpvec(m1.s.z-m2.s.z).norm());
        if(!comp_maxrecs)continue;
        checked_maximize(ret,std::abs(m1.min_distance-m2.min_distance)/checked_min(m1.min_distance,m2.min_distance));
        checked_maximize(ret,std::abs(m1.max_influence-m2.max_influence)/checked_min(m1.max_influence,m2.max_influence));
    }
    ret/=ref_maxrec_err;
    const fast_real errfactor=TEST_TOTAL_T/ms1.ephemeris_time();
    for(int i=0;i<5;++i)
        checked_maximize(ret,merr[i]*errfactor/ref_err[i]);
    return ret;
}

static uint32_t msystem_crc32(const msystem &ms,uint32_t crc32val){
    MFILE fchpt;
    ms.save_checkpoint(&fchpt);
    return crc32(fchpt.data(),fchpt.size(),crc32val);
}

int test_integrator(){
    msystem mcombine=get_test_msystem();
    const size_t n_mass=mcombine.size();
    msystem mhalfdt=mcombine;
    msystem mcpu=mcombine;
    mcpu.clear_accel();
    mcpu.accel();
    msystem mgpu=mcombine;
    mgpu.clear_accel();
    mgpu.accel(1);

    double scpu=CalcTime();
    mcpu.integrate(TEST_DELTA_T,1,0);
    double sgpu=CalcTime();
    scpu=sgpu-scpu;
    struct thread_pool_guard{
        thread_pool_guard(){ ThreadPool::thread_local_pool_alloc(); }
        ~thread_pool_guard(){ ThreadPool::thread_local_pool_free(); }
    } _;
    sgpu=CalcTime();
    mgpu.integrate(TEST_DELTA_T,1,1);
    sgpu=CalcTime()-sgpu;

    //compare cpu and gpu
    fast_real max_err=msystem_compare(ref_cgpu_err,mcpu,mgpu,true);
    if(!(max_err<1)){
        LogError(
            "\nMax Reduced Difference between CPU & GPU Integrator %.16le Too Large",
            max_err);
        return 1;
    }

    double sgpu2=CalcTime();
    mgpu.integrate(TEST_DELTA_T,(TEST_COMBINED_T/TEST_DELTA_T)-1,1);
    double scombined=CalcTime();
    sgpu2=sgpu+scombined-sgpu2;
    mcombine.combined_integrate(TEST_DELTA_T,TEST_COMBINED_T/TEST_DELTA_T,1);
    scombined=CalcTime()-scombined;

    //compare gpu and combined
    checked_maximize(max_err,msystem_compare(ref_combined_err,mgpu,mcombine));
    if(!(max_err<1)){
        LogError(
            "\nMax Reduced Difference between C/GPU & Combined Integrator %.16le Too Large",
            max_err);
        return 2;
    }

    std::thread thhalf([&](){
        ThreadPool::thread_local_pool_alloc();
        mhalfdt.combined_integrate(TEST_DELTA_T/2,TEST_COMBINED_T/TEST_DELTA_T,2*TEST_TOTAL_T/TEST_COMBINED_T,0);
        ThreadPool::thread_local_pool_free();
        });

    mcombine.combined_integrate(TEST_DELTA_T,TEST_COMBINED_T/TEST_DELTA_T,TEST_TOTAL_T/TEST_COMBINED_T-1);
    thhalf.join();

    //compare combined and halfdt
    checked_maximize(max_err,msystem_compare(ref_cutoff_err,mcombine,mhalfdt));
    if(!(max_err<1)){
        LogError(
            "\nMax Cutoff Error of Combined Integrator %.16le Too Large",
            max_err);
        return 3;
    }

    uint32_t mcrc;
    mcrc=msystem_crc32(mcpu,0);
    mcrc=msystem_crc32(mgpu,mcrc);
    mcrc=msystem_crc32(mcombine,mcrc);
    mcrc=msystem_crc32(mhalfdt,mcrc);
    sgpu=scpu/sgpu;
    scombined=sgpu2/scombined;
    LogInfo("\n      Passed(%f)[%.2f * %.2f = %.2f] @ %08x, ",max_err,sgpu,scombined,sgpu*scombined,mcrc);
    return 0;
}
