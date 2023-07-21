#define INLINE __device__
#include"ephemeris.h"
#include<stdlib.h>
#include<algorithm>
#include<cuda_runtime.h>
#include<cooperative_groups.h>

#define WARP_SIZE 32
#define CUDA_CORES 1920
#define MAXBLOCKS (1920/32)

#define cuda_max(a,b) ((a)>(b)?(a):(b))

//geopotential data, mlist[first].gpmodel==second
typedef std::vector<std::pair<int_t,geopotential*>> gpdata_t;
//ring data, mlist[first].ringmodel==second
typedef std::vector<std::pair<int_t,ring*>> ringdata_t;

#if 0
#define mycudaMalloc cudaMalloc
#define mycudaFree cudaFree
#define mymalloc malloc
#define myfree myfree
#else
void mycudaMalloc(void *devPtr,size_t size){
    static void *mem=0;
    static size_t memsize=0;
    if(memsize<size){
        cudaFree(mem);
        cudaMalloc(&mem,size);
        memsize=size;
    }
    *(void**)devPtr=mem;
}
void *mymalloc(size_t size){
    static void *mem=0;
    static size_t memsize=0;
    if(memsize<size){
        free(mem);
        mem=malloc(size);
        memsize=size;
    }
    return mem;
}
void myfree(void *){

}
void mycudaFree(void *){

}
#endif

struct cuda_rungekutta_kernel_config{
    //number of mass
    int nmass;
    int nblocks,mass_per_block;
    int nthreads,mass_per_thread;
    int_t n_step;
    fast_real dt;
    mass *dmlist;
    mass_state *x0,*f;
    real t_eph;

    void load(std::vector<mass> &mlist,gpdata_t &mgp,ringdata_t &mrg,fast_real _dt,int_t _nstep){
        dt=_dt;
        n_step=_nstep;
        int mn=mlist.size();
        nmass=mn;
        mass_per_block=(mn+MAXBLOCKS-1)/MAXBLOCKS;
        nblocks=(mn+mass_per_block-1)/mass_per_block;
        nthreads=WARP_SIZE*std::min(MAXBLOCKS/nblocks,(mn+WARP_SIZE-1)/WARP_SIZE);
        //Make sure nthreads is power of 2
        int new_nth;
        while(new_nth=nthreads&nthreads-1)nthreads=new_nth;
        mass_per_thread=(mn+nthreads-1)/nthreads;
        int_t grsize=0;
        for(int_t i=0;i<mn;++i){
            mass &m=mlist[i];
            if(m.gpmodel){
                int_t thissize=m.gpmodel->size();
                mgp.push_back({i,m.gpmodel});
                m.gpmodel=(geopotential*)grsize;
                grsize+=thissize;
            }
            if(m.ringmodel){
                int_t thissize=m.ringmodel->size();
                mrg.push_back({i,m.ringmodel});
                m.ringmodel=(ring*)grsize;
                grsize+=thissize;
            }
        }
        mycudaMalloc(&x0,nmass*(sizeof(mass)+26*sizeof(mass_state))+grsize);
        f=x0+nmass;
        dmlist=(mass*)(x0+26*nmass);
        void *grdata=mymalloc(grsize);
        for(auto &mgpi:mgp){
            int_t gpoffset=(int_t)mlist[mgpi.first].gpmodel;
            mlist[mgpi.first].gpmodel=(geopotential*)(gpoffset+(int_t)(dmlist+nmass));
            memcpy((geopotential*)(gpoffset+(int_t)grdata),mgpi.second,mgpi.second->size());
        }
        for(auto &mrgi:mrg){
            int_t rgoffset=(int_t)mlist[mrgi.first].ringmodel;
            mlist[mrgi.first].ringmodel=(ring*)(rgoffset+(int_t)(dmlist+nmass));
            memcpy((ring*)(rgoffset+(int_t)grdata),mrgi.second,mrgi.second->size());
        }
        cudaMemcpy(dmlist,mlist.data(),nmass*sizeof(mass),cudaMemcpyHostToDevice);
        cudaMemcpy(dmlist+nmass,grdata,grsize,cudaMemcpyHostToDevice);
        myfree(grdata);
    }

    void save(std::vector<mass> &mlist,gpdata_t &mgp,ringdata_t &mrg){
        cudaMemcpy(mlist.data(),dmlist,nmass*sizeof(mass),cudaMemcpyDeviceToHost);
        for(auto &mgpi:mgp){
            mlist[mgpi.first].gpmodel=mgpi.second;
        }
        for(auto &mrgi:mrg){
            mlist[mrgi.first].ringmodel=mrgi.second;
        }
        mycudaFree((mass_state*)dmlist-26*nmass);
    }
};

__constant__ cuda_rungekutta_kernel_config dkf;

struct maccel_1{
    fast_mpmat C_potential;
    fast_mpvec naccel;
    fast_real phi;
};
struct maccel_2{
    fast_mpvec gaccel,daccel,dtorque;
    fast_real min_distance;
    fast_real max_influence;
};

//copied from geopotential.cpp, subroutines are added ``static __device__''
#if 1
//static const int_t Precompute_Table_size[1+Max_N]={0, 0, 10, 25, 46, 74, 110, 155, 210};
static __device__ int_t Precompute_Table_size(int_t n){
    return (n-1)*(36+n*(10+n))/6;
}

#if 0

template<const int_t n>
fast_mpvec sum_n(const fast_mpvec &r,const fast_mpvec *cn){
    int_t np=n+1,cp=0;
    fast_real xbase(1);
    fast_mpvec an(0);
    for(int_t px=0;px<=np;++px){
        fast_real ybase(xbase);
        for(int_t py=0;py<=np-px;++py){
            fast_real zbase(ybase);
            for(int_t pz=0;pz<np-px-py;++pz)zbase*=r.z;
            an+=zbase*cn[cp];
            ++cp;

            ybase*=r.y;
        }
        xbase*=r.x;
    }
    return an;
}

static const uint8_t j_table[]={1,2,1,2,4,0,4,1,2,4,1,2,1,2,1,4,0,4,0,1,2,1,4,0,1,1,2,1,2,1,2,4,0,4,0,4,1,2,1,2,4,0,4,1,2,4,1,2,1,2,1,2,1,4,0,4,0,4,0,1,2,1,2,1,4,0,4,0,1,2,1,4,0,1,1,2,1,2,1,2,1,2,4,0,4,0,4,0,4,1,2,1,2,1,2,4,0,4,0,4,1,2,1,2,4,0,4,1,2,4,1,2,1,2,1,2,1,2,1,4,0,4,0,4,0,4,0,1,2,1,2,1,2,1,4,0,4,0,4,0,1,2,1,2,1,4,0,4,0,1,2,1,4,0,1,1,2,1,2,1,2,1,2,1,2,4,0,4,0,4,0,4,0,4,1,2,1,2,1,2,1,2,4,0,4,0,4,0,4,1,2,1,2,1,2,4,0,4,0,4,1,2,1,2,4,0,4,1,2,4};
//summation over only terms include Jn coefficients
template<const int_t n>
fast_mpvec sum_Jn(const fast_mpvec &r,const fast_mpvec *cn){
    int_t np=n+1,cp=0;
    fast_real xbase(1);
    fast_mpvec an(0);
    const uint8_t *jt=j_table+Precompute_Table_size(n-1);
    for(int_t px=0;px<=np;++px){
        fast_real ybase(xbase);
        for(int_t py=0;py<=np-px;++py){
            fast_real zbase(ybase);
            for(int_t pz=0;pz<np-px-py;++pz)zbase*=r.z;
            if(jt[cp]&4)an.x+=zbase*cn[cp].x;
            if(jt[cp]&2)an.y+=zbase*cn[cp].y;
            if(jt[cp]&1)an.z+=zbase*cn[cp].z;
            ++cp;

            ybase*=r.y;
        }
        xbase*=r.x;
    }
    return an;
}

static fast_mpvec (*const sum_Jn_funlist[])(const fast_mpvec &,const fast_mpvec *)={
    nullptr,nullptr,
    sum_Jn<int_t(2)>,
    sum_Jn<int_t(3)>,
    sum_Jn<int_t(4)>,
    sum_Jn<int_t(5)>,
    sum_Jn<int_t(6)>,
    sum_Jn<int_t(7)>,
    sum_Jn<int_t(8)>
};

static fast_mpvec (*const sum_n_funlist[])(const fast_mpvec &,const fast_mpvec *)={
    nullptr,nullptr,
    sum_n<int_t(2)>,
    sum_n<int_t(3)>,
    sum_n<int_t(4)>,
    sum_n<int_t(5)>,
    sum_n<int_t(6)>,
    sum_n<int_t(7)>,
    sum_n<int_t(8)>
};

#else
static __device__ fast_mpvec sum_2(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);an+=base2.z*cn[0];an+=base0.y*base1.z*cn[1];an+=base1.y*base0.z*cn[2];an+=base2.y*cn[3];an+=base0.x*base1.z*cn[4];an+=base0.x*base0.y*base0.z*cn[5];an+=base0.x*base1.y*cn[6];an+=base1.x*base0.z*cn[7];an+=base1.x*base0.y*cn[8];an+=base2.x*cn[9];return an;
}

static __device__ fast_mpvec sum_3(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);an+=base3.z*cn[0];an+=base0.y*base2.z*cn[1];an+=base1.y*base1.z*cn[2];an+=base2.y*base0.z*cn[3];an+=base3.y*cn[4];an+=base0.x*base2.z*cn[5];an+=base0.x*base0.y*base1.z*cn[6];an+=base0.x*base1.y*base0.z*cn[7];an+=base0.x*base2.y*cn[8];an+=base1.x*base1.z*cn[9];an+=base1.x*base0.y*base0.z*cn[10];an+=base1.x*base1.y*cn[11];an+=base2.x*base0.z*cn[12];an+=base2.x*base0.y*cn[13];an+=base3.x*cn[14];return an;
}

static __device__ fast_mpvec sum_4(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);an+=base4.z*cn[0];an+=base0.y*base3.z*cn[1];an+=base1.y*base2.z*cn[2];an+=base2.y*base1.z*cn[3];an+=base3.y*base0.z*cn[4];an+=base4.y*cn[5];an+=base0.x*base3.z*cn[6];an+=base0.x*base0.y*base2.z*cn[7];an+=base0.x*base1.y*base1.z*cn[8];an+=base0.x*base2.y*base0.z*cn[9];an+=base0.x*base3.y*cn[10];an+=base1.x*base2.z*cn[11];an+=base1.x*base0.y*base1.z*cn[12];an+=base1.x*base1.y*base0.z*cn[13];an+=base1.x*base2.y*cn[14];an+=base2.x*base1.z*cn[15];an+=base2.x*base0.y*base0.z*cn[16];an+=base2.x*base1.y*cn[17];an+=base3.x*base0.z*cn[18];an+=base3.x*base0.y*cn[19];an+=base4.x*cn[20];return an;
}

static __device__ fast_mpvec sum_5(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);an+=base5.z*cn[0];an+=base0.y*base4.z*cn[1];an+=base1.y*base3.z*cn[2];an+=base2.y*base2.z*cn[3];an+=base3.y*base1.z*cn[4];an+=base4.y*base0.z*cn[5];an+=base5.y*cn[6];an+=base0.x*base4.z*cn[7];an+=base0.x*base0.y*base3.z*cn[8];an+=base0.x*base1.y*base2.z*cn[9];an+=base0.x*base2.y*base1.z*cn[10];an+=base0.x*base3.y*base0.z*cn[11];an+=base0.x*base4.y*cn[12];an+=base1.x*base3.z*cn[13];an+=base1.x*base0.y*base2.z*cn[14];an+=base1.x*base1.y*base1.z*cn[15];an+=base1.x*base2.y*base0.z*cn[16];an+=base1.x*base3.y*cn[17];an+=base2.x*base2.z*cn[18];an+=base2.x*base0.y*base1.z*cn[19];an+=base2.x*base1.y*base0.z*cn[20];an+=base2.x*base2.y*cn[21];an+=base3.x*base1.z*cn[22];an+=base3.x*base0.y*base0.z*cn[23];an+=base3.x*base1.y*cn[24];an+=base4.x*base0.z*cn[25];an+=base4.x*base0.y*cn[26];an+=base5.x*cn[27];return an;
}

static __device__ fast_mpvec sum_6(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);fast_mpvec base6(base5.x*r.x,base5.y*r.y,base5.z*r.z);an+=base6.z*cn[0];an+=base0.y*base5.z*cn[1];an+=base1.y*base4.z*cn[2];an+=base2.y*base3.z*cn[3];an+=base3.y*base2.z*cn[4];an+=base4.y*base1.z*cn[5];an+=base5.y*base0.z*cn[6];an+=base6.y*cn[7];an+=base0.x*base5.z*cn[8];an+=base0.x*base0.y*base4.z*cn[9];an+=base0.x*base1.y*base3.z*cn[10];an+=base0.x*base2.y*base2.z*cn[11];an+=base0.x*base3.y*base1.z*cn[12];an+=base0.x*base4.y*base0.z*cn[13];an+=base0.x*base5.y*cn[14];an+=base1.x*base4.z*cn[15];an+=base1.x*base0.y*base3.z*cn[16];an+=base1.x*base1.y*base2.z*cn[17];an+=base1.x*base2.y*base1.z*cn[18];an+=base1.x*base3.y*base0.z*cn[19];an+=base1.x*base4.y*cn[20];an+=base2.x*base3.z*cn[21];an+=base2.x*base0.y*base2.z*cn[22];an+=base2.x*base1.y*base1.z*cn[23];an+=base2.x*base2.y*base0.z*cn[24];an+=base2.x*base3.y*cn[25];an+=base3.x*base2.z*cn[26];an+=base3.x*base0.y*base1.z*cn[27];an+=base3.x*base1.y*base0.z*cn[28];an+=base3.x*base2.y*cn[29];an+=base4.x*base1.z*cn[30];an+=base4.x*base0.y*base0.z*cn[31];an+=base4.x*base1.y*cn[32];an+=base5.x*base0.z*cn[33];an+=base5.x*base0.y*cn[34];an+=base6.x*cn[35];return an;
}

static __device__ fast_mpvec sum_7(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);fast_mpvec base6(base5.x*r.x,base5.y*r.y,base5.z*r.z);fast_mpvec base7(base6.x*r.x,base6.y*r.y,base6.z*r.z);an+=base7.z*cn[0];an+=base0.y*base6.z*cn[1];an+=base1.y*base5.z*cn[2];an+=base2.y*base4.z*cn[3];an+=base3.y*base3.z*cn[4];an+=base4.y*base2.z*cn[5];an+=base5.y*base1.z*cn[6];an+=base6.y*base0.z*cn[7];an+=base7.y*cn[8];an+=base0.x*base6.z*cn[9];an+=base0.x*base0.y*base5.z*cn[10];an+=base0.x*base1.y*base4.z*cn[11];an+=base0.x*base2.y*base3.z*cn[12];an+=base0.x*base3.y*base2.z*cn[13];an+=base0.x*base4.y*base1.z*cn[14];an+=base0.x*base5.y*base0.z*cn[15];an+=base0.x*base6.y*cn[16];an+=base1.x*base5.z*cn[17];an+=base1.x*base0.y*base4.z*cn[18];an+=base1.x*base1.y*base3.z*cn[19];an+=base1.x*base2.y*base2.z*cn[20];an+=base1.x*base3.y*base1.z*cn[21];an+=base1.x*base4.y*base0.z*cn[22];an+=base1.x*base5.y*cn[23];an+=base2.x*base4.z*cn[24];an+=base2.x*base0.y*base3.z*cn[25];an+=base2.x*base1.y*base2.z*cn[26];an+=base2.x*base2.y*base1.z*cn[27];an+=base2.x*base3.y*base0.z*cn[28];an+=base2.x*base4.y*cn[29];an+=base3.x*base3.z*cn[30];an+=base3.x*base0.y*base2.z*cn[31];an+=base3.x*base1.y*base1.z*cn[32];an+=base3.x*base2.y*base0.z*cn[33];an+=base3.x*base3.y*cn[34];an+=base4.x*base2.z*cn[35];an+=base4.x*base0.y*base1.z*cn[36];an+=base4.x*base1.y*base0.z*cn[37];an+=base4.x*base2.y*cn[38];an+=base5.x*base1.z*cn[39];an+=base5.x*base0.y*base0.z*cn[40];an+=base5.x*base1.y*cn[41];an+=base6.x*base0.z*cn[42];an+=base6.x*base0.y*cn[43];an+=base7.x*cn[44];return an;
}

static __device__ fast_mpvec sum_8(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);fast_mpvec base6(base5.x*r.x,base5.y*r.y,base5.z*r.z);fast_mpvec base7(base6.x*r.x,base6.y*r.y,base6.z*r.z);fast_mpvec base8(base7.x*r.x,base7.y*r.y,base7.z*r.z);an+=base8.z*cn[0];an+=base0.y*base7.z*cn[1];an+=base1.y*base6.z*cn[2];an+=base2.y*base5.z*cn[3];an+=base3.y*base4.z*cn[4];an+=base4.y*base3.z*cn[5];an+=base5.y*base2.z*cn[6];an+=base6.y*base1.z*cn[7];an+=base7.y*base0.z*cn[8];an+=base8.y*cn[9];an+=base0.x*base7.z*cn[10];an+=base0.x*base0.y*base6.z*cn[11];an+=base0.x*base1.y*base5.z*cn[12];an+=base0.x*base2.y*base4.z*cn[13];an+=base0.x*base3.y*base3.z*cn[14];an+=base0.x*base4.y*base2.z*cn[15];an+=base0.x*base5.y*base1.z*cn[16];an+=base0.x*base6.y*base0.z*cn[17];an+=base0.x*base7.y*cn[18];an+=base1.x*base6.z*cn[19];an+=base1.x*base0.y*base5.z*cn[20];an+=base1.x*base1.y*base4.z*cn[21];an+=base1.x*base2.y*base3.z*cn[22];an+=base1.x*base3.y*base2.z*cn[23];an+=base1.x*base4.y*base1.z*cn[24];an+=base1.x*base5.y*base0.z*cn[25];an+=base1.x*base6.y*cn[26];an+=base2.x*base5.z*cn[27];an+=base2.x*base0.y*base4.z*cn[28];an+=base2.x*base1.y*base3.z*cn[29];an+=base2.x*base2.y*base2.z*cn[30];an+=base2.x*base3.y*base1.z*cn[31];an+=base2.x*base4.y*base0.z*cn[32];an+=base2.x*base5.y*cn[33];an+=base3.x*base4.z*cn[34];an+=base3.x*base0.y*base3.z*cn[35];an+=base3.x*base1.y*base2.z*cn[36];an+=base3.x*base2.y*base1.z*cn[37];an+=base3.x*base3.y*base0.z*cn[38];an+=base3.x*base4.y*cn[39];an+=base4.x*base3.z*cn[40];an+=base4.x*base0.y*base2.z*cn[41];an+=base4.x*base1.y*base1.z*cn[42];an+=base4.x*base2.y*base0.z*cn[43];an+=base4.x*base3.y*cn[44];an+=base5.x*base2.z*cn[45];an+=base5.x*base0.y*base1.z*cn[46];an+=base5.x*base1.y*base0.z*cn[47];an+=base5.x*base2.y*cn[48];an+=base6.x*base1.z*cn[49];an+=base6.x*base0.y*base0.z*cn[50];an+=base6.x*base1.y*cn[51];an+=base7.x*base0.z*cn[52];an+=base7.x*base0.y*cn[53];an+=base8.x*cn[54];return an;
}

static __device__ fast_mpvec sum_J2(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_real zbase;an.z+=base2.z*cn[0].z;an.y+=base0.y*base1.z*cn[1].y;an.z+=base1.y*base0.z*cn[2].z;an.y+=base2.y*cn[3].y;an.x+=base0.x*base1.z*cn[4].x;an.x+=base0.x*base1.y*cn[6].x;an.z+=base1.x*base0.z*cn[7].z;an.y+=base1.x*base0.y*cn[8].y;an.x+=base2.x*cn[9].x;return an;
}

static __device__ fast_mpvec sum_J3(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_real zbase;an.z+=base3.z*cn[0].z;an.y+=base0.y*base2.z*cn[1].y;an.z+=base1.y*base1.z*cn[2].z;an.y+=base2.y*base0.z*cn[3].y;an.z+=base3.y*cn[4].z;an.x+=base0.x*base2.z*cn[5].x;an.x+=base0.x*base1.y*base0.z*cn[7].x;an.z+=base1.x*base1.z*cn[9].z;an.y+=base1.x*base0.y*base0.z*cn[10].y;an.z+=base1.x*base1.y*cn[11].z;an.x+=base2.x*base0.z*cn[12].x;an.z+=base3.x*cn[14].z;return an;
}

static __device__ fast_mpvec sum_J4(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_real zbase;an.z+=base4.z*cn[0].z;an.y+=base0.y*base3.z*cn[1].y;an.z+=base1.y*base2.z*cn[2].z;an.y+=base2.y*base1.z*cn[3].y;an.z+=base3.y*base0.z*cn[4].z;an.y+=base4.y*cn[5].y;an.x+=base0.x*base3.z*cn[6].x;an.x+=base0.x*base1.y*base1.z*cn[8].x;an.x+=base0.x*base3.y*cn[10].x;an.z+=base1.x*base2.z*cn[11].z;an.y+=base1.x*base0.y*base1.z*cn[12].y;an.z+=base1.x*base1.y*base0.z*cn[13].z;an.y+=base1.x*base2.y*cn[14].y;an.x+=base2.x*base1.z*cn[15].x;an.x+=base2.x*base1.y*cn[17].x;an.z+=base3.x*base0.z*cn[18].z;an.y+=base3.x*base0.y*cn[19].y;an.x+=base4.x*cn[20].x;return an;
}

static __device__ fast_mpvec sum_J5(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);fast_real zbase;an.z+=base5.z*cn[0].z;an.y+=base0.y*base4.z*cn[1].y;an.z+=base1.y*base3.z*cn[2].z;an.y+=base2.y*base2.z*cn[3].y;an.z+=base3.y*base1.z*cn[4].z;an.y+=base4.y*base0.z*cn[5].y;an.z+=base5.y*cn[6].z;an.x+=base0.x*base4.z*cn[7].x;an.x+=base0.x*base1.y*base2.z*cn[9].x;an.x+=base0.x*base3.y*base0.z*cn[11].x;an.z+=base1.x*base3.z*cn[13].z;an.y+=base1.x*base0.y*base2.z*cn[14].y;an.z+=base1.x*base1.y*base1.z*cn[15].z;an.y+=base1.x*base2.y*base0.z*cn[16].y;an.z+=base1.x*base3.y*cn[17].z;an.x+=base2.x*base2.z*cn[18].x;an.x+=base2.x*base1.y*base0.z*cn[20].x;an.z+=base3.x*base1.z*cn[22].z;an.y+=base3.x*base0.y*base0.z*cn[23].y;an.z+=base3.x*base1.y*cn[24].z;an.x+=base4.x*base0.z*cn[25].x;an.z+=base5.x*cn[27].z;return an;
}

static __device__ fast_mpvec sum_J6(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);fast_mpvec base6(base5.x*r.x,base5.y*r.y,base5.z*r.z);fast_real zbase;an.z+=base6.z*cn[0].z;an.y+=base0.y*base5.z*cn[1].y;an.z+=base1.y*base4.z*cn[2].z;an.y+=base2.y*base3.z*cn[3].y;an.z+=base3.y*base2.z*cn[4].z;an.y+=base4.y*base1.z*cn[5].y;an.z+=base5.y*base0.z*cn[6].z;an.y+=base6.y*cn[7].y;an.x+=base0.x*base5.z*cn[8].x;an.x+=base0.x*base1.y*base3.z*cn[10].x;an.x+=base0.x*base3.y*base1.z*cn[12].x;an.x+=base0.x*base5.y*cn[14].x;an.z+=base1.x*base4.z*cn[15].z;an.y+=base1.x*base0.y*base3.z*cn[16].y;an.z+=base1.x*base1.y*base2.z*cn[17].z;an.y+=base1.x*base2.y*base1.z*cn[18].y;an.z+=base1.x*base3.y*base0.z*cn[19].z;an.y+=base1.x*base4.y*cn[20].y;an.x+=base2.x*base3.z*cn[21].x;an.x+=base2.x*base1.y*base1.z*cn[23].x;an.x+=base2.x*base3.y*cn[25].x;an.z+=base3.x*base2.z*cn[26].z;an.y+=base3.x*base0.y*base1.z*cn[27].y;an.z+=base3.x*base1.y*base0.z*cn[28].z;an.y+=base3.x*base2.y*cn[29].y;an.x+=base4.x*base1.z*cn[30].x;an.x+=base4.x*base1.y*cn[32].x;an.z+=base5.x*base0.z*cn[33].z;an.y+=base5.x*base0.y*cn[34].y;an.x+=base6.x*cn[35].x;return an;
}

static __device__ fast_mpvec sum_J7(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);fast_mpvec base6(base5.x*r.x,base5.y*r.y,base5.z*r.z);fast_mpvec base7(base6.x*r.x,base6.y*r.y,base6.z*r.z);fast_real zbase;an.z+=base7.z*cn[0].z;an.y+=base0.y*base6.z*cn[1].y;an.z+=base1.y*base5.z*cn[2].z;an.y+=base2.y*base4.z*cn[3].y;an.z+=base3.y*base3.z*cn[4].z;an.y+=base4.y*base2.z*cn[5].y;an.z+=base5.y*base1.z*cn[6].z;an.y+=base6.y*base0.z*cn[7].y;an.z+=base7.y*cn[8].z;an.x+=base0.x*base6.z*cn[9].x;an.x+=base0.x*base1.y*base4.z*cn[11].x;an.x+=base0.x*base3.y*base2.z*cn[13].x;an.x+=base0.x*base5.y*base0.z*cn[15].x;an.z+=base1.x*base5.z*cn[17].z;an.y+=base1.x*base0.y*base4.z*cn[18].y;an.z+=base1.x*base1.y*base3.z*cn[19].z;an.y+=base1.x*base2.y*base2.z*cn[20].y;an.z+=base1.x*base3.y*base1.z*cn[21].z;an.y+=base1.x*base4.y*base0.z*cn[22].y;an.z+=base1.x*base5.y*cn[23].z;an.x+=base2.x*base4.z*cn[24].x;an.x+=base2.x*base1.y*base2.z*cn[26].x;an.x+=base2.x*base3.y*base0.z*cn[28].x;an.z+=base3.x*base3.z*cn[30].z;an.y+=base3.x*base0.y*base2.z*cn[31].y;an.z+=base3.x*base1.y*base1.z*cn[32].z;an.y+=base3.x*base2.y*base0.z*cn[33].y;an.z+=base3.x*base3.y*cn[34].z;an.x+=base4.x*base2.z*cn[35].x;an.x+=base4.x*base1.y*base0.z*cn[37].x;an.z+=base5.x*base1.z*cn[39].z;an.y+=base5.x*base0.y*base0.z*cn[40].y;an.z+=base5.x*base1.y*cn[41].z;an.x+=base6.x*base0.z*cn[42].x;an.z+=base7.x*cn[44].z;return an;
}

static __device__ fast_mpvec sum_J8(const fast_mpvec &r,const fast_mpvec *cn){
    fast_mpvec base0(r),an(0);fast_mpvec base1(base0.x*r.x,base0.y*r.y,base0.z*r.z);fast_mpvec base2(base1.x*r.x,base1.y*r.y,base1.z*r.z);fast_mpvec base3(base2.x*r.x,base2.y*r.y,base2.z*r.z);fast_mpvec base4(base3.x*r.x,base3.y*r.y,base3.z*r.z);fast_mpvec base5(base4.x*r.x,base4.y*r.y,base4.z*r.z);fast_mpvec base6(base5.x*r.x,base5.y*r.y,base5.z*r.z);fast_mpvec base7(base6.x*r.x,base6.y*r.y,base6.z*r.z);fast_mpvec base8(base7.x*r.x,base7.y*r.y,base7.z*r.z);fast_real zbase;an.z+=base8.z*cn[0].z;an.y+=base0.y*base7.z*cn[1].y;an.z+=base1.y*base6.z*cn[2].z;an.y+=base2.y*base5.z*cn[3].y;an.z+=base3.y*base4.z*cn[4].z;an.y+=base4.y*base3.z*cn[5].y;an.z+=base5.y*base2.z*cn[6].z;an.y+=base6.y*base1.z*cn[7].y;an.z+=base7.y*base0.z*cn[8].z;an.y+=base8.y*cn[9].y;an.x+=base0.x*base7.z*cn[10].x;an.x+=base0.x*base1.y*base5.z*cn[12].x;an.x+=base0.x*base3.y*base3.z*cn[14].x;an.x+=base0.x*base5.y*base1.z*cn[16].x;an.x+=base0.x*base7.y*cn[18].x;an.z+=base1.x*base6.z*cn[19].z;an.y+=base1.x*base0.y*base5.z*cn[20].y;an.z+=base1.x*base1.y*base4.z*cn[21].z;an.y+=base1.x*base2.y*base3.z*cn[22].y;an.z+=base1.x*base3.y*base2.z*cn[23].z;an.y+=base1.x*base4.y*base1.z*cn[24].y;an.z+=base1.x*base5.y*base0.z*cn[25].z;an.y+=base1.x*base6.y*cn[26].y;an.x+=base2.x*base5.z*cn[27].x;an.x+=base2.x*base1.y*base3.z*cn[29].x;an.x+=base2.x*base3.y*base1.z*cn[31].x;an.x+=base2.x*base5.y*cn[33].x;an.z+=base3.x*base4.z*cn[34].z;an.y+=base3.x*base0.y*base3.z*cn[35].y;an.z+=base3.x*base1.y*base2.z*cn[36].z;an.y+=base3.x*base2.y*base1.z*cn[37].y;an.z+=base3.x*base3.y*base0.z*cn[38].z;an.y+=base3.x*base4.y*cn[39].y;an.x+=base4.x*base3.z*cn[40].x;an.x+=base4.x*base1.y*base1.z*cn[42].x;an.x+=base4.x*base3.y*cn[44].x;an.z+=base5.x*base2.z*cn[45].z;an.y+=base5.x*base0.y*base1.z*cn[46].y;an.z+=base5.x*base1.y*base0.z*cn[47].z;an.y+=base5.x*base2.y*cn[48].y;an.x+=base6.x*base1.z*cn[49].x;an.x+=base6.x*base1.y*cn[51].x;an.z+=base7.x*base0.z*cn[52].z;an.y+=base7.x*base0.y*cn[53].y;an.x+=base8.x*cn[54].x;return an;
}


static __device__ fast_mpvec (*const sum_Jn_funlist[])(const fast_mpvec &,const fast_mpvec *)={
    nullptr,nullptr,
    sum_J2,
    sum_J3,
    sum_J4,
    sum_J5,
    sum_J6,
    sum_J7,
    sum_J8
};

static __device__ fast_mpvec (*const sum_n_funlist[])(const fast_mpvec &,const fast_mpvec *)={
    nullptr,nullptr,
    sum_2,
    sum_3,
    sum_4,
    sum_5,
    sum_6,
    sum_7,
    sum_8
};
#endif


__device__ fast_mpvec geopotential::cuda_sum(fast_real R,fast_mpvec r,int_t N_start,int_t N_end) const{
    fast_real x=r.x,y=r.y,z=r.z;
    fast_real rr2=1/(r%r);
    fast_real rr=sqrt(rr2);
    fast_real R_r2=R*rr2;

    if(N_start<2)N_start=2;
    if(N_end<0||N_end>N)N_end=N;

    fast_mpvec a(0);
    for(int_t n=N_end;n>=1;--n){
        if(n>=N_start){
            a+=(n>Nt?sum_Jn_funlist:sum_n_funlist)[n](r,c_table+Precompute_Table_size(n-1));
        }
        a*=R_r2;
    }
    a*=-rr2*rr;

    return a;
}
#endif

//copied from ring.cpp, subroutines are added ``static __device__''
#if 1

#define CONST_TABLE static __device__ const
#include"disk_approx.impl"
static __device__ float fzcorr(float x,float y,int kh){
    const size_t N=32;

    x*=N;
    y*=N;
    int i=floor(x),j=floor(y);
    j-=kh*N;
    if(i<0)i=0;
    if(j<0)j=0;
    if(i>=N)i=N-1;
    if(j>=N)j=N-1;
    j+=kh*N;
    x-=i;
    y-=j;
    const float *p=(const float *)disk_approx_table+16*(2*N*i+j);
    return  p[0]+x*(p[4]+x*(p[8]+x*p[12]))
        +y*(p[1]+x*(p[5]+x*(p[9]+x*p[13]))
            +y*(p[2]+x*(p[6]+x*(p[10]+x*p[14]))
                +y*(p[3]+x*(p[7]+x*(p[11]+x*p[15])))));
}
static __device__ double padesum(double x,const double coef[][2],int n){
    double resn=0,resd=0;
    while(n>0){
        --n;
        resn=x*resn+coef[n][0];
        resd=x*resd+coef[n][1];
    }
    return resn/resd;
}
#define PadeSum(x,coef) padesum(x,coef,sizeof(coef)/sizeof(coef[0]))

__device__ fast_mpvec ring::cuda_sum(fast_mpvec r) const{
    const double pi=3.1415926535897932;
    const double xserkk[][2]={{-4228.1846570321674587,6387.9414439779279738},{-13276.372802884782407,16984.953742612298093},{-11316.549465243054011,8997.7319379726254196},{-3242.2867356286923985,403.15182320396046619},{-704.60633921130372399,-5.7789477668119521439}};
    const double xserkm[][2]={{-2530.8255093002657219,7950.8733392382639504},{-5471.7803365810451727,17189.101860053144376},{-1624.0389725420726306,7100.5562444641275465},{1164.6058458813108946,532.54915854738706820},{270.03897254207263065,-5.0806023029229406554}};
    const double xser0[][2]={{3216.9908772759482762,32768.000000000000000},{-6727.1375380948220662,-93098.060290997395032},{4534.8257724052420971,96814.901233532309690},{-1004.5251977686921953,-43973.540867249017404},{17.915064837735955511,7754.0323173314244943},{-0.039157108395559676838,-224.00490789882410101}};

    const double kserkk[][2]={{-2373.3522804139070577,1124.9866194216978637},{-8300.7524798211745830,17215.217000406477760},{15279.068333441193395,26590.060415743810476},{-4354.4701057632348105,-9776.5071563707548874},{-250.49346744287694398,-2385.7568792012312128}};
    const double kserkm[][2]={{-523.02599873176918199,290.10526315789473684},{-5485.2200380469245403,5336.9790741915028535},{-8338.1052631578947368,15206.361445783132530},{-2017.9378566899175650,10515.743817374762207},{-19.710843373493975904,1418.8103994927076728}};
    const double kser0[][2]={{51471.854036415172419,32768.000000000000000},{-119621.18573838632696,-84345.212035111671234},{96754.491735189724987,78074.125512775430927},{-31810.735585315211544,-31108.829182924968396},{3678.1383813624327949,4926.4461273804873303},{-74.040670509817904983,-206.66608791871105442}};

    bool psign;
    if(psign=r.z<0)r.z=-r.z;

    const double re2(r.x*r.x+r.y*r.y);
    const double re=sqrt(re2);
    const double rz2(r.z*r.z);
    const double r2=re2+rz2;
    const double rz(r.z);

    fast_mpvec ret(0);

    for(int_t i=0;i<N;++i){
        double R=c_table[i].R;
        double H=c_table[i].H;

        double reR=re*R;
        double R2=R*R,H2,fzh=1;
        double k2dre=4/(r2+R2+2*reR);
        if(H>0){
            H2=H*H;
            k2dre-=H2*(k2dre*k2dre)*(1-k2dre*rz2)/48;
            fzh=(re>R?re-R:0);
            fzh=4*(rz2+fzh*fzh);
            fzh=sqrt(fzh/(fzh+H2));
        }
        double k2=reR*k2dre;
        k2dre*=R2*sqrt(k2dre);
        double phi=atan2(rz,re-R);
        double k,lk,fx;
        bool lkh=k2>=0.45;
        if(lkh){
            k=sqrt(k2);lk=log((1-k2)/(8*(1+k)));
            fx=(PadeSum(k,xserkk)+lk*PadeSum(k,xserkm))/(k2*k2);
        }
        else{
            fx=PadeSum(k2,xser0);
        }
        fx*=k2dre;

        //approx for fz=diskfz_approx(phi,k2)
        double fz;
        double s=sin(phi),c=cos(phi);
        bool kh=k2>=.5;
        bool exphi;

        if(kh){
            const double phimax=(pi/2);
            exphi=c>0;
            double ck=1-k2,lk=log(ck/16),sk=sqrt(ck);
            double cc=exphi?-c:c;
            double cphi=exphi?pi-phi:phi;
            fz=cphi-(s/8)*((lk*(k2-5)-2*ck)*sk+2*cc*(1+lk)*ck);
            fz*=fzcorr((exphi?phi:pi-phi)*(1/phimax),k2*2,kh);
        }
        else{
            const double phimax=(809./512);
            exphi=phi>phimax;
            double cr2=exphi?R2+rz2:r2;
            double cRdr2=(exphi?re2:R2)/cr2;
            fz=(pi/8)*(4-3*cRdr2)*(rz/sqrt(cr2)*cRdr2);
            fz*=fzcorr((exphi?pi-phi:phi)*(1/phimax),k2*2,kh);
        }

        if(exphi){
            double ek;
            if(lkh){
                ek=PadeSum(k,kserkk)+lk*PadeSum(k,kserkm);
            }
            else{
                ek=PadeSum(k2,kser0);
            }
            fz=pi-2*sqrt(1-k2)*ek*s-fz;
        }

        double Gs=c_table[i].Gs;
        fx*=Gs;
        fz*=Gs;

        ret+=fast_mpvec(-4*fx*r.x,-4*fx*r.y,-2*fz*fzh);
    }

    if(psign)ret.z=-ret.z;

    return ret;
}
#endif

#include"RungeKutta.impl"


extern __shared__ char sharedMem[];
void __device__ accel_0(){//deform

    const fast_real c=299792458;
    const fast_real c2=c*c;

    int i0=blockIdx.x*dkf.mass_per_block;
    for(int di=0;di<dkf.mass_per_block;++di){
        int i=di+i0;
        if(i<dkf.nmass){
            mass &mi=dkf.dmlist[i];

            int_t max_iter=4;
            do{
                maccel_1 *tpmi=(maccel_1 *)sharedMem+threadIdx.x;
                tpmi[0].phi=0;
                tpmi[0].naccel=0;
                tpmi[0].C_potential=0;
                for(int dj=0;dj<dkf.nmass;dj+=blockDim.x){
                    int j=dj+threadIdx.x;
                    if(j<dkf.nmass&&i!=j){
                        mass &mj=dkf.dmlist[j];
                        //---accel-prepare
                        fast_mpvec r=mj.r-mi.r;
                        fast_mpvec v=mj.v-mi.v;
                        fast_real rr2=1/(r%r);
                        fast_real rr=sqrt(rr2);
                        fast_real tp_dphi=rr*mj.GM;
                        fast_real tp_dg=rr2*tp_dphi;

                        tpmi[0].phi-=tp_dphi;
                        tpmi[0].naccel+=tp_dg*r;

                        //damped tidal deformation matrix
                        fast_real fmj=-3*tp_dg;
                        fast_mpvec dw=r*(mi.w*r-v)*rr2;
                        fast_real dw2=dw%dw,dw1=sqrt(dw2),dwt=dw1*mi.tide_delay*mj.tide_delay_factor;
                        fast_real
                            ecwt2=1+4*dwt*dwt,
                            secwt2=sqrt(ecwt2),
                            x0=sqrt((1+secwt2)/(2*ecwt2)),
                            y0=dwt*sqrt(2/(ecwt2*(1+secwt2))),
                            z02=dwt*dwt*(2/(ecwt2+secwt2));
                        fast_mpvec
                            dr=r*x0+dw*r*(y0/dw1);
                        tpmi[0].C_potential+=(
                            fast_mpmat(fast_real(1)/3-z02)
                            -fast_mpmat(dr*rr2,dr)
                            +fast_mpmat(dw*(z02/dw2),dw)
                            )*(fmj*mi.k2);
                    }
                }

                for(int wing=blockDim.x>>1;wing>1;wing>>=1){
                    __syncthreads();
                    if(threadIdx.x<wing){
                        tpmi[0].phi+=tpmi[wing].phi;
                        tpmi[0].naccel+=tpmi[wing].naccel;
                        tpmi[0].C_potential+=tpmi[wing].C_potential;
                    }
                }
                __syncthreads();
                if(threadIdx.x<1){
                    mi.phi=tpmi[0].phi+tpmi[1].phi;
                    mi.naccel=tpmi[0].naccel+tpmi[1].naccel;
                    mi.C_potential=tpmi[0].C_potential+tpmi[1].C_potential;
                }
                __syncthreads();

                __shared__ bool should_break;
                //angular accelerate
                if(threadIdx.x==0){
                    fast_real w2=mi.w%mi.w;
                    mi.C_potential+=fast_mpmat(w2*mi.k2r/3)-fast_mpmat(mi.w*mi.k2r,mi.w);
                    mi.C_potential*=mi.R*mi.R2/(2*mi.GM);

                    fast_mpmat fmis(mi.s);
                    fast_mpmat mc(mi.C_static);
                    mc.x.x+=mi.exJ2/2;
                    mc.y.y+=mi.exJ2/2;
                    mc.z.z-=mi.exJ2;
                    mc=fmis.toworld(mc);
                    mi.C_potential+=mc;
                    mi.GI=2*mi.R2/3*(fast_mpmat(mi.A)-mi.C_potential);
                    fast_mpvec oldw=mi.w;
                    mi.w=mi.GI.inverse()%(fast_mpvec(mi.GL));
                    should_break=((mi.w-oldw).norm()<1e-9*oldw.norm());
                }
                __syncthreads();
                if(should_break)break;
            } while(--max_iter);
            if(threadIdx.x==0){
                mi.beta=fast_mpvec(mi.v)/c;
                mi.beta2=mi.beta%mi.beta;
                mi.phi/=c2;
                mi.naccel/=c2;
            }
        }
    }
}

void __device__ accel_1(){//accel
    const fast_real c=299792458;
    const fast_real c2=c*c;

    int i0=blockIdx.x*dkf.mass_per_block;
    for(int di=0;di<dkf.mass_per_block;++di){
        int i=di+i0;
        if(i<dkf.nmass){
            mass &mi=dkf.dmlist[i];
            maccel_2 *tpmi=(maccel_2*)sharedMem+threadIdx.x;
            tpmi[0].gaccel=0;
            tpmi[0].daccel=0;
            tpmi[0].dtorque=0;
            tpmi[0].min_distance=0;
            tpmi[0].max_influence=0;
            for(int dj=0;dj<dkf.nmass;dj+=blockDim.x){
                int j=dj+threadIdx.x;
                if(j<dkf.nmass&&i!=j){
                    mass &mj=dkf.dmlist[j];
                    //---accel
                    fast_mpvec r=mj.r-mi.r;
                    fast_real rr2=1/(r%r);
                    fast_real rr=sqrt(rr2);
                    fast_real rr3=rr*rr2;

                    //start post-newtonian correction
                    fast_real tp_dphi=rr*mj.GM;
                    fast_real tp_dg=rr2*tp_dphi;
                    fast_real rbj=r%mj.beta;
                    fast_real tp_rbjrr2=rbj*rr;
                    tp_rbjrr2*=tp_rbjrr2;
                    fast_real delta1=4*mi.phi+mj.phi+mi.beta2+2*mj.beta2-4*(mi.beta%mj.beta)+(r%mj.naccel-3*tp_rbjrr2)/2;
                    fast_mpvec b=mj.beta-mi.beta;
                    tpmi[0].gaccel+=7*tp_dphi/2*mj.naccel+tp_dg*delta1*r+tp_dg*(rbj-r%b*4)*b;
                    //end post-newtonian correction
                    tpmi[0].min_distance=cuda_max(tpmi[0].min_distance,rr);
                    tpmi[0].max_influence=cuda_max(tpmi[0].max_influence,tp_dg);
                    //rotational & tidal deformation: gravity, torque
                    fast_mpvec Cr=mi.C_potential%r;
                    fast_mpvec dg=mj.GM*rr3*rr2*mi.R2*(r%Cr*5*rr2*r-(Cr+Cr));
                    tpmi[0].daccel+=dg;
                    tpmi[0].dtorque+=r*dg;
                    //to avoid reduce cross thread block, re-calculate daccel instead of using anti-force
                    fast_mpvec Cjr=mj.C_potential%r;
                    tpmi[0].daccel+=mj.GM*rr3*rr2*mj.R2*(r%Cjr*5*rr2*r-(Cjr+Cjr));
                    //start lense thirring
                    fast_real rcr32=2*mj.GM/c*rr3;
                    fast_mpvec GL=mi.GL;
                    fast_mpvec GLb=GL*b;
                    tpmi[0].daccel-=rcr32*(GLb-3*(GLb%r)*rr2*r);
                    tpmi[0].dtorque-=rcr32*(GL*(r*b));
                    fast_mpvec GLj=mj.GL;
                    tpmi[0].daccel+=(rcr32*(GLj-3*(GLj%r)*rr2*r))*b;
                    //end lense thirring

                    //start radiation pressure
                    tpmi[0].daccel+=mj.lum*mi.rR2_4Mc*rr3*r;
                    //end radiation pressure

                    //---accel
                    //higher harmonics will be done later
                }
            }

            for(int wing=blockDim.x>>1;wing>1;wing>>=1){
                __syncthreads();
                if(threadIdx.x<wing){
                    tpmi[0].gaccel+=tpmi[wing].gaccel;
                    tpmi[0].daccel+=tpmi[wing].daccel;
                    tpmi[0].dtorque+=tpmi[wing].dtorque;
                    tpmi[0].min_distance=cuda_max(tpmi[0].min_distance,tpmi[wing].min_distance);
                    tpmi[0].max_influence=cuda_max(tpmi[0].max_influence,tpmi[wing].max_influence);
                }
            }
            __syncthreads();
            if(threadIdx.x<1){
                mi.gaccel=tpmi[0].gaccel+tpmi[1].gaccel;
                mi.daccel=tpmi[0].daccel+tpmi[1].daccel;
                mi.dtorque=tpmi[0].dtorque+tpmi[1].dtorque;
                tpmi[0].min_distance=cuda_max(tpmi[0].min_distance,tpmi[1].min_distance);
                mi.min_distance=cuda_max(mi.min_distance,tpmi[0].min_distance);
                tpmi[0].max_influence=cuda_max(tpmi[0].max_influence,tpmi[1].max_influence);
                mi.max_influence=cuda_max(mi.max_influence,tpmi[0].max_influence);
            }
            __syncthreads();
        }
    }
}

void __device__ accel_2(){//higher harmonics
    int mn=dkf.nmass;
    int i0=blockIdx.x*dkf.mass_per_block;

    for(int j=0;j<dkf.nmass;++j){
        mass &mi=dkf.dmlist[j];
        if(mi.gpmodel||mi.ringmodel){
            maccel_2 *tpmi=(maccel_2*)sharedMem+threadIdx.x;
            tpmi[0].daccel=0;
            tpmi[0].dtorque=0;

            for(int di=0;di<dkf.mass_per_block;di+=blockDim.x){
                int i=i0+di+threadIdx.x;
                if(i<mn&&di+threadIdx.x<dkf.mass_per_block&&i!=j){
                    mass &mj=dkf.dmlist[i];
                    fast_mpvec an(0);
                    fast_mpvec r=mj.r-mi.r;
                    if(mi.gpmodel){
                        fast_mpmat fmis(mi.s);
                        fast_mpvec lr=fmis.tolocal(r);
                        an+=fmis.toworld(mi.gpmodel->cuda_sum(mi.R,lr));
                    }
                    if(mi.ringmodel){
                        fast_mpvec migl=mi.GL;
                        fast_mpmat fgls(migl.perpunit(),0,migl/migl.norm());
                        fast_mpvec lr=fgls.tolocal(r);
                        an+=fgls.toworld(mi.ringmodel->cuda_sum(lr));
                    }

                    mj.daccel+=mi.GM*an;
                    tpmi[0].daccel-=mj.GM*an;
                    tpmi[0].dtorque-=mj.GM*(r*an);


                }
            }
            //the initial wing should be lift to 2's power and the mask should be adjusted.
            for(int wing=blockDim.x>>1;wing>1;wing>>=1){
                __syncthreads();
                if(threadIdx.x<wing){
                    tpmi[0].daccel+=tpmi[wing].daccel;
                    tpmi[0].dtorque+=tpmi[wing].dtorque;
                }
            }
            __syncthreads();
            if(threadIdx.x<1){
                //Note: by construction, nblocks <= nmass
                //      so we store partial acceleration per block in global mlist memory
                mass &mbi=dkf.dmlist[blockIdx.x];
                mbi.idaccel=tpmi[0].daccel+tpmi[1].daccel;
                mbi.idtorque=tpmi[0].dtorque+tpmi[1].dtorque;
            }
            __syncthreads();


            cooperative_groups::grid_group grid=cooperative_groups::this_grid();
            grid.sync();
            if(blockIdx.x==0){
                // block 0
                tpmi[0].daccel=0;
                tpmi[0].dtorque=0;
                for(int i=threadIdx.x;i<gridDim.x;i+=blockDim.x){
                    mass &mj=dkf.dmlist[i];
                    tpmi[0].daccel+=mj.idaccel;
                    tpmi[0].dtorque+=mj.idtorque;
                }

                for(int wing=blockDim.x>>1;wing>1;wing>>=1){
                    __syncthreads();
                    if(threadIdx.x<wing){
                        tpmi[0].daccel+=tpmi[wing].daccel;
                        tpmi[0].dtorque+=tpmi[wing].dtorque;
                    }
                }
                __syncthreads();
                if(threadIdx.x<1){
                    mi.daccel+=tpmi[0].daccel+tpmi[1].daccel;
                    mi.dtorque+=tpmi[0].dtorque+tpmi[1].dtorque;
                }
                __syncthreads();

            }
            grid.sync();
        }
    }
}

void __device__ Cuda_accel(){
    const fast_real c=299792458;
    const fast_real c2=c*c;
    cooperative_groups::grid_group grid=cooperative_groups::this_grid();

    mass *x=dkf.dmlist;
    int mn=dkf.nmass;
    int i0=blockIdx.x*dkf.mass_per_block;

    grid.sync();
    accel_0();
    grid.sync();
    accel_1();
    grid.sync();
    accel_2();
    grid.sync();
    for(int di=0;di<dkf.mass_per_block;di+=blockDim.x){
        int i=i0+di+threadIdx.x;
        if(i<mn&&di+threadIdx.x<dkf.mass_per_block){
            mass &mi=x[i];
            mi.phi*=c2;
            mi.naccel*=c2;

            //ring correction
            if(mi.ringmodel){
                mass &m=mi;
                ring &mr=*m.ringmodel;
                //ring has inertia that prevents extra accelerations
                m.daccel-=(m.gaccel+m.daccel+m.naccel)*mr.GM_ratio;
                //ring has angular momentum that prevents extra angular accelerations
                fast_mpvec mGL=m.GL;
                fast_real rmGL2=1/(mGL%mGL),rmGL=sqrt(rmGL2);
                fast_mpvec ptorque=m.dtorque;
                ptorque=ptorque%mGL*rmGL2*mGL;
                fast_mpvec dtorque=m.dtorque-ptorque;
                //ptorque(parallel to GL) will not change
                //  since this part has nothing to do with the ring
                //dtorque(perpendicular to GL) will decrease
                //  due to ring's angular momentum
                dtorque*=1/(1+rmGL*mr.GL);
                m.dtorque=dtorque+ptorque;
            }
        }
    }
}

void __global__ Cuda_RungeKutta_Kernel(){
    const real *clist=(const real *)rk12_coefs;

    mass_state *x0=dkf.x0,*f=dkf.f;
    mass *x=dkf.dmlist;
    fast_real dt=dkf.dt;

    int mn=dkf.nmass;
    int i0=blockIdx.x*dkf.mass_per_block;

    for(int_t i_step=0;i_step<dkf.n_step;++i_step){
        for(int di=0;di<dkf.mass_per_block;di+=blockDim.x){
            int i=i0+di+threadIdx.x;
            if(i<mn&&di+threadIdx.x<dkf.mass_per_block){
                x0[i]=x[i];
            }
        }

        fast_real dt_k;
        for(int_t k=1;k<=25;++k){
            dt_k=0;
            for(int di=0;di<dkf.mass_per_block;di+=blockDim.x){
                int i=i0+di+threadIdx.x;
                if(i<mn&&di+threadIdx.x<dkf.mass_per_block){
                    mass_state &fi=f[mn*(k-1)+i];
                    mass &xi=x[i];
                    fi.v=xi.v;
                    fi.naccel=xi.gaccel+xi.daccel+xi.naccel;
                    fi.w=xi.w;
                    fi.dtorque=xi.dtorque;

                    if(k==25){
                        fast_mpvec j(xi.GL),&w=xi.w;
                        fast_mpvec wxj=w*j;
                        fast_real wxj2=wxj%wxj;
                        if(wxj2!=0){
                            fast_mpmat ir=xi.GI.inverse();
                            fast_real dt2=dt*dt;
                            fast_real e=w%j*(fast_real(1)/2);
                            fast_real de=(ir%w)%((ir%wxj)*j+w*wxj);
                            de*=e*dt2*dt2/(36*wxj2);
                            xi.Erot=de;
                            xi.Egrad=wxj;
                        }
                        else{
                            xi.Erot=0;
                            xi.Egrad=0;
                        }
                    }

                    xi.v=0;
                    xi.naccel=0;
                    xi.w=0;
                    xi.dtorque=0;
                }
            }

            for(int_t j=0;j<k;++j){
                const real &ckj=clist[k*(k-1)/2+j];
                const fast_real fckj=(fast_real)ckj;
                if(ckj.hi){
                    dt_k+=fckj;
                    for(int di=0;di<dkf.mass_per_block;di+=blockDim.x){
                        int i=i0+di+threadIdx.x;
                        if(i<mn&&di+threadIdx.x<dkf.mass_per_block){
                            mass_state &fi=f[mn*j+i];
                            mass &xi=x[i];
                            xi.v+=ckj*fi.v;
                            xi.naccel+=fckj*fi.naccel;
                            xi.w+=fckj*fi.w;
                            xi.dtorque+=fckj*fi.dtorque;
                        }
                    }
                }
            }
            
            fast_real t=fast_real(dkf.t_eph)+dt*(dt_k+i_step);
            for(int di=0;di<dkf.mass_per_block;di+=blockDim.x){
                int i=i0+di+threadIdx.x;
                if(i<mn&&di+threadIdx.x<dkf.mass_per_block){
                    const mass_state &xi0=x0[i];
                    mass &xi=x[i];

                    xi.r=xi0.r+real(dt)*xi.v;
                    xi.v=xi0.v+mpvec(dt*xi.naccel);
                    xi.w/=dt_k;
                    xi.s=xi0.s;
                    fast_mpvec dw=(fast_real(2)/3)*(xi.w-xi0.w);
                    xi.s+=rotation_matrix(xi.w-dw,dt*dt_k/2)%fast_mpmat(xi.s);
                    xi.s+=rotation_matrix(xi.w+dw,dt*dt_k/2)%fast_mpmat(xi.s);
                    xi.GL=xi0.GL+mpvec(dt*xi.dtorque);

                    if(k==25){
                        xi.s+=rotation_matrix(xi.Egrad,xi.Erot)%fast_mpmat(xi.s);
                    }
                    //update
                    xi.GM=xi.GM0+xi.dGM*t;
            
                    xi.exJ2=xi.dJ2*t;
                }
            }

            Cuda_accel();
        }
    }
}

void msystem::Cuda_RungeKutta12(fast_real dt,int_t n_step){
    if(n_step<=0)return;
    gpdata_t mgp;
    ringdata_t mrg;
    cuda_rungekutta_kernel_config kf;
    kf.load(mlist,mgp,mrg,dt,n_step);
    kf.t_eph=t_eph;
    cudaMemcpyToSymbol(dkf,&kf,sizeof(kf));
    //Cuda_Kernel<<<kf.nblocks,kf.nthreads>>>();
    cudaLaunchCooperativeKernel(
        (void*)Cuda_RungeKutta_Kernel,
        dim3(kf.nblocks),
        dim3(kf.nthreads),
        nullptr,
        kf.nthreads*std::max(sizeof(maccel_1),sizeof(maccel_2))
    );
    cudaDeviceSynchronize();
    kf.save(mlist,mgp,mrg);
}

void __global__ Cuda_accel_Kernel(){
    Cuda_accel();
}

void msystem::Cuda_accel(){
    gpdata_t mgp;
    ringdata_t mrg;
    cuda_rungekutta_kernel_config kf;
    kf.load(mlist,mgp,mrg,0,0);
    kf.t_eph=t_eph;
    cudaMemcpyToSymbol(dkf,&kf,sizeof(kf));
    //Cuda_Kernel<<<kf.nblocks,kf.nthreads>>>();
    cudaLaunchCooperativeKernel(
        (void*)Cuda_RungeKutta_Kernel,
        dim3(kf.nblocks),
        dim3(kf.nthreads),
        nullptr,
        kf.nthreads*std::max(sizeof(maccel_1),sizeof(maccel_2))
    );
    cudaDeviceSynchronize();
    kf.save(mlist,mgp,mrg);
}
