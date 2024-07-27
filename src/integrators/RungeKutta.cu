#define INLINE __device__
#include"physics/mass.h"
#include"physics/geopotential.h"
#include"physics/ring.h"
#include"physics/mass.impl"
#include<stdlib.h>
#include<algorithm>
#include<cuda_runtime.h>
#include<cooperative_groups.h>
#include<mutex>

#define CUDA_IMPL

#define WARP_SIZE 32
#define CUDA_CORES 1920
#define MAXBLOCKS (1920/32)

#define cuda_max(a,b) ((a)>(b)?(a):(b))

std::mutex cuda_mutex;

//geopotential data, mlist[first].gpmodel==second
typedef std::vector<std::pair<int_t,const geopotential*>> gpdata_t;
//ring data, mlist[first].ringmodel==second
typedef std::vector<std::pair<int_t,const ring*>> ringdata_t;

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

#include"physics/geopotential.impl"
#include"physics/ring.impl"

extern __shared__ char sharedMem[];
void __device__ accel_0(){//deform

    const fast_real c=CONSTANT_VALUE_C;
    const fast_real c2=c*c;

    int i0=blockIdx.x*dkf.mass_per_block;
    for(int di=0;di<dkf.mass_per_block;++di){
        int i=di+i0;
        if(i<dkf.nmass){
            mass &mi=dkf.dmlist[i];

            int_t max_iter=MAX_ANGULAR_VELOCITY_ITER;
            do{
                maccel_1 *tpmi=(maccel_1 *)sharedMem+threadIdx.x;
                tpmi[0].phi=0;
                tpmi[0].naccel=0;
                tpmi[0].C_potential=0;
                for(int dj=0;dj<dkf.nmass;dj+=blockDim.x){
                    int j=dj+threadIdx.x;
                    if(j<dkf.nmass&&i!=j){
                        mass &mj=dkf.dmlist[j];
                        DAMPED_TIDAL_DEFORMATION_MATRIX(tpmi[0]);
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
                    UPDATE_HARMONICS;
                    UPDATE_ANGULAR_VELOCITY;
                }
                __syncthreads();
                if(should_break)break;
            } while(--max_iter);
            if(threadIdx.x==0){
                PREPARE_RELATIVITY;
            }
        }
    }
}

void __device__ accel_1(){//accel
    const fast_real c=CONSTANT_VALUE_C;

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
                    RELATIVITY(tpmi[0]);
                    tpmi[0].min_distance=cuda_max(tpmi[0].min_distance,rr);
                    tpmi[0].max_influence=cuda_max(tpmi[0].max_influence,tp_dg);
                    //to avoid reduce cross thread block, re-calculate daccel instead of using anti-force
                    ROTATIONAL_TIDAL_DEFORMATION_NANTI_FORCE(tpmi[0]);
                    LENSE_THIRRING(tpmi[0]);
                    RADIATION_PRESSURE(tpmi[0]);

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
                        an+=fgls.toworld(mi.ringmodel->cuda_sum(mi.R,lr));
                    }

                    APPLY_NONPOINT_FORCE(tpmi[0]);
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
    const fast_real c=CONSTANT_VALUE_C;
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
            FINALIZE_RELATIVITY;

            if(mi.ringmodel){
                RING_CORRECTION;
            }
        }
    }
}

#include"RungeKutta.impl"

void msystem::Cuda_RungeKutta12(fast_real dt,int_t n_step){
    if(n_step<=0)return;
    gpdata_t mgp;
    ringdata_t mrg;
    cuda_rungekutta_kernel_config kf;
    cuda_mutex.lock();
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
    cuda_mutex.unlock();
}

void __global__ Cuda_accel_Kernel(){
    Cuda_accel();
}

void msystem::Cuda_accel(){
    gpdata_t mgp;
    ringdata_t mrg;
    cuda_rungekutta_kernel_config kf;
    cuda_mutex.lock();
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
    cuda_mutex.unlock();
}
