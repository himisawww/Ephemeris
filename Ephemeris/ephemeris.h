#pragma once
#include<cstdint>
#include<cstdio>
#include<cmath>
#include<vector>
#include<map>
#include<string>
#include"dfloat_t.h"
typedef dfloat_t<double> real;

#define strtoreal atof

typedef double fast_real;

#include"vec_t.h"

typedef vec_t<double> vec;
typedef mat_t<double> mat;
typedef vec_t<real> mpvec;
typedef mat_t<real> mpmat;
typedef vec_t<fast_real> fast_mpvec;
typedef mat_t<fast_real> fast_mpmat;

const double pi=3.1415926535897932;
const double degree=pi/180;

typedef int64_t int_t;

//parameters for orbital ephemerides
//see https://bridge.kamine.cloud/archives/421
class ephem_orb{
public:
    //t = (time relative to epoch) * sqrt(GM)
    double t;
    //j = (angular momentum per mass) / sqrt(GM)
    vec j;
    //q = e-1 , defined to reduce numeric error at e~0 and e~1
    double q;
    //earg: argument of eccentricity vector in j-frame
    double earg;
    //m = (mean anomaly) / sqrt(1-e^2)^3
    //dm/dt = 1/j^3
    double m;

    ephem_orb()=default;
    //convert r and v to orbital parameters
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    ephem_orb(double t,const vec &r,const vec &v);

    //convert orbital parameters to r and v
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    void rv(double t,vec &r,vec &v) const;
};
//parameters for rotational ephemerides
class ephem_rot{
public:
    //t = time relative to epoch
    double t;
    //rotational angular velocity
    vec w;
    //polar wandering in surface-frame
    double ptheta,pphi;
    //rotation angle
    double angle;

    ephem_rot()=default;
    //convert surface-frame and w to rotational parameters
    ephem_rot(double t,const mat &s,const vec &w);

    //convert rotational parameters to surface-frame and w
    void sw(double t,mat &s,vec &w) const;
};

class barycen{
public:
    mpvec r,v;
    real GM;
    /*
    *   GMrv_sys is total/mean GMrv counting all children,
    *   and is used for total momentum of system hosted by barycen.
    *   GMrv is used for orbital center,
    *   for barycen, GMrv is total/mean of center binary's GMrv_sys
    *   for mass, GMrv is actual GMrv of center mass itself.
    */
    mpvec r_sys,v_sys;
    real GM_sys;

    //index of parent in blist
    //if no parent, -1
    int_t pid;
    //for barycen, index of host and guest in blist
    //if not barycen, -1
    int_t hid,gid;
    //id of barycen which represents the orbital motion
    //i.e. (recursive-)parent through hosts of binary
    int_t tid;
    //for mass, index of mass in mlist
    //for barycen, index of (recursive-)host mass
    int_t mid;
    //indices of direct children barycens in blist
    std::vector<int_t> children;

    barycen();
};


//geopotential summer
class geopotential{
    //Nz: max zonal degree
    //Nt: max tesseral degree
    // N: max(Nz,Nt)
    int_t Nz,Nt,N;
//    fast_real *J,*C,*S;
//    fast_mpvec *c_table;
    fast_mpvec c_table[];

public:
    geopotential()=delete;
    //ref_radius_factor:
    // treat the coefficients as their reference radius is R*ref_radius_factor,
    // degree-n coefficients will be multiplied by ref_radius_factor^n,
    // default to 1, geopotential is disabled when equal to 0.
    static geopotential *load(const char *file,fast_real ref_radius_factor=1);
    static void unload(geopotential *);
    //sizeof this structure
    int_t size() const;
    //R: radius of the Extended Body
    //r: point mass position from the Extended Body,
    //      represented in body-fixed coordinate system
    //      in which the harmonics are expressed
    //N_start: starting degree of harmonics
    //N_end:   ending degree of harmonics
    //      default is -1, means to sum to maximum degree
    //return acceleration of point mass per Extended Body GM,
    //      represented in body-fixed coordinate system
    fast_mpvec sum(fast_real R,fast_mpvec r,int_t N_start=3,int_t N_end=-1) const;
    INLINE fast_mpvec cuda_sum(fast_real R,fast_mpvec r,int_t N_start=3,int_t N_end=-1) const;
};

//ring attractor
class ring{
public:
    //number of disks to sum
    int_t N;
    //GM_ratio: GM of ring system / GM of host body system(include ring)
    //GM,A,J2,GL: parameters of the ring system
    fast_real GM_ratio,A,J2,GL;

    struct {
        fast_real Gs,R,H;
    } c_table[];
    ring()=delete;
    //ref_GM, ref_R2 : GM, R2 of host for normalization
    //direction_mass_factor:
    // if negative, the angular momentum of the ring is reversed,
    // mass of the ring system will be multiplied by its absolute value,
    // default to 1, ring is disabled when equal to 0.
    static ring *load(const char *file,fast_real ref_GM,fast_real ref_R2,fast_real direction_mass_factor=1);
    static void unload(ring *);
    //sizeof this structure
    int_t size() const;
    //r: point mass position from the body which hosts the ring,
    //      represented in body-fixed coordinate system
    //      in which the ring lies on xy plane
    //return acceleration of point mass per host body GM,
    //      represented in body-fixed coordinate system
    fast_mpvec sum(fast_mpvec r) const;
    INLINE fast_mpvec cuda_sum(fast_mpvec r) const;
};

//celestial object
class mass{
public:
    //evolving parameters:
    //r: position
    //v: velocity
    //GL: G * angular momentum / GM
    mpvec r,v,GL;
    //s: surface frame
    mpmat s;
    //w: angular velocity
    fast_mpvec w;

    //time-variable parameters:
    //linear calculated from time and rate
    //GM: gravitational parameter
    //exJ2: change of J2 since t=0
    fast_real GM,exJ2;

    //constant parameters:
    //sid: unique object id, actually is char[8]
    uint64_t sid;
    //GM0: GM at epoch 0.
    //dGM: dGM/dt, mass loss rate, should be <0 for stars
    //dJ2: dJ2/dt, J2 rate
    //R: mean radius
    //inertia: mean moment of inertia factor, default to 0.4
    //k2: degree 2 tidal love number, default to 0
    //k2r: degree 2 rotational love number, default to k2
    //tide_delay: tidal relaxation time scale, default to 0
    //tide_delay_factor: the tide_delay raised on others by this object
    //                   will be multiplied by this factor, default to 1
    //in https://bridge.kamine.cloud/archives/284 ,
    //{inertia, k2, k2r, tide_delay} are {2A/3,(k-1)kt,k-1,tau}, respectively
    fast_real GM0,dGM,dJ2,R,inertia,k2,k2r,tide_delay,tide_delay_factor;
    //lum: luminosity, which causes radiation pressure
    //recpt: receptance factor of radiation pressure, 1 for total absorb
    fast_real lum,recpt;
    //auxiliary constants:
    //A = 3*inertia/2
    //R2 = GM*R^2
    //rR2_4Mc = recpt*R^2/(4*GM*c/G)
    fast_real A,R2,rR2_4Mc;
    //static deformation of potential in surface frame
    fast_mpmat C_static;

    //geopoteintial model:
    geopotential *gpmodel;
    //ring model:
    ring *ringmodel;

    //auxiliary variables:
    //phi: newtonian gravitational potential
    //beta2: (v/c)^2
    fast_real phi,beta2;
    //naccel: point mass newtonian acceleration
    //gaccel: point mass gravitoelectric acceleration
    //daccel: acceleration due to extended body deformation and lense thirring effect
    //dtorque: G * torque / GM due to extended body deformation and lense thirring effect
    //beta: v/c
    fast_mpvec naccel,daccel,dtorque,gaccel,beta;
    //GI: G * moment of inertia / GM
    //C_potential = k2 * C + C_static
    fast_mpmat GI,C_potential;

    //temporary variables used by integrators
    fast_real Erot;
    fast_mpvec Egrad;
    fast_mpvec idaccel,idtorque;
    //used by msystem::integrate to test if a collision occurred and whether time step is appropriate
    // and by msystem::combined_integrate to test if a capture occurred
    fast_real min_distance;
    fast_real max_influence;

    //calculate deformation matrix(C_potential) and inertia matrix(GI)
    //calculate Newtonian acceleration(naccel) & potential(phi)
    void deform_this(const std::vector<mass> &mlist);
    //calculate deformation matrices(C_potential) and inertia matrices(GI)
    //calculate Newtonian acceleration(naccel) & potential(phi)
    static void deform_all(std::vector<mass> &mlist);
    //resize the mass by a factor
    //updates all relevant parameters
    void scale(fast_real factor);
};

//short mass used in Runge-Kutta-integrators
struct mass_state{
    //evolving parameters:
    //r: position
    //v: velocity
    //GL: G * angular momentum / GM
    mpvec r,v,GL;
    //s: surface frame
    mpmat s;
    //w: angular velocity
    fast_mpvec w;

    //auxiliary variables:
    //naccel: point mass acceleration
    //dtorque: G * torque / GM
    fast_mpvec naccel,dtorque;

    INLINE mass_state &operator =(const mass &src){
        r=src.r;
        v=src.v;
        GL=src.GL;
        s=src.s;
        w=src.w;
        /*naccel=0;
        dtorque=0;
        Prot=0;*/
        return *this;
    }
};

//stellar system
class msystem{
public:
    //relativistic coordinate time
    real t_eph;
    //frame time, dt for integration
    fast_real delta_t;
    
    //  0: use CPU RungeKutta
    //  1: use GPU RungeKutta
    //  2: use CPU/GPU Combined RungeKutta
    enum {
        CPU_RK12        =0,
        GPU_RK12        =1,
        COMBINED_RK12   =2
    };
    int_t integrator;
    
    //cadence of data_output
    fast_real data_cadence;

    //max length of a single output file, in seconds
    fast_real max_ephm_length;

    // if integrator==COMBINED_RK12, following should be specified
    fast_real combined_delta_t;
    fast_real GM_max_child,GM_max_parent,GM_max_tiny,Period_max_child;

    //list of barycens
    std::vector<barycen> blist;
    //list of masses
    std::vector<mass> mlist;
    //index of masses
    std::map<uint64_t,int_t> midx;

    //tidal corrections
    int_t tidal_parent;
    std::vector<int_t> tidal_childlist;
    fast_mpmat tidal_matrix;

    //update time-variables of mass list to epoch t
    void update(fast_real t);
    //calculate deformation matrices(C_potential) and inertia matrices(GI)
    //calculate Newtonian acceleration(naccel) & potential(phi)
    void deform();
    //calculate acceleration and solve for angular velocity
    void accel();
    //same as accel, use GPU
    void Cuda_accel();

    //Runge-Kutta-12 integrator
    void RungeKutta12(fast_real dt,int_t n_step);
    //Runge-Kutta-12 integrator
    void Cuda_RungeKutta12(fast_real dt,int_t n_step);

    //integrate ephemerides
    //USE_GPU:   0: CPU Runge Kutta 12
    //           1: GPU Runge Kutta 12
    void integrate(fast_real dt,int_t n_step,int USE_GPU=0);
    //integrate ephemerides
    //for full system with un-parented tiny masses, use (dt*n_combine) as time step
    //for subsystem with tiny children masses, use dt as time step, always use CPU
    //USE_GPU:   use CPU(0)/GPU(1) for single-step integration of full system
    void combined_integrate(fast_real dt,int_t n_combine,int_t n_step,int USE_GPU=1);

    //analyse position of masses to build barycen list
    //reconstruct: if true, analyse from scratch;
    //             if false, update existing structure when necessary, this is faster.
    bool analyse(bool reconstruct=false);
    //update barycens' GM & rv
    void update_barycens();

    //load system from directory
    //     dir : directory
    // fconfig : config file name
    bool load_dir(std::map<std::string,std::string> &config,const char *dir,const char *fconfig);

    //load system from files
    //   fbase : basic parameters and initial states
    //    fext : extra parameters
    //  gppath : path for geopotential models
    //ringpath : path for ring models
    bool load(const char *fbase,
              const char *fext=nullptr,
              const char *gppath=nullptr,
              const char *ringpath=nullptr);

    //load system from checkpoint file
    bool load_checkpoint(class mem_file *fcp);
    //save system as checkpoint
    bool save_checkpoint(class mem_file *fcp);

    //get index of mass from sid, return -1 if not found
    int_t get_mid(const char *sid);
    int_t get_mid(uint64_t sid);
};
