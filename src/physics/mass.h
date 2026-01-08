#pragma once
#include"definitions.h"

#include<cstdint>
#include<string>
#include<vector>
#include<map>

class msystem;

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

class bsystem:public std::vector<barycen>{
    int_t decompose(int_t bid);
    int_t compose(int_t bid);
public:
    //fill tid and mid
    void fill_tid();
    //update state vectors of blist by ms.mlist
    void update_barycens(const msystem &ms);
    //convert state vectors in blist as relative to direct parent
    int_t decompose(){ return decompose(root_id()); }
    //restore state vectors in blist to absolute
    int_t compose(){ return compose(root_id()); }
    //get index of root, -1 if non-exist;
    int_t root_id() const;

    bool load_barycen_structure(MFILE *fin,size_t bsize);
    void save_barycen_structure(MFILE *fout) const;
    //return number of mass in this bsystem, or 0 if not self-compatible
    size_t compatible_size() const;

    void print_structure(MFILE *mf,const msystem &) const;
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
    //R: volumetric mean radius, also reference radius for geopotential/ring
    //   if negative, use soft gravity(prop. to r when r < |R|)
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
    //R2 = R^2
    //sR2 = R^2*(R<0), present for soft mass
    //rR2G_4c = recpt*R^2*G/(4*c)
    fast_real A,R2,sR2,rR2G_4c;
    //static deformation of potential in surface frame
    fast_mpmat C_static;

    //geopotential model:
    const geopotential *gpmodel;
    //ring model:
    const ring *ringmodel;

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
    fast_mpvec Egrad;
    fast_mpvec idaccel,idtorque;
    fast_mpvec jdaccel,jdtorque;
    int_t task_index,task_jndex,task_count;
    void *pmlist,*qmlist;

    //used by msystem::integrate to test if a collision occurred and whether time step is appropriate
    // and by msystem::combined_integrate to test if a capture occurred
    fast_real min_distance;
    fast_real max_influence;

    const char *get_ssid() const{ return (const char*)&sid; }

    INLINE void orthogonalize(){
        fast_mpmat fmis(s);
        s+=0.5*fmis%(1-fmis.transpose()%fmis);
    }

    //calculate deformation matrix(C_potential) and inertia matrix(GI)
    //calculate Newtonian acceleration(naccel) & potential(phi)
    void deform_by(const std::vector<mass> &mlist);
    //resize the radius by a factor
    //updates all relevant parameters
    void scale(fast_real factor);
    //check sanity of states & params
    bool sanity(bool alert=false) const;
};

constexpr auto mass_constant_parameters_begin=&mass::sid;
constexpr auto mass_constant_parameters_end=&mass::gpmodel;
constexpr auto mass_auxiliary_variables=&mass::phi;
constexpr auto mass_temporary_variables=&mass::Egrad;

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

class ephemeris_generator;
class ephemeris_substeper;
class ephemeris_reader;

//stellar system
class msystem{
public:
    //  0: use CPU RungeKutta
    //  1: use GPU RungeKutta
    //  2: use CPU/GPU Combined RungeKutta
    enum {
        CPU_RK12        =   0,
        GPU_RK12        =   1,
        COMBINED_RK12   =   2
    };

    friend class ephemeris_generator;
    friend class ephemeris_substeper;
    friend class ephemeris_reader;
private:
    //relativistic coordinate time
    real t_eph;
    //frame time, dt for integration
    fast_real delta_t;

    int_t integrator;

    //cadence of data_output
    fast_real data_cadence;

    //max length of a single output file, in seconds
    fast_real max_ephm_length;

    // if integrator==COMBINED_RK12, following should be specified
    fast_real combined_delta_t;
    fast_real GM_max_child,GM_max_parent,GM_max_tiny,Period_max_child;

    //tidal corrections
    int_t tidal_parent;
    fast_mpmat tidal_matrix;
    //above two is undefined if tidal_childlist.empty()
    std::vector<int_t> tidal_childlist;
    ephemeris_substeper *p_substeper;

    //t_eph when blist is analysed
    real t_barycen;
    //t_eph when analyse updates blist
    real t_update;
    //list of barycens
    bsystem blist;
    //list of masses
    std::vector<mass> mlist;
    //index of masses
    std::map<uint64_t,int_t> midx;

    //loaded components used by masses, will be freed on destruction
    std::vector<const geopotential *> gp_components;
    std::vector<const ring *> ring_components;
private:
    //load system from files
    //   fbase : basic parameters and initial states
    //    fext : extra parameters
    //  gppath : path for geopotential models
    //ringpath : path for ring models
    bool load(const char *fbase,
        const char *fext,
        const char *gppath,
        const char *ringpath);
    //Runge-Kutta-12 integrator
    void RungeKutta12(fast_real dt,int_t n_step);
    //Runge-Kutta-12 integrator
    void Cuda_RungeKutta12(fast_real dt,int_t n_step);
    //same as accel, use GPU
    void Cuda_accel();
    //calculate deformation matrices(C_potential) and inertia matrices(GI)
    //calculate Newtonian acceleration(naccel) & potential(phi)
    void deform();
    //used by combined_integrate with p_substeper
    void record_substeps(fast_real dt,bool initialize=false);

public:
    //update time-variables of mass list to epoch t
    //if pblist present, also update GM/GM_sys in *pblist
    void update(fast_real t,bsystem *pblist=nullptr);
    //calculate acceleration and solve for angular velocity
    //parallel_option:
    //-1: CPU single thread;
    // 0: CPU threadpool when available
    // 1: GPU
    void accel(int parallel_option=0);
    //reset all the parameters calculated by accel() to a initial state
    void clear_accel();

    //integrate ephemerides
    //USE_GPU:   0: CPU Runge Kutta 12
    //           1: GPU Runge Kutta 12
    void integrate(fast_real dt,int_t n_step,int USE_GPU=0);
    //integrate ephemerides
    //for full system with un-parented tiny masses, use (dt*n_combine) as time step
    //for subsystem with tiny children masses, use dt as time step, always use CPU
    //USE_GPU:   use CPU(0)/GPU(1) for combined integration of full system
    //ps: if !nullptr, collects substeps of children in subsystems
    void combined_integrate(fast_real dt,int_t n_combine,int_t n_step,int USE_GPU=1,
        ephemeris_substeper *ps=nullptr);

    //analyse position of masses to build barycen list
    //reconstruct: if true, analyse from scratch;
    //             if false, update existing structure when necessary, this is faster.
    //return latest t_eph when analyse updates blist
    real analyse(bool reconstruct=false);

    const bsystem &get_barycens() const{ return blist; }

    //load internal system, optional save to a checkpoint
    bool load_internal(const char *fcheckpoint=nullptr);
    //load system from config file, optional save to a checkpoint
    bool load(const char *fconfig,const char *fcheckpoint=nullptr);
    //load system from checkpoint file
    bool load_checkpoint(MFILE *fcp);
    //save system as checkpoint
    bool save_checkpoint(MFILE *fcp) const;

    void build_mid();
    //get index of mass from sid, return -1 if not found
    int_t get_mid(const char *ssid) const;
    int_t get_mid(uint64_t sid) const;

    void reset_params();
    void copy_params(const msystem &);
    void clear();
    bool push_back(const mass &m);

    msystem &operator =(const msystem &other);
    //compare time-constant parameters of masses and msystem
    //that is, msystem is considered same after integrate() for any amount of time
    bool is_same(const msystem &other);

    msystem(){}
    msystem(msystem &&)=default;
    msystem(const msystem &other){ *this=other; }
    ~msystem(){ clear(); }

    auto begin() const{ return mlist.begin(); }
    auto end() const{ return mlist.end(); }
    auto size() const{ return mlist.size(); }
    auto empty() const{ return mlist.empty(); }
    const mass &operator [](int_t mid) const{ return mlist[mid]; }
    const mass &operator [](const char *ssid) const{ return mlist[get_mid(ssid)]; }

    //functions modifying mass of the system; once called,
    //  accel() must be called to restore correct internal state of msystem.
    void scale(int_t mid,fast_real factor);
    void scale_geopotential(int_t mid,fast_real factor);
    void scale_ring(int_t mid,fast_real factor);
    auto begin(){ return mlist.begin(); }
    auto end(){ return mlist.end(); }
    mass &operator [](int_t mid){ return mlist[mid]; }
    mass &operator [](const char *ssid){ return mlist[get_mid(ssid)]; }

    // t_eph should always be integers
    int_t ephemeris_time() const;
};
