#pragma once
#include<mutex>
#include<set>
#include<map>
#include"physics/mass.h"
#include"utils/memio.h"
#include"math/Keplerians.h"

class ephemeris_compressor{
public:
    //compressed data format
    enum Format:int_t{
        UNKNOWN             =0,
        //orbital                  type[channel], data components , error function
        STATE_VECTORS       =1,     // vec[1]   , use x,y,z,      , absolute_state_error
        KEPLERIAN_VECTORS   =2,     // vec[1]   , use x,y,z,      , relative_state_error
        KEPLERIAN_CIRCULAR  =3,     // double[6], use j, ex,ey, ml, circular_kepler_error
        KEPLERIAN_RAW       =4,     // double[6], use j, q, el, m , raw_kepler_error
        //rotational
        AXIAL_OFFSET        =5,     // double[ ], use w, ???      , axial_rotation_error
        QUATERNION          =6,     // quat[1]  , use q(s)        , quaterion_rotation_error
    };

    class keplerian:public ephem_orb{
    public:
        //longitude of ascending node + earg
        double el;
        //eccentricity vector
        //e * sincos(el)
        double ex,ey;
        //mean longitude
        //el + mean anomaly
        double ml;

        keplerian()=default;
        keplerian(const ephem_orb &other);
        keplerian(const double *k,bool circular);

        // test if two states are interpolatable
        static bool interpolatable(bool circular,const keplerian &k1,const keplerian &k2);

        //blend multiple keplerians by weights
        void blend_initialize(bool circular);
        void blend_add(bool circular,const keplerian &component,double weight);
        void blend_finalize(bool circular);

    };

    struct orbital_state_t{
        vec r,v;
    };
    struct rotational_state_t{
        vec w,x,z;
    };
    struct data_header{
        uint64_t uformat;  // = ~uint64_t(enum Format)
        int_t degree;  // of bspline basis
        int_t n;       // # of segments. total size of data is (n+degree)*channel*sizeof(type).
        double relative_error;
    };
private:
    template<typename T,size_t N_Channel,int_t format>
    static MFILE compress_data(const T *pdata,int_t N,int_t d);

    static double infer_GM_from_data(const orbital_state_t *pdata,int_t N);

    static double relative_state_error(const vec *r,const vec *rp);
    static double absolute_state_error(const vec *r,const vec *rp);
    static double circular_kepler_error(const double *k,const double *kp);
    static double raw_kepler_error(const double *k,const double *kp);
    static double axial_rotation_error(const double *a,const double *ap);
    static double quaterion_rotation_error(const quat *q,const quat *qp);

    //minimizing this will achieve a balance between relative error and file size
    static double compression_score(double relative_error,double compressed_size);
public:
    //max possible degree of bspline fitting, must be odd
    //11 is maximum odd number not exceed the degree of RungeKutta integrator
    //do not change this
    static constexpr int_t max_bspline_degree=11;
    //criterion below which error is thought mainly due to fp-imprecision, not model.
    static constexpr double epsilon_relative_error=1e-12;
    //same, but for absolute positional error(m)
    static constexpr double epsilon_absolute_error=1e-4;

    // mf: contains raw orbital_state_t data
    // dt: time cadence between data points
    static bool compress_orbital_data(MFILE &mf,double dt);
    // mf: contains raw rotational_state_t data
    // dt: time cadence between data points
    static bool compress_rotational_data(MFILE &mf,double dt);
};

class ephemeris_collector{
    struct orbital_datapack_t{
        int_t tid;
        int_t t_start;
        int_t t_end;
        int_t parent_barycen_id;
        MFILE data;
    };
    struct rotational_datapack_t{
        int_t t_start;
        int_t t_end;
        MFILE data;
    };

    struct index_entry_t{
        //+: orbital data file, from 1...
        //-: rotational data file, from -1...
        //0: barycen structure
        int_t fid;
        //for data file: sid of related mass
        //for barycen structure: vector<barycen>.size
        uint64_t sid;
        int_t t_start;
        int_t t_end;
    };

private:
    msystem &ms;
    std::vector<barycen> blist;

    // { { mids of parent barycen, mids of child barycen }, index of pair }
    std::map<std::pair<std::set<int_t>,std::set<int_t>>,int_t> barycen_ids;
    std::vector<orbital_datapack_t> orbital_data;
    std::vector<rotational_datapack_t> rotational_data;
    int_t t_start;
public:
    ephemeris_collector(msystem &_ms);

    //update state vectors of blist by ms.mlist
    void update_barycens();
    const std::vector<barycen> &get_barycens() const{ return blist; }

    //convert state vectors in blist as relative to direct parent
    int_t decompose(int_t bid=-1);
    //restore state vectors in blist to absolute
    int_t compose(int_t bid=-1);

    //record state vectors
    void record();

    //rebind ephemeris data with parent barycen
    //update datapacks::tid & pbarycen
    void rebind();

    void extract(std::vector<MFILE> &ephm_files,bool force);

    //convert zips to old data pack
    static int convert_format(const char *path);
};

class ephemeris_generator{
    static std::mutex io_mutex;
public:
    const char *ip;
    const char *op;

    double t_years;
    //  1: only do forward integration
    // -1: only do backward integration
    //  0: do both, default
    int fix_dir=0;

    int make_ephemeris(int dir);
};
