#pragma once
#include"math/Keplerians.h"
#include"utils/memio.h"

class ephemeris_compressor{
public:
    struct orbital_state_t{
        vec r,v;
    };
    struct rotational_state_t{
        vec w,x,z;
    };
private:
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
    //compressed data format
    enum Format:int_t{
        NONE                =0,
        //orbital                  type[channel], data components , error function
        STATE_VECTORS       =1,     // vec[1]   , use x,y,z,      , absolute_state_error
        KEPLERIAN_VECTORS   =2,     // vec[1]   , use x,y,z,      , relative_state_error
        KEPLERIAN_CIRCULAR  =3,     // double[6], use j, ex,ey, ml, circular_kepler_error
        KEPLERIAN_RAW       =4,     // double[6], use j, q, el, m , raw_kepler_error
        //rotational
        AXIAL_OFFSET        =5,     // double[6], use w, py,pz,ang, axial_rotation_error
        QUATERNION          =6,     // quat[1]  , use q(s)        , quaterion_rotation_error
        TIDAL_LOCK          =7,     // quat[1]  , use q(s in -r0j), quaterion_rotation_error
        UNKNOWN                     // last enum
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

    class axial:public ephem_rot{
    public:
        axial()=default;
        axial(const ephem_rot &other):ephem_rot(other){}
        axial(const double *a);
    };

    struct header_base{
        float relative_error;
        int16_t degree;     // of bspline basis
        uint16_t uformat;   // = ~uint16_t(enum Format)
                            // here most 12 bits of uformat is 0xfff for valid compression methods.
                            // i.e. first 8 bytes of file is double(NAN) if data is compressed.
                            //      otherwise data was left unchanged as a raw sequence of double.
        int_t n;            // # of segments. total size of data is (n+degree)*channel*sizeof(type).
    };
    template<Format format>struct header_t;
    template<>struct header_t<STATE_VECTORS>:public header_base{
        typedef vec type;
        static constexpr size_t channel=1;
        static constexpr auto error_function=absolute_state_error;

        vec r0,v0,r1,v1;
    };
    template<>struct header_t<KEPLERIAN_VECTORS>:public header_base{
        typedef vec type;
        static constexpr size_t channel=1;
        static constexpr auto error_function=relative_state_error;

        vec r0,v0,r1,v1;
    };
    template<>struct header_t<KEPLERIAN_CIRCULAR>:public header_base{
        typedef double type;
        static constexpr size_t channel=6;
        static constexpr auto error_function=circular_kepler_error;

        typedef double kep_ct[6];
        kep_ct k0,k1;
    };
    template<>struct header_t<KEPLERIAN_RAW>:public header_base{
        typedef double type;
        static constexpr size_t channel=6;
        static constexpr auto error_function=raw_kepler_error;

        typedef double kep_rt[6];
        kep_rt k0,k1;
    };
    template<>struct header_t<AXIAL_OFFSET>:public header_base{
        typedef double type;
        static constexpr size_t channel=6;
        static constexpr auto error_function=axial_rotation_error;

        double local_axis_theta,local_axis_phi;
        quat q0;
        vec w0;
        quat q1;
        vec w1;
    };
    template<>struct header_t<QUATERNION>:public header_base{
        typedef quat type;
        static constexpr size_t channel=1;
        static constexpr auto error_function=quaterion_rotation_error;

        quat q0;
        vec w0;
        quat q1;
        vec w1;
    };
    template<>struct header_t<TIDAL_LOCK>:public header_base{
        typedef quat type;
        static constexpr size_t channel=1;
        static constexpr auto error_function=quaterion_rotation_error;

        quat q0;
        vec w0;
        quat q1;
        vec w1;
    };
private:
    template<Format format,typename T=typename header_t<format>::type,size_t N_Channel=header_t<format>::channel>
    static MFILE compress_data(const T *pdata,int_t N,int_t d);
    static void select_best(MFILE &mf,std::vector<MFILE> &compressed_results);
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
    // ephemeris_data: full datapack produced by ephemeris_generator and collected by ephemeris_collector,
    //      see ephemeris_generator::make_ephemeris
    // return number of failures
    static int_t compress(std::vector<MFILE> &ephemeris_data);
};
