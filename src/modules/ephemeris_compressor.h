#pragma once
#include"definitions.h"
#include"utils/memio.h"

struct ephemeris_entry;

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
    static double quaternion_rotation_error(const quat *q,const quat *qp);

    //minimizing this will achieve a balance between relative error and file size
    static double compression_score(double relative_error,double compressed_size);
public:
    //compressed data format
    enum format_t:uint16_t{
        NONE                =0,     // uncompressed orbital_state_t/rotational_state_t
        //orbital                  data_type[data_channel], data components , error function
        STATE_VECTORS       =1,     // vec[1]   , use x,y,z,      , absolute_state_error
        KEPLERIAN_VECTORS   =2,     // vec[1]   , use x,y,z,      , relative_state_error
        KEPLERIAN_CIRCULAR  =3,     // double[6], use j, ex,ey, ml, circular_kepler_error
        KEPLERIAN_RAW       =4,     // double[6], use j, q, el, m , raw_kepler_error
        //rotational
        AXIAL_OFFSET        =5,     // double[6], use w, py,pz,ang, axial_rotation_error
        QUATERNION          =6,     // quat[1]  , use q(s)        , quaternion_rotation_error
        TIDAL_LOCK          =7,     // quat[1]  , use q(s in -r0j), quaternion_rotation_error
        UNKNOWN               ,     // last enum
        INVALID             =255
    };
    static constexpr bool is_valid_format(format_t f){ return NONE<f&&f<UNKNOWN; }
    static constexpr bool is_orbital_format(format_t f){ return STATE_VECTORS<=f&&f<=KEPLERIAN_RAW; }
    static constexpr bool is_rotational_format(format_t f){ return AXIAL_OFFSET<=f&&f<=TIDAL_LOCK; }
    static constexpr const char *format_name(format_t f){
        switch(f){
#define FORMAT_NAME_ENTRY(F) case F:return #F
            FORMAT_NAME_ENTRY(NONE);

            FORMAT_NAME_ENTRY(STATE_VECTORS);
            FORMAT_NAME_ENTRY(KEPLERIAN_VECTORS);
            FORMAT_NAME_ENTRY(KEPLERIAN_CIRCULAR);
            FORMAT_NAME_ENTRY(KEPLERIAN_RAW);

            FORMAT_NAME_ENTRY(AXIAL_OFFSET);
            FORMAT_NAME_ENTRY(QUATERNION);
            FORMAT_NAME_ENTRY(TIDAL_LOCK);

            FORMAT_NAME_ENTRY(INVALID);
#undef FORMAT_NAME_ENTRY
        }
        return "UNKNOWN";
    }

    struct header_base{
        float relative_error;
        int16_t degree;     // of bspline basis
        uint16_t uformat;   // = ~uint16_t(enum format_t)
                            // here most 12 bits of uformat is 0xfff for valid compression methods.
                            // i.e. first 8 bytes of file is double(NAN) if data is compressed.
                            //      otherwise data was left unchanged as a raw sequence of double.
        int_t n;            // # of segments. total size of data is (n+degree)*data_channel*sizeof(data_type).
        float smooth_range[2];  //{left, right}, in (0,1], in unit of segments
    };
    template<format_t F>struct header_t;
    template<>struct header_t<STATE_VECTORS>:public header_base{
        typedef vec data_type;
        static constexpr size_t data_channel=1;
        static constexpr auto error_function=absolute_state_error;
        typedef orbital_state_t state_type;
        
        //{ left_end_data - fit, derivative, right_end_data - fit, derivative }
        data_type fix[4][data_channel];
    };
    template<>struct header_t<KEPLERIAN_VECTORS>:public header_base{
        typedef vec data_type;
        static constexpr size_t data_channel=1;
        static constexpr auto error_function=relative_state_error;
        typedef orbital_state_t state_type;

        data_type fix[4][data_channel];
        double GM;
    };
    template<>struct header_t<KEPLERIAN_CIRCULAR>:public header_base{
        typedef double data_type;
        static constexpr size_t data_channel=6;
        static constexpr auto error_function=circular_kepler_error;
        typedef orbital_state_t state_type;

        data_type fix[4][data_channel];
        double GM;
    };
    template<>struct header_t<KEPLERIAN_RAW>:public header_base{
        typedef double data_type;
        static constexpr size_t data_channel=6;
        static constexpr auto error_function=raw_kepler_error;
        typedef orbital_state_t state_type;

        data_type fix[4][data_channel];
        double GM;
    };
    template<>struct header_t<AXIAL_OFFSET>:public header_base{
        typedef double data_type;
        static constexpr size_t data_channel=6;
        static constexpr auto error_function=axial_rotation_error;
        typedef rotational_state_t state_type;

        data_type fix[4][data_channel];
        double local_axis_theta,local_axis_phi;
    };
    template<>struct header_t<QUATERNION>:public header_base{
        typedef quat data_type;
        static constexpr size_t data_channel=1;
        static constexpr auto error_function=quaternion_rotation_error;
        typedef rotational_state_t state_type;

        data_type fix[4][data_channel];
    };
    template<>struct header_t<TIDAL_LOCK>:public header_base{
        typedef rotational_state_t state_type;
        typedef quat data_type;
        static constexpr size_t data_channel=1;
        static constexpr auto error_function=quaternion_rotation_error;

        data_type fix[4][data_channel];
    };

    class interpolator:private header_base{
        typedef header_base base_t;

        //auxiliary data
        char fix_data[192];
        union{
            double sGM;
            mat local_axis;
            orbital_state_t orb;
        };

        void *pfitter;
        double t_range;

        template<format_t F>
        //convert data to state
        void convert(
            typename header_t<F>::state_type *state,
            const typename header_t<F>::data_type *x,
            const typename header_t<F>::data_type *v) const;

        template<typename T,size_t N_Channel>
        //assume 0<=x<smooth_range
        void smooth(T (&r)[N_Channel],T (&v)[N_Channel],
            const T (&fix_r)[N_Channel],const T (&fix_v)[N_Channel],
            double x,double smooth_range) const;
    public:
        interpolator(MFILE *fin,double _range);
        interpolator(interpolator &&);
        ~interpolator();

        void expand();

        format_t data_format() const{ return format_t(~base_t::uformat); }
        double relative_error() const{ return base_t::relative_error; }
        operator bool() const{ return pfitter; }
        bool is_orbital() const{ return is_orbital_format(data_format()); }
        bool is_rotational() const{ return is_rotational_format(data_format()); }

        // if data_format()==TIDAL_LOCK, orbital state at the same instant
        // should be set by this function before interpolating rotational states.
        void set_orbital_state(const vec &r,const vec &v);
        // state = &orbital_state_t or &rotational_state_t
        void operator()(double t,void *state) const;
    };
private:
    template<format_t F,typename T=typename header_t<F>::data_type,size_t N_Channel=header_t<F>::data_channel>
    static MFILE compress_data(const T *pdata,int_t N,int_t d,const T fix_derivative[2][N_Channel]=nullptr);
    static int_t select_best(MFILE &mf,std::vector<MFILE> &compressed_results,int_t N,int_t d);
    // for number of datapoints and degree of bspline fitting,
    // fill possible choices of segments of compressed bsplines in decreasing order.
    static void segment_choices(std::vector<int_t> &result,int_t N,int_t d);

    struct compress_work{
        MFILE *morb,*mrot,*msuborb,*msubrot;
        const ephemeris_entry *pindex;

        // compress informations
        header_base *pheaders[2];
        int_t clevels[2];
        size_t newsize[2],oldsize[2];
        double max_r,max_v,max_xz,max_w;
        double end_r,end_v,end_xz,end_w;
        void run();
        int_t priority() const;
    };
    static void do_compress_work(void *pcw,size_t thread_id){ ((compress_work*)pcw)->run(); }
public:
    //max possible degree of bspline fitting, must be odd
    //11 is maximum odd number not exceed the degree of RungeKutta integrator
    //do not change this
    static constexpr int_t max_bspline_degree=11;
    //criterion below which error is thought mainly due to fp-imprecision, not model.
    static constexpr double epsilon_relative_error=1e-12;
    //same, but for absolute positional error(m)
    static constexpr double epsilon_absolute_error=1e-4;
    //emit warning if relative fit error not less than this
    static constexpr double relative_error_warning_threshold=1e-5;

    // mf: contains raw orbital_state_t data
    // time_span: time between first & last data point, i.e. delta_t*(N-1)
    // return level of compression, 0 means failed and mf is untouched.
    // if successed, 0 <= relative_error < 1
    static int_t compress_orbital_data(MFILE &mf,double time_span);
    // mf: contains raw rotational_state_t data
    // if(morb) also try to use TIDAL_LOCK method with orbital_data
    // others same as compress_orbital_data
    static int_t compress_rotational_data(MFILE &mf,double time_span,MFILE *morb=nullptr);

    // ephemeris_data: full datapack produced by ephemeris_generator and collected by ephemeris_collector,
    //      see ephemeris_generator::make_ephemeris
    // return -1 if critical failure occurs, input is untouched in that case.
    // otherwise return number of failures (entries that cannot be compressed).
    static int_t compress(std::vector<MFILE> &ephemeris_data);
};

typedef ephemeris_compressor::format_t ephemeris_format;
typedef ephemeris_compressor::interpolator ephemeris_interpolator;
