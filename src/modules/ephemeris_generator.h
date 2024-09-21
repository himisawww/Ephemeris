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
    enum class Format{
        //orbital
        STATE_VECTORS,
        KEPLERIAN_CIRCULAR, //use j, ex,ey, ml
        KEPLERIAN_RAW,      //use j, q, el, m
        //rotational
        POLAR_OFFSET,
        QUATERNION
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
private:
    static double infer_GM_from_data(orbital_state_t *pdata,int_t N);
public:
    //max possible degree of bspline fitting, must be odd
    //11 is maximum odd number not exceed the degree of RungeKutta integrator
    //do not change this
    static constexpr int_t max_bspline_degree=11;

    // mf: contains raw orbital_state_t data
    // dt: time cadence between data points
    static void compress_orbital_data(MFILE &mf,double dt);
    // mf: contains raw rotational_state_t data
    // dt: time cadence between data points
    static void compress_rotational_data(MFILE &mf,double dt);
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
