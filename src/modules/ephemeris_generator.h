#pragma once
#include<mutex>
#include<set>
#include<map>
#include<string>
#include"physics/mass.h"
#include"utils/memio.h"
#include"integrators/mass_combined.h"

class ephemeris_collector{
    struct datapack_t{
        int_t tid;
        int_t t_start;
        int_t t_end;
        int_t parent_barycen_id;
        MFILE orbital_data;
        MFILE rotational_data;
    };

    struct index_entry_t{
        //0: barycen structure
        //+: orbital&rotational data file, from 1...
        int_t fid;
        //for data file: sid of related mass
        //for barycen structure: vector<barycen>.size
        uint64_t sid;
        int_t t_start;
        int_t t_end;

        std::string entry_name(bool rotational,bool substep);
    };

private:
    msystem &ms;
    std::vector<barycen> blist;
    //t_eph when blist is bind with ms
    real t_bind;

    // { { mids of parent barycen, mids of child barycen }, index of pair }
    std::map<std::pair<std::set<int_t>,std::set<int_t>>,int_t> barycen_ids;
    std::vector<datapack_t> data;
    int_t t_start;

    std::map<uint64_t,int_t> file_index;
public:
    ephemeris_substeper m_substeper;

    ephemeris_collector(msystem &_ms);

    const std::vector<barycen> &get_barycens() const{ return blist; }

    //record state vectors
    void record();

    //return if structure of blist is up to date
    bool synchronized();
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
