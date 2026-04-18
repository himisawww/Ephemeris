#pragma once
#include<mutex>
#include"htl/set.h"
#include"htl/map.h"
#include<string>
#include"physics/mass.h"
#include"utils/memio.h"

struct ephemeris_entry{
    //0: barycen structure
    //+: orbital&rotational data file, from 1...
    //for orbital barycentric offset cache file: +index from 1...
    int_t fid;
    //for data file: sid of related mass
    //for barycen structure: vector<barycen>.size
    //for orbital barycentric offset cache file: 0
    uint64_t sid;
    int_t t_start;
    int_t t_end;

    std::string entry_name(bool rotational,bool substep) const;
    std::string offset_name(bool substep) const;
};

class ephemeris_collector{
    struct datapack_t{
        int_t tid;
        int_t t_start;
        int_t t_end;
        int_t parent_barycen_id;
        MFILE orbital_data;
        MFILE rotational_data;
    };
    struct cachepack_t{
        int_t t_start;
        int_t t_end;
        MFILE offset_data;
        MFILE offset_subdata;
    };
    msystem &ms;
    bsystem blist;
    //t_eph when blist is bind with ms
    real t_bind;

    // { { mids of parent barycen, mids of child barycen }, index of pair }
    htl::map<std::pair<htl::set<int_t>,htl::set<int_t>>,int_t> barycen_ids,cache_ids;
    htl::vector<datapack_t> data;
    htl::vector<cachepack_t> cache;
    // [i-th barycen in ephm_index file]={{bid in ms::blist, cid in cache}...}
    htl::vector<htl::map<int_t,int_t>> bcache_maps;
    int_t t_start;

    htl::map<uint64_t,int_t> file_index;

    struct subdatapack_t{
        MFILE orbital_data;
        MFILE rotational_data;
    };
    struct subsystem_t{
        bsystem sublist;
        //{ bid in mssub, cid in cache}
        htl::map<int64_t,int64_t> subcache_map;
    };

    fast_real t_substep;
    real t_link;

    htl::map<uint64_t,subsystem_t> sublists;
    htl::map<uint64_t,subdatapack_t> subdata;

    friend class msystem;
    // returns corresponding bid in subsys::sublist for bid in ms::blist
    int_t link(subsystem_t &subsys,msystem &mssub,int_t bid);
public:

    //build blist for a subsystem mssub of c::ms
    //  s.t. for tidal_childlist of mssub,
    //  decomposed states in mssub is the same as in c::ms
    //  i.e. extract a sub-blist rooted from tidal_parent
    int_t link(msystem &mssub);

    ephemeris_collector(msystem &_ms);

    //record state vectors
    void record();

    //return if structure of blist is up to date
    bool synchronized();
    //rebind ephemeris data with parent barycen
    //update datapacks::tid & pbarycen
    void rebind();

    //return error message
    std::string extract(htl::vector<MFILE> &ephm_files,bool force);

    //convert zips to old data pack
    //for compressed format only:
    //  fix_interval: set output data interval
    //   psid_subset: if present, only output objects in *psid_subset
    static int convert_format(const char *path,int_t fix_interval=0,htl::vector<const char*> *psid_subset=nullptr);
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
