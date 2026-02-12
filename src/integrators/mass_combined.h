#pragma once
#include"physics/mass.h"
#include"utils/memio.h"

class ephemeris_substeper{
    struct subdatapack_t{
        MFILE orbital_data;
        MFILE rotational_data;
    };


    msystem *pms;
    fast_real t_substep;
    real t_link;

    std::map<uint64_t,bsystem> sublists;
    std::map<uint64_t,subdatapack_t> subdata;

    friend class msystem;
    friend class ephemeris_collector;
    int_t link(msystem &mssub,int_t tid);
public:
    ephemeris_substeper():pms(nullptr){}

    //build blist for a subsystem mssub of *pms
    //  s.t. for tidal_childlist of mssub,
    //  decomposed states in mssub is the same as in *pms
    //  i.e. extract a sub-blist rooted from tidal_parent
    int_t link(msystem &mssub){ return link(mssub,-1); }
};
