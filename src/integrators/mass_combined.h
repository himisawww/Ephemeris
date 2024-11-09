#pragma once
#include"physics/mass.h"
#include"utils/memio.h"

class msystem_combinator{
    msystem *pms;
    fast_real t_substep;
    real t_split;

    std::map<uint64_t,std::vector<barycen>> sublists;
    std::map<uint64_t,MFILE> orbital_subdata;
    //std::map<uint64_t,MFILE> rotational_subdata;

    friend class msystem;
    friend class ephemeris_collector;
public:
    msystem_combinator():pms(nullptr){}

    //build blist for a subsystem mssub of *pms
    //  s.t. for tidal_childlist of mssub,
    //  decomposed states in mssub is the same as in *pms
    //  i.e. extract a sub-blist rooted from tidal_parent
    int_t link(msystem &mssub,int_t tid=-1);
};
