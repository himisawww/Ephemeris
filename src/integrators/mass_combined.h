#pragma once
#include"physics/mass.h"


class msystem_combinator{
    msystem *pms;
    real t_split;

    std::map<uint64_t,std::vector<barycen>> sublists;
    friend class msystem;
public:
    msystem_combinator():pms(nullptr){}

    //build blist for a subsystem mssub of *pms
    //  s.t. for tidal_childlist of mssub,
    //  decomposed states in mssub is the same as in *pms
    //  i.e. extract a sub-blist rooted from tidal_parent
    int_t link(msystem &mssub,int_t tid=-1);
};
