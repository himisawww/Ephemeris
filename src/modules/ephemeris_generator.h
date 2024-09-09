#pragma once
#include<mutex>
#include"physics/mass.h"
#include"utils/memio.h"

class ephemeris_collector{
private:
    msystem &ms;
    std::vector<barycen> blist;
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

    void extract(std::vector<MFILE> &ephm_files,bool force);
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
