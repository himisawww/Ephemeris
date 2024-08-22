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
    void update_barycens();
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
