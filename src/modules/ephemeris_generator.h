#pragma once
#include"physics/mass.h"
#include<mutex>

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
