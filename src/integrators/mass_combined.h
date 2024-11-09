#pragma once
#include"physics/mass.h"


class msystem_combinator{
    msystem *pms;
    real t_split;

    friend class msystem;
public:
    msystem_combinator():pms(nullptr){}

    void split_system(msystem &ms);
};
