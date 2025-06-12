#pragma once
#include"physics/mass.h"

class ephemeris_reader{
public:
    ephemeris_reader(const char *ephemeris_path);

    operator bool() const{ return false; }
};
