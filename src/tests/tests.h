#pragma once
#include"physics/mass.h"
#include<vector>

const msystem &get_test_msystem();
msystem get_test_subsystem(std::vector<const char *> sids);
int test_all(uint64_t seed=0);
