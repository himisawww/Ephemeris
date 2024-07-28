#pragma once
#include<vector>
#include"physics/mass.h"

const msystem &get_test_msystem();
msystem get_test_subsystem(std::vector<const char *> sids);
int test_all(uint64_t seed=0);
