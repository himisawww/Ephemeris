#pragma once
#include"htl/vector.h"
#include"physics/mass.h"

const msystem &get_test_msystem();
msystem get_test_subsystem(htl::vector<const char *> sids);
int test_all(uint64_t seed=0);
