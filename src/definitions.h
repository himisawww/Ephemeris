#pragma once
#include<cstdint>

#include"math/dfloat_t.h"
typedef dfloat_t<double> real;

#define strtoreal atof

typedef double fast_real;

#include"math/vec_t.h"

typedef vec_t<double> vec;
typedef mat_t<double> mat;
typedef vec_t<real> mpvec;
typedef mat_t<real> mpmat;
typedef vec_t<fast_real> fast_mpvec;
typedef mat_t<fast_real> fast_mpmat;

typedef int64_t int_t;

//compile switches
#define USE_NEW_GEOPOTENTIAL

//consts
const double pi=3.1415926535897932;
const double degree=pi/180;

// pre-declarations
class mem_file;

class geopotential;
class ring;
