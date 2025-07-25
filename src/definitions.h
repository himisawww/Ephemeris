#pragma once
#pragma fp_contract (off)
#include<cstdint>
#include"math/basic_math.h"
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

typedef quat_t<double> quat;
typedef quat_t<real> mpquat;
typedef quat_t<fast_real> fast_mpquat;

typedef int64_t int_t;

//compile switches


// pre-declarations
class MFILE;

class geopotential;
class ring;
