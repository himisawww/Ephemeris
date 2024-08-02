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


// pre-declarations
class MFILE;

class geopotential;
class ring;

//consts
namespace Constants{

#define CONSTANT_VALUE_PI    3.1415926535897932
#define CONSTANT_VALUE_C     299792458

constexpr double pi         = CONSTANT_VALUE_PI;
constexpr double pi_mul2    = 2*pi;
constexpr double pi_div2    = pi/2;
constexpr double pi_div4    = pi/4;
constexpr double degree     = pi/180;

constexpr double c      = CONSTANT_VALUE_C;
constexpr double c2     = c*c;

constexpr double G      = 6.67430E-11;

//  1 year = 8766 h, Julian
constexpr double year   = 365.25*86400;

}

// return nan if one of parameters is nan
template<typename T>
const T &checked_max(const T &a,const T &b){
    if(b!=b)return b;
    if(a!=a)return a;
    return (a<b?b:a);
}
template<typename T>
const T &checked_min(const T &a,const T &b){
    if(b!=b)return b;
    if(a!=a)return a;
    return (b<a?b:a);
}
template<typename T>
void checked_maximize(T &a,const T &b){
    a=(b<a?a:b);
}
template<typename T>
void checked_minimize(T &a,const T &b){
    a=(a<b?a:b);
}
