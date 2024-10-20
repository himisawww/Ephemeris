#pragma once
#ifndef INLINE
#define INLINE inline
#endif

//consts
namespace Constants{

#define CONSTANT_VALUE_PI    3.1415926535897932
#define CONSTANT_VALUE_C     299792458

constexpr double pi=CONSTANT_VALUE_PI;
constexpr double pi_mul2=2*pi;
constexpr double pi_div2=pi/2;
constexpr double pi_div4=pi/4;
constexpr double degree=pi/180;

constexpr double c=CONSTANT_VALUE_C;
constexpr double c2=c*c;

constexpr double G=6.67430E-11;

//  1 year = 8766 h, Julian
constexpr double year=365.25*86400;

}

// return nan if one of parameters is nan
template<typename T>
INLINE T checked_max(const T &a,const T &b){
    return (b!=b||a<b)?b:a;
}
template<typename T>
INLINE T checked_min(const T &a,const T &b){
    return (b!=b||b<a)?b:a;
}
template<typename T>
INLINE void checked_maximize(T &a,const T &b){
    a=checked_max(b,a);
}
template<typename T>
INLINE void checked_minimize(T &a,const T &b){
    a=checked_min(b,a);
}

// return non-nan if one of parameters is nan
template<typename T>
INLINE T filtered_max(const T &a,const T &b){
    return (a!=a||a<b)?b:a;
}
template<typename T>
INLINE T filtered_min(const T &a,const T &b){
    return (a!=a||b<a)?b:a;
}
template<typename T>
INLINE void filtered_maximize(T &a,const T &b){
    a=filtered_max(b,a);
}
template<typename T>
INLINE void filtered_minimize(T &a,const T &b){
    a=filtered_min(b,a);
}

// reduce angle to [-pi,pi] by factors of 2pi
double angle_reduce(double x);
