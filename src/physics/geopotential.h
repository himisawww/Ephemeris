#pragma once
#include"definitions.h"

class geopotential{
    //Nz: max zonal degree
    //Nt: max tesseral degree
    int_t Nz,Nt;
    //    fast_real *J,*C,*S;
    fast_real c_table[];
public:

    geopotential()=delete;
    geopotential(geopotential &&)=delete;
    //ref_radius_factor:
    // treat the coefficients as their reference radius is R*ref_radius_factor,
    // degree-n coefficients will be multiplied by ref_radius_factor^n,
    // default to 1, geopotential is disabled when equal to 0.
    //N_start: starting degree of harmonics
    static const geopotential *load(const char *file,fast_real ref_radius_factor=1,int_t N_start=3);
    static void unload(const geopotential *);
    //multiplier: scale the intensity of original geopotential proportionally
    static const geopotential *copy(const geopotential *,fast_real multiplier=1);
    static bool is_same(const geopotential *,const geopotential *);
    //sizeof this structure
    int_t size() const;
    int_t max_degree() const{ return Nz; }
    double get_J(int_t n) const;
    double get_C(int_t n,int_t m) const;
    double get_S(int_t n,int_t m) const;
    //R: radius of the Extended Body
    //r: point mass position from the Extended Body,
    //      represented in body-fixed coordinate system
    //      in which the harmonics are expressed
    //return acceleration of point mass per Extended Body GM,
    //      represented in body-fixed coordinate system
    fast_mpvec sum(fast_real R,fast_mpvec r) const;
    INLINE fast_mpvec cuda_sum(fast_real R,fast_mpvec r) const;
};
