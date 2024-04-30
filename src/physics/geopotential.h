#pragma once
#include"definitions.h"
#include"utils/memfile.h"

#ifndef USE_NEW_GEOPOTENTIAL
//geopotential summer
class geopotential{
    //Nz: max zonal degree
    //Nt: max tesseral degree
    // N: max(Nz,Nt)
    int_t Nz,Nt,N;
    //    fast_real *J,*C,*S;
    //    fast_mpvec *c_table;
    fast_mpvec c_table[];

public:
    geopotential()=delete;
    //ref_radius_factor:
    // treat the coefficients as their reference radius is R*ref_radius_factor,
    // degree-n coefficients will be multiplied by ref_radius_factor^n,
    // default to 1, geopotential is disabled when equal to 0.
    static geopotential *load(mem_file &file,fast_real ref_radius_factor=1);
    static void unload(geopotential *);
    //sizeof this structure
    int_t size() const;
    //R: radius of the Extended Body
    //r: point mass position from the Extended Body,
    //      represented in body-fixed coordinate system
    //      in which the harmonics are expressed
    //N_start: starting degree of harmonics
    //N_end:   ending degree of harmonics
    //      default is -1, means to sum to maximum degree
    //return acceleration of point mass per Extended Body GM,
    //      represented in body-fixed coordinate system
    fast_mpvec sum(fast_real R,fast_mpvec r,int_t N_start=3,int_t N_end=-1) const;
    INLINE fast_mpvec cuda_sum(fast_real R,fast_mpvec r,int_t N_start=3,int_t N_end=-1) const;
};
#else
class geopotential{
    //Nz: max zonal degree
    //Nt: max tesseral degree
    int_t Nz,Nt;
    //    fast_real *J,*C,*S;
    fast_real c_table[];
public:

    geopotential()=delete;
    //ref_radius_factor:
    // treat the coefficients as their reference radius is R*ref_radius_factor,
    // degree-n coefficients will be multiplied by ref_radius_factor^n,
    // default to 1, geopotential is disabled when equal to 0.
    //N_start: starting degree of harmonics
    static geopotential *load(mem_file &file,fast_real ref_radius_factor=1,int_t N_start=3);
    static void unload(geopotential *);
    //sizeof this structure
    int_t size() const;
    //R: radius of the Extended Body
    //r: point mass position from the Extended Body,
    //      represented in body-fixed coordinate system
    //      in which the harmonics are expressed
    //return acceleration of point mass per Extended Body GM,
    //      represented in body-fixed coordinate system
    fast_mpvec sum(fast_real R,fast_mpvec r) const;
    INLINE fast_mpvec cuda_sum(fast_real R,fast_mpvec r) const;
};
#endif
