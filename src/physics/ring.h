#pragma once
#include"definitions.h"

//ring attractor
class ring{
public:
    //number of disks to sum
    int_t N;
    //GM_ratio: GM of ring system / GM of host body system(include ring)
    //GM,A,J2,GL: parameters of the ring system
    fast_real GM_ratio,A,J2,GL;

    struct {
        fast_real Gs,R,H;
    } c_table[];
    ring()=delete;
    //ref_GM, ref_R2 : GM, R2 of host for normalization
    //direction_mass_factor:
    // if negative, the angular momentum of the ring is reversed,
    // mass of the ring system will be multiplied by its absolute value,
    // default to 1, ring is disabled when equal to 0.
    static ring *load(const char *file,fast_real ref_GM,fast_real ref_R2,fast_real direction_mass_factor=1);
    static void unload(ring *);
    //sizeof this structure
    int_t size() const;
    //r: point mass position from the body which hosts the ring,
    //      represented in body-fixed coordinate system
    //      in which the ring lies on xy plane
    //return acceleration of point mass per host body GM,
    //      represented in body-fixed coordinate system
    fast_mpvec sum(fast_mpvec r) const;
    INLINE fast_mpvec cuda_sum(fast_mpvec r) const;
};
