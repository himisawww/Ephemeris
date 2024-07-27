#pragma once
#include"definitions.h"

//ring attractor
class ring{
public:
    //number of disks to sum
    int_t N;
    //GM_ratio: GM of ring system / GM of host body system(include ring)
    //A,J2: parameters of the ring system, same as in mass, using host R as reference radius
    //GL_R2: G * angular momentum of ring / GMR2 of host (exclude ring)
    fast_real GM_ratio,A,J2,GL_R2;

    struct {
        //GsR2: surface GM density of disk * R2/GM of host (exclude ring)
        //R_R: radius of disk / R of host
        //H_R: thickness of disk / R of host
        fast_real GsR2,R_R,H_R;
    } c_table[];
    ring()=delete;
    ring(ring &&)=delete;
    //ref_GM, ref_R : GM, R of host for normalization
    //direction_mass_factor:
    // if negative, the angular momentum of the ring is reversed,
    // mass of the ring system will be multiplied by its absolute value,
    // default to 1, ring is disabled when equal to 0.
    static const ring *load(const char *file,fast_real ref_GM,fast_real ref_R,fast_real direction_mass_factor=1);
    static void unload(const ring *);
    //multiplier: scale the intensity of original ring proportionally
    static const ring *copy(const ring *rp,fast_real multiplier=1);
    //sizeof this structure
    int_t size() const;
    //R: radius of the Extended Body
    //r: point mass position from the body which hosts the ring,
    //      represented in body-fixed coordinate system
    //      in which the ring lies on xy plane
    //return acceleration of point mass per host body GM,
    //      represented in body-fixed coordinate system
    fast_mpvec sum(fast_real R,fast_mpvec r) const;
    INLINE fast_mpvec cuda_sum(fast_real R,fast_mpvec r) const;
};
