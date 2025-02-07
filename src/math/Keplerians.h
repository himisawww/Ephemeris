#pragma once
#include"definitions.h"

//parameters for orbital ephemerides
//see https://bridge.kamine.cloud/archives/421
class ephem_orb{
public:
    //t = (time relative to epoch) * sqrt(GM)
    //double t;
    //j = (angular momentum per mass) / sqrt(GM)
    vec j;
    //q = e-1 , defined to reduce numeric error at e~0 and e~1
    double q;
    //earg: argument of eccentricity vector in j-frame
    double earg;
    //m = (mean anomaly) / sqrt(1-e^2)^3
    //dm/dt = 1/j^3
    double m;

    ephem_orb()=default;
    //convert r and v to orbital parameters
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    ephem_orb(const vec &r,const vec &v);

    //convert orbital parameters to r and v
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    void rv(double t,vec &r,vec &v) const;
};
//parameters for rotational ephemerides
class ephem_rot{
public:
    //t = time relative to epoch
    //double t;
    //rotational angular velocity
    vec w;
    //polar wandering in surface-frame
    double ptheta,pphi;
    //rotation angle
    double angle;

    ephem_rot()=default;
    //convert surface-frame and w to rotational parameters
    ephem_rot(const mat &s,const vec &w);

    //convert rotational parameters to surface-frame and w
    void sw(double t,mat &s,vec &w) const;
};
