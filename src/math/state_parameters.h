#pragma once
#include"math/keplerian.h"

//parameters for orbital ephemerides
class orbital_param_t:public keplerian{
public:
    //longitude of ascending node + earg
    double el;
    //eccentricity vector
    //e * sincos(el)
    double ex,ey;
    //mean longitude
    //el + mean anomaly
    double ml;

    orbital_param_t()=default;
    orbital_param_t(const keplerian &other);
    orbital_param_t(const double *k,bool circular);

    // test if two states are interpolatable
    static bool interpolatable(bool circular,const orbital_param_t &k1,const orbital_param_t &k2);

    //blend multiple keplerians by weights
    void blend_initialize(bool circular);
    void blend_add(bool circular,const orbital_param_t &component,double weight);
    void blend_finalize(bool circular);

};

//parameters for rotational ephemerides
class rotational_param_t{
public:
    //t = time relative to epoch
    //double t;
    //rotational angular velocity
    vec w;
    //polar wandering in surface-frame
    double ptheta,pphi;
    //rotation angle
    double angle;

    rotational_param_t()=default;
    //convert surface-frame and w to rotational parameters
    rotational_param_t(const mat &s,const vec &w);

    //convert rotational parameters to surface-frame and w
    void sw(double t,mat &s,vec &w) const;

    rotational_param_t(const double *a);
};
