#pragma once
#include"definitions.h"

//parameters for orbital ephemerides
//see https://bridge.kamine.cloud/archives/421
class keplerian{
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

    keplerian()=default;
    //convert r and v to orbital parameters
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    keplerian(const vec &r,const vec &v);

    //convert orbital parameters to r and v
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    void rv(double t,vec &r,vec &v) const;

    //kep(x)
    static double kep(double);
    //dkep(x)/dx
    static double dkep(double);
};

