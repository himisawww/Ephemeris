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

    //adaptor to simpify usage with standard gravitional parameter
    class mu{
    public:
        double sqrt_mu;
        mu(double GM):sqrt_mu(std::sqrt(GM)){}
    };

    keplerian()=default;
    //convert r and v to orbital parameters
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    keplerian(const vec &r,const vec &v);
    keplerian(const mu &_mu,const vec &r,const vec &v):keplerian(r,v/_mu.sqrt_mu){}

    //transform keplerian to local axes
    static keplerian tolocal(const mat &frame,const keplerian &kep_world);
    //transform keplerian from local axes
    static keplerian toworld(const mat &frame,const keplerian &kep_local);

    //convert orbital parameters to r and v
    //note v = velocity / sqrt(GM)
    //     t = (time - t_epoch) * sqrt(GM)
    void rv(double t,vec &r,vec &v) const;
    void rv(const mu &_mu,double t,vec &r,vec &v){ rv(t*_mu.sqrt_mu,r,v);v*=_mu.sqrt_mu; }

    //note mean_motion() * sqrt(GM) for actual value, nan for q>0
    double mean_motion() const;
    double mean_motion(const mu &_mu) const{ return mean_motion()*_mu.sqrt_mu; }

    //kep(x)
    static double kep(double);
    //dkep(x)/dx
    static double dkep(double);
};

