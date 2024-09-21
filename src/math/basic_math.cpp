#include"basic_math.h"
#include<cmath>

double angle_reduce(double x){
    using Constants::pi;
    constexpr double pi_mul2l=2.44929359829470641e-16;
    if(x*0!=0)return NAN;
    while(x>pi||x<-pi){
        double rx=x/(2*pi);
        // faster than std::round by 7x under AVX
        // more robust and avoid infinite loop due to potential precision loss (force round up at 0.5, no round to even)
        // however note 0.49...95 may be rounded to 1. This will not cause error since we tested range of x.
        rx=rx>=0?std::floor(rx+0.5):std::ceil(rx-0.5);
        x=std::fma(-2*pi,rx,x)-rx*pi_mul2l;
    }
    return x;
}
