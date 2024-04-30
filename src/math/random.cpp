#include"random.h"
#include<random>

static thread_local std::mt19937_64 g;

//Re-init with a given seed
void seedrandom(const uint64_t seed){
    g.seed(seed);
}

//a 64-bit random number
uint64_t random64(){
    return g();
}
//a random real in (0,1)
double randomreal(){
    const double rbase=1.0/18446744073709551616.;
    return ((double)random64()+0.5)*rbase;
}
//a random real from normal distribution
double randomnormal(){
    const double dpi=6.2831853071795864;
    return sqrt(-2*log(randomreal()))*cos(dpi*randomreal());
}

vec randomdirection(){
    double cost=randomreal()*2-1,sint=sqrt(1-cost*cost),phi=2*pi*randomreal();
    return vec(cos(phi)*sint,sin(phi)*sint,cost);
}
mat randommatrix(){
    vec rd=randomdirection();
    mat res(rd.perpunit(),0,rd);
    res.rotz(2*pi*randomreal());
    return res;
}
