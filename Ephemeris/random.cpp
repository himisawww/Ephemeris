#include"utils.h"
#include<random>

double randomreal(){
    static std::mt19937_64 g64(0);
    static std::uniform_real_distribution<double> dist(0.0,1.0);
    return dist(g64);
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
