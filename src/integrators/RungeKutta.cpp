#include"physics/mass.h"
#include"RungeKutta.impl"

void msystem::Cuda_RungeKutta12(fast_real dt,int_t n_step){
    RungeKutta12(dt,n_step);
}

void msystem::Cuda_accel(){
    accel();
}
