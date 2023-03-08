#include"ephemeris.h"

#define CONST_TABLE static const
#include"RungeKutta.impl"

void msystem::RungeKutta12(fast_real dt,int_t n_step){
    static const real *clist=(const real *)rk12_coefs;

    //prepare space for states
    std::vector<mass> &x=mlist;
    int_t mn=x.size();
    static thread_local std::vector<mass_state> x0,f;
    x0.resize(mn);
    f.resize(25*mn);
    
    for(int_t i_step=0;i_step<n_step;++i_step){
        for(int_t i=0;i<mn;++i)x0[i]=x[i];

        fast_real dt_k;
        for(int_t k=1;k<=25;++k){
            dt_k=0;
            for(int_t i=0;i<mn;++i){
                mass_state &fi=f[mn*(k-1)+i];
                mass &xi=x[i];
                fi.v=xi.v;
                fi.naccel=xi.gaccel+xi.daccel+xi.naccel;
                fi.w=xi.w;
                fi.dtorque=xi.dtorque;

                if(k==25){
                    fast_mpvec j(xi.GL),&w=xi.w;
                    fast_mpvec wxj=w*j;
                    fast_real wxj2=wxj%wxj;
                    if(wxj2!=0){
                        fast_mpmat ir=xi.GI.inverse();
                        fast_real dt2=dt*dt;
                        fast_real e=w%j*(fast_real(1)/2);
                        fast_real de=(ir%w)%((ir%wxj)*j+w*wxj);
                        de*=e*dt2*dt2/(36*wxj2);
                        xi.Erot=de;
                        xi.Egrad=wxj;
                    }
                    else{
                        xi.Erot=0;
                        xi.Egrad=0;
                    }
                }

                xi.v=0;
                xi.naccel=0;
                xi.w=0;
                xi.dtorque=0;
            }

            for(int_t j=0;j<k;++j){
                const real &ckj=clist[k*(k-1)/2+j];
                const fast_real fckj=(fast_real)ckj;
                if(ckj.hi){
                    dt_k+=fckj;
                    for(int_t i=0;i<mn;++i){
                        mass_state &fi=f[mn*j+i];
                        mass &xi=x[i];
                        xi.v+=ckj*fi.v;
                        xi.naccel+=fckj*fi.naccel;
                        xi.w+=fckj*fi.w;
                        xi.dtorque+=fckj*fi.dtorque;
                    }
                }
            }
            for(int_t i=0;i<mn;++i){
                const mass_state &xi0=x0[i];
                mass &xi=x[i];

                xi.r=xi0.r+real(dt)*xi.v;
                xi.v=xi0.v+mpvec(dt*xi.naccel);
                xi.w/=dt_k;
                xi.s=xi0.s;
                fast_mpvec dw=(fast_real(2)/3)*(xi.w-xi0.w);
                xi.s+=rotation_matrix(xi.w-dw,dt*dt_k/2)%fast_mpmat(xi.s);
                xi.s+=rotation_matrix(xi.w+dw,dt*dt_k/2)%fast_mpmat(xi.s);
                xi.GL=xi0.GL+mpvec(dt*xi.dtorque);
            
                if(k==25){
                    xi.s+=rotation_matrix(xi.Egrad,xi.Erot)%fast_mpmat(xi.s);
                }
            }
            
            update(fast_real(t_eph)+dt*(dt_k+i_step));
            accel();
        }
    }
}