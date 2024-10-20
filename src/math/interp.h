#include"definitions.h"
#include<vector>

// maximum implemented degree d for bspline_* series functions
int_t bspline_basis_max_degree();

/*
    Basis Functions of BSpline(B).

    defined as:
        B[0,x] = x is in [0,1) ? 1 : 0;
        B[d,x] = ((x+d)*B[d-1,x+1]+(1-x)*B[d-1,x])/d;
    for 0<=d<=bspline_basis_max_degree().

    B[d,x] is a C^{d-1} continuous non-negative function of x on Reals,
        it is positive when x is in (-d,1), is 0 when x is not in [-d,1).

    if(pdb) also store derivative dB[d,x]/dx to *pdb.
*/
double bspline_basis(int_t d,double x,double *pdb=nullptr);

/*
    The Chebyshev Coefficients(c) of BSpline Basis Functions.
    for real x in [0,1) and integer k in [0,d],
        B[d,x-k] = sum_{j=0}^{d} c(d,k,j)*chebyshev_T(j,2*x-1);
*/
double bspline_basis_chebyshev_coef(int_t d,int_t k,int_t j);

//same as bspline_basis but use Chebyshev expansion for calculation,
//faster but loses relative accuracy when x is near boundaries of [-d,1).
double bspline_basis_chebyshev(int_t d,double x);
// also store derivative to *pdb
double bspline_basis_chebyshev(int_t d,double x,double *pdb);


template<typename T,size_t N_Channel>
class bspline_fitter{
    const int_t d,n,mn;
    std::vector<T> bsplines;
    std::vector<T> chebyshevs;
    bool is_valid;

    typedef dfloat_t<double> ddouble;
    typedef dfloat_t<T> dprec_t;
    static_assert(N_Channel>0,"Number of Channels should be positive.");

public:

    //Fit T source_data[mn+1][N_Channel] with n+d d-th-degree bspline basis functions.
    //Coefficients will be stored internally as T bspline_coefs[n+d][N_Channel].
    //need n>=1 && d>0 && n+d<=mn+1
    //Note with high degree d, the problem may be ill-conditioned when n+d is close to mn+1,
    // thence the fitting may be very poor or result to nan.
    bspline_fitter(int_t _d,int_t _n,const T *source_data,int_t _mn):d(_d),n(_n),mn(_mn){
        is_valid=n>0&&d>0&&mn+1>=n+d;
        if(!is_valid)return;

        std::vector<dprec_t> yvec(N_Channel*(n+d),dprec_t(0));

        //(l,u)=(d,d) diagonal ordered form of matrix
        std::vector<ddouble> mb((n+d)*(2*d+1),ddouble(0));

#define all_channels    int_t ich=0;ich<N_Channel;++ich
#define index(i)        (i)*N_Channel+ich
#define mat(i,j)        mb[(i)+d*(2*(j)+1)]
        dprec_t *y=yvec.data();

        for(int_t i=0;i<n+d;++i){
            int_t jmin=std::max(int_t(0),mn*(i-d)+n)/n;
            int_t jmax=std::min(mn,(mn*(i+1)-1)/n);
            const double rmn=double(1)/mn;
            ddouble w=0;
            for(int_t j=jmin;j<=jmax;++j){
                ddouble b=bspline_basis_chebyshev(d,j*n*rmn-i);
                int_t kmin=j*n/mn;
                int_t kmax=(j*n+mn-1)/mn+d-1;
                for(int_t k=kmin;k<=kmax;++k)
                    mat(i,k)+=bspline_basis_chebyshev(d,j*n*rmn-k)*b;
                for(all_channels)
                    y[index(i)]+=dprec_t(source_data[index(j)])*b;
                w+=b;
            }
            //normalize by w
            w=ddouble(1)/w;
            int_t kmin=std::max(int_t(0),i-d);
            int_t kmax=std::min(n+d-1,i+d);
            for(int_t k=kmin;k<=kmax;++k)
                mat(i,k)*=w;
            for(all_channels)
                y[index(i)]*=w;
        }

        //banded solver
        for(int_t i=0;i<n+d;++i){
            int_t jmin=std::max(int_t(0),i-d);
            int_t jmax=std::min(n+d-1,i+d);
            for(int_t j=jmin;j<i;++j){
                ddouble sum=0;
                for(int_t k=std::max(int_t(0),i-d);k<j;++k)
                    sum+=mat(i,k)*mat(k,j);
                mat(i,j)-=sum;
                mat(i,j)/=mat(j,j);
            }
            for(int_t j=i;j<=jmax;++j){
                ddouble sum=0;
                for(int_t k=std::max(int_t(0),j-d);k<i;++k)
                    sum+=mat(i,k)*mat(k,j);
                mat(i,j)-=sum;
            }
            for(int_t k=std::max(int_t(0),i-d);k<i;++k){
                for(all_channels)
                    y[index(i)]-=mat(i,k)*y[index(k)];
            }
        }
        for(int_t i=n+d-1;i>=0;--i){
            int_t jmax=std::min(n+d-1,i+d);
            for(int_t j=i+1;j<=jmax;++j){
                for(all_channels)
                    y[index(i)]-=mat(i,j)*y[index(j)];
            }
            ddouble w=ddouble(1)/mat(i,i);
            for(all_channels)
                y[index(i)]*=w;
        }

        bsplines.reserve(N_Channel*(n+d));
        for(int_t i=0;i<n+d;++i){
            for(all_channels)
                bsplines.push_back(T(y[index(i)]));
        }
#undef mat
        //must expand for d==1
        //otherwise discrete derivative will lead to errors on edge cases
        if(d==1)expand();
    }

    //Construct from fitted coefficients T bspline_coefs[n+d][N_Channel].
    //need n>=1 && d>0 && n+d<=mn+1
    bspline_fitter(const T *bspline_coefs,int_t _d,int_t _n,int_t _mn):d(_d),n(_n),mn(_mn){
        is_valid=n>0&&d>0&&mn+1>=n+d;
        if(!is_valid)return;
        std::vector<T>(bspline_coefs,bspline_coefs+N_Channel*(n+d)).swap(bsplines);
        //must expand for d==1
        //otherwise discrete derivative will lead to errors on edge cases
        if(d==1)expand();
    }

    //convert bspline coefs to chebyshev coefs to accelerate future samplings
    void expand(){
        if(!is_valid)return;
        chebyshevs.resize(0);
        chebyshevs.resize(N_Channel*(d+1)*n,T(0));
        T *cdata=chebyshevs.data();
        const T *bdata=bsplines.data();
        for(int_t xi=0;xi<n;++xi){
            for(int_t k=0;k<=d;++k)for(int_t j=0;j<=d;++j){
                const double c=bspline_basis_chebyshev_coef(d,k,j);
                for(all_channels)
                    cdata[index(j)]+=c*bdata[index(xi+k)];
            }
            cdata+=N_Channel*(d+1);
        }
    }

    //return fitted T bspline_coefs[n+d][N_Channel]
    const T *get_fitted_data() const{
        return is_valid?bsplines.data():nullptr;
    }

    //fit original data[mn] at a position in [0,mn]
    void operator ()(double position,T result[N_Channel]) const{
        double x=position/mn*n;
        if(!(is_valid&&0<=x&&x<=n)){
            for(all_channels)
                result[index(0)]=T(NAN);
            return;
        }

        if(chebyshevs.empty()){
            int_t i_min=std::max<int_t>(0,(int_t)std::floor(x));
            int_t i_max=std::min<int_t>(n,(int_t)std::ceil(x))+d-1;
            for(all_channels)
                result[index(0)]=T(0);
            const T *bdata=bsplines.data();
            for(int_t i=i_min;i<=i_max;++i){
                const double b_coef=bspline_basis_chebyshev(d,x-i);
                for(all_channels)
                    result[index(0)]+=b_coef*bdata[index(i)];
            }
        }
        else{
            int_t xi=(int_t)std::floor(x);
            xi-=xi==n;
            const T *cdata=chebyshevs.data()+N_Channel*(d+1)*xi;
            x-=xi;
            double xp=2*x-1;
            double x2=2*xp;
            double tm,t=1,tp=xp;

            for(all_channels)
                result[index(0)]=cdata[index(0)];
            for(int_t j=1;j<=d;++j){
                for(all_channels)
                    result[index(0)]+=cdata[index(j)]*tp;
                tm=t;
                t=tp;
                tp=x2*t-tm;
            }
        }
    }

    //fit original data[mn] at a position in [0,mn], with derivative
    void operator ()(double position,T result[N_Channel],T derivative[N_Channel]) const{
        double x=position/mn*n;
        double dx_factor=double(n)/mn;
        if(!(is_valid&&0<=x&&x<=n)){
            for(all_channels)
                result[index(0)]=T(NAN);
            for(all_channels)
                derivative[index(0)]=T(NAN);
            return;
        }

        if(chebyshevs.empty()){
            int_t i_min=std::max<int_t>(0,(int_t)std::floor(x));
            int_t i_max=std::min<int_t>(n,(int_t)std::ceil(x))+d-1;
            for(all_channels)
                result[index(0)]=T(0);
            for(all_channels)
                derivative[index(0)]=T(0);
            const T *bdata=bsplines.data();
            for(int_t i=i_min;i<=i_max;++i){
                double db_coef;
                const double b_coef=bspline_basis_chebyshev(d,x-i,&db_coef);
                for(all_channels)
                    result[index(0)]+=b_coef*bdata[index(i)];
                db_coef*=dx_factor;
                for(all_channels)
                    derivative[index(0)]+=db_coef*bdata[index(i)];
            }
        }
        else{
            int_t xi=(int_t)std::floor(x);
            xi-=xi==n;
            const T *cdata=chebyshevs.data()+N_Channel*(d+1)*xi;
            x-=xi;
            double xp=2*x-1;
            double x2=2*xp;
            double tm,t=1,tp=xp;
            double vm,v=0,vp=2;
            for(all_channels)
                result[index(0)]=cdata[index(0)];
            for(all_channels)
                derivative[index(0)]=T(0);
            for(int_t j=1;j<=d;++j){
                for(all_channels)
                    result[index(0)]+=cdata[index(j)]*tp;
                double dc_coef=vp*dx_factor;
                for(all_channels)
                    derivative[index(0)]+=cdata[index(j)]*dc_coef;
                tm=t;
                vm=v;
                t=tp;
                v=vp;
                tp=x2*t-tm;
                vp=x2*v+4*t-vm;
            }
        }
#undef index
#undef all_channels
    }
};
