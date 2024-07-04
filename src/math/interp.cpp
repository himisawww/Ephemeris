#pragma once
#include"interp.h"
#include<algorithm>
#include<vector>

template<typename T,int_t d>
class _bsp_coef_eval{
    T _data[d+1][d+1]{};
    static constexpr T ijj(int_t j){ return T(j==0?2:1); }
public:
    constexpr _bsp_coef_eval(){ _data[0][0]=T(1); }
    constexpr explicit _bsp_coef_eval(const _bsp_coef_eval<T,d-1> &c){
        for(int_t k=0;k<=d;++k){
            T bjm{};
            T bj=c(k-1,1)-c(k,1);
            T bjp=c(k-1,0)-c(k,0);
            for(int_t j=0;j<=d;++j){
                T aj=(d-k)*c(k-1,j)+(k+1)*c(k,j);
                bjm=bj;
                bj=bjp;
                bjp=c(k-1,j+1)-c(k,j+1);
                aj+=bj/T(2);
                _data[k][j]=(ijj(j)*aj+(ijj(j-1)*bjm+ijj(j+1)*bjp)/T(4))/(T(d)*ijj(j));
            }
        }
    }
    constexpr T operator()(int_t k,int_t j) const{
        if(k<0||k>d||j<0||j>d)
            return T(0);
        return _data[k][j];
    }
    const T *get_data() const{ return _data[0]; }
};

template<typename T,int_t d>
class _bsp_coef_table{
    static_assert(d>=0,"Degree of BSpline function should be non-negative.");
    static constexpr auto _make_table(){
        if constexpr(d<=0)
            return _bsp_coef_eval<T,0>{};
        else
            return _bsp_coef_eval<T,d>{_bsp_coef_table<T,d-1>::coef_table};
    }
public:
    static constexpr auto coef_table=_make_table();
    static const T *table_data(){ return coef_table.get_data(); }
};

static const double *_bsp_coef_tables[]={
    _bsp_coef_table<double,0>::table_data(),
    _bsp_coef_table<double,1>::table_data(),
    _bsp_coef_table<double,2>::table_data(),
    _bsp_coef_table<double,3>::table_data(),
    _bsp_coef_table<double,4>::table_data(),
    _bsp_coef_table<double,5>::table_data(),
    _bsp_coef_table<double,6>::table_data(),
    _bsp_coef_table<double,7>::table_data(),
    _bsp_coef_table<double,8>::table_data(),
    _bsp_coef_table<double,9>::table_data(),
    _bsp_coef_table<double,10>::table_data(),
    _bsp_coef_table<double,11>::table_data(),
    _bsp_coef_table<double,12>::table_data(),
    _bsp_coef_table<double,13>::table_data(),
    _bsp_coef_table<double,14>::table_data(),
    _bsp_coef_table<double,15>::table_data(),
    _bsp_coef_table<double,16>::table_data(),
    _bsp_coef_table<double,17>::table_data(),
    _bsp_coef_table<double,18>::table_data(),
    _bsp_coef_table<double,19>::table_data(),
    _bsp_coef_table<double,20>::table_data(),
    _bsp_coef_table<double,21>::table_data(),
    _bsp_coef_table<double,22>::table_data(),
    _bsp_coef_table<double,23>::table_data(),
    _bsp_coef_table<double,24>::table_data(),
    _bsp_coef_table<double,25>::table_data(),
    _bsp_coef_table<double,26>::table_data(),
    _bsp_coef_table<double,27>::table_data(),
    _bsp_coef_table<double,28>::table_data(),
    _bsp_coef_table<double,29>::table_data(),
    _bsp_coef_table<double,30>::table_data(),
    _bsp_coef_table<double,31>::table_data(),
    _bsp_coef_table<double,32>::table_data(),
    _bsp_coef_table<double,33>::table_data(),
    _bsp_coef_table<double,34>::table_data(),
    _bsp_coef_table<double,35>::table_data(),
    _bsp_coef_table<double,36>::table_data(),
    _bsp_coef_table<double,37>::table_data(),
    _bsp_coef_table<double,38>::table_data(),
    _bsp_coef_table<double,39>::table_data()
};

static constexpr int_t _bsp_maxdeg=sizeof(_bsp_coef_tables)/sizeof(_bsp_coef_tables[0])-1;

int_t bspline_basis_max_degree(){ return _bsp_maxdeg; }

#define _coef_unchecked(d,k,j) _bsp_coef_tables[d][(k)*((d)+1)+(j)]
double bspline_basis_chebyshev_coef(int_t d,int_t k,int_t j){
    if(d<0||d>_bsp_maxdeg||k<0||k>d||j<0||j>d)
        return 0;
    return _coef_unchecked(d,k,j);
}

double bspline_basis_chebyshev(int_t d,double x){
    if(d<0||d>_bsp_maxdeg||!(x<1&&x>=-d))
        return 0;
    int_t k=-(int_t)std::floor(x);
    double s=_coef_unchecked(d,k,0);
    if(d>0){
        x+=k;
        double xp=2*x-1;
        s+=_coef_unchecked(d,k,1)*xp;
        if(d>1){
            double x2=2*xp;
            double tm,t=1,tp=xp;
            for(int_t j=2;j<=d;++j){
                tm=t;
                t=tp;
                tp=x2*t-tm;
                s+=_coef_unchecked(d,k,j)*tp;
            }
        }
    }
    return s;
}

double bspline_basis(int_t d,double x){
    if(d<0||d>_bsp_maxdeg||!(x<1&&x>=-d))
        return 0;
    int_t k=-(int_t)std::floor(x);
    x+=k;
    double bspdxmk[_bsp_maxdeg+2];
    double *const b=&bspdxmk[1];
    b[-1]=0;
    b[0]=1;
    for(int_t i=1;i<=d;++i){
        b[i]=0;
        int_t jmax=std::max<int_t>(k,i);
        int_t jmin=std::max<int_t>(0,k+i-d);
        for(int_t j=jmax;j>=jmin;--j)
            b[j]=((i-j+x)*b[j-1]+(1+j-x)*b[j])/i;
    }
    return b[k];
}

template<typename T,size_t N_Channel>
bool interp_bspline(T *result,const int_t d,const int_t n,const T *data,const int_t mn){
    static_assert(N_Channel>0,"Number of Channels should be positive.");
    if(n<1||d<1||!(n+d<=mn+1))return false;

    //(l,u)=(d,d) diagonal ordered form of matrix
    std::vector<double> mb((n+d)*(2*d+1),double(0));

#define mat(i,j)        mb[(i)+d*(2*(j)+1)]
#define all_channels    int_t ich=0;ich<N_Channel;++ich
#define index(i)        (i)*N_Channel+ich
    T *y=result;

    for(int_t i=0;i<n+d;++i){
        int_t jmin=std::max(int_t(0),mn*(i-d)+n)/n;
        int_t jmax=std::min(mn,(mn*(i+1)-1)/n);
        const double rmn=double(1)/mn;
        double w=0;
        for(int_t j=jmin;j<=jmax;++j){
            double b=bspline_basis_chebyshev(d,j*n*rmn-i);
            int_t kmin=j*n/mn;
            int_t kmax=(j*n+mn-1)/mn+d-1;
            for(int_t k=kmin;k<=kmax;++k)
                mat(i,k)+=bspline_basis_chebyshev(d,j*n*rmn-k)*b;
            for(all_channels)
                y[index(i)]+=data[index(j)]*b;
            w+=b;
        }
        //normalize by w
        w=1/w;
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
            double sum=0;
            for(int_t k=std::max(int_t(0),i-d);k<j;++k)
                sum+=mat(i,k)*mat(k,j);
            mat(i,j)-=sum;
            mat(i,j)/=mat(j,j);
        }
        for(int_t j=i;j<=jmax;++j){
            double sum=0;
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
        double w=1/mat(i,i);
        for(all_channels)
            y[index(i)]*=w;
    }

#undef index
#undef all_channels
#undef mat

    return true;
}
