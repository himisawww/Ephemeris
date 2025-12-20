#pragma once
#include"interp.h"
#include<algorithm>
#include"htl/vector.h"

template<typename T,int_t d>
class bsp_coef_eval{
    T _data[d+1][d+1]{};
    static constexpr T ijj(int_t j){ return T(j==0?2:1); }
public:
    constexpr bsp_coef_eval(){ _data[0][0]=T(1); }
    constexpr explicit bsp_coef_eval(const bsp_coef_eval<T,d-1> &c){
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
class bsp_coef_table{
    static_assert(d>=0,"Degree of BSpline function should be non-negative.");
    static constexpr auto _make_table(){
        if constexpr(d<=0)
            return bsp_coef_eval<T,0>{};
        else
            return bsp_coef_eval<T,d>{bsp_coef_table<T,d-1>::coef_table};
    }
public:
    static constexpr auto coef_table=_make_table();
    static const T *table_data(){ return coef_table.get_data(); }
};

static const double *bsp_coef_tables[]={
    bsp_coef_table<double,0>::table_data(),
    bsp_coef_table<double,1>::table_data(),
    bsp_coef_table<double,2>::table_data(),
    bsp_coef_table<double,3>::table_data(),
    bsp_coef_table<double,4>::table_data(),
    bsp_coef_table<double,5>::table_data(),
    bsp_coef_table<double,6>::table_data(),
    bsp_coef_table<double,7>::table_data(),
    bsp_coef_table<double,8>::table_data(),
    bsp_coef_table<double,9>::table_data(),
    bsp_coef_table<double,10>::table_data(),
    bsp_coef_table<double,11>::table_data(),
    bsp_coef_table<double,12>::table_data(),
    bsp_coef_table<double,13>::table_data(),
    bsp_coef_table<double,14>::table_data(),
    bsp_coef_table<double,15>::table_data(),
    bsp_coef_table<double,16>::table_data(),
    bsp_coef_table<double,17>::table_data(),
    bsp_coef_table<double,18>::table_data(),
    bsp_coef_table<double,19>::table_data(),
    bsp_coef_table<double,20>::table_data(),
    bsp_coef_table<double,21>::table_data(),
    bsp_coef_table<double,22>::table_data(),
    bsp_coef_table<double,23>::table_data(),
    bsp_coef_table<double,24>::table_data(),
    bsp_coef_table<double,25>::table_data(),
    bsp_coef_table<double,26>::table_data(),
    bsp_coef_table<double,27>::table_data(),
    bsp_coef_table<double,28>::table_data(),
    bsp_coef_table<double,29>::table_data(),
    bsp_coef_table<double,30>::table_data(),
    bsp_coef_table<double,31>::table_data(),
    bsp_coef_table<double,32>::table_data(),
    bsp_coef_table<double,33>::table_data(),
    bsp_coef_table<double,34>::table_data(),
    bsp_coef_table<double,35>::table_data(),
    bsp_coef_table<double,36>::table_data(),
    bsp_coef_table<double,37>::table_data(),
    bsp_coef_table<double,38>::table_data(),
    bsp_coef_table<double,39>::table_data()
};

static constexpr int_t bsp_maxdeg=sizeof(bsp_coef_tables)/sizeof(bsp_coef_tables[0])-1;

int_t bspline_basis_max_degree(){ return bsp_maxdeg; }

#define coef_unchecked(d,k,j) bsp_coef_tables[d][(k)*((d)+1)+(j)]
double bspline_basis_chebyshev_coef(int_t d,int_t k,int_t j){
    if(d<0||d>bsp_maxdeg||k<0||k>d||j<0||j>d)
        return 0;
    return coef_unchecked(d,k,j);
}

double bspline_basis_chebyshev(int_t d,double x){
    if(d<0||d>bsp_maxdeg||!(x<1&&x>=-d))
        return 0;
    int_t k=-(int_t)std::floor(x);
    double s=coef_unchecked(d,k,0);
    if(d>0){
        x+=k;
        double xp=2*x-1;
        s+=coef_unchecked(d,k,1)*xp;
        if(d>1){
            double x2=2*xp;
            double tm,t=1,tp=xp;
            for(int_t j=2;j<=d;++j){
                tm=t;
                t=tp;
                tp=x2*t-tm;
                s+=coef_unchecked(d,k,j)*tp;
            }
        }
    }
    return s;
}

double bspline_basis_chebyshev(int_t d,double x,double *pdb){
    double &ds=*pdb;
    ds=0;
    if(d<0||d>bsp_maxdeg||!(x<1&&x>=-d))
        return 0;
    int_t k=-(int_t)std::floor(x);
    double s=coef_unchecked(d,k,0);
    if(d>0){
        x+=k;
        double xp=2*x-1;
        const double c1=coef_unchecked(d,k,1);
        s+=c1*xp;
        ds=c1*2;
        if(d>1){
            double x2=2*xp;
            double tm,t=1,tp=xp;
            double vm,v=0,vp=2;
            for(int_t j=2;j<=d;++j){
                tm=t;
                vm=v;
                t=tp;
                v=vp;
                tp=x2*t-tm;
                vp=x2*v+4*t-vm;
                const double cj=coef_unchecked(d,k,j);
                s+=cj*tp;
                ds+=cj*vp;
            }
        }
    }
    return s;
}

double bspline_basis(int_t d,double x,double *pdb){
    if(pdb)*pdb=0;
    if(d<0||d>bsp_maxdeg||!(x<1&&x>=-d))
        return 0;
    int_t k=-(int_t)std::floor(x);
    x+=k;
    double bspdxmk[bsp_maxdeg+2];
    double *const b=&bspdxmk[1];
    b[-1]=0;
    b[0]=1;
    for(int_t i=1;i<=d;++i){
        b[i]=0;
        int_t jmax=std::min<int_t>(k,i);
        int_t jmin=std::max<int_t>(0,k+i-d);
        if(i==d&&pdb)*pdb=b[k-1]-b[k];
        for(int_t j=jmax;j>=jmin;--j)
            b[j]=((i-j+x)*b[j-1]+(1+j-x)*b[j])/i;
    }
    return b[k];
}

