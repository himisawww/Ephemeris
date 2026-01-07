#pragma once
#include"definitions.h"
#include<vector>
#include"utils/memio.h"

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
    const int_t d,n;
    const double range;
    //offset of bspline_coefs in mdata in bytes
    size_t mdata_offset;
    MFILE *mdata;
    //offset of bsplines.data() in bspline_coefs in data_size
    int_t bsplines_offset;
    //load() will at least expand cache to this threshold(in data_size)
    int_t min_cache_pts;
    std::vector<T> bsplines;
    //if present, full bspline_coefs are dumped to bsplines and expanded to chebyshevs
    std::vector<T> chebyshevs;
    bool is_valid;

    typedef dfloat_t<double> ddouble;
    typedef dfloat_t<T> dprec_t;
    static_assert(N_Channel>0,"Number of Channels should be positive.");
public:
    typedef T data_type;
    static constexpr size_t data_channel=N_Channel;
    //sizeof(data_type)*data_count
    static constexpr size_t data_size=sizeof(T)*N_Channel;

    //Fit T source_data[mn+1][N_Channel] with n+d d-th-degree bspline basis functions.
    //Coefficients will be stored internally as T bspline_coefs[n+d][N_Channel].
    //need n>=1 && d>0 && n+d<=mn+1
    //Note with high degree d, the problem may be ill-conditioned when n+d is close to mn+1,
    // thence the fitting may be very poor or result to nan.
    bspline_fitter(int_t _d,int_t _n,const T *source_data,const int_t mn)
        :d(_d),n(_n),range(double(mn)),mdata(nullptr){
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
    //need n>=1 && d>0 && range!=0,+-inf
    bspline_fitter(const T *bspline_coefs,int_t _d,int_t _n,double _range)
        :d(_d),n(_n),range(_range),mdata(nullptr){
        is_valid=n>0&&d>0&&(range!=0&&range*0==0);
        if(!is_valid)return;
        std::vector<T>(bspline_coefs,bspline_coefs+N_Channel*(n+d)).swap(bsplines);
        //must expand for d==1
        //otherwise discrete derivative will lead to errors on edge cases
        if(d==1)expand();
    }

    //Construct from fitted coefficients T bspline_coefs[n+d][N_Channel] that stored in _file at _offset.
    // caller must ensure _file is valid during lifetime of bspline_fitter
    //need n>=1 && d>0 && range!=0,+-inf
    bspline_fitter(MFILE *_file,size_t _offset,int_t _d,int_t _n,double _range,int_t _cache_bytes=0)
        :d(_d),n(_n),range(_range),mdata_offset(_offset),mdata(_file),min_cache_pts(_cache_bytes/data_size){
        is_valid=false;
        do{
            if(!(n>0&&d>0&&(range!=0&&range*0==0)))return;
            if(!(mdata&&mdata->is_read()))return;
            fseek(mdata,0,SEEK_END);
            size_t fsize=ftell(mdata);
            size_t max_size=data_size*(n+d);
            if(_offset>=fsize||_offset+max_size>fsize)return;
        } while(0);
        is_valid=true;
        bsplines_offset=0;
        //must expand for d==1
        //otherwise discrete derivative will lead to errors on edge cases
        if(d==1)expand();
    }

private:
    //assume valid, 0<=i_begin<i_end<=n+d
    const T *load_data(int_t i_begin,int_t i_end){
        MFILE *fdata=mdata;
        if(!fdata)return bsplines.data();
        int_t old_size=bsplines.size()/N_Channel;
        if(bsplines_offset<=i_begin&&i_end<=bsplines_offset+old_size)
            return bsplines.data()-bsplines_offset*N_Channel;
        int_t load_size=std::max((i_end-i_begin)*3+6*d,std::max(old_size,min_cache_pts));
        int_t load_begin=(i_begin+i_end-load_size)/2;
        int_t load_excess=load_begin+load_size-(n+d);
        if(load_excess>0)
            load_begin-=load_excess;
        if(load_begin<=0){
            load_begin=0;
            if(load_size>=n+d){
                load_size=n+d;
                mdata=nullptr;
            }
        }
        bsplines.resize(load_size*N_Channel);
        bsplines_offset=load_begin;
        fseek(fdata,mdata_offset+load_begin*data_size,SEEK_SET);
        fread(bsplines.data(),data_size,load_size,fdata);
        return bsplines.data()-load_begin*N_Channel;
    }
    const T *load_data(){
        return load_data(0,n+d);
    }
public:
    explicit operator bool() const{ return is_valid; }

    //convert bspline coefs to chebyshev coefs to accelerate future samplings
    void expand(){
        if(!is_valid)return;
        chebyshevs.resize(0);
        chebyshevs.resize(N_Channel*(d+1)*n,T(0));
        T *cdata=chebyshevs.data();
        const T *bdata=load_data();
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
    const T *get_fitted_data(){
        return is_valid?load_data():nullptr;
    }

    //return heap size of fitter
    int_t memory_size() const{
        return bsplines.size()*sizeof(T)+chebyshevs.size()*sizeof(T);
    }

    //fit original data[mn] at a position in [0,mn]
    void operator ()(double position,T result[N_Channel]){
        double x=position/range*n;
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
            const T *bdata=load_data(i_min,i_max+1);
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
    void operator ()(double position,T result[N_Channel],T derivative[N_Channel]){
        double x=position/range*n;
        double dx_factor=double(n)/range;
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
            const T *bdata=load_data(i_min,i_max+1);
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
