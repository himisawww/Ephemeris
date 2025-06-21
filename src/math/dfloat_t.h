#pragma once
#include<cmath>//for fma
#include<type_traits>
#ifndef INLINE
#define INLINE inline
#endif

//minimal implementation of double-float
template<typename T>
class dfloat_impl_t{
    static_assert(std::is_same_v<T,float>||std::is_same_v<T,double>,"T must be float/double");
public:
    T hi,lo;

    typedef dfloat_impl_t<T> dfloat_type;

    INLINE dfloat_impl_t(){}
    INLINE dfloat_impl_t(T _hi,T _lo=0):hi(_hi),lo(_lo){}
    INLINE explicit operator T() const{ return hi+lo; }
    template<typename I,typename=std::enable_if_t<std::is_integral_v<I>>>
    INLINE explicit dfloat_impl_t(I i):hi(i){
        if constexpr(sizeof(I)>=sizeof(T))
            lo=i-I(hi);
        else
            lo=0;
    }
    //actually, floor
    template<typename I,typename=std::enable_if_t<std::is_integral_v<I>>>
    INLINE explicit operator I() const{
        I iapprox=std::round(hi+lo);
        iapprox+=(I)std::floor(T(*this-dfloat_type(iapprox)));
        return iapprox;
    }

    INLINE dfloat_type &operator+=(const dfloat_type &x){ return (*this=*this+x); }
    INLINE dfloat_type &operator-=(const dfloat_type &x){ return (*this=*this-x); }
    INLINE dfloat_type &operator*=(const dfloat_type &x){ return (*this=*this*x); }
    INLINE dfloat_type &operator/=(const dfloat_type &x){ return (*this=*this/x); }

    friend INLINE bool operator> (const dfloat_type &t,const dfloat_type &x){ return T(t-x)> 0; }
    friend INLINE bool operator< (const dfloat_type &t,const dfloat_type &x){ return T(t-x)< 0; }
    friend INLINE bool operator==(const dfloat_type &t,const dfloat_type &x){ return T(t-x)==0; }
    friend INLINE bool operator>=(const dfloat_type &t,const dfloat_type &x){ return T(t-x)>=0; }
    friend INLINE bool operator<=(const dfloat_type &t,const dfloat_type &x){ return T(t-x)<=0; }
    friend INLINE bool operator!=(const dfloat_type &t,const dfloat_type &x){ return T(t-x)!=0; }

    friend INLINE dfloat_type operator+(dfloat_type x,dfloat_type y){
        dfloat_type re;
        twosum(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.lo+y.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_type operator-(dfloat_type x,dfloat_type y){
        dfloat_type re;
        twodiff(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.lo-y.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_type operator*(dfloat_type x,dfloat_type y){
        dfloat_type re;
        twoprod(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.hi*y.lo+y.hi*x.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_type operator/(dfloat_type x,dfloat_type y){
        dfloat_type re,sf;
        re.hi=x.hi/y.hi;
        twoprod(re.hi,y.hi,sf.hi,sf.lo);
        re.lo=(x.hi-sf.hi-sf.lo+x.lo-re.hi*y.lo)/y.hi;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }


    friend INLINE dfloat_type operator +(dfloat_type x){
        return dfloat_type(+x.hi,+x.lo);
    }

    friend INLINE dfloat_type operator -(dfloat_type x){
        return dfloat_type(-x.hi,-x.lo);
    }

    friend INLINE dfloat_type abs(dfloat_type x){
        return T(x)>=0?+x:-x;
    }

    friend INLINE dfloat_type sqrt(dfloat_type x){
        dfloat_type re=sqrt(x.hi);
        if(re.hi>0)
            re-=(re-x/re)/2;
        return re;
    }

    static INLINE void twosum(T x,T y,T &r,T &e){
        r=x+y;
        T t=r-x;
        e=(x-(r-t))+(y-t);
    }
    static INLINE void twosumq(T x,T y,T &r,T &e){
        r=x+y;
        e=y-(r-x);
    }
    static INLINE void twodiff(T x,T y,T &r,T &e){
        r=x-y;
        T t=r-x;
        e=(x-(r-t))-(y+t);
    }
    static INLINE void twoprod(T x,T y,T &r,T &e){
        r=x*y;
        e=std::fma(x,y,-r);
    }
};

// make dfloat_t<vec_t<T>> same as vec_t<dfloat_t<T>>
template<typename T,template<typename> typename V>
struct dfloat_impl_t<V<T>>{
    typedef V<dfloat_impl_t<T>> dfloat_type;
};

template<typename T>
using dfloat_t=typename dfloat_impl_t<T>::dfloat_type;
