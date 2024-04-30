#pragma once
#include<cmath>//for fma
#ifndef INLINE
#define INLINE inline
#endif

template<typename T>
INLINE void twosum(T x,T y,T &r,T &e){
    r=x+y;
    T t=r-x;
    e=(x-(r-t))+(y-t);
}
template<typename T>
INLINE void twosumq(T x,T y,T &r,T &e){
    r=x+y;
    e=y-(r-x);
}
template<typename T>
INLINE void twodiff(T x,T y,T &r,T &e){
    r=x-y;
    T t=r-x;
    e=(x-(r-t))-(y+t);
}
template<typename T>
INLINE void twoprod(T x,T y,T &r,T &e){
    r=x*y;
    e=fma(x,y,-r);
}

//minimal implementation of double-float
template<typename T>
class dfloat_t{
public:
    T hi,lo;

    INLINE dfloat_t(){}
    INLINE dfloat_t(T _hi,T _lo=0):hi(_hi),lo(_lo){}
    INLINE explicit operator T() const{ return hi+lo; }

    INLINE dfloat_t<T> &operator+=(const dfloat_t<T> &x){ return (*this=*this+x); }
    INLINE dfloat_t<T> &operator-=(const dfloat_t<T> &x){ return (*this=*this-x); }
    INLINE dfloat_t<T> &operator*=(const dfloat_t<T> &x){ return (*this=*this*x); }
    INLINE dfloat_t<T> &operator/=(const dfloat_t<T> &x){ return (*this=*this/x); }

    friend INLINE dfloat_t<T> operator+(dfloat_t<T> x,dfloat_t<T> y){
        dfloat_t<T> re;
        twosum(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.lo+y.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_t<T> operator-(dfloat_t<T> x,dfloat_t<T> y){
        dfloat_t<T> re;
        twodiff(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.lo-y.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_t<T> operator*(dfloat_t<T> x,dfloat_t<T> y){
        dfloat_t<T> re;
        twoprod(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.hi*y.lo+y.hi*x.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_t<T> operator/(dfloat_t<T> x,dfloat_t<T> y){
        dfloat_t<T> re,sf;
        re.hi=x.hi/y.hi;
        twoprod(re.hi,y.hi,sf.hi,sf.lo);
        re.lo=(x.hi-sf.hi-sf.lo+x.lo-re.hi*y.lo)/y.hi;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }


    friend INLINE dfloat_t<T> operator +(dfloat_t<T> x){
        return 0+x;
    }

    friend INLINE dfloat_t<T> operator -(dfloat_t<T> x){
        return 0-x;
    }

    friend INLINE dfloat_t<T> sqrt(dfloat_t<T> x){
        dfloat_t<T> re=sqrt(x.hi);
        re-=(re-x/re)/2;
        return re;
    }
};