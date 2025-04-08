#pragma once
#include<cmath>//for fma
#ifndef INLINE
#define INLINE inline
#endif

//minimal implementation of double-float
template<typename T>
class dfloat_impl_t{
public:
    T hi,lo;

    typedef dfloat_impl_t<T> dfloat_type;

    INLINE dfloat_impl_t(){}
    INLINE dfloat_impl_t(T _hi,T _lo=0):hi(_hi),lo(_lo){}
    INLINE explicit operator T() const{ return hi+lo; }

    INLINE dfloat_impl_t<T> &operator+=(const dfloat_impl_t<T> &x){ return (*this=*this+x); }
    INLINE dfloat_impl_t<T> &operator-=(const dfloat_impl_t<T> &x){ return (*this=*this-x); }
    INLINE dfloat_impl_t<T> &operator*=(const dfloat_impl_t<T> &x){ return (*this=*this*x); }
    INLINE dfloat_impl_t<T> &operator/=(const dfloat_impl_t<T> &x){ return (*this=*this/x); }

    INLINE bool operator> (const dfloat_impl_t<T> &x) const{ return T(*this-x)> 0; }
    INLINE bool operator< (const dfloat_impl_t<T> &x) const{ return T(*this-x)< 0; }
    INLINE bool operator==(const dfloat_impl_t<T> &x) const{ return T(*this-x)==0; }
    INLINE bool operator>=(const dfloat_impl_t<T> &x) const{ return T(*this-x)>=0; }
    INLINE bool operator<=(const dfloat_impl_t<T> &x) const{ return T(*this-x)<=0; }
    INLINE bool operator!=(const dfloat_impl_t<T> &x) const{ return T(*this-x)!=0; }

    friend INLINE dfloat_impl_t<T> operator+(dfloat_impl_t<T> x,dfloat_impl_t<T> y){
        dfloat_impl_t<T> re;
        twosum(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.lo+y.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_impl_t<T> operator-(dfloat_impl_t<T> x,dfloat_impl_t<T> y){
        dfloat_impl_t<T> re;
        twodiff(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.lo-y.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_impl_t<T> operator*(dfloat_impl_t<T> x,dfloat_impl_t<T> y){
        dfloat_impl_t<T> re;
        twoprod(x.hi,y.hi,re.hi,re.lo);
        re.lo+=x.hi*y.lo+y.hi*x.lo;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }

    friend INLINE dfloat_impl_t<T> operator/(dfloat_impl_t<T> x,dfloat_impl_t<T> y){
        dfloat_impl_t<T> re,sf;
        re.hi=x.hi/y.hi;
        twoprod(re.hi,y.hi,sf.hi,sf.lo);
        re.lo=(x.hi-sf.hi-sf.lo+x.lo-re.hi*y.lo)/y.hi;
        twosumq(re.hi,re.lo,re.hi,re.lo);
        return re;
    }


    friend INLINE dfloat_impl_t<T> operator +(dfloat_impl_t<T> x){
        return 0+x;
    }

    friend INLINE dfloat_impl_t<T> operator -(dfloat_impl_t<T> x){
        return 0-x;
    }

    friend INLINE dfloat_impl_t<T> sqrt(dfloat_impl_t<T> x){
        dfloat_impl_t<T> re=sqrt(x.hi);
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
