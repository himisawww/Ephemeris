#pragma once
#include<cmath>//for trigs
#include<cfloat>
#include<tuple>
#ifndef INLINE
#define INLINE inline
#endif

template<typename T>
struct hypot_traits;

template<>
struct hypot_traits<float>{
    typedef float float_t;
    static constexpr float_t high=0x1p100;
    static constexpr float_t mul=0x1p96;
};

template<>
struct hypot_traits<double>{
    typedef double float_t;
    static constexpr float_t high=0x1p960;
    static constexpr float_t mul=0x1p768;
};

// make hypot_traits<dfloat_t<T>> same as hypot_traits<T>
template<typename T,template<typename> typename V>
struct hypot_traits<V<T>>:public hypot_traits<T>{};

template<typename T>
class hypot_t:private hypot_traits<T>{
    typedef hypot_traits<T> base_t;
    using typename base_t::float_t;
    using base_t::high;
    using base_t::mul;
    static constexpr float_t div=1/mul;
    static constexpr float_t low=1/high;
    template<typename ...Args>
    static INLINE T _length(Args ...x){
        T r=(...+(x*x));
        if(r<low){
            ((x*=mul),...);
            r=(...+(x*x));
            return div*sqrt(r);
        }
        if(!(r<high)){
            ((x*=div),...);
            r=(...+(x*x));
            return mul*sqrt(r);
        }
        return sqrt(r);
    }
    template<typename ...Args,size_t ...Is>
    static INLINE T _normalize(std::integer_sequence<size_t,Is...>,Args &..._x){
        std::tuple<Args...> x{_x...};
        T r=(...+(std::get<Is>(x)*std::get<Is>(x)));
        if(r<low){
            ((std::get<Is>(x)*=mul),...);
            r=(...+(std::get<Is>(x)*std::get<Is>(x)));
            if(r==0){
                ((_x=Is+1==sizeof...(Args)),...);
                return r;
            }
            T s=sqrt(r),rs=1/s;
            ((_x=std::get<Is>(x)*rs),...);
            return div*s;
        }
        if(!(r<high)){
            ((std::get<Is>(x)*=div),...);
            r=(...+(std::get<Is>(x)*std::get<Is>(x)));
            T s=sqrt(r),rs=1/s;
            ((_x=std::get<Is>(x)*rs),...);
            return mul*s;
        }
        T s=sqrt(r),rs=1/s;
        ((_x=std::get<Is>(x)*rs),...);
        return s;
    }
public:
    template<typename ...Args>
    static INLINE T length(const Args &...x){
        return _length(T(x)...);
    }
    template<typename ...Args>
    static INLINE T normalize(Args &...x){
        static_assert((...&&(std::is_same_v<T,Args>)),"All parameters should be T &");
        return _normalize(std::index_sequence_for<Args...>{},x...);
    }
};

template<typename T>
class mat_t;

template<typename T>
class quat_t;

template<typename T>
class vec_t{
public:
    T x,y,z;

    INLINE vec_t(){}
    explicit INLINE vec_t(const T &x):x(x),y(x),z(x){}
    INLINE vec_t(const T &x,const T &y,const T &z):x(x),y(y),z(z){}
    INLINE vec_t(const T &theta,const T &phi){
        T c=sin(theta);
        x=c*cos(phi);
        y=c*sin(phi);
        z=cos(theta);
    }
    template<typename T2>
    INLINE vec_t(const vec_t<T2> &a):x(a.x),y(a.y),z(a.z){}


    INLINE vec_t<T> &rotx(const T &b){
        T s=sin(b),c=cos(b),ny=y*c-z*s;
        z=z*c+y*s;
        y=ny;
        return *this;
    }
    INLINE vec_t<T> &roty(const T &b){
        T s=sin(b),c=cos(b),nz=z*c-x*s;
        x=x*c+z*s;
        z=nz;
        return *this;
    }
    INLINE vec_t<T> &rotz(const T &b){
        T s=sin(b),c=cos(b),nx=x*c-y*s;
        y=y*c+x*s;
        x=nx;
        return *this;
    }
    INLINE T theta() const{
        return atan2(hypot_t<T>::length(x,y),z);
    }
    INLINE T phi() const{
        return atan2(y,x);
    }
    INLINE T norm() const{
        return hypot_t<T>::length(x,y,z);
    }
    INLINE T normsqr() const{
        return x*x+y*y+z*z;
    }
    //normalize *this by norm()
    //  *this = (0,0,1) for (0,0,0)
    //  return norm()
    INLINE T normalize(){
        return hypot_t<T>::normalize(x,y,z);
    }
    INLINE vec_t<T> unit() const{
        vec_t<T> result(x,y,z);
        result.normalize();
        return result;
    }

    INLINE vec_t<T> &operator =(const T &a){
        x=a;
        y=a;
        z=a;
        return *this;
    }
    template<typename T2>
    INLINE vec_t<T> &operator =(const vec_t<T2> &a){
        x=a.x;
        y=a.y;
        z=a.z;
        return *this;
    }

    INLINE vec_t<T> &operator *=(const T &a){
        x*=a;
        y*=a;
        z*=a;
        return *this;
    }
    INLINE vec_t<T> &operator /=(const T &a){
        T reca=1/a;
        x*=reca;
        y*=reca;
        z*=reca;
        return *this;
    }
    INLINE vec_t<T> &operator +=(const vec_t<T> &b){
        x+=b.x;
        y+=b.y;
        z+=b.z;
        return *this;
    }
    INLINE vec_t<T> &operator -=(const vec_t<T> &b){
        x-=b.x;
        y-=b.y;
        z-=b.z;
        return *this;
    }
    INLINE vec_t<T> &operator *=(const vec_t<T> &b){
        T _x=y*b.z-z*b.y,_y=z*b.x-x*b.z;
        z=x*b.y-y*b.x;
        x=_x;
        y=_y;
        return *this;
    }
    //ascending node: a perpendicular unit vector on xy-plane
    //  return (1,0,0) for (0,0,*)
    //  return exactly same vectors for multiple calls from exactly same vectors
    INLINE vec_t<T> asc_node() const{
        vec_t<T> p(-y,x,0);
        hypot_t<T>::normalize(p.y,p.x);
        return p;
    }

    //return matrix R s.t.
    //(I+R)%x rotates vector x along axis w counterclockwise by angle norm()*t
    INLINE mat_t<T> rotation_matrix(const T &t) const;

    friend INLINE vec_t<T> operator +(const vec_t<T> &a){
        return vec_t<T>(+a.x,+a.y,+a.z);
    }
    friend INLINE vec_t<T> operator -(const vec_t<T> &a){
        return vec_t<T>(-a.x,-a.y,-a.z);
    }
    friend INLINE vec_t<T> operator +(const vec_t<T> &a,const vec_t<T> &b){
        return vec_t<T>(a.x+b.x,a.y+b.y,a.z+b.z);
    }
    friend INLINE vec_t<T> operator -(const vec_t<T> &a,const vec_t<T> &b){
        return vec_t<T>(a.x-b.x,a.y-b.y,a.z-b.z);
    }
    friend INLINE vec_t<T> operator *(const vec_t<T> &a,const vec_t<T> &b){
        return vec_t<T>(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
    }
    friend INLINE T operator %(const vec_t<T> &a,const vec_t<T> &b){
        return a.x*b.x+a.y*b.y+a.z*b.z;
    }

    friend INLINE vec_t<T> operator *(const T &a,const vec_t<T> &b){
        return vec_t<T>(a*b.x,a*b.y,a*b.z);
    }
    friend INLINE vec_t<T> operator *(const vec_t<T> &b,const T &a){
        return vec_t<T>(a*b.x,a*b.y,a*b.z);
    }

    friend INLINE vec_t<T> operator /(const vec_t<T> &b,const T &a){
        T reca=1/a;
        return vec_t<T>(b.x*reca,b.y*reca,b.z*reca);
    }
};

//mat_t: a coordinate represented by its 3 axes xyz.
template<typename T>
class mat_t{
public:
    vec_t<T> x,y,z;
    INLINE mat_t(){}
    INLINE mat_t(const T &a):x(a,0,0),y(0,a,0),z(0,0,a){}
    INLINE mat_t(const vec_t<T> &a,const vec_t<T> &b):x(a*b.x),y(a*b.y),z(a*b.z){}
    INLINE mat_t(const vec_t<T> &x0,const vec_t<T> &y0,const vec_t<T> &z0):x(x0),y(y0),z(z0){}
    INLINE mat_t(const vec_t<T> &x0,const vec_t<T> &y0,int):x(x0),y(y0){ z=x*y; }
    INLINE mat_t(const vec_t<T> &x0,int,const vec_t<T> &z0):x(x0),z(z0){ y=z*x; }
    INLINE mat_t(int,const vec_t<T> &y0,const vec_t<T> &z0):y(y0),z(z0){ x=y*z; }
    template<typename T2>
    INLINE mat_t(const mat_t<T2> &a):x(a.x),y(a.y),z(a.z){}

    //converts a world vec_t to local
    INLINE vec_t<T> tolocal(const vec_t<T> &a) const{ return a%(*this); }
    //converts a local vec_t to world
    INLINE vec_t<T> toworld(const vec_t<T> &a) const{ return (*this)%a; }
    //converts a world mat_t to local
    INLINE mat_t<T> tolocal(const mat_t<T> &a) const{
        mat_t<T> sTa(a.x%(*this),a.y%(*this),a.z%(*this));
        return sTa%(*this);
    }
    //converts a local mat_t to world
    INLINE mat_t<T> toworld(const mat_t<T> &a) const{
        mat_t<T> sa((*this)%a);
        return mat_t<T>(sa.x*x.x+sa.y*y.x+sa.z*z.x,sa.x*x.y+sa.y*y.y+sa.z*z.y,sa.x*x.z+sa.y*y.z+sa.z*z.z);
    }
    INLINE void thetaphi(const vec_t<T> &a,T &theta,T &phi) const{
        vec_t<T> na(tolocal(a));
        phi=na.phi();
        theta=na.theta();
    }
    INLINE mat_t<T> &rotx(const T &b){
        T s=sin(b),c=cos(b);
        vec_t<T> ny(y*c+z*s);
        z=z*c-y*s;
        y=ny;
        return *this;
    }
    INLINE mat_t<T> &roty(const T &b){
        T s=sin(b),c=cos(b);
        vec_t<T> nz(z*c+x*s);
        x=x*c-z*s;
        z=nz;
        return *this;
    }
    INLINE mat_t<T> &rotz(const T &b){
        T s=sin(b),c=cos(b);
        vec_t<T> nx(x*c+y*s);
        y=y*c-x*s;
        x=nx;
        return *this;
    }
    INLINE T theta(const vec_t<T> &a) const{ return tolocal(a).theta(); }
    INLINE T phi(const vec_t<T> &a) const{
        return atan2(a%y,a%x);
    }

    INLINE mat_t<T> &operator =(const T &a){
        x=vec_t<T>(a,0,0);
        y=vec_t<T>(0,a,0);
        z=vec_t<T>(0,0,a);
        return *this;
    }
    template<typename T2>
    INLINE mat_t<T> &operator =(const mat_t<T2> &a){
        x=a.x;
        y=a.y;
        z=a.z;
        return *this;
    }
    INLINE mat_t<T> &operator +=(const T &a){
        x.x+=a;
        y.y+=a;
        z.z+=a;
        return *this;
    }
    INLINE mat_t<T> &operator -=(const T &a){
        x.x-=a;
        y.y-=a;
        z.z-=a;
        return *this;
    }
    INLINE mat_t<T> &operator *=(const T &a){
        x*=a;
        y*=a;
        z*=a;
        return *this;
    }
    INLINE mat_t<T> &operator /=(const T &a){
        T reca=1/a;
        x*=reca;
        y*=reca;
        z*=reca;
        return *this;
    }
    INLINE mat_t<T> &operator +=(const mat_t<T> &b){
        x+=b.x;
        y+=b.y;
        z+=b.z;
        return *this;
    }
    INLINE mat_t<T> &operator -=(const mat_t<T> &b){
        x-=b.x;
        y-=b.y;
        z-=b.z;
        return *this;
    }
    INLINE mat_t<T> &operator %=(const mat_t<T> &b){
        *this=*this%b;
        return *this;
    }
    INLINE T tr() const{
        return x.x+y.y+z.z;
    }
    INLINE T det() const{
        return x*y%z;
    }
    //square of Frobenius norm
    INLINE T normsqr() const{
        return x.normsqr()+y.normsqr()+z.normsqr();
    }
    //Frobenius norm
    INLINE T norm() const{
        return hypot_t<T>::length(x.x,x.y,x.z,y.x,y.y,y.z,z.x,z.y,z.z);
    }
    INLINE mat_t<T> adjoint() const{
        return mat_t<T>(
            vec_t<T>(y.y*z.z-y.z*z.y,z.y*x.z-x.y*z.z,x.y*y.z-x.z*y.y),
            vec_t<T>(y.z*z.x-y.x*z.z,z.z*x.x-x.z*z.x,x.z*y.x-x.x*y.z),
            vec_t<T>(y.x*z.y-y.y*z.x,z.x*x.y-x.x*z.y,x.x*y.y-x.y*y.x)
        );
    }
    INLINE mat_t<T> inverse() const{
        mat_t<T> madj=adjoint();
        T rdet=1/(madj.x.z*z.x+madj.y.z*z.y+madj.z.z*z.z);
        return madj*rdet;
    }
    INLINE mat_t<T> transpose() const{
        mat_t<T> trans(
            vec_t<T>(x.x,y.x,z.x),
            vec_t<T>(x.y,y.y,z.y),
            vec_t<T>(x.z,y.z,z.z)
        );
        return trans;
    }

    friend INLINE mat_t<T> operator +(const mat_t<T> &a){
        return mat_t<T>(+a.x,+a.y,+a.z);
    }
    friend INLINE mat_t<T> operator -(const mat_t<T> &a){
        return mat_t<T>(-a.x,-a.y,-a.z);
    }
    friend INLINE mat_t<T> operator +(const mat_t<T> &a,const mat_t<T> &b){
        return mat_t<T>(a.x+b.x,a.y+b.y,a.z+b.z);
    }
    friend INLINE mat_t<T> operator -(const mat_t<T> &a,const mat_t<T> &b){
        return mat_t<T>(a.x-b.x,a.y-b.y,a.z-b.z);
    }
    friend INLINE vec_t<T> operator %(const mat_t<T> &a,const vec_t<T> &b){
        return a.x*b.x+a.y*b.y+a.z*b.z;
    }
    friend INLINE vec_t<T> operator %(const vec_t<T> &a,const mat_t<T> &b){
        return vec_t<T>(a%b.x,a%b.y,a%b.z);
    }
    friend INLINE mat_t<T> operator %(const mat_t<T> &a,const mat_t<T> &b){
        return mat_t<T>(a%b.x,a%b.y,a%b.z);
    }

    friend INLINE mat_t<T> operator +(const T &a,const mat_t<T> &b){
        mat_t<T> result(b);
        result+=a;
        return result;
    }
    friend INLINE mat_t<T> operator +(const mat_t<T> &b,const T &a){
        mat_t<T> result(b);
        result+=a;
        return result;
    }
    friend INLINE mat_t<T> operator -(const T &a,const mat_t<T> &b){
        mat_t<T> result(-b);
        result+=a;
        return result;
    }
    friend INLINE mat_t<T> operator -(const mat_t<T> &b,const T &a){
        mat_t<T> result(b);
        result-=a;
        return result;
    }
    friend INLINE mat_t<T> operator *(const T &a,const mat_t<T> &b){
        return mat_t<T>(a*b.x,a*b.y,a*b.z);
    }
    friend INLINE mat_t<T> operator *(const mat_t<T> &b,const T &a){
        return mat_t<T>(a*b.x,a*b.y,a*b.z);
    }
    friend INLINE mat_t<T> operator /(const T &a,const mat_t<T> &b){
        mat_t<T> madj=b.adjoint();
        T rdet=a/(madj.x.z*b.z.x+madj.y.z*b.z.y+madj.z.z*b.z.z);
        return madj*rdet;
    }
    friend INLINE mat_t<T> operator /(const mat_t<T> &b,const T &a){
        T reca=1/a;
        return mat_t<T>(b.x*reca,b.y*reca,b.z*reca);
    }
};

template<typename T>
INLINE mat_t<T> vec_t<T>::rotation_matrix(const T &t) const{
    T _x=x,_y=y,_z=z;
    T wt=hypot_t<T>::normalize(_x,_y,_z)*t;
    T st=sin(wt),st222=sin(wt/2);
    st222*=2*st222;
    mat_t<T> dr(
        vec_t<T>(_x*_x-1,_x*_y,_x*_z),
        vec_t<T>(_x*_y,_y*_y-1,_y*_z),
        vec_t<T>(_x*_z,_y*_z,_z*_z-1)
    );
    dr*=st222;
    dr.x.y+=_z*st;
    dr.y.x-=_z*st;
    dr.x.z-=_y*st;
    dr.z.x+=_y*st;
    dr.y.z+=_x*st;
    dr.z.y-=_x*st;
    return dr;
}

template<typename T>
class quat_t{
public:
    T x,y,z,w;

    INLINE quat_t(){}
    INLINE quat_t(const T &w):x(0),y(0),z(0),w(w){}
    INLINE quat_t(const T &x,const T &y,const T &z,const T &w=0):x(x),y(y),z(z),w(w){}
    template<typename T2>
    INLINE quat_t(const vec_t<T2> &a,const T &w=0):x(a.x),y(a.y),z(a.z),w(w){}
    template<typename T2>
    INLINE quat_t(const quat_t<T2> &a):x(a.x),y(a.y),z(a.z),w(a.w){}
    template<typename T2>
    explicit INLINE operator vec_t<T2>() const{ return vec_t<T2>(T2(x),T2(y),T2(z)); }
    
    INLINE quat_t(const mat_t<T> &m){
        do{
            T mx=m.y.y+m.z.z;
            if(mx<=0&&m.y.y<=m.x.x&&m.z.z<=m.x.x){
                x=sqrt(m.x.x-mx+T(1));
                w=T(1)/x;
                y=w*(m.x.y+m.y.x);
                z=w*(m.x.z+m.z.x);
                w=w*(m.y.z-m.z.y);
                break;
            }
            T my=m.z.z+m.x.x;
            if(my<=0&&m.z.z<=m.y.y&&m.x.x<=m.y.y){
                y=sqrt(m.y.y-my+T(1));
                w=T(1)/y;
                z=w*(m.y.z+m.z.y);
                x=w*(m.y.x+m.x.y);
                w=w*(m.z.x-m.x.z);
                break;
            }
            T mz=m.x.x+m.y.y;
            if(mz<=0&&m.x.x<=m.z.z&&m.y.y<=m.z.z){
                z=sqrt(m.z.z-mz+T(1));
                w=T(1)/z;
                x=w*(m.z.x+m.x.z);
                y=w*(m.z.y+m.y.z);
                w=w*(m.x.y-m.y.x);
            }
            else{
                w=sqrt(m.z.z+mz+T(1));
                z=T(1)/w;
                x=z*(m.y.z-m.z.y);
                y=z*(m.z.x-m.x.z);
                z=z*(m.x.y-m.y.x);
            }
        } while(0);
        normalize();
    }
    explicit INLINE operator mat_t<T>() const{
        T xx=x*x,yy=y*y,zz=z*z,ww=w*w;
        T q1=w*x,q2=w*y,q3=w*z;
        T qx=y*z,qy=z*x,qz=x*y;
        q1+=q1;q2+=q2;q3+=q3;
        qx+=qx;qy+=qy;qz+=qz;
        return mat_t<T>(
            vec_t<T>(ww+xx-(yy+zz),qz+q3,qy-q2),
            vec_t<T>(qz-q3,ww+yy-(xx+zz),qx+q1),
            vec_t<T>(qy+q2,qx-q1,ww+zz-(xx+yy))
            );
    }

    INLINE T norm() const{
        return hypot_t<T>::length(x,y,z,w);
    }
    INLINE T normsqr() const{
        return x*x+y*y+z*z+w*w;
    }
    //normalize *this by norm()
    //  *this = (0,0,0,1) for (0,0,0,0)
    //  return norm()
    INLINE T normalize(){
        return hypot_t<T>::normalize(x,y,z,w);
    }
    INLINE quat_t<T> unit() const{
        quat_t<T> result(x,y,z,w);
        result.normalize();
        return result;
    }

    INLINE quat_t<T> &operator =(const T &a){
        x=0;
        y=0;
        z=0;
        w=a;
        return *this;
    }
    template<typename T2>
    INLINE quat_t<T> &operator =(const vec_t<T2> &a){
        x=a.x;
        y=a.y;
        z=a.z;
        w=0;
        return *this;
    }
    template<typename T2>
    INLINE quat_t<T> &operator =(const quat_t<T2> &a){
        x=a.x;
        y=a.y;
        z=a.z;
        w=a.w;
        return *this;
    }
    INLINE quat_t<T> &operator +=(const T &a){
        w+=a;
        return *this;
    }
    INLINE quat_t<T> &operator -=(const T &a){
        w-=a;
        return *this;
    }
    INLINE quat_t<T> &operator *=(const T &a){
        x*=a;
        y*=a;
        z*=a;
        w*=a;
        return *this;
    }
    INLINE quat_t<T> &operator /=(const T &a){
        T reca=1/a;
        x*=reca;
        y*=reca;
        z*=reca;
        w*=reca;
        return *this;
    }
    INLINE quat_t<T> &operator +=(const vec_t<T> &b){
        x+=b.x;
        y+=b.y;
        z+=b.z;
        return *this;
    }
    INLINE quat_t<T> &operator -=(const vec_t<T> &b){
        x-=b.x;
        y-=b.y;
        z-=b.z;
        return *this;
    }
    INLINE quat_t<T> &operator *=(const vec_t<T> &b){
        *this=*this*b;
        return *this;
    }
    INLINE quat_t<T> &operator /=(const vec_t<T> &b){
        *this=*this/b;
        return *this;
    }
    INLINE quat_t<T> &operator +=(const quat_t<T> &b){
        x+=b.x;
        y+=b.y;
        z+=b.z;
        w+=b.w;
        return *this;
    }
    INLINE quat_t<T> &operator -=(const quat_t<T> &b){
        x-=b.x;
        y-=b.y;
        z-=b.z;
        w-=b.w;
        return *this;
    }
    INLINE quat_t<T> &operator *=(const quat_t<T> &b){
        *this=*this*b;
        return *this;
    }
    INLINE quat_t<T> &operator /=(const quat_t<T> &b){
        *this=*this/b;
        return *this;
    }
    INLINE quat_t<T> inverse() const{
        const T r=1/normsqr();
        return quat_t<T>(-x,-y,-z,w)*r;
    }
    INLINE quat_t<T> conjugate() const{
        return quat_t<T>(-x,-y,-z,w);
    }

    friend INLINE quat_t<T> operator +(const quat_t<T> &a){
        return quat_t<T>(+a.x,+a.y,+a.z,+a.w);
    }
    friend INLINE quat_t<T> operator -(const quat_t<T> &a){
        return quat_t<T>(-a.x,-a.y,-a.z,-a.w);
    }
    friend INLINE quat_t<T> operator +(const quat_t<T> &a,const quat_t<T> &b){
        return quat_t<T>(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w);
    }
    friend INLINE quat_t<T> operator -(const quat_t<T> &a,const quat_t<T> &b){
        return quat_t<T>(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w);
    }
    friend INLINE quat_t<T> operator *(const quat_t<T> &a,const quat_t<T> &b){
        return quat_t<T>(
          a.w*b.x+(a.x*b.w+a.y*b.z-a.z*b.y),
          a.w*b.y+(a.y*b.w+a.z*b.x-a.x*b.z),
          a.w*b.z+(a.z*b.w+a.x*b.y-a.y*b.x),
          a.w*b.w-(a.x*b.x+a.y*b.y+a.z*b.z));
    }
    friend INLINE quat_t<T> operator /(const quat_t<T> &a,const quat_t<T> &b){
        return a*b.inverse();
    }
    friend INLINE T operator %(const quat_t<T> &a,const quat_t<T> &b){
        return a.x*b.x+a.y*b.y+a.z*b.z+a.w*b.w;
    }

    friend INLINE quat_t<T> operator +(const T &a,const quat_t<T> &b){
        quat_t<T> result(b);
        result+=a;
        return result;
    }
    friend INLINE quat_t<T> operator +(const quat_t<T> &b,const T &a){
        quat_t<T> result(b);
        result+=a;
        return result;
    }
    friend INLINE quat_t<T> operator -(const T &a,const quat_t<T> &b){
        quat_t<T> result(-b);
        result+=a;
        return result;
    }
    friend INLINE quat_t<T> operator -(const quat_t<T> &b,const T &a){
        quat_t<T> result(b);
        result-=a;
        return result;
    }
    friend INLINE quat_t<T> operator *(const T &a,const quat_t<T> &b){
        return quat_t<T>(a*b.x,a*b.y,a*b.z,a*b.w);
    }
    friend INLINE quat_t<T> operator *(const quat_t<T> &b,const T &a){
        return quat_t<T>(a*b.x,a*b.y,a*b.z,a*b.w);
    }
    friend INLINE quat_t<T> operator /(const T &a,const quat_t<T> &b){
        return (a/b.normsqr())*b.conjugate();
    }
    friend INLINE quat_t<T> operator /(const quat_t<T> &b,const T &a){
        T reca=1/a;
        return quat_t<T>(b.x*reca,b.y*reca,b.z*reca,b.w*reca);
    }

    friend INLINE quat_t<T> operator +(const vec_t<T> &a,const quat_t<T> &b){
        quat_t<T> result(b);
        result+=a;
        return result;
    }
    friend INLINE quat_t<T> operator +(const quat_t<T> &b,const vec_t<T> &a){
        quat_t<T> result(b);
        result+=a;
        return result;
    }
    friend INLINE quat_t<T> operator -(const vec_t<T> &a,const quat_t<T> &b){
        quat_t<T> result(-b);
        result+=a;
        return result;
    }
    friend INLINE quat_t<T> operator -(const quat_t<T> &b,const vec_t<T> &a){
        quat_t<T> result(b);
        result-=a;
        return result;
    }
    friend INLINE quat_t<T> operator *(const vec_t<T> &a,const quat_t<T> &b){
        return quat_t<T>(
            a.x*b.w+a.y*b.z-a.z*b.y,
            a.y*b.w+a.z*b.x-a.x*b.z,
            a.z*b.w+a.x*b.y-a.y*b.x,
          -(a.x*b.x+a.y*b.y+a.z*b.z));
    }
    friend INLINE quat_t<T> operator *(const quat_t<T> &b,const vec_t<T> &a){
        return quat_t<T>(
            b.w*a.x+b.y*a.z-b.z*a.y,
            b.w*a.y+b.z*a.x-b.x*a.z,
            b.w*a.z+b.x*a.y-b.y*a.x,
          -(b.x*a.x+b.y*a.y+b.z*a.z));
    }
    friend INLINE quat_t<T> operator /(const vec_t<T> &a,const quat_t<T> &b){
        T recb=1/b.normsqr();
        return (a*recb)*b.conjugate();
    }
    friend INLINE quat_t<T> operator /(const quat_t<T> &b,const vec_t<T> &a){
        T reca=1/(-a.normsqr());
        return b*(a*reca);
    }

    friend INLINE quat_t<T> exp(const quat_t<T> &q){
        T qv=hypot_t<T>::length(q.x,q.y,q.z);
        T ew=exp(q.w);
        T sincqv=(qv==0?1:sin(qv)/qv)*ew;
        return quat_t<T>(q.x*sincqv,q.y*sincqv,q.z*sincqv,cos(qv)*ew);
    }
    friend INLINE quat_t<T> log(const quat_t<T> &q){
        vec_t<T> v=vec_t<T>(q);
        T qv=hypot_t<T>::normalize(v.z,v.y,v.x);
        T qn=hypot_t<T>::length(qv,q.w);
        return quat_t<T>(v*atan2(qv,q.w),log(qn));
    }
};
