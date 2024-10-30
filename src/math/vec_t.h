#pragma once
#include<cmath>//for trigs
#ifndef INLINE
#define INLINE inline
#endif

template<typename T>
class vec_t{
public:
    T x,y,z;

    INLINE vec_t(){}
    INLINE vec_t(const T &x):x(x),y(x),z(x){}
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
        return atan2(sqrt(x*x+y*y),z);
    }
    INLINE T phi() const{
        return atan2(y,x);
    }
    INLINE T norm() const{
        return sqrt(x*x+y*y+z*z);
    }
    INLINE T normsqr() const{
        return x*x+y*y+z*z;
    }
    //normalize *this by norm()
    //  *this = (0,0,1) for (0,0,0)
    //  return norm()
    INLINE T normalize(){
        const T r=norm();
        if(r==0){
            x=0;
            y=0;
            z=1;
        }
        else{
            T rr=1/r;
            x*=rr;
            y*=rr;
            z*=rr;
        }
        return r;
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
    INLINE vec_t<T> &operator +=(const T &a){
        x+=a;
        y+=a;
        z+=a;
        return *this;
    }
    INLINE vec_t<T> &operator -=(const T &a){
        x-=a;
        y-=a;
        z-=a;
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
        T r=x*x+y*y;
        if(r==0)p.x=1;
        else{
            r=1/sqrt(r);
            p.x*=r;
            p.y*=r;
        }
        return p;
    }
};

template<typename T>
INLINE vec_t<T> operator +(const vec_t<T> &a){
    return vec_t<T>(+a.x,+a.y,+a.z);
}
template<typename T>
INLINE vec_t<T> operator -(const vec_t<T> &a){
    return vec_t<T>(-a.x,-a.y,-a.z);
}
template<typename T>
INLINE vec_t<T> operator +(const vec_t<T> &a,const vec_t<T> &b){
    return vec_t<T>(a.x+b.x,a.y+b.y,a.z+b.z);
}
template<typename T>
INLINE vec_t<T> operator -(const vec_t<T> &a,const vec_t<T> &b){
    return vec_t<T>(a.x-b.x,a.y-b.y,a.z-b.z);
}
template<typename T>
INLINE vec_t<T> operator *(const vec_t<T> &a,const vec_t<T> &b){
    return vec_t<T>(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}
template<typename T>
INLINE T operator %(const vec_t<T> &a,const vec_t<T> &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

template<typename T>
INLINE vec_t<T> operator +(const T &a,const vec_t<T> &b){
    return vec_t<T>(a+b.x,a+b.y,a+b.z);
}
template<typename T>
INLINE vec_t<T> operator +(const vec_t<T> &b,const T &a){
    return vec_t<T>(a+b.x,a+b.y,a+b.z);
}
template<typename T>
INLINE vec_t<T> operator -(const T &a,const vec_t<T> &b){
    return vec_t<T>(a-b.x,a-b.y,a-b.z);
}
template<typename T>
INLINE vec_t<T> operator -(const vec_t<T> &b,const T &a){
    return vec_t<T>(b.x-a,b.y-a,b.z-a);
}
template<typename T>
INLINE vec_t<T> operator *(const T &a,const vec_t<T> &b){
    return vec_t<T>(a*b.x,a*b.y,a*b.z);
}
template<typename T>
INLINE vec_t<T> operator *(const vec_t<T> &b,const T &a){
    return vec_t<T>(a*b.x,a*b.y,a*b.z);
}
//template<typename T>
//INLINE vec_t<T> operator /(const T &a,const vec_t<T> &b){
//    return vec_t<T>(a/b.x,a/b.y,a/b.z);
//}
template<typename T>
INLINE vec_t<T> operator /(const vec_t<T> &b,const T &a){
    T reca=1/a;
    return vec_t<T>(b.x*reca,b.y*reca,b.z*reca);
}

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
    INLINE mat_t<T> inverse() const{
        T rdet=1/(x*y%z);
        mat_t<T> inv(
            vec_t<T>(y.y*z.z-y.z*z.y,z.y*x.z-x.y*z.z,x.y*y.z-x.z*y.y),
            vec_t<T>(y.z*z.x-y.x*z.z,z.z*x.x-x.z*z.x,x.z*y.x-x.x*y.z),
            vec_t<T>(y.x*z.y-y.y*z.x,z.x*x.y-x.x*z.y,x.x*y.y-x.y*y.x)
        );
        return inv*rdet;
    }
    INLINE mat_t<T> transpose() const{
        mat_t<T> trans(
            vec_t<T>(x.x,y.x,z.x),
            vec_t<T>(x.y,y.y,z.y),
            vec_t<T>(x.z,y.z,z.z)
        );
        return trans;
    }
};

template<typename T>
INLINE mat_t<T> operator +(const mat_t<T> &a){
    return mat_t<T>(+a.x,+a.y,+a.z);
}
template<typename T>
INLINE mat_t<T> operator -(const mat_t<T> &a){
    return mat_t<T>(-a.x,-a.y,-a.z);
}
template<typename T>
INLINE mat_t<T> operator +(const mat_t<T> &a,const mat_t<T> &b){
    return mat_t<T>(a.x+b.x,a.y+b.y,a.z+b.z);
}
template<typename T>
INLINE mat_t<T> operator -(const mat_t<T> &a,const mat_t<T> &b){
    return mat_t<T>(a.x-b.x,a.y-b.y,a.z-b.z);
}
template<typename T>
INLINE vec_t<T> operator %(const mat_t<T> &a,const vec_t<T> &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
template<typename T>
INLINE vec_t<T> operator %(const vec_t<T> &a,const mat_t<T> &b){
    return vec_t<T>(a%b.x,a%b.y,a%b.z);
}
template<typename T>
INLINE mat_t<T> operator %(const mat_t<T> &a,const mat_t<T> &b){
    return mat_t<T>(a%b.x,a%b.y,a%b.z);
}

template<typename T>
INLINE mat_t<T> operator +(const T &a,const mat_t<T> &b){
    mat_t<T> result(b);
    result+=a;
    return result;
}
template<typename T>
INLINE mat_t<T> operator +(const mat_t<T> &b,const T &a){
    mat_t<T> result(b);
    result+=a;
    return result;
}
template<typename T>
INLINE mat_t<T> operator -(const T &a,const mat_t<T> &b){
    mat_t<T> result(-b);
    result+=a;
    return result;
}
template<typename T>
INLINE mat_t<T> operator -(const mat_t<T> &b,const T &a){
    mat_t<T> result(b);
    result-=a;
    return result;
}
template<typename T>
INLINE mat_t<T> operator *(const T &a,const mat_t<T> &b){
    return mat_t<T>(a*b.x,a*b.y,a*b.z);
}
template<typename T>
INLINE mat_t<T> operator *(const mat_t<T> &b,const T &a){
    return mat_t<T>(a*b.x,a*b.y,a*b.z);
}
//template<typename T>
//INLINE mat_t<T> operator /(const T &a,const mat_t<T> &b){
//    return mat_t<T>(a/b.x,a/b.y,a/b.z);
//}
template<typename T>
INLINE mat_t<T> operator /(const mat_t<T> &b,const T &a){
    T reca=1/a;
    return mat_t<T>(b.x*reca,b.y*reca,b.z*reca);
}

//return matrix R s.t.
//(I+R)%x rotates vector x along axis w counterclockwise by angle w.norm()*t
template<typename T>
INLINE mat_t<T> rotation_matrix(const vec_t<T> &w,const T &t){
    T w2=w%w;
    if(w2==0)return 0;
    
    T w1=sqrt(w2),wt=w1*t,wr=1/w1;
    T st=sin(wt),st222=sin(wt/2);
    st222*=2*st222;
    T _x=w.x*wr,_y=w.y*wr,_z=w.z*wr;
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