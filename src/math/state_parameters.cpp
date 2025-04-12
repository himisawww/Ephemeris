#include"state_parameters.h"

using Constants::pi;
using Constants::epsilon;

rotational_param_t::rotational_param_t(const mat &s,const vec &_w):w(_w){
    vec wz=w;
    wz.normalize();
    mat wm(w.asc_node(),0,wz);

    s.thetaphi(wz,ptheta,pphi);
    wm.roty(-ptheta).rotz(-pphi);

    vec za,zb;
    if(std::abs(ptheta-pi/2)<pi/4){
        za=wz*wm.z;
        zb=wz*s.z;
    }
    else{
        za=wz*wm.x;
        zb=wz*s.x;
    }
    angle=atan2(za*zb%wz,za%zb);
}
void rotational_param_t::sw(double t,mat &s,vec &w) const{
    w=this->w;
    vec wz=w;
    double sww=wz.normalize();
    s=mat(w.asc_node(),0,wz);
    s.roty(-ptheta).rotz(-pphi);

    //current rotation angle
    double theta=angle+t*sww;

    double cth=cos(theta),sth=sin(theta);
    s.x+=(s.x-wz%s.x*wz)*(cth-1)+wz*s.x*sth;
    s.y+=(s.y-wz%s.y*wz)*(cth-1)+wz*s.y*sth;
    s.z+=(s.z-wz%s.z*wz)*(cth-1)+wz*s.z*sth;
}

rotational_param_t::rotational_param_t(const double *a){
    w.x=a[0];
    w.y=a[1];
    w.z=a[2];
    pphi=std::asin(a[3]/std::sqrt(1-a[4]*a[4]));
    ptheta=std::acos(a[4]);
    angle=a[5];
}

orbital_param_t::orbital_param_t(const keplerian &other):keplerian(other){
    el=earg+j.asc_node().phi();
    ex=(q+1)*std::cos(el);
    ey=(q+1)*std::sin(el);
    if(q>=0)
        ml=NAN;
    else{
        double s=-q*(q+2);
        s*=std::sqrt(s);
        ml=el+m*s;
        //ml=angle_reduce(ml);
    }
}
orbital_param_t::orbital_param_t(const double *p,bool circular){
    double &w=circular?q:ml;
    w=1;
    j.x=p[0];
    j.y=p[1];
    j.z=p[2];
    if(circular){
        ex=p[3];
        ey=p[4];
        ml=p[5];
    }
    else{
        q=p[3];
        el=p[4];
        m=p[5];
    }
    blend_finalize(circular);
}

bool orbital_param_t::interpolatable(bool circular,const orbital_param_t &k1,const orbital_param_t &k2){
    if(circular){
        double mq=checked_max(k1.q,k2.q);
        if(!(mq<0))return false;
        double rdq=std::abs(k1.q-k2.q)/mq;
        if(!(rdq<1))return false;
    }
    else{
        if(!(k2.m>k1.m))return false;
        double del=angle_reduce(k2.el-k1.el);
        if(!(del<1))return false;
    }
    return (k2.j-k1.j).normsqr()/checked_min(k1.j.normsqr(),k2.j.normsqr())<1;
}

void orbital_param_t::blend_initialize(bool circular){
    double &w=circular?q:ml;
    w=0;
    j=0;
    if(circular){
        ex=0;
        ey=0;
        ml=0;
        earg=NAN;//last_ml;
    }
    else{
        q=0;
        el=0;
        m=0;
        earg=NAN;//last_el;
    }
}
void orbital_param_t::blend_add(bool circular,const orbital_param_t &ki,double ws){
    double &w=circular?q:ml;
    if(circular){
        double &last_ml=earg;
        w+=ws;
        j+=ws*ki.j;
        ex+=ws*ki.ex;
        ey+=ws*ki.ey;
        double mli=ki.ml;
        if(last_ml==last_ml)mli=last_ml+angle_reduce(mli-last_ml);
        ml+=ws*mli;
        last_ml=mli;
    }
    else{
        double &last_el=earg;
        w+=ws;
        j+=ws*ki.j;
        q+=ws*ki.q;
        m+=ws*ki.m;
        double eli=ki.el;
        if(last_el==last_el)eli=last_el+angle_reduce(eli-last_el);
        el+=ws*eli;
        last_el=eli;
    }
}
void orbital_param_t::blend_finalize(bool circular){
    double &w=circular?q:ml;
    w=1/w;
    j*=w;
    if(circular){
        ex*=w;
        ey*=w;
        ml*=w;
        q=std::sqrt(ex*ex+ey*ey)-1;
        el=std::atan2(ey,ex);
        double s=-q*(q+2);
        s*=std::sqrt(s);
        m=(ml-el)/s;
        earg=el-j.asc_node().phi();
    }
    else{
        q*=w;
        m*=w;
        el*=w;
        earg=el-j.asc_node().phi();
        ex=(q+1)*std::cos(el);
        ey=(q+1)*std::sin(el);
        if(q>=0)
            ml=NAN;
        else{
            double s=-q*(q+2);
            s*=std::sqrt(s);
            ml=el+m*s;
        }
    }
}
