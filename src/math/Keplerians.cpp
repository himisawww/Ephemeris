#include"Keplerians.h"
#define USE_NEW_ALGORITHM
//12(2^-53)^(1/4)
static const double parabolic_threshold=0.00123;
static const double epsilon=1e-16;

double kep(double x){
    double s=x+1;
    double t=x-1;

    if(abs(t)>0.5){
        double f=1/(s*t);
        f-=(t>0?acosh(x)*(f*sqrt(f)):acos(x)*(f*sqrt(-f)));
        return f;
    }

#define SUM_SERIES             \
    int i=1;                   \
    double df=t/35;            \
    double f=4.0/5+df;         \
    do{                        \
        ++i;                   \
        df*=-t*i/(2*i+5);      \
        f+=df;                 \
    } while(abs(df)>=2e-16)
    SUM_SERIES;

    return (4.0/3+t*f)/(s*s*s);
}
double dkep(double x){
    double s=x+1;
    double t=x-1;

    if(abs(t)>0.5){
        double f=1/(s*t);
        f=f/s-3*x*f*f*(1-(t>0?acosh(x)*sqrt(f):acos(x)*sqrt(-f)));
        return f;
    }

    SUM_SERIES;
#undef SUM_SERIES
    
    s*=s;
    return (t-3*x*f)/(s*s);
}

ephem_orb::ephem_orb(double t,const vec &r,const vec &v):t(t){
    j=r*v;
    double rr=r%r,r0=sqrt(rr);
    //adjust j-direction if j not perpendicular to r due to numerical errors
    j-=(j%r/rr)*r;
    double jj=j%j;
    double jjeps=epsilon*epsilon*r0;
    if(jj<jjeps){
        //if jj==0, use finite-pseudo-j to avoid singularity
        j=r.perpunit();
        j*=sqrt(jjeps);
        jj=j%j;
    }

    vec e=v*j-r/r0;
    double ee=e%e,e0=sqrt(ee);
    double rjj=1/jj,srjj=sqrt(rjj);
    vec jx=j.perpunit(),jy=j*jx*srjj;
    double ex=e%jx,ey=e%jy;

    earg=atan2(ey,ex);

    double rx=(r%jx)*rjj,ry=(r%jy)*rjj;
    double ce=cos(earg),se=sin(earg);
    r0*=rjj;
    double x=ce*rx+se*ry,xx=x*x;
    double y,yy;

    //adjust y for higher precision (when e->1)
    if(e0>0.5)y=v%r/e0*srjj;
    else y=ce*ry-se*rx;
    yy=y*y;

    //adjust q for higher precision (when e->1)
    if(e0<2&&ee*xx>1+yy)q=(1+yy-2*r0)/xx;
    else q=ee-1;
#ifndef USE_NEW_ALGORITHM
    if(abs((2-x)*q)<parabolic_threshold/2){
        //parabolic orbit
        double y4=yy*yy;
        m=y*((3+yy)/6+q*(q*(7+5*y4*yy)/112-(5+3*y4)/40));
    }
    else if(q>0){
        //hyperbolic orbit
        double sq=sqrt(q);
        m=(e0*y-asinh(y*sq)/sq)/q;
    }
    else{
        //elliptic orbit
        double sq=sqrt(-q);
        m=(e0*y-atan2(y*sq,e0-q*x)/sq)/q;
    }
#else
    double c=q<0?e0-q*x:sqrt(1+q*yy);
    if(c<0){
        double sq=sqrt(-q);
        m=(e0*y-atan2(y*sq,c)/sq)/q;
    }
    else{
        m=y/(1+e0)+y*yy*kep(c);
    }
#endif
    //adjust definition of q to avoid numeric error at e=0
    q/=e0+1;
}
void ephem_orb::rv(double t,vec &r,vec &v) const{
    double jj=j%j,rjj=1/jj,srjj=sqrt(rjj);
    vec jx=j.perpunit(),jy=j*jx*srjj;

    //current modified mean anomaly
    double mp=m+(t-this->t)*srjj*rjj;
    
    //restore e and q=e^2-1 from modified q-parameter
    double e0=q;
    double q=e0*(e0+2);
    e0+=1;

    double sq3=abs(q),sq=sqrt(sq3);
    sq3*=sq;
    if(q<0){
        //period folding for elliptic orbit
        double p=2*pi/sq3;
        mp-=p*round(mp/p);
    }

    bool msign=mp<0;
    if(msign)mp=-mp;

    double x=NAN,y=NAN;

#ifndef USE_NEW_ALGORITHM
    y=3*mp;
    y=pow(y+sqrt(1+y*y),1/3.);
    y-=1/y;

    if(abs((3+y*y)*q)<parabolic_threshold){
        //parabolic orbit

        //solve y
        double oldy,yy,y4;
        double olddiff=HUGE_VAL,diff;
        do{
            yy=y*y;
            y4=yy*yy;
            diff=y*((3+yy)/6+q*(q*(7+5*y4*yy)/112-(5+3*y4)/40))-mp;
            if(!(abs(diff)<abs(olddiff)))break;
            oldy=y;
            y-=diff/((1+yy)/2+q*(q-2+(5*q*yy-6)*y4)/16);
            olddiff=diff;
        } while(1);

        y=oldy;
        yy=y*y;
        x=(1-yy)/(e0+sqrt(1+q*yy));
    }
    else if(q>0){
        //hyperbolic orbit
        mp*=sq3;

        //solve e0*sh(ea)-ea==mp
        double ea=mp<pi?pow(6*mp,1/3.):asinh(mp/e0);
        double olddiff=HUGE_VAL,diff;
        do{
            double shea=sinh(ea);
            diff=e0*shea-ea-mp;
            if(!(abs(diff)<abs(olddiff)))break;
            x=cosh(ea);
            y=shea;
            ea-=diff/(e0*x-1);
            olddiff=diff;
        } while(1);

        x=(e0-x)/q;
        y=y/sq;
    }
    else{
        //elliptic orbit
        mp*=sq3;

        //solve ea-e0*sin(ea)==mp
        double ea=pow(6*mp,1/3.);
        double olddiff=HUGE_VAL,diff;
        do{
            double sinea=sin(ea);
            diff=ea-e0*sinea-mp;
            if(!(abs(diff)<abs(olddiff)))break;
            x=cos(ea);
            y=sinea;
            ea-=diff/(1-e0*x);
            olddiff=diff;
        } while(1);

        x=(e0-x)/q;
        y=y/sq;
    }
#else
    y=3*mp;
    y=pow(y+sqrt(1+y*y),1/3.);
    y-=1/y;

    if(abs((3+y*y)*q)<parabolic_threshold){
        //parabolic orbit

        //solve y
        double oldx=NAN,oldy=NAN,yy,y4;
        double olddiff=HUGE_VAL,diff;
        do{
            yy=y*y;
            y4=yy*yy;
            x=sqrt(1+q*yy);
            diff=y/(1+e0)+y*yy*kep(x)-mp;
            if(!(abs(diff)<abs(olddiff)))break;
            oldx=x;
            oldy=y;
            y-=diff/((1+yy)/2+q*(q-2+(5*q*yy-6)*y4)/16);
            olddiff=diff;
        } while(1);

        y=oldy;
        x=(1-y*y)/(e0+oldx);
    }
    else if(q>0){
        //hyperbolic orbit
        mp*=sq3;

        //solve e0*sh(ea)-ea==mp
        double ea=mp<pi?pow(6*mp,1/3.):asinh(mp/e0);
        double olddiff=HUGE_VAL,diff;
        do{
            double c=cosh(ea);
            double s=sinh(ea);
            diff=s*(q/(1+e0)+s*s*kep(c))-mp;
            if(!(abs(diff)<abs(olddiff)))break;
            x=c;
            y=s;
            ea-=diff/(e0*c-1);
            olddiff=diff;
        } while(1);

        y=y/sq;
        x=(1-y*y)/(e0+x);
    }
    else{
        //elliptic orbit
        mp*=sq3;

        //solve ea-e0*sin(ea)==mp
        double ea=pow(6*mp,1/3.);
        double olddiff=HUGE_VAL,diff;
        do{
            double c=cos(ea);
            double s=sin(ea);
            diff=(c<0?ea-e0*s:s*(-q/(1+e0)+s*s*kep(c)))-mp;
            if(!(abs(diff)<abs(olddiff)))break;
            x=c;
            y=s;
            ea-=diff/(1-e0*c);
            olddiff=diff;
        } while(1);

        y=y/sq;
        if(x+q<0)x=(e0-x)/q;
        else x=(1-y*y)/(e0+x);
    }
#endif
    if(msign)y=-y;

    double rr=x*x+y*y;
    double vx=-y,vy=e0-q*x,vv=1/sqrt(jj*rr);
    //adjust vy for higher precision (when e->+inf)
    if(q>0)vy=sqrt(1+q*y*y);
    vx*=vv;
    vy*=vv;
    x*=jj;
    y*=jj;

    double ce=cos(earg),se=sin(earg);
    r=jx*( x*ce- y*se)+jy*( x*se+ y*ce);
    v=jx*(vx*ce-vy*se)+jy*(vx*se+vy*ce);
}

ephem_rot::ephem_rot(double t,const mat &s,const vec &_w):t(t),w(_w){
    double ww=w%w;
    if(ww==0){
        //if ww==0, use finite-pseudo-w to avoid singularity
        w=s.z*epsilon;
        ww=w%w;
    }

    double sww=sqrt(ww),srww=1/sww;
    mat wm(w.perpunit(),0,w*srww);
    vec wz=wm.z;

    s.thetaphi(wz,ptheta,pphi);
    wm.roty(-ptheta).rotz(-pphi);

    vec za,zb;
    if(abs(ptheta-pi/2)<pi/4){
        za=wz*wm.z;
        zb=wz*s.z;
    }
    else{
        za=wz*wm.x;
        zb=wz*s.x;
    }
    angle=atan2(za*zb%wz,za%zb);
}
void ephem_rot::sw(double t,mat &s,vec &w) const{
    w=this->w;
    double sww=sqrt(w%w),srww=1/sww;
    s=mat(w.perpunit(),0,w*srww);
    vec wz=s.z;
    s.roty(-ptheta).rotz(-pphi);

    //current rotation angle
    double theta=angle+(t-this->t)*sww;

    double cth=cos(theta),sth=sin(theta);
    s.x+=(s.x-wz%s.x*wz)*(cth-1)+wz*s.x*sth;
    s.y+=(s.y-wz%s.y*wz)*(cth-1)+wz*s.y*sth;
    s.z+=(s.z-wz%s.z*wz)*(cth-1)+wz*s.z*sth;
}
