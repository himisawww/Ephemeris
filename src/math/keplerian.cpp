#include"keplerian.h"

using Constants::pi;
using Constants::epsilon;

//12(2^-53)^(1/4)
constexpr double parabolic_threshold=0x.5p-8;

double keplerian::kep(double x){
    double s=x+1;
    double t=x-1;
    if(std::abs(t)<0.5){
#define SUM_SERIES             \
    int i=1;                   \
    double df=t/35;            \
    double f=4.0/5+df;         \
    do{                        \
        ++i;                   \
        df*=-t*i/(2*i+5);      \
        f+=df;                 \
    } while(std::abs(df)>=epsilon)
        SUM_SERIES;
        return (4.0/3+t*f)/(s*s*s);
    }
    if(x<128/epsilon){
        double f=1/(s*t);
        return f*(1-(t>0?acosh(x)*sqrt(f):acos(x)*sqrt(-f)));
    }
    return 1/x/x;
}
double keplerian::dkep(double x){
    double s=x+1;
    double t=x-1;
    if(std::abs(t)<0.5){
        SUM_SERIES;
#undef SUM_SERIES
        s*=s;
        return (t-3*x*f)/(s*s);
    }
    if(x<128/epsilon){
        double f=1/(s*t);
        return f/s-3*x*f*f*(1-(t>0?acosh(x)*sqrt(f):acos(x)*sqrt(-f)));
    }
    return -2/x/(x*x);
}

keplerian::keplerian(const vec &r,const vec &v){
    //precision of j is crutial for further calculation...
    j=vec(mpvec(r)*mpvec(v));
    double rr=r%r,r0=sqrt(rr);
    //adjust j-direction if j not perpendicular to r due to numerical errors
    j-=(j%r/rr)*r;
    double jj=j%j;
    double jjeps=epsilon*epsilon*r0;
    if(jj<jjeps){
        //if jj==0, use finite-pseudo-j to avoid singularity
        j=r.asc_node();
        j*=sqrt(jjeps);
        jj=j%j;
    }

    vec e=v*j-r/r0;
    double ee=e%e,e0=sqrt(ee);
    double rjj=1/jj,srjj=sqrt(rjj);
    vec jx=j.asc_node(),jy=j*jx*srjj;
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

    double c=q<0?e0-q*x:sqrt(1+q*yy);
    if(c<0){
        double sq=sqrt(-q);
        m=(e0*y-atan2(y*sq,c)/sq)/q;
    }
    else{
        m=y/(1+e0)+y*yy*kep(c);
    }

    //adjust definition of q to avoid numeric error at e=0
    q/=e0+1;
}
void keplerian::rv(double t,vec &r,vec &v) const{
    double jj=j%j,rjj=1/jj,srjj=sqrt(rjj);
    vec jx=j.asc_node(),jy=j*jx*srjj;

    //current modified mean anomaly
    double mm=m+srjj*rjj*t;
    
    //restore e and q=e^2-1 from modified q-parameter
    double e0=q;
    double q=e0*(e0+2);
    e0+=1;

    double sq3=std::abs(q),sq=sqrt(sq3);
    sq3*=sq;
    double mp=mm*sq3;
    if(q<0&&std::abs(mp)>pi){
        //period folding for elliptic orbit
        mp=angle_reduce(mp);
        mm=mp/sq3;
    }

    bool msign=mm<0;
    if(msign){
        mm=-mm;
        mp=-mp;
    }

    double x=NAN,y=3*mm;
    y=cbrt(y+sqrt(1+y*y));
    y-=1/y;

    if(std::abs((3+y*y)*q)<parabolic_threshold){
        //parabolic orbit

        //solve y
        double oldx=NAN,oldy=NAN,yy,y4;
        double olddiff=HUGE_VAL,diff;
        do{
            yy=y*y;
            y4=yy*yy;
            x=sqrt(1+q*yy);
            diff=y/(1+e0)+y*yy*kep(x)-mm;
            if(!(std::abs(diff)<std::abs(olddiff)))break;
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

        //solve e0*sh(ea)-ea==mp
        double ea=mp<pi?cbrt(6*mp):asinh(mp/e0);
        double olddiff=HUGE_VAL,diff;
        do{
            double c=cosh(ea);
            double s=sinh(ea);
            diff=s*(q/(1+e0)+s*s*kep(c))-mp;
            if(!(std::abs(diff)<std::abs(olddiff)))break;
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

        //solve ea-e0*sin(ea)==mp
        double ea=cbrt(6*mp);
        double olddiff=HUGE_VAL,diff;
        do{
            double c=cos(ea);
            double s=sin(ea);
            diff=(c<0?ea-e0*s:s*(-q/(1+e0)+s*s*kep(c)))-mp;
            if(!(std::abs(diff)<std::abs(olddiff)))break;
            x=c;
            y=s;
            ea-=diff/(1-e0*c);
            olddiff=diff;
        } while(1);

        y=y/sq;
        if(x+q<0)x=(e0-x)/q;
        else x=(1-y*y)/(e0+x);
    }

    if(msign)y=-y;

    double vv=1/sqrt(jj*(x*x+y*y));
    //adjust vy for higher precision (when e->+inf)
    double vx=-y,vy=q<0?e0-q*x:sqrt(1+q*y*y);
    vx*=vv;
    vy*=vv;
    x*=jj;
    y*=jj;

    double ce=cos(earg),se=sin(earg);
    r=jx*( x*ce- y*se)+jy*( x*se+ y*ce);
    v=jx*(vx*ce-vy*se)+jy*(vx*se+vy*ce);
}
