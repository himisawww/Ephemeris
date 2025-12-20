#include"definitions.h"
#include"utils/memio.h"
#include"physics/geopotential.h"
#include"utils/calctime.h"
#include"utils/logger.h"

#define TEST_MEMPATH "__gptest.gp"

using Constants::pi;

// test parameter
#define TEST_JVAL 0.1919810
#define TEST_CVAL 0.114
#define TEST_SVAL 0.0514
#define REFERENCE_RADIUS 1415.9265
#define REFERENCE_RADIUS_FACTOR 1.0123456
#define TEST_RADIUS 1828.1828
// test size
#define TEST_N 17
#define TEST_RADIUS_FACTOR (-2.2)
// max allowed relative error
#define TEST_EPSILON 1e-13

class geopotential_reference{
public:
	struct component{
		int_t n,m;
		fast_real c,s;
	};
	htl::vector<component> data;

	fast_mpvec sum(fast_real R,fast_mpvec r){
		fast_real r2=r%r;
		fast_real r1=sqrt(r2);
		fast_real rbar=R/r1;
		// longitude
		fast_real lambda=r.phi();
		// latitude
		fast_real phi=pi/2-r.theta();
		fast_real c=cos(phi),s=sin(phi);
		fast_mpvec g(0);
		for(const auto &gc:data){
			int_t n=gc.n,m=gc.m;
			fast_real rbarn=std::pow(rbar,n);
			fast_real csign=std::copysign(fast_real(1),c);
			if(m==0){
				// Legendre P_n(sin(phi))
				fast_real pn=std::legendre(n,s);
				// Legendre cos(phi)P'_n(sin(phi))
				fast_real cppn=csign*std::assoc_legendre(n,1,s);
				g+=rbarn*gc.c*fast_mpvec((n+1)*pn,0,-cppn);
			}
			else{
				fast_real cm=cos(m*lambda),sm=sin(m*lambda);
				fast_real tx=gc.c*cm+gc.s*sm,ty=-gc.c*sm+gc.s*cm;
				// Legendre sec(phi)P_n^m(sin(phi))
				fast_real sp;
				if(std::abs(c)>std::abs(s))
					sp=std::assoc_legendre(n,m,s)/c;
				else
					sp=csign*(
						std::assoc_legendre(n,m+1,s)
						+(n+m)*(n-m+1)*std::assoc_legendre(n,m-1,s)
						)/(2*m*s);
				// Legendre cos(phi)P'_n^m(sin(phi))
				fast_real cp=(m*s*sp-(n+m)*(n-m+1)*csign*std::assoc_legendre(n,m-1,s));
				g+=rbarn*fast_mpvec(-(n+1)*c*sp*tx,m*sp*ty,cp*tx);
			}
		}
		g/=r2;
		fast_mpvec x,y,z=fast_mpvec::from_theta_phi(-phi,lambda);
		x=r/r1;
		y=z*x;
		return g.x*x+g.y*y+g.z*z;
	}
};

int test_geopotential(){
	fast_real max_rerr(0);
	//test every components
	for(int_t n=2;n<=8;++n)for(int_t m=0;m<=n;++m)for(int_t cs=0;cs<=1;++cs){
		fast_real c,s;
		if(m==0){
			if(cs==1)continue;
			c=TEST_JVAL/n;
			s=0;
		}
		else{
			if(cs==0){
				c=TEST_CVAL/n;
				s=0;
			}
			else{
				c=0;
				s=TEST_SVAL/n;
			}
		}

		geopotential_reference gpref;
		gpref.data.push_back({n,m,c,s});
		MFILE gpmf;
		fprintf(&gpmf,"%lld\t%lld\n",n,m?n:0);
		fprintf(&gpmf,"%lld\t%lld\t%.16le\t%.16le\n",n,m,c,s);
		gpmf.publish(TEST_MEMPATH);
		const geopotential *gp=geopotential::load(TEST_MEMPATH,REFERENCE_RADIUS_FACTOR,2);
		if(!gp)return 2;
		fast_real rn=TEST_RADIUS;
		for(int_t ir=0;ir<TEST_N;++ir){
			fast_real max_norm(0);
			fast_real max_diff(0);
			rn*=TEST_RADIUS_FACTOR;
			for(int_t itheta=0;itheta<=TEST_N;++itheta)for(int_t iphi=0;iphi<=2*TEST_N;++iphi){
				fast_mpvec r=fast_mpvec::from_theta_phi(itheta*pi/TEST_N,iphi*pi/TEST_N);
				r*=rn;
				fast_mpvec gtest=gp->sum(REFERENCE_RADIUS,r);
				fast_mpvec gref=gpref.sum(REFERENCE_RADIUS*REFERENCE_RADIUS_FACTOR,r);
				checked_maximize(max_norm,gtest.norm());
				checked_maximize(max_diff,(gtest-gref).norm());
			}
			for(int_t iz=-1;iz<=1;iz+=2)for(int_t iphi=0;iphi<=2*TEST_N;++iphi){
				fast_real xynorm=std::sqrt(DBL_MIN)*(3*ir+5)/(2*ir+4);
				fast_real xyphi=iphi*pi/TEST_N;
				fast_mpvec r(xynorm*cos(xyphi),xynorm*sin(xyphi),iz*rn);
				fast_mpvec grz=gp->sum(REFERENCE_RADIUS,vec(0,0,r.z));
				do{
					fast_mpvec gsubr=gp->sum(REFERENCE_RADIUS,r);
					checked_maximize(max_diff,(gsubr-grz).norm());
					r.x/=2;
					r.y/=2;
				} while(r.x*r.x||r.y*r.y);
			}
			checked_maximize(max_rerr,max_diff/max_norm);
		}
		geopotential::unload(gp);

		if(!(max_rerr<TEST_EPSILON)){
			LogError(
				"\nMax Rel. Error %.16le Too Large at n=%lld, m=%lld",
				max_rerr,n,m);
			return 1;
		}
	}

	LogInfo("\n      Passed(%f), ",max_rerr/TEST_EPSILON);
    return 0;
}
