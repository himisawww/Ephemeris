#include"definitions.h"
#include"utils/memio.h"
#include"physics/geopotential.h"
#include"utils/calctime.h"

#define TEST_MEMPATH "__gptest.gp"

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

// c=cos(phi), s=sin(phi)
// where phi is latitude in surface frame of extended body

// Legendre P_n(sin(phi))
static fast_real pn(int_t n,fast_real s){
	fast_real s2=s*s;
	switch(n){
	case 2: return (3*s2-1)/2;
	case 3: return (5*s2-3)*s/2;
	case 4: return (3+s2*(-30+s2*35))/8;
	case 5: return (15+s2*(-70+s2*63))*s/8;
	case 6: return (-5+21*s2*(5+s2*(-15+11*s2)))/16;
	case 7: return (-35+s2*(315+s2*(-693+s2*429)))*s/16;
	case 8: return (35+s2*(-1260+s2*(6930+s2*(-12012+s2*6435))))/128;
	default: return NAN;
	}
}
// Legendre P'_n(sin(phi))
static fast_real ppn(int_t n,fast_real s){
	fast_real s2=s*s;
	switch(n){
	case 2: return 3*s;
	case 3: return (5*s2-1)*3/2;
	case 4: return (7*s2-3)*s*5/2;
	case 5: return (1+s2*(-14+s2*21))*15/8;
	case 6: return (5+s2*(-30+s2*33))*s*21/8;
	case 7: return (-5+s2*(135+s2*(-495+s2*429)))*7/16;
	case 8: return (-35+11*s2*(35+s2*(-91+s2*65)))*s*9/16;
	default: return NAN;
	}
}
// Legendre sec(phi)P_n^m(sin(phi))
static fast_real spnm(int_t n,int_t m,fast_real c,fast_real s){
	int_t degord=(m<1||m>n)?-1:n*10+m;
	fast_real c2=c*c,s2=s*s;
	switch(degord){
	case 21: return 3*s;
	case 22: return 3*c;
	case 31: return 3*(-1+s2*5)/2;
	case 32: return 15*s*c;
	case 33: return 15*c2;
	case 41: return 5*s*(-3+s2*7)/2;
	case 42: return -15*(c*(-3+c2*7)+s2*c*-21)/8;
	case 43: return 105*s*c2;
	case 44: return 105*c*c2;
	case 51: return 15*((15+c2*(-28+c2*21))+s2*((28+c2*-126)+s2*21))/64;
	case 52: return -105*s*(c*(-1+c2*3)+s2*c*-3)/4;
	case 53: return -105*((-5+c2*(4+c2*9))+s2*((-4+c2*-54)+s2*9))/16;
	case 54: return 945*s*c*c2;
	case 55: return 945*c2*c2;
	case 61: return 21*s*((50+c2*(-135+c2*165))+s2*((45+c2*-330)+s2*33))/128;
	case 62: return 105*(c*(10+c2*(-27+c2*33))+s2*(c*(81+c2*-330)+s2*c*165))/128;
	case 63: return -315*s*((-10+c2*(3+c2*55))+s2*((-1+c2*-110)+s2*11))/32;
	case 64: return -945*(c*(-10+c2*(15+c2*11))+s2*(c*(-45+c2*-110)+s2*c*55))/32;
	case 65: return 10395*s*c2*c2;
	case 66: return 10395*c*c2*c2;
	case 71: return 7*((350+c2*(-675+c2*(594+c2*-429)))+s2*((675+c2*(-3564+c2*6435))+s2*((594+c2*-6435)+s2*429)))/512;
	case 72: return 63*s*(c*(75+c2*(-264+c2*429))+s2*(c*(264+c2*-1430)+s2*c*429))/128;
	case 73: return 315*((70+c2*(-95+c2*(-22+c2*143)))+s2*((95+c2*(132+c2*-2145))+s2*((-22+c2*2145)+s2*-143)))/256;
	case 74: return -3465*s*(c*(-15+c2*(24+c2*39))+s2*(c*(-24+c2*-130)+s2*c*39))/32;
	case 75: return -10395*((-14+c2*(3+c2*(30+c2*13)))+s2*((-3+c2*(-180+c2*-195))+s2*((30+c2*195)+s2*-13)))/64;
	case 76: return 135135*s*c*c2*c2;
	case 77: return 135135*c2*c2*c2;
	case 81: return 9*s*((1225+c2*(-3465+c2*(5005+c2*-5005)))+s2*((1155+c2*(-10010+c2*25025))+s2*((1001+c2*-15015)+s2*715)))/1024;
	case 82: return -315*(c*(-35+c2*(99+c2*(-143+c2*143)))+s2*(c*(-297+c2*(1430+c2*-3003))+s2*(c*(-715+c2*5005)+s2*c*-1001)))/1024;
	case 83: return 3465*s*((35+c2*(-51+c2*(-65+c2*273)))+s2*((17+c2*(130+c2*-1365))+s2*((-13+c2*819)+s2*-39)))/512;
	case 84: return 10395*(c*(35+c2*(-75+c2*(39+c2*65)))+s2*(c*(225+c2*(-390+c2*-1365))+s2*(c*(195+c2*2275)+s2*c*-455)))/512;
	case 85: return 135135*s*((7+c2*(9+c2*(-45+c2*-35)))+s2*((-3+c2*(90+c2*175))+s2*((-9+c2*-105)+s2*5)))/128;
	case 86: return -135135*(c*(-35+c2*(35+c2*(49+c2*15)))+s2*(c*(-105+c2*(-490+c2*-315))+s2*(c*(245+c2*525)+s2*c*-105)))/128;
	case 87: return 2027025*s*c2*c2*c2;
	case 88: return 2027025*c*c2*c2*c2;
	default: return NAN;
	}
}
// Legendre cos(phi)P'_n^m(sin(phi))
static fast_real cppnm(int_t n,int_t m,fast_real c,fast_real s){
	int_t degord=(m<1||m>n)?-1:n*10+m;
	fast_real c2=c*c,s2=s*s;
	switch(degord){
	case 21: return 3*(c2+s2*-1);
	case 22: return -6*s*c;
	case 31: return -3*s*((1+c2*-45)+s2*15)/8;
	case 32: return 15*(c+s2*c*-3);
	case 33: return -45*s*c2;
	case 41: return -5*(c2*(-1+c2*7)+s2*((1+c2*-42)+s2*7))/4;
	case 42: return 15*s*(c*(1+c2*7)+s2*c*-7);
	case 43: return 105*(c2*(1+c2)+s2*((-1+c2*-6)+s2))/2;
	case 44: return -420*s*c*c2;
	case 51: return -15*s*((2+c2*(-63+c2*525))+s2*((21+c2*-1050)+s2*105))/128;
	case 52: return -105*(c*(-2+c2*(3+c2*15))+s2*(c*(-9+c2*-150)+s2*c*75))/32;
	case 53: return 315*s*((-2+c2*(39+c2*75))+s2*((-13+c2*-150)+s2*15))/32;
	case 54: return 945*(c*(2+c2*(9+c2*5))+s2*(c*(-27+c2*-50)+s2*c*25))/16;
	case 55: return -4725*s*c2*c2;
	case 61: return 21*(c2*(5+c2*(-24+c2*99))+s2*((-5+c2*(144+c2*-1485))+s2*((-24+c2*1485)+s2*-99)))/128;
	case 62: return -105*s*(c*(-17+c2*(24+c2*297))+s2*(c*(-24+c2*-990)+s2*c*297))/64;
	case 63: return -945*(c2*(-3+c2*(8+c2*11))+s2*((3+c2*(-48+c2*-165))+s2*((8+c2*165)+s2*-11)))/32;
	case 64: return 945*s*(c*(5+c2*(104+c2*99))+s2*(c*(-104+c2*-330)+s2*c*99))/16;
	case 65: return 10395*(c2*(5+c2*(8+c2*3))+s2*((-5+c2*(-48+c2*-45))+s2*((8+c2*45)+s2*-3)))/16;
	case 66: return -62370*s*c*c2*c2;
	case 71: return -7*s*((25+c2*(-729+c2*(4125+c2*-21021)))+s2*((243+c2*(-8250+c2*105105))+s2*((825+c2*-63063)+s2*3003)))/1024;
	case 72: return 63*(c*(75+c2*(-171+c2*(55+c2*1001)))+s2*(c*(513+c2*(-550+c2*-21021))+s2*(c*(275+c2*35035)+s2*c*-7007)))/512;
	case 73: return -315*s*((45+c2*(-1053+c2*(3025+c2*7007)))+s2*((351+c2*(-6050+c2*-35035))+s2*((605+c2*21021)+s2*-1001)))/512;
	case 74: return -3465*(c*(-15+c2*(-9+c2*(125+c2*91)))+s2*(c*(27+c2*(-1250+c2*-1911))+s2*(c*(625+c2*3185)+s2*c*-637)))/128;
	case 75: return 10395*s*((-25+c2*(297+c2*(1075+c2*637)))+s2*((-99+c2*(-2150+c2*-3185))+s2*((215+c2*1911)+s2*-91)))/128;
	case 76: return 135135*(c*(5+c2*(27+c2*(25+c2*7)))+s2*(c*(-81+c2*(-250+c2*-147))+s2*(c*(125+c2*245)+s2*c*-49)))/64;
	case 77: return -945945*s*c2*c2*c2;
	case 81: return 9*(c2*(35+c2*(-154+c2*(429+c2*-1430)))+s2*((-35+c2*(924+c2*(-6435+c2*40040)))+s2*((-154+c2*(6435+c2*-100100))+s2*((-429+c2*40040)+s2*-1430))))/512;
	case 82: return 315*s*(c*(4+c2*(-11+c2*c2*143))+s2*(c*(11+c2*c2*-1001)+s2*(c*c2*1001+s2*c*-143)))/32;
	case 83: return 10395*(c2*(3+c2*(-10+c2*(13+c2*26)))+s2*((-3+c2*(60+c2*(-195+c2*-728)))+s2*((-10+c2*(195+c2*1820))+s2*((-13+c2*-728)+s2*26))))/256;
	case 84: return -10395*s*(c*(-5+c2*(-18+c2*(117+c2*130)))+s2*(c*(18+c2*(-390+c2*-910))+s2*(c*(117+c2*910)+s2*c*-130)))/32;
	case 85: return -135135*(c2*(-5+c2*(6+c2*(21+c2*10)))+s2*((5+c2*(-36+c2*(-315+c2*-280)))+s2*((6+c2*(315+c2*700))+s2*((-21+c2*-280)+s2*10))))/64;
	case 86: return 405405*s*(c*c2*(7+c2*(12+c2*5))+s2*(c*(-7+c2*(-40+c2*-35))+s2*(c*(12+c2*35)+s2*c*-5)))/4;
	case 87: return 2027025*(c2*(7+c2*(14+c2*(9+c2*2)))+s2*((-7+c2*(-84+c2*(-135+c2*-56)))+s2*((14+c2*(135+c2*140))+s2*((-9+c2*-56)+s2*2))))/32;
	case 88: return -16216200*s*c*c2*c2*c2;
	default: return NAN;
	}
}

class geopotential_reference{
public:
	struct component{
		int_t n,m;
		fast_real c,s;
	};
	std::vector<component> data;

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
			if(m==0){
				g+=rbarn*gc.c*fast_mpvec((n+1)*pn(n,s),0,-c*ppn(n,s));
			}
			else{
				fast_real cm=cos(m*lambda),sm=sin(m*lambda);
				fast_real tx=gc.c*cm+gc.s*sm,ty=-gc.c*sm+gc.s*cm;
				fast_real sp=spnm(n,m,c,s);
				fast_real cp=cppnm(n,m,c,s);
				g+=rbarn*fast_mpvec(-(n+1)*c*sp*tx,m*sp*ty,cp*tx);
			}
		}
		g/=r2;
		fast_mpvec x,y,z(-phi,lambda);
		x=r/r1;
		y=z*x;
		return g.x*x+g.y*y+g.z*z;
	}
};

int test_geopotential(){
	double start_time=CalcTime();
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
#ifndef USE_NEW_GEOPOTENTIAL
		geopotential *gp=geopotential::load(TEST_MEMPATH,REFERENCE_RADIUS_FACTOR);
#else
		geopotential *gp=geopotential::load(TEST_MEMPATH,REFERENCE_RADIUS_FACTOR,2);
#endif
		if(!gp)return 2;
		fast_real rn=TEST_RADIUS;
		for(int_t ir=0;ir<TEST_N;++ir){
			fast_real max_norm(0);
			fast_real max_diff(0);
			rn*=TEST_RADIUS_FACTOR;
			for(int_t itheta=0;itheta<=TEST_N;++itheta)for(int_t iphi=0;iphi<=2*TEST_N;++iphi){
				fast_mpvec r(itheta*pi/TEST_N,iphi*pi/TEST_N);
				r*=rn;
#ifndef USE_NEW_GEOPOTENTIAL
				fast_mpvec gtest=gp->sum(REFERENCE_RADIUS,r,2);
#else
				fast_mpvec gtest=gp->sum(REFERENCE_RADIUS,r);
#endif
				fast_mpvec gref=gpref.sum(REFERENCE_RADIUS*REFERENCE_RADIUS_FACTOR,r);
				max_norm=std::max(max_norm,gtest.norm());
				max_diff=std::max(max_diff,(gtest-gref).norm());
			}
			max_rerr=std::max(max_rerr,max_diff/max_norm);
		}
		geopotential::unload(gp);

		if(!(max_rerr<TEST_EPSILON)){
			fprintf(stderr,
				"\nMax Rel. Error %.16le Too Large at n=%lld, m=%lld",
				max_rerr,n,m);
			return 1;
		}
	}

	printf("Passed(%f), Done in %fs",max_rerr/TEST_EPSILON,CalcTime()-start_time);
    return 0;
}
