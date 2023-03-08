#include"ephemeris.h"
#include"utils.h"

#define MAX_LINESIZE 1024

struct disk{
    fast_real Gs;
    fast_real R;
    fast_real H;
};

void adddisk(std::vector<disk> &disks,fast_real Gs,fast_real R,fast_real H){
    for(auto &d:disks){
        if(d.R==R&&d.H==H){
            d.Gs+=Gs;
            return;
        }
    }
    disks.push_back({Gs,R,H});
}

ring *ring::load(const char *file,fast_real ref_GM,fast_real ref_R2,fast_real direction_mass_factor){
    std::string linebuf;
    FILE *fin=fopen(file,"r");
    if(!fin){
        fprintf(stderr,"%s : Error opening ring model\n",file);
        return nullptr;
    }

    bool failed=false;
    char rname[MAX_LINESIZE];
    std::vector<disk> disks;
    while((linebuf=readline(fin)).size()){
        const char *chbuf=linebuf.c_str();
        if(linebuf.size()>=MAX_LINESIZE){
            fprintf(stderr,
                "Loading %s\n Line Too Long: %s\n",
                file,chbuf);
            failed=true;
            break;
        }

        double rin,rout,thick,GM;
        int n=sscanf(chbuf,
            "%[^\t]%lf%lf%lf%lf",
            rname,&rin,&rout,&thick,&GM
            );

        if(n!=5){
            fprintf(stderr,
                "Loading %s\n Line Length Error: %s\n",
                file,chbuf);
            failed=true;
            break;
        }
        GM*=std::abs(direction_mass_factor);
        //convert km to m
        rin*=1000;
        rout*=1000;
        thick*=1000;

        //surface grav density
        double Gs=GM/(pi*(rout-rin)*(rout+rin));

        adddisk(disks,-Gs, rin,thick);
        adddisk(disks, Gs,rout,thick);
    }
    fclose(fin);
    if(failed)return nullptr;

    int_t N=disks.size();
    int_t csize=3*sizeof(fast_real)*(N);//radius, thickness, Gsigma
    ring *ret=(ring*)malloc(sizeof(ring)+csize);

    ret->N=N;
    ret->GM_ratio=0;
    ret->A=0;
    ret->J2=0;
    ret->GL=0;

    fast_real ref_GMR2=ref_GM*ref_R2;

    for(int_t i=0;i<N;++i){
        const auto &d=disks[i];
        fast_real R2=d.R*d.R;
        fast_real H2=d.H*d.H;
        fast_real GM=pi*R2*d.Gs;
        ret->GM_ratio+=GM;
        ret->A+=GM*(6*R2+H2)/(12*ref_GMR2);
        ret->J2+=GM*(3*R2-H2)/(12*ref_GMR2);
        ret->GL+=4*GM*sqrt(ref_GM*d.R)/(5*ref_GM);
        //since host GM&R not known for now
        //GL is now G * angular momentum / sqrt(host GM)
        //and J2, A are not normalized by (host GM*R^2)

        ret->c_table[i].Gs=d.Gs/ref_GM;
        ret->c_table[i].R=d.R;
        ret->c_table[i].H=d.H;
    }

    ret->GM_ratio=1/(1+ref_GM/ret->GM_ratio);
    //retrograde ring, is this possible?
    if(direction_mass_factor<0)ret->GL*=-1;
    return ret;
}

void ring::unload(ring *rp){
    free(rp);
}

int_t ring::size() const{
    return sizeof(ring)+3*sizeof(fast_real)*N;
}

#define CONST_TABLE static const
#include"disk_approx.impl"
float fzcorr(float x,float y,int kh){
    const size_t N=32;

    x*=N;
    y*=N;
    int i=floor(x),j=floor(y);
    j-=kh*N;
    if(i<0)i=0;
    if(j<0)j=0;
    if(i>=N)i=N-1;
    if(j>=N)j=N-1;
    j+=kh*N;
    x-=i;
    y-=j;
    const float *p=(const float *)disk_approx_table+16*(2*N*i+j);
    return  p[0]+x*(p[4]+x*(p[ 8]+x*p[12]))
        +y*(p[1]+x*(p[5]+x*(p[ 9]+x*p[13]))
        +y*(p[2]+x*(p[6]+x*(p[10]+x*p[14]))
        +y*(p[3]+x*(p[7]+x*(p[11]+x*p[15])))));
}
static double padesum(double x,const double coef[][2],int n){
    double resn=0,resd=0;
    while(n>0){
        --n;
        resn=x*resn+coef[n][0];
        resd=x*resd+coef[n][1];
    }
    return resn/resd;
}
#define PadeSum(x,coef) padesum(x,coef,sizeof(coef)/sizeof(coef[0]))

fast_mpvec ring::sum(fast_mpvec r) const{
    const double xserkk[][2]={{-4228.1846570321674587,6387.9414439779279738},{-13276.372802884782407,16984.953742612298093},{-11316.549465243054011,8997.7319379726254196},{-3242.2867356286923985,403.15182320396046619},{-704.60633921130372399,-5.7789477668119521439}};
    const double xserkm[][2]={{-2530.8255093002657219,7950.8733392382639504},{-5471.7803365810451727,17189.101860053144376},{-1624.0389725420726306,7100.5562444641275465},{1164.6058458813108946,532.54915854738706820},{270.03897254207263065,-5.0806023029229406554}};
    const double xser0[][2]={{3216.9908772759482762,32768.000000000000000},{-6727.1375380948220662,-93098.060290997395032},{4534.8257724052420971,96814.901233532309690},{-1004.5251977686921953,-43973.540867249017404},{17.915064837735955511,7754.0323173314244943},{-0.039157108395559676838,-224.00490789882410101}};

    const double kserkk[][2]={{-2373.3522804139070577,1124.9866194216978637},{-8300.7524798211745830,17215.217000406477760},{15279.068333441193395,26590.060415743810476},{-4354.4701057632348105,-9776.5071563707548874},{-250.49346744287694398,-2385.7568792012312128}};
    const double kserkm[][2]={{-523.02599873176918199,290.10526315789473684},{-5485.2200380469245403,5336.9790741915028535},{-8338.1052631578947368,15206.361445783132530},{-2017.9378566899175650,10515.743817374762207},{-19.710843373493975904,1418.8103994927076728}};
    const double kser0[][2]={{51471.854036415172419,32768.000000000000000},{-119621.18573838632696,-84345.212035111671234},{96754.491735189724987,78074.125512775430927},{-31810.735585315211544,-31108.829182924968396},{3678.1383813624327949,4926.4461273804873303},{-74.040670509817904983,-206.66608791871105442}};

    bool psign;
    if(psign=r.z<0)r.z=-r.z;

    const double re2(r.x*r.x+r.y*r.y);
    const double re=sqrt(re2);
    const double rz2(r.z*r.z);
    const double r2=re2+rz2;
    const double rz(r.z);

    fast_mpvec ret(0);

    for(int_t i=0;i<N;++i){
        double R=c_table[i].R;
        double H=c_table[i].H;

        double reR=re*R;
        double R2=R*R,H2,fzh=1;
        double k2dre=4/(r2+R2+2*reR);
        if(H>0){
            H2=H*H;
            k2dre-=H2*(k2dre*k2dre)*(1-k2dre*rz2)/48;
            fzh=(re>R?re-R:0);
            fzh=4*(rz2+fzh*fzh);
            fzh=sqrt(fzh/(fzh+H2));
        }
        double k2=reR*k2dre;
        k2dre*=R2*sqrt(k2dre);
        double phi=atan2(rz,re-R);
        double k,lk,fx;
        bool lkh=k2>=0.45;
        if(lkh){
            k=sqrt(k2);lk=log((1-k2)/(8*(1+k)));
            fx=(PadeSum(k,xserkk)+lk*PadeSum(k,xserkm))/(k2*k2);
        }
        else{
            fx=PadeSum(k2,xser0);
        }
        fx*=k2dre;

        //approx for fz=diskfz_approx(phi,k2)
        double fz;
        double s=sin(phi),c=cos(phi);
        bool kh=k2>=.5;
        bool exphi;

        if(kh){
            const double phimax=(pi/2);
            exphi=c>0;
            double ck=1-k2,lk=log(ck/16),sk=sqrt(ck);
            double cc=exphi?-c:c;
            double cphi=exphi?pi-phi:phi;
            fz=cphi-(s/8)*((lk*(k2-5)-2*ck)*sk+2*cc*(1+lk)*ck);
            fz*=fzcorr((exphi?phi:pi-phi)*(1/phimax),k2*2,kh);
        }
        else{
            const double phimax=(809./512);
            exphi=phi>phimax;
            double cr2=exphi?R2+rz2:r2;
            double cRdr2=(exphi?re2:R2)/cr2;
            fz=(pi/8)*(4-3*cRdr2)*(rz/sqrt(cr2)*cRdr2);
            fz*=fzcorr((exphi?pi-phi:phi)*(1/phimax),k2*2,kh);
        }

        if(exphi){
            double ek;
            if(lkh){
                ek=PadeSum(k,kserkk)+lk*PadeSum(k,kserkm);
            }
            else{
                ek=PadeSum(k2,kser0);
            }
            fz=pi-2*sqrt(1-k2)*ek*s-fz;
        }

        double Gs=c_table[i].Gs;
        fx*=Gs;
        fz*=Gs;

        ret+=fast_mpvec(-4*fx*r.x,-4*fx*r.y,-2*fz*fzh);
    }

    if(psign)ret.z=-ret.z;

    return ret;
}
