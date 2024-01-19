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
