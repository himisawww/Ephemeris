#include"ring.h"
#include"htl/vector.h"
#include<string>
#include"utils/memio.h"
#include"utils/logger.h"

#define MAX_LINESIZE 1024

struct disk{
    fast_real Gs;
    fast_real R;
    fast_real H;
};

void adddisk(htl::vector<disk> &disks,fast_real Gs,fast_real R,fast_real H){
    for(auto &d:disks){
        if(d.R==R&&d.H==H){
            d.Gs+=Gs;
            return;
        }
    }
    disks.push_back({Gs,R,H});
}

const ring *ring::load(const char *file,fast_real ref_GM,fast_real ref_R,fast_real direction_mass_factor){
    std::string linebuf;
    MFILE *fin=mopen(file);
    if(!fin){
        LogError("%s : Error opening ring model\n",file);
        return nullptr;
    }

    using Constants::pi;

    bool failed=false;
    char rname[MAX_LINESIZE];
    htl::vector<disk> disks;
    while((linebuf=readline(fin)).size()){
        const char *chbuf=linebuf.c_str();
        if(linebuf.size()>=MAX_LINESIZE){
            LogError(
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
            LogError(
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
    ret->GL_R2=0;

    fast_real ref_R2=ref_R*ref_R;
    fast_real ref_GMR2=ref_GM*ref_R2;

    for(int_t i=0;i<N;++i){
        const auto &d=disks[i];
        fast_real R2=d.R*d.R;
        fast_real H2=d.H*d.H;
        fast_real GM=pi*R2*d.Gs;
        ret->GM_ratio+=GM;
        ret->A+=GM*(6*R2+H2)/(12*ref_GMR2);
        ret->J2+=GM*(3*R2-H2)/(12*ref_GMR2);
        ret->GL_R2+=4*GM*sqrt(ref_GM*d.R)/(5*ref_GMR2);

        ret->c_table[i].GsR2=d.Gs/ref_GM*ref_R2;
        ret->c_table[i].R_R=d.R/ref_R;
        ret->c_table[i].H_R=d.H/ref_R;
    }

    ret->GM_ratio=1/(1+ref_GM/ret->GM_ratio);
    //retrograde ring, is this possible?
    if(direction_mass_factor<0)ret->GL_R2*=-1;
    return ret;
}

void ring::unload(const ring *rp){
    free((void*)rp);
}

int_t ring::size() const{
    return sizeof(ring)+3*sizeof(fast_real)*N;
}

const ring *ring::copy(const ring *rp,fast_real multiplier){
    if(!rp||!multiplier)return nullptr;
    size_t rpsize=rp->size();
    ring *ret=(ring *)malloc(rpsize);
    memcpy(ret,rp,rpsize);
    if(multiplier!=1){
        fast_real host_ratio=1-ret->GM_ratio;
        fast_real new_ratio=ret->GM_ratio*multiplier;
        ret->GM_ratio=new_ratio/(host_ratio+new_ratio);
        ret->A*=multiplier;
        ret->J2*=multiplier;
        ret->GL_R2*=multiplier;
        for(int_t i=0;i<ret->N;++i)
            ret->c_table[i].GsR2*=multiplier;
    }
    return ret;
}

bool ring::is_same(const ring *rp,const ring *rq){
    if(rp==rq)return true;
    if(rp&&rq){
        size_t rsize=rp->size();
        if(rq->size()!=rsize)return false;
        return memcmp(rp,rq,rsize)==0;
    }
    return false;
}

#define CONST_TABLE static const
#include"ring.impl"
