#include"geopotential.h"
#include<string>
#include<cstdio>
#include<map>
#include"utils/memio.h"
#include"utils/logger.h"
#include"geopotential.impl"

#define Precompute_Table_size(nz,nt) ((nt)*((nt)+1)+(nz)-3)

const geopotential *geopotential::load(const char *file,fast_real ref_radius_factor,int_t N_start){
    //we support at most order 8 harmonics...
    static constexpr int_t Max_N=8;
    std::string linebuf;
    MFILE *fin=mopen(file);
    if(!fin){
        LogError("%s : Error opening harmonics model\n",file);
        return nullptr;
    }
    int_t Nz=-1,Nt;
    int_t lskip=0,hskip=0,dwarning=0;
    int_t Nzmax=0,Ntmax=0;
    std::map<std::pair<int_t,int_t>,std::pair<fast_real,fast_real>> coefs;
    while((linebuf=readline(fin)).size()){
        const char *chbuf=linebuf.c_str();
        if(Nz<0){
            if(2!=sscanf(chbuf,"%lld%lld",&Nz,&Nt)){
                LogError("%s : Error reading max numbers of zonal/tesseral degrees\n",file);
                break;
            }
            int_t N=std::max(Nz,Nt);
            if(N>Max_N){
                LogWarning("%s : Geopotential degree[%lld] is too high, truncated to supported maximum[%lld].\n",
                    file,N,Max_N);
                Nz=std::min(Nz,Max_N);
                Nt=std::min(Nt,Max_N);
                N=Max_N;
            }
            if(N<N_start){
                LogWarning("%s : Geopotential degree[%lld] is too low. Ignored.\n",
                    file,N);
                break;
            }
            continue;
        }
        int_t n=0,m=0;
        fast_real Cnm=0,Snm=0;
        int_t nr=sscanf(chbuf,"%lld%lld%lf%lf",&n,&m,&Cnm,&Snm);
        if(nr<3||n<2||m<0||m>n||nr==3&&m!=0||m==0&&Snm!=0){
            LogError("%s : %s :\n"
                "\tshould in format of \n"
                "\tdegree(n>=2)  order(m==0)     Jn   [0]\n"
                "\t\tor\n"
                "\tdegree(n>=2)  order(1<=m<=n)  Cnm  Snm\n",file,chbuf);
            break;
        }
        if(n<N_start){
            ++lskip;
            continue;
        }
        if(n>(m>0?Nt:Nz)){
            ++hskip;
            continue;
        }
        double rrffactor=std::pow(ref_radius_factor,n);
        Cnm*=rrffactor;
        Snm*=rrffactor;
        auto ir=coefs.insert({{n,m},{Cnm,Snm}});
        if(ir.second){
            if(m==0)
                Nzmax=std::max(Nzmax,n);
            else
                Ntmax=std::max(Ntmax,n);
        }
        else{
            LogWarning("%s : %s :\n"
                "Duplicated terms will be superposed. Is this intentional?\n",
                file,chbuf);
            auto &CSnm=ir.first->second;
            CSnm.first+=Cnm;
            CSnm.second+=Snm;
        }
    }
    fclose(fin);

    if(lskip)
        LogWarning(
            "%s : Ignored [%lld] lines with degree < [%lld].\n"
            "\tNote: 2nd degree harmonics are handled separately in this program.\n",
            file,lskip,N_start);
    if(hskip)
        LogWarning("%s : Ignored [%lld] lines with degree higher than maximum\n",
            file,hskip);

    geopotential *ret=nullptr;
    do{
        if(!std::max(Nzmax,Ntmax)){
            LogWarning("%s : No valid geopotential terms. Ignored.\n",file);
            break;
        }
        if(Nzmax<Nz)
            LogWarning(
                "%s : Maximum degree of zonal terms less than specified [%lld < %lld].\n",
                file,Nzmax,Nz);
        if(Ntmax<Nt)
            LogWarning(
                "%s : Maximum degree of tesseral terms less than specified [%lld < %lld].\n",
                file,Ntmax,Nt);
        //avoid negative size
        Nz=std::max(Nzmax,int_t(1));
        Nt=std::max(Ntmax,int_t(1));

        int_t csize=sizeof(fast_real)*Precompute_Table_size(Nz,Nt);
        ret=(geopotential *)malloc(sizeof(geopotential)+csize);
        ret->Nz=Nz;
        ret->Nt=Nt;
        fast_real *c_table=(fast_real *)ret->c_table;
        memset(c_table,0,csize);
        for(auto &c:coefs){
            int_t n=c.first.first;
            int_t m=c.first.second;
            auto &cs=c.second;
            if(m==0){
                fast_real *JCSn=c_table+Nt*(Nt+1)+n-4;
                JCSn[0]=cs.first;
            }
            else{
                fast_real *JCSn=c_table+(Nt-m+1)*(Nt-m)+2*(n-(m+(m==1)));
                JCSn[0]=cs.first;
                JCSn[1]=cs.second;
            }
        }
    } while(0);

    return ret;
}

void geopotential::unload(const geopotential *gp){
    free((void*)gp);
}

int_t geopotential::size() const{
    return sizeof(geopotential)+sizeof(fast_real)*Precompute_Table_size(Nz,Nt);
}

const geopotential *geopotential::copy(const geopotential *gp,fast_real multiplier){
    if(!gp)return nullptr;
    size_t gpsize=gp->size();
    geopotential *ret=(geopotential *)malloc(gpsize);
    memcpy(ret,gp,gpsize);
    if(multiplier!=1){
        size_t n=Precompute_Table_size(ret->Nz,ret->Nt);
        for(int_t i=0;i<n;++i)
            ret->c_table[i]*=multiplier;
    }
    return ret;
}
