#include"geopotential.h"
#include<string>
#include<cstdio>
#include"utils/memio.h"

#include"geopotential.impl"

geopotential *geopotential::load(const char *file,fast_real ref_radius_factor,int_t N_start){
    std::string linebuf;
    MFILE *fin=mopen(file);
    if(!fin){
        fprintf(stderr,"%s : Error opening harmonics model\n",file);
        return nullptr;
    }
    int_t Nz=-1,Nt;
    fast_real *c_table;
    geopotential *ret=nullptr;
    int_t lskip=0,hskip=0,dwarning=0;
    int_t Nzmax=0,Ntmax=0;
    while((linebuf=readline(fin)).size()){
        const char *chbuf=linebuf.c_str();
        if(Nz<0){
            if(2!=sscanf(chbuf,"%lld%lld",&Nz,&Nt)){
                fprintf(stderr,"%s : Error reading max numbers of zonal/tesseral degrees\n",file);
                break;
            }
            if(Nt<0)Nt=0;
            if(Nz<Nt)Nz=Nt;
            if(Nz>Max_N){
                fprintf(stderr,"%s : Geopotential degree[%lld] is too high, truncated to supported maximum[%lld].\n",
                    file,Nz,Max_N);
                Nz=Max_N;
                if(Nt>Max_N)Nt=Max_N;
            }
            if(Nz<N_start){
                fprintf(stderr,"%s : Geopotential degree[%lld] is too low. Ignored.\n",
                    file,Nz);
                break;
            }
            int_t csize=sizeof(fast_real)*Precompute_Table_size(Nz);
            ret=(geopotential *)malloc(sizeof(geopotential)+csize);
            ret->Nz=Nz;
            ret->Nt=Nt;
            c_table=(fast_real *)ret->c_table;
            memset(c_table,0,csize);
            continue;
        }
        int_t n=0,m=0;
        fast_real Cnm=0,Snm=0;
        int_t nr=sscanf(chbuf,"%lld%lld%lf%lf",&n,&m,&Cnm,&Snm);
        if(nr<3||n<2||m<0||m>n||nr==3&&m!=0||m==0&&Snm!=0){
            fprintf(stderr,"%s : %s :\n"
                "\tshould in format of \n"
                "\tdegree(n>=2)  order(m==0)     Jn   [0]\n"
                "\t\tor\n"
                "\tdegree(n>=2)  order(1<=m<=n)  Cnm  Snm\n",file,chbuf);
            free(ret);
            ret=nullptr;
            break;
        }
        if(n<N_start){
            ++lskip;
            continue;
        }
        if(n>Nz||n>Nt&&m>0){
            ++hskip;
            continue;
        }
        double rrffactor=std::pow(ref_radius_factor,n);
        Cnm*=rrffactor;
        Snm*=rrffactor;
        bool isdup=false;
        if(m==0){
            fast_real &jn=c_table[J(n)];
            if(jn)isdup=true;
            jn+=Cnm;
        }
        else{
            fast_real &cnm=c_table[C(n,m)];
            fast_real &snm=c_table[S(n,m)];
            if(cnm||snm)isdup=true;
            cnm+=Cnm;
            snm+=Snm;
            Ntmax=std::max(Ntmax,n);
        }
        Nzmax=std::max(Nzmax,n);
        if(isdup){
            fprintf(stderr,"%s : %s :\n"
                "Duplicated terms will be superposed. Is this intentional?\n",
                file,chbuf);
        }
    }
    fclose(fin);
    if(!ret)return nullptr;

    if(lskip)
        fprintf(stderr,
            "%s : Ignored [%lld] lines with degree < [%lld].\n"
            "\tNote: 2nd degree harmonics are handled separately in this program.\n",
            file,lskip,N_start);
    if(hskip)
        fprintf(stderr,"%s : Ignored [%lld] lines with degree higher than maximum\n",
            file,hskip);

    if(!Nzmax){
        fprintf(stderr,"%s : No valid geopotential terms. Ignored.\n",file);
        free(ret);
        ret=nullptr;
    }
    else{
        if(Nzmax<Nz){
            fprintf(stderr,
                "%s : Maximum degree of zonal terms less than specified [%lld < %lld].\n",
                file,Nzmax,Nz);
            ret->Nz=Nzmax;
        }
        if(Ntmax<Nt){
            fprintf(stderr,
                "%s : Maximum degree of tesseral terms less than specified [%lld < %lld].\n",
                file,Ntmax,Nt);
            ret->Nt=Ntmax;
        }
    }

    return ret;
}

void geopotential::unload(geopotential *gp){
    free(gp);
}

int_t geopotential::size() const{
    return sizeof(geopotential)+sizeof(fast_real)*Precompute_Table_size(Nz);
}

