#include"ephemeris.h"
#include"utils.h"
#include<map>
#include<thread>

//used to freeze rotation while keep angular momentum
//also used to reduce the mass of a system to make it a test mass
#define VERY_BIG_VALUE  (1.e16)
//minimum steps per rotation period for non-parent objects
#define ROTATION_STEP_LIMIT (10)

const fast_real                 //combine threshold (only for Solar System)
    GM_max_child    =   1E13,   //( Ganymede < * < Mercury )
    GM_max_parent   =   2E17,   //( Jupiter < * < Sun )
//  GM_max_tiny     =   5E11;   //( Haumea < * < Pluto, Eris )
    GM_max_tiny     =  14E11;   //( Pluto, Eris < * < Triton )

#define IS_MINOR(GM)        ((GM)<GM_max_tiny)

void mass_copy(mass &dst,const mass &src){
    dst.r=src.r;
    dst.v=src.v;
    dst.GL=src.GL;
    dst.s=src.s;
    dst.w=src.w;
}

void thread_work(msystem *ms,fast_real dt,int_t n_combine){
    ms->integrate(dt,n_combine);
}

void msystem::combined_integrate(fast_real dt,int_t n_combine,int_t n_step,int USE_GPU){
    int_t bn=blist.size();
    std::map<int_t,int_t> clist;
    std::map<int_t,std::vector<int_t>> cvecs;
    for(int_t id=0;id<bn;++id){
        const auto &b=blist[id];
        if(b.hid>=0)continue;

        const auto &t=blist[b.tid];
        if(t.pid<0)continue;
        const auto &tp=blist[t.pid];
        //combine target = tp.mid;

        //combine threshold (only for Solar System)
        //GM < 1E13 ( Ganymede < * < Mercury )
        //target GM < 2E17 ( Jupiter < * < Sun )
        if((fast_real)t.GM_sys>GM_max_child
            ||(fast_real)tp.GM>GM_max_parent)continue;

        clist.insert({b.mid,tp.mid});
    }

    do{
        bool changed=false;
        for(auto &c:clist){
            auto fs=clist.find(c.second);
            if(fs==clist.end())continue;

            c.second=fs->second;
            changed=true;
        }
        if(!changed)break;
    } while(1);

    for(const auto &c:clist){
        auto is=cvecs.insert({c.second,std::vector<int_t>()});
        is.first->second.push_back(c.first);
    }

    msystem Sc,Sx;
    std::map<int_t,msystem> Sn;
    //major subsystems
    std::vector<int_t> Smajsub;
    int_t mn=mlist.size();

    for(const auto &p:cvecs)
        Sn.insert({p.first,msystem()});
    
    fast_real dt_long=dt*n_combine;
    for(int_t i_step=0;i_step<n_step;++i_step){
        //initialize
        double s;
        double sp=CalcTime();

        Sc.mlist.clear();

        Sx.mlist.clear();

        for(auto &sns:Sn){
            sns.second.tidal_childlist.clear();
            sns.second.mlist.clear();
        }

        Smajsub.clear();

        for(int_t i=0;i<mn;++i){
            const mass &mi=mlist[i];

            // if is a child
            auto clp=clist.find(i);
            if(clp!=clist.end()){
                // it will only appear in its parents' subsystem
                msystem &sn=Sn[clp->second];
                sn.tidal_childlist.push_back(sn.mlist.size());
                sn.mlist.push_back(mi);
                continue;
            }

            // if is a parent
            auto cvp=cvecs.find(i);
            if(cvp!=cvecs.end()){
                auto &chs=cvp->second;//children
                //combined parent
                mass mci=mi;
                mci.gpmodel=nullptr;
                
                if(mci.ringmodel){
                    ring &mr=*mci.ringmodel;
                    mci.ringmodel=nullptr;
                    fast_real rGM=mr.GM_ratio;
                    rGM/=1-rGM;
                    mci.GM+=mci.GM*rGM;
                    mci.GM0+=mci.GM*rGM;

                    fast_mpvec mgl=mci.GL;
                    fast_mpmat fgls(mgl.perpunit(),0,mgl/mgl.norm());
                    fast_mpmat rh(0);
                    rh.x.x=mr.J2/2;
                    rh.y.y=mr.J2/2;
                    rh.z.z=-mr.J2;
                    mci.C_potential+=fgls.toworld(rh);
                    mci.C_potential-=mci.C_potential*mr.GM_ratio;
                }

                // clean parent rotational informations
                mci.k2=0;
                mci.k2r=0;
                mci.dJ2=0;
                mci.C_static=mci.C_potential;
                mci.s=1;
                mci.A*=VERY_BIG_VALUE;
                mci.inertia*=VERY_BIG_VALUE;
                mci.w/=VERY_BIG_VALUE;
                
                //newtonian combination
                real GM=mci.GM;
                real GM0=mci.GM0;
                mpvec r=mi.r*GM,v=mi.v*GM;
                for(auto &j:chs){
                    const mass &mj=mlist[j];
                    real GMj=mj.GM;
                    r+=mj.r*GMj;
                    v+=mj.v*GMj;
                    GM+=GMj;
                    GM0+=mj.GM0;

                    mci.dGM+=mj.dGM;
                    mci.lum+=mj.lum;
                    
                    //CombineCorrector
                    //assume dr.dv==0
                    fast_mpvec dr=mj.r-mci.r;
                    fast_mpvec dv=mj.v-mci.v;
                    fast_mpvec dj=dr*dv;
                    mci.GL+=(fast_mpvec(mj.GL)+dj)*(mj.GM/mi.GM);

                    fast_real r2=dr%dr,rr=1/sqrt(r2);
                    fast_real j2=dj%dj;
                    fast_mpmat dh;
                    if(j2==0){
                        dh=fast_mpmat(1)-fast_mpmat((3/r2)*dr,dr);
                    }
                    else{
                        fast_real rj=1/sqrt(j2);
                        fast_mpmat mrot(dr*rr,0,dj*rj);
                        fast_real dw=1/(r2*rj);
                        fast_real t=dt_long*dw;
                        fast_real st=sin(t),ct=cos(t),sit=3*st/t;
                        dh=fast_mpmat(
                            fast_mpvec((-1-sit*ct)/2,-sit*st/2,0),
                            fast_mpvec(-sit*st/2,(-1+sit*ct)/2,0),
                            fast_mpvec(0,0,1)
                        );
                        dh=mrot.toworld(dh);
                    }
                    fast_real dhnorm=mj.GM+mci.GM;
                    dhnorm=-r2*mj.GM*mci.GM/(dhnorm*dhnorm*2*mci.R2);
                    dh*=dhnorm;
                    mci.C_static+=dh;
                }
                r/=GM;
                v/=GM;
                mci.r=r;
                mci.v=v;
                mci.GM0=(fast_real)GM0;
                mci.GM=(fast_real)GM;

                //parent will present in Sc
                Sc.mlist.push_back(mci);

                //parent will always in its own subsystem
                msystem &sn=Sn[i];
                sn.tidal_parent=sn.mlist.size();
                sn.mlist.push_back(mi);

                if(!IS_MINOR(mci.GM0)){
                    Sx.mlist.push_back(mci);
                    Smajsub.push_back(i);

                    //parent will present in all Sn as combined form
                    //if not in its own subsystem
                    for(auto &sns:Sn)if(sns.first!=i)
                        sns.second.mlist.push_back(mci);
                }
                else{
                    //tiny will not present in Sn, and only as a test mass in Sx
                    mci.scale(1/VERY_BIG_VALUE);
                    Sx.mlist.push_back(mci);
                }
                continue;
            }

            // not in parent-child relationship
            // will present in Sc
            Sc.mlist.push_back(mi);

            // tiny will not present in Sx and Sn
            if(!IS_MINOR(mi.GM0)){
                Sx.mlist.push_back(mi);

                for(auto &sns:Sn)
                    sns.second.mlist.push_back(mi);
            }
        }

        Sc.t_eph=t_eph;
        Sc.accel();

        Sx.t_eph=t_eph;
        Sx.accel();

        for(auto &sns:Sn){
            msystem &sn=sns.second;
            sn.tidal_matrix=0;
            const mass &mp=Sc.mlist[Sc.get_mid(mlist[sns.first].sid)];
            for(const auto &mj:Sc.mlist){
                if(IS_MINOR(mj.GM0)&&mj.sid!=mp.sid){
                    fast_mpvec rm=mj.r-mp.r;
                    fast_mpvec vm=mj.v-mp.v;
                    rm+=vm*(dt_long/2);
                    fast_real rr2=1/(rm%rm),rr=sqrt(rr2);
                    fast_real tmatn=-mj.GM0*rr2*rr;
                    sn.tidal_matrix+=fast_mpmat(tmatn)-fast_mpmat(
                        (3*tmatn*rr2)*rm,rm
                    );
                }
            }
            sn.t_eph=t_eph;
            sn.accel();
        }

        s=CalcTime();
        sp=s-sp;

        // integrate
#ifndef NDEBUG
        double sc=CalcTime();
        Sc.integrate(dt_long,1,USE_GPU);
        s=CalcTime();
        sc=s-sc;

        double sx=s;
        Sx.integrate(dt_long,1,0);
        s=CalcTime();
        sx=s-sx;

        double sn=s;
        for(auto &sns:Sn){
            sns.second.integrate(dt,n_combine,0);
        }
        s=CalcTime();
        sn=s-sn;
#else
        std::vector<std::thread> threads;
        threads.reserve(Sn.size());
        for(auto &sns:Sn){
            threads.push_back(std::thread(thread_work,&sns.second,dt,n_combine));
        }

        double sx=s;
        Sx.integrate(dt_long,1,0);
        s=CalcTime();
        sx=s-sx;

        double sc=CalcTime();
        Sc.integrate(dt_long,1,USE_GPU);
        s=CalcTime();
        sc=s-sc;

        double sn=s;
        for(auto &th:threads){
            th.join();
        }
        s=CalcTime();
        sn=s-sn;
#endif
        //finalize
        double sf=CalcTime();

        for(int_t i=0;i<mn;++i){
            mass &mi=mlist[i];

            // if is a child
            auto clp=clist.find(i);
            if(clp!=clist.end()){
                continue;
            }

            // if is a parent
            auto cvp=cvecs.find(i);
            if(cvp!=cvecs.end()){
                auto &chs=cvp->second;
                mass &Sci=Sc.mlist[Sc.get_mid(mi.sid)];
                mass &Sxi=Sx.mlist[Sx.get_mid(mi.sid)];

                //sum for perturbation of children of other subsystems
                mpvec dr,dv;
                dr=Sci.r-Sxi.r;
                dv=Sci.v-Sxi.v;

                if(!IS_MINOR(Sci.GM0)){
                    for(auto &si:Smajsub)if(si!=i){
                        msystem &Sj=Sn[si];
                        mass &Sji=Sj.mlist[Sj.get_mid(mi.sid)];
                        dr+=Sji.r-Sxi.r;
                        dv+=Sji.v-Sxi.v;
                    }
                }
                msystem &Si=Sn[i];
                mass &Sni=Si.mlist[Si.get_mid(mi.sid)];

                mass_copy(mi,Sni);
                mi.r+=dr;
                mi.v+=dv;

                for(auto &j:chs){
                    mass &mj=mlist[j];
                    mass &Snj=Si.mlist[Si.get_mid(mj.sid)];
                    mass_copy(mj,Snj);
                    mj.r+=dr;
                    mj.v+=dv;
                }
                continue;
            }

            // not in parent-child relationship
            // will present in Sc
            mass &Sci=Sc.mlist[Sc.get_mid(mi.sid)];
            mass_copy(mi,Sci);

            // tiny will not present in Sx and Sn
            if(!IS_MINOR(mi.GM0)){
                mpvec dr(0),dv(0);
                mass &Sxi=Sx.mlist[Sx.get_mid(mi.sid)];

                for(auto &si:Smajsub){
                    msystem &Sj=Sn[si];
                    mass &Sji=Sj.mlist[Sj.get_mid(mi.sid)];
                    dr+=Sji.r-Sxi.r;
                    dv+=Sji.v-Sxi.v;
                }
                mi.r+=dr;
                mi.v+=dv;
            }

            //fix rotations
            if(mi.w.norm()*std::abs(dt_long)>2*pi/ROTATION_STEP_LIMIT){
                fast_mpvec zl=mi.s.z*mi.GL;
                fast_real l=fast_mpvec(mi.GL).norm();
                if(l>0){
                    zl/=l;
                    mi.s+=rotation_matrix(zl,fast_real(1))%fast_mpmat(mi.s);
                }
            }
        }

        t_eph=Sc.t_eph;
        update((fast_real)t_eph);
        accel();

        s=CalcTime();
        sf=s-sf;
        if(0)printf("\n"
            "prepare: %lfs\n"
            "int(Sc): %lfs\n"
            "int(Sx): %lfs\n"
            "int(Sn): %lfs\n"
            "  final: %lfs\n",
            sp,sc,sx,sn,sf);
    }

    return;
}
