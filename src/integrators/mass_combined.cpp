#include"mass_combined.h"
#include<map>
#include"math/distribute.h"
#include"physics/geopotential.h"
#include"physics/ring.h"
#include"utils/threadpool.h"

//used to freeze rotation while keep angular momentum
//also used to reduce the mass of a system to make it a test mass
#define VERY_BIG_VALUE  (1.e16)
//minimum steps per rotation period for non-parent objects
#define ROTATION_STEP_LIMIT (10)

#define IS_MINOR(GM)        ((GM)<GM_max_tiny)

void mass_copy(mass &dst,const mass &src){
    dst.r=src.r;
    dst.v=src.v;
    dst.GL=src.GL;
    dst.s=src.s;
    dst.w=src.w;
}

struct thread_works{
    const std::vector<msystem*> *pm;
    const std::vector<std::vector<size_t>> *widxs;
    fast_real dt;
    int_t n_combine;
};

void do_thread_works(void *pworks,size_t thread_id){
    const thread_works &w=*(const thread_works*)pworks;
    const auto &idxs=(*w.widxs)[thread_id];
    const auto &mss=*w.pm;
    for(const auto idx:idxs)
        mss[idx]->integrate(w.dt,w.n_combine,0);
}


void msystem::combined_integrate(fast_real dt,int_t n_combine,int_t n_step,int USE_GPU,msystem_combinator *pmc){
    real t_latest=analyse();
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

    for(const auto &c:clist)
        cvecs[c.second].push_back(c.first);

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
                    const ring &mr=*mci.ringmodel;
                    mci.ringmodel=nullptr;
                    fast_real rGM=mr.GM_ratio;
                    rGM/=1-rGM;
                    mci.GM+=mci.GM*rGM;
                    mci.GM0+=mci.GM*rGM;

                    fast_mpvec mgl=mci.GL;
                    fast_mpmat fgls(mgl.asc_node(),0,mgl/mgl.norm());
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
        Sc.build_mid();
        Sc.accel();

        Sx.t_eph=t_eph;
        Sx.build_mid();
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
            sn.build_mid();
            sn.accel();
            sn.p_substep_recorder=pmc;
        }

        if(pmc&&(pmc->pms!=this||pmc->t_split!=t_latest)){
            pmc->pms=this;
            pmc->t_substep=dt;
            pmc->t_split=t_latest;
            pmc->sublists.clear();
            for(auto &sns:Sn)
                pmc->link(sns.second);
        }

        // integrate
#ifndef NDEBUG
        Sc.integrate(dt_long,1,USE_GPU);
        
        Sx.integrate(dt_long,1,0);
        
        for(auto &sns:Sn){
            sns.second.integrate(dt,n_combine,0);
        }
#else
        int_t cost=0,max_cost=0,n_threads;
        std::vector<int_t> costs;
        std::vector<msystem*> msys;
        std::vector<std::vector<size_t>> distributed;
        for(auto &sns:Sn){
            int_t cur_cost=sns.second.mlist.size();
            cur_cost*=cur_cost;
            max_cost=std::max(max_cost,cur_cost);
            cost+=cur_cost;
            msys.push_back(&sns.second);
            costs.push_back(cur_cost);
        }
        n_threads=(cost+max_cost-1)/max_cost;
        distribute(distributed,costs,n_threads);

        thread_works tasks;
        tasks.pm=&msys;
        tasks.widxs=&distributed;
        tasks.dt=dt;
        tasks.n_combine=n_combine;
        ThreadPool *pthread_pool=ThreadPool::get_thread_pool();
        std::vector<std::thread> threads;
        ThreadPool::TaskGroup task_group;
        if(pthread_pool)
            pthread_pool->distribute_tasks(n_threads,do_thread_works,&tasks,&task_group);
        else for(int_t i=0;i<n_threads;++i)
            threads.push_back(std::thread(do_thread_works,&tasks,i));

        Sx.integrate(dt_long,1,0);

        Sc.integrate(dt_long,1,USE_GPU);

        if(pthread_pool)
            pthread_pool->wait_for_all(&task_group);
        else for(int_t i=0;i<n_threads;++i)
            threads[i].join();
#endif
        //finalize

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
            if(mi.w.norm()*std::abs(dt_long)>2*Constants::pi/ROTATION_STEP_LIMIT){
                fast_mpvec zl=mi.s.z*mi.GL;
                fast_real l=fast_mpvec(mi.GL).norm();
                if(l>0){
                    zl/=l;
                    mi.s+=zl.rotation_matrix(fast_real(1))%fast_mpmat(mi.s);
                }
            }
        }

        t_eph=Sc.t_eph;
        update((fast_real)t_eph);
        if(USE_GPU)Cuda_accel();
        else accel();
    }
}


int_t msystem_combinator::link(msystem &mssub,int_t bid){
    const msystem &ms=*pms;
    const auto &blist=ms.get_barycens();

    uint64_t psid=mssub[mssub.tidal_parent].sid;
    std::vector<barycen> &sublist=sublists[psid];
    if(bid<0){
        int_t pmid=ms.get_mid(psid);
        //assume for ms.mlist[i], ms.blist[i].mid==i
        int_t pbid=link(mssub,blist[pmid].tid);
        sublist[pbid].pid=-1;
        barycen::fill_tid(sublist);
        return sublist.size();
    }

    barycen sbi;
    const barycen &bi=blist[bid];
    if(bi.hid>=0){
        const auto &h=blist[bi.hid];
        const auto &g=blist[bi.gid];
        int_t hmid=mssub.get_mid(ms[h.mid].sid);
        int_t gmid=mssub.get_mid(ms[g.mid].sid);
        if(gmid>=0){
            sbi.hid=link(mssub,bi.hid);
            sbi.gid=link(mssub,bi.gid);
        }
    }
    else{
        sbi.mid=mssub.get_mid(ms[bi.mid].sid);
    }

    for(auto cid:bi.children){
        const auto &ci=blist[cid];
        if(mssub.get_mid(ms[ci.mid].sid)>=0)
            sbi.children.push_back(link(mssub,cid));
    }
    int_t sbid=sublist.size();
    if(sbi.hid>=0){
        sublist[sbi.hid].pid=sbid;
        sublist[sbi.gid].pid=sbid;
    }
    for(auto scid:sbi.children)
        sublist[scid].pid=sbid;
    sublist.push_back(sbi);
    return sbid;
}
