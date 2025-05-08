#include"physics/mass.h"
#include"physics/geopotential.h"
#include"physics/ring.h"
#include"physics/mass.impl"
#include"utils/logger.h"
#include"utils/threadpool.h"

void mass::scale(fast_real factor){

    fast_real f2=factor*factor,f3=factor*f2;

    R*=factor;

    R2*=f2;
    sR2*=f2;
    GL*=f2;
    GI*=f2;

    GM*=f3;
    GM0*=f3;
    dGM*=f3;

    lum*=f2;        //luminosity scales with surface
    recpt*=factor;  //receptance scales to not influence acceleration
    rR2G_4c*=f3;
}

void msystem::update(fast_real t){
    int_t mn=mlist.size();
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];

        mi.GM=mi.GM0+mi.dGM*t;

        /*if(i==3){//earth rotation adjust
            fast_real targetw=2*pi/(86164.10013880364+0.00172/(86400*365.25*100)*t);
            fast_real wn=mi.w.norm();
            fast_real dw=targetw-wn;

            fast_mpmat fmis(mi.s);
            fast_mpmat AmC=fmis.tolocal(mi.GI/(2*mi.R2/3));
            fast_mpmat dc(-1);
            dc.z.z=2;
            mi.C_static+=dc*AmC.z.z*dw/(2*wn);
        }*/
        mi.exJ2=mi.dJ2*t;

    }
}
void mass::deform_by(const std::vector<mass> &mlist){
    int_t mn=mlist.size();
    mass &mi=*this;

    mi.phi=0;
    mi.naccel=0;
    mi.C_potential=0;
    if(mi.k2!=0){
        if(mi.tide_delay!=0){
    for(int_t j=0;j<mn;++j)if(this!=&mlist[j]){
        const mass &mj=mlist[j];
        DAMPED_TIDAL_DEFORMATION_MATRIX(mi);
    }
    }else{//mi.tide_delay==0
    for(int_t j=0;j<mn;++j)if(this!=&mlist[j]){
        const mass &mj=mlist[j];
        DAMPED_TIDAL_DEFORMATION_MATRIX_NDELAY(mi);
    }
    }}else{//mi.k2==0
    for(int_t j=0;j<mn;++j)if(this!=&mlist[j]){
        const mass &mj=mlist[j];
        DAMPED_TIDAL_DEFORMATION_MATRIX_NK2(mi);
    }
    }
    UPDATE_HARMONICS;
}
void msystem::deform(){
    int_t mn=mlist.size();
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];
        mi.deform_by(mlist);
    }
}

static void accel_deform(void *param,size_t){
    using Constants::c;
    using Constants::c2;

    mass &mi=*(mass*)param;
    auto &mlist=*(std::vector<mass>*)mi.pmlist;
    int_t i_start=&mi-mlist.data();
    int_t i_end=mi.task_index;
    for(int_t i=i_start;i<i_end;++i){
        mass &mi=mlist[i];

        int_t max_iter=MAX_ANGULAR_VELOCITY_ITER;
        do{
            mi.deform_by(mlist);
            bool should_break;
            UPDATE_ANGULAR_VELOCITY;
            if(should_break)break;
        } while(--max_iter);

        PREPARE_RELATIVITY;

        mi.daccel=mi.dtorque=mi.gaccel=0;
    }
}
static void accel_mainforce(void *param,size_t){
    using Constants::c;

    mass &mi=*(mass*)param;
    auto &mlist=*(std::vector<mass>*)mi.pmlist;
    int_t mn=mlist.size();
    int_t i_start=&mi-mlist.data();
    int_t i_end=mi.task_index;
    for(int_t i=i_start;i<i_end;++i){
        mass &mi=mlist[i];

        for(int_t j=0;j<mn;++j)if(i!=j){
            mass &mj=mlist[j];
            RELATIVITY(mi);
            checked_maximize(mi.min_distance,rr);
            checked_maximize(mi.max_influence,tp_dg);
            ROTATIONAL_TIDAL_DEFORMATION_NANTI_FORCE(mi);
            LENSE_THIRRING(mi);
            RADIATION_PRESSURE(mi);
        }
    }
}
static void accel_nonpoint(void *param,size_t){
    auto end_task=*(uintptr_t*)param;
    do{
        if(end_task&1)break;
        end_task=*(uintptr_t*)end_task;
    } while(1);
    auto &mlist=*(std::vector<mass>*)(end_task-1);
    int_t mn=mlist.size();
    int_t i=((char*)param-(char*)mlist.data())/sizeof(mass);
    mass &minit=mlist[i];
    bool is_geopotential=param==&minit.pmlist;
    int_t k=is_geopotential?minit.task_index:minit.task_jndex;
    int_t n_task=minit.task_count;
    if(k>=n_task)k-=n_task;
    int_t j_start=k*mn/n_task;
    int_t j_end=(k+1)*mn/n_task;
    void *next_task;
    do{
        mass &mi=mlist[i];
        struct{
            fast_mpvec daccel,dtorque;
        } tpmi;
        tpmi.daccel=tpmi.dtorque=0;
        if(is_geopotential){
            for(int_t j=j_start;j<j_end;++j)if(i!=j){
                mass &mj=mlist[j];
                fast_mpvec r=mj.r-mi.r;

                fast_mpmat fmis(mi.s);
                fast_mpvec lr=fmis.tolocal(r);
                fast_mpvec an=fmis.toworld(mi.gpmodel->sum(mi.R,lr));
                APPLY_NONPOINT_FORCE(tpmi);
            }
            mi.idaccel+=tpmi.daccel;
            mi.idtorque+=tpmi.dtorque;
            next_task=mi.pmlist;
        }
        else{
            for(int_t j=j_start;j<j_end;++j)if(i!=j){
                mass &mj=mlist[j];
                fast_mpvec r=mj.r-mi.r;

                fast_mpvec migl=mi.GL;
                fast_mpmat fgls(migl.asc_node(),0,migl.unit());
                fast_mpvec lr=fgls.tolocal(r);
                fast_mpvec an=fgls.toworld(mi.ringmodel->sum(mi.R,lr));
                APPLY_NONPOINT_FORCE(tpmi);
            }
            mi.jdaccel+=tpmi.daccel;
            mi.jdtorque+=tpmi.dtorque;
            next_task=mi.qmlist;
        }
        if(uintptr_t(next_task)==end_task)break;
        i=((char*)next_task-(char*)mlist.data())/sizeof(mass);
        is_geopotential=next_task==&mlist[i].pmlist;
    } while(1);
}

void msystem::accel(int parallel_option){
    int_t mn=mlist.size();
    if(parallel_option>0){
        Cuda_accel();
        return;
    }

    using Constants::c;
    using Constants::c2;

    constexpr int_t n_tasklimit=64;
    constexpr int_t n_tasksize=48*48;
    int_t n_task=parallel_option|std::min(mn,(mn*mn+(n_tasksize-1))/n_tasksize);
    ThreadPool *pthreadpool;
  if(n_task>1&&(pthreadpool=ThreadPool::get_thread_pool())){
    n_task=std::min(n_task,n_tasklimit);
    ThreadPool::TaskGroup g;
    for(int_t i_task=0;i_task<n_task;++i_task){
        int_t i_start=i_task*mn/n_task;
        int_t i_end=(i_task+1)*mn/n_task;
        mass &mi=mlist[i_start];
        mi.pmlist=&mlist;
        mi.task_index=i_end;
        pthreadpool->add_task(accel_deform,&mi,&g);
    }
    pthreadpool->wait_for_all(&g);

    for(int_t i_task=0;i_task<n_task;++i_task){
        int_t i_start=i_task*mn/n_task;
        int_t i_end=(i_task+1)*mn/n_task;
        mass &mi=mlist[i_start];
        mi.pmlist=&mlist;
        mi.task_index=i_end;
        pthreadpool->add_task(accel_mainforce,&mi,&g);
    }
    pthreadpool->wait_for_all(&g);

    int_t n_nonpoint=0;
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];
        if(mi.gpmodel){
            ++n_nonpoint;
            mi.idaccel=mi.idtorque=0;
        }
        if(mi.ringmodel){
            ++n_nonpoint;
            mi.jdaccel=mi.jdtorque=0;
        }
    }

    n_task=std::min(n_nonpoint,n_task);
    void *const null_task=(char*)&mlist+1;
    void *next_task=null_task;
    for(int_t k=0;k<n_task;++k){
        int_t j=0,j_start=0;
        for(int_t i=mn;i>0;){--i;
            mass &mi=mlist[i];
            if(!mi.gpmodel&&!mi.ringmodel)continue;
            mi.task_count=n_task;
            if(mi.gpmodel){
                mi.pmlist=next_task;
                next_task=&mi.pmlist;
                if(++j==(j_start+1)*n_nonpoint/n_task){
                    mi.task_index=j_start+k;
                    pthreadpool->add_task(accel_nonpoint,next_task,&g);
                    ++j_start;
                    next_task=null_task;
                }
            }
            if(mi.ringmodel){
                mi.qmlist=next_task;
                next_task=&mi.qmlist;
                if(++j==(j_start+1)*n_nonpoint/n_task){
                    mi.task_jndex=j_start+k;
                    pthreadpool->add_task(accel_nonpoint,next_task,&g);
                    ++j_start;
                    next_task=null_task;
                }
            }
        }
        pthreadpool->wait_for_all(&g);
    }
    if(n_nonpoint)for(auto &mi:mlist){
        if(mi.gpmodel){
            mi.daccel+=mi.idaccel;
            mi.dtorque+=mi.idtorque;
        }
        if(mi.ringmodel){
            mi.daccel+=mi.jdaccel;
            mi.dtorque+=mi.jdtorque;
        }
    }
  }
  else{
    //prepare
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];

        int_t max_iter=MAX_ANGULAR_VELOCITY_ITER;
        do{
            mi.deform_by(mlist);
            bool should_break;
            UPDATE_ANGULAR_VELOCITY;
            if(should_break)break;
        } while(--max_iter);

        PREPARE_RELATIVITY;

        mi.daccel=mi.dtorque=mi.gaccel=0;

    }


    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];

        for(int_t j=0;j<mn;++j)if(i!=j){
            mass &mj=mlist[j];
            RELATIVITY(mi);
            checked_maximize(mi.min_distance,rr);
            checked_maximize(mi.max_influence,tp_dg);
            ROTATIONAL_TIDAL_DEFORMATION;
            LENSE_THIRRING(mi);
            RADIATION_PRESSURE(mi);

            /*start shape-shape torque, 17cm for lunar libration
            if(i==10&&j==3){
                fast_mpmat es(mj.s);
                fast_mpmat eGI(es.tolocal(mj.GI));
                fast_real J2=(eGI.z.z-(eGI.x.x+eGI.y.y)/2)/(mj.R2);
                //printf("%le\n",J2);
                fast_mpvec rem=r*rr,pe=es.z;
                fast_mpvec irem=mi.GI%rem,ipe=mi.GI%pe;
                fast_real sinphi=rem%pe;
                fast_mpvec sstorque=(fast_real)15/2*mj.GM*mj.R2*J2*(rr3*rr2)*(
                    (1-7*sinphi*sinphi)*(rem*(irem))
                    +2*sinphi*(rem*(ipe)+pe*(irem))
                    -(fast_real)2/5*(pe*(ipe))
                    );
                //printf("%le ",sstorque.norm());
                mi.dtorque+=sstorque;*/
            //start higher harmonics
            if(mi.gpmodel){
                fast_mpmat fmis(mi.s);
                fast_mpvec lr=fmis.tolocal(r);
                fast_mpvec an=fmis.toworld(mi.gpmodel->sum(mi.R,lr));
                APPLY_NONPOINT_FORCE(mi);
            }
            //start ring gravity model
            if(mi.ringmodel){
                fast_mpvec migl=mi.GL;
                fast_mpmat fgls(migl.asc_node(),0,migl.unit());
                fast_mpvec lr=fgls.tolocal(r);
                fast_mpvec an=fgls.toworld(mi.ringmodel->sum(mi.R,lr));
                APPLY_NONPOINT_FORCE(mi);
            }
        }
    }
  }

    //tidal matrix correction
    //Not implemented in GPU version
    if(tidal_childlist.size()){
        mass &mp=mlist[tidal_parent];
        for(int_t cidx:tidal_childlist){
            mass &mc=mlist[cidx];
            fast_mpvec r=mc.r-mp.r;
            fast_mpvec dg=(tidal_matrix%r)/(mc.GM+mp.GM);
            mc.daccel+=mp.GM*dg;
            mp.daccel-=mc.GM*dg;
        }
    }
    
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];
        FINALIZE_RELATIVITY;

        if(mi.ringmodel){
            RING_CORRECTION;
        }
    }
}

void msystem::integrate(fast_real dt,int_t n_step,int USE_GPU){
    if(n_step==0)return;
    if(n_step<0){
        n_step=-n_step;
        dt=-dt;
    }
    
    for(auto &m:mlist)m.max_influence=m.min_distance=0;
    
    if(!USE_GPU)RungeKutta12(dt,n_step);
    else Cuda_RungeKutta12(dt,n_step);
    
    //square of maximum angular displacement per step for circular orbit
    const double TIMESTEP_THRESHOLD=0.01;
    double dt2=dt*dt;

    for(auto &m:mlist){
        m.min_distance=1/m.min_distance;
        if(m.R>m.min_distance){
            LogError(
                "\n\n************************************************************************\n"
                "Catastrophy: At Ephemeris Time: %lld s\n"
                "   Collision of Celestial Bodies Detected !!!\n"
                "   Another body entered radius of <%s>:\n"
                "       %.16e m <  %.16le m\n"
                "   Note the program is not designed to handle body collisions.\n"
                "   Ephemeris further than this may be unreliable.\n"
                "************************************************************************\n\n",
                ephemeris_time(),(char*)&m.sid,m.min_distance,m.R);
        }
        if(m.max_influence*dt2>TIMESTEP_THRESHOLD){
            LogWarning(
                "\n\n************************************************************************\n"
                "Warning: At Ephemeris Time: %lld s\n"
                "   Time step is too large for accurate integration of <%s>.\n"
                "       %.16e s >  %.16le s\n"
                "   Ephemeris further than this may be inaccurate.\n"
                "************************************************************************\n\n",
                ephemeris_time(),(char*)&m.sid,dt,std::sqrt(TIMESTEP_THRESHOLD/m.max_influence));
        }
    }

    t_eph+=dt*n_step;
}

void msystem::clear_accel(){
    for(mass &m:mlist){
        char *pstart=(char*)&(m.*mass_auxiliary_variables);
        char *pend=(char*)&m+sizeof(m);
        memset(pstart,-1,pend-pstart);
    }
}

int_t msystem::ephemeris_time() const{
    return int_t(t_eph.hi)+int_t(t_eph.lo);
}
