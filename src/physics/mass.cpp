#include"physics/CelestialSystem.h"
#include"physics/geopotential.h"
#include"physics/ring.h"
#include"physics/mass.h"

void mass::scale(fast_real factor){

    fast_real f2=factor*factor,f3=factor*f2;

    R*=factor;

    R2*=f2;
    GL*=f2;
    GI*=f2;

    GM*=f3;
    GM0*=f3;
    dGM*=f3;

    lum*=f2;        //luminosity scales with surface
    recpt*=factor;  //receptance scales to not influence acceleration

    if(ringmodel){
        ringmodel->GL*=f2;
        for(int_t i=0;i<ringmodel->N;++i){
            ringmodel->c_table[i].Gs/=f2;
            ringmodel->c_table[i].R*=factor;
            ringmodel->c_table[i].H*=factor;
        }
    }
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
void mass::deform_this(const std::vector<mass> &mlist){
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
void mass::deform_all(std::vector<mass> &mlist){
    int_t mn=mlist.size();
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];
        mi.deform_this(mlist);
    }
}
void msystem::deform(){
    mass::deform_all(mlist);
}
void msystem::accel(){
    int_t mn=mlist.size();

    using Constants::c;
    using Constants::c2;

    //prepare
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];

        int_t max_iter=MAX_ANGULAR_VELOCITY_ITER;
        do{
            mi.deform_this(mlist);
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
            mi.min_distance=std::max(mi.min_distance,rr);
            mi.max_influence=std::max(mi.max_influence,tp_dg);
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
                fast_mpmat fgls(migl.perpunit(),0,migl/migl.norm());
                fast_mpvec lr=fgls.tolocal(r);
                fast_mpvec an=fgls.toworld(mi.ringmodel->sum(lr));
                APPLY_NONPOINT_FORCE(mi);
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

int_t msystem::get_mid(const char *sid){
    size_t slen=strlen(sid);
    uint64_t isid=0;
    const size_t maxslen=sizeof(isid)-1;
    if(slen>maxslen)return -1;
    memcpy(&isid,sid,slen);
    return get_mid(isid);
}
int_t msystem::get_mid(uint64_t sid){
    auto it=midx.find(sid);

    do{
        if(it==midx.end())break;
        size_t result=it->second;
        if(result>=mlist.size())break;
        if(mlist[result].sid!=sid)break;
        //correct
        return result;
    } while(0);

    //not correct, reconstruct midx
    midx.clear();
    size_t mn=mlist.size();
    for(size_t i=0;i<mn;++i){
        const mass &mi=mlist[i];
        midx.insert({mi.sid,i});
    }

    it=midx.find(sid);

    if(it==midx.end())
        return -1;

    return it->second;
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
            fprintf(stderr,
                "\n\n************************************************************************\n"
                "Catastrophy: At Ephemeris Time: %lld s\n"
                "   Collision of Celestial Bodies Detected !!!\n"
                "   Another body entered radius of <%s>:\n"
                "       %.16e m <  %.16le m\n"
                "   Note the program is not designed to handle body collisions.\n"
                "   Ephemeris further than this may be unreliable.\n"
                "************************************************************************\n\n",
               int_t(t_eph.hi)+int_t(t_eph.lo),(char*)&m.sid,m.min_distance,m.R);
        }
        if(m.max_influence*dt2>TIMESTEP_THRESHOLD){
            fprintf(stderr,
                "\n\n************************************************************************\n"
                "Warning: At Ephemeris Time: %lld s\n"
                "   Time step is too large for accurate integration of <%s>.\n"
                "       %.16e s >  %.16le s\n"
                "   Ephemeris further than this may be inaccurate.\n"
                "************************************************************************\n\n",
                int_t(t_eph.hi)+int_t(t_eph.lo),(char*)&m.sid,dt,std::sqrt(TIMESTEP_THRESHOLD/m.max_influence));
        }
    }

    t_eph+=dt*n_step;
}
