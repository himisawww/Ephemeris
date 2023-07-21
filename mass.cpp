#include"ephemeris.h"

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

    fast_mpmat fmis(mi.s);
    fast_mpmat mc(mi.C_static);
    mc.x.x+=mi.exJ2/2;
    mc.y.y+=mi.exJ2/2;
    mc.z.z-=mi.exJ2;
    mc=fmis.toworld(mc);

    mi.phi=0;
    mi.naccel=0;
    mi.C_potential=0;
    if(mi.k2!=0){
        if(mi.tide_delay!=0){
    for(int_t j=0;j<mn;++j)if(this!=&mlist[j]){
        const mass &mj=mlist[j];
        fast_mpvec r=mj.r-mi.r;
        fast_mpvec v=mj.v-mi.v;
        fast_real rr2=1/(r%r);
        fast_real rr=sqrt(rr2);
        fast_real tp_dphi=rr*mj.GM;
        mi.phi-=tp_dphi;
        fast_real tp_dg=rr2*tp_dphi;
        mi.naccel+=tp_dg*r;

        //damped tidal deformation matrix
        fast_real fmj=-3*tp_dg;
        fast_mpvec dw=r*(mi.w*r-v)*rr2;
        fast_real dw2=dw%dw,dw1=sqrt(dw2),dwt=dw1*mi.tide_delay*mj.tide_delay_factor;
        fast_real
            ecwt2=1+4*dwt*dwt,
            secwt2=sqrt(ecwt2),
            x0=sqrt((1+secwt2)/(2*ecwt2)),
            y0=dwt*sqrt(2/(ecwt2*(1+secwt2))),
            z02=dwt*dwt*(2/(ecwt2+secwt2));
        fast_mpvec
            dr=r*x0+dw*r*(y0/dw1);
        mi.C_potential+=(
            fast_mpmat(fast_real(1)/3-z02)
            -fast_mpmat(dr*rr2,dr)
            +fast_mpmat(dw*(z02/dw2),dw)
            )*(fmj*mi.k2);
    }
    }else{//mi.tide_delay==0
    for(int_t j=0;j<mn;++j)if(this!=&mlist[j]){
        const mass &mj=mlist[j];
        fast_mpvec r=mj.r-mi.r;
        fast_mpvec v=mj.v-mi.v;
        fast_real rr2=1/(r%r);
        fast_real rr=sqrt(rr2);
        fast_real tp_dphi=rr*mj.GM;
        mi.phi-=tp_dphi;
        fast_real tp_dg=rr2*tp_dphi;
        mi.naccel+=tp_dg*r;

        //damped tidal deformation matrix
        fast_real fmj=-3*tp_dg;
        fast_mpvec dw=r*(mi.w*r-v)*rr2;
        fast_real dw2=dw%dw,dw1=sqrt(dw2);
        mi.C_potential+=(
            fast_mpmat(fast_real(1)/3)
            -fast_mpmat(r*rr2,r)
            )*(fmj*mi.k2);
    }
    }}else{//mi.k2==0
    for(int_t j=0;j<mn;++j)if(this!=&mlist[j]){
        const mass &mj=mlist[j];
        fast_mpvec r=mj.r-mi.r;
        fast_mpvec v=mj.v-mi.v;
        fast_real rr2=1/(r%r);
        fast_real rr=sqrt(rr2);
        fast_real tp_dphi=rr*mj.GM;
        mi.phi-=tp_dphi;
        fast_real tp_dg=rr2*tp_dphi;
        mi.naccel+=tp_dg*r;
    }
    }
    fast_real w2=mi.w%mi.w;
    mi.C_potential+=fast_mpmat(w2*mi.k2r/3)-fast_mpmat(mi.w*mi.k2r,mi.w);
    mi.C_potential*=mi.R*mi.R2/(2*mi.GM);
    mi.C_potential+=mc;
    mi.GI=2*mi.R2/3*(fast_mpmat(mi.A)-mi.C_potential);
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

    const fast_real c=299792458;
    const fast_real c2=c*c;

    //prepare
    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];

        int_t max_iter=4;
        do{
            mi.deform_this(mlist);
            //angular accelerate
            fast_mpvec oldw=mi.w;
            mi.w=mi.GI.inverse()%(fast_mpvec(mi.GL));
            if((mi.w-oldw).norm()<1e-9*oldw.norm())break;
        } while(--max_iter);

        mi.beta=fast_mpvec(mi.v)/c;
        mi.beta2=mi.beta%mi.beta;
        mi.phi/=c2;
        mi.naccel/=c2;

        mi.daccel=mi.dtorque=mi.gaccel=0;

    }


    for(int_t i=0;i<mn;++i){
        mass &mi=mlist[i];

        for(int_t j=0;j<mn;++j)if(i!=j){
            mass &mj=mlist[j];
            fast_mpvec r=mj.r-mi.r;
            fast_real rr2=1/(r%r);
            fast_real rr=sqrt(rr2);
            fast_real rr3=rr*rr2;

            //start post-newtonian correction
            fast_real tp_dphi=rr*mj.GM;
            fast_real tp_dg=rr2*tp_dphi;
            fast_real rbj=r%mj.beta;
            fast_real tp_rbjrr2=rbj*rr;
            tp_rbjrr2*=tp_rbjrr2;
            fast_real delta1=4*mi.phi+mj.phi+mi.beta2+2*mj.beta2-4*(mi.beta%mj.beta)+(r%mj.naccel-3*tp_rbjrr2)/2;
            fast_mpvec b=mj.beta-mi.beta;
            mi.gaccel+=7*tp_dphi/2*mj.naccel+tp_dg*delta1*r+tp_dg*(rbj-r%b*4)*b;
            //end post-newtonian correction
            mi.min_distance=std::max(mi.min_distance,rr);
            mi.max_influence=std::max(mi.max_influence,tp_dg);
            //rotational & tidal deformation: gravity, torque
            fast_mpvec Cr=mi.C_potential%r;
            fast_mpvec dg=rr3*rr2*mi.R2*(r%Cr*5*rr2*r-(Cr+Cr));
            mi.daccel+=dg*mj.GM;
            mj.daccel-=dg*mi.GM;
            mi.dtorque+=mj.GM*(r*dg);
            //start lense thirring
            fast_real rcr32=2*mj.GM/c*rr3;
            fast_mpvec GL=mi.GL;
            fast_mpvec GLb=GL*b;
            mi.daccel-=rcr32*(GLb-3*(GLb%r)*rr2*r);
            mi.dtorque-=rcr32*(GL*(r*b));
            fast_mpvec GLj=mj.GL;
            mi.daccel+=(rcr32*(GLj-3*(GLj%r)*rr2*r))*b;
            //end lense thirring

            //start radiation pressure
            mi.daccel+=mj.lum*mi.rR2_4Mc*rr3*r;
            //end radiation pressure

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
#if 0
            if(i==10){// lunar geopotential model
                const int_t maxN=6;

                // Jn = J[n-3]
                const fast_real J[]={
                    8.4597026974594570E-06, //J3
                   -9.7044138365700000E-06, //J4
                    7.4221608384052890E-07, //J5
                   -1.3767531350969900E-05  //J6
                };

                // Cij = C[i*(i-1)/2+j-4]
                const fast_real C[]={
                    2.8480741195592860E-05, //C31
                    4.8449420619770600E-06, //C32
                    1.6756178134114570E-06, //C33
                   -5.7048697319733210E-06, //C41
                   -1.5912271792977430E-06, //C42
                   -8.0678881596778210E-08, //C43
                   -1.2692158612216040E-07, //C44
                   -8.6629769308983560E-07, //C51
                    7.1199537967353330E-07, //C52
                    1.5399750424904520E-08, //C53
                    2.1444704319218450E-08, //C54
                    7.6596153884006140E-09, //C55
                    1.2024363601545920E-06, //C61
                   -5.4703897324156850E-07, //C62
                   -6.8785612757292010E-08, //C63
                    1.2915580402925160E-09, //C64
                    1.1737698784460500E-09, //C65
                   -1.0913395178881540E-09  //C66
                };

                // Sij = S[i*(i-1)/2+j-4]
                const fast_real S[]={
                    5.8915551555318640E-06, //S31
                    1.6844743962783900E-06, //S32
                   -2.4742714379805760E-07, //S33
                    1.5789202789245720E-06, //S41
                   -1.5153915796731720E-06, //S42
                   -8.0349266627431070E-07, //S43
                    8.2964257754075220E-08, //S44
                   -3.5272289393243820E-06, //S51
                    1.7107886673430380E-07, //S52
                    2.8736257616334340E-07, //S53
                    5.2652110720146800E-10, //S54
                   -6.7824035473995330E-09, //S55
                   -2.0453507141252220E-06, //S61
                   -2.6966834353574270E-07, //S62
                   -7.1063745295915780E-08, //S63
                   -1.5361616966632300E-08, //S64
                   -8.3465073195142520E-09, //S65
                    1.6844213702632920E-09  //S66
                };
                
                fast_real phi,lambda;
                fast_mpmat fmis(mi.s);
                fast_mpvec lr=fmis.tolocal(r);
                phi=atan2(lr.z,sqrt(lr.x*lr.x+lr.y*lr.y));
                lambda=atan2(lr.y,lr.x);
                fast_real
                    cl=cos(lambda),cp=cos(phi),
                    sl=sin(lambda),sp=sin(phi);

                fast_mpvec
                    ax( cl*cp, sl*cp,sp),
                    ay(-   sl,    cl, 0),
                    az(-cl*sp,-sl*sp,cp);

                fast_real R_r=mi.R*rr;
                fast_mpvec an(0);
                for(int_t n=maxN;n>=3;--n){
                    const fast_real &Jn=J[n-3];
                    fast_real sgncp=std::copysign(fast_real(1),cp);

                    an.x+=Jn*(n+1)*std::legendre(n,sp);
                    an.z-=Jn*sgncp*std::assoc_legendre(n,1,sp);
                    for(int_t m=1;m<=n;++m){
                        const fast_real &Cnm=C[n*(n-1)/2+m-4];
                        const fast_real &Snm=S[n*(n-1)/2+m-4];
                        fast_real cml=cos(m*lambda),sml=sin(m*lambda);
                        an.x-=(n+1)*std::assoc_legendre(n,m,sp)*(Cnm*cml+Snm*sml);
                        fast_real anyp;
                        if(abs(cp)>abs(sp)){
                            anyp=std::assoc_legendre(n,m,sp)/cp;
                        }
                        else{
                            anyp=sgncp*(
                                std::assoc_legendre(n,m+1,sp)
                               +(n+m)*(n-m+1)*std::assoc_legendre(n,m-1,sp)
                                )/(2*m*sp);
                        }
                        an.y+=m*anyp*(-Cnm*sml+Snm*cml);
                        an.z+=(m*sp*anyp-(n+m)*(n-m+1)*sgncp*std::assoc_legendre(n,m-1,sp))*(Cnm*cml+Snm*sml);
                    }
                    an*=R_r;
                }
                an*=rr2*R_r*R_r;
                an=fmis.toworld(an.x*ax+an.y*ay+an.z*az);
                mi.daccel-=mj.GM*an;
                mj.daccel+=mi.GM*an;
                mi.dtorque-=mj.GM*(r*an);
            }
#else
            if(mi.gpmodel){
                fast_mpmat fmis(mi.s);
                fast_mpvec lr=fmis.tolocal(r);
                fast_mpvec an=fmis.toworld(mi.gpmodel->sum(mi.R,lr));
                mi.daccel-=mj.GM*an;
                mj.daccel+=mi.GM*an;
                mi.dtorque-=mj.GM*(r*an);
            }
#endif
            //start ring gravity model
            if(mi.ringmodel){
                fast_mpvec migl=mi.GL;
                fast_mpmat fgls(migl.perpunit(),0,migl/migl.norm());
                fast_mpvec lr=fgls.tolocal(r);
                fast_mpvec an=fgls.toworld(mi.ringmodel->sum(lr));
                mi.daccel-=mj.GM*an;
                mj.daccel+=mi.GM*an;
                mi.dtorque-=mj.GM*(r*an);
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
        mi.phi*=c2;
        mi.naccel*=c2;

        if(mi.ringmodel){
            mass &m=mi;
            ring &mr=*m.ringmodel;
            //ring has inertia that prevents extra accelerations
            m.daccel-=(m.gaccel+m.daccel+m.naccel)*mr.GM_ratio;
            //ring has angular momentum that prevents extra angular accelerations
            fast_mpvec mGL=m.GL;
            fast_real rmGL2=1/(mGL%mGL),rmGL=sqrt(rmGL2);
            fast_mpvec ptorque=m.dtorque;
            ptorque=ptorque%mGL*rmGL2*mGL;
            fast_mpvec dtorque=m.dtorque-ptorque;
            //ptorque(parallel to GL) will not change
            //  since this part has nothing to do with the ring
            //dtorque(perpendicular to GL) will decrease
            //  due to ring's angular momentum
            dtorque*=1/(1+rmGL*mr.GL);
            m.dtorque=dtorque+ptorque;
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
