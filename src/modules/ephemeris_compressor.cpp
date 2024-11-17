#include"ephemeris_compressor.h"
#include<map>
#include"math/interp.h"
#include"utils/logger.h"

#define HALF 0.5
typedef ephemeris_compressor::keplerian kep_t;

kep_t::keplerian(const ephem_orb &other):ephem_orb(other){
    el=earg+j.asc_node().phi();
    ex=(q+1)*std::cos(el);
    ey=(q+1)*std::sin(el);
    if(q>=0)
        ml=NAN;
    else{
        double s=-q*(q+2);
        s*=std::sqrt(s);
        ml=el+m*s;
        //ml=angle_reduce(ml);
    }
}
kep_t::keplerian(const double *p,bool circular){
    t=1;
    j.x=p[0];
    j.y=p[1];
    j.z=p[2];
    if(circular){
        ex=p[3];
        ey=p[4];
        ml=p[5];
    }
    else{
        q=p[3];
        el=p[4];
        m=p[5];
    }
    blend_finalize(circular);
}

bool kep_t::interpolatable(bool circular,const kep_t &k1,const kep_t &k2){
    if(circular){
        double mq=checked_max(k1.q,k2.q);
        if(!(mq<0))return false;
        double rdq=std::abs(k1.q-k2.q)/mq;
        if(!(rdq<1))return false;
    }
    else{
        if(!(k2.m>k1.m))return false;
        double del=angle_reduce(k2.el-k1.el);
        if(!(del<1))return false;
    }
    return (k2.j-k1.j).normsqr()/checked_min(k1.j.normsqr(),k2.j.normsqr())<1;
}

void kep_t::blend_initialize(bool circular){
    t=0;//w
    j=0;
    if(circular){
        ex=0;
        ey=0;
        ml=0;
        earg=NAN;//last_ml;
    }
    else{
        q=0;
        el=0;
        m=0;
        earg=NAN;//last_el;
    }
}
void kep_t::blend_add(bool circular,const keplerian &ki,double ws){
    double &w=t;
    if(circular){
        double &last_ml=earg;
        w+=ws;
        j+=ws*ki.j;
        ex+=ws*ki.ex;
        ey+=ws*ki.ey;
        double mli=ki.ml;
        if(last_ml==last_ml)mli=last_ml+angle_reduce(mli-last_ml);
        ml+=ws*mli;
        last_ml=mli;
    }
    else{
        double &last_el=earg;
        w+=ws;
        j+=ws*ki.j;
        q+=ws*ki.q;
        m+=ws*ki.m;
        double eli=ki.el;
        if(last_el==last_el)eli=last_el+angle_reduce(eli-last_el);
        el+=ws*eli;
        last_el=eli;
    }
}
void kep_t::blend_finalize(bool circular){
    double &w=t;
    w=1/w;
    j*=w;
    if(circular){
        ex*=w;
        ey*=w;
        ml*=w;
        q=std::sqrt(ex*ex+ey*ey)-1;
        el=std::atan2(ey,ex);
        double s=-q*(q+2);
        s*=std::sqrt(s);
        m=(ml-el)/s;
        earg=el-j.asc_node().phi();
    }
    else{
        q*=w;
        m*=w;
        el*=w;
        earg=el-j.asc_node().phi();
        ex=(q+1)*std::cos(el);
        ey=(q+1)*std::sin(el);
        if(q>=0)
            ml=NAN;
        else{
            double s=-q*(q+2);
            s*=std::sqrt(s);
            ml=el+m*s;
        }
    }
    w=0;
}

ephemeris_compressor::axial::axial(const double *a){
    t=0;
    w.x=a[0];
    w.y=a[1];
    w.z=a[2];
    pphi=std::asin(a[3]/std::sqrt(1-a[4]*a[4]));
    ptheta=std::acos(a[4]);
    angle=a[5];
}

double ephemeris_compressor::infer_GM_from_data(const orbital_state_t *pdata,int_t N){
    std::vector<double> angles(N);
    angles[0]=0;
    for(int_t i=1;i<N;++i){
        const vec &rm=pdata[i-1].r;
        const vec &rp=pdata[i].r;
        vec rdiff=rp-rm;
        angles[i]=angles[i-1]+std::sqrt(rdiff%rdiff/std::min(rm%rm,rp%rp));
    }

    constexpr double min_angle=Constants::degree*15;
    double musum=0;
    double muweight=0;
    int_t muerr=0;
    for(int_t l=0,r=0;l<N&&r<N;){
        int_t lold=l,rold=r;
        while(r+1<N&&angles[r]-angles[l]<min_angle)++r;
        if(l<r){
            const vec &r0=pdata[l].r;
            const vec &v0=pdata[l].v;
            const vec &r1=pdata[r].r;
            const vec &v1=pdata[r].v;
            vec drn=r0/r0.norm()-r1/r1.norm();
            vec dvj=v0*(r0*v0)-v1*(r1*v1);
            double dweight=drn.normalize();
            double dmu=dvj%drn;
            if(dmu*2<dvj.norm())
                ++muerr;

            musum+=dmu;
            muweight+=dweight;
        }
        if(rold==r)++r;
        if(lold==l)++l;
    }

    return muerr?0:musum/muweight;
}

static constexpr double midinterp_coefficients[][8]={
    {1./2},
    {2./3,-1./6},
    {3./7,3./14,-1./7},
    {27./50,3./25,-11./50,3./50},
    {60./143,30./143,-5./143,-45./286,9./143},
    {1000./2009,325./2009,-250./2009,-515./4018,242./2009,-55./2009},
    {1750./4199,875./4199,-70./4199,-1085./8398,-14./247,35./323,-10./323},
    {245./514,91./514,-133./1542,-35./257,7./514,175./1542,-37./514,7./514}
};
constexpr int_t midinterp_max_wing=sizeof(midinterp_coefficients)/sizeof(midinterp_coefficients[0]);

double ephemeris_compressor::compression_score(double relative_error,double compressed_size){
    //score acts like (relative error)^[this index] * size at epsilon_relative_error
    static constexpr double index_compress=0.25;
    //make index_compress vary as (relative error)^[this index] to suppress large relative error
    static constexpr double index_accuracy=0.3;

    const double a=index_compress/std::pow(epsilon_relative_error,index_accuracy);
    return std::log(compressed_size)/a+std::pow(relative_error,index_accuracy)/index_accuracy;
}

template<ephemeris_format F,typename T,size_t N_Channel>
MFILE ephemeris_compressor::compress_data(const T *pdata,int_t N,int_t d){
    typedef header_t<F> data_header;

    struct scored_compress_data{
        MFILE data;
        data_header *pheader;
        double max_error;
        double reduced_error;
        double score;
    };
    std::map<int_t,scored_compress_data> mfmap;
    
    auto make_fit=[&](int_t n){
        auto &target=mfmap[n];
        MFILE &mf=target.data;
        const size_t bsp_size=N_Channel*sizeof(T)*(n+d);
        const size_t data_size=bsp_size+sizeof(data_header);
        bspline_fitter<T,N_Channel> bf(d,n,pdata,N-1);
        target.pheader=(data_header*)mf.prepare(data_size);
        memcpy(target.pheader+1,bf.get_fitted_data(),bsp_size);
        bf.expand();
        T result[N_Channel];
        double max_error=0;
        bf(double(0),target.pheader->fix[0],target.pheader->fix[1]);
        for(int_t i=0;i<N;++i){
            bf(double(i),result);
            double ierror=data_header::error_function(pdata+N_Channel*i,result);
            checked_maximize(max_error,ierror);
        }
        bf(double(N-1),target.pheader->fix[2],target.pheader->fix[3]);
        target.max_error=max_error;
        double score=filtered_min<double>(max_error,INFINITY)+epsilon_relative_error;
        target.reduced_error=score;
        target.score=compression_score(score,double(data_size));
        return &target;
    };

    std::vector<int_t> n_all,i_chs;
    segment_choices(n_all,N,d);
    const int_t n_min=n_all.back();
    const int_t n_max=n_all.front();

    make_fit(n_max);
    int_t i_current=0,i_left=0,i_right=n_all.size()-1;
    int_t step=1,dir=0,n_optimal=n_max;
    n_all.front()=-1;

    for(int_t k=0;k<2;++k){
        int_t next_dir=dir^k;
        i_chs.clear();
        int_t i_next=i_current;
        do{
            i_next+=next_dir==0?-1:1;
            if(i_next<i_left||i_next>i_right)
                break;
            if(n_all[i_next]<n_min)
                continue;
            i_chs.push_back(i_next);
        } while(1);
        if(i_chs.empty())continue;
        step=dir==next_dir?step*2:std::max<int_t>(1,step/2);
        step=std::min<int_t>(step,(i_chs.size()+1)/2);
        dir=next_dir;
        i_current=i_chs[step-1];
        int_t &n_chs=n_all[i_current];
        if(make_fit(n_chs)->score<mfmap[n_optimal].score){
            n_optimal=n_chs;
            i_left=i_current;
            while(i_left>0&&n_all[i_left-1]>=n_min)--i_left;
            int_t i_max=n_all.size();
            i_right=i_current;
            while(i_right+1<i_max&&n_all[i_right+1]>=n_min)++i_right;
        }
        else (n_chs<n_optimal?i_right:i_left)=i_current;

        n_chs=-1;
        k=-1;
    }

    scored_compress_data &result=mfmap[n_optimal];
    //fill header
    data_header &header=*result.pheader;
    header.relative_error=(float)result.max_error;
    header.degree=(int16_t)d;
    header.uformat=~(uint16_t)F;
    header.n=n_optimal;
    if(!(header.relative_error<1))
        result.data.reset();
    return std::move(result.data);
}

double ephemeris_compressor::relative_state_error(const vec *r,const vec *rp){
    return std::sqrt((*rp-*r).normsqr()/(*r).normsqr());
};
double ephemeris_compressor::absolute_state_error(const vec *r,const vec *rp){
    constexpr double minref=epsilon_absolute_error/epsilon_relative_error;
    return std::sqrt((*rp-*r).normsqr()/std::max(minref*minref,(*r).normsqr()));
};
double ephemeris_compressor::circular_kepler_error(const double *k,const double *kp){
    vec r,rp,v;
    kep_t(k,true).rv(0,r,v);
    kep_t(kp,true).rv(0,rp,v);
    return std::sqrt((rp-r).normsqr()/r.normsqr());
}
double ephemeris_compressor::raw_kepler_error(const double *k,const double *kp){
    vec r,rp,v;
    kep_t(k,false).rv(0,r,v);
    kep_t(kp,false).rv(0,rp,v);
    return std::sqrt((rp-r).normsqr()/r.normsqr());
}
double ephemeris_compressor::axial_rotation_error(const double *a,const double *ap){
    mat s,sp;
    vec w;
    axial(a).sw(0,s,w);
    axial(ap).sw(0,sp,w);
    quat q(s),qp(sp);
    return quaterion_rotation_error(&q,&qp);
}
double ephemeris_compressor::quaterion_rotation_error(const quat *q,const quat *qp){
    return std::sqrt(checked_min((*q-*qp).normsqr(),(*q+*qp).normsqr()));
}

int_t ephemeris_compressor::compress_orbital_data(MFILE &mf,double time_span){
    mf.load_data();
    int_t N=mf.size();
    const orbital_state_t *pdata=(const orbital_state_t*)mf.prepare(N);
    constexpr int_t state_size=sizeof(orbital_state_t);
    if(N<2*state_size||N%state_size){
        LogError("compress_orbital_data::invalid file size.\n");
        return 0;
    }
    N/=state_size;

    double delta_t=time_span/(N-1);
    //DEBUG
    //N=2;

    int_t d=std::min(max_bspline_degree,N-1);
    d-=d+1&1;
    int_t dtest=std::min(d,N-2);
    dtest-=dtest+1&1;

    struct fit_err_t{
        //relative errors of orbital formats
        double serr,kcerr,krerr;
    };
    //fill NAN
    std::vector<fit_err_t> fiterrs(N);
    memset(fiterrs.data(),-1,N*sizeof(fit_err_t));

    std::vector<kep_t> ostates;
    double GM=infer_GM_from_data(pdata,N);
    const double tfac=std::sqrt(GM);
    const double vfac=1/tfac;
    const double h=delta_t*tfac;
    if(GM>=Constants::G){
        ostates.reserve(N);
        for(int_t i=0;i<N;++i)
            ostates.emplace_back(ephem_orb(0,pdata[i].r,pdata[i].v*vfac));
    }
    bool can_kepler=!ostates.empty();

    auto *state_error_fun=can_kepler?relative_state_error:absolute_state_error;

    if(dtest<1){//N==2 here
        if(can_kepler){
            vec r0,v0,r1,v1;
            kep_t kmix;
            for(int k=0;k<2;++k){
                const bool circ=k==0;
                if(!kep_t::interpolatable(circ,ostates[0],ostates[1]))
                    continue;
                kmix.blend_initialize(circ);
                kmix.blend_add(circ,ostates[0],HALF);
                kmix.blend_add(circ,ostates[1],HALF);
                kmix.blend_finalize(circ);
                kmix.rv(-HALF*h,r0,v0);
                kmix.rv(+HALF*h,r1,v1);
                double kfac=1/(circ?-kmix.q:1+kmix.q);
                const auto perr=circ?&fit_err_t::kcerr:&fit_err_t::krerr;
                fiterrs[0].*perr=kfac*state_error_fun(&pdata[0].r,&r0);
                fiterrs[1].*perr=kfac*state_error_fun(&pdata[1].r,&r1);
            }
        }
        vec dr=pdata[1].r-pdata[0].r;
        dr-=HALF*delta_t*(pdata[0].v+pdata[1].v);
        double rerr=HALF*std::sqrt(dr.normsqr()/checked_min(pdata[0].r.normsqr(),pdata[1].r.normsqr()));
        fiterrs[1].serr=fiterrs[0].serr=rerr;
    }
    else{
        int_t wing=(dtest+1)/2;
        if(wing>midinterp_max_wing){
            LogError("Error: Interpolation Order %d Not Implemented.\n",dtest);
            return 0;
        }
        const double *midinterp_coefs=midinterp_coefficients[wing-1];

        std::vector<bool> can_interp_c(N,can_kepler),can_interp_r(N,can_kepler);
        if(can_kepler){
            for(int_t i=0;i+1<N;++i){
                can_interp_c[i]=kep_t::interpolatable(1,ostates[i],ostates[i+1]);
                can_interp_r[i]=kep_t::interpolatable(0,ostates[i],ostates[i+1]);
            }
            for(int_t i=N-1;i>0;--i){
                can_interp_c[i]=can_interp_c[i]&&can_interp_c[i-1];
                can_interp_r[i]=can_interp_r[i]&&can_interp_r[i-1];
            }
        }

        for(int_t i=wing;i+wing<N;++i){
            bool do_interp_c=can_kepler,do_interp_r=can_kepler;
            for(int_t j=i-wing;do_interp_c&&j<=i+wing;++j)
                do_interp_c=do_interp_c&&can_interp_c[j];
            for(int_t j=i-wing;do_interp_r&&j<=i+wing;++j)
                do_interp_r=do_interp_r&&can_interp_r[j];

            vec r,v;
            kep_t kmix;
            for(int k=0;k<2;++k){
                const bool circ=k==0;
                if(!(circ?do_interp_c:do_interp_r))continue;
                kmix.blend_initialize(circ);
                for(int_t j=-wing;j<=wing;++j){
                    if(j==0)continue;
                    const double wj=midinterp_coefs[std::abs(j)-1];
                    kmix.blend_add(circ,ostates[i+j],wj);
                }
                kmix.blend_finalize(circ);
                kmix.rv(0,r,v);
                double kfac=1/(circ?-kmix.q:1+kmix.q);
                const auto perr=circ?&fit_err_t::kcerr:&fit_err_t::krerr;
                fiterrs[i].*perr=kfac*state_error_fun(&pdata[i].r,&r);
            }
            r=0;
            v=0;
            for(int_t j=wing;j>0;--j){
                const double wj=midinterp_coefs[j-1];
                r+=wj*(pdata[i+j].r+pdata[i-j].r);
                v+=wj*(pdata[i+j].v+pdata[i-j].v);
            }
            fiterrs[i].serr=state_error_fun(&pdata[i].r,&r);
        }

        for(int_t i=0;i<wing;++i){
            fiterrs[i]=fiterrs[wing];
            fiterrs[N-1-i]=fiterrs[N-1-wing];
        }
    }

    fit_err_t max_errs,mean_errs;
    memset(&max_errs,0,sizeof(fit_err_t));
    memset(&mean_errs,0,sizeof(fit_err_t));
    for(int_t i=0;i<N;++i){
        const auto &ei=fiterrs[i];
        checked_maximize(max_errs.serr,ei.serr);
        checked_maximize(max_errs.kcerr,ei.kcerr);
        checked_maximize(max_errs.krerr,ei.krerr);
        mean_errs.serr+=ei.serr*ei.serr;
        mean_errs.kcerr+=ei.kcerr*ei.kcerr;
        mean_errs.krerr+=ei.krerr*ei.krerr;
    }
    mean_errs.serr=std::sqrt(mean_errs.serr/N);
    mean_errs.kcerr=std::sqrt(mean_errs.kcerr/N);
    mean_errs.krerr=std::sqrt(mean_errs.krerr/N);
    checked_minimize(max_errs.serr,3*mean_errs.serr);
    checked_minimize(max_errs.kcerr,3*mean_errs.kcerr);
    checked_minimize(max_errs.krerr,3*mean_errs.krerr);

    //try fit state
    double min_kep_err=filtered_min(max_errs.kcerr,max_errs.krerr);
    double use_state=max_errs.serr<epsilon_relative_error?
        max_errs.serr:filtered_min(max_errs.serr,min_kep_err*(1<<d));
    double use_kepler=min_kep_err<epsilon_relative_error?
        min_kep_err:filtered_min(max_errs.serr,min_kep_err);
    std::vector<MFILE> compressed_results;
    if(max_errs.serr==use_state){
        //{x,y,z}=vec[1]
        std::vector<vec> fdata;
        fdata.reserve(N);
        for(int_t i=0;i<N;++i)fdata.push_back(pdata[i].r);
        compressed_results.push_back(
            (can_kepler?compress_data<KEPLERIAN_VECTORS>:compress_data<STATE_VECTORS>)
            (fdata.data(),N,d));
        size_t header_size=(can_kepler?sizeof(header_t<KEPLERIAN_VECTORS>):sizeof(header_t<STATE_VECTORS>));
        if(compressed_results.back().size()<=header_size)compressed_results.pop_back();
        else{
            vec left_r=pdata[0].r,left_v=pdata[0].v,right_r=pdata[N-1].r,right_v=pdata[N-1].v;
            if(can_kepler){
                auto *pheader=(header_t<KEPLERIAN_VECTORS>*)compressed_results.back().data();
                double v_factor=time_span/pheader->n;
                double compress_factor=double(N-1)/pheader->n;
                pheader->fix[0][0]= left_r-pheader->fix[0][0];
                pheader->fix[1][0]=( left_v*v_factor-pheader->fix[1][0]*compress_factor);
                pheader->fix[2][0]=right_r-pheader->fix[2][0];
                pheader->fix[3][0]=(right_v*v_factor-pheader->fix[3][0]*compress_factor);
            }
            else{
                auto *pheader=(header_t<STATE_VECTORS>*)compressed_results.back().data();
                double v_factor=time_span/pheader->n;
                double compress_factor=double(N-1)/pheader->n;
                pheader->fix[0][0]= left_r-pheader->fix[0][0];
                pheader->fix[1][0]=( left_v*v_factor-pheader->fix[1][0]*compress_factor);
                pheader->fix[2][0]=right_r-pheader->fix[2][0];
                pheader->fix[3][0]=(right_v*v_factor-pheader->fix[3][0]*compress_factor);
            }
        }
    }
    if(max_errs.kcerr==use_kepler){
        //kep=double[6]
        std::vector<double> fdata;
        fdata.reserve(6*N);
        double last_ml=0;
        for(const auto &oi:ostates){
            fdata.push_back(oi.j.x);
            fdata.push_back(oi.j.y);
            fdata.push_back(oi.j.z);
            fdata.push_back(oi.ex);
            fdata.push_back(oi.ey);
            double ml=oi.ml;
            ml=last_ml+angle_reduce(ml-last_ml);
            fdata.push_back(ml);
            last_ml=ml;
        }
        compressed_results.push_back(compress_data<KEPLERIAN_CIRCULAR>(fdata.data(),N,d));
        if(compressed_results.back().size()<=sizeof(header_t<KEPLERIAN_CIRCULAR>))compressed_results.pop_back();
        else{
            auto *pheader=(header_t<KEPLERIAN_CIRCULAR>*)compressed_results.back().data();
            const double *pleft=fdata.data();
            const double *pright=fdata.data()+6*(N-1);
            for(int_t ich=0;ich<6;++ich)
                pheader->fix[0][ich]= pleft[ich]-pheader->fix[0][ich];
            memset(pheader->fix[1],0,sizeof(pheader->fix[1]));
            for(int_t ich=0;ich<6;++ich)
                pheader->fix[2][ich]=pright[ich]-pheader->fix[2][ich];
            memset(pheader->fix[3],0,sizeof(pheader->fix[3]));
        }
    }
    if(max_errs.krerr==use_kepler){
        //kep=double[6]
        std::vector<double> fdata;
        fdata.reserve(6*N);
        for(const auto &oi:ostates){
            fdata.push_back(oi.j.x);
            fdata.push_back(oi.j.y);
            fdata.push_back(oi.j.z);
            fdata.push_back(oi.q);
            fdata.push_back(oi.el);
            fdata.push_back(oi.m);
        }
        compressed_results.push_back(compress_data<KEPLERIAN_RAW>(fdata.data(),N,d));
        if(compressed_results.back().size()<=sizeof(header_t<KEPLERIAN_RAW>))compressed_results.pop_back();
        else{
            auto *pheader=(header_t<KEPLERIAN_RAW>*)compressed_results.back().data();
            const double *pleft=fdata.data();
            const double *pright=fdata.data()+6*(N-1);
            for(int_t ich=0;ich<6;++ich)
                pheader->fix[0][ich]= pleft[ich]-pheader->fix[0][ich];
            memset(pheader->fix[1],0,sizeof(pheader->fix[1]));
            for(int_t ich=0;ich<6;++ich)
                pheader->fix[2][ich]=pright[ich]-pheader->fix[2][ich];
            memset(pheader->fix[3],0,sizeof(pheader->fix[3]));
        }
    }
    if(compressed_results.empty()){
        LogWarning("compress_orbital_data::no available method.\n");
        return 0;
    }

    return select_best(mf,compressed_results,N,d);
}
int_t ephemeris_compressor::compress_rotational_data(MFILE &mf,double time_span,MFILE *morb){
    mf.load_data();
    int_t N=mf.size();
    const rotational_state_t *pdata=(const rotational_state_t*)mf.prepare(N);
    constexpr int_t state_size=sizeof(rotational_state_t);
    if(N<2*state_size||N%state_size){
        LogError("compress_rotational_data::invalid file size.\n");
        return 0;
    }
    N/=state_size;

    int_t d=std::min(max_bspline_degree,N-1);
    d-=d+1&1;

    std::vector<MFILE> compressed_results;
    do{
        vec wlocal(0);
        double wmax=0;
        for(int_t i=0;i<N;++i){
            mat s=mat(pdata[i].x,0,pdata[i].z);
            vec wi=s.tolocal(pdata[i].w);
            wlocal+=wi;
            checked_maximize(wmax,wi.normsqr());
        }
        wmax=std::sqrt(wmax);
        double local_axis_theta=wlocal.theta(),local_axis_phi=wlocal.phi();
        mat local_axis(vec(local_axis_theta,local_axis_phi),vec(local_axis_theta+Constants::pi_div2,local_axis_phi),0);
        std::vector<double> fdata;
        fdata.reserve(6*N);
        double last_angle=0;
        for(int_t i=0;i<N;++i){
            mat s=mat(pdata[i].x,0,pdata[i].z);
            mat saxis(s.toworld(local_axis.x),s.toworld(local_axis.y),s.toworld(local_axis.z));
            ephem_rot ri(0,saxis,pdata[i].w);
            double py=std::sin(ri.ptheta);
            double pz=std::cos(ri.ptheta);
            py*=std::sin(ri.pphi);
            if(!(-HALF<=py&&py<=HALF)||!(-HALF<=pz&&pz<=HALF)||!(wmax*HALF<=saxis.x%pdata[i].w)){
                fdata.clear();
                break;
            }
            fdata.push_back(ri.w.x);
            fdata.push_back(ri.w.y);
            fdata.push_back(ri.w.z);
            fdata.push_back(py);
            fdata.push_back(pz);
            double angle=ri.angle;
            angle=last_angle+angle_reduce(angle-last_angle);
            fdata.push_back(angle);
            last_angle=angle;
        }

        if(fdata.empty())
            break;
        compressed_results.push_back(compress_data<AXIAL_OFFSET>(fdata.data(),N,d));
        if(compressed_results.back().size()<=sizeof(header_t<AXIAL_OFFSET>))compressed_results.pop_back();
        else{
            auto *pheader=(header_t<AXIAL_OFFSET>*)compressed_results.back().data();
            pheader->local_axis_theta=local_axis_theta;
            pheader->local_axis_phi=local_axis_phi;
            const double *pleft=fdata.data();
            const double *pright=fdata.data()+6*(N-1);
            for(int_t ich=0;ich<6;++ich)
                pheader->fix[0][ich]= pleft[ich]-pheader->fix[0][ich];
            memset(pheader->fix[1],0,sizeof(pheader->fix[1]));
            for(int_t ich=0;ich<6;++ich)
                pheader->fix[2][ich]=pright[ich]-pheader->fix[2][ich];
            memset(pheader->fix[3],0,sizeof(pheader->fix[3]));
        }
    } while(0);

    bool tidal_lockable=false;
    do{
        interpolator eorb(morb,time_span);
        if(!eorb.is_orbital())
            break;
        eorb.expand();
        vec left_w,right_w;
        std::vector<quat> fdata;
        fdata.reserve(N);
        for(int_t i=0;;++i){
            orbital_state_t os;
            eorb(double(i)/(N-1)*time_span,&os);
            vec &r=os.r,&j=os.v;
            double r2=r.normsqr();
            vec w=r*j/r2;
            r/=-std::sqrt(r2);
            mat jframe(r,0,w/w.norm());
            mat js(jframe.tolocal(pdata[i].x),0,jframe.tolocal(pdata[i].z));
            if(!(js.x.x>HALF)){
                fdata.clear();
                break;
            }
            fdata.emplace_back(js);
            if(i==0||i+1==N){
                (i==0?left_w:right_w)=jframe.tolocal(pdata[i].w-w);
                if(i)break;
            }
        }
        if(fdata.empty())
            break;
        for(int_t i=1;i<N;++i)if(fdata[i-1]%fdata[i]<0)fdata[i]=-fdata[i];
        compressed_results.push_back(compress_data<TIDAL_LOCK>(fdata.data(),N,d));
        if(compressed_results.back().size()<=sizeof(header_t<TIDAL_LOCK>))compressed_results.pop_back();
        else{
            tidal_lockable=true;
            quat left_q=fdata[0],right_q=fdata[N-1];
            quat left_dq=HALF*left_w*left_q;
            quat right_dq=HALF*right_w*right_q;
            auto *pheader=(header_t<TIDAL_LOCK>*)compressed_results.back().data();
            double v_factor=time_span/pheader->n;
            double compress_factor=double(N-1)/pheader->n;
            pheader->fix[0][0]= left_q-pheader->fix[0][0];
            pheader->fix[1][0]=( left_dq*v_factor-pheader->fix[1][0]*compress_factor);
            pheader->fix[2][0]=right_q-pheader->fix[2][0];
            pheader->fix[3][0]=(right_dq*v_factor-pheader->fix[3][0]*compress_factor);
        }
    } while(0);

    if(!tidal_lockable){
        std::vector<quat> fdata;
        fdata.reserve(N);
        for(int_t i=0;i<N;++i)fdata.emplace_back(mat(pdata[i].x,0,pdata[i].z));
        for(int_t i=1;i<N;++i)if(fdata[i-1]%fdata[i]<0)fdata[i]=-fdata[i];
        compressed_results.push_back(compress_data<QUATERNION>(fdata.data(),N,d));
        if(compressed_results.back().size()<=sizeof(header_t<QUATERNION>))compressed_results.pop_back();
        else{
            quat left_q=fdata[0],right_q=fdata[N-1];
            quat left_dq=HALF*pdata[0].w*left_q;
            quat right_dq=HALF*pdata[N-1].w*right_q;
            auto *pheader=(header_t<QUATERNION>*)compressed_results.back().data();
            double v_factor=time_span/pheader->n;
            double compress_factor=double(N-1)/pheader->n;
            pheader->fix[0][0]= left_q-pheader->fix[0][0];
            pheader->fix[1][0]=( left_dq*v_factor-pheader->fix[1][0]*compress_factor);
            pheader->fix[2][0]=right_q-pheader->fix[2][0];
            pheader->fix[3][0]=(right_dq*v_factor-pheader->fix[3][0]*compress_factor);
        }
    }
    if(compressed_results.empty()){
        LogWarning("compress_rotational_data::no available method.\n");
        return 0;
    }

    return select_best(mf,compressed_results,N,d);
}

int_t ephemeris_compressor::select_best(MFILE &mf,std::vector<MFILE> &compressed_results,int_t N,int_t d){
    MFILE *best_result=nullptr;
    double best_score=INFINITY;
    int_t n_optimal=0;
    for(auto &cmf:compressed_results){
        const auto *header=(const header_base*)cmf.data();
        double cscore=filtered_min<double>(header->relative_error,INFINITY)+epsilon_relative_error;
        cscore=compression_score(cscore,double(cmf.size()));
        if(cscore<best_score){
            best_score=cscore;
            best_result=&cmf;
            n_optimal=header->n;
        }
    }

    int_t ret=0;
    std::vector<int_t> n_all;
    segment_choices(n_all,N,d);
    int_t ns=n_all.size();
    for(int_t i=0;i<ns;++i){
        if(n_all[i]<=n_optimal){
            ret=i+1;
            break;
        }
    }
    if(!best_result||!ret)
        return 0;

    memcpy(mf.prepare(best_result->size()),best_result->data(),best_result->size());
    return ret;
}

void ephemeris_compressor::segment_choices(std::vector<int_t> &result,int_t N,int_t d){
    const int_t n_min=1;
    int_t n=std::max(n_min,(N-d)/(d+1));
    result.resize(1,n);
    while(n>n_min){
        n=std::max(n_min,n*3/4);
        result.push_back(n);
    }
}

typedef ephemeris_compressor::interpolator interp_t;

#define CONVERT_IMPLEMENT(F)  template<>        \
void interp_t::convert<ephemeris_format::F>(    \
    typename header_t<F>::state_type *pstate,   \
    const typename header_t<F>::data_type *x,   \
    const typename header_t<F>::data_type *v) const

CONVERT_IMPLEMENT(STATE_VECTORS){
    pstate->r=*x;
    pstate->v=*v;
}
CONVERT_IMPLEMENT(KEPLERIAN_VECTORS){
    pstate->r=*x;
    pstate->v=*v;
}
CONVERT_IMPLEMENT(KEPLERIAN_CIRCULAR){
    kep_t k(x,true);
    k.rv(0,pstate->r,pstate->v);
}
CONVERT_IMPLEMENT(KEPLERIAN_RAW){
    kep_t k(x,false);
    k.rv(0,pstate->r,pstate->v);
}
CONVERT_IMPLEMENT(AXIAL_OFFSET){
    axial a(x);
    mat saxis;
    a.sw(0,saxis,pstate->w);
    saxis%=local_axis;
    pstate->x=saxis.x;
    pstate->z=saxis.z;
}
CONVERT_IMPLEMENT(QUATERNION){
    mat s(*x);
    pstate->w=vec(2*(*v)/(*x));
    pstate->x=s.x;
    pstate->z=s.z;
}
CONVERT_IMPLEMENT(TIDAL_LOCK){
    mat s(*x);
    pstate->w=vec(2*(*v)/(*x));
    pstate->x=s.x;
    pstate->z=s.z;
    // correct by orbital state
    vec r=orb.r;
    double r2=r.normsqr();
    vec w=r*orb.v/r2;
    r/=-std::sqrt(r2);
    mat jframe(r,0,w/w.norm());
    pstate->x=jframe.toworld(pstate->x);
    pstate->z=jframe.toworld(pstate->z);
    pstate->w=jframe.toworld(pstate->w)+w;
}
#undef CONVERT_IMPLEMENT


#define CASE(F,...) case F:{            \
    constexpr format_t FORMAT=F;        \
    typedef header_t<F> header_type;    \
    typedef bspline_fitter<typename header_type::data_type,header_type::data_channel> bspline_t;  \
    bspline_t *&pfitter=(bspline_t*&)(this->pfitter);   \
    __VA_ARGS__                                         \
}break

#define SWITCH(F,...)   \
        switch(F){      \
            CASE(STATE_VECTORS      ,__VA_ARGS__);  \
            CASE(KEPLERIAN_VECTORS  ,__VA_ARGS__);  \
            CASE(KEPLERIAN_CIRCULAR ,__VA_ARGS__);  \
            CASE(KEPLERIAN_RAW      ,__VA_ARGS__);  \
            CASE(AXIAL_OFFSET       ,__VA_ARGS__);  \
            CASE(QUATERNION         ,__VA_ARGS__);  \
            CASE(TIDAL_LOCK         ,__VA_ARGS__);  \
        }


interp_t::interpolator(MFILE *fin,double _range):t_range(_range){
    pfitter=nullptr;
    do{
        if(!fin||!fin->publish()||fin->size()<=sizeof(base_t))
            break;
        const uint8_t *const fdata=fin->data();
        memcpy((base_t*)this,fdata,sizeof(base_t));
        double r=relative_error();
        if(!(0<=r&&r<1))
            break;
        ephemeris_format f=data_format();
        if(!is_valid_format(f))
            break;

        const size_t fsize=fin->size();

        SWITCH(f,
            if(fsize!=sizeof(header_type)+(n+degree)*bspline_t::data_size)break;
            pfitter=new bspline_t((typename bspline_t::data_type*)(fdata+sizeof(header_type)),
                degree,n,t_range);
            if(!*pfitter){
                delete pfitter;
                pfitter=nullptr;
            }
        );

        if(pfitter){
            //set auxiliary data
            if(f==AXIAL_OFFSET){
                auto *pheader=(header_t<AXIAL_OFFSET>*)fdata;
                double local_axis_theta=pheader->local_axis_theta;
                double local_axis_phi=pheader->local_axis_phi;
                local_axis=mat(
                    vec(local_axis_theta,local_axis_phi),
                    vec(local_axis_theta+Constants::pi_div2,local_axis_phi),
                    0).transpose();
            }
            return;
        }
    } while(0);
    //failed
    base_t::uformat=~uint16_t(ephemeris_format::INVALID);
    base_t::relative_error=NAN;
}

interp_t::interpolator(interpolator &&other){
    memcpy(this,&other,sizeof(interpolator));
    memset(&other,0,sizeof(interpolator));
}

interp_t::~interpolator(){
    if(!pfitter)
        return;
    ephemeris_format f=data_format();
    SWITCH(f,
        delete pfitter;
    );
}

void interp_t::expand(){
    if(!pfitter)
        return;
    ephemeris_format f=data_format();
    SWITCH(f,
        pfitter->expand();
    );
}

void interp_t::set_orbital_state(const vec &_r,const vec &_v){
    if(data_format()==TIDAL_LOCK){
        orb.r=_r;
        orb.v=_v;
    }
}

void interp_t::operator()(double t,void *state) const{
    if(!pfitter)
        return;

    double x=t/t_range*n;
    bool fix_left=std::floor(x)==0;
    bool fix_right=std::ceil(x)==n;

    ephemeris_format f=data_format();
    SWITCH(f,
        typename header_type::data_type r[header_type::data_channel],v[header_type::data_channel];
        (*pfitter)(t,r,v);
        if(fix_left){

        }
        if(fix_right){
            x=n-x;

        }
        convert<FORMAT>((typename header_type::state_type *)state,r,v);
    );
}

#undef SWITCH
#undef CASE
