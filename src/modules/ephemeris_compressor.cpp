#include"ephemeris_compressor.h"
#include<map>
#include"math/interp.h"
#include"utils/logger.h"

#define HALF 0.5
typedef ephemeris_compressor::keplerian kep_t;

kep_t::keplerian(const ephem_orb &other){
    static_cast<ephem_orb&>(*this)=other;
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

template<ephemeris_compressor::Format format,typename T,size_t N_Channel>
MFILE ephemeris_compressor::compress_data(const T *pdata,int_t N,int_t d){
    typedef data_header<format> header_t;

    double (*error_fun)(const T*ref,const T*val);
    if constexpr(format==Format::STATE_VECTORS)
        error_fun=absolute_state_error;
    else if constexpr(format==Format::KEPLERIAN_VECTORS)
        error_fun=relative_state_error;
    else if constexpr(format==Format::KEPLERIAN_CIRCULAR)
        error_fun=circular_kepler_error;
    else if constexpr(format==Format::KEPLERIAN_RAW)
        error_fun=raw_kepler_error;
    else if constexpr(format==Format::AXIAL_OFFSET)
        error_fun=axial_rotation_error;
    else if constexpr(format==Format::QUATERNION)
        error_fun=quaterion_rotation_error;
    else return MFILE();

    struct scored_compress_data{
        MFILE data;
        header_t *pheader;
        double max_error;
        double reduced_error;
        double score;
    };
    std::map<int_t,scored_compress_data> mfmap;
    const int_t n_min=1;
    const int_t n_max=std::max(n_min,(N-d)/(d+1));
    
    auto make_fit=[&](int_t n){
        auto &target=mfmap[n];
        MFILE &mf=target.data;
        const size_t bsp_size=N_Channel*sizeof(T)*(n+d);
        const size_t data_size=bsp_size+sizeof(header_t);
        bspline_fitter<T,N_Channel> bf(d,n,pdata,N-1);
        target.pheader=(header_t*)mf.prepare(data_size);
        memcpy(target.pheader+1,bf.get_fitted_data(),bsp_size);
        bf.expand();
        T result[N_Channel];
        double max_error=0;
        for(int_t i=0;i<N;++i){
            bf(double(i),result);
            double ierror=error_fun(pdata+N_Channel*i,result);
            checked_maximize(max_error,ierror);
        }
        target.max_error=max_error;
        double score=filtered_min<double>(max_error,INFINITY)+epsilon_relative_error;
        target.reduced_error=score;
        target.score=compression_score(score,double(data_size));
        return &target;
    };

    std::vector<int_t> n_all(1,-1),i_chs;
    int_t n=n_max;
    while(n>n_min){
        n=std::max(n_min,n*3/4);
        n_all.push_back(n);
    }

    make_fit(n_max);
    int_t i_current=0,i_left=0,i_right=n_all.size()-1;
    int_t step=1,dir=0,n_optimal=n_max;

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
    header_t &header=*result.pheader;
    header.relative_error=(float)result.max_error;
    header.degree=(int16_t)d;
    header.uformat=~(uint16_t)format;
    header.n=n_optimal;
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
    return 0;
}
double ephemeris_compressor::quaterion_rotation_error(const quat *q,const quat *qp){
    return 0;
}

bool ephemeris_compressor::compress_orbital_data(MFILE &mf,double delta_t){
    mf.load_data();
    int_t N=mf.size();
    const orbital_state_t *pdata=(const orbital_state_t*)mf.prepare(N);
    constexpr int_t state_size=sizeof(orbital_state_t);
    if(N<2*state_size||N%state_size){
        LogError("compress_orbital_data::invalid file size.\n");
        return false;
    }
    N/=state_size;

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
            return false;
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
        size_t header_size=(can_kepler?sizeof(data_header<KEPLERIAN_VECTORS>):sizeof(data_header<STATE_VECTORS>));
        if(compressed_results.back().size()<=header_size)compressed_results.pop_back();
    }
    if(max_errs.kcerr==use_kepler){
        //kep=double[6]
        std::vector<double> fdata;
        fdata.reserve(6*N);
        double last_ml=ostates.front().ml;
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
        compressed_results.push_back(compress_data<Format::KEPLERIAN_CIRCULAR>(fdata.data(),N,d));
        if(compressed_results.back().size()<=sizeof(data_header<KEPLERIAN_CIRCULAR>))compressed_results.pop_back();
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
        compressed_results.push_back(compress_data<Format::KEPLERIAN_RAW>(fdata.data(),N,d));
        if(compressed_results.back().size()<=sizeof(data_header<KEPLERIAN_RAW>))compressed_results.pop_back();
    }
    if(compressed_results.empty()){
        LogWarning("compress_orbital_data::no available method.\n");
        return false;
    }

    MFILE *best_result=nullptr;
    double best_score=INFINITY;
    for(auto &cmf:compressed_results){
        const auto *header=(const data_header_base*)cmf.data();
        double cscore=filtered_min<double>(header->relative_error,INFINITY)+epsilon_relative_error;
        cscore=compression_score(cscore,double(cmf.size()));
        if(cscore<best_score){
            best_score=cscore;
            best_result=&cmf;
        }
    }
    memcpy(mf.prepare(best_result->size()),best_result->data(),best_result->size());
    return true;
}
bool ephemeris_compressor::compress_rotational_data(MFILE &mf,double dt){

    return false;
}
