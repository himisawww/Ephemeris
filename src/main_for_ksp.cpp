#include"physics/mass.h"
#include"modules/ephemeris_reader.h"
#include"modules/ephemeris_generator.h"
#include"math/state_parameters.h"
#include"utils/logger.h"
#include"physics/geopotential.h"

class principia_model{
    typedef const std::string &strref;
    static bool is_starts_with(strref s,const char *start){
        return !strncmp(s.c_str(),start,strlen(start));
    }
    static std::string trim(strref s,const char *space=" \r\n\t"){
        size_t b=s.find_first_not_of(space);
        size_t e=s.find_last_not_of(space)+1;
        return e>b?s.substr(b,e-b):"";
    }
    static bool expect(double &result,strref s,const char *unit){
        char *pend;
        result=strtod(s.c_str(),&pend);
        return s.c_str()!=pend&&trim(pend)==unit;
    }

public:
    struct body_model{
        std::string name;
        int_t mid,pmid;
        //m^3/s^2
        double gravitational_parameter;
        //deg & deg / d
        double ra,dec,W,angfreq;
        //m
        double ref_radius;
        //[n][m]=[J,S], normalized as principia's convention
        std::map<int,std::vector<std::pair<double,double>>> gpmodel;

        body_model(){
            mid=pmid=-1;
            gravitational_parameter=ra=dec=W=angfreq=ref_radius=NAN;
        }
        bool sanity() const{
            return !name.empty()&&mid>=0&&!(0*gravitational_parameter*ra*dec*W*angfreq)&&(gpmodel.empty()||!(0*ref_radius));
        }

        std::string to_string() const{
            std::string result=strprintf(
                "  body {\n"
                "    name                    = %s\n"
                "    gravitational_parameter = %.17e km^3/s^2\n"
                "    reference_instant       = JD2451545.000000000\n"
                "    axis_right_ascension    = %.6f deg\n"
                "    axis_declination        = %.6f deg\n"
                "    reference_angle         = %.6f deg\n"
                "    angular_frequency       = %.9f deg / d\n"
                ,name.c_str(),
                gravitational_parameter/1e9,ra,dec,W,angfreq);
            if(!gpmodel.empty()){
                result+=strprintf(
                    "    reference_radius        = %.9f km\n"
                    ,ref_radius/1000
                );
                auto use_jn=[](const std::vector<std::pair<double,double>> &gpn){
                    int n=gpn.size()-1;
                    bool result=false;
                    for(int m=0;m<=n;++m){
                        if(m&&gpn[m].first||gpn[m].second)
                            break;
                        result=m==n;
                    }
                    return result;
                };
                if(gpmodel.count(2)==gpmodel.size()&&use_jn(gpmodel.at(2)))
                    result+=strprintf(
                        "    j2                      = %.17e\n"
                        ,gpmodel.at(2)[0].first/geopotential_factor(2,0)
                    );
                else for(const auto &gp:gpmodel){
                    int n=gp.first;
                    const auto &gpn=gp.second;
                    result+=strprintf(
                        "    geopotential_row {\n"
                        "      degree = %d\n"
                        ,n);
                    if(use_jn(gpn))
                        result+=strprintf(
                            "      geopotential_column {\n"
                            "        order = 0\n"
                            "        j     = %.17e\n"
                            "        sin   = 0.00000000000000000e+00\n"
                            "      }\n"
                            ,gpn[0].first/geopotential_factor(n,0));
                    else for(int m=0;m<=n;++m){
                        if(!gpn[m].first&&!gpn[m].second)
                            continue;
                        result+=strprintf(
                            "      geopotential_column {\n"
                            "        order = %d\n"
                            "        cos   = %.17e\n"
                            "        sin   = %.17e\n"
                            "      }\n"
                            ,m,gpn[m].first,gpn[m].second);
                    }
                    result+="    }\n";
                }
            }
            result+="  }";
            return result;
        }
    };
    std::string header;
    std::vector<body_model> models;
    bool valid;

    explicit operator bool() const{ return valid; }
    principia_model(const ephemeris_reader &ereader,const char *gravity_model_path){
        valid=false;
        MFILE *fin=mopen(gravity_model_path,MFILE_STATE::READ_FILE);
        if(!fin)
            return;

        //state
        int level=0,m=-1,n=-1;

        std::string line;
        while(!(line=readline(fin)).empty()){
            line=trim(line);
            if(line.back()=='{'){
                const char *levels[]={"principia_gravity_model","body","geopotential_row","geopotential_column"};
                if(level>=4||!is_starts_with(line,levels[level]))
                    return;
                if(level==0)
                    header=line;
                else if(level==1)
                    models.emplace_back();
                ++level;
            }
            else if(line.back()=='}'){
                if(--level<0)
                    return;
                else if(level==1){
                    if(models.empty())
                        return;
                    body_model &cur_model=models.back();
                    if(!cur_model.sanity())
                        return;
                }
                else if(level==2)
                    n=-1;
                else if(level==3)
                    m=-1;
            }
            else if(level<2||4<level||models.empty())
                return;
            else{
                body_model &cur_model=models.back();
                size_t eqpos=line.find('=');
                if(eqpos==std::string::npos||line.find('=',eqpos+1)!=std::string::npos)
                    return;
                std::string eql=trim(line.substr(0,eqpos-1)),eqr=trim(line.substr(eqpos+1));
                if(level==2){
                    bool success=false;
                    if(eql=="name"){
                        cur_model.name=eqr;
                        const msystem &ms=ereader.get_msystem();
                        for(int_t mid=0,mn=ms.size();mid<mn;++mid){
                            if(ereader.get_massinfo(mid).name==cur_model.name){
                                cur_model.mid=mid;
                                success=true;
                                break;
                            }
                        }
                    }
                    else if(eql=="gravitational_parameter"){
                        if(success=expect(cur_model.gravitational_parameter,eqr,"km^3/s^2"))
                            cur_model.gravitational_parameter*=1e9;
                    }
                    else if(eql=="reference_instant")
                        success=eqr=="JD2451545.000000000";
                    else if(eql=="axis_right_ascension")
                        success=expect(cur_model.ra,eqr,"deg");
                    else if(eql=="axis_declination")
                        success=expect(cur_model.dec,eqr,"deg");
                    else if(eql=="reference_angle")
                        success=expect(cur_model.W,eqr,"deg");
                    else if(eql=="angular_frequency")
                        success=expect(cur_model.angfreq,eqr,"deg / d");
                    else if(eql=="reference_radius"){
                        if(success=expect(cur_model.ref_radius,eqr,"km"))
                            cur_model.ref_radius*=1000;
                        else
                            success=expect(cur_model.ref_radius,eqr,"m");
                    }
                    else if(eql=="j2"&&cur_model.gpmodel.empty()){
                        double &j2=cur_model.gpmodel.try_emplace(2,3).first->second[0].first;
                        if(success=expect(j2,eqr,""))
                            j2*=geopotential_factor(2,0);
                    }
                    if(!success)
                        return;
                }
                else if(level==3){
                    double dn;
                    if(n!=-1||eql!="degree"||!expect(dn,eqr,"")||(n=std::floor(dn))!=dn)
                        return;
                }
                else if(eql=="order"){
                    double dm;
                    if(m!=-1||!expect(dm,eqr,"")||(m=std::floor(dm))!=dm)
                        return;
                }
                else if(n<2||m<0||n<m)
                    return;
                else{
                    auto &gp=cur_model.gpmodel.try_emplace(n,n+1).first->second[m];
                    bool success=false;
                    if(eql=="j"){
                        if(success=m==0&&expect(gp.first,eqr,""))
                            gp.first*=geopotential_factor(n,0);
                    }
                    else if(eql=="cos")
                        success=expect(gp.first,eqr,"");
                    else if(eql=="sin")
                        success=expect(gp.second,eqr,"")&&(m!=0||!gp.second);
                    if(!success)
                        return;
                }
            }
        }

        if(level==0)
            valid=true;
    }
    bool append(int_t mid,const std::string &name,double GM,double R){
        if(std::find_if(models.begin(),models.end(),[mid](const body_model &bm){return mid==bm.mid;})!=models.end())
            return false;
        if(std::find_if(models.begin(),models.end(),[&name](const body_model &bm){return name==bm.name;})!=models.end()){
            LogError("Duplicate name %s in Principia config.\n",name.c_str());
            return false;
        }
        //LogInfo("Appending %s to Principia...\n",name.c_str());
        auto &mnew=models.emplace_back();
        mnew.mid=mid;
        mnew.name=name;
        mnew.gravitational_parameter=GM;
        mnew.ref_radius=R;
        return true;
    }

    std::string to_string() const{
        std::string result(header);
        result+='\n';
        for(auto &model:models){
            result+=model.to_string();
            result+='\n';
        }
        result+="}\n";
        std::string result_rn;
        result_rn.reserve(result.size());
        for(char ch:result){
            if(ch=='\n')
                result_rn+='\r';
            result_rn+=ch;
        }
        return result_rn;
    }

    //the normalization convention used by principia is differed from us...
    static double geopotential_factor(int n,int m){
        double factor=1/double(2*n+1);
        for(int k=n-m;k<n+m;)factor*=++k;
        return m==0?-std::sqrt(factor):std::sqrt(factor/2);
    }
};

//J2000=TDB 0
static double JDToTDB(double JD){
    return 86400*(JD-2451545);
}

int main_for_ksp(const char *eph_path,const char *principia_model_path,const char *export_path){
    ephemeris_reader ereader(eph_path);
    if(!ereader)
        return __LINE__;
    const msystem &ms=ereader.get_msystem();

    std::set<int_t> mids;
    //initial state for Principia
    mat ICRS_Equator(1);
    ICRS_Equator.rotx(-Constants::J2000_obliquity);
    constexpr double JD_Pricipia=2433282.5;
    if(ereader.checkout(JDToTDB(JD_Pricipia))){
        //Principia's sun, as origin
        vec rsun(+1.309126697236264e+05,+3.443856610385113e+05,+1.364602296561306e+05),
            vsun(-7.799754996220354e-03,-5.561927893069310e-03,-2.253148239533338e-03);
        rsun=ICRS_Equator.tolocal(ms["10"].r)/1000-rsun;
        vsun=ICRS_Equator.tolocal(ms["10"].v)/1000-vsun;
        MFILE *fout=mopen(strprintf("%sICRS_Equator_State[JD.%f].txt",export_path,JD_Pricipia),MFILE_STATE::WRITE_FILE);
        fprintf(fout,
            "Note: Instant state vectors at TDB JD %f, represented in ICRS-Equator.\n"
            "Note: Rotational states of highly librated objects (including most of moons, Luna, or even Earth) vary drastically "
            "even over a short period, hence the listed states shall not be used as fixed parameters for long term prediction.\n"
            "Note: Rotational states of irregular/chaotic satellites may be not reliable.\n\n",
            JD_Pricipia
        );
        for(const mass &m:ms){
            const char *ssid=m.get_ssid();
            if(*ssid=='A')
                break;
            vec r=ICRS_Equator.tolocal(m.r)/1000;
            vec v=ICRS_Equator.tolocal(m.v)/1000;
            vec w=ICRS_Equator.tolocal(m.w)*(86400/Constants::degree);
            vec x=ICRS_Equator.tolocal(m.s.x);
            vec z=ICRS_Equator.tolocal(m.s.z);
            r-=rsun;
            v-=vsun;
            fprintf(fout,"[%16s] %s:\n",ereader.get_massinfo(ssid).name.c_str(),ssid);
            fprintf(fout,"Position    { %+.17e, %+.17e, %+.17e }  km\n"  ,r.x,r.y,r.z);
            fprintf(fout,"Velocity    { %+.17e, %+.17e, %+.17e }  km/s\n",v.x,v.y,v.z);
            fprintf(fout,"Rotation    { %+.17e, %+.17e, %+.17e } deg/d\n",w.x,w.y,w.z);
            fprintf(fout,"       X    { %+.17e, %+.17e, %+.17e } unit\n" ,x.x,x.y,x.z);
            fprintf(fout,"       Z    { %+.17e, %+.17e, %+.17e } unit\n" ,z.x,z.y,z.z);
            fprintf(fout,"\n");
            mids.insert(ms.get_mid(ssid));
        }
        fclose(fout);
    }

    principia_model pmodel(ereader,principia_model_path);
    if(!pmodel){
        LogError("Cannot load %s as Principia gravity_model.\n",principia_model_path);
        return __LINE__;
    }

    for(int_t mid:mids){
        const mass &mi=ms[mid];
        const auto &minfo=ereader.get_massinfo(mid);
        pmodel.append(mid,minfo.name,mi.GM0,mi.R);
    }

    //record data over 1951~2049
    constexpr double JD_Elements=2433647.5;
    constexpr double interval=3600;
    const auto &bs=ms.get_barycens();
    int_t n_points=0,i_j2000=-1;
    struct statistics{
        struct stat_param{
            //rotation
            vec w,x,z;
            //if has orbit:
            vec r,v;
            double GM;
            //harmonics 2
            mat h2;
        };
        int_t pmid;
        mat eq_frame;
        rotational_param_t rotparam;
        double GM_nominal;
        bool rotation_reverse;
        std::vector<stat_param> data;

        statistics(int_t _pmid):pmid(_pmid){}
    };
    std::map<int_t,statistics> mstats;
    ereader.update_physics_parallel_option=1;
    ereader.update_bsystem=true;
    ereader.update_physics=true;
    ereader.update_orbits=true;
    const double t_start=JDToTDB(JD_Elements);
    for(double t_eph=t_start;t_eph<=-t_start;t_eph+=interval){
        if(!ereader.checkout(t_eph)){
            LogError("Checkout Failed at %llds\n",(int_t)t_eph);
            return __LINE__;
        }
        if(n_points%256==0)LogInfo(" %lld\r",n_points);
        if(t_eph==0)i_j2000=n_points;

        for(auto &pm:pmodel.models){
            const mass &mi=ms[pm.mid];
            const auto &minfo=ereader.get_massinfo(pm.mid);
            const auto &bi=bs[pm.mid];
            if(bi.mid!=pm.mid||bi.hid>=0){
                LogError("Weired bsystem.\n");
                return __LINE__;
            }
            const int_t pbid=bs[bi.tid].pid,pmid=pbid<0?-1:bs[pbid].mid;
            if(n_points==0)
                pm.pmid=pmid;
            else if(pm.pmid!=pmid){
                LogError("Parent changed, not supported.\n");
                return __LINE__;
            }
            auto &mrec=mstats.try_emplace(pm.mid,pmid).first->second.data.emplace_back();
            //rotations
            mrec.w=ICRS_Equator.tolocal(mi.w);
            mrec.x=ICRS_Equator.tolocal(mi.s.x);
            mrec.z=ICRS_Equator.tolocal(mi.s.z);
            //orbits
            if(pmid>=0){
                mrec.GM=minfo.keplerian_GM;
                mrec.r=ICRS_Equator.tolocal(minfo.state_vectors.r);
                mrec.v=ICRS_Equator.tolocal(minfo.state_vectors.v);
            }
            //harmonics
            mrec.h2=fast_mpmat(mi.s).tolocal(mi.C_potential);
        }
        ++n_points;
    }

    //process statistics
    const double total_time=(n_points-1)*interval;
    auto make_ra=[](double alpha){
        alpha=angle_reduce(alpha);
        return (alpha<0?alpha+Constants::pi_mul2:alpha)/Constants::degree;
    };
    auto make_dec=[](double delta){
        return delta/Constants::degree;
    };
    //build rotation & physical properties
    for(auto &pm:pmodel.models){
        int_t mid=pm.mid;
        auto &mstat=mstats.at(mid);
        int_t pmid=mstat.pmid;
        const auto &mdata=mstat.data;
        if(mdata.size()!=n_points){
            LogError("Missing data.\n");
            return __LINE__;
        }

        vec jrot(0);
        double ra2000(NAN),dec2000(NAN),W2000(NAN);
        double rot_sum=0,last_rot;
        bool rotation_failed=false/*,rotation_reverse=pm.angfreq<0*/;
        mat h2sum(0);
        for(int_t i=0;i<n_points;++i){
            const auto &mdi=mdata[i];

            h2sum+=mdi.h2;
            while(!rotation_failed){
                jrot+=mdi.w;
                vec mz=/*rotation_reverse?-mdi.z:*/mdi.z;
                if(i==i_j2000&&pm.name=="Earth"){
                    //fix earth's ra/dec to 0/90 at J2000, this may lead to some error... on the order of nutation & leap seconds.
                    mz=vec(0,0,1);
                }
                double alpha=mz.lon();
                double delta=mz.lat();
                mat asc_node(vec::from_lon_lat(alpha+Constants::pi_div2,0),0,mz);
                double W=asc_node.lon(mdi.x);
                double rot_angle=alpha+Constants::pi_div2+W;
                if(i){
                    double rot_expect=last_rot+mdi.w.norm()*(/*rotation_reverse?-interval:*/interval);
                    double delta_rot=angle_reduce(rot_angle-rot_expect);
                    if(std::abs(delta_rot)>Constants::pi_div4){
                        LogError("Rotation Error for %s.\n",pm.name.c_str());
                        rotation_failed=true;
                        break;
                    }
                    rot_angle=rot_expect+delta_rot;
                    rot_sum+=rot_angle-last_rot;
                }
                last_rot=rot_angle;
                if(i==i_j2000){
                    ra2000=alpha;
                    dec2000=delta;
                    W2000=W;
                }
                break;
            }
        }
        {
            pm.gravitational_parameter=ms[mid].GM0;

            double c2[3],s2[2];
            h2sum/=n_points;
            h2sum.to_harmonics(c2[0],c2[1],c2[2],s2[0],s2[1]);

            if(pm.gpmodel.empty()){
                bool add_s21=std::abs(s2[0])>0.005;
                bool add_s22=add_s21||std::abs(s2[1])>0.005;
                bool add_c21=add_s21||std::abs(c2[1])>0.005;
                bool add_c22=add_s22||std::abs(c2[2])>0.005;
                bool add_j2=add_c22||std::abs(c2[0])>0.01;
                if(add_j2||add_c21||add_c22||add_s21||add_s22){
                    auto &hp=pm.gpmodel.try_emplace(2,3).first->second;
                    if(add_j2)hp[0].first=NAN;
                    if(add_c21)hp[1].first=NAN;
                    if(add_s21)hp[1].second=NAN;
                    if(add_c22)hp[2].first=NAN;
                    if(add_s22)hp[2].second=NAN;
                    pm.ref_radius=ms[mid].R;
                }
            }

            double Rfac=ms[mid].R/pm.ref_radius;
            const geopotential *mgp=ms[mid].gpmodel;
            for(auto &hp:pm.gpmodel){
                int n=hp.first;
                for(int m=0;m<=n;++m){
                    auto &cs=hp.second[m];
                    bool has_cos=cs.first,has_sin=cs.second;
                    if(!has_cos&&!has_sin)
                        continue;
                    double pfac=principia_model::geopotential_factor(n,m)*std::pow(Rfac,n);
                    if(cs.first){
                        if(n==2)cs.first=pfac*c2[m];
                        else if(mgp){
                            double cnm=mgp->get_C(n,m);
                            if(cnm)cs.first=cnm*pfac;
                        }
                    }
                    if(cs.second){
                        if(n==2)cs.second=pfac*s2[m-1];
                        else if(mgp){
                            double snm=mgp->get_S(n,m);
                            if(snm)cs.second=snm*pfac;
                        }
                    }
                }
            }
        }
        if(rotation_failed)
            mstat.rotparam.w=NAN;
        else{
            rot_sum/=total_time;
            mat sj2000(1);
            sj2000.rotz(Constants::pi_div2+ra2000).rotx(Constants::pi_div2-dec2000).rotz(W2000);
            mstat.rotparam=rotational_param_t(sj2000,sj2000.z*rot_sum);
            mstat.rotation_reverse=!(pm.angfreq>0)&&ICRS_Equator.toworld(jrot).z<0;
            if(mstat.rotation_reverse){
                jrot=-jrot;
                rot_sum=-rot_sum;
                W2000=Constants::pi-W2000;
                ra2000+=Constants::pi;
                dec2000=-dec2000;
            }
            mstat.eq_frame=mat(1).rotz(Constants::pi_div2+ra2000).rotx(Constants::pi_div2-dec2000);
            mstat.GM_nominal=pm.gravitational_parameter;
            double mean_alpha=make_ra(jrot.lon()),mean_delta=make_dec(jrot.lat());
            ra2000=make_ra(ra2000);
            dec2000=make_dec(dec2000);
            W2000=make_ra(W2000);
            if(W2000>359.9999||W2000<0.0001)W2000=0;

            pm.ra=ra2000;
            pm.dec=dec2000;
            pm.W=W2000;
            pm.angfreq=rot_sum*(86400/Constants::degree);
        }
    }

    //build orbit
    double err_k=0;
    int err_flag=0;
    MFILE *forbit=mopen(strprintf("%sorbit_parameter[JD.%f].txt",export_path,JD_Elements),MFILE_STATE::WRITE_FILE);
    for(auto &pm:pmodel.models){
        int_t mid=pm.mid;
        auto &mstat=mstats.at(mid);
        int_t pmid=mstat.pmid;
        const auto &mdata=mstat.data;
        if(pmid<0)
            continue;
        const auto &pstat=mstats.at(pmid);
        const auto &pdata=pstat.data;
        const mat &ps=pstat.eq_frame;
        if(pdata.size()!=n_points||!pstat.rotparam.w.is_finite()){
            LogError("Missing orbit parent data.\n");
            err_flag=__LINE__;
            break;
        }

        double GM_sum=0;
        double mean_motion_sum=0,last_mean_longitude;
        double inclsum=0,esum=0;
        double lan0(NAN),argp0(NAN),ml2000(NAN);
        bool orbit_failed=false;
        for(int_t i=0;i<n_points;++i){
            const auto &mdi=mdata[i];

            while(!orbit_failed){
                const auto &pdi=pdata[i];
                mat psi(pdi.x,0,pstat.rotation_reverse?-pdi.z:pdi.z);
                GM_sum+=mdi.GM;
                keplerian::mu mu(mdi.GM);
                keplerian keq(mu,mdi.r,mdi.v),kprot(mu,psi.tolocal(mdi.r),psi.tolocal(mdi.v));
                orbital_param_t orbeq(keq);
                inclsum+=kprot.j.theta();
                esum+=keq.q+1;
                if(i){
                    double ml_expect=last_mean_longitude+keq.mean_motion(mu)*interval;
                    double delta_ml=angle_reduce(orbeq.ml-ml_expect);
                    if(std::abs(delta_ml)>Constants::pi_div4){
                        LogError("Orbit Error for %s.\n",pm.name.c_str());
                        orbit_failed=true;
                        break;
                    }
                    mean_motion_sum+=ml_expect+delta_ml-last_mean_longitude;
                }
                if(i==i_j2000)
                    ml2000=orbeq.ml;
                if(i==0){
                    keplerian k0(mu,ps.tolocal(mdi.r),ps.tolocal(mdi.v));
                    lan0=k0.j.asc_node().lon();
                    argp0=k0.earg;
                }
                last_mean_longitude=orbeq.ml;
                break;
            }
        }
        if(!orbit_failed){
            GM_sum/=n_points;
            inclsum/=n_points;
            esum/=n_points;
            mean_motion_sum/=total_time;

            //rebuild orbit using
            //longitude of ascending node & argument of periapsis, in parent-body's equator frame, at Epoch 0;
            //inclination about parent's equator & eccentricity & mean sidereal motion, averaged over [1951,2049];
            //mean anomaly at Epoch 0, obtained by instant mean anomaly at J2000 and mean motion.
            keplerian kprot;
            keplerian::mu mu(GM_sum);
            kprot.j=vec(0,0,1).rotx(inclsum).rotz(lan0);
            kprot.q=esum-1;
            kprot.earg=argp0;
            kprot.m=0;
            kprot.j*=std::cbrt(kprot.mean_motion(mu)/mean_motion_sum);

            //convert to ICRS Equator
            vec r,v,ri,vi;
            keplerian keq=keplerian::toworld(ps,kprot);
            double mfac=-keq.q*(keq.q+2);
            mfac*=std::sqrt(mfac);
            keq.m+=(angle_reduce(ml2000+mean_motion_sum*t_start)-orbital_param_t(keq).ml)/mfac;
            kprot=keplerian::tolocal(ps,keq);

            //orbital error analysis
            double err_m=0,err_p=0;
            for(int_t i=0;i<n_points;++i){
                const auto &mdi=mdata[i];
                kprot.rv(mu,i*interval,r,v);
                keq.rv(mu,i*interval,ri,vi);
                checked_maximize(err_p,std::atan2((mdi.r*ri).norm(),mdi.r%ri));
                ri=ps.tolocal(ri);
                checked_maximize(err_k,std::atan2((ri*r).norm(),ri%r)/(1+mean_motion_sum*i*interval/Constants::pi_mul2));

                double cur_ml=orbital_param_t(keplerian(mu,mdi.r,mdi.v)).ml;
                double ml_expect=ml2000+(t_start+i*interval)*mean_motion_sum;
                checked_maximize(err_m,std::abs(angle_reduce(cur_ml-ml_expect)));
            }

            //rotational error analysis
            double err_x=0,err_z=0;
            if(mstat.rotparam.w.is_finite())for(int_t i=0;i<n_points;++i){
                const auto &mdi=mdata[i];
                mat st;
                vec wt;
                mstat.rotparam.sw(t_start+i*interval,st,wt);
                checked_maximize(err_x,(st.x-mdi.x).norm());
                checked_maximize(err_z,(st.z-mdi.z).norm());
            }
            else err_x=err_z=NAN;
            fprintf(forbit,
                "%s[%s]:\n"
                "   # maxerr of orbital mean longitude, position: [%f, %f] rad\n"
                "   # maxerr of rotation meridian, pole: [%f, %f] chord\n"
                "       meanMotion = %.17e\n"
                "       semiMajorAxis = %.17e\n"
                "       eccentricity = %.17f\n"
                "       meanAnomalyAtEpochD = %.17f\n"
                "   # in body mean equator frame of parent\n"
                "       inclination = %.17f\n"
                "       longitudeOfAscendingNode = %.17f\n"
                "       argumentOfPeriapsis = %.17f\n"
                "   # in ICRS Equator\n"
                "       inclination = %.17f\n"
                "       longitudeOfAscendingNode = %.17f\n"
                "       argumentOfPeriapsis = %.17f\n"
                "       epoch = 0\n"
                ,pm.name.c_str(),ms[mid].get_ssid(),
                err_m,err_p,err_x,err_z,
                mean_motion_sum,
                std::cbrt(pstat.GM_nominal/(mean_motion_sum*mean_motion_sum)),
                esum,
                make_ra(kprot.m*mfac),
                make_dec(kprot.j.theta()),
                make_ra(kprot.j.asc_node().lon()),
                make_ra(kprot.earg),
                make_dec(keq.j.theta()),
                make_ra(keq.j.asc_node().lon()),
                make_ra(keq.earg));

        }
    }
    fclose(forbit);
    if(err_flag)
        return err_flag;
    printf("max keplerian error %.17e\n",err_k);

    MFILE *fpcfg=mopen(strprintf("%snew_gravity_model.cfg",export_path),MFILE_STATE::WRITE_FILE);
    std::string pcfg=pmodel.to_string();
    fwrite(pcfg.c_str(),1,pcfg.size(),fpcfg);
    fclose(fpcfg);

    return 0;
}
