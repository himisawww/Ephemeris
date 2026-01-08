#include"physics/mass.h"
#include"modules/ephemeris_reader.h"
#include"modules/ephemeris_generator.h"
#include"math/state_parameters.h"
#include"utils/logger.h"

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
        int_t mid;
        double gravitational_parameter;
        double ra,dec,W,angfreq;

        double ref_radius;
        //[n][m]=[J,S], normalized as principia
        std::map<int,std::vector<std::pair<double,double>>> gpmodel;

        body_model(){
            mid=-1;
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
                        ,gpmodel.at(2)[0].first
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
                            ,gpn[0].first);
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
                        success=expect(cur_model.gpmodel.try_emplace(2,3).first->second[0].first,eqr,"");
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
                    if(eql=="j")
                        success=m==0&&expect(gp.first,eqr,"");
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
};

//J2000=TDB 0
static double JDToTDB(double JD){
    return 86400*(JD-2451545);
}

int main_for_ksp(const char *eph_path){
    ephemeris_reader ereader(eph_path);
    if(!ereader)
        return -1;
    const msystem &ms=ereader.get_msystem();

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
        MFILE *fout=mopen(strprintf("r:\\ICRS_Equator_State[JD.%f].txt",JD_Pricipia),MFILE_STATE::WRITE_FILE);
        fprintf(fout,
            "Note: Instant state vectors at TDB JD %f, represented in ICRS-Equator.\n"
            "Note: Rotational states of highly librated objects (including most of moons, Luna, or even Earth) "
            "vary drastically even over a short period, hence the listed states shall not be used as fixed parameters for long term prediction.\n"
            "Note: Rotational states of irregular/chaotic satellites may be not reliable.\n\n",
            JD_Pricipia
        );
        for(const mass &m:ms){
            const char *ssid=(const char *)&m.sid;
            if(*ssid=='A')
                break;
            vec r=ICRS_Equator.tolocal(m.r)/1000;
            vec v=ICRS_Equator.tolocal(m.v)/1000;
            vec w=ICRS_Equator.tolocal(m.w);
            vec x=ICRS_Equator.tolocal(m.s.x);
            vec z=ICRS_Equator.tolocal(m.s.z);
            r-=rsun;
            v-=vsun;
            fprintf(fout,"[%16s] %s:\n",ereader.get_massinfo(ssid).name.c_str(),ssid);
            fprintf(fout,"Position    { %+.17e, %+.17e, %+.17e }  km\n"  ,r.x,r.y,r.z);
            fprintf(fout,"Velocity    { %+.17e, %+.17e, %+.17e }  km/s\n",v.x,v.y,v.z);
            fprintf(fout,"Rotation    { %+.17e, %+.17e, %+.17e } rad/s\n",w.x,w.y,w.z);
            fprintf(fout,"       X    { %+.17e, %+.17e, %+.17e } unit\n" ,x.x,x.y,x.z);
            fprintf(fout,"       Z    { %+.17e, %+.17e, %+.17e } unit\n" ,z.x,z.y,z.z);
            fprintf(fout,"\n");
        }
        fclose(fout);
    }

    const char *principia_model_path="E:\\SteamLibrary\\steamapps\\KSP_RSS\\GameData\\Principia\\real_solar_system\\gravity_model.cfg";
    principia_model pmodel(ereader,principia_model_path);
    if(!pmodel){
        LogError("Cannot load %s as Principia gravity_model.\n",principia_model_path);
        return 1;
    }

    MFILE *fpcfg=mopen("R:\\gravity_model.cfg",MFILE_STATE::WRITE_FILE);
    std::string pcfg=pmodel.to_string();
    fwrite(pcfg.c_str(),1,pcfg.size(),fpcfg);
    fclose(fpcfg);

    int_t n_points=0;
    constexpr double JD_Elements=2433647.5;
    for(double t_start=JDToTDB(JD_Elements),t_eph=t_start;t_eph<=-t_start;t_eph+=3600){
        if(!ereader.checkout(t_eph)){
            LogError("Checkout Failed at %llds\n",(int_t)t_eph);
            break;
        }
        ++n_points;
        if(n_points%256==0)printf(" %lld\r",n_points);
    }

    return 0;
}