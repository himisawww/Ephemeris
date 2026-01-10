#include"physics/mass.h"
#include<deque>
#include"physics/geopotential.h"
#include"physics/ring.h"
#include"utils/zipio.h"
#include"configs.h"
#include"utils/logger.h"

using Configs::MAX_LINESIZE;
using Configs::MAX_PATHSIZE;

const char *default_path=".";

static bool sanity(double dt){
    return dt>0&&dt!=INFINITY;
}

bool mass::sanity(bool alert) const{
    std::vector<const char*> errs;
    constexpr fast_real eps2=Constants::epsilon*Constants::epsilon;
    constexpr fast_real large_finite=1/eps2;

    if(!(fast_mpvec(r).normsqr()<large_finite*large_finite)||!(fast_mpvec(v).normsqr()<Constants::c2))
        errs.push_back("position/velocity not finite\n");

    fast_mpmat fmis(s);
    if(!((1-fmis.transpose()%fmis).normsqr()<64*eps2)||!(fmis.det()>0))
        errs.push_back("rotation not normalized\n");

    if(!(w.normsqr()*R2<Constants::c2))
        errs.push_back("rotation too fast\n");

    if(!(eps2<GM0&&GM0<large_finite))
        errs.push_back("invalid gravititional parameter\n");

    if(!(2*GM0<Constants::c2*std::abs(R)))
        errs.push_back("invalid radius\n");

    if(!(eps2<inertia&&inertia<=1))
        errs.push_back("invalid inertia factor\n");

    if(!(0<=k2&&k2<=1.5&&0<=k2r&&k2r<=1.5))
        errs.push_back("invalid tidal/rotational deformation factor\n");

    if(!(0<=tide_delay&&tide_delay<large_finite&&0<=tide_delay_factor&&tide_delay_factor<=large_finite))
        errs.push_back("invalid tidal delay time/factor\n");

    fast_real j2,c21,c22,s21,s22;
    C_static.to_harmonics(j2,c21,c22,s21,s22);
    if(!(0<=j2&&j2+6*std::abs(c22)<3*inertia&&std::abs(c21)<inertia&&std::abs(s21)<inertia&&std::abs(s22)*2<inertia))
        errs.push_back("invalid static harmonics\n");

    if(R<0&&(k2||k2r||j2||dJ2||c21||c22||s21||s22||gpmodel||ringmodel||tide_delay))
        errs.push_back("soft mass shall not have non-point components(harmonics/deformations/geopotential/ring)\n");

    if(errs.empty())
        return true;
    
    if(alert){
        LogError("Error when loading <%s>:\n",(char*)&sid);
        for(const char *estr:errs)
            LogError("    %s",estr);
    }
    return false;
}

bool msystem::load(const char *fconfig,const char *fcheckpoint){
    bool success=false;
  do{
    const std::string sip=fconfig;
    size_t spos=1+sip.find_last_of("/\\");
    std::map<std::string,std::string> config;
    std::string dirstr(spos?sip.substr(0,spos):".\\");
    std::string fconfigstr(sip.substr(spos));

    const char *dir=dirstr.c_str();
    fconfig=fconfigstr.c_str();

    if(!dir||!*dir)dir=default_path;
    if(!fconfig||!*fconfig)break;
    if(strlen(dir)+strlen(fconfig)+1>=MAX_PATHSIZE){
        LogError("Path too long: %s\\%s\n",dir,fconfig);
        break;
    }

    const char *integrators_list[]={
        "CPU_RK12","GPU_RK12","COMBINED_RK12"
    };
    const char *config_keys_list[]={
        "Initial",
        "Extra",
        "Geopotentials",
        "Rings",
        "TDB",
        "Integrator",
        "Delta_t",
        "Cadence",
        "Max_Ephemeris_Length",
        "Combined_Delta_t_Max",
        "Combined_GM_Max_Parent",
        "Combined_GM_Max_Child",
        "Combined_GM_Max_Tiny",
        "Combined_Period_Max_Child"
    };
    const int config_keys_size=sizeof(config_keys_list)/sizeof(char *);
    config.clear();
    for(int i=0;i<config_keys_size;++i)
        config.insert({config_keys_list[i],""});

    char fname[MAX_LINESIZE],sname[MAX_LINESIZE],sval[MAX_LINESIZE];
    sprintf(fname,"%s\\%s",dir,fconfig);

    bool failed=false;
    MFILE *fin=mopen(fname);
    if(!fin){
        LogError("File Not Exist: %s\n",fname);
        break;
    }

    while(1){
        std::string chbuf=readline(fin);
        if(chbuf.size()==0)break;
        if(chbuf.size()>=MAX_LINESIZE){
            LogError(
                "Loading %s\n Line Too Long: %s\n",
                fname,chbuf.c_str());
            failed=true;
            break;
        }

        int n=sscanf(chbuf.c_str(),"%s%s",sname,sval);
        if(n!=2){
            LogError(
                "Loading %s\n Line Length Error: %s\n",
                fname,chbuf.c_str());
            failed=true;
            break;
        }

        auto it=config.find(sname);
        if(it==config.end()){
            LogError(
                "Loading %s\n Unknown Parameter: %s\n",
                fname,sname);
            failed=true;
            break;
        }

        if(it->second.size()){
            LogError(
                "Loading %s\n Duplicate Parameter: %s\n",
                fname,sname);
            failed=true;
            break;
        }
        
        it->second=sval;
    }
    fclose(fin);
    if(failed)break;

    std::string
        &fbase      =config["Initial"],
        &fext       =config["Extra"],
        &gppath     =config["Geopotentials"],
        &ringpath   =config["Rings"];

    if(!fbase.size()){
        LogError(
            "Loading %s\n Missing Parameter: Initial\n",
            fname);
        break;
    }
    
    sprintf(sname,"%s\\",dir);

    if(!load((sname+fbase).c_str(),
        fext.size()?(sname+fext).c_str():nullptr,
        gppath.size()?(sname+gppath).c_str():nullptr,
        ringpath.size()?(sname+ringpath).c_str():nullptr
        ))break;

    const char *pval;

    pval=config["Integrator"].c_str();
    if(*pval){
        integrator=-1;
        const int integrators_size=sizeof(integrators_list)/sizeof(char*);
        for(int i=0;i<integrators_size;++i){
            if(0==strcmp(pval,integrators_list[i]))integrator=i;
        }
        if(integrator<0){
            LogError(
                "Loading %s\n Invalid Integrator: %s\n",
                fname,pval);
            break;
        }
    }
    else LogWarning(
        "Warning: Using Default Integrator: %s\n",
        integrators_list[integrator]);

#define LOAD_CONFIG(NAME,PARAM) do{                \
    pval=config[NAME].c_str();                     \
    if(*pval)PARAM=atof(pval);                     \
    else LogWarning(                               \
        "Warning: Using Default " NAME ": %e\n",   \
        double(PARAM));                     }while(0)

    LOAD_CONFIG("TDB",t_eph);
    LOAD_CONFIG("Delta_t",delta_t);
    LOAD_CONFIG("Cadence",data_cadence);
    LOAD_CONFIG("Max_Ephemeris_Length",max_ephm_length);

    if(integrator==int_t(COMBINED_RK12)){
        LOAD_CONFIG("Combined_Delta_t_Max",combined_delta_t);
        LOAD_CONFIG("Combined_GM_Max_Parent",GM_max_parent);
        LOAD_CONFIG("Combined_GM_Max_Child",GM_max_child);
        LOAD_CONFIG("Combined_GM_Max_Tiny",GM_max_tiny);
        LOAD_CONFIG("Combined_Period_Max_Child",Period_max_child);
    }
#undef LOAD_CONFIG

    if(!mlist.size())
        break;

    if(int_t(t_eph.hi)!=t_eph.hi||t_eph.lo!=0){
        LogError("TDB (Initial Ephemeris Time, seconds) should be an integer.");
        break;
    }

    //correct time variables to t_eph 0 (see update())
    for(auto &m:mlist){
        m.GM0=m.GM-m.dGM*t_eph.hi;
        m.exJ2=m.dJ2*t_eph.hi;
        m.C_static.x.x-=m.exJ2/2;
        m.C_static.y.y-=m.exJ2/2;
        m.C_static.z.z+=m.exJ2;
    }

    //check config
    if(!sanity(delta_t)||!sanity(data_cadence)||!sanity(max_ephm_length)||!sanity(combined_delta_t)){
        LogError("Delta_t/Cadence/Max_Ephemeris_Length/Combined_Delta_t_Max should be finite positive real numbers.\n");
        break;
    }

    LogInfo("Loaded %lld bodies from initial %s.\n",mlist.size(),sip.c_str());
    if(fcheckpoint&&*fcheckpoint){
        ozippack zp(fcheckpoint);
        if(!zp){
            LogError("Cannot open output file. Is its directory exist?\n");
            break;
        }
        zp.resize(1);
        MFILE &mf=zp[0];
        save_checkpoint(&mf);
        mf.set_name(Configs::SaveNameCheckpoint);

        std::string
            &fbase=config["Initial"],
            &fext=config["Extra"],
            &gppath=config["Geopotentials"],
            &ringpath=config["Rings"];
        std::string initdir=std::string(Configs::SaveNameInitialDirectory)+"/";
        zp.push_back(initdir,{});

        zp.push_back(initdir+fconfig,MFILE(dirstr+fconfig));
        zp.push_back(initdir+fbase,MFILE(dirstr+fbase));
        if(fext.size())
            zp.push_back(initdir+fext,MFILE(dirstr+fext));
        if(gppath.size())
            zp.push_back(initdir+gppath+"/",{});
        if(ringpath.size())
            zp.push_back(initdir+ringpath+"/",{});
        for(const auto &m:mlist){
            if(m.gpmodel){
                std::string gname=gppath+"/"+m.get_ssid()+".txt";
                zp.push_back(initdir+gname,MFILE(dirstr+gname));
            }
            if(m.ringmodel){
                std::string rname=ringpath+"/"+m.get_ssid()+".txt";
                zp.push_back(initdir+rname,MFILE(dirstr+rname));
            }
        }
        LogInfo("Saving checkpoint %s\n",fcheckpoint);
    }
    success=true;
  }while(0);
    if(!success){
        clear();
        return false;
    }
    return true;
}

void msystem::reset_params(){
    //default value for solar system
    t_eph=0;
    delta_t=300;//5min
    integrator=int_t(COMBINED_RK12);
    data_cadence=7200;//2h
    max_ephm_length=631152000;//20yr
    combined_delta_t=28800;//8h
    GM_max_parent=2E17;
    GM_max_child=1E13;
    GM_max_tiny=14E11;
    Period_max_child=0;
}

void msystem::clear(){
    t_barycen=NAN;
    t_update=NAN;
    blist.clear();
    mlist.clear();
    midx.clear();
    tidal_childlist.clear();

    for(auto *p:gp_components)geopotential::unload(p);
    gp_components.clear();
    for(auto *p:ring_components)ring::unload(p);
    ring_components.clear();
}

bool msystem::load(
    const char *fbase,
    const char *fext,
    const char *gppath,
    const char *ringpath){
    reset_params();
    clear();

    if(!fbase)return false;
    if(!gppath||!*gppath)gppath=default_path;
    if(!ringpath||!*ringpath)ringpath=default_path;
    if(strlen(gppath)>MAX_PATHSIZE){
        LogError("Geopotentials Path too long: %s\n",gppath);
        return false;
    }
    if(strlen(ringpath)>MAX_PATHSIZE){
        LogError("Rings Path too long: %s\n",ringpath);
        return false;
    }

    bool failed=false;
    MFILE *fin=mopen(fbase);
    if(!fin){
        LogError("File Not Exist: %s\n",fbase);
        return false;
    }
    char sname[MAX_LINESIZE],sid[MAX_LINESIZE];
    while(1){
        std::string chbuf=readline(fin);
        if(chbuf.size()==0)break;
        if(chbuf.size()>=MAX_LINESIZE){
            LogError(
                "Loading %s\n Line Too Long: %s\n",
                fbase,chbuf.c_str());
            failed=true;
            break;
        }
        mass m;
        double j2,c21,c22,s21,s22;
        double ra,dec,W,p;
        double gprf,ringmf;
        fast_mpvec r,v,x,z,w;
        int n=sscanf(chbuf.c_str(),
            "%[^\t]%s%lf%lf"//Name ID GM Radius
            "%lf%lf%lf%lf"//inertia k2 k2r td
            "%lf%lf%lf%lf%lf"//Harmonics
            "%lf%lf%lf"//tdf,gprf,ringmf
            "%lf%lf%lf%lf%lf%lf"//r, v
            "%lf%lf%lf%lf"//alpha delta W period
            "%lf%lf%lf%lf%lf%lf%lf%lf%lf",//x, z, w
            sname,sid,&m.GM,&m.R,
            &m.inertia,&m.k2,&m.k2r,&m.tide_delay,
            &j2,&c21,&c22,&s21,&s22,
            &m.tide_delay_factor,&gprf,&ringmf,
            &r.x,&r.y,&r.z,&v.x,&v.y,&v.z,
            &ra,&dec,&W,&p,
            &x.x,&x.y,&x.z,&z.x,&z.y,&z.z,&w.x,&w.y,&w.z
        );
        if(!(n==35||n==26)){
            LogError(
                "Loading %s\n Line Length Error: %s\n",
                fbase,chbuf.c_str());
            failed=true;
            break;
        }

        m.sid=0;
        size_t sidlen=strlen(sid);
        if(sidlen>=8){
            LogError("Loading %s\n sid Too Long: len(%s) >= 8\n",
                fbase,sid);
            failed=true;
            break;
        }
        memcpy(&m.sid,sid,sidlen);

        if(!midx.insert({m.sid,mlist.size()}).second){
            LogError("Loading %s\n Duplicate sid: %s\n",
                fbase,sid);
            failed=true;
            break;
        }

        m.R*=1000;
        m.GM0=m.GM;
        m.dGM=0;
        m.lum=0;
        m.recpt=0;
        m.exJ2=0;
        m.dJ2=0;

        m.gpmodel=nullptr;
        m.ringmodel=nullptr;
        if(gprf!=0){
            sprintf(sname,"%s\\%s.txt",gppath,sid);
            m.gpmodel=geopotential::load(sname,gprf);
            if(!m.gpmodel){
                failed=true;
                break;
            }
            gp_components.push_back(m.gpmodel);
        }
        if(ringmf!=0){
            sprintf(sname,"%s\\%s.txt",ringpath,sid);
            m.ringmodel=ring::load(sname,m.GM,m.R,ringmf);
            if(!m.ringmodel){
                failed=true;
                break;
            }
            if(m.ringmodel->N>0)
                ring_components.push_back(m.ringmodel);
            else{
                ring::unload(m.ringmodel);
                m.ringmodel=nullptr;
            }
        }

        if(n==26){
            using Constants::pi;
            using Constants::degree;
            using Constants::J2000_obliquity;
            z=vec::from_theta_phi((90-dec)*degree,ra*degree);
            x=vec::from_theta_phi(90*degree,(ra+90)*degree);
            x+=z.rotation_matrix(W*degree)%x;
            x.rotx(-J2000_obliquity);
            z.rotx(-J2000_obliquity);
            w=z*(2*pi/(p*3600));
        }
        m.s.x=x;
        m.s.z=z;
        m.s.y=z*x;
        m.orthogonalize();
        m.r=r;
        m.v=v;
        m.w=w;
        m.A=3*m.inertia/2;
        m.R2=m.R*m.R;
        m.sR2=m.R<0?m.R2:0;
        m.rR2G_4c=0;
        m.C_static.from_harmonics(j2,c21,c22,s21,s22);

        //Calculate initial GI and GL
        fast_mpmat fmis(m.s);
        m.GI=fast_real(-2)/3*m.R2*(fmis.toworld(m.C_static)-m.A);
        m.GL=m.GI%m.w;

        mlist.push_back(m);
    }
    fclose(fin);
    if(failed)return false;

    if(fext){
        MFILE *finex=mopen(fext);
        if(!finex){
            LogError("File Not Exist: %s\n",fext);
            return false;
        }
        while(1){
            std::string chbuf=readline(finex);
            if(chbuf.size()==0)break;
            if(chbuf.size()>=MAX_LINESIZE){
                LogError(
                    "Loading %s\n Line Too Long: %s\n",
                    fext,chbuf.c_str());
                failed=true;
                break;
            }
            double param;
            int n=sscanf(chbuf.c_str(),
                "%s%s%lf",//sid Param value
                sid,sname,&param
            );

            if(n!=3){
                LogError(
                    "Loading %s\n Line Length Error: %s\n",
                    fext,chbuf.c_str());
                failed=true;
                break;
            }

            int_t mid=get_mid(sid);
            if(mid<0){
                LogWarning(
                    "Loading %s\n Ignorning %s of None-Existing sid %s\n",
                    fext,sname,sid);
                continue;
            }
            mass &m=mlist[mid];
            if     (0==strcmp(sname,"MassRate")){
                m.dGM=param;
            }
            else if(0==strcmp(sname,"Luminosity")){
                m.lum=param;
            }
            else if(0==strcmp(sname,"Receptance")){
                m.recpt=param;
                using Constants::c;
                using Constants::G;
                m.rR2G_4c=m.recpt*m.R2*(G/(4*c));
            }
            else if(0==strcmp(sname,"J2Rate")){
                m.dJ2=param;
            }
            else {
                LogError(
                    "Loading %s\n Unknown Parameter: %s\n",
                    fext,sname);
                failed=true;
                break;
            }
        }
        fclose(finex);
        if(failed)return false;
    }

    //Calculate initial deform
    deform();

    //remove deform from initial harmonics
    for(auto &m:mlist){
        fast_mpmat fmis(m.s);
        fast_mpmat ncp=fmis.toworld(m.C_static);
        m.C_static-=fmis.tolocal(m.C_potential)-m.C_static;
        m.C_potential=ncp;
        if(!m.sanity(true))
            failed=true;
    }
    if(failed)return false;

    //Calculate initial accel
    accel();

    //analyse orbit structure
    analyse(true);

    return true;
}

using Configs::CheckPointMagic;
using Configs::CheckPointVersion;

struct fcp_header{
    int_t magic;//0x53484c486d687045
    int_t version;
    real t_eph;
    fast_real delta_t;
    int_t integrator;
    fast_real data_cadence,max_ephm_length;
    fast_real combined_delta_t;
    fast_real GM_max_child,GM_max_parent,GM_max_tiny,Period_max_child;
    int_t nmass,nbarycen;
    int_t masssize;
};
struct barycen_ids{
    int_t pid,hid,gid,tid,mid;
    int_t nch;
};

bool bsystem::load_barycen_structure(MFILE *fin,size_t bsize){
    bsystem &blist=*this;
    bool failed=false;
    blist.resize(0);
    for(size_t i=0;i<bsize;++i){
        barycen &b=blist.emplace_back();
        barycen_ids bids;
        if(1!=fread(&bids,sizeof(bids),1,fin)){
            failed=true;
            break;
        }
        b.pid=bids.pid;
        b.hid=bids.hid;
        b.gid=bids.gid;
        b.tid=bids.tid;
        b.mid=bids.mid;
        for(int_t ic=0;ic<bids.nch;++ic)
            if(1!=fread(&b.children.emplace_back(),sizeof(int_t),1,fin))
                bids.nch=-1;
        if(bids.nch<0){
            failed=true;
            break;
        }
    }
    if(failed){
        blist.clear();
        return false;
    }
    return true;
}

bool msystem::load_checkpoint(MFILE *fin){
    clear();
    
    mass mwrite;
    const size_t masssize=(char*)&(mwrite.*mass_temporary_variables)-(char*)&mwrite;

    fcp_header h;

    bool failed=false;
    do{
        if(1!=fread(&h,sizeof(h),1,fin)
         ||h.magic!=CheckPointMagic
         ||h.nmass<=0
         ||h.nbarycen<=0){
            failed=true;
            break;
        }
        if(h.version!=CheckPointVersion
         ||h.masssize!=masssize){
            LogError(
                "Error: The version of checkpoint file is inconsistent with the program.\n"
                "Please use the same version as produced this checkpoint to continue.\n");
            failed=true;
            break;
        }

        t_eph=h.t_eph;
        delta_t=h.delta_t;
        integrator=h.integrator;
        data_cadence=h.data_cadence;
        max_ephm_length=h.max_ephm_length;
        combined_delta_t=h.combined_delta_t;
        GM_max_child=h.GM_max_child;
        GM_max_parent=h.GM_max_parent;
        GM_max_tiny=h.GM_max_tiny;
        Period_max_child=h.Period_max_child;

        failed=!blist.load_barycen_structure(fin,h.nbarycen);
        if(failed)break;
        
        for(int_t i=0;i<h.nmass;++i){
            mass &m=mlist.emplace_back();
            if(1!=fread(&m,masssize,1,fin)){
                failed=true;
                break;
            }
            if(m.gpmodel){
                size_t gps=(size_t)m.gpmodel;
                geopotential *gpmodel=(geopotential*)malloc(gps);
                gp_components.push_back(gpmodel);
                if(1!=fread(gpmodel,gps,1,fin)){
                    m.gpmodel=nullptr;
                    failed=true;
                    break;
                }
                m.gpmodel=gpmodel;
            }
            if(m.ringmodel){
                size_t rms=(size_t)m.ringmodel;
                ring *ringmodel=(ring*)malloc(rms);
                ring_components.push_back(ringmodel);
                if(1!=fread(ringmodel,rms,1,fin)){
                    m.ringmodel=nullptr;
                    failed=true;
                    break;
                }
                m.ringmodel=ringmodel;
            }
        }
        if(failed)break;

        if(blist.compatible_size()==h.nmass){
            blist.update_barycens(*this);
            t_update=t_eph;
        }
    } while(false);

    if(failed){
        clear();
        return false;
    }

    build_mid();
    return !mlist.empty();
}

void bsystem::save_barycen_structure(MFILE *fout) const{
    const bsystem &blist=*this;
    for(const auto &b:blist){
        barycen_ids bids;
        bids.pid=b.pid;
        bids.hid=b.hid;
        bids.gid=b.gid;
        bids.tid=b.tid;
        bids.mid=b.mid;
        bids.nch=b.children.size();
        fwrite(&bids,sizeof(bids),1,fout);
        fwrite(b.children.data(),sizeof(int_t),bids.nch,fout);
    }
}

bool msystem::save_checkpoint(MFILE *fout) const{

    mass mwrite;
    const size_t masssize=(char*)&(mwrite.*mass_temporary_variables)-(char*)&mwrite;
    
    fcp_header h;
    h.magic=CheckPointMagic;
    h.version=CheckPointVersion;
    h.t_eph=t_eph;
    h.delta_t=delta_t;
    h.integrator=integrator;
    h.data_cadence=data_cadence;
    h.max_ephm_length=max_ephm_length;
    h.combined_delta_t=combined_delta_t;
    h.GM_max_child=GM_max_child;
    h.GM_max_parent=GM_max_parent;
    h.GM_max_tiny=GM_max_tiny;
    h.Period_max_child=Period_max_child;
    h.nmass=mlist.size();
    h.nbarycen=blist.size();
    h.masssize=masssize;
    
    fwrite(&h,sizeof(h),1,fout);

    blist.save_barycen_structure(fout);

    for(const auto &m:mlist){        
        memcpy(&mwrite,&m,masssize);
        if(m.gpmodel){
            mwrite.gpmodel=(geopotential*)(m.gpmodel->size());
        }
        if(m.ringmodel){
            mwrite.ringmodel=(ring*)(m.ringmodel->size());
        }
        fwrite(&mwrite,masssize,1,fout);
        if(m.gpmodel){
            fwrite(m.gpmodel,(size_t)mwrite.gpmodel,1,fout);
        }
        if(m.ringmodel){
            fwrite(m.ringmodel,(size_t)mwrite.ringmodel,1,fout);
        }
    }

    return true;
}

void bsystem::print_structure(MFILE *mf,const msystem &ms) const{
    class barycen_structure_printer{
        const msystem &ms;
        const bsystem &blist;
    public:
        void print_structure(MFILE *mf,int_t root,int_t level) const{
            if(root<0)return;

            const barycen& br=blist[root];
            int_t cn=br.children.size();
            std::string barycen_name;

            auto print_string=[mf,level](const std::string &str){
                if(mf->tell()==0||mf->data()[mf->tell()-1]=='\n'){
                    //indent
                    std::string indent(level,' ');
                    fwrite(indent.data(),indent.size(),1,mf);
                }
                fwrite(str.data(),str.size(),1,mf);
            };

            if(br.hid<0){
                barycen_name=ms[br.mid].get_ssid();
                if(cn==0){
                    print_string("\""+barycen_name+"\"");
                    return;
                }
                print_string("{\n");
                print_string("       \"ID\":\""+barycen_name+"\"");
            }
            else{
                barycen_name="Barycenter[";
                barycen_name+=ms[blist[br.hid].mid].get_ssid();
                if(blist[br.hid].hid>=0)barycen_name+=" System";
                barycen_name+=", ";
                barycen_name+=ms[blist[br.gid].mid].get_ssid();
                if(blist[br.gid].hid>=0)barycen_name+=" System";
                barycen_name+="]";
                print_string("{\n");
                print_string("       \"ID\":\""+barycen_name+"\",\n");
                print_string("  \"Primary\":");
                print_structure(mf,br.hid,level+12);
                print_string(",\n");
                print_string("\"Secondary\":");
                print_structure(mf,br.gid,level+12);
            }

            if(cn){
                print_string(",\n");
                print_string(" \"Children\":[\n");
                for(int_t i=0;i<cn;){
                    print_structure(mf,br.children[i],level+16);
                    print_string((++i==cn)+",\n");
                }

                print_string("            ]\n");
            }
            else{
                print_string("\n");
            }
            print_string("}");
        }
        barycen_structure_printer(const msystem &_ms,const bsystem &_blist):ms(_ms),blist(_blist){}
    };
    barycen_structure_printer(ms,*this).print_structure(mf,root_id(),0);
}

size_t bsystem::compatible_size() const{
    const bsystem &blist=*this;
    int_t root=root_id(),bn=blist.size(),mn=bn;
    if(root<0)return false;
    for(const barycen &b:blist)
        mn-=b.hid>=0;

    std::vector<bool> has_done(bn,false);
    std::vector<bool> has_mass(mn,false);
    std::vector<bool> mid_done(bn,false);
    std::deque<std::pair<int_t,int_t>> active;
    active.emplace_back(root,-1);
    while(!active.empty()){
        auto cur=active.front();
        int_t bid=cur.first;
        active.pop_front();
        if(bid<0||bid>=bn||has_done[bid])return false;
        has_done[bid]=true;
        const barycen &b=blist[bid];
        if(cur.second!=b.pid)return false;
        bool is_barycen=b.hid>=0;
        if(is_barycen!=(b.gid>=0))return false;
        if(is_barycen){
            active.emplace_back(b.hid,bid);
            active.emplace_back(b.gid,bid);
        }
        else{
            if(b.mid<0||b.mid>=mn||has_mass[b.mid]||mid_done[bid])return false;
            has_mass[b.mid]=true;
            mid_done[bid]=true;
        }
        int_t ctid=bid;
        if(b.pid>=0){
            const barycen &p=blist[b.pid];
            bool is_host=p.hid==bid;
            if(is_host){
                ctid=p.tid;
                if(b.mid!=p.mid||mid_done[b.pid])return false;
                mid_done[b.pid]=true;
            }
        }
        if(b.tid!=ctid)return false;
        for(int_t c:b.children)active.emplace_back(c,bid);
    }

    for(bool _:has_done)if(!_)return false;
    for(bool _:mid_done)if(!_)return false;
    for(bool _:has_mass)if(!_)return false;
    return mn;
}

int_t bsystem::root_id() const{
    const bsystem &blist=*this;
    int_t bn=blist.size();
    int_t rootid=-1;
    int_t nroot=0;
    for(int_t i=0;i<bn;++i){
        if(blist[i].pid<0){
            rootid=i;
            ++nroot;
        }
    }
    return nroot==1?rootid:-1;
}

void msystem::build_mid(){
    midx.clear();
    size_t mn=mlist.size();
    for(size_t i=0;i<mn;++i){
        const mass &mi=mlist[i];
        midx.insert({mi.sid,i});
    }
}
int_t msystem::get_mid(const char *ssid) const{
    size_t slen=strlen(ssid);
    uint64_t isid=0;
    const size_t maxslen=sizeof(isid)-1;
    if(slen>maxslen)return -1;
    memcpy(&isid,ssid,slen);
    return get_mid(isid);
}
int_t msystem::get_mid(uint64_t sid) const{
    auto it=midx.find(sid);

    do{
        if(it==midx.end())break;
        size_t result=it->second;
        if(result>=mlist.size())break;
        if(mlist[result].sid!=sid)break;
        //correct
        return result;
    } while(0);

    return -1;
}

void msystem::copy_params(const msystem &other){
    if(this==&other)return;
#define copy_member(_m) _m=other._m
    copy_member(t_eph);
    copy_member(delta_t);
    copy_member(integrator);
    copy_member(data_cadence);
    copy_member(max_ephm_length);
    copy_member(combined_delta_t);
    copy_member(GM_max_child);
    copy_member(GM_max_parent);
    copy_member(GM_max_tiny);
    copy_member(Period_max_child);
}
msystem &msystem::operator =(const msystem &other){
    if(this==&other)return *this;
    copy_params(other);
    clear();
    copy_member(tidal_parent);
    copy_member(tidal_matrix);
    copy_member(tidal_childlist);
    copy_member(p_substeper);
    copy_member(t_barycen);
    copy_member(t_update);
    copy_member(blist);
    copy_member(mlist);
#undef copy_member
    build_mid();
    for(mass &m:mlist){
        if(m.gpmodel)
            gp_components.push_back(m.gpmodel=geopotential::copy(m.gpmodel));
        if(m.ringmodel)
            ring_components.push_back(m.ringmodel=ring::copy(m.ringmodel));
    }
    return *this;
}
bool msystem::is_same(const msystem &other){
    if(this==&other)return true;
#define compare_member(_m) if(memcmp(&_m,&other._m,sizeof(_m)))return false
    compare_member(delta_t);
    compare_member(integrator);
    compare_member(data_cadence);
    compare_member(max_ephm_length);
    compare_member(combined_delta_t);
    compare_member(GM_max_child);
    compare_member(GM_max_parent);
    compare_member(GM_max_tiny);
    compare_member(Period_max_child);
    if(tidal_childlist!=other.tidal_childlist)return false;
    if(!tidal_childlist.empty()){
        compare_member(tidal_parent);
        compare_member(tidal_matrix);
#undef compare_member
    }
    int_t mn=mlist.size();
    if(other.mlist.size()!=mn)return false;
    for(int_t i=0;i<mn;++i){
        const mass &mi=mlist[i];
        const mass &mj=other.mlist[i];
        const char *pbegin=(const char*)&(mi.*mass_constant_parameters_begin);
        const char *pother=(const char*)&(mj.*mass_constant_parameters_begin);
        const char *pend=(const char*)&(mi.*mass_constant_parameters_end);
        if(memcmp(pbegin,pother,pend-pbegin)
            ||!geopotential::is_same(mi.gpmodel,mj.gpmodel)
            ||!ring::is_same(mi.ringmodel,mj.ringmodel))
            return false;
    }
    return true;
}
bool msystem::push_back(const mass &msrc){
    auto result=midx.insert({msrc.sid,mlist.size()});
    if(!result.second)return false;
    mlist.push_back(msrc);
    mass &m=mlist.back();
    if(m.gpmodel)
        gp_components.push_back(m.gpmodel=geopotential::copy(m.gpmodel));
    if(m.ringmodel)
        ring_components.push_back(m.ringmodel=ring::copy(m.ringmodel));
    return true;
}

void msystem::scale(int_t mid,fast_real factor){
    mlist[mid].scale(factor);
}
void msystem::scale_geopotential(int_t mid,fast_real factor){
    mass &mi=mlist[mid];
    auto &pmodel=mi.gpmodel;
    if(pmodel=geopotential::copy(pmodel,factor))
        gp_components.push_back(pmodel);
    mi.C_potential*=factor;
    mi.C_static*=factor;
    mi.k2*=factor;
    mi.k2r*=factor;
    mi.exJ2*=factor;
    mi.dJ2*=factor;
    mi.GI=fast_real(-2)/3*mi.R2*(mi.C_potential-mi.A);
    mi.GL=mi.GI%mi.w;
}
void msystem::scale_ring(int_t mid,fast_real factor){
    auto &pmodel=mlist[mid].ringmodel;
    if(pmodel=ring::copy(pmodel,factor))
        ring_components.push_back(pmodel);
}
