#include"physics/mass.h"
#include"physics/geopotential.h"
#include"physics/ring.h"
#include"utils/zipio.h"
#include"configs.h"
#include"utils/logger.h"

#define MAX_LINESIZE 1024
#define MAX_PATHSIZE 260

const char *default_path=".";

bool sanity(double dt){
    return dt>0&&dt!=INFINITY;
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
                zp.push_back(
                    initdir+gppath+"/"+(const char *)&m.sid+".txt",
                    MFILE(dirstr+gppath+"\\"+(const char *)&m.sid+".txt"));
            }
            if(m.ringmodel){
                zp.push_back(
                    initdir+ringpath+"/"+(const char *)&m.sid+".txt",
                    MFILE(dirstr+ringpath+"\\"+(const char *)&m.sid+".txt"));
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
    data_cadence=3600;//1h
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
            z=vec((90-dec)*degree,ra*degree);
            x=vec(90*degree,(ra+90)*degree);
            x+=z.rotation_matrix(W*degree)%x;
            x.rotx(-84381.448/3600*degree);
            z.rotx(-84381.448/3600*degree);
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
        using Constants::c;
        using Constants::G;
        m.rR2_4Mc=m.recpt*m.R2/(4*m.GM*c/G);
        m.C_static=fast_mpmat(
            fast_mpvec(3*c22+j2/2,3*s22,3*c21/2),
            fast_mpvec(3*s22,-3*c22+j2/2,3*s21/2),
            fast_mpvec(3*c21/2,3*s21/2,-j2)
        );

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
    }

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

bool msystem::load_barycen_structure(MFILE *fin,std::vector<barycen> &blist){
    bool failed=false;
    for(auto &b:blist){
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
        b.children.resize(bids.nch);
        if(bids.nch!=fread(b.children.data(),sizeof(int_t),bids.nch,fin)){
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
         ||h.magic!=CheckPointMagic){
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
        mlist.resize(h.nmass);
        blist.resize(h.nbarycen);
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

        failed=!load_barycen_structure(fin,blist);
        if(failed)break;
        
        for(auto &m:mlist){
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

    } while(false);

    if(failed){
        clear();
        return false;
    }

    build_mid();
    return true;
}

void msystem::save_barycen_structure(MFILE *fout,const std::vector<barycen> &blist){
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

    save_barycen_structure(fout,blist);

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

void print_string(MFILE *mf,int_t level,const std::string &str){
    if(mf->tell()==0||mf->data()[mf->tell()-1]=='\n'){
        //indent
        std::string indent(level,' ');
        fwrite(indent.data(),indent.size(),1,mf);
    }
    fwrite(str.data(),str.size(),1,mf);
}
void msystem::print_structure(MFILE *mf,int_t root,int_t level) const{
    const msystem &ms=*this;
    if(root<0){
        root=ms.blist.size();
        do{
            if(--root<0)return;
        } while(!(ms.blist[root].pid<0));
    }

    const barycen& br=ms.blist[root];
    int_t cn=br.children.size();
    std::string barycen_name;

    if(br.hid<0){
        barycen_name=(const char *)&ms.mlist[br.mid].sid;
        if(cn==0){
            print_string(mf,level,"\""+barycen_name+"\"");
            return;
        }
        print_string(mf,level,"{\n");
        print_string(mf,level,"       \"ID\":\""+barycen_name+"\"");
    }
    else{
        barycen_name="Barycenter[";
        barycen_name+=(const char *)&ms.mlist[ms.blist[br.hid].mid].sid;
        if(ms.blist[br.hid].hid>=0)barycen_name+=" System";
        barycen_name+=", ";
        barycen_name+=(const char *)&ms.mlist[ms.blist[br.gid].mid].sid;
        if(ms.blist[br.gid].hid>=0)barycen_name+=" System";
        barycen_name+="]";
        print_string(mf,level,"{\n");
        print_string(mf,level,"       \"ID\":\""+barycen_name+"\",\n");
        print_string(mf,level,"  \"Primary\":");
        print_structure(mf,br.hid,level+12);
        print_string(mf,level,",\n");
        print_string(mf,level,"\"Secondary\":");
        print_structure(mf,br.gid,level+12);
    }

    if(cn){
        print_string(mf,level,",\n");
        print_string(mf,level," \"Children\":[\n");
        for(int_t i=0;i<cn;){
            print_structure(mf,br.children[i],level+16);
            print_string(mf,level,(++i==cn)+",\n");
        }

        print_string(mf,level,"            ]\n");
    }
    else{
        print_string(mf,level,"\n");
    }
    print_string(mf,level,"}");
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
    if(pmodel)
        gp_components.push_back(pmodel=geopotential::copy(pmodel,factor));
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
    if(pmodel)
        ring_components.push_back(pmodel=ring::copy(pmodel,factor));
}
