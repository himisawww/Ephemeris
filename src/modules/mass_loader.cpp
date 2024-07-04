#include"physics/CelestialSystem.h"
#include"physics/geopotential.h"
#include"physics/ring.h"
#include"utils/zipio.h"
#include"configs.h"

#define MAX_LINESIZE 1024
#define MAX_PATHSIZE 260

const char *default_path=".";

bool sanity(double dt){
    return dt>0&&dt!=INFINITY;
}

bool msystem::load(const char *fconfig,const char *fcheckpoint){
    const std::string sip=fconfig;
    size_t spos=1+sip.find_last_of("/\\");
    std::map<std::string,std::string> config;
    std::string dirstr(spos?sip.substr(0,spos):".\\");
    std::string fconfigstr(sip.substr(spos));

    const char *dir=dirstr.c_str();
    fconfig=fconfigstr.c_str();

    if(!dir||!*dir)dir=default_path;
    if(!fconfig||!*fconfig)return false;
    if(strlen(dir)+strlen(fconfig)+1>=MAX_PATHSIZE){
        fprintf(stderr,"Path too long: %s\\%s\n",dir,fconfig);
        return false;
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
        fprintf(stderr,"File Not Exist: %s\n",fname);
        return false;
    }

    while(1){
        std::string chbuf=readline(fin);
        if(chbuf.size()==0)break;
        if(chbuf.size()>=MAX_LINESIZE){
            fprintf(stderr,
                "Loading %s\n Line Too Long: %s\n",
                fname,chbuf.c_str());
            failed=true;
            break;
        }

        int n=sscanf(chbuf.c_str(),"%s%s",sname,sval);
        if(n!=2){
            fprintf(stderr,
                "Loading %s\n Line Length Error: %s\n",
                fname,chbuf.c_str());
            failed=true;
            break;
        }

        auto it=config.find(sname);
        if(it==config.end()){
            fprintf(stderr,
                "Loading %s\n Unknown Parameter: %s\n",
                fname,sname);
            failed=true;
            break;
        }

        if(it->second.size()){
            fprintf(stderr,
                "Loading %s\n Duplicate Parameter: %s\n",
                fname,sname);
            failed=true;
            break;
        }
        
        it->second=sval;
    }
    fclose(fin);
    if(failed)return false;

    std::string
        &fbase      =config["Initial"],
        &fext       =config["Extra"],
        &gppath     =config["Geopotentials"],
        &ringpath   =config["Rings"];

    if(!fbase.size()){
        fprintf(stderr,
            "Loading %s\n Missing Parameter: Initial\n",
            fname);
        return false;
    }
    
    sprintf(sname,"%s\\",dir);

    if(!load((sname+fbase).c_str(),
        fext.size()?(sname+fext).c_str():nullptr,
        gppath.size()?(sname+gppath).c_str():nullptr,
        ringpath.size()?(sname+ringpath).c_str():nullptr
        ))return false;

    const char *pval;

    pval=config["Integrator"].c_str();
    if(*pval){
        integrator=-1;
        const int integrators_size=sizeof(integrators_list)/sizeof(char*);
        for(int i=0;i<integrators_size;++i){
            if(0==strcmp(pval,integrators_list[i]))integrator=i;
        }
        if(integrator<0){
            fprintf(stderr,
                "Loading %s\n Invalid Integrator: %s\n",
                fname,pval);
            return false;
        }
    }
    else fprintf(stderr,
        "Warning: Using Default Integrator: %s\n",
        integrators_list[integrator]);

#define LOAD_CONFIG(NAME,PARAM) do{                \
    pval=config[NAME].c_str();                     \
    if(*pval)PARAM=atof(pval);                     \
    else fprintf(stderr,                           \
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
        return false;

    //check config
    if(!sanity(delta_t)||!sanity(data_cadence)||!sanity(max_ephm_length)||!sanity(combined_delta_t)){
        fprintf(stderr,"Delta_t/Cadence/Max_Ephemeris_Length/Combined_Delta_t_Max should be finity positive real numbers.\n");
        return false;
    }

    printf("Loaded %lld bodies from initial %s.\n",mlist.size(),sip.c_str());
    if(fcheckpoint&&*fcheckpoint){
        ozippack zp(fcheckpoint);
        if(!zp){
            fprintf(stderr,"Cannot open output file. Is its directory exist?\n");
            return false;
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
        printf("Saving checkpoint %s\n",fcheckpoint);
    }
    
    return true;
}

bool msystem::load(
    const char *fbase,
    const char *fext,
    const char *gppath,
    const char *ringpath){

    //defalt value for solar system
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

    blist.clear();
    mlist.clear();
    midx.clear();
    tidal_childlist.clear();

    if(!fbase)return false;
    if(!gppath||!*gppath)gppath=default_path;
    if(!ringpath||!*ringpath)ringpath=default_path;
    if(strlen(gppath)>MAX_PATHSIZE){
        fprintf(stderr,"Geopotentials Path too long: %s\n",gppath);
        return false;
    }
    if(strlen(ringpath)>MAX_PATHSIZE){
        fprintf(stderr,"Rings Path too long: %s\n",ringpath);
        return false;
    }

    bool failed=false;
    MFILE *fin=mopen(fbase);
    if(!fin){
        fprintf(stderr,"File Not Exist: %s\n",fbase);
        return false;
    }
    char sname[MAX_LINESIZE],sid[MAX_LINESIZE];
    while(1){
        std::string chbuf=readline(fin);
        if(chbuf.size()==0)break;
        if(chbuf.size()>=MAX_LINESIZE){
            fprintf(stderr,
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
            fprintf(stderr,
                "Loading %s\n Line Length Error: %s\n",
                fbase,chbuf.c_str());
            failed=true;
            break;
        }

        m.sid=0;
        size_t sidlen=strlen(sid);
        if(sidlen>=8){
            fprintf(stderr,"Loading %s\n sid Too Long: len(%s) >= 8\n",
                fbase,sid);
            failed=true;
            break;
        }
        memcpy(&m.sid,sid,sidlen);

        if(!midx.insert({m.sid,mlist.size()}).second){
            fprintf(stderr,"Loading %s\n Duplicate sid: %s\n",
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
        }
        if(ringmf!=0){
            sprintf(sname,"%s\\%s.txt",ringpath,sid);
            m.ringmodel=ring::load(sname,m.GM,m.R2,ringmf);
            if(!m.ringmodel){
                failed=true;
                break;
            }
            if(m.ringmodel->N==0){
                ring::unload(m.ringmodel);
                m.ringmodel=nullptr;
            }
        }

        if(n==26){
            using Constants::pi;
            using Constants::degree;
            z=vec((90-dec)*degree,ra*degree);
            x=vec(90*degree,(ra+90)*degree);
            x+=rotation_matrix(z,W*degree)%x;
            x.rotx(-84381.448/3600*degree);
            z.rotx(-84381.448/3600*degree);
            w=z*(2*pi/(p*3600));
        }
        m.s.x=x;
        m.s.z=z;
        m.s.y=z*x;
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
        m.GI=2*m.R2/3*(fast_mpmat(m.A)-fmis.toworld(m.C_static));
        m.GL=m.GI%m.w;

        mlist.push_back(m);
    }
    fclose(fin);
    if(failed)return false;

    if(fext){
        MFILE *finex=mopen(fext);
        if(!finex){
            fprintf(stderr,"File Not Exist: %s\n",fext);
            return false;
        }
        while(1){
            std::string chbuf=readline(finex);
            if(chbuf.size()==0)break;
            if(chbuf.size()>=MAX_LINESIZE){
                fprintf(stderr,
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
                fprintf(stderr,
                    "Loading %s\n Line Length Error: %s\n",
                    fext,chbuf.c_str());
                failed=true;
                break;
            }

            int_t mid=get_mid(sid);
            if(mid<0){
                fprintf(stderr,
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
                fprintf(stderr,
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
    analyse();

    return true;
}

static const int_t CheckPoint_Magic=0x53484c486d687045;
static const int_t Version_Number=2;
//contents in mass is irrelevant after this parameter
#define Mass_Auxiliary_Head Erot

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

bool msystem::load_checkpoint(MFILE *fin){
    blist.clear();
    mlist.clear();
    midx.clear();
    tidal_childlist.clear();
    
    mass mwrite;
    const size_t masssize=(char*)&mwrite.Mass_Auxiliary_Head-(char*)&mwrite;

    fcp_header h;

    bool failed=false;
    do{
        if(1!=fread(&h,sizeof(h),1,fin)
         ||h.magic!=CheckPoint_Magic){
            failed=true;
            break;
        }
        if(h.version!=Version_Number
         ||h.masssize!=masssize){
            fprintf(stderr,
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
        if(failed)break;
        
        for(auto &m:mlist){
            if(1!=fread(&m,masssize,1,fin)){
                failed=true;
                break;
            }
            if(m.gpmodel){
                size_t gps=(size_t)m.gpmodel;
                m.gpmodel=(geopotential*)malloc(gps);
                if(1!=fread(m.gpmodel,gps,1,fin)){
                    failed=true;
                    break;
                }
            }
            if(m.ringmodel){
                size_t rms=(size_t)m.ringmodel;
                m.ringmodel=(ring*)malloc(rms);
                if(1!=fread(m.ringmodel,rms,1,fin)){
                    failed=true;
                    break;
                }
            }
        }
        if(failed)break;

        //Calculate initial accel
        //accel();
        update_barycens();
    } while(false);

    return !failed;
}

bool msystem::save_checkpoint(MFILE *fout){

    mass mwrite;
    const size_t masssize=(char*)&mwrite.Mass_Auxiliary_Head-(char*)&mwrite;
    
    fcp_header h;
    h.magic=CheckPoint_Magic;
    h.version=Version_Number;
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
