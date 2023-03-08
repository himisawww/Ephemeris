#include"ephemeris.h"
#include"utils.h"

#define MAX_LINESIZE 1024

std::string readline(FILE *fin){
    std::string result;
    char chbuf[16];
    char *fret;
    do{
        while(fret=fgets(chbuf,sizeof(chbuf),fin)){
            result+=chbuf;
            if(result.back()=='\n'){
                result.pop_back();
                break;
            }
        }

        auto cpos=result.find('#');
        if(cpos==0)
            result.resize(cpos);
        
        if(result.size()==0&&fret)
            continue;

        return result;
    } while(1);
}

bool msystem::load(
    const char *fbase,
    const char *fext,
    const char *gppath,
    const char *ringpath){

    t_eph=0;
    delta_t=300;//5min
    node=96;//28800s = 8h
    page=4383;//4yr

    blist.clear();
    mlist.clear();
    midx.clear();
    tidal_childlist.clear();

    //load solar system
    bool failed=false;
    FILE *fin=fopen(fbase,"r");
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
        const double c=299792458;
        const double G=6.67430E-11;
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
        fin=fopen(fext,"r");
        if(!fin){
            fprintf(stderr,"File Not Exist: %s\n",fext);
            return false;
        }
        while(1){
            std::string chbuf=readline(fin);
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
        fclose(fin);
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
//contents in mass is irrelevant after this parameter
#define Mass_Auxiliary_Head Erot

struct fcp_header{
    int_t magic;//0x53484c486d687045
    int_t masssize;
    int_t nmass,nbarycen;
    real t_eph,delta_t;
    int_t node,page;
};
struct barycen_ids{
    int_t pid,hid,gid,tid,mid;
    int_t nch;
};

bool msystem::load_checkpoint(const char *fcp){
    blist.clear();
    mlist.clear();
    midx.clear();
    tidal_childlist.clear();
    
    FILE *fin=fopen(fcp,"rb");
    if(!fin){
        fprintf(stderr,"File Not Exist: %s\n",fcp);
        return false;
    }

    mass mwrite;
    const size_t masssize=(char*)&mwrite.Mass_Auxiliary_Head-(char*)&mwrite;

    fcp_header h;

    bool failed=false;
    do{
        if(1!=fread(&h,sizeof(h),1,fin)
         ||h.magic!=CheckPoint_Magic
         ||h.masssize!=masssize){
            failed=true;
            break;
        }
        mlist.resize(h.nmass);
        blist.resize(h.nbarycen);
        t_eph=h.t_eph;
        delta_t=h.delta_t;
        node=h.node;
        page=h.page;

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

    fclose(fin);
    return !failed;
}

bool msystem::save_checkpoint(const char *fcp){
    FILE *fout=fopen(fcp,"wb");
    if(!fout){
        fprintf(stderr,"Cannot open: %s\n",fcp);
        return false;
    }

    mass mwrite;
    const size_t masssize=(char*)&mwrite.Mass_Auxiliary_Head-(char*)&mwrite;
    
    fcp_header h;
    h.magic=CheckPoint_Magic;
    h.masssize=masssize;
    h.nmass=mlist.size();
    h.nbarycen=blist.size();
    h.t_eph=t_eph;
    h.delta_t=delta_t;
    h.node=node;
    h.page=page;
    
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
    fclose(fout);
    return true;
}
