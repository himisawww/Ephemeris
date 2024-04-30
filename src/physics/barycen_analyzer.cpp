#include"physics/CelestialSystem.h"
#include<map>

barycen::barycen(){
    pid=-1;
    hid=-1;
    gid=-1;
    tid=-1;
    mid=-1;
}

struct bdata{
    fast_real minr,factor;
    fast_real maxinfl;
    bool del;
};

bool msystem::analyse(bool reconstruct){
    int_t old_bn=reconstruct?0:blist.size();
    int_t mn=mlist.size();
    std::vector<barycen> bl;
    std::vector<bdata> dl;
    if(old_bn){
        bl.resize(old_bn);
        dl.resize(old_bn);
        for(int_t i=0;i<old_bn;++i){
            barycen &bi=bl[i];
            bdata &di=dl[i];
            bi.pid=blist[i].pid;
            bi.hid=blist[i].hid;
            bi.gid=blist[i].gid;
            bi.mid=bi.hid<0?blist[i].mid:-1;
            di.del=false;
        }
    }
    else{
        bl.resize(mn);
        dl.resize(mn);
        for(int_t i=0;i<mn;++i){
            barycen &bi=bl[i];
            bdata &di=dl[i];
            bi.pid=-1;
            bi.hid=-1;
            bi.gid=-1;
            bi.mid=i;
            di.del=false;
        }
    }

    bool diff=false;
    do{
        if(diff)old_bn=0;
        int_t bn=bl.size();

        //update barycens: bl.rvGM,rvGM_sys;
        for(int_t i=0;i<bn;++i)if(!dl[i].del){//initialize
            bl[i].GM_sys=0;
            bl[i].r_sys=0;
            bl[i].v_sys=0;
        }
        for(int_t i=0;i<bn;++i)if(bl[i].mid>=0){//is mass
            barycen &bi=bl[i];
            bdata &di=dl[i];
            const mass &mi=mlist[bi.mid];
            const real mGM=mi.GM;
            const mpvec rGM=mi.r*mGM;
            const mpvec vGM=mi.v*mGM;
            bi.r=mi.r;
            bi.v=mi.v;
            bi.GM=mGM;
            bi.r_sys+=rGM;
            bi.v_sys+=vGM;
            bi.GM_sys+=mGM;
            int_t bp=bi.pid;
            int_t ncount=0;
            while(bp>=0&&!dl[bp].del){
                if(++ncount>=bl.size())return false;
                bl[bp].GM_sys+=mGM;
                bl[bp].r_sys+=rGM;
                bl[bp].v_sys+=vGM;
                bp=bl[bp].pid;
            }
        }
        for(int_t i=0;i<bn;++i)if(bl[i].mid<0&&!dl[i].del){//is barycen
            barycen &bi=bl[i];
            barycen &h=bl[bi.hid];
            barycen &g=bl[bi.gid];
            bi.GM=h.GM_sys+g.GM_sys;
            bi.r=(h.r_sys+g.r_sys)/bi.GM;
            bi.v=(h.v_sys+g.v_sys)/bi.GM;
        }
        for(int_t i=0;i<bn;++i)if(!dl[i].del){//finalize
            barycen &bi=bl[i];
            bi.r_sys/=bi.GM_sys;
            bi.v_sys/=bi.GM_sys;
        }

        //update dl.minr&factor
        for(int_t i=0;i<bn;++i)if(!dl[i].del){
            barycen &bi=bl[i];
            bdata &di=dl[i];
            if(bi.mid>=0){//is mass
                di.minr=0;
            }
            else{//is barycen
                fast_mpvec r=bl[bi.gid].r_sys-bl[bi.hid].r_sys;
                di.minr=r.norm();
                fast_real gmh=fast_real(bl[bi.hid].GM_sys);
                fast_real gmg=fast_real(bl[bi.gid].GM_sys);
                fast_real q=gmh>gmg?gmg/gmh:gmh/gmg;
                di.minr/=1-1.0001*(std::cbrt(q)-q);
                q+=1;q*=q;
                di.factor=1/q;
            }
        }

        int_t n_root=0,root;
        //update influences: dl.maxinfl
        for(int_t i=0;i<bn;++i)if(!dl[i].del){
            barycen &bi=bl[i];
            bdata &di=dl[i];
            if(bi.pid<0||dl[bi.pid].del){//no parent
                bi.pid=-1;
                di.maxinfl=0;
                ++n_root;
                root=i;
            }
            else{//have parent
                barycen &p=bl[bi.pid];
                if(p.hid==i||p.gid==i){//is one of binary
                    fast_real r1=dl[bi.pid].minr,r2=r1*r1;
                    if(r1<dl[p.hid].minr||r1<dl[p.gid].minr){//too close to be a companion of binary
                        bi.pid=-1;
                        di.maxinfl=0;
                    }
                    else{
                        di.maxinfl=fast_real(p.GM)/(r2*r1);
                        di.maxinfl*=2;//make binaries more stable. does it make sense?
                    }
                }
                else{
                    fast_mpvec r=bi.r_sys-p.r;
                    fast_real r2=r%r,r1=sqrt(r2);
                    if(r1<dl[bi.pid].minr||r1<dl[i].minr){//too close to be a child of binary
                        bi.pid=-1;
                        di.maxinfl=0;
                    }
                    else{
                        di.maxinfl=fast_real(p.GM)/(r2*r1);
                    }
                }
            }
            if(old_bn)di.maxinfl*=2;//make existing stucture more stable.
        }
        if(n_root!=1)root=-1;

        diff=false;
        //find maximum influences
        for(int_t i=0;i<bn;++i)if(i!=root){
            barycen &bi=bl[i];
            bdata &di=dl[i];
            if(di.del)continue;
            int_t pid=bi.pid;
            int_t cid=-1;
            if(pid>=0&&bl[pid].mid<0){//parent is a barycen
                if(bl[pid].hid==i)cid=bl[pid].gid;
                if(bl[pid].gid==i)cid=bl[pid].hid;
            }
            for(int_t j=0;j<bn;++j)if(i!=j&&pid!=j&&cid!=j){//pass self, parent and companion
                barycen &bj=bl[j];
                bdata &dj=dl[j];
                if(dj.del)continue;
                fast_mpvec r=bj.r-bi.r;
                fast_real r2=r%r,r1=sqrt(r2);
                if(r1<di.minr||r1<dj.minr)continue;
                fast_real infl=fast_real(bj.GM)/(r2*r1);
                int_t tj=j,tp=bj.pid;
                // if j is one of a binary-system, and i is distant from the system-center
                // j's influence should be reduced by a factor
                while(tp>=0){
                    barycen &bp=bl[tp];
                    if(!(bp.mid<0&&(bp.hid==tj||bp.gid==tj)))break;
                    fast_mpvec dr(bp.r_sys-bi.r);
                    fast_real minr=dl[tp].minr;
                    if(dr%dr>minr*minr)
                        infl*=dl[tp].factor;
                    tj=tp;
                    tp=bp.pid;
                }
                if(infl>di.maxinfl){
                    di.maxinfl=infl;
                    bi.pid=j;
                    diff=true;
                }
            }
        }
        //delete unbounded barycens
        for(int_t i=mn;i<bn;++i){
            barycen &bi=bl[i];
            bdata &di=dl[i];
            if(di.del)continue;
            if(bl[bi.hid].pid!=i||bl[bi.gid].pid!=i){
                di.del=true;
                diff=true;
                int_t ti=i,tp=bi.pid;
                while(tp>=0&&(bl[tp].hid==ti||bl[tp].gid==ti)){
                    dl[tp].del=true;
                    ti=tp;
                    tp=bl[ti].pid;
                }
            }
        }
        //create new barycens
        for(int_t i=0;i<bn;++i)if(!dl[i].del){
            int_t j=bl[i].pid;
            if(j<0||dl[j].del||j>=bn)continue;
            barycen &bi=bl[i];
            barycen &bj=bl[j];
            if(bj.pid!=i){
                if((fast_real)bj.GM_sys<(fast_real)bi.GM_sys){
                    bi.pid=-1;
                    diff=true;
                }
                continue;
            }
            if(j<i)continue;
            int_t k=bl.size();
            bi.pid=k;
            bj.pid=k;
            bl.push_back(barycen());
            dl.push_back(bdata());
            barycen &bk=bl[k];
            bdata &dk=dl[k];
            bk.pid=-1;
            bk.hid=i;
            bk.gid=j;
            bk.mid=-1;
            dk.del=false;
            diff=true;
        }
    } while(diff);

    if(old_bn)return false;

    //remove deleted barycens
    std::vector<int_t> newindex(bl.size()+1,-1);
    int_t bsize=0;
    for(int_t i=0;i<bl.size();++i)if(!dl[i].del){
        newindex[i+1]=bsize;
        if(i!=bsize)bl[bsize]=bl[i];
        ++bsize;
    }
    bl.resize(bsize);
    for(int_t i=0;i<bsize;++i){
        barycen &bi=bl[i];
        bi.pid=newindex[bi.pid+1];
        bi.hid=newindex[bi.hid+1];
        bi.gid=newindex[bi.gid+1];
    }

    //adjust host & guest, fill children list
    for(int_t i=0;i<bsize;++i){
        barycen &bi=bl[i];
        if(bi.pid>=0){
            barycen &p=bl[bi.pid];
            if(!(p.hid==i||p.gid==i))
                p.children.push_back(i);
        }
        if(bi.hid>=0){
            if((fast_real)bl[bi.hid].GM_sys<(fast_real)bl[bi.gid].GM_sys){
                int_t temp=bi.hid;
                bi.hid=bi.gid;
                bi.gid=temp;
            }
        }
    }

    //fill tid and mid
    for(int_t i=0;i<bsize;++i){
        barycen &bi=bl[i];

        int_t oid=i,pid=bi.pid;
        while(pid>=0&&bl[pid].hid==oid){
            oid=pid;
            pid=bl[pid].pid;
        }
        bi.tid=oid;

        if(bi.mid<0){
            int_t hid=i;
            while(bl[hid].hid>=0){
                hid=bl[hid].hid;
            }
            bi.mid=bl[hid].mid;
        }
    }

    blist.swap(bl);
    return true;
}
/*
//combine minor bodies to majors
int_t msystem::combine(int method){
    if(method<0)return 0;

    int_t bn=blist.size();
    std::map<int_t,int_t> clist;
    std::map<int_t,std::vector<int_t>> cvecs;
    for(int_t id=0;id<bn;++id){
        const auto &b=blist[id];
        if(b.hid>=0)continue;

        const auto &t=blist[b.tid];
        if(t.pid<0)continue;
        const auto &tp=blist[t.pid];
        //combine target = tp.mid;

        //combine threshold (only for Solar System)
        //GM < 1E13 ( Ganymede < * < Mercury )
        //target GM < 2E17 ( Jupiter < * < Sun )
        if((fast_real)t.GM_sys>1E13
         ||(fast_real)tp.GM>2E17)continue;

        clist.insert({b.mid,tp.mid});
    }

    do{
        bool changed=false;
        for(auto &c:clist){
            auto fs=clist.find(c.second);
            if(fs==clist.end())continue;

            c.second=fs->second;
            changed=true;
        }
        if(!changed)break;
    } while(1);

    for(const auto &c:clist){
        auto is=cvecs.insert({c.second,std::vector<int_t>()});
        is.first->second.push_back(c.first);
    }

    std::vector<mass> ml;
    int_t mn=mlist.size();
    for(int_t i=0;i<mn;++i){
        if(clist.find(i)!=clist.end())continue;
        mass &mi=mlist[i];
        auto fs=cvecs.find(i);
        if(fs!=cvecs.end()){
            //start combine process
            if(method==0){
                //do nothing
            }
            else if(method==1){
                //newtonian combination
                real GMi=mi.GM;
                real GM=mi.GM0;
                mpvec r=mi.r*GM,v=mi.v*GM;
                for(auto &j:fs->second){
                    mass &mj=mlist[j];
                    real GMj=mj.GM;
                    r+=mj.r*GMj;
                    v+=mj.v*GMj;
                    GMi+=GMj;
                    GM+=mj.GM0;

                    mi.dGM+=mj.dGM;
                    mi.lum+=mj.lum;
                }
                r/=GM;
                v/=GM;
                mi.r=r;
                mi.v=v;
                mi.GM0=(fast_real)GM;
                mi.GM=(fast_real)GMi;
            }
            else if(method==2){
                //more complex considerations
                //relativity: 20cm for earth, not too large ...
                //GL, w, J2, C_static, GI, ...


            }
        }
        ml.push_back(mi);
    }
    mlist.swap(ml);
    accel();
    analyse();
    return clist.size();
}
*/
void msystem::update_barycens(){
    std::vector<barycen> &bl=blist;
    int_t bn=bl.size();

    //update barycens: bl.rvGM,rvGM_sys;
    for(int_t i=0;i<bn;++i){//initialize
        bl[i].GM_sys=0;
        bl[i].r_sys=0;
        bl[i].v_sys=0;
    }
    for(int_t i=0;i<bn;++i)if(bl[i].mid>=0){//is mass
        barycen &bi=bl[i];
        const mass &mi=mlist[bi.mid];
        const real mGM=mi.GM;
        const mpvec rGM=mi.r*mGM;
        const mpvec vGM=mi.v*mGM;
        bi.r=mi.r;
        bi.v=mi.v;
        bi.GM=mGM;
        bi.r_sys+=rGM;
        bi.v_sys+=vGM;
        bi.GM_sys+=mGM;
        int_t bp=bi.pid;
        while(bp>=0){
            bl[bp].GM_sys+=mGM;
            bl[bp].r_sys+=rGM;
            bl[bp].v_sys+=vGM;
            bp=bl[bp].pid;
        }
    }
    for(int_t i=0;i<bn;++i)if(bl[i].mid<0){//is barycen
        barycen &bi=bl[i];
        barycen &h=bl[bi.hid];
        barycen &g=bl[bi.gid];
        bi.GM=h.GM_sys+g.GM_sys;
        bi.r=(h.r_sys+g.r_sys)/bi.GM;
        bi.v=(h.v_sys+g.v_sys)/bi.GM;
    }
    for(int_t i=0;i<bn;++i){//finalize
        barycen &bi=bl[i];
        bi.r_sys/=bi.GM_sys;
        bi.v_sys/=bi.GM_sys;
    }
}
