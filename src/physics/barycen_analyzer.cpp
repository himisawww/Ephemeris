#include<map>
#include"physics/mass.h"

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

real msystem::analyse(bool reconstruct){
    int_t old_bn=reconstruct?0:blist.size();
    int_t mn=mlist.size();
    if(old_bn&&t_barycen==t_eph)return t_update;
    t_barycen=t_eph;
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

            while(bp>=0&&!dl[bp].del){

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
                    fast_real rmax=std::max(dl[p.hid].minr,dl[p.gid].minr);
                    if(old_bn)rmax*=std::cbrt(.5);//make existing stucture more stable.
                    if(r1<rmax){//too close to be a companion of binary
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
                    fast_real rmax=std::max(dl[bi.pid].minr,dl[i].minr);
                    if(old_bn)rmax*=std::cbrt(.5);//make existing stucture more stable.
                    if(r1<rmax){//too close to be a child of binary
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

    if(old_bn){
        if(t_update!=t_update){
            for(int_t i=0;i<old_bn;++i){
                const barycen &bsrc=bl[i];
                barycen &bdst=blist[i];
                bdst.r=bsrc.r;
                bdst.v=bsrc.v;
                bdst.GM=bsrc.GM;
                bdst.r_sys=bsrc.r_sys;
                bdst.v_sys=bsrc.v_sys;
                bdst.GM_sys=bsrc.GM_sys;
            }
            t_update=t_eph;
        }
        return t_update;
    }

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
    barycen::fill_tid(bl);

    blist.swap(bl);
    t_update=t_eph;
    return t_update;
}

void barycen::fill_tid(std::vector<barycen> &bl){
    int_t bsize=bl.size();
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
}

void barycen::update_barycens(const msystem &ms,std::vector<barycen> &bl){
    int_t bn=bl.size();

    //update barycens: bl.rvGM,rvGM_sys;
    for(int_t i=0;i<bn;++i){//initialize
        bl[i].GM_sys=0;
        bl[i].r_sys=0;
        bl[i].v_sys=0;
    }
    for(int_t i=0;i<bn;++i)if(bl[i].hid<0){//is mass
        barycen &bi=bl[i];
        const mass &mi=ms[bi.mid];
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
    for(int_t i=0;i<bn;++i)if(bl[i].hid>=0){//is barycen
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

int_t barycen::decompose(std::vector<barycen> &blist,int_t bid){
    if(bid<0){
        int_t bn=blist.size();
        int_t rootid=-1;
        int_t nroot=0;
        for(int_t i=0;i<bn;++i){
            if(blist[i].pid<0){
                rootid=i;
                ++nroot;
            }
        }
        return nroot==1?decompose(blist,rootid):-1;
    }

    barycen &b=blist[bid];
    int_t nret=0;
    for(const auto cid:b.children)
        nret+=decompose(blist,cid);

    if(b.gid>=0){
        nret+=decompose(blist,b.gid);
        nret+=decompose(blist,b.hid);
    }

    if(b.pid>=0){
        barycen &p=blist[b.pid];
        if(bid==p.hid){
            b.r=NAN;
            b.v=NAN;
        }
        else if(bid==p.gid){
            b.r=b.r_sys-blist[p.hid].r_sys;
            b.v=b.v_sys-blist[p.hid].v_sys;
        }
        else{
            b.r=b.r_sys-p.r;
            b.v=b.v_sys-p.v;
        }
    }
    else{
        b.r=b.r_sys;
        b.v=b.v_sys;
    }

    b.r_sys=NAN;
    b.v_sys=NAN;

    return nret+1;
}

int_t barycen::compose(std::vector<barycen> &blist,int_t bid){
    if(bid<0){
        int_t bn=blist.size();
        int_t rootid=-1;
        int_t nroot=0;
        for(int_t i=0;i<bn;++i){
            if(blist[i].pid<0){
                rootid=i;
                ++nroot;
            }
        }
        return nroot==1?compose(blist,rootid):-1;
    }

    barycen &b=blist[bid];
    int_t nret=1;

    if(b.pid<0){
        b.r_sys=b.r;
        b.v_sys=b.v;
    }
    else{
        barycen &p=blist[b.pid];
        if(bid==p.gid){
            b.r_sys=b.r+blist[p.hid].r_sys;
            b.v_sys=b.v+blist[p.hid].v_sys;
        }
        else if(bid==p.hid){
            barycen &g=blist[p.gid];
            real gdm=g.GM_sys/p.GM;
            b.r_sys=p.r-gdm*g.r;
            b.v_sys=p.v-gdm*g.v;
        }
        else{
            b.r_sys=b.r+p.r;
            b.v_sys=b.v+p.v;
        }
    }

    mpvec cracc(0),cvacc(0);
    for(const auto cid:b.children){
        barycen &c=blist[cid];
        cracc+=c.r*c.GM_sys;
        cvacc+=c.v*c.GM_sys;
    }
    b.r=b.r_sys-cracc/b.GM_sys;
    b.v=b.v_sys-cvacc/b.GM_sys;

    if(b.gid>=0){
        nret+=compose(blist,b.hid);
        nret+=compose(blist,b.gid);
    }

    for(const auto cid:b.children)
        nret+=compose(blist,cid);

    return nret;
}
