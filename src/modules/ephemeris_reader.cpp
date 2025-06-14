#include"ephemeris_reader.h"
#include<algorithm>
#include"ephemeris_generator.h"
#include"configs.h"
#include"utils/logger.h"

ephemeris_reader::chapter::chapter(msystem &ms,const std::string &_):izippack(_){
    std::string failure;
    do{
        if(!izippack::operator bool())
            break;
        auto mf_files=load_central_directory();
        MFILE mf_index;
        bool index_loaded=false,chpt_loaded=false;
        for(izipfile zf:mf_files){
            const std::string &zfn=zf.name();
            if(zfn==Configs::SaveNameIndex){
                if(index_loaded)failure+="    Duplicated index;\n";
                else{
                    index_loaded=true;
                    zf.fetch();
                    zf.dumpfile(mf_index);
                }
            }
            else if(zfn==Configs::SaveNameCheckpoint){
                if(chpt_loaded)failure+="    Duplicated checkpoint;\n";
                else{
                    chpt_loaded=true;
                    MFILE mf_chpt;
                    zf.fetch();
                    zf.dumpfile(mf_chpt);
                    if(ms.empty()){
                        if(!ms.load_checkpoint(&mf_chpt))
                            failure+="    Invalid checkpoint;\n";
                    }
                    else{
                        msystem msc;
                        msc.load_checkpoint(&mf_chpt);
                        if(!msc.is_same(ms))
                            failure+="    Incompatible celestial system between ephemerides;\n";
                    }
                }
            }
            if(failure.size()||index_loaded&&chpt_loaded)break;
        }
        if(failure.size()||!(index_loaded&&chpt_loaded)){
            if(failure.empty())failure+="    Missing index/checkpoint;\n";
            break;
        }

        std::map<std::string,izipfile> allfiles;
        for(const auto &zf:mf_files){
            if(!allfiles.try_emplace(get_file_name(zf.name()),zf).second){
                failure+="    Duplicated file;\n";
                break;
            }
        }
        if(failure.size())break;

        dir=0;
        size_t mn=ms.size();
        ephm_index.resize(mn);
        while(failure.empty()){
            ephemeris_entry index;
            if(1!=fread(&index,sizeof(index),1,&mf_index))
                break;
            bool fwd=index.t_start<index.t_end;
            bool bak=index.t_start>index.t_end;
            int _dir=fwd?1:-1;
            if(fwd==bak||dir&&dir!=_dir){
                failure+="    Incompatible direction;\n";
                break;
            }
            int_t cur_min=dir>0?index.t_start:index.t_end;
            int_t cur_max=dir>0?index.t_end:index.t_start;
            if(dir){
                t_start=std::min(t_start,cur_min);
                t_end=std::max(t_end,cur_max);
            }
            else{
                dir=_dir;
                t_start=cur_min;
                t_end=cur_max;
            }
            if(index.fid==0){
                index.fid=blists.size();
                auto it=blist_index.try_emplace(dir*index.t_end,index);
                if(!it.second||!barycen::load_barycen_structure(blists.emplace_back(),&mf_index,index.sid))
                    failure+="    Error loading system structure;\n";
            }
            else{
                auto itorb=allfiles.find(index.entry_name(false,false));
                auto itrot=allfiles.find(index.entry_name(true,false));
                if(itorb==allfiles.end()||itrot==allfiles.end()){
                    failure+="    Missing file;\n";
                    break;
                }
                size_t mid=ms.get_mid(index.sid);
                if(mid>=mn){
                    failure+="    Invalid sid;\n";
                    break;
                }
                index.fid=ephm_files.size();
                if(!ephm_index[mid].try_emplace(dir*index.t_end,index).second){
                    failure+="    Duplicated entry;\n";
                    break;
                }
                ephm_files.push_back(itorb->second);
                ephm_files.push_back(itrot->second);
            }
        }
        if(failure.size()||!dir){
            if(failure.empty())failure+="    Missing direction;\n";
            break;
        }

        std::set<izipfile> dedup_file(ephm_files.begin(),ephm_files.end());
        dedup_file.insert(izippack::end());
        if(dedup_file.size()!=ephm_files.size()+1){
            failure+="    Mismatch ephemeris;\n";
            break;
        }

        if(dir<0)
            std::swap(t_start,t_end);

        for(size_t i=0;i<=mn;++i){
            const auto &range=i==mn?blist_index:ephm_index[i];
            int_t t_min,t_max,n_count=0;
            for(const auto &p:range){
                const auto &eid=p.second;
                if(n_count==0)
                    t_min=eid.t_start;
                else if(t_max!=eid.t_start)
                    break;
                t_max=eid.t_end;
                ++n_count;
            }
            if(!n_count||n_count!=range.size()||t_min!=t_start||t_max!=t_end){
                failure+="    Invalid coverage;\n";
                break;
            }
        }
        if(failure.size())break;

        return;
    } while(0);

    if(failure.size())
        LogError("Error loading %s:\n%s",_.c_str(),failure.c_str());
    close();
}

ephemeris_reader::ephemeris_reader(const char *ephemeris_path){
    int_t middle_fwd,middle_bak;
    int_t t_fwd,t_bak;
    int_t n_fwd=0,n_bak=0;
    for(int dir=1;dir>=-1;dir-=2){
        bool fwd=dir>0;
        const char *fwdbak=fwd?"fwd":"bak";
        size_t cur_index=0;
        int_t &t_mid=fwd?middle_fwd:middle_bak;
        int_t &t_last=fwd?t_fwd:t_bak;
        int_t &n_chpt=fwd?n_fwd:n_bak;
        std::string failure;
        do{
            std::string chname=strprintf("%s.%llu.%s.zip",ephemeris_path,++cur_index,fwdbak);
            chapter curchpt(ms,chname);
            if(!curchpt)
                break;
            if(curchpt.dir!=dir||(dir>0?curchpt.t_start>=curchpt.t_end:curchpt.t_start<=curchpt.t_end)){
                failure+="    Incompatible direction;\n";
                break;
            }
            if(cur_index==1){
                t_mid=curchpt.t_start;
                if(!fwd&&n_fwd&&t_mid!=middle_fwd){
                    failure+="    Mismatch epoch;\n";
                    break;
                }
            }
            else if(curchpt.t_start!=t_last){
                failure+="    Uncontinuous ephemeris;\n";
                break;
            }
            t_last=curchpt.t_end;
            n_chpt+=1;
            chapters.emplace_back(std::move(curchpt));
        } while(1);
        if(failure.size())
            LogWarning("Warning: %s ephemerides of %s is interrupted at %llu.%s.zip:\n%s",
                fwd?"Forward":"Backward",ephemeris_path,cur_index,fwdbak,failure.c_str());
        std::vector<chapter> revchpts;
        for(auto it=chapters.rbegin();it!=chapters.rend();++it)
            revchpts.emplace_back(std::move(*it));
        revchpts.swap(chapters);
    }
    if(!n_fwd&&!n_bak){
        ms.clear();
        return;
    }
    if(!n_fwd)t_fwd=middle_bak;
    if(!n_bak)t_bak=middle_fwd;
    LogInfo(
        "Loaded %llu bodies from ephemeris %s\n"
        "    Time Span: [%lld, %lld]\n",
        ms.size(),ephemeris_path,t_bak,t_fwd);
}

