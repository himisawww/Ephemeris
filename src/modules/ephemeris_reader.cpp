#include"ephemeris_reader.h"
#include"utils/zipio.h"
#include"configs.h"
#include<string>

ephemeris_reader::ephemeris_reader(const char *ephemeris_path){
    bool is_broken=false;
    for(int dir=1;dir>=-1;dir-=2){
        const char *fwdbak=dir>0?"fwd":"bak";
        size_t cur_index=0;
        do{
            ++cur_index;
            izippack zp(strprintf("%s.%llu.%s.zip",ephemeris_path,cur_index,fwdbak));
            if(!zp)break;
            MFILE mf_index;
            auto mf_files=zp.load_central_directory();
            for(izipfile zf:mf_files){
                const std::string &zfn=zf.name();
                if(zfn==Configs::SaveNameIndex){
                    zf.fetch();
                    zf.dumpfile(mf_index);
                    break;
                }
            }
            printf("%s.%3d[%llu,%llu]\n",fwdbak,cur_index,mf_files.size(),mf_index.size());
        } while(1);
    }
}
