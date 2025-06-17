#pragma once
#include"physics/mass.h"
#include"utils/zipio.h"
#include"modules/ephemeris_compressor.h"

class ephemeris_reader{
    msystem ms;
    class chapter:public izippack{
        //{dir*t_end, blist index over [t_start,t_end]}
        std::map<int_t,ephemeris_entry> blist_index;
        //blists[blist index.fid] = system structure
        std::vector<bsystem> blists;

        //[mid]={dir*t_end, ephemeris index over [t_start,t_end]}
        std::vector<std::map<int_t,ephemeris_entry>> ephm_index;
        //ephm_files[ephemeris index.fid+(0/1)] = {orbital,rotational} ephemerides
        std::vector<izipfile> ephm_files;

        int dir;
        //if dir<0, t_end < t_start
        int_t t_start,t_end;

        friend class ephemeris_reader;
    public:
        int_t t_min() const{ return dir>0?t_start:t_end; }
        int_t t_max() const{ return dir>0?t_end:t_start; }
        chapter(msystem &,const std::string &);
    };
    std::vector<chapter> chapters;
    std::vector<std::string> massnames;
public:
    ephemeris_reader(const char *ephemeris_path);

    operator bool() const{ return !ms.empty(); }
};
