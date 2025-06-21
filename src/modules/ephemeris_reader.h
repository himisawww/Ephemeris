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
        std::vector<ephemeris_interpolator> ephm_interps;
        std::vector<bool> ephm_expand;

        //if dir<0, t_end < t_start
        int_t t_start,t_end;
        int_t _interp_size,_cache_bytes;

        friend class ephemeris_reader;
    public:
        int_t t_min() const{ return std::min(t_start,t_end); }
        int_t t_max() const{ return std::max(t_end,t_start); }
        int_t interpolator_size() const{ return _interp_size; }
        //free interpolator cache
        void unload();

        chapter(msystem &,const std::string &,int_t memory_budget);
        bool checkout(msystem &,real t_eph);
    };
    std::vector<chapter> chapters;
    std::vector<std::string> massnames;
    int_t cur_chid;
    int_t mem_budget;
    std::vector<int_t> active_chapters;

    //assume chapters[cur_chid valid]
    void lru();
public:
    //memory_budget: a suggestion(not mandatory) of memory limit in bytes for this reader
    ephemeris_reader(const char *ephemeris_path,int_t memory_budget=768*int_t(1024*1024));

    operator bool() const{ return !ms.empty(); }
    int_t t_min() const{ return chapters.empty()?0:chapters.front().t_min(); }
    int_t t_max() const{ return chapters.empty()?0:chapters.back().t_max(); }
    int_t interpolator_size() const;
    //free interpolator cache
    void unload();

    const msystem &get_msystem(){ return ms; }
    bool checkout(real t_eph);
};
