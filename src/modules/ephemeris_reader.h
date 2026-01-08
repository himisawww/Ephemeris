#pragma once
#include"math/keplerian.h"
#include"physics/mass.h"
#include"utils/zipio.h"
#include"modules/ephemeris_compressor.h"

class ephemeris_reader{
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
        bool checkout(ephemeris_reader &,real t_eph);
    };
    struct massinfo{
        std::string name;
        //followings are about/relative to orbital center:
        //and is updated by checkout only-if update_orbits is true
        keplerian parameters;
        orbital_state_t state_vectors;
        double keplerian_GM;
    };
    //data and states
    msystem ms;
    std::vector<chapter> chapters;
    std::vector<int_t> active_chapters;
    std::vector<massinfo> minfos;
    int_t cur_chid;
    //configs
    int_t mem_budget;
public:
    //if update_physics, passed to accel()
    //default 0
    int update_physics_parallel_option;
    //if true, checkout() also calls accel() for get_msystem();
    //default: false
    bool update_physics;
    //if true, checkout() also updates get_msystem().get_barycens();
    //default: false
    bool update_bsystem;
    //if true, checkout() also updates massinfo::GM,parameter,r,v;
    //default: false
    bool update_orbits;
public:
    //memory_budget: a suggestion(not mandatory) of memory limit in bytes for this reader
    ephemeris_reader(const char *ephemeris_path,int_t memory_budget=768*int_t(1024*1024));

    explicit operator bool() const{ return !ms.empty(); }
    int_t t_min() const{ return chapters.empty()?0:chapters.front().t_min(); }
    int_t t_max() const{ return chapters.empty()?0:chapters.back().t_max(); }
    int_t interpolator_size() const;
    //free interpolator cache
    void unload();

    const msystem &get_msystem() const{ return ms; }
    const massinfo &get_massinfo(int_t mid) const{ return minfos[mid]; }
    const massinfo &get_massinfo(const char *ssid) const{ return minfos[ms.get_mid(ssid)]; }
    bool checkout(real t_eph);
private:
    //assume chapters[cur_chid valid]
    void lru();
};
