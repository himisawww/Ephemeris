#pragma once
#include"math/keplerian.h"
#include"physics/mass.h"
#include"utils/zipio.h"
#include"modules/ephemeris_compressor.h"

class ephemeris_reader{
public:
    enum select_type:uint8_t{
        NONE     = 0,
        ORBIT    = 1,
        ROTATION = 2,
        BOTH     = 3 //however, rotation requires orbit
    };
    class massinfo{
        const mass *const mref;
        uint8_t config;
        friend class ephemeris_reader;
    public:
        std::string name;
        //followings are about/relative to orbital center:
        //and is updated by checkout only-if update_orbits is true
        keplerian parameters;
        orbital_state_t state_vectors;
        double keplerian_GM;

        massinfo(const mass &_mref):mref(&_mref),config(NONE){}
        const mass &operator *() const{ return *mref; }
        const mass *operator->() const{ return  mref; }
    };
private:
    enum active_type:uint8_t{
        REQUIRE_RVSYSC  = 1,
        REQUIRE_RVC     = 2,
        REQUIRE_RVD     = 4,
        REQUIRE_OFFSET  = 8,
        LOADED_RVD      = 16,
        LOADED_OFFSET   = 32,
        LOADED_ROTATION = 64
    };
    struct active_info{
        active_type type;
        int_t bid,fid,t_start,t_end;
        active_info(active_type,int_t,ephemeris_entry);
    };
    struct GM_info{
        real GM0,GM0_sys,dGM,dGM_sys;
        real GM    (fast_real t) const{ return GM0    +t*dGM    ; }
        real GM_sys(fast_real t) const{ return GM0_sys+t*dGM_sys; }
    };
    class chapter:public izippack{
        //file name of chapter
        std::string chname;
        //{dir*t_end, blist index over [t_start,t_end]}
        htl::map<int_t,ephemeris_entry> blist_index;
        //blists[blist index.fid] = system structure
        htl::vector<bsystem> blists;

        //[mid]={dir*t_end, ephemeris index over [t_start,t_end]}
        htl::vector<htl::map<int_t,ephemeris_entry>> ephm_index;
        //[blist index.fid]={bid, offset index}
        htl::vector<htl::map<int_t,ephemeris_entry>> offset_maps;
        //[ephemeris index.fid+(0/1)] = {orbital,rotational} ephemerides
        //[   offset index.fid      ] = {offset            } ephemerides
        htl::vector<izipfile> ephm_files;
        //[same as ephm_files]
        htl::vector<ephemeris_interpolator> ephm_interps;

        //if dir<0, t_end < t_start
        int_t t_start,t_end;
        //preparation for partial checkout
        int_t _bselect,_brootid,_interp_size;
        fast_real _ft_eph;
        htl::vector<int_t> active_composed;
        htl::vector<active_info> active_files;
        htl::vector<uint8_t> active_map;
        htl::vector<GM_info> GM_map;

        friend class ephemeris_reader;
        int_t compose_active(int_t bid);
    public:
        int_t t_min() const{ return std::min(t_start,t_end); }
        int_t t_max() const{ return std::max(t_end,t_start); }
        int_t interpolator_size() const{ return _interp_size; }
        //free interpolator cache
        void unload();

        chapter(msystem &,const std::string &);
        bool checkout(ephemeris_reader &,real t_eph);
    };
    //data and states
    msystem ms;
    htl::vector<chapter> chapters;
    htl::vector<int_t> active_chapters;
    htl::vector<massinfo> minfos;
    int_t cur_chid;
    //configs
    int_t memory_limit;
    bool update_orbits;
    bool update_bsystem;
    bool update_physics;
    int update_physics_parallel_option;
public:
    ephemeris_reader(const char *ephemeris_path);

    explicit operator bool() const{ return !ms.empty(); }
    int_t t_min() const{ return chapters.empty()?0:chapters.front().t_min(); }
    int_t t_max() const{ return chapters.empty()?0:chapters.back().t_max(); }
    int_t interpolator_size() const;
    //free interpolator cache
    void unload();

    size_t size() const{ return ms.size(); }
    bool empty() const{ return ms.empty(); }
    auto begin() const{ return minfos.begin(); }
    auto end() const{ return minfos.end(); }
    const massinfo &operator[](int_t mid) const{ return minfos[mid]; }
    const massinfo &operator[](const char *ssid) const{ return minfos[ms.get_mid(ssid)]; }

//vvvv configurators vvvv
    //memory_limit: a suggestion(not mandatory) of memory limit in bytes for this reader
    //default: see parameter
    int_t set_memory_limit(int_t new_bytes=768*int_t(1024*1024)){ return std::exchange(memory_limit,new_bytes); }
    //if true, checkout() also updates massinfo::GM,parameter,r,v;
    //default: false
    bool set_update_orbits(bool new_setting){ return std::exchange(update_orbits,new_setting); }
    //if true, bsystem of internal msystem will be filled, including structure & all GM/GM_sys,
    // but not rv/rv_sys of unselected objects.
    //default: false
    bool set_update_bsystem(bool new_setting){ return std::exchange(update_bsystem,new_setting); }
    //note: performance warning;
    //if true, implies set_update_bsystem() and select_all();
    // all fields of internal msystem & bsystem will be filled, ready for further integrating.
    //default: false
    bool set_update_physics(bool new_setting);
    //if update_physics, passed to accel()
    //default 0
    int set_physics_parallel_option(int new_option){ return std::exchange(update_physics_parallel_option,new_option); }

    //only checkout some part of system ephemerides, as an optimization
    //returns config before select/deselect
    select_type select(const massinfo &,select_type _s=BOTH);
    select_type deselect(const massinfo &,select_type _s=BOTH);
    const massinfo &select(int_t mid,select_type _s=BOTH){
        const massinfo &minfo=minfos[mid];
        select(minfo,_s);
        return minfo;
    }
    const massinfo &select(const char *ssid,select_type _s=BOTH){
        const massinfo &minfo=minfos[ms.get_mid(ssid)];
        select(minfo,_s);
        return minfo;
    }
    //note: performance warning
    size_t select_all(select_type _s=BOTH);
    size_t deselect_all(select_type _s=BOTH);
//^^^^ configurators ^^^^

    //must select something to checkout before call this
    bool checkout(real t_eph);
private:

    void reset_selection();
    //assume chapters[cur_chid valid]
    void lru();
};
