#pragma once
#include"memio.h"

/*
    Simple .zip IO utilities,
        support Store only, no (de)compression/en(/de)cryption/multi-voluming.
        support utf8, zip64 extension.

    open & traverse a .zip file:

        izippack zp(path_to_zip);

        // first method:
        for(const izipfile &zf:zp){
            //  izipfile zf contains file information
            //  and can be dumped to MFILE memory by zf.dumpfile(memory)
        }

        //second method:
        for(izipfile zf:zp.load_central_directory()){
            //  avoid traversing full izippack by looking up central directory
            //  in this case, zf.fetch() is required before zf.dumpfile(memory)
        }

    create and save a .zip file:

        MFILE zm;
        // zm contains content of a file
        ozippack zp(path_to_zip);
        zp.push_back("filename",zm)...
        // .zip file will be written on destruction of ozippack zp.

*/

class izippack;
class izipfile{
    const izippack *pzip;
    size_t locoffset,fileoffset,filesize;
    std::string filename;

    friend class izippack;
    uint16_t load_zip64(uint16_t exremain,bool is_central);
    //load from local header at _offset
    izipfile(const izippack *,size_t _offset);
    //load from central header at ftell
    izipfile(const izippack *);
public:
    static constexpr size_t npos=-1;
    operator bool() const{ return locoffset!=npos; }
    //ready for dumpfile(), operator++(), or offset()
    bool is_ready() const{ return fileoffset!=npos; }
    //make izipfile is_ready(), return success
    bool fetch();
    //get the next file by scanning izippack from beginning
    //requires is_ready()
    izipfile &operator++();
    izipfile operator++(int);
    //return success
    //requires is_ready()
    bool dumpfile(MFILE &mf) const;

    bool operator==(const izipfile &ft2) const{
        return pzip==ft2.pzip&&locoffset==ft2.locoffset;
    }
    bool operator!=(const izipfile &ft2) const{
        return pzip!=ft2.pzip||locoffset!=ft2.locoffset;
    }
    bool operator<(const izipfile &ft2) const{
        return locoffset<ft2.locoffset;
    }
    izipfile &operator*(){
        return *this;
    }
    izipfile *operator->(){
        return this;
    }

    const std::string &name() const{ return filename; }
    std::string fullname() const;
    size_t size() const{ return filesize; }
    //requires is_ready()
    size_t offset() const{ return fileoffset; }
};

class izippack{
    MFILE *fzip;
    friend class izipfile;
public:
    operator bool() const{ return fzip; }
    //izipfile can be located by offset() and size() on this MFILE
    //for read-only access, resource managed by izippack, do not close()
    MFILE *get_file() const{ return fzip; }

    izippack(izippack &&_i);
    //name of the zip file
    explicit izippack(const std::string &filename);
    int close();
    ~izippack(){ close(); }
    //scan from beginning of zippack, is_ready() is ensured
    izipfile begin() const;
    izipfile end() const;
    //load central directory at end of zippack, avoid full scan
    //in this case fetch() is required to make izipfile is_ready()
    std::vector<izipfile> load_central_directory() const;
};

class ozippack:private std::vector<MFILE>{
    MFILE *fzip;
    typedef std::vector<MFILE> base_t;
public:
    using base_t::begin;
    using base_t::end;
    using base_t::size;
    using base_t::operator[];
    using base_t::swap;
    using base_t::resize;
    using base_t::reserve;

    operator bool() const{ return fzip; }
    ozippack(ozippack &&_o);
    //name of the zip file
    explicit ozippack(const std::string &filename);
    int close();
    ~ozippack(){ close(); }
    void push_back(std::string fullname,MFILE &&mf);
    void push_back(const izipfile &zf);
    void push_back(const izippack &izip);
};
