#include"memio.h"

/*
    Simple .zip IO utilities,
        support Store only, no (de)compression.
        support utf8, zip64 extension.

    open & traverse a .zip file:

        izippack zp(path_to_zip);
        for(const izipfile &zf:zp){
            //  izipfile zf contains file information
            //  and can be dumped to memory by zf.dumpfile(&memory)
        }

    create and save a .zip file:

        zipmem zm;
        // zipmem zm contains name and content of a file
        ozippack zp(path_to_zip, true);
        zp.zipmems.push_back(zm...)
        // .zip file will be written on destruction of ozippack zp.

*/

class izippack{
    MFILE *fzip;
    friend class izipfile;
public:
    operator bool() const{ return fzip; }

    izippack(izippack &&_i);
    //name of the zip file
    explicit izippack(const std::string &filename);
    ~izippack();
    izipfile begin() const;
    izipfile end() const;
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
    ~ozippack();
    void push_back(std::string fullname,MFILE &&mf);
    void push_back(const izipfile &zf);
    void push_back(const izippack &izip);
};

class izipfile{
    const izippack *pzip;
    size_t offset,fileoffset,filesize;
    std::string filename;
    friend class izippack;
public:
    bool operator==(const izipfile &ft2) const{
        return pzip==ft2.pzip&&offset==ft2.offset;
    }
    bool operator!=(const izipfile &ft2) const{
        return pzip!=ft2.pzip||offset!=ft2.offset;
    }
    void validate();
    izipfile &operator*(){
        return *this;
    }
    izipfile *operator->(){
        return this;
    }

    izipfile &operator++();
    izipfile operator++(int);
    const std::string &name() const{ return filename; }
    std::string fullname() const;
    size_t size() const{ return filesize; }
    void dumpfile(MFILE &mf) const;
};
