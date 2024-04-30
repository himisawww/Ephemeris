#pragma once
#include<string>
#include<vector>

/*
    Simple .zip IO utilities,
        support Store only, no (de)compression.
        support utf8, zip64 extension.
     
    open & traverse a .zip file:

        zippack zp(path_to_zip);
        for(const zipfile &zf:zp){
            //  zipfile zf contains file information
            //  and can be dumped to memory by zf.dumpfile(&memory)
        }

    create and save a .zip file:
        
        zipmem zm;
        // zipmem zm contains name and content of a file
        zippack zp(path_to_zip, true);
        zp.zipmems.push_back(zm...)
        // .zip file will be written on destruction of zippack zp.
    
*/
class zippack;

class zipfile{
public:
    zippack *pzip;
    size_t offset,fileoffset,filesize;
    std::wstring filename;
public:
    bool operator==(const zipfile &ft2) const{
        return pzip==ft2.pzip&&offset==ft2.offset;
    }
    bool operator!=(const zipfile &ft2) const{
        return pzip!=ft2.pzip||offset!=ft2.offset;
    }
    void validate();
    zipfile &operator*(){
        return *this;
    }
    zipfile *operator->(){
        return this;
    }

    zipfile &operator++();
    zipfile operator++(int);
    std::wstring fullname() const;
    void dumpfile(void *dst) const;
};

//file content to be written to zip
class zipmem{
public:
    std::wstring filename;
    std::vector<std::uint8_t> data;
};

class zippack{
public:
    std::wstring zippath;
    FILE *fzip;
    std::vector<zipmem> zipmems;
    bool is_construct;
public:
    //name of the zip file
    zippack(const std::wstring &filename,bool construct=false);
    ~zippack();
    zipfile begin();
    zipfile end();
};

