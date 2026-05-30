#pragma once
#include<cstdint>
#include<cstdarg>
#include<cstdio>
#include<string>
#include"htl/vector.h"

/* usage:
    
    1. Open file for read/write
    2. Open memory for read/write
    3. Convert memory-write to memory-read
    4. Write content to file

*/
class MFILE;

enum class MFILE_STATE{
    INVALID,
    READ_FILE,READ_CACHE,
    WRITE_FILE,WRITE_CACHE
};

// open file
// fname should be utf-8 encoded
FILE *fopen(const std::string &fname,bool is_read);

// open memory to write (WRITE_CACHE)
MFILE *mopen();
// open memory to read (READ_CACHE)
MFILE *mopen(const void *_mem,size_t _size);
// open file to read/write (read?READ:WRITE)_(cache?CACHE:FILE)
MFILE *mopen(const std::string &_fname,MFILE_STATE _state=MFILE_STATE::READ_FILE);

size_t fread(void *buffer,size_t e_size,size_t e_count,MFILE *mem);
size_t fwrite(const void *buffer,size_t e_size,size_t e_count,MFILE *mem);
std::string vstrprintf(const char *format,va_list _arg_list);
std::string strprintf(const char *format,...);
int64_t fprintf(MFILE *mem,const char *format,...);
int fseek(MFILE *_stream,int64_t _offset,int _origin);
int64_t ftell(MFILE *_stream);
int fclose(MFILE *_stream);

// if _enable_skip, skip lines empty or starts with #
std::string fgetstr(MFILE *mem,bool _enable_skip=false);
bool file_exist(const std::string &path);
std::string get_file_name(const std::string &path);
std::string get_file_extension(const std::string &path);

class MFILE{
    // specified when open a file.
    // WRITE_FILE && !fp : cache will be published in this name when close;
    // [otherwise]: auxiliary;
    std::string filename;
    // local cache used for memory read/write
    // [READ_CACHE] : store data to be read;
    // WRITE_CACHE || WRITE_FILE && !fp : used to store wrote data;
    htl::vector<char> cached_data;
    // used for actuall disk io
    // READ_FILE || WRITE_FILE && fp : the file reading/writing;
    // WRITE_CACHE : if fp, cache will flush to this file when close;
    FILE *fp;
    // used for memory read
    // READ_CACHE : the begining address and byte size of reading cache
    const char *idata;
    int64_t isize;
    // used for memory read/write
    // READ_CACHE || WRITE_CACHE || WRITE_FILE &&!fp : current position of io
    int64_t offset;
    // if fp valid && fp, fp is owned by MFILE and will be closed at destruction
    bool own;

    /* usage of members under possible states:
        READ_CACHE:          [cached_data],         idata, isize, offset
         READ_FILE:                        *fp, own
       WRITE_CACHE:           cached_data,  fp, own,              offset
    WRITE_FILE:
            &&  fp:                        *fp, own
            && !fp: filename, cached_data,                        offset
    */
    MFILE_STATE state;

public:
    MFILE(MFILE &&mf);
    // open memory to write (WRITE_CACHE)
    MFILE();
    // open memory to read (READ_CACHE)
    MFILE(const void *_mem,size_t _size);
    // open file to read/write (read?READ:WRITE)_(cache?CACHE:FILE)
    MFILE(const std::string &_fname,MFILE_STATE _state=MFILE_STATE::READ_FILE);
    // FILE * wrapper. caller should specify correct READ_FILE/WRITE_FILE _state. 
    // if take_own, fclose(_fp) is called when MFILE destructs/closes.
    MFILE(FILE *_fp,MFILE_STATE _state=MFILE_STATE::READ_FILE,bool take_own=false);

    ~MFILE(){ close(); }

    int close();
    //resize local cache, set reading cache pointer to it, and set state to READ_CACHE.
    //return address of prepared local cache.
    char *prepare(size_t new_cache_size);
    //empty local cache, and set state to WRITE_CACHE
    void reset();

    const std::string &get_name() const{ return filename; }
    void set_name(const std::string &_name){ filename=_name; }

    // (READ_FILE/READ_CACHE) read file/memory to local cache and converts to READ_CACHE
    void load_data();
    // (WRITE_CACHE/READ_CACHE) get data
    const char *data() const{ return state==MFILE_STATE::READ_CACHE?idata:cached_data.data(); }
    size_t size() const{ return state==MFILE_STATE::READ_CACHE?isize:cached_data.size(); }

    bool is_valid() const;
    bool is_read() const;
    bool is_write() const;
    bool is_cache() const;

    void reserve(size_t rsize);

    int64_t tell() const;
    int seek(int64_t fpos,int forg);

    size_t read(void *buffer,size_t e_size,size_t e_count);
    size_t write(const void *buffer,size_t e_size,size_t e_count);
    // if _enable_skip, skip lines empty or starts with #
    std::string fgetstr(bool _enable_skip=false);

    // (WRITE_CACHE/READ_CACHE/READ_FILE)
    //  convert *this to READ_CACHE, load_data to local cache,
    //      and set offset to 0
    //  if(filename):
    //      put cache to a memory library with a file name;
    //      no writing to disk;
    //      can be further read by MFILE(filename);
    bool publish();
    bool publish(const std::string &filename);

    // When open failed in WRITE_FILE/WRITE_CACHE mode, that is, cannot create FILE for writing,
    // if true: MFILE is still a valid WRITE_FILE stream, all wrote data will be
    //          cached in memory and published in its filename when close.
    // if false: MFILE will be INVALID.
    static bool set_wcache_onfail(bool do_publish);
    static bool get_wcache_onfail();
};
