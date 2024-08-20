#pragma once
#include<cstdint>
#include<cstdarg>
#include<cstdio>
#include<string>
#include<vector>

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
int fprintf(MFILE *mem,const char *format,...);
int fseek(MFILE *_Stream,int64_t _Offset,int _Origin);
int64_t ftell(MFILE *_Stream);
int fclose(MFILE *_Stream);

std::string readline(MFILE *mem);
bool file_exist(const std::string &path);

class MFILE{
public:
    typedef uint8_t byte_t;
    typedef int64_t index_t;
private:
    // specified when open a file.
    // WRITE_FILE && !fp : cache will be published in this name when close;
    // [otherwise]: auxiliary;
    std::string filename;
    // local cache used for memory read/write
    // [READ_CACHE] : store data to be read;
    // WRITE_CACHE || WRITE_FILE && !fp : used to store wrote data;
    std::vector<byte_t> cached_data;
    // used for actuall disk io
    // READ_FILE || WRITE_FILE && fp : the file reading/writing;
    // WRITE_CACHE : if fp, cache will flush to this file when close;
    FILE *fp;
    // used for memory read
    // READ_CACHE : the begining address and byte size of reading cache
    const byte_t *idata;
    index_t isize;
    // used for memory read/write
    // READ_CACHE || WRITE_CACHE || WRITE_FILE &&!fp : current position of io
    index_t offset;

    /* usage of members under possible states:
        READ_CACHE:          [cached_data],   idata, isize, offset
         READ_FILE:                        fp
       WRITE_CACHE:           cached_data, fp,              offset
    WRITE_FILE:          
            &&  fp:                        fp
            && !fp: filename, cached_data,                  offset
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

    ~MFILE(){ close(); }

    int close();
    //resize local cache, set reading cache pointer to it, and set state to READ_CACHE.
    //return address of prepared local cache.
    byte_t *reset(size_t new_cache_size);

    const std::string &get_name() const{ return filename; }
    void set_name(const std::string &_name){ filename=_name; }

    // (READ_FILE/READ_CACHE) read file/memory to local cache and converts to READ_CACHE
    void load_data();
    // (WRITE_CACHE/READ_CACHE) get data
    const byte_t *data() const{ return state==MFILE_STATE::READ_CACHE?idata:cached_data.data(); }
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
    std::string readline();

    // (WRITE_CACHE/READ_CACHE)
    //  convert *this to READ_CACHE, load_data to local cache,
    //      and set offset to 0
    //  if(fname):
    //      put cache to a memory library with a file name;
    //      no writing to disk;
    //      can be further read by MFILE(filename);
    bool publish(const std::string &fname="");

    // When open failed in WRITE_FILE mode, that is, cannot create FILE for writing,
    // if true: MFILE is still a valid WRITE_FILE stream, all wrote data will be
    //          cached in memory and published in its filename when close.
    // if false: MFILE will be INVALID.
    static bool set_wcache_onfail(bool do_publish);
    static bool get_wcache_onfail();
};
