#pragma once
#include<vector>
#include<cstdint>
#include<cstdio>
#include<string>
#include"ziptree.h"

/* usage:
    
    1. Open file for read
    2. Open memory for read
    3. Open memory for write
    4. Convert memory-write to memory-read
    5. Write memory to file

*/
class mem_file{
    typedef uint8_t byte_t;
    typedef int64_t index_t;
    enum class mem_file_state{
        INVALID,
        READ_FILE,READ_CACHE,
        WRITE_FILE,WRITE_CACHE
    };

    // specified when open a file
    std::string fstr;
    // used for memory read/write
    // [READ_CACHE]/WRITE_CACHE
    std::vector<byte_t> cache;
    // used for file read/write
    // when in WRITE_CACHE state, cache will flush to this file when close
    // READ_FILE/WRITE_FILE/WRITE_CACHE
    FILE *fp;
    // used for memory read
    // READ_CACHE
    const byte_t *idata;
    index_t isize;
    // used for memory read/write
    // WRITE_CACHE/READ_CACHE
    index_t offset;

    mem_file_state state;

public:
    // open memory to write (WRITE_CACHE)
    mem_file();
    // open memory to read (READ_CACHE)
    mem_file(const void *_mem,size_t _size);
    // open file to read/write (read?READ:WRITE)_(cache?CACHE:FILE)
    mem_file(const char *fname,bool is_read=true,bool is_cache=false);

    // (READ_CACHE/WRITE_CACHE) save cache to file
    size_t save(const char *fname) const;
    // (READ_CACHE/WRITE_CACHE) swap cache with e_data
    void swap(std::vector<byte_t> &e_data);

    void close();

    // (READ_FILE) read file to cache
    void cache_read();
    // (WRITE_CACHE/READ_CACHE/READ_FILE) get data
    const std::vector<byte_t> &data();

    const char *c_str() const;

    // =is_valid()
    operator bool() const;
    bool is_valid() const;
    bool is_read() const;
    bool is_write() const;
    bool is_cache() const;

    void reserve(size_t rsize);

    int64_t tell() const;
    int seek(int64_t fpos,int forg);

    void read(const zipfile &zf);
    size_t read(void *buffer,size_t e_size,size_t e_count);
    size_t write(const void *buffer,size_t e_size,size_t e_count);
    std::string readline();

    // (WRITE_CACHE)
    //  if(fname):
    //      put cache to a memory library with a file name;
    //      no writing to disk;
    //      can be further read by mem_file(filename);
    //      *this will be invalidated;
    //  else :
    //      convert *this to READ_CACHE
    bool publish(const char *fname=nullptr);
};
size_t fread(void *buffer,size_t e_size,size_t e_count,mem_file *mem);
size_t fwrite(const void *buffer,size_t e_size,size_t e_count,mem_file *mem);
std::string readline(mem_file &mem);
int fprintf(mem_file *mem,const char *format,...);
