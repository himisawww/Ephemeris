#pragma once
#include<cstdint>
#include<string>
#include<vector>
#include<map>
#include<cstdio>

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
    // specified when open a file
    std::string filename;
    // used for memory read/write
    // [READ_CACHE]/WRITE_CACHE
    std::vector<byte_t> cached_data;
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

    MFILE_STATE state;

    static std::map<std::string,std::vector<byte_t>> mem_library;

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
    //resize cache and reset state to READ_CACHE,
    //return address of prepared cache.
    byte_t *reset(size_t new_cache_size);

    const std::string &get_name() const{ return filename; }
    void set_name(const std::string &_name){ filename=_name; }

    // (READ_FILE) read file to cache and converts to READ_CACHE
    void load_data();
    // (WRITE_CACHE/READ_CACHE) get data
    const byte_t *data() const{ return state==MFILE_STATE::READ_CACHE?idata:cached_data.data(); }
    const size_t size() const{ return state==MFILE_STATE::READ_CACHE?isize:cached_data.size(); }

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

    // (WRITE_CACHE)
    //  convert *this to READ_CACHE
    //  if(fname):
    //      put cache to a memory library with a file name;
    //      no writing to disk;
    //      can be further read by MFILE(filename);
    bool publish(const std::string &fname="");
};
