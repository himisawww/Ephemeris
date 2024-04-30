#include"memfile.h"
#include<cstdarg>
#include<unordered_map>

static std::unordered_map<std::string,mem_file *> mem_library;

mem_file::mem_file(){
    fp=nullptr;
    offset=0;
    state=mem_file_state::WRITE_CACHE;
}
mem_file::mem_file(const void *_mem,size_t _size){
    idata=(const byte_t *)_mem;
    isize=_size;
    offset=0;
    state=mem_file_state::READ_CACHE;
}
mem_file::mem_file(const char *fname,bool is_read,bool is_cache){
    fstr=fname;
    state=mem_file_state::INVALID;

    if(is_read){
        auto it=mem_library.find(fstr);
        if(it!=mem_library.end()){
            mem_file *pmem=it->second;
            idata=pmem->cache.data();
            isize=pmem->cache.size();
            offset=0;
            state=mem_file_state::READ_CACHE;
            return;
        }
    }

    fp=fopen(fname,is_read?"rb":"wb");
    if(!fp)return;
    
    if(is_read){
        state=mem_file_state::READ_FILE;
        if(is_cache)cache_read();
    }
    else{
        state=is_cache?mem_file_state::WRITE_CACHE:mem_file_state::WRITE_FILE;
        offset=0;
    }
}

size_t mem_file::save(const char *fname) const{
    if(!is_cache())return 0;
    FILE *fout=fopen(fname,"wb");
    if(!fout){
        return 0;
    }

    size_t result=fwrite(cache.data(),1,cache.size(),fout);

    fclose(fout);
    return result;
}
void mem_file::swap(std::vector<byte_t> &e_data){
    if(is_cache()){
        cache.swap(e_data);
        idata=cache.data();
        isize=cache.size();
        offset=0;
    }
}

void mem_file::close(){
    if(fp)switch(state){
    case mem_file_state::WRITE_CACHE:
        fwrite(cache.data(),1,cache.size(),fp);
    case mem_file_state::READ_FILE:
    case mem_file_state::WRITE_FILE:
        fclose(fp);
    }
    state=mem_file_state::INVALID;
}

void mem_file::cache_read(){
    if(state==mem_file_state::READ_FILE){
        fseek(fp,0,SEEK_END);
        isize=_ftelli64(fp);
        fseek(fp,0,SEEK_SET);
        cache.resize(isize);
        fread(cache.data(),1,isize,fp);
        idata=cache.data();
        fclose(fp);
        offset=0;
        state=mem_file_state::READ_CACHE;
    }
}
const std::vector<mem_file::byte_t> &mem_file::data(){
    cache_read();
    return cache;
}

const char *mem_file::c_str() const{
    return fstr.c_str();
}

mem_file::operator bool() const{
    return is_valid();
}
bool mem_file::is_valid() const{
    return state!=mem_file_state::INVALID;
}
bool mem_file::is_read() const{
    return state==mem_file_state::READ_CACHE||state==mem_file_state::READ_FILE;
}
bool mem_file::is_write() const{
    return state==mem_file_state::WRITE_CACHE||state==mem_file_state::WRITE_FILE;
}
bool mem_file::is_cache() const{
    return state==mem_file_state::READ_CACHE||state==mem_file_state::WRITE_CACHE;
}

void mem_file::reserve(size_t rsize){
    cache.reserve(rsize);
}

int64_t mem_file::tell() const{
    if(!is_valid())return EOF;
    if(is_cache())return offset;
    return _ftelli64(fp);
}
int mem_file::seek(int64_t fpos,int forg){
    if(!is_valid())return EOF;
    if(is_cache()){
        if(forg==SEEK_SET)offset=fpos;
        else if(forg==SEEK_CUR)offset+=fpos;
        else if(forg==SEEK_END)offset=(int64_t)cache.size()+fpos;
        else return EOF;
        return size_t(offset)>cache.size()?EOF:0;
    }
    return _fseeki64(fp,fpos,forg);
}

size_t mem_file::read(void *buffer,size_t e_size,size_t e_count){
    if(state==mem_file_state::READ_FILE){
        return fread(buffer,e_size,e_count,fp);
    }
    else if(state==mem_file_state::READ_CACHE){
        if(e_size==0)return e_count;
        if(offset>=isize)return 0;
        size_t max_count=isize-offset;
        max_count/=e_size;
        if(e_count>max_count)e_count=max_count;

        size_t wsize=e_size*e_count;
        memcpy(buffer,idata+offset,wsize);
        offset+=wsize;
        return e_count;
    }
    return 0;
}
size_t mem_file::write(const void *buffer,size_t e_size,size_t e_count){
    if(state==mem_file_state::WRITE_FILE){
        return fwrite(buffer,e_size,e_count,fp);
    }
    else if(state==mem_file_state::WRITE_CACHE){
        if(e_size==0)return e_count;
        size_t max_count=-1-offset;
        max_count/=e_size;
        if(e_count>max_count)e_count=max_count;

        size_t wsize=e_size*e_count;
        size_t endpos=offset+wsize;
        if(endpos>cache.size())cache.resize(endpos);
        memcpy(cache.data()+offset,buffer,wsize);
        offset+=wsize;
        return e_count;
    }
    return 0;
}

std::string mem_file::readline(){
    std::string result;
    if(!is_read())return result;
    bool use_file=state==mem_file_state::READ_FILE;
    const int bufsize=16;
    char chbuf[bufsize];
    char *fret;
    do{
        do{
            fret=chbuf;
            int c;
            int n=bufsize;
            while(--n){
                if(use_file)c=fgetc(fp);
                else{
                    if(offset>=isize)c=EOF;
                    else {
                        c=idata[offset];
                        ++offset;
                    }
                }
                if(c==EOF)break;
                if(c=='\r')c='\n';
                *fret=c;
                ++fret;
                if(c=='\n')break;
            }
            *fret=0;
            if(c==EOF&&fret==chbuf){
                fret=nullptr;
                break;
            }
            result+=chbuf;
            if(result.back()=='\n'){
                result.pop_back();
                break;
            }
        } while(1);

        auto cpos=result.find('#');
        if(cpos==0)
            result.resize(cpos);

    } while(result.size()==0&&fret);
    return result;
}
bool mem_file::publish(const char *fname){
    if(state!=mem_file_state::WRITE_CACHE)return false;
    close();
    if(!fname){
        idata=cache.data();
        isize=cache.size();
        offset=0;
        state=mem_file_state::READ_CACHE;
    }
    else{
        auto it=mem_library.insert({std::string(fname),nullptr});
        mem_file *&pmem=it.first->second;
        if(!it.second)delete pmem;
        pmem=new mem_file;
        pmem->idata=cache.data();
        pmem->isize=cache.size();
        pmem->offset=0;
        pmem->state=mem_file_state::READ_CACHE;
        pmem->cache.swap(cache);
    }
    return true;
}

size_t fread(void *buffer,size_t e_size,size_t e_count,mem_file *mem){
    return mem->read(buffer,e_size,e_count);
}
size_t fwrite(const void *buffer,size_t e_size,size_t e_count,mem_file *mem){
    return mem->write(buffer,e_size,e_count);
}
std::string readline(mem_file &mem){
    return mem.readline();
}
int fprintf(mem_file *mem,const char *format,...){
    const int max_size=4096;
    char buffer[max_size];
    va_list args;
    va_start(args,format);
    int ret=vsnprintf(buffer,max_size,format,args);
    va_end(args);
    return fwrite(buffer,1,ret,mem);
}

void mem_file::read(const zipfile &zf){
    close();
    cache.resize(zf.filesize);
    zf.dumpfile(cache.data());
    idata=cache.data();
    isize=zf.filesize;
    offset=0;
    state=mem_file_state::READ_CACHE;
}
