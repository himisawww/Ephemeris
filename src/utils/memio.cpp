#include"memio.h"
#include<cstdarg>

#include"wcs_convert.h"

FILE *fopen(const std::string &fname,bool is_read){
    std::wstring wfname=strtowcs(fname,true);
    if(!wfname.size())return nullptr;
    return _wfopen(wfname.c_str(),is_read?L"rb":L"wb");
}
bool file_exist(const std::string &path){
    return MFILE(path).is_valid();
}

std::map<std::string,std::vector<MFILE::byte_t>> MFILE::mem_library;

MFILE::MFILE(MFILE &&mf):
    filename(std::move(mf.filename)),
    cached_data(std::move(mf.cached_data))
{
    fp=mf.fp;
    mf.fp=nullptr;
    idata=mf.idata;
    isize=mf.isize;
    offset=mf.offset;
    state=mf.state;
    mf.state=MFILE_STATE::INVALID;
}
MFILE::MFILE(){
    fp=nullptr;
    offset=0;
    state=MFILE_STATE::WRITE_CACHE;
}
MFILE::MFILE(const void *_mem,size_t _size){
    idata=(const byte_t *)_mem;
    isize=_size;
    offset=0;
    state=MFILE_STATE::READ_CACHE;
}
MFILE::MFILE(const std::string &_fname,MFILE_STATE _state){
    state=_state;
    if(state==MFILE_STATE::INVALID)return;
    bool is_read_=is_read();
    bool is_cache_=is_cache();
    filename=_fname;
    if(is_read_){
        auto it=mem_library.find(filename);
        if(it!=mem_library.end()){
            const auto &pmem=it->second;
            idata=pmem.data();
            isize=pmem.size();
            offset=0;
            state=MFILE_STATE::READ_CACHE;
            return;
        }
    }

    fp=fopen(_fname,is_read_);
    if(!fp)
        state=MFILE_STATE::INVALID;
    else if(is_read_){
        state=MFILE_STATE::READ_FILE;
        if(is_cache_)load_data();
    }
    else
        offset=0;
}

// open memory to write (WRITE_CACHE)
MFILE *mopen(){
    return new MFILE();
}
// open memory to read (READ_CACHE)
MFILE *mopen(const void *_mem,size_t _size){
    return new MFILE(_mem,_size);
}
// open file to read/write (read?READ:WRITE)_(cached_data?CACHE:FILE)
MFILE *mopen(const std::string &_fname,MFILE_STATE _state){
    MFILE *mf=new MFILE(_fname,_state);
    if(!mf->is_valid()){
        delete mf;
        return nullptr;
    }
    return mf;
}

int MFILE::close(){
    MFILE_STATE ostate=state;
    state=MFILE_STATE::INVALID;
    if(ostate==MFILE_STATE::INVALID)
        return EOF;
    if(ostate==MFILE_STATE::READ_CACHE||!fp)
        return 0;
    if(ostate==MFILE_STATE::WRITE_CACHE)
        fwrite(cached_data.data(),1,cached_data.size(),fp);
    return fclose(fp);
}
MFILE::byte_t *MFILE::reset(size_t new_cache_size){
    close();
    cached_data.resize(new_cache_size);
    state=MFILE_STATE::READ_CACHE;
    idata=cached_data.data();
    isize=new_cache_size;
    offset=0;
    return cached_data.data();
}

void MFILE::load_data(){
    if(state==MFILE_STATE::READ_FILE){
        offset=tell();
        seek(0,SEEK_END);
        isize=tell();
        seek(0,SEEK_SET);
        cached_data.resize(isize);
        fread(cached_data.data(),1,isize,fp);
        idata=cached_data.data();
        fclose(fp);
        state=MFILE_STATE::READ_CACHE;
    }
}

bool MFILE::is_valid() const{
    return this&&state!=MFILE_STATE::INVALID;
}
bool MFILE::is_read() const{
    return state==MFILE_STATE::READ_CACHE||state==MFILE_STATE::READ_FILE;
}
bool MFILE::is_write() const{
    return state==MFILE_STATE::WRITE_CACHE||state==MFILE_STATE::WRITE_FILE;
}
bool MFILE::is_cache() const{
    return state==MFILE_STATE::READ_CACHE||state==MFILE_STATE::WRITE_CACHE;
}

void MFILE::reserve(size_t rsize){
    cached_data.reserve(rsize);
}

int64_t MFILE::tell() const{
    if(is_cache())return offset;
    return _ftelli64(fp);
}
int MFILE::seek(int64_t fpos,int forg){
    if(is_cache()){
        size_t s=size();
        if(forg==SEEK_SET)offset=fpos;
        else if(forg==SEEK_CUR)offset+=fpos;
        else if(forg==SEEK_END)offset=(int64_t)s+fpos;
        else return EOF;
        return size_t(offset)>s?EOF:0;
    }
    return _fseeki64(fp,fpos,forg);
}
int fseek(MFILE *_Stream,int64_t _Offset,int _Origin){
    return _Stream->seek(_Offset,_Origin);
}
int64_t ftell(MFILE *_Stream){
    return _Stream->tell();
}
int fclose(MFILE *_Stream){
    int ret=_Stream->close();
    delete _Stream;
    return ret;
}

size_t MFILE::read(void *buffer,size_t e_size,size_t e_count){
    if(state==MFILE_STATE::READ_FILE){
        return fread(buffer,e_size,e_count,fp);
    }
    else if(state==MFILE_STATE::READ_CACHE){
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
size_t MFILE::write(const void *buffer,size_t e_size,size_t e_count){
    if(state==MFILE_STATE::WRITE_FILE){
        return fwrite(buffer,e_size,e_count,fp);
    }
    else if(state==MFILE_STATE::WRITE_CACHE){
        if(e_size==0)return e_count;
        size_t max_count=-1-offset;
        max_count/=e_size;
        if(e_count>max_count)e_count=max_count;

        size_t wsize=e_size*e_count;
        size_t endpos=offset+wsize;
        if(endpos>cached_data.size())cached_data.resize(endpos);
        memcpy(cached_data.data()+offset,buffer,wsize);
        offset+=wsize;
        return e_count;
    }
    return 0;
}

std::string MFILE::readline(){
    std::string result;
    if(!is_read())return result;
    bool use_file=state==MFILE_STATE::READ_FILE;
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
bool MFILE::publish(const std::string &fname){
    if(state!=MFILE_STATE::WRITE_CACHE)return false;
    reset(size());
    if(fname.size()){
        auto &pmem=mem_library[fname];
        pmem.swap(cached_data);
        idata=pmem.data();
        cached_data.clear();
    }
    return true;
}

size_t fread(void *buffer,size_t e_size,size_t e_count,MFILE *mem){
    return mem->read(buffer,e_size,e_count);
}
size_t fwrite(const void *buffer,size_t e_size,size_t e_count,MFILE *mem){
    return mem->write(buffer,e_size,e_count);
}
std::string readline(MFILE *mem){
    return mem->readline();
}
int fprintf(MFILE *mem,const char *format,...){
    const int max_size=4096;
    char buffer[max_size];
    va_list args;
    va_start(args,format);
    int ret=vsnprintf(buffer,max_size,format,args);
    va_end(args);
    return fwrite(buffer,1,ret,mem);
}