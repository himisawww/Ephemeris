#include"memio.h"
#include"htl/map.h"
#include"wcs_convert.h"

FILE *fopen(const std::string &fname,bool is_read){
    std::wstring wfname=strtowcs(fname,true);
    if(!wfname.size())return nullptr;
    return _wfopen(wfname.c_str(),is_read?L"rb":L"wb");
}
bool file_exist(const std::string &path){
    return MFILE(path).is_valid();
}
std::string get_file_name(const std::string &path){
    return path.substr(1+path.find_last_of("/\\"));
}
std::string get_file_extension(const std::string &path){
    auto subpos=path.rfind('.')+1;
    return path.substr(subpos&&!(path.find_first_of("/\\",subpos)+1)?subpos:path.size());
}

static htl::map<std::string,htl::vector<MFILE::byte_t>> mem_library;
static bool s_publish_invalid_ofile=false;
static std::string filepath_normalize(const std::string &fpath){
    std::string result;
    result.reserve(fpath.size());
    for(const char c:fpath){
        if(c=='/'||c=='\\'){
            if(!result.empty()&&result.back()=='/')
                continue;
            result.push_back('/');
        }
        else
            result.push_back(c);
    }
    return result;
}

bool MFILE::set_wcache_onfail(bool do_publish){
    bool ret=s_publish_invalid_ofile;
    s_publish_invalid_ofile=do_publish;
    return ret;
}
bool MFILE::get_wcache_onfail(){
    return s_publish_invalid_ofile;
}

MFILE::MFILE(MFILE &&mf):
    filename(std::move(mf.filename)),
    cached_data(std::move(mf.cached_data))
{
    fp=mf.fp;
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
    bool is_read_=is_read();
    bool is_cache_=is_cache();
    if(is_read_){
        auto it=mem_library.find(filepath_normalize(_fname));
        if(it!=mem_library.end()){
            const auto &pmem=it->second;
            idata=pmem.data();
            isize=pmem.size();
            offset=0;
            state=MFILE_STATE::READ_CACHE;
            return;
        }
    }
    else if(!is_write()){
        state=MFILE_STATE::INVALID;
        return;
    }

    fp=fopen(_fname,is_read_);
    if(!fp&&(is_read_||!s_publish_invalid_ofile))
        state=MFILE_STATE::INVALID;
    else if(is_read_){
        state=MFILE_STATE::READ_FILE;
        if(is_cache_)load_data();
    }
    else{
        if(!fp||is_cache_)offset=0;
        if(!fp&&!is_cache_)filename=_fname;
    }
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
    if(state==MFILE_STATE::INVALID)
        return EOF;
    if(state==MFILE_STATE::WRITE_CACHE&&fp)
        fwrite(cached_data.data(),1,cached_data.size(),fp);
    if(state==MFILE_STATE::WRITE_FILE&&!fp){
        state=MFILE_STATE::WRITE_CACHE;
        publish(filename);
    }
    MFILE_STATE ostate=state;
    state=MFILE_STATE::INVALID;
    if(ostate==MFILE_STATE::READ_CACHE||!fp)
        return 0;
    return fclose(fp);
}
MFILE::byte_t *MFILE::prepare(size_t new_cache_size){
    close();
    cached_data.resize(new_cache_size);
    if(cached_data.capacity()-new_cache_size>(new_cache_size>>1))
        cached_data.shrink_to_fit();
    state=MFILE_STATE::READ_CACHE;
    idata=cached_data.data();
    isize=new_cache_size;
    offset=0;
    return cached_data.data();
}
void MFILE::reset(){
    close();
    std::string().swap(filename);
    htl::vector<byte_t>().swap(cached_data);
    fp=0;
    offset=0;
    state=MFILE_STATE::WRITE_CACHE;
}

void MFILE::load_data(){
    if(state==MFILE_STATE::READ_CACHE){
        if(idata!=cached_data.data()||isize!=cached_data.size()){
            htl::vector<byte_t>(idata,idata+isize).swap(cached_data);
            idata=cached_data.data();
        }
    }
    else if(state==MFILE_STATE::READ_FILE){
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
    return state!=MFILE_STATE::INVALID;
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
    if(is_cache()||state==MFILE_STATE::WRITE_FILE&&!fp)return offset;
    return _ftelli64(fp);
}
int MFILE::seek(int64_t fpos,int forg){
    if(is_cache()||state==MFILE_STATE::WRITE_FILE&&!fp){
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
        if(fp)return fwrite(buffer,e_size,e_count,fp);
    }
    else if(state!=MFILE_STATE::WRITE_CACHE)
        return 0;

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
bool MFILE::publish(){
    if(state==MFILE_STATE::WRITE_CACHE)
        prepare(size());
    else if(is_read()){
        load_data();
        offset=0;
    }
    else
        return false;
    return true;
}
bool MFILE::publish(const std::string &fname){
    if(!publish())return false;
    if(fname.size()){
        auto &pmem=mem_library[filepath_normalize(fname)];
        pmem.swap(cached_data);
        idata=pmem.data();
        htl::vector<byte_t>().swap(cached_data);
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
    va_list args;
    va_start(args,format);
    std::string result=vstrprintf(format,args);
    va_end(args);
    return fwrite(result.data(),1,result.size(),mem);
}
std::string vstrprintf(const char *format,va_list _arg_list){
    char buffer[4096];
    size_t max_size=sizeof(buffer);
    va_list vcopy;
    va_copy(vcopy,_arg_list);
    int ret=vsnprintf(buffer,max_size,format,_arg_list);
    std::string result;
    if(ret>=0){
        size_t req_size=size_t(ret)+1;
        if(req_size>max_size){
            result.resize(req_size);
            int ret2=vsnprintf(result.data(),req_size,format,vcopy);
            result.resize(ret2!=ret?0:ret);
        }
        else result=buffer;
    }
    va_end(vcopy);
    return result;
}
std::string strprintf(const char *format,...){
    va_list args;
    va_start(args,format);
    std::string result=vstrprintf(format,args);
    va_end(args);
    return result;
}
