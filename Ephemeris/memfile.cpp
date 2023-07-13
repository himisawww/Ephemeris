#include"memfile.h"

mem_file::mem_file(){
    idata=nullptr;
    isize=0;
    offset=0;
}
mem_file::mem_file(const void *mem,size_t size){
    idata=(const uint8_t*)mem;
    isize=size;
    offset=0;
}
mem_file::mem_file(const char *fname){
    idata=nullptr;
    isize=0;
    offset=0;

    FILE *fin=fopen(fname,"rb");
    if(!fin){
        return;
    }
    
    fseek(fin,0,SEEK_END);
    isize=_ftelli64(fin);
    fseek(fin,0,SEEK_SET);
    wdata.resize(isize);
    fread(wdata.data(),1,isize,fin);
    idata=wdata.data();
    fclose(fin);
}
bool mem_file::save(const char *fname){
    FILE *fout=fopen(fname,"wb");
    if(!fout){
        return false;
    }

    fwrite(wdata.data(),1,wdata.size(),fout);

    fclose(fout);
    return true;
}

size_t fread(void *buffer,size_t e_size,size_t e_count,mem_file *mem){
    if(e_size==0)return e_count;
    if(mem->offset>=mem->isize)return 0;
    size_t max_count=mem->isize-mem->offset;
    max_count/=e_size;
    if(e_count>max_count)e_count=max_count;

    size_t wsize=e_size*e_count;
    memcpy(buffer,mem->idata+mem->offset,wsize);
    mem->offset+=wsize;
    return e_count;
}
size_t fwrite(const void *buffer,size_t e_size,size_t e_count,mem_file *mem){
    if(e_size==0)return e_count;
    size_t max_count=-1-mem->offset;
    max_count/=e_size;
    if(e_count>max_count)e_count=max_count;

    size_t wsize=e_size*e_count;
    size_t endpos=mem->offset+wsize;
    if(endpos>mem->wdata.size())mem->wdata.resize(endpos);
    memcpy(mem->wdata.data()+mem->offset,buffer,wsize);
    mem->offset+=wsize;
    return e_count;
}
