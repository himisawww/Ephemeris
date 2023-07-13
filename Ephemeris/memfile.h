#pragma once
#include<vector>
#include<cstdint>

class mem_file{
public:
    std::vector<uint8_t> wdata;
    const uint8_t *idata;
    size_t isize,offset;

    mem_file();
    mem_file(const void *mem,size_t size);
    mem_file(const char *fname);

    bool save(const char *fname);
};
size_t fread(void *buffer,size_t e_size,size_t e_count,mem_file *mem);
size_t fwrite(const void *buffer,size_t e_size,size_t e_count,mem_file *mem);
