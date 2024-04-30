#pragma once
#include<stdint.h>

//for example: dat = dat1 + dat2
//crc32 of dat1 = crc32(dat1, length of dat1, 0)
//crc32 of dat = crc32(dat2, length of dat2, crc32 of dat1)
uint32_t crc32(const void *dat,size_t bytes,uint32_t crc32val=0);
