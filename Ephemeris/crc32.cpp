#include"crc32.h"
#include<stdint.h>

uint32_t crc32(const void *dat,size_t bytes,uint32_t crc32val){
	static uint32_t table[0x100];
	if(!table[0]){
		for(size_t i=0;i<0x100;++i){
			uint32_t crci=i;
			for(int j=0;j<8;++j){
				crci=(crci&1?0:(uint32_t)0xEDB88320)^crci>>1;
			}
			table[i]=crci^(uint32_t)0xFF000000;
		}
	}
	for(size_t i=0;i<bytes;++i){
		crc32val=table[(uint8_t)crc32val^((uint8_t*)dat)[i]]^crc32val>>8;
	}
	return crc32val;
}
