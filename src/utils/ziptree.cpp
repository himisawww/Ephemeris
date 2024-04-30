#include"ziptree.h"
#include"crc32.h"
#include<vector>
#include"wcs_convert.h"

#pragma pack(push,2)
struct zipheader{
	uint32_t header;
	uint16_t version,gpbits,compress;
	uint32_t lasttime,crc32,compsize,uncompsize;
	uint16_t lenfname,lenexfield;
};
struct zip64header{
	uint16_t header,exsize;
	uint64_t uncompsize,compsize;
	uint64_t lhroffset;
	uint32_t disknum;
};
struct zipcentralheader{
	uint32_t header;
	uint16_t version_vendor,version,gpbits,compress;
	uint32_t lasttime,crc32,compsize,uncompsize;
	uint16_t lenfname,lenexfield,lencomment,disknum,inattrib;
	uint32_t exattrib,lhroffset;
};
struct zip64end{
	uint32_t header;
	uint64_t endsize;
	uint16_t version_vendor,version;
	uint32_t disknum,diskchs;
	uint64_t diskchsnum,totalchsnum,chssize,chsroffset;

	uint32_t headerloc,diskend;
	uint64_t endroffset;
	uint32_t diskcount;

	uint32_t header32;
	uint64_t dummy1,dummy2;
	uint16_t commentsize;
};
#pragma pack(pop)

void zipfile::validate(){
	bool success=false;

	do{
		FILE *fzip=pzip->fzip;
		if(!fzip)break;
		zipheader zh;
		_fseeki64(fzip,offset,SEEK_SET);
		if(fread(&zh,sizeof(zipheader),1,fzip)!=1)break;
		if(zh.header!=0x04034b50||zh.compress!=0)break;
		if(zh.compsize!=zh.uncompsize)break;
		if(zh.gpbits!=0&&zh.gpbits!=2048)break;
		fileoffset=offset+sizeof(zipheader)+zh.lenfname+zh.lenexfield;
		filesize=zh.compsize;

		std::string str;
		str.resize(zh.lenfname);
		if(fread(str.data(),zh.lenfname,1,fzip)!=1)break;
		filename=strtowcs(str,zh.gpbits);
		if(filename.size()==0)break;

		if(zh.lenexfield>=0x14){
			zip64header z64h;
			if(fread(&z64h,zh.lenexfield,1,fzip)!=1)break;
			if(z64h.header!=1||z64h.compsize!=z64h.uncompsize)break;
			filesize=z64h.compsize;
		}

		success=true;
	} while(0);
	if(!success)offset=-1;
	return;
}

zipfile &zipfile::operator++(){
	offset=fileoffset+filesize;
	validate();
	return *this;
}
zipfile zipfile::operator++(int){
	zipfile res(*this);
	this->operator++();
	return res;
}

std::wstring zipfile::fullname() const{
	return pzip->zippath+L"\\"+filename;
}

void zipfile::dumpfile(void *dst) const{
	FILE *fzip=pzip->fzip;
	_fseeki64(fzip,fileoffset,SEEK_SET);
	fread(dst,1,filesize,fzip);
}

zippack::zippack(const std::wstring &filename,bool construct){
	zippath=filename;
	fzip=_wfopen(filename.c_str(),construct?L"wb":L"rb");
	is_construct=construct;
}
zippack::~zippack(){
	if(!fzip)return;

	if(is_construct){
		//create zip file from zipmem
		_fseeki64(fzip,0,SEEK_SET);

		const uint32_t max_uint32=-1;
		std::vector<zipcentralheader> zchs(zipmems.size());
		std::vector<zip64header> z64hs(zipmems.size());
		std::vector<std::string> namestrs(zipmems.size());
		memset(zchs.data(),0,sizeof(zipcentralheader)*zchs.size());
		memset(z64hs.data(),0,sizeof(zip64header)*z64hs.size());

		for(size_t i=0;i<zipmems.size();++i){
			const auto &zm=zipmems[i];
			std::string &namestr=namestrs[i];
			namestr=wcstostr(zm.filename,true);
			if(!namestr.size())continue;

			zipheader zh;
			zip64header &z64h=z64hs[i];
			zipcentralheader &zch=zchs[i];

			memset(&zh,0,sizeof(zh));
			
			zh.header=0x04034b50;
			zch.header=0x02014b50;
			zch.gpbits=zh.gpbits=2048;
			zch.crc32=zh.crc32=crc32(zm.data.data(),zm.data.size(),0);
			size_t zsize=zm.data.size();

			z64h.header=1;
			z64h.compsize=z64h.uncompsize=zsize;
			z64h.lhroffset=_ftelli64(fzip);

			zch.lenfname=zh.lenfname=namestr.size();

			if(zsize>max_uint32){
				z64h.exsize=16;
				zh.lenexfield=z64h.exsize+4;
				zh.compsize=zh.uncompsize=max_uint32;
			}
			else{
				zh.compsize=zh.uncompsize=z64h.uncompsize;
			}

			fwrite(&zh,sizeof(zh),1,fzip);
			fwrite(namestr.data(),zh.lenfname,1,fzip);
			fwrite(&z64h,zh.lenexfield,1,fzip);
			fwrite(zm.data.data(),1,zm.data.size(),fzip);
		}

		zip64end ze;
		memset(&ze,0,sizeof(ze));
		ze.chsroffset=_ftelli64(fzip);

		for(size_t i=0;i<zipmems.size();++i){
			std::string &namestr=namestrs[i];
			zip64header &z64h=z64hs[i];
			zipcentralheader &zch=zchs[i];
			if(!namestr.size())continue;

			if(z64h.uncompsize>max_uint32||z64h.lhroffset>max_uint32){
				z64h.exsize=28;
				zch.lenexfield=z64h.exsize+4;
				zch.compsize=zch.uncompsize=max_uint32;
				zch.lhroffset=max_uint32;
				zch.disknum=-1;
			}
			else{
				zch.compsize=zch.uncompsize=z64h.uncompsize;
				zch.lhroffset=z64h.lhroffset;
			}

			fwrite(&zch,sizeof(zch),1,fzip);
			fwrite(namestr.data(),zch.lenfname,1,fzip);
			fwrite(&z64h,zch.lenexfield,1,fzip);
			++ze.totalchsnum;
		}

		if(ze.totalchsnum>0){
			ze.headerloc=0x07064b50;
			ze.endroffset=_ftelli64(fzip);
			ze.diskcount=1;

			ze.header=0x06064b50;
			ze.header32=0x06054b50;
			ze.dummy1=ze.dummy2=-1;
			ze.endsize=44;
			ze.chssize=ze.endroffset-ze.chsroffset;
			ze.diskchsnum=ze.totalchsnum;

			fwrite(&ze,sizeof(ze),1,fzip);
		}
		else{
			ze.header32=0x06054b50;
			ze.dummy1=ze.dummy2=0;
			fwrite(&ze.header32,22,1,fzip);
		}
	}
	fclose(fzip);
}
zipfile zippack::begin(){
	zipfile res;
	res.pzip=this;
	res.offset=0;

	//validate first file
	res.validate();
	return res;
}
zipfile zippack::end(){
	zipfile res;
	res.pzip=this;
	res.offset=-1;
	return res;
}
