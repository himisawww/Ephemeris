#include"zipio.h"
#include<algorithm>
#include"crc32.h"
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
	//zip64 end of central directory, 56 bytes
	uint32_t header;
	uint64_t endsize;
	uint16_t version_vendor,version;
	uint32_t disknum,diskchs;
	uint64_t diskchsnum,totalchsnum,chssize,chsroffset;

	//zip64 locator, 20 bytes
	uint32_t headerloc,diskend;
	uint64_t endroffset;
	uint32_t diskcount;

	//zip end of central directory, 22 bytes
	uint32_t header32;
	union{
		struct{
			uint64_t dummy1,dummy2;
		};
		struct{
			uint16_t disknum32,diskchs32,diskchsnum32,totalchsnum32;
			uint32_t chssize32,chsroffset32;
		};
	};
	uint16_t commentsize;
private:
	template<typename T,typename U>
	static bool check_valid(const T &v32,U v64){ return v32==T(-1)||U(v32)==v64; }
public:
	bool check_valid() const{
		return check_valid(disknum32,disknum)
			 &&check_valid(diskchs32,diskchs)
			 &&check_valid(diskchsnum32,diskchsnum)
			 &&check_valid(totalchsnum32,totalchsnum)
			 &&check_valid(chssize32,chssize)
			 &&check_valid(chsroffset32,chsroffset);
	}
};
#pragma pack(pop)

izipfile::izipfile(const izippack *_pzip,size_t _offset):pzip(_pzip),locoffset(_offset){
	filesize=fileoffset=npos;
	if(locoffset==npos)return;
	bool success=false;

	do{
		MFILE *fzip=pzip->fzip;
		if(!fzip)break;
		zipheader zh;
		if(fseek(fzip,locoffset,SEEK_SET))break;
		if(fread(&zh,sizeof(zipheader),1,fzip)!=1)break;
		if(zh.header!=0x04034b50||zh.compress!=0)break;
		if(zh.compsize!=zh.uncompsize)break;
		if(zh.gpbits!=0&&zh.gpbits!=2048)break;
		fileoffset=locoffset+sizeof(zipheader)+zh.lenfname+zh.lenexfield;
		filesize=zh.compsize;

		std::string str;
		str.resize(zh.lenfname);
		if(fread(str.data(),zh.lenfname,1,fzip)!=1)break;
		if(!zh.gpbits)
			str=wcstostr(strtowcs(str,false),true);
		filename=str;
		if(filename.size()==0)break;

		if(load_zip64(zh.lenexfield,false))break;

		size_t nextoffset=fileoffset+filesize;
		if(nextoffset<fileoffset||fseek(fzip,nextoffset,SEEK_SET))break;

		success=true;
	} while(0);
	if(!success)filesize=fileoffset=locoffset=npos;
}
izipfile::izipfile(const izippack *_pzip):pzip(_pzip),fileoffset(npos){
	filesize=npos;
	bool success=false;

	do{
		MFILE *fzip=pzip->fzip;
		if(!fzip)break;
		zipcentralheader ch;
		if(fread(&ch,sizeof(zipcentralheader),1,fzip)!=1)break;
		if(ch.header!=0x02014b50||ch.compress!=0)break;
		if(ch.compsize!=ch.uncompsize)break;
		if(ch.gpbits!=0&&ch.gpbits!=2048)break;
		locoffset=ch.lhroffset;
		filesize=ch.compsize;

		std::string str;
		str.resize(ch.lenfname);
		if(fread(str.data(),ch.lenfname,1,fzip)!=1)break;
		if(!ch.gpbits)
			str=wcstostr(strtowcs(str,false),true);
		filename=str;
		if(filename.size()==0)break;

		int64_t nextoffset=ftell(fzip);
		if(nextoffset<0)break;

		if(load_zip64(ch.lenexfield,true))break;

		if(fseek(fzip,nextoffset+ch.lenexfield+ch.lencomment,SEEK_SET))break;

		success=true;
	} while(0);
	if(!success)filesize=locoffset=npos;
}
uint16_t izipfile::load_zip64(uint16_t exremain,bool is_central){
	bool size64=filesize==uint32_t(-1);
	bool offset64=is_central&&locoffset==uint32_t(-1);
	if(!(size64||offset64))return 0;

	//maybe zip64
	MFILE *fzip=pzip->fzip;
	zip64header z64h;
	while(exremain>=4){
		if(fread(&z64h,4,1,fzip)!=1)break;
		if(exremain<z64h.exsize+4)
			break;
		if(z64h.header==1)
		{
			auto loadentry=[&](bool may_exist,uint64_t &entry,uint64_t def){
				if(may_exist&&z64h.exsize>=8){
					if(fread(&entry,8,1,fzip)!=1)return false;
					z64h.exsize-=8;
				}
				else entry=def;
				return true;
			};
			if(!loadentry(size64,z64h.uncompsize,filesize))break;
			if(!loadentry(size64,z64h.compsize,filesize))break;
			if(!loadentry(offset64,z64h.lhroffset,locoffset))break;
			if(z64h.compsize!=z64h.uncompsize)break;
			filesize=z64h.compsize;
			locoffset=z64h.lhroffset;
			exremain=0;
			break;
		}
		if(fseek(fzip,z64h.exsize,SEEK_CUR))break;
		exremain-=z64h.exsize+4;
	}
	return exremain;
}

izipfile &izipfile::operator++(){
	*this=izipfile(pzip,is_ready()?fileoffset+filesize:npos);
	return *this;
}
izipfile izipfile::operator++(int){
	izipfile res(*this);
	this->operator++();
	return res;
}

std::string izipfile::fullname() const{
	return pzip->fzip->get_name()+"/"+filename;
}

bool izipfile::dumpfile(MFILE &mf) const{
	bool success=false;
	do{
		if(!is_ready())break;
		void *dst=mf.prepare(size());
		MFILE *fzip=pzip->fzip;
		if(fseek(fzip,fileoffset,SEEK_SET))break;
		if(fread(dst,1,filesize,fzip)!=filesize)break;
		mf.set_name(filename);
		success=true;
	} while(0);
	if(!success)
		mf.prepare(0);
	return success;
}
bool izipfile::fetch(){
	if(fileoffset==npos&&locoffset!=npos)
		*this=izipfile(pzip,locoffset);
	return is_ready();
}

izippack::izippack(izippack &&_i){
	fzip=_i.fzip;
	_i.fzip=nullptr;
}
izippack::izippack(const std::string &filename){
	fzip=mopen(filename);
}
izippack::~izippack(){
	if(fzip)fclose(fzip);
}
izipfile izippack::begin() const{
	return izipfile(this,0);
}
izipfile izippack::end() const{
	return izipfile(this,izipfile::npos);
}
std::vector<izipfile> izippack::load_central_directory(){
	std::vector<izipfile> result;
	if(!fzip||fseek(fzip,0,SEEK_END))return result;
	int64_t fremain=ftell(fzip);
	size_t fsize=fremain;

	constexpr size_t buffer_size=4096;
	constexpr size_t min_offset=sizeof(zip64end);
	char buf[buffer_size];
	char *const data=buf+buffer_size;

	zip64end ze;
	ze.header32=0;
	size_t boffset=0;
	while(fremain>0/*&&fremain+0x10000>fsize*/){
		size_t readsize=std::min(size_t(fremain),buffer_size-boffset);
		fremain-=readsize;
		if(fseek(fzip,fremain,SEEK_SET))break;
		size_t bend=boffset+readsize;
		if(fread(data-bend,1,readsize,fzip)!=readsize)break;
		for(size_t i=boffset;i<bend;){
			char *ph=data-++i;
			if(i<22||ph[0]!=0x50||ph[1]!=0x4b||ph[2]!=0x05||ph[3]!=0x06)
				continue;
			memcpy(&ze.header32,ph,22);
			//try load zip64end
			bool has_zip64=false;
			size_t zipend=fremain+bend-i;
			do{
				if(zipend<76||fseek(fzip,zipend-20,SEEK_SET))break;
				if(fread(&ze.headerloc,20,1,fzip)!=1)break;
				if(ze.headerloc!=0x07064b50)break;
				if(ze.diskcount!=1||ze.diskend!=0)break;
				if(ze.endroffset>zipend-76||fseek(fzip,ze.endroffset,SEEK_SET))break;
				if(fread(&ze.header,56,1,fzip)!=1)break;
				if(ze.header!=0x06064b50||ze.endroffset+ze.endsize+32!=zipend)break;
				has_zip64=true;
			} while(0);
			if(!has_zip64){
				ze.disknum=ze.disknum32;
				ze.diskchs=ze.diskchs32;
				ze.diskchsnum=ze.diskchsnum32;
				ze.totalchsnum=ze.totalchsnum32;
				ze.chssize=ze.chssize32;
				ze.chsroffset=ze.chsroffset32;
			}
			fremain=0;
			break;
		}
		if(bend>min_offset){
			size_t clipsize=bend-min_offset;
			std::memmove(data-min_offset,data-bend,min_offset);
			bend=min_offset;
		}
		boffset=bend;
	}

	do{
		if(ze.header32!=0x06054b50||!ze.check_valid())
			break;
		if(ze.diskchsnum!=ze.totalchsnum)
			break;
		if(fsize<ze.chsroffset||fsize<ze.chssize||fsize<ze.chsroffset+ze.chssize)
			break;
		if(ze.chssize/sizeof(zipcentralheader)<ze.totalchsnum)
			break;
		if(fseek(fzip,ze.chsroffset,SEEK_SET))
			break;

		result.reserve(ze.totalchsnum);
		size_t i=0;
		while(i!=ze.totalchsnum){
			if(!result.emplace_back(izipfile(this)))
				break;
			++i;
		}
		if(i!=ze.totalchsnum)
			std::vector<izipfile>().swap(result);
		std::sort(result.begin(),result.end());
	} while(0);
	return result;
}

ozippack::ozippack(ozippack &&_o){
	swap(_o);
	fzip=_o.fzip;
	_o.fzip=nullptr;
}
ozippack::ozippack(const std::string &filename){
	fzip=mopen(filename,MFILE_STATE::WRITE_FILE);
}
ozippack::~ozippack(){
	if(!fzip)return;

	std::vector<MFILE> &zipmems=static_cast<std::vector<MFILE> &>(*this);
	const size_t nzips=zipmems.size();
	//create zip file from zipmem
	fseek(fzip,0,SEEK_SET);

	const uint32_t max_uint32=-1;
	std::vector<zipcentralheader> zchs(nzips);
	std::vector<zip64header> z64hs(nzips);
	std::vector<std::string> zipnames(nzips);

	memset(zchs.data(),0,sizeof(zipcentralheader)*nzips);
	memset(z64hs.data(),0,sizeof(zip64header)*nzips);
	
	for(size_t i=0;i<nzips;++i){
		zipmems[i].load_data();
		const auto &zm=zipmems[i];
		if(zm.is_cache())zipnames[i]=zm.get_name();
		const std::string &namestr=zipnames[i];
		if(!namestr.size())continue;

		zipheader zh;
		zip64header &z64h=z64hs[i];
		zipcentralheader &zch=zchs[i];
		size_t zsize=zm.size();
	
		memset(&zh,0,sizeof(zh));
		
		zh.header=0x04034b50;
		zch.header=0x02014b50;
		zch.gpbits=zh.gpbits=2048;
		zch.crc32=zh.crc32=crc32(zm.data(),zsize,0);
	
		z64h.header=1;
		z64h.compsize=z64h.uncompsize=zsize;
		z64h.lhroffset=ftell(fzip);
	
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
		fwrite(zm.data(),1,zsize,fzip);
		zipmems[i].reset();
	}
	
	zipmems.clear();

	zip64end ze;
	memset(&ze,0,sizeof(ze));
	ze.chsroffset=ftell(fzip);
	
	for(size_t i=0;i<nzips;++i){
		const std::string &namestr=zipnames[i];
		if(!namestr.size())continue;
		zip64header &z64h=z64hs[i];
		zipcentralheader &zch=zchs[i];
	
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
		ze.endroffset=ftell(fzip);
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
	
	fclose(fzip);
}
void ozippack::push_back(std::string fullname,MFILE &&mf){
	emplace_back(std::move(mf)).set_name(fullname);
}
void ozippack::push_back(const izipfile &zf){
	MFILE mf;
	zf.dumpfile(mf);
	emplace_back(std::move(mf));
}
void ozippack::push_back(const izippack &izip){
	for(const izipfile &zf:izip)
		push_back(zf);
}
