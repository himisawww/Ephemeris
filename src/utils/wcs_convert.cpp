#include"wcs_convert.h"
#define NOMINMAX
#include<Windows.h>

std::wstring strtowcs(std::string str,bool from_utf8,bool allow_error){
	std::wstring result;
	int codepage=from_utf8?CP_UTF8:CP_ACP;//UTF-8|GBK to UTF-16
	int dwflags=allow_error?0:MB_ERR_INVALID_CHARS;
	int reqsize=MultiByteToWideChar(codepage,dwflags,
		(LPCSTR)str.data(),str.size(),
		(LPWSTR)result.data(),0);
	if(reqsize>0){
		result.resize(reqsize);
		MultiByteToWideChar(codepage,dwflags,
			(LPCSTR)str.data(),str.size(),
			(LPWSTR)result.data(),result.size());
	}
	return result;
}
std::string wcstostr(std::wstring str,bool to_utf8,bool allow_error){
	std::string result;
	int codepage=to_utf8?CP_UTF8:CP_ACP;//UTF-8|GBK to UTF-16
	int dwflags=to_utf8?(allow_error?0:WC_ERR_INVALID_CHARS):0;
	BOOL haserr;
	LPBOOL useddef=to_utf8||allow_error?NULL:&haserr;
	int reqsize=WideCharToMultiByte(codepage,dwflags,
		(LPCWSTR)str.data(),str.size(),
		(LPSTR)result.data(),0,NULL,useddef);
	if(reqsize>0&&!(useddef&&haserr)){
		result.resize(reqsize);
		WideCharToMultiByte(codepage,dwflags,
			(LPCWSTR)str.data(),str.size(),
			(LPSTR)result.data(),result.size(),NULL,useddef);
	}
	return result;
}
