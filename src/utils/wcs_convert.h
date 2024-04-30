#pragma once
#include<string>

std::wstring strtowcs(std::string,bool from_utf8=false,bool allow_error=false);
std::string wcstostr(std::wstring,bool to_utf8=false,bool allow_error=false);
