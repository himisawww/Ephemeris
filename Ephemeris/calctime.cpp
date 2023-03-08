#include"utils.h"
#ifdef _MSC_VER
#include<windows.h>
double CalcTime(){
	static LARGE_INTEGER time_li;
	static double time_w=(QueryPerformanceFrequency(&time_li),(double)time_li.QuadPart);
	return (QueryPerformanceCounter(&time_li),(double)time_li.QuadPart/time_w);
}
#else
#include<time.h>
double CalcTime(){
	timespec ts;
	clock_gettime(CLOCK_MONOTONIC,&ts);
	return ts.tv_nsec*1e-9+ts.tv_sec;
}
#endif