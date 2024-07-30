#include"logger.h"
#include<cstdarg>

Logger global_logger;

Logger::Logger(LogTo _target,const char *_logfile,MFILE_STATE _filestate){
    target=_target;
    if(target&LogTo::FILE)
        mbind=mopen(_logfile,_filestate);
    else
        mbind=nullptr;
    if(!mbind)
        target&=~LogTo::FILE;
    //default to verbose
    all_levels=0;
    //higher than debug goes to stdout
    levels[0]=LogLevel::DEBUG;
    //higher than warning goes to stderr
    levels[1]=LogLevel::WARNING;
}

Logger::~Logger(){
    if(mbind)fclose(mbind);
}

uint64_t Logger::set_log_levels(LogLevel lv,LogTo select){
    uint64_t ret=all_levels;
    for(int i=0;i<8;++i)if(target&select&1<<i)
        levels[i]=lv;
    return ret;
}
uint64_t Logger::restore_log_levels(uint64_t _all_levels){
    uint64_t ret=all_levels;
    all_levels=_all_levels;
    return ret;
}

size_t Logger::logformat(LogLevel lv,const char *format,...){
    bool log_out  = lv>=levels[0]&&(target&1<<0);
    bool log_err  = lv>=levels[1]&&(target&1<<1);
    bool log_file = lv>=levels[2]&&(target&1<<2);
    if(!(log_out||log_err||log_file))return 0;
    va_list args;
    va_start(args,format);
    std::string result=vstrprintf(format,args);
    va_end(args);
    const char *buffer=result.data();
    const size_t ret=result.size();
    if(log_err)
        fwrite(buffer,1,ret,stderr);
    else if(log_out)
        fwrite(buffer,1,ret,stdout);
    if(mbind&&log_file)
        fwrite(buffer,1,ret,mbind);
    return ret;
}
