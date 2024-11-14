#include<cstdint>
class MFILE;

//like enum class, but can be used as uint8_t without static_cast
class LogTo{
    uint8_t e;
public:
    enum:uint8_t{
        VOID        =   0,
        STDOUT      =   1,  //  1<<0
        STDERR      =   2,  //  1<<1
        FILE        =   4,  //  1<<2
        ALL         =   uint8_t(-1)
    };
    LogTo(uint8_t _e=0):e(_e){}
    operator uint8_t &(){ return e; }
};

class LogLevel{
    uint8_t e;
public:
    enum:uint8_t{
        VERBOSE     =   0,
        DEBUG       =   10,
        INFO        =   20,
        WARNING     =   30,
        ERROR       =   40,
        CRITICAL    =   50,
        ANNOUNCEMENT=   uint8_t(-1)
    };
    LogLevel(uint8_t _e=0):e(_e){}
    operator uint8_t &(){ return e; }
};

class Logger{
public:
    typedef uint64_t LogLevelSettings;
private:
    LogTo target;
    MFILE *mbind;
    union{
        LogLevelSettings all_levels;
        LogLevel levels[8];
    };
public:
    // if target&LogTo::FILE and !_logfile, log to memory.
    Logger(LogTo _target=LogTo::STDOUT|LogTo::STDERR,
           const char *_logfile=nullptr,bool _use_cache=false);
    Logger(Logger &&)=delete;
    ~Logger();

    LogLevelSettings get_log_levels(){ return all_levels; }
    //return old levels for restore_log_levels
    LogLevelSettings set_log_levels(LogLevel,LogTo select=LogTo::ALL);
    //return old levels for restore_log_levels
    LogLevelSettings restore_log_levels(LogLevelSettings _all_levels);

    size_t logformat(LogLevel lv,const char *format,...);
};

extern Logger global_logger;
#define LogVerbose(...)  do{ global_logger.logformat(LogLevel::VERBOSE , __VA_ARGS__); } while(0)
#define LogDebug(...)    do{ global_logger.logformat(LogLevel::DEBUG   , __VA_ARGS__); } while(0)
#define LogInfo(...)     do{ global_logger.logformat(LogLevel::INFO    , __VA_ARGS__); } while(0)
#define LogWarning(...)  do{ global_logger.logformat(LogLevel::WARNING , __VA_ARGS__); } while(0)
#define LogError(...)    do{ global_logger.logformat(LogLevel::ERROR   , __VA_ARGS__); } while(0)
#define LogCritical(...) do{ global_logger.logformat(LogLevel::CRITICAL, __VA_ARGS__); } while(0)

#define LogAnnouncement  printf

class ScopedLogLevelSettings{
    Logger &logger;
    Logger::LogLevelSettings push_log_levels;
public:
    ScopedLogLevelSettings(Logger &_logger,LogLevel lv,LogTo select=LogTo::ALL)
        :logger(_logger){
        push_log_levels=logger.set_log_levels(lv,select);
    }
    ~ScopedLogLevelSettings(){
        logger.restore_log_levels(push_log_levels);
    }
};
