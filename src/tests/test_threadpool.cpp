#include"utils/threadpool.h"
#include"utils/logger.h"
#include"utils/calctime.h"

#define TEST_N 16
#define TEST_RANGE 10000000
#define TEST_RECURSE_LIMIT 2000
#define TEST_MOD 0xC0000001

int test_threadpool(){
    class test_task{
        uint64_t _begin,_end,_result;
    public:
        test_task()=default;
        test_task(uint64_t b,uint64_t e):_begin(b),_end(e),_result(1){}
        uint64_t get_result() const{ return _result; }
        void work_single(){
            for(uint64_t i=_begin;i!=_end;++i)
                _result=_result*i%TEST_MOD;
        }
        static void work(void *p,size_t){
            test_task &_this=*(test_task*)p;
            uint64_t range=_this._end-_this._begin;
            if(range>TEST_RECURSE_LIMIT){
                auto *p=ThreadPool::get_thread_pool();
                uint64_t _mid=_this._begin+range/2;
                test_task l(_this._begin,_mid),r(_mid,_this._end);
                ThreadPool::TaskGroup g;
                p->add_task(work,&l,&g);
                p->add_task(work,&r,&g);
                p->wait_for_all(&g);
                _this._result=_this._result*l._result%TEST_MOD;
                _this._result=_this._result*r._result%TEST_MOD;
            }
            else _this.work_single();
        }
    };

    test_task base(1,1+TEST_RANGE);
    double s=CalcTime();
    ThreadPool::LocalGuard _;
    auto *p=ThreadPool::get_thread_pool();
    test_task param[TEST_N];
    for(int i=0;i<TEST_N;++i){
        param[i]=base;
        p->add_task(test_task::work,param+i);
    }
    p->wait_for_all();
    s=CalcTime()-s;
    base.work_single();
    for(int i=0;i<TEST_N;++i)if(param[i].get_result()!=base.get_result()){
        LogError(
            "\nResult Mismatch: 0x%llx != 0x%llx",param[i].get_result(),base.get_result());
        return 1;
    }
    LogInfo("\n      Passed(0x%llx, %fs), ",base.get_result(),s);
    return 0;
}
