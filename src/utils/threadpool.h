#pragma once
#include<vector>
#include<mutex>
#include<condition_variable>
#include<deque>

class ThreadPool{
    typedef void (*TaskFunction)(void*);
    typedef void (*IndexedTaskFunction)(void*,size_t);
    struct task_t{
        IndexedTaskFunction task;
        task_t(IndexedTaskFunction f):task(f){}
        task_t(TaskFunction f):task((IndexedTaskFunction)f){}
    };
    typedef void *TaskParameter;
    typedef std::pair<IndexedTaskFunction,TaskParameter> ThreadTask;
private:
    std::vector<std::thread> m_threads;
    std::deque<ThreadTask> m_tasks;
    std::deque<std::deque<ThreadTask>> m_thread_tasks;
    std::mutex m_mutex;
    std::condition_variable m_distribute;
    std::condition_variable m_collect;
    std::atomic_size_t m_size;
    std::atomic_size_t m_busy;
    static void thread_loop(ThreadPool *pPool,size_t thread_id);

    ThreadPool(ThreadPool&&)=delete;
public:
    ThreadPool(size_t n_threads=0);
    ~ThreadPool();

    // task can be consumed by any thread
    void add_task(task_t,void *param=nullptr);

    // n tasks will be distributed to first n threads
    // if size()<n, a resize(n) will occur
    void distribute_tasks(size_t n,task_t,void *param=nullptr);

    // assign a task to a specific thread
    // assigned tasks have higher priority than those added by add_task()
    // if size()<thread_id+1, a resize(thread_id+1) will occur
    // must call run() after all assignments are done,
    //  if not, dead sleep can occur.
    void assign_task(size_t thread_id,task_t,void *param=nullptr);
    // wake up all threads to consume tasks
    void run();

    // number of worker threads
    size_t size();
    // return size() before call
    size_t resize(size_t n_threads);
    // return number of unfinished tasks
    size_t busy();

    // block until all tasks are done
    void wait_for_all();
};
