#pragma once
#include<vector>
#include<mutex>
#include<condition_variable>
#include<deque>

class ThreadPool{
public:
    // task function
    // the prototype of functions that will be executed(called/run) by worker threads of the pool
    typedef void (*TaskFunction)(void *task_paramter,size_t thread_index);
    // counter of subtasks
    // usage:
    //  declare a local variable TaskGroup group;
    //  call add_task/distribute_tasks/assign_task with &group to add tasks to this group,
    //  then wait_for_all(&group) will block until all the tasks of this group are done.
    // combined with get_thread_pool(), this mechanism can also be used in worker threads of a thread pool, recursively.
    typedef std::atomic_size_t TaskGroup;

    static constexpr size_t npos_tid=-1;
private:
    struct ThreadTask{
        TaskFunction task_function;
        void *parameter;
        TaskGroup *p_counter;
    };

    class TaskQueue{
        std::vector<std::deque<ThreadTask>> v;
        size_t max_priority;
    public:
        TaskQueue();
        bool empty() const;
        size_t priority() const;
        void push(const ThreadTask &,size_t priority);
        ThreadTask pop();
    };

    std::vector<std::thread> m_threads;
    TaskQueue m_tasks;
    std::deque<TaskQueue> m_thread_tasks;
    std::mutex m_mutex;
    std::condition_variable m_distribute;
    std::condition_variable m_collect;
    std::atomic_size_t m_size;
    std::atomic_size_t m_busy;
    
    static void thread_loop(ThreadPool *pPool,size_t thread_id,TaskGroup *p_group);
    static inline void set_thread_pool(ThreadPool *);
    static inline void set_thread_id(size_t);

    // return thread_id if current thread is in this pool, npos_tid(-1) otherwise
    inline size_t get_thread_id();
    // return depth in stack if current thread is in this pool, 0 otherwise
    inline size_t get_stack_depth();
    // return size() before call
    inline size_t resize_unchecked(size_t old_size,size_t n_threads);
    inline size_t expand_unchecked(size_t n_threads);

    ThreadPool(ThreadPool&&)=delete;
public:
    ThreadPool(size_t n_threads=0);
    ~ThreadPool();

    //if called, subsequent functions in the same thread can use a thread pool
    // returned by get_thread_pool() to do tasks in parallel.
    //if current thread is a worker of a thread pool, calling this is unnecessary (and effectless),
    //  as get_thread_pool() will return its parent thread pool.
    static void thread_local_pool_alloc();
    //must call this in the same thread if thread_local_pool_alloc() is called,
    //  otherwise, dead wait will occur when the thread joins/exits.
    static void thread_local_pool_free();
    //return either the pool created by thread_local_pool_alloc(), or the pool to which current thread is in;
    //nullptr if none.
    static ThreadPool *get_thread_pool();

    // task can be consumed by any thread
    // if size()<1, a resize(1) will occur
    void add_task(TaskFunction,void *param,TaskGroup *p_group=nullptr);

    // n tasks will be distributed to first n threads
    // if size()<n, a resize(n) will occur
    void distribute_tasks(size_t n,TaskFunction,void *param,TaskGroup *p_group=nullptr);

    // assign a task to a specific thread
    // assigned tasks have higher priority than those added by add_task()
    // if size()<thread_id+1, a resize(thread_id+1) will occur
    // must call run() after all assignments are done,
    //  if not, dead sleep can occur.
    void assign_task(size_t thread_id,TaskFunction,void *param,TaskGroup *p_group=nullptr);
    // wake up all threads to consume tasks
    // must call this after assign_task() calls
    void run();

    // number of worker threads
    size_t size();
    // return npos_tid (and fails) if called from a worker thread of this
    // otherwise return size() before call
    // internally calls wait_for_all() to ensure no resize while busy
    size_t resize(size_t n_threads);
    // return number of unfinished tasks
    size_t busy();

    // block until all tasks are done
    // return true when success
    // return false if called from a worker thread of this without specifying a group (to avoid deadwait)
    bool wait_for_all(TaskGroup *p_group=nullptr);
};
