#include"threadpool.h"

void ThreadPool::thread_loop(ThreadPool *pPool,size_t thread_id){
    //printf("++++++++ thread id %llu start\n",thread_id);
    do{
        ThreadPool::ThreadTask task;
        auto &assigned_queue=pPool->m_thread_tasks[thread_id];
        auto &task_queue=pPool->m_tasks;
        auto &busy_counter=pPool->m_busy;
        {
            std::unique_lock<std::mutex> lock(pPool->m_mutex);
            do{
                if(!assigned_queue.empty()){
                    task=assigned_queue.front();
                    assigned_queue.pop_front();
                    break;
                }
                if(!task_queue.empty()){
                    task=task_queue.front();
                    task_queue.pop_front();
                    break;
                }
                pPool->m_collect.notify_all();
                if(thread_id>=pPool->m_size.load()){
                    //printf("-------- thread id %llu exit\n",thread_id);
                    return;
                }
                pPool->m_distribute.wait(lock);
            } while(1);
            lock.unlock();
            pPool->m_distribute.notify_one();
        }
        task.first(task.second,thread_id);
        --busy_counter;
    } while(1);
}

ThreadPool::ThreadPool(size_t n_threads){
    m_busy.store(0);
    m_size.store(0);
    resize(n_threads);
}
ThreadPool::~ThreadPool(){
    resize(0);
}

size_t ThreadPool::size(){
    return m_size.load();
}
size_t ThreadPool::resize(size_t n_threads){
    size_t old_size=m_size.exchange(n_threads);
    if(n_threads<old_size){
        m_distribute.notify_all();
        for(size_t i=n_threads;i<old_size;++i)
            if(m_threads[i].joinable())m_threads[i].join();
    }
    m_threads.resize(n_threads);
    m_thread_tasks.resize(n_threads);
    for(size_t i=old_size;i<n_threads;++i)
        m_threads[i]=std::thread(thread_loop,this,i);
    return old_size;
}

void ThreadPool::add_task(task_t f,void *param){
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_busy.fetch_add(1);
        m_tasks.push_back({f.task,param});
    }
    m_distribute.notify_one();
}
void ThreadPool::distribute_tasks(size_t n,task_t f,void *param){
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        if(size()<n)resize(n);
        m_busy.fetch_add(n);
        for(size_t i=0;i<n;++i)
            m_thread_tasks[i].push_back({f.task,param});
    }
    m_distribute.notify_all();
}
void ThreadPool::assign_task(size_t thread_id,task_t f,void *param){
    std::lock_guard<std::mutex> lock(m_mutex);
    if(size()<thread_id+1)resize(thread_id+1);
    m_busy.fetch_add(1);
    m_thread_tasks[thread_id].push_back({f.task,param});
}
void ThreadPool::run(){
    m_distribute.notify_all();
}

void ThreadPool::wait_for_all(){
    std::unique_lock<std::mutex> lock(m_mutex);
    do{
        if(!busy())return;
        m_collect.wait(lock);
    } while(1);
}

size_t ThreadPool::busy(){
    return m_busy.load();
}
