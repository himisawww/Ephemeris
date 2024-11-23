#include"threadpool.h"

static thread_local ThreadPool *this_thread_pool=nullptr;
static thread_local size_t this_thread_id=ThreadPool::npos_tid;
static thread_local size_t this_stack_depth=0;
void ThreadPool::set_thread_pool(ThreadPool *pPool){
    this_thread_pool=pPool;
}
void ThreadPool::set_thread_id(size_t thread_id){
    this_thread_id=thread_id;
}

void ThreadPool::thread_local_pool_alloc(){
    if(!this_thread_pool){
        thread_local ThreadPool thread_pool;
        this_thread_pool=&thread_pool;
    }
}
void ThreadPool::thread_local_pool_free(){
    if(this_thread_pool&&(this_thread_id==npos_tid)){
        this_thread_pool->resize(0);
        this_thread_pool=nullptr;
    }
}

ThreadPool *ThreadPool::get_thread_pool(){
    return this_thread_pool;
}
size_t ThreadPool::get_thread_id(){
    return this==this_thread_pool?this_thread_id:npos_tid;
}
size_t ThreadPool::get_stack_depth(){
    return this==this_thread_pool?this_stack_depth:0;
}

ThreadPool::TaskQueue::TaskQueue():max_priority(0){

}
bool ThreadPool::TaskQueue::empty() const{
    return v.empty()||v[max_priority].empty();
}
size_t ThreadPool::TaskQueue::priority() const{
    return max_priority;
}
void ThreadPool::TaskQueue::push(const ThreadTask &t,size_t priority){
    if(priority>=v.size())
        v.resize(priority+1);
    if(priority>max_priority)
        max_priority=priority;
    v[priority].push_back(t);
}
ThreadPool::ThreadTask ThreadPool::TaskQueue::pop(){
    auto *pq=&v[max_priority];
    auto ret=pq->front();
    pq->pop_front();
    while(max_priority&&pq->empty())
        pq=&v[--max_priority];
    return ret;
}

void ThreadPool::thread_loop(ThreadPool *pPool,size_t thread_id,TaskGroup *pc){
    if(!pc){
        set_thread_pool(pPool);
        set_thread_id(thread_id);
    }
    do{
        if(pc&&!pc->load())
            return;
        bool has_task=false;
        size_t task_priority;
        ThreadPool::ThreadTask task;
        auto &assigned_queue=pPool->m_thread_tasks[thread_id];
        auto &task_queue=pPool->m_tasks;
        auto &busy_counter=pPool->m_busy;
        do{
            std::unique_lock<std::mutex> lock(pPool->m_mutex_distribute);
            size_t depth=pPool->get_stack_depth();
            do{
                if(!assigned_queue.empty()){
                    task_priority=1+assigned_queue.priority();
                    if(has_task=task_priority>=depth){
                        task=assigned_queue.pop();
                        break;
                    }
                }
                if(!task_queue.empty()){
                    task_priority=1+task_queue.priority();
                    if(has_task=task_priority>=depth){
                        task=task_queue.pop();
                        break;
                    }
                }
                if(pc){
                    if(!pc->load())return;
                    break;
                }
                else if(thread_id>=pPool->m_size.load()){
                    //printf("-------- thread id %llu exit\n",thread_id);
                    return;
                }
                pPool->m_distribute.wait(lock);
            } while(1);
            lock.unlock();
            pPool->m_distribute.notify_one();
        } while(!has_task);
        std::swap(task_priority,this_stack_depth);
        task.task_function(task.parameter,thread_id);
        std::swap(task_priority,this_stack_depth);
        bool do_notify=false;
        if(task.p_counter){
            if(!--*task.p_counter)do_notify=true;
        }
        if(!--busy_counter)do_notify=true;
        if(do_notify){
            std::lock_guard<std::mutex>(pPool->m_mutex_collect);
            pPool->m_collect.notify_all();
        }
    } while(1);
}

ThreadPool::ThreadPool(size_t n_threads){
    n_threads=std::min<size_t>(n_threads,std::thread::hardware_concurrency());
    m_busy.store(0);
    m_size.store(n_threads);
    resize_unchecked(0,n_threads);
}
ThreadPool::~ThreadPool(){
    resize(0);
}

size_t ThreadPool::size(){
    return m_size.load();
}
size_t ThreadPool::resize_unchecked(size_t old_size,size_t n_threads){
    m_threads.resize(n_threads);
    m_thread_tasks.resize(n_threads);
    for(size_t i=old_size;i<n_threads;++i)
        m_threads[i]=std::thread(thread_loop,this,i,nullptr);
    return old_size;
}
size_t ThreadPool::expand_unchecked(size_t n_threads){
    size_t old_size=m_size.load();
    if(n_threads>old_size){
        m_size.store(n_threads);
        resize_unchecked(old_size,n_threads);
    }
    return old_size;
}
size_t ThreadPool::resize(size_t n_threads){
    if(!wait_for_all())return npos_tid;
    size_t old_size=m_size.exchange(n_threads);
    if(n_threads<old_size){
        std::lock_guard<std::mutex>(this->m_mutex_distribute);
        m_distribute.notify_all();
        for(size_t i=n_threads;i<old_size;++i)
            if(m_threads[i].joinable())m_threads[i].join();
    }
    return resize_unchecked(old_size,n_threads);
}

void ThreadPool::add_task(TaskFunction f,void *param,TaskGroup *pc){
    {
        std::lock_guard<std::mutex> lock(m_mutex_distribute);
        expand_unchecked(1);
        m_busy.fetch_add(1);
        if(pc)
            pc->fetch_add(1);
        m_tasks.push({f,param,pc},get_stack_depth());
    }
    m_distribute.notify_one();
}
void ThreadPool::distribute_tasks(size_t n,TaskFunction f,void *param,TaskGroup *pc){
    {
        std::lock_guard<std::mutex> lock(m_mutex_distribute);
        expand_unchecked(n);
        m_busy.fetch_add(n);
        if(pc)
            pc->fetch_add(n);
        size_t depth=get_stack_depth();
        for(size_t i=0;i<n;++i)
            m_thread_tasks[i].push({f,param,pc},depth);
    }
    m_distribute.notify_all();
}
void ThreadPool::assign_task(size_t thread_id,TaskFunction f,void *param,TaskGroup *pc){
    std::lock_guard<std::mutex> lock(m_mutex_distribute);
    expand_unchecked(thread_id+1);
    m_busy.fetch_add(1);
    if(pc)
        pc->fetch_add(1);
    m_thread_tasks[thread_id].push({f,param,pc},get_stack_depth());
}
void ThreadPool::run(){
    m_distribute.notify_all();
}

bool ThreadPool::wait_for_all(TaskGroup *pc,std::function<void()> callback,double wakeup_seconds){
    size_t thread_id=get_thread_id();
    if(thread_id!=npos_tid){
        if(!pc)return false;
        ++this_stack_depth;
        thread_loop(this,thread_id,pc);
        --this_stack_depth;
        return true;//assert(!pc.load());
    }

    std::unique_lock<std::mutex> lock(m_mutex_collect);
    do{
        if(pc?!pc->load():!busy())break;
        if(callback){
            m_collect.wait_for(lock,std::chrono::duration<double>(std::max(minimum_wait_for,wakeup_seconds)));
            callback();
        }
        else
            m_collect.wait(lock);
    } while(1);
    return true;
}

size_t ThreadPool::busy(){
    return m_busy.load();
}
