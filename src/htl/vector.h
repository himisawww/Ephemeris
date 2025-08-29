#pragma once
#include"base.h"
#include<memory>

namespace htl{

template<typename T,bool=::std::is_scalar_v<T>>
struct _vector_inplace_trait;
template<typename T>
struct _vector_inplace_trait<T,true>{
    static constexpr size_t inplace_size=16/sizeof(T);
    static constexpr size_t type_align=alignof(T);
    static constexpr size_t type_size(size_t count){ return sizeof(T)*count; }
};
template<typename T>
struct _vector_inplace_trait<T,false>{
    static constexpr size_t inplace_size=0;
    static constexpr size_t type_align=1;
    static constexpr size_t type_size(size_t count){ return 1; }
};

template<typename T,int InplaceSize=-1,bool SizeFixed=false,typename Allocator=::std::allocator<T>>
class vector{
    typedef ::std::allocator_traits<Allocator> alloc_traits;
public:
    typedef Allocator   allocator_type;
    typedef         T       value_type;
    typedef         T       &reference;
    typedef   const T &const_reference;
    typedef typename alloc_traits::pointer         pointer;
    typedef typename alloc_traits::const_pointer   const_pointer;
    typedef typename alloc_traits::size_type       size_type;
    typedef typename alloc_traits::difference_type difference_type;
private:
    template<bool Const>
    class _iterator{
        friend class _iterator<!Const>;
        friend class vector;
        vector::pointer _ptr;
    public:
        typedef vector::difference_type difference_type;
        typedef vector::value_type      value_type;
        typedef ::std::conditional_t<Const,vector::const_pointer,  vector::pointer>   pointer;
        typedef ::std::conditional_t<Const,vector::const_reference,vector::reference> reference;
        typedef ::std::random_access_iterator_tag iterator_category;
    public:
        template<bool C2=!Const,typename=::std::enable_if_t<C2>>
        operator _iterator<C2>() const{
            _iterator<true> it;
            it._ptr=_ptr;
            return it;
        }

        bool operator==(const _iterator &_other) const{ return _ptr==_other._ptr; }
        bool operator!=(const _iterator &_other) const{ return _ptr!=_other._ptr; }
        bool operator< (const _iterator &_other) const{ return _ptr< _other._ptr; }
        bool operator<=(const _iterator &_other) const{ return _ptr<=_other._ptr; }
        bool operator> (const _iterator &_other) const{ return _ptr> _other._ptr; }
        bool operator>=(const _iterator &_other) const{ return _ptr>=_other._ptr; }

        _iterator &operator++(){
            ++_ptr;
            return *this;
        }
        _iterator &operator--(){
            --_ptr;
            return *this;
        }
        _iterator operator++(int){
            _iterator _tmp=*this;
            ++*this;
            return _tmp;
        }
        _iterator operator--(int){
            _iterator _tmp=*this;
            --*this;
            return _tmp;
        }

        _iterator &operator+=(const difference_type _offset){
            _ptr+=_offset;
            return *this;
        }
        _iterator operator+(const difference_type _offset) const{
            _iterator _tmp=*this;
            _tmp+=_offset;
            return _tmp;
        }
        friend _iterator operator+(const difference_type _offset,_iterator _it){
            _it+=_offset;
            return _it;
        }
        _iterator &operator-=(const difference_type _offset){
            _ptr-=_offset;
            return *this;
        }
        _iterator operator-(const difference_type _offset) const{
            _iterator _tmp=*this;
            _tmp-=_offset;
            return _tmp;
        }
        template<bool C2>
        difference_type operator-(const _iterator<C2> &_other) const{
            return _ptr-_other._ptr;
        }

        reference operator *() const{
            return *_ptr;
        }
        pointer operator ->() const{
            return  _ptr;
        }
        reference operator[](const difference_type _offset) const{
            return *(_ptr+_offset);
        }
    };
public:
    typedef _iterator<false>       iterator;
    typedef _iterator< true> const_iterator;
    typedef ::std::reverse_iterator<      iterator>       reverse_iterator;
    typedef ::std::reverse_iterator<const_iterator> const_reverse_iterator;

private:
    static_assert(size_type(-1)>size_type(0)&&sizeof(size_type)==8,"size should be 64-bit unsigned type.");
    static_assert(::std::is_same_v<T,typename Allocator::value_type>,"value_type of vector and allocator mismatch.");

    typedef struct{
        size_type _size:63;
        size_type _alloc:1;
    } _size_flag;
    compressed_pair<Allocator,_size_flag> _a;

    static constexpr size_type inplace_size=InplaceSize>=0?size_t(InplaceSize):_vector_inplace_trait<T>::inplace_size;
    union{
        alignas(_vector_inplace_trait<T,inplace_size!=0>::type_align) unsigned char
        _buffer[_vector_inplace_trait<T,inplace_size!=0>::type_size(inplace_size)];
        struct{
            pointer _b;
            size_type _c;
        };
    } _u;

    template<typename V=value_type,typename A=alloc_traits,bool I=bool(inplace_size)>
    struct value_traits{
        static constexpr bool just_copy_constructible=!::std::is_move_constructible_v<V>&&::std::is_copy_constructible_v<V>;
        static constexpr bool nothrow_constructible=just_copy_constructible?::std::is_nothrow_copy_constructible_v<V>: ::std::is_nothrow_move_constructible_v<V>;
        static constexpr bool just_copy_assignable=!::std::is_move_assignable_v<V>&&::std::is_copy_assignable_v<V>;
        static constexpr bool nothrow_assignable=just_copy_assignable?::std::is_nothrow_copy_assignable_v<V>: ::std::is_nothrow_move_assignable_v<V>;
        static constexpr bool is_trivially_swappable=::std::is_trivially_destructible_v<V>&&::std::is_trivially_move_constructible_v<V>&&::std::is_trivially_move_assignable_v<V>;

        static constexpr bool nothrow_move_constructible=!I||nothrow_constructible;
        static constexpr bool nothrow_move_constructible_with_allocator=!I&&A::is_always_equal::value||nothrow_constructible;
        static constexpr bool nothrow_move_assignable=!I&&(A::propagate_on_container_move_assignment::value||A::is_always_equal::value)||nothrow_constructible&&nothrow_assignable;
        static constexpr bool nothrow_swappable=!I&&(A::propagate_on_container_swap::value||A::is_always_equal::value)||nothrow_constructible&&nothrow_assignable;
    };

    static decltype(auto) _move_if_constructible(T &val){ return ::std::conditional_t<value_traits<>::just_copy_constructible,const T &,T &&>(val); }
    static decltype(auto) _move_if_assignable(T &val){ return ::std::conditional_t<value_traits<>::just_copy_assignable,const T &,T &&>(val); }

    template<typename ...Args>
    pointer _construct(const pointer ptr,Args &&...args){
        alloc_traits::construct(_a.get_first(),ptr,::std::forward<Args>(args)...);
        return ptr;
    }
    void _destroy(pointer ptr){ alloc_traits::destroy(_a.get_first(),ptr); }
    pointer _allocate(const size_type n){
        if constexpr(SizeFixed){
            HTL_ASSERT(!SizeFixed);
            return nullptr;
        }
        else return alloc_traits::allocate(_a.get_first(),n);
    }
    void _deallocate(const pointer ptr,const size_type n){
        alloc_traits::deallocate(_a.get_first(),ptr,n);
    }
    void _move_destroy(pointer pdst,pointer psrc,size_type _size){
        pointer psrcend=psrc+_size;
        while(psrc<psrcend){
            _construct(pdst,_move_if_constructible(*psrc));
            _destroy(psrc);
            ++pdst;++psrc;
        }
    }
    void _reverse(pointer pbegin,pointer pend){
        if constexpr(value_traits<>::is_trivially_swappable)
            ::std::reverse(pbegin,pend);
        else{
            alignas(value_type) unsigned char _tmp_buffer[sizeof(value_type)];
            pointer ptmp=pointer(_tmp_buffer);
            while(pend>pbegin+1){
                _construct(ptmp,_move_if_constructible(*pbegin));
                *pbegin=_move_if_assignable(*--pend);
                *pend=_move_if_assignable(*ptmp);
                ++pbegin;
                _destroy(ptmp);
            }
        }
    }
    void _rotate(pointer pbegin,pointer pmiddle,pointer pend){
        if(pbegin!=pmiddle&&pmiddle!=pend){
            _reverse(pbegin,pmiddle);
            _reverse(pmiddle,pend);
            _reverse(pbegin,pend);
        }
    }
    size_type _calculate_growth(size_type old_capacity,size_type new_size) const{
        constexpr size_type _capacity_init=32/sizeof(value_type);
        size_type new_capacity=old_capacity+(old_capacity>>1)+_capacity_init;
        HTL_ASSERT(old_capacity<=new_capacity&&new_capacity<=max_size());
        return new_capacity<new_size?new_size:new_capacity;
    }
    void _change_array(const size_type _size,const pointer _begin,const size_type _capacity){
        _a.second._size=_size;
        _a.second._alloc=true;
        _u._b=_begin;
        _u._c=_capacity;
    }
    void _init(){
        _a.second._size=0;
        _a.second._alloc=0;
        if constexpr(inplace_size==0){
            _u._b=nullptr;
            _u._c=0;
        }
    }
    pointer _begin() const{
        if constexpr(inplace_size)return _a.second._alloc?_u._b:pointer(_u._buffer);
        else return _u._b;
    }
    //initialize {_size,_alloc,_u._b,_u._c} & return _begin()
    pointer _init_capacity(size_type count){
        pointer pbegin;
        if(count>inplace_size){
            pbegin=_allocate(count);
            _change_array(count,pbegin,count);
        }
        else{
            pbegin=pointer(_u._buffer);
            _a.second._size=count;
            _a.second._alloc=false;
            if constexpr(inplace_size==0){
                _u._b=nullptr;
                _u._c=0;
            }
        }
        return pbegin;
    }
    template<typename I>
    void _init_range_counted(I ibegin,size_type count){
        pointer pbegin=_init_capacity(count);
        if(count>0){
            const pointer pend=pbegin+count;
            do _construct(pbegin,*ibegin);
            while(++ibegin,++pbegin<pend);
        }
    }
    //assume this->get_allocator == (_other.get_allocator(), before move)
    void _init_move(vector &_other){
        if constexpr(::std::is_trivially_destructible_v<value_type>
                   &&::std::is_trivially_move_constructible_v<value_type>
                   ||!inplace_size){
            _a.second=_other._a.second;
            ::std::memcpy(&_u,&_other._u,sizeof(_u));
        }
        else{
            if(_other._a.second._alloc){
                _a.second=_other._a.second;
                ::std::memcpy(&_u,&_other._u,sizeof(_u));
            }
            else{
                const size_type _size=_other.size();
                _move_destroy((pointer)_u._buffer,(pointer)_other._u._buffer,_size);
                _a.second._size=_size;
                _a.second._alloc=false;
            }
        }
        _other._init();
    }
    //caller should ensure if count>1, args are const l-reference
    template<typename ...Args>
    void _resize(size_type new_size,Args &&...args){
        const pointer pbegin=_begin();
        const size_type old_size=size();
        if(new_size<old_size){
            pointer pend=pbegin+old_size;
            const pointer pnewend=pbegin+new_size;
            do _destroy(--pend);
            while(pend>pnewend);
            _a.second._size=new_size;
        }
        else if(new_size>old_size){
            const size_type old_capacity=capacity();
            if(new_size>old_capacity){
                const bool _old_inplace=inplace();
                const size_type new_capacity=_calculate_growth(old_capacity,new_size);
                const pointer pnewbegin=_allocate(new_capacity);
                const pointer pnewend=pnewbegin+new_size;
                pointer pend=pnewbegin+old_size;
                do _construct(pend,::std::forward<Args>(args)...);
                while(++pend<pnewend);
                _move_destroy(pnewbegin,pbegin,old_size);
                if(!_old_inplace)_deallocate(pbegin,old_capacity);
                _change_array(new_size,pnewbegin,new_capacity);
            }
            else{
                const pointer pnewend=pbegin+new_size;
                pointer pend=pbegin+old_size;
                do _construct(pend,::std::forward<Args>(args)...);
                while(++pend<pnewend);
                _a.second._size=new_size;
            }
        }
    }
    //caller should ensure if count>1, args are const l-reference
    template<typename ...Args>
    iterator _insert(const_iterator pos,size_type count,Args &&...args){
        const pointer pbegin=_begin();
        const size_type old_size=size();
        const pointer pend=pbegin+old_size;
        const pointer psrc=pos._ptr;
        HTL_ASSERT(pbegin<=psrc&&psrc<=pend);
        pointer pdst;
        if(count==0)pdst=psrc;
        else{
            const size_type old_capacity=capacity();
            const size_type new_size=old_size+count;
            if(new_size>old_capacity){
                const bool _old_inplace=inplace();
                const size_type new_capacity=_calculate_growth(old_capacity,new_size);
                const pointer pnewbegin=_allocate(new_capacity);
                const size_type _offset=psrc-pbegin;
                pdst=pnewbegin+_offset;
                const pointer pdstend=pdst+count;
                pointer pcur=pdst;
                do _construct(pcur,::std::forward<Args>(args)...);
                while(++pcur<pdstend);
                _move_destroy(pnewbegin,pbegin,_offset);
                _move_destroy(pdstend,psrc,old_size-_offset);
                if(!_old_inplace)_deallocate(pbegin,old_capacity);
                _change_array(new_size,pnewbegin,new_capacity);
            }
            else{
                pdst=psrc;
                const pointer pdstend=pdst+count;
                if(psrc==pend){
                    pointer pcur=pdst;
                    do _construct(pcur,::std::forward<Args>(args)...);
                    while(++pcur<pdstend);
                }
                else{
                    alignas(value_type) unsigned char _tmp_buffer[sizeof(value_type)];
                    pointer ptmp=pointer(_tmp_buffer);
                    _construct(ptmp,::std::forward<Args>(args)...);
                    pointer pcursrc=pend,pcur=pend+count;
                    if(psrc+count>pend){
                        do _construct(--pcur,_move_if_constructible(*--pcursrc));
                        while(pcursrc>psrc);
                        do _construct(--pcur,const_reference(*ptmp));
                        while(pcur>pend);
                    }
                    else{
                        do _construct(--pcur,_move_if_constructible(*--pcursrc));
                        while(pcur>pend);
                        while(pcursrc>psrc)
                            *--pcur=_move_if_assignable(*--pcursrc);
                    }
                    while(--pcur>psrc)
                        *pcur=const_reference(*ptmp);
                    *pdst=_move_if_assignable(*ptmp);
                    _destroy(ptmp);
                }
                _a.second._size=new_size;
            }
        }
        return _make<iterator>(pdst);
    }
    template<typename I>
    iterator _insert_range_counted(const_iterator pos,I ibegin,size_type count){
        const pointer pbegin=_begin();
        const size_type old_size=size();
        const pointer pend=pbegin+old_size;
        const pointer psrc=pos._ptr;
        HTL_ASSERT(pbegin<=psrc&&psrc<=pend);
        pointer pdst;
        if(count==0)pdst=psrc;
        else{
            const size_type old_capacity=capacity();
            const size_type new_size=old_size+count;
            if(new_size>old_capacity){
                const bool _old_inplace=inplace();
                const size_type new_capacity=_calculate_growth(old_capacity,new_size);
                const pointer pnewbegin=_allocate(new_capacity);
                const size_type _offset=psrc-pbegin;
                pdst=pnewbegin+_offset;
                const pointer pdstend=pdst+count;
                pointer pcur=pdst;
                do _construct(pcur,*ibegin);
                while(++ibegin,++pcur<pdstend);
                _move_destroy(pnewbegin,pbegin,_offset);
                _move_destroy(pdstend,psrc,old_size-_offset);
                if(!_old_inplace)_deallocate(pbegin,old_capacity);
                _change_array(new_size,pnewbegin,new_capacity);
            }
            else{
                pdst=psrc;
                pointer pcur=pend;
                const pointer pnewend=pcur+count;
                do _construct(pcur,*ibegin);
                while(++ibegin,++pcur<pnewend);
                _rotate(pdst,pend,pnewend);
                _a.second._size=new_size;
            }
        }
        return _make<iterator>(pdst);
    }
    template<typename I>
    iterator _insert_range_uncounted(const_iterator pos,I ibegin,I iend){
        const pointer pbegin=_begin();
        const size_type old_size=size();
        const pointer pend=pbegin+old_size;
        const pointer psrc=pos._ptr;
        HTL_ASSERT(pbegin<=psrc&&psrc<=pend);
        pointer pdst;
        if(ibegin!=iend){
            const size_type old_capacity=capacity();
            const pointer ptail=pbegin+old_capacity;
            pointer pcur=pend;
            bool remain=true;
            while(pcur<ptail){
                _construct(pcur,*ibegin);
                ++pcur;
                if(!(++ibegin!=iend)){
                    remain=false;
                    break;
                }
            }
            const size_type add_size=pcur-pend;
            if(remain){
                vector _tmp(_a.get_first());
                do _tmp.emplace_back(*ibegin);
                while(++ibegin!=iend);
                const size_type insert_size=add_size+_tmp.size();
                const size_type new_size=old_size+insert_size;
                const bool _old_inplace=inplace();
                const size_type new_capacity=_calculate_growth(old_capacity,new_size);
                const pointer pnewbegin=_allocate(new_capacity);
                const size_type _offset=psrc-pbegin;
                pdst=pnewbegin+_offset;
                _move_destroy(pnewbegin,pbegin,_offset);
                _move_destroy(pdst,pend,add_size);
                _move_destroy(pdst+add_size,_tmp._begin(),insert_size-add_size);
                _move_destroy(pdst+insert_size,psrc,old_size-_offset);
                if(!_old_inplace)_deallocate(pbegin,old_capacity);
                _tmp._a.second._size=0;
                _change_array(new_size,pnewbegin,new_capacity);
            }
            else{
                pdst=psrc;
                _rotate(pdst,pend,pcur);
                _a.second._size=old_size+add_size;
            }
        }
        else pdst=psrc;
        return _make<iterator>(pdst);
    }
    template<bool Move=false,typename I>
    void _assign_range_counted(I ibegin,const size_type new_size){
        static_assert(!Move||::std::is_same_v<I,pointer>);
        const pointer pbegin=_begin();
        const size_type old_size=size();
        const pointer pend=pbegin+old_size;
        const size_type old_capacity=capacity();
        if(new_size>old_capacity){
            const bool _old_inplace=inplace();
            const size_type new_capacity=_calculate_growth(old_capacity,new_size);
            const pointer pnewbegin=_allocate(new_capacity);
            const pointer pnewend=pnewbegin+new_size;
            pointer pcur=pnewbegin;
            do{
                if constexpr(Move)
                    _construct(pcur,_move_if_constructible(*ibegin));
                else
                    _construct(pcur,*ibegin);
                ++ibegin;
            } while(++pcur<pnewend);
            pcur=pend;
            while(pbegin<pcur)
                _destroy(--pcur);
            if(!_old_inplace)_deallocate(pbegin,old_capacity);
            _change_array(new_size,pnewbegin,new_capacity);
        }
        else{
            pointer pcur=pbegin;
            const pointer pcurend=pbegin+(::std::min)(new_size,old_size);
            while(pcur<pcurend){
                if constexpr(Move)
                    *pcur=_move_if_assignable(*ibegin);
                else
                    *pcur=*ibegin;
                ++ibegin;
                ++pcur;
            }
            if(new_size<old_size){
                do _destroy(pcur);
                while(++pcur<pend);
                _a.second._size=new_size;
            }
            else if(new_size>old_size){
                const pointer pnewend=pbegin+new_size;
                do{
                    if constexpr(Move)
                        _construct(pcur,_move_if_constructible(*ibegin));
                    else
                        _construct(pcur,*ibegin);
                    ++ibegin;
                } while(++pcur<pnewend);
                _a.second._size=new_size;
            }
        }
    }
    template<typename I>
    void _assign_range_uncounted(I ibegin,I iend){
        const pointer pbegin=_begin();
        const size_type old_size=size();
        const pointer pend=pbegin+old_size;
        pointer pcur=pbegin;
        while(ibegin!=iend){
            bool remain=true;
            while(pcur<pend){
                *pcur=*ibegin;
                ++pcur;
                if(!(++ibegin!=iend)){
                    remain=false;
                    break;
                }
            }
            if(!remain)break;
            const size_type old_capacity=capacity();
            pointer ptail=pbegin+old_capacity;
            while(pcur<ptail){
                _construct(pcur,*ibegin);
                ++pcur;
                if(!(++ibegin!=iend)){
                    remain=false;
                    break;
                }
            }
            if(!remain)break;
            vector _tmp(_a.get_first());
            do _tmp.emplace_back(*ibegin);
            while(++ibegin!=iend);
            const size_type new_size=old_capacity+_tmp.size();
            const bool _old_inplace=inplace();
            const size_type new_capacity=_calculate_growth(old_capacity,new_size);
            const pointer pnewbegin=_allocate(new_capacity);
            _move_destroy(pnewbegin,pbegin,old_capacity);
            _move_destroy(pnewbegin+old_capacity,_tmp._begin(),new_size-old_capacity);
            if(!_old_inplace)_deallocate(pbegin,old_capacity);
            _tmp._a.second._size=0;
            _change_array(new_size,pnewbegin,new_capacity);
            return;
        }
        pointer pexend=pend;
        while(pcur<pexend)
            _destroy(--pexend);
        _a.second._size=pcur-pbegin;
    }
    template<typename I>
    auto _make(pointer ptr) const{
        I it;
        it._ptr=ptr;
        return it;
    }
    pointer _dereference_index(size_type index) const{
        HTL_ASSERT(index<size());
        return _begin()+index;
    }
public:
    vector(){ _init(); }
    explicit vector(const allocator_type &alloc):_a(alloc){ _init(); }
    explicit vector(size_type _size,const allocator_type &alloc=allocator_type()):_a(alloc){
        _init();
        _resize(_size);
    }
    vector(size_type _size,const value_type &_val,const allocator_type &alloc=allocator_type()):_a(alloc){
        _init();
        _resize(_size,_val);
    }
    template<typename I,typename C=typename ::std::iterator_traits<I>::iterator_category>
    vector(I ibegin,I iend,const allocator_type &alloc=allocator_type()):_a(alloc){
        if constexpr(::std::is_convertible_v<C,::std::forward_iterator_tag>)
            _init_range_counted(ibegin,::std::distance(ibegin,iend));
        else{
            _init();
            for(;ibegin!=iend;++ibegin)
                emplace_back(*ibegin);
        }
    }
    vector(const vector &_other):_a(alloc_traits::select_on_container_copy_construction(_other.get_allocator())){
        _init_range_counted(_other._begin(),_other.size());
    }
    vector(const vector &_other,const allocator_type &alloc):_a(alloc){
        _init_range_counted(_other._begin(),_other.size());
    }
    vector(vector &&_other) noexcept(value_traits<>::nothrow_move_constructible)
        :_a(::std::move(_other.get_allocator())){
        _init_move(_other);
    }
    vector(vector &&_other,const allocator_type &alloc) noexcept(value_traits<>::nothrow_move_constructible_with_allocator)
        :_a(alloc){
        if constexpr(alloc_traits::is_always_equal::value)
            _init_move(_other);
        else{
            if(_other._a.get_first()==_a.get_first())
                _init_move(_other);
            else{
                const size_type count=_other.size();
                pointer pbegin=_init_capacity(count);
                if(count>0){
                    const pointer pend=pbegin+count;
                    const pointer pother=_other._begin();
                    do{
                        _construct(pbegin,_move_if_constructible(*pother));
                        _other._destroy(pother);
                    } while(++pother,++pbegin<pend);
                    _other._a.second._size=0;
                }
            }
        }
    }
    vector(::std::initializer_list<value_type> _array,const allocator_type &alloc=allocator_type()):_a(alloc){
        _init_range_counted(_array.begin(),_array.size());
    }
    ~vector(){ free(); }

    vector &operator=(const vector &_other){
        if(this!=&_other){
            if constexpr(alloc_traits::propagate_on_container_copy_assignment::value){
                allocator_type &this_alloc=_a.get_first();
                allocator_type &other_alloc=_other._a.get_first();
                if constexpr(!alloc_traits::is_always_equal::value){
                    if(this_alloc!=other_alloc)
                        free();
                }
                this_alloc=other_alloc;
            }
            _assign_range_counted(_other._begin(),_other.size());
        }
        return *this;
    }
    vector &operator=(vector &&_other) noexcept(value_traits<>::nothrow_move_assignable){
        if(this!=&_other){
            allocator_type &this_alloc=_a.get_first();
            allocator_type &other_alloc=_other._a.get_first();
            if constexpr(!alloc_traits::propagate_on_container_move_assignment::value){
                if constexpr(!alloc_traits::is_always_equal::value){
                    if(this_alloc!=other_alloc){
                        _assign_range_counted<true>(_other._begin(),_other.size());
                        _other.clear();
                        return *this;
                    }
                }
                free();
            }
            else{
                free();
                this_alloc=::std::move(other_alloc);
            }
            _init_move(_other);
        }
        return *this;
    }
    vector &operator=(::std::initializer_list<value_type> _array){
        _assign_range_counted(_array.begin(),_array.size());
        return *this;
    }
    void assign(size_type count,const value_type &val){
        class dummy_iterator{
            const value_type &ref;
        public:
            dummy_iterator(const value_type &_r):ref(_r){}
            const value_type &operator *() const{ return ref; }
            void operator++() const{}
        };
        _assign_range_counted(dummy_iterator(val),count);
    }
    template<typename I,typename C=typename ::std::iterator_traits<I>::iterator_category>
    void assign(I ibegin,I iend){
        if constexpr(::std::is_convertible_v<C,::std::forward_iterator_tag>)
            _assign_range_counted(ibegin,::std::distance(ibegin,iend));
        else
            _assign_range_uncounted(ibegin,iend);
    }
    void assign(::std::initializer_list<value_type> _array){
        _assign_range_counted(_array.begin(),_array.size());
    }
    void swap(vector &_other) noexcept(value_traits<>::nothrow_swappable){
        if(this==&_other)return;

        allocator_type &this_alloc=_a.get_first();
        allocator_type &other_alloc=_other._a.get_first();
        auto _swap_trivially=[&](){
            if constexpr(alloc_traits::propagate_on_container_swap::value)
                swap(this_alloc,other_alloc);
            ::std::swap(_a.second,_other._a.second);
            ::std::swap(_u,_other._u);
        };

        if constexpr(!inplace_size||value_traits<>::is_trivially_swappable)
            _swap_trivially();
        else{
            size_type this_cost=inplace()?size():0;
            size_type other_cost=_other.inplace()?_other.size():0;
            if(this_cost||other_cost){
                if constexpr(!alloc_traits::is_always_equal::value){
                    if(this_alloc!=other_alloc){
                        if constexpr(alloc_traits::propagate_on_container_swap::value){
                            vector _tmp(::std::move(_other));
                            other_alloc=::std::move(this_alloc);
                            _other._init_move(*this);
                            this_alloc=::std::move(_tmp._a.get_first());
                            _init_move(_tmp);
                        }
                        else{
                            // according to C++ standard, it is undefined behavior 
                            // calling swap on vectors with unequal allocators when propagate_on_container_swap is false.
                            // however, here defined as is non-specialized templated ::std::swap:
                            vector _tmp(::std::move(_other));
                            _other=::std::move(*this);
                            *this=::std::move(_tmp);
                        }
                        return;
                    }
                }
                const bool keep_this=this_cost>other_cost;
                vector &heavy=keep_this?*this:_other;
                vector &light=keep_this?_other:*this;
                vector _tmp(this_alloc);
                _tmp._init_move(light);
                light._init_move(heavy);
                heavy._init_move(_tmp);
            }
            else
                _swap_trivially();
        }
    }

    allocator_type get_allocator() const{ return _a.get_first(); }
    size_type size() const{ return _a.second._size; }
    size_type max_size() const{
        return ::std::min<size_type>(alloc_traits::max_size(_a.get_first()),(::std::numeric_limits<difference_type>::max)());
    }
    size_type capacity() const{
        if constexpr(inplace_size)return _a.second._alloc?_u._c:inplace_size;
        else return _u._c;
    }
    bool empty() const{ return size()==0; }
    bool inplace() const{
        if constexpr(inplace_size)return !_a.second._alloc;
        else return false;
    }
    size_type inplace_capacity() const{ return inplace_size; }

    void resize(size_type new_size){ _resize(new_size); }
    void resize(size_type new_size,const value_type &val){ _resize(new_size,val); }
    void clear(){
        const pointer pbegin=_begin();
        pointer pend=pbegin+size();
        while(pbegin<pend)
            _destroy(--pend);
        _a.second._size=0;
    }
    void free(){
        const pointer pbegin=_begin();
        pointer pend=pbegin+size();
        const bool _old_inplace=inplace();
        while(pbegin<pend)
            _destroy(--pend);
        if(!_old_inplace)_deallocate(pbegin,_u._c);
        _init();
    }
    void shrink_to_fit(){
        if(inplace())
            return;
        const size_type _size=size();
        const size_type old_capacity=_u._c;
        if(_size==old_capacity)
            return;

        const pointer pbegin=_u._b;
        const pointer pnewbegin=_init_capacity(_size);
        _move_destroy(pnewbegin,pbegin,_size);
        _deallocate(pbegin,old_capacity);
    }
    void reserve(size_type new_capacity){
        size_type old_capacity=capacity();
        if(new_capacity>old_capacity){
            const size_type _size=size();
            const bool _old_inplace=inplace();
            const pointer pbegin=_begin();
            const pointer pnewbegin=_allocate(new_capacity);
            _move_destroy(pnewbegin,pbegin,_size);
            if(!_old_inplace)_deallocate(pbegin,old_capacity);
            _change_array(_size,pnewbegin,new_capacity);
        }
    }

    iterator begin(){ return _make<iterator>(_begin()); }
    iterator end(){ return _make<iterator>(_begin()+size()); }
    const_iterator begin() const{ return _make<const_iterator>(_begin()); }
    const_iterator end() const{ return _make<const_iterator>(_begin()+size()); }
    const_iterator cbegin() const{ return _make<const_iterator>(_begin()); }
    const_iterator cend() const{ return _make<const_iterator>(_begin()+size()); }
    reverse_iterator rbegin(){ return reverse_iterator(end()); }
    reverse_iterator rend(){ return reverse_iterator(begin()); }
    const_reverse_iterator rbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator rend() const{ return const_reverse_iterator(begin()); }
    const_reverse_iterator crbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator crend() const{ return const_reverse_iterator(begin()); }

    pointer data(){ return _begin(); }
    const_pointer data() const{ return _begin(); }
    reference front(){ return *_dereference_index(0); }
    reference back(){ return *_dereference_index(size()-1); }
    const_reference front() const{ return *_dereference_index(0); }
    const_reference back() const{ return *_dereference_index(size()-1); }
    reference at(size_type index){
        HTL_ASSERT(index<size());
        return *(_begin()+index);
    }
    reference operator[](size_type index){ return *_dereference_index(index); }
    const_reference at(size_type index) const{
        HTL_ASSERT(index<size());
        return *(_begin()+index);
    }
    const_reference operator[](size_type index) const{ return *_dereference_index(index); }

    template<int I2,bool S2,typename A2>
    bool operator==(const vector<T,I2,S2,A2> &_other) const{
        const size_type _size=size();
        if(_size!=_other.size())return false;
        const const_pointer pbegin=_begin();
        return ::std::equal(pbegin,pbegin+_size,_other._begin());
    }
    template<int I2,bool S2,typename A2>
    bool operator!=(const vector<T,I2,S2,A2> &_other) const{ return !(*this==_other); }
    template<int I2,bool S2,typename A2>
    bool operator< (const vector<T,I2,S2,A2> &_other) const{
        return ::std::lexicographical_compare(begin(),end(),_other.begin(),_other.end());
    }
    template<int I2,bool S2,typename A2>
    bool operator<=(const vector<T,I2,S2,A2> &_other) const{ return !(_other<*this); }
    template<int I2,bool S2,typename A2>
    bool operator> (const vector<T,I2,S2,A2> &_other) const{ return _other<*this; }
    template<int I2,bool S2,typename A2>
    bool operator>=(const vector<T,I2,S2,A2> &_other) const{ return !(*this<_other); }

    template<typename ...Args>
    reference emplace_back(Args &&...args){
        const pointer pbegin=_begin();
        const size_type old_size=size();
        const size_type old_capacity=capacity();
        const size_type new_size=old_size+1;
        pointer ptr;
        if(new_size>old_capacity){
            const bool _old_inplace=inplace();
            const size_type new_capacity=_calculate_growth(old_capacity,new_size);
            const pointer pnewbegin=_allocate(new_capacity);
            ptr=_construct(pnewbegin+old_size,::std::forward<Args>(args)...);
            _move_destroy(pnewbegin,pbegin,old_size);
            if(!_old_inplace)_deallocate(pbegin,old_capacity);
            _change_array(new_size,pnewbegin,new_capacity);
        }
        else{
            ptr=_construct(pbegin+old_size,::std::forward<Args>(args)...);
            _a.second._size=new_size;
        }
        return *ptr;
    }
    void push_back(const value_type &val){ emplace_back(val); }
    void push_back(value_type &&val){ emplace_back(::std::move(val)); }
    void pop_back(){
        size_type _size=size();
        HTL_ASSERT(_size>0);
        const pointer pbegin=_begin();
        _destroy(pbegin+--_size);
        _a.second._size=_size;
    }

    template<typename ...Args>
    iterator emplace(const_iterator pos,Args &&...args){
        return _insert(pos,1,::std::forward<Args>(args)...);
    }
    iterator insert(const_iterator pos,const value_type &val){
        return _insert(pos,1,val);
    }
    iterator insert(const_iterator pos,value_type &&val){
        return _insert(pos,1,::std::move(val));
    }
    iterator insert(const_iterator pos,size_type count,const value_type &val){
        return _insert(pos,count,val);
    }
    template<typename I,typename C=typename ::std::iterator_traits<I>::iterator_category>
    iterator insert(const_iterator pos,I ibegin,I iend){
        if constexpr(::std::is_convertible_v<C,::std::forward_iterator_tag>)
            return _insert_range_counted(pos,ibegin,::std::distance(ibegin,iend));
        else
            return _insert_range_uncounted(pos,ibegin,iend);
    }
    iterator insert(const_iterator pos,::std::initializer_list<value_type> _array){
        insert(pos,_array.begin(),_array.end());
    }

    iterator erase(const_iterator pos){
        const pointer pbegin=_begin();
        const size_type _size=size();
        const pointer pend=pbegin+_size;
        const pointer perase=pos._ptr;
        HTL_ASSERT(pbegin<=perase&&perase<pend);
        pointer pcur=perase;
        while(++pcur<pend)
            *(pcur-1)=_move_if_assignable(*pcur);
        _destroy(pcur-1);
        _a.second._size=_size-1;
        return _make<iterator>(perase);
    }
    iterator erase(const_iterator ibegin,const_iterator iend){
        const pointer pbegin=_begin();
        size_type _size=size();
        const pointer pend=pbegin+_size;
        const pointer perase=ibegin._ptr;
        pointer pcursrc=iend._ptr;
        HTL_ASSERT(pbegin<=perase&&perase<=pcursrc&&pcursrc<=pend);
        pointer pcurdst=perase;
        while(pcursrc<pend){
            *pcurdst=_move_if_assignable(*pcursrc);
            ++pcursrc;
            ++pcurdst;
        }
        while(pcurdst<pend){
            _destroy(pcurdst);
            ++pcurdst;
            --_size;
        }
        _a.second._size=_size;
        return _make<iterator>(perase);
    }

};

template<typename T,int N>
using inplace_vector=vector<T,N,true>;

template<typename T,int I,bool S,typename A>
void swap(vector<T,I,S,A> &_lhs,vector<T,I,S,A> &_rhs){
    _lhs.swap(_rhs);
}

}//namespace htl
