#pragma once
#include"rbtree.h"

namespace htl{

template<typename K,typename ...Args>
struct _map_key_from_args{ static constexpr bool extractable=false; };
template<typename K,typename V>
struct _map_key_from_args<K,K,V>{
    static constexpr bool extractable=true;
    static const K &get_key(const K &k,const V &v){ return k; }
};
template<typename K,typename U,typename V>
struct _map_key_from_args<K,::std::pair<U,V>>{
    static constexpr bool extractable=::std::is_same_v<K,remove_cvref_t<U>>;
    static const K &get_key(const ::std::pair<U,V> &pair){ return pair.first; }
};

template<typename NodeHandle,typename K,typename V>
class _map_node_handle_base{
    auto &_value() const{ return static_cast<const NodeHandle &>(*this)._ptr->_d; }
public:
    typedef K key_type;
    typedef V mapped_type;

    key_type &key() const{ return const_cast<key_type &>(_value().first); }
    mapped_type &mapped() const{ return _value().second; }
};

template<bool M,bool L,typename K,typename V,typename C,typename A>
struct _map_traits{
    typedef K key_type;
    typedef ::std::pair<const K,V> value_type;
    typedef V mapped_type;
    typedef C key_compare;
    typedef A allocator_type;
    static constexpr bool Multi=M,Linked=L,Set=false;
    typedef ::std::conditional_t<Linked,_linked_rbnode,_rbnode> node_base;
    typedef _rbnode_data<node_base,value_type> node_data;
    typedef _rbnode_handle<node_data,allocator_type,_map_node_handle_base,key_type,mapped_type> node_type;

    static const key_type &get_key(const value_type &data){ return data.first; }

    template<typename ...Args>
    using get_key_from_args=_map_key_from_args<K,remove_cvref_t<Args>...>;

    class value_compare{
    protected:
        friend class _rbtree<_map_traits>;
        key_compare comp;
        value_compare(key_compare _comp):comp(_comp){}
    public:
        typedef value_type first_argument_type;
        typedef value_type second_argument_type;
        typedef bool result_type;

        bool operator()(const value_type &_left,const value_type &_right) const{
            return comp(get_key(_left),get_key(_right));
        }
    };
};

template<typename K,typename V,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<::std::pair<const K,V>>>
class linked_map;

template<typename K,typename V,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<::std::pair<const K,V>>>
class linked_multimap;

template<typename K,typename V,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<::std::pair<const K,V>>,bool Linked=false>
class map:public _rbtree<_map_traits<false,Linked,K,V,Comparer,Allocator>>{
protected:
    typedef _map_traits<false,Linked,K,V,Comparer,Allocator> TreeTraits;
    typedef _rbtree<TreeTraits> rbtree_type;
public:
    using typename rbtree_type::              key_type;
    using typename rbtree_type::            value_type;
    using typename rbtree_type::           key_compare;
    using typename rbtree_type::        allocator_type;
    using typename rbtree_type::         value_compare;

    using typename rbtree_type::             size_type;
    using typename rbtree_type::       difference_type;
    using typename rbtree_type::               pointer;
    using typename rbtree_type::         const_pointer;
    using typename rbtree_type::             reference;
    using typename rbtree_type::       const_reference;

    using typename rbtree_type::              iterator;
    using typename rbtree_type::        const_iterator;
    using typename rbtree_type::      reverse_iterator;
    using typename rbtree_type::const_reverse_iterator;

    using typename rbtree_type::node_type;
    typedef _rbinsert_return_type<iterator,node_type> insert_return_type;
    typedef typename TreeTraits::mapped_type mapped_type;

    using rbtree_type::rbtree_type;
    using rbtree_type::operator=;

    explicit map(const map<K,V,Comparer,Allocator,!Linked> &_other)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other)){
    }
    map(const map<K,V,Comparer,Allocator,!Linked> &_other,const allocator_type &alloc)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other),alloc){
    }
    map &operator=(const map<K,V,Comparer,Allocator,!Linked> &_other){
        rbtree_type::template _copy_base<false>(_other);
        rbtree_type::insert(_other.begin(),_other.end());
        return *this;
    }
    explicit map(const linked_map<K,V,Comparer,Allocator> &)=delete;
    map(const linked_map<K,V,Comparer,Allocator> &,const allocator_type &)=delete;
    map &operator=(const linked_map<K,V,Comparer,Allocator> &)=delete;

    mapped_type &at(const key_type &key){
        return rbtree_type::_at(key)->second;
    }
    const mapped_type &at(const key_type &key) const{
        return rbtree_type::_at(key)->second;
    }

    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator==(const C &_other) const{
        return ::std::equal(rbtree_type::begin(),rbtree_type::end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator!=(const C &_other) const{ return !(*this==_other); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator< (const C &_other) const{
        return ::std::lexicographical_compare(rbtree_type::begin(),rbtree_type::end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator<=(const C &_other) const{ return !(_other<*this); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator> (const C &_other) const{ return _other<*this; }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator>=(const C &_other) const{ return !(*this<_other); }

    using rbtree_type::try_emplace;
    mapped_type &operator[](const key_type &key){
        return rbtree_type::try_emplace(key).first->second;
    }
    mapped_type &operator[](key_type &&key){
        return rbtree_type::try_emplace(::std::move(key)).first->second;
    }
    using rbtree_type::insert;
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    ::std::pair<iterator,bool> insert(V2 &&v){
        return rbtree_type::emplace(::std::forward<V2>(v));
    }
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    iterator insert(const_iterator hint,V2 &&v){
        return rbtree_type::emplace_hint(hint,::std::forward<V2>(v));
    }
    template<typename M>
    ::std::pair<iterator,bool> insert_or_assign(const key_type &key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::try_emplace(key,::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result;
    }
    template<typename M>
    ::std::pair<iterator,bool> insert_or_assign(key_type &&key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::try_emplace(::std::move(key),::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result;
    }
    template<typename M>
    iterator insert_or_assign(const_iterator hint,const key_type &key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::_try_emplace_hint(hint,key,::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result.first;
    }
    template<typename M>
    iterator insert_or_assign(const_iterator hint,key_type &&key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::_try_emplace_hint(hint,::std::move(key),::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result.first;
    }
};

template<typename K,typename V,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<::std::pair<const K,V>>,bool Linked=false>
class multimap:public _rbtree<_map_traits<true,Linked,K,V,Comparer,Allocator>>{
protected:
    typedef _map_traits<true,Linked,K,V,Comparer,Allocator> TreeTraits;
    typedef _rbtree<TreeTraits> rbtree_type;
public:
    using typename rbtree_type::              key_type;
    using typename rbtree_type::            value_type;
    using typename rbtree_type::           key_compare;
    using typename rbtree_type::        allocator_type;
    using typename rbtree_type::         value_compare;

    using typename rbtree_type::             size_type;
    using typename rbtree_type::       difference_type;
    using typename rbtree_type::               pointer;
    using typename rbtree_type::         const_pointer;
    using typename rbtree_type::             reference;
    using typename rbtree_type::       const_reference;

    using typename rbtree_type::              iterator;
    using typename rbtree_type::        const_iterator;
    using typename rbtree_type::      reverse_iterator;
    using typename rbtree_type::const_reverse_iterator;

    using typename rbtree_type::node_type;
    typedef typename TreeTraits::mapped_type mapped_type;

    using rbtree_type::rbtree_type;
    using rbtree_type::operator=;

    explicit multimap(const multimap<K,V,Comparer,Allocator,!Linked> &_other)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other)){
    }
    multimap(const multimap<K,V,Comparer,Allocator,!Linked> &_other,const allocator_type &alloc)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other),alloc){
    }
    multimap &operator=(const multimap<K,V,Comparer,Allocator,!Linked> &_other){
        rbtree_type::template _copy_base<false>(_other);
        rbtree_type::insert(_other.begin(),_other.end());
        return *this;
    }
    explicit multimap(const linked_multimap<K,V,Comparer,Allocator> &)=delete;
    multimap(const linked_multimap<K,V,Comparer,Allocator> &,const allocator_type &)=delete;
    multimap &operator=(const linked_multimap<K,V,Comparer,Allocator> &)=delete;

    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator==(const C &_other) const{
        return ::std::equal(rbtree_type::begin(),rbtree_type::end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator!=(const C &_other) const{ return !(*this==_other); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator< (const C &_other) const{
        return ::std::lexicographical_compare(rbtree_type::begin(),rbtree_type::end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator<=(const C &_other) const{ return !(_other<*this); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator> (const C &_other) const{ return _other<*this; }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator>=(const C &_other) const{ return !(*this<_other); }

    using rbtree_type::insert;
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    iterator insert(V2 &&v){
        return rbtree_type::emplace(::std::forward<V2>(v));
    }
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    iterator insert(const_iterator hint,V2 &&v){
        return rbtree_type::emplace_hint(hint,::std::forward<V2>(v));
    }
};

template<typename K,typename V,typename Comparer,typename Allocator>
class linked_map:public map<K,V,Comparer,Allocator,true>{
public:
    typedef map<K,V,Comparer,Allocator,true> unlinked_type;
protected:
    using typename unlinked_type::TreeTraits;
    typedef _rbtree<TreeTraits> rbtree_type;
    using rbtree_type::lower_bound;
    using rbtree_type::upper_bound;
    using rbtree_type::equal_range;
    using rbtree_type::emplace_hint;
public:
    using typename rbtree_type::       key_type;
    using typename rbtree_type::     value_type;
    using typename rbtree_type::    key_compare;
    using typename rbtree_type:: allocator_type;
    using typename rbtree_type::  value_compare;

    using typename rbtree_type::      size_type;
    using typename rbtree_type::difference_type;
    using typename rbtree_type::        pointer;
    using typename rbtree_type::  const_pointer;
    using typename rbtree_type::      reference;
    using typename rbtree_type::const_reference;

    typedef typename rbtree_type::               linked_iterator               iterator;
    typedef typename rbtree_type::         const_linked_iterator         const_iterator;
    typedef typename rbtree_type::       reverse_linked_iterator       reverse_iterator;
    typedef typename rbtree_type:: const_reverse_linked_iterator const_reverse_iterator;

    using typename rbtree_type::node_type;
    typedef _rbinsert_return_type<iterator,node_type> insert_return_type;
    typedef typename TreeTraits::mapped_type mapped_type;

    const unlinked_type &as_unlinked() const{ return *this; }
    unlinked_type &as_unlinked(){ return *this; }

    using unlinked_type::unlinked_type;
    using rbtree_type::operator=;
    using rbtree_type::sort;

    linked_map(const linked_map &_other,const allocator_type &alloc):unlinked_type(_other.as_unlinked(),alloc){}
    linked_map(linked_map &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other.as_unlinked()),alloc){}
    explicit linked_map(const unlinked_type &_other):unlinked_type(_other){
        rbtree_type::sort();
    }
    linked_map(const unlinked_type &_other,const allocator_type &alloc):unlinked_type(_other,alloc){
        rbtree_type::sort();
    }
    linked_map &operator=(const unlinked_type &_other){
        if(this!=::std::addressof(_other)){
            rbtree_type::template _copy_base<>(_other);
            rbtree_type::template _copy_init<false,false>(_other);
        }
        rbtree_type::sort();
        return *this;
    }
    explicit linked_map(unlinked_type &&_other):unlinked_type(::std::move(_other)){
        rbtree_type::sort();
    }
    linked_map(unlinked_type &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other),alloc){
        rbtree_type::sort();
    }
    linked_map &operator=(unlinked_type &&_other){
        rbtree_type::operator=(::std::move(_other));
        rbtree_type::sort();
        return *this;
    }
    explicit linked_map(const map<K,V,Comparer,Allocator> &)=delete;
    explicit linked_map(const map<K,V,Comparer,Allocator> &,const allocator_type &)=delete;

    iterator begin(){ return rbtree_type::template _make<iterator>(rbtree_type::_lbegin()); }
    iterator end(){ return rbtree_type::template _make<iterator>(rbtree_type::_end()); }
    const_iterator begin() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_lbegin()); }
    const_iterator end() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_end()); }
    const_iterator cbegin() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_lbegin()); }
    const_iterator cend() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_end()); }
    reverse_iterator rbegin(){ return reverse_iterator(end()); }
    reverse_iterator rend(){ return reverse_iterator(begin()); }
    const_reverse_iterator rbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator rend() const{ return const_reverse_iterator(begin()); }
    const_reverse_iterator crbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator crend() const{ return const_reverse_iterator(begin()); }

    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator==(const C &_other) const{
        return ::std::equal(begin(),end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator!=(const C &_other) const{ return !(*this==_other); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator< (const C &_other) const{
        return ::std::lexicographical_compare(begin(),end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator<=(const C &_other) const{ return !(_other<*this); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator> (const C &_other) const{ return _other<*this; }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator>=(const C &_other) const{ return !(*this<_other); }

    //try emplace directly before pos
    template<typename ...Args>
    ::std::pair<iterator,bool> try_emplace(const_iterator pos,const key_type &key,Args &&...args){
        return rbtree_type::_ltry_emplace(pos,key,::std::forward<Args>(args)...);
    }
    //try emplace directly before pos
    template<typename ...Args>
    ::std::pair<iterator,bool> try_emplace(const_iterator pos,key_type &&key,Args &&...args){
        return rbtree_type::_ltry_emplace(pos,::std::move(key),::std::forward<Args>(args)...);
    }
    //try emplace directly before end()
    template<typename ...Args>
    ::std::pair<iterator,bool> try_emplace(const key_type &key,Args &&...args){
        return rbtree_type::_ltry_emplace(cend(),key,::std::forward<Args>(args)...);
    }
    //try emplace directly before end()
    template<typename ...Args>
    ::std::pair<iterator,bool> try_emplace(key_type &&key,Args &&...args){
        return rbtree_type::_ltry_emplace(cend(),::std::move(key),::std::forward<Args>(args)...);
    }
    //emplace directly before pos
    template<typename ...Args>
    ::std::pair<iterator,bool> emplace(const_iterator pos,Args &&...args){
        return rbtree_type::_lemplace(pos,::std::forward<Args>(args)...);
    }
    //emplace directly before end()
    template<typename ...Args>
    ::std::pair<iterator,bool> emplace(Args &&...args){
        return rbtree_type::_lemplace(cend(),::std::forward<Args>(args)...);
    }
    mapped_type &operator[](const key_type &key){
        return rbtree_type::_ltry_emplace(cend(),key).first->second;
    }
    mapped_type &operator[](key_type &&key){
        return rbtree_type::_ltry_emplace(cend(),::std::move(key)).first->second;
    }

    ::std::pair<iterator,bool> insert(const value_type &v){
        return rbtree_type::_lemplace(cend(),v);
    }
    ::std::pair<iterator,bool> insert(value_type &&v){
        return rbtree_type::_lemplace(cend(),::std::move(v));
    }
    ::std::pair<iterator,bool> insert(const_iterator pos,const value_type &v){
        return rbtree_type::_lemplace(pos,v);
    }
    ::std::pair<iterator,bool> insert(const_iterator pos,value_type &&v){
        return rbtree_type::_lemplace(pos,::std::move(v));
    }
    insert_return_type insert(node_type &&n){
        return rbtree_type::_linsert_node(cend(),n);
    }
    insert_return_type insert(const_iterator pos,node_type &&n){
        return rbtree_type::_linsert_node(pos,n);
    }
    template<typename I>
    void insert(I ibegin,I iend){
        rbtree_type::insert(ibegin,iend);
    }
    void insert(::std::initializer_list<value_type> ilist){
        rbtree_type::insert(ilist.begin(),ilist.end());
    }
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    ::std::pair<iterator,bool> insert(V2 &&v){
        return rbtree_type::_lemplace(cend(),::std::forward<V2>(v));
    }
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    ::std::pair<iterator,bool> insert(const_iterator pos,V2 &&v){
        return rbtree_type::_lemplace(pos,::std::forward<V2>(v));
    }
    template<typename M>
    ::std::pair<iterator,bool> insert_or_assign(const key_type &key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::_ltry_emplace(cend(),key,::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result;
    }
    template<typename M>
    ::std::pair<iterator,bool> insert_or_assign(key_type &&key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::_ltry_emplace(cend(),::std::move(key),::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result;
    }
    template<typename M>
    ::std::pair<iterator,bool> insert_or_assign(const_iterator pos,const key_type &key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::_ltry_emplace(pos,key,::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result;
    }
    template<typename M>
    ::std::pair<iterator,bool> insert_or_assign(const_iterator pos,key_type &&key,M &&mapped){
        ::std::pair<iterator,bool> result=rbtree_type::_ltry_emplace(pos,::std::move(key),::std::forward<M>(mapped));
        if(!result.second)result.first->second=::std::forward<M>(mapped);
        return result;
    }

    iterator find(const key_type &key){
        return rbtree_type::template _make<iterator>(rbtree_type::_find(key));
    }
    const_iterator find(const key_type &key) const{
        return rbtree_type::template _make<const_iterator>(rbtree_type::_find(key));
    }
};

template<typename K,typename V,typename Comparer,typename Allocator>
class linked_multimap:public multimap<K,V,Comparer,Allocator,true>{
public:
    typedef multimap<K,V,Comparer,Allocator,true> unlinked_type;
protected:
    using typename unlinked_type::TreeTraits;
    typedef _rbtree<TreeTraits> rbtree_type;
    using rbtree_type::lower_bound;
    using rbtree_type::upper_bound;
    using rbtree_type::equal_range;
    using rbtree_type::emplace_hint;
public:
    using typename rbtree_type::       key_type;
    using typename rbtree_type::     value_type;
    using typename rbtree_type::    key_compare;
    using typename rbtree_type:: allocator_type;
    using typename rbtree_type::  value_compare;

    using typename rbtree_type::      size_type;
    using typename rbtree_type::difference_type;
    using typename rbtree_type::        pointer;
    using typename rbtree_type::  const_pointer;
    using typename rbtree_type::      reference;
    using typename rbtree_type::const_reference;

    typedef typename rbtree_type::               linked_iterator               iterator;
    typedef typename rbtree_type::         const_linked_iterator         const_iterator;
    typedef typename rbtree_type::       reverse_linked_iterator       reverse_iterator;
    typedef typename rbtree_type:: const_reverse_linked_iterator const_reverse_iterator;

    using typename rbtree_type::node_type;
    typedef typename TreeTraits::mapped_type mapped_type;

    const unlinked_type &as_unlinked() const{ return *this; }
    unlinked_type &as_unlinked(){ return *this; }

    using unlinked_type::unlinked_type;
    using rbtree_type::operator=;
    using rbtree_type::sort;

    linked_multimap(const linked_multimap &_other,const allocator_type &alloc):unlinked_type(_other.as_unlinked(),alloc){}
    linked_multimap(linked_multimap &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other.as_unlinked()),alloc){}
    explicit linked_multimap(const unlinked_type &_other):unlinked_type(_other){
        rbtree_type::sort();
    }
    linked_multimap(const unlinked_type &_other,const allocator_type &alloc):unlinked_type(_other,alloc){
        rbtree_type::sort();
    }
    linked_multimap &operator=(const unlinked_type &_other){
        if(this!=::std::addressof(_other)){
            rbtree_type::template _copy_base<>(_other);
            rbtree_type::template _copy_init<false,false>(_other);
        }
        rbtree_type::sort();
        return *this;
    }
    explicit linked_multimap(unlinked_type &&_other):unlinked_type(::std::move(_other)){
        rbtree_type::sort();
    }
    linked_multimap(unlinked_type &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other),alloc){
        rbtree_type::sort();
    }
    linked_multimap &operator=(unlinked_type &&_other){
        rbtree_type::operator=(::std::move(_other));
        rbtree_type::sort();
        return *this;
    }
    explicit linked_multimap(const multimap<K,V,Comparer,Allocator> &)=delete;
    explicit linked_multimap(const multimap<K,V,Comparer,Allocator> &,const allocator_type &)=delete;

    iterator begin(){ return rbtree_type::template _make<iterator>(rbtree_type::_lbegin()); }
    iterator end(){ return rbtree_type::template _make<iterator>(rbtree_type::_end()); }
    const_iterator begin() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_lbegin()); }
    const_iterator end() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_end()); }
    const_iterator cbegin() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_lbegin()); }
    const_iterator cend() const{ return rbtree_type::template _make<const_iterator>(rbtree_type::_end()); }
    reverse_iterator rbegin(){ return reverse_iterator(end()); }
    reverse_iterator rend(){ return reverse_iterator(begin()); }
    const_reverse_iterator rbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator rend() const{ return const_reverse_iterator(begin()); }
    const_reverse_iterator crbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator crend() const{ return const_reverse_iterator(begin()); }

    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator==(const C &_other) const{
        return ::std::equal(begin(),end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator!=(const C &_other) const{ return !(*this==_other); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator< (const C &_other) const{
        return ::std::lexicographical_compare(begin(),end(),_other.begin(),_other.end());
    }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator<=(const C &_other) const{ return !(_other<*this); }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator> (const C &_other) const{ return _other<*this; }
    template<typename C,typename=::std::enable_if_t<::std::is_same_v<typename C::value_type,value_type>>>
    bool operator>=(const C &_other) const{ return !(*this<_other); }

    //emplace directly before pos
    template<typename ...Args>
    iterator emplace(const_iterator pos,Args &&...args){
        return rbtree_type::_lemplace(pos,::std::forward<Args>(args)...);
    }
    //emplace directly before end()
    template<typename ...Args>
    iterator emplace(Args &&...args){
        return rbtree_type::_lemplace(cend(),::std::forward<Args>(args)...);
    }

    iterator insert(const value_type &v){
        return rbtree_type::_lemplace(cend(),v);
    }
    iterator insert(value_type &&v){
        return rbtree_type::_lemplace(cend(),::std::move(v));
    }
    iterator insert(const_iterator pos,const value_type &v){
        return rbtree_type::_lemplace(pos,v);
    }
    iterator insert(const_iterator pos,value_type &&v){
        return rbtree_type::_lemplace(pos,::std::move(v));
    }
    iterator insert(node_type &&n){
        return rbtree_type::_linsert_node(cend(),n);
    }
    iterator insert(const_iterator pos,node_type &&n){
        return rbtree_type::_linsert_node(pos,n);
    }
    template<typename I>
    void insert(I ibegin,I iend){
        rbtree_type::insert(ibegin,iend);
    }
    void insert(::std::initializer_list<value_type> ilist){
        rbtree_type::insert(ilist.begin(),ilist.end());
    }
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    iterator insert(V2 &&v){
        return rbtree_type::_lemplace(cend(),::std::forward<V2>(v));
    }
    template<typename V2,::std::enable_if_t<::std::is_constructible_v<value_type,V2 &&>,int> E=0>
    iterator insert(const_iterator pos,V2 &&v){
        return rbtree_type::_lemplace(pos,::std::forward<V2>(v));
    }

    iterator find(const key_type &key){
        return rbtree_type::template _make<iterator>(rbtree_type::_find(key));
    }
    const_iterator find(const key_type &key) const{
        return rbtree_type::template _make<const_iterator>(rbtree_type::_find(key));
    }
};

}//namespace htl
