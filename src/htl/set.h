#pragma once
#include"rbtree.h"

namespace htl{

template<typename K,typename ...Args>
struct _set_key_from_args{ static constexpr bool extractable=false; };
template<typename K>
struct _set_key_from_args<K,K>{
    static constexpr bool extractable=true;
    static const K &get_key(const K &k){ return k; }
};

template<typename NodeHandle,typename K>
class _set_node_handle_base{
public:
    typedef K value_type;

    value_type &value() const{ return static_cast<const NodeHandle &>(*this)._ptr->_d; }
};

template<bool M,bool L,typename K,typename C,typename A>
struct _set_traits{
    typedef K key_type;
    typedef K value_type;
    typedef C key_compare;
    typedef A allocator_type;
    static constexpr bool Multi=M,Linked=L,Set=true;
    typedef ::std::conditional_t<Linked,_linked_rbnode,_rbnode> node_base;
    typedef _rbnode_data<node_base,value_type> node_data;
    typedef _rbnode_handle<node_data,allocator_type,_set_node_handle_base,value_type> node_type;

    static const key_type &get_key(const value_type &data){ return data; }

    template<typename ...Args>
    using get_key_from_args=_set_key_from_args<K,remove_cvref_t<Args>...>;

    typedef C value_compare;
};

template<typename K,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<K>>
class linked_set;

template<typename K,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<K>>
class linked_multiset;

template<typename K,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<K>,bool Linked=false>
class set:public _rbtree<_set_traits<false,Linked,K,Comparer,Allocator>>{
protected:
    typedef _set_traits<false,Linked,K,Comparer,Allocator> TreeTraits;
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

    using rbtree_type::rbtree_type;
    using rbtree_type::operator=;

    explicit set(const set<K,Comparer,Allocator,!Linked> &_other)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other)){
    }
    set(const set<K,Comparer,Allocator,!Linked> &_other,const allocator_type &alloc)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other),alloc){
    }
    set &operator=(const set<K,Comparer,Allocator,!Linked> &_other){
        rbtree_type::template _copy_base<false>(_other);
        rbtree_type::insert(_other.begin(),_other.end());
        return *this;
    }
    explicit set(const linked_set<K,Comparer,Allocator> &)=delete;
    set(const linked_set<K,Comparer,Allocator> &,const allocator_type &)=delete;
    set &operator=(const linked_set<K,Comparer,Allocator> &)=delete;

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
};

template<typename K,typename Comparer=::std::less<K>,typename Allocator=::std::allocator<K>,bool Linked=false>
class multiset:public _rbtree<_set_traits<true,Linked,K,Comparer,Allocator>>{
protected:
    typedef _set_traits<true,Linked,K,Comparer,Allocator> TreeTraits;
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

    using rbtree_type::rbtree_type;
    using rbtree_type::operator=;

    explicit multiset(const multiset<K,Comparer,Allocator,!Linked> &_other)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other)){
    }
    multiset(const multiset<K,Comparer,Allocator,!Linked> &_other,const allocator_type &alloc)
        :rbtree_type(_other.begin(),_other.end(),rbtree_type::_get_compare(_other),alloc){
    }
    multiset &operator=(const multiset<K,Comparer,Allocator,!Linked> &_other){
        rbtree_type::template _copy_base<false>(_other);
        rbtree_type::insert(_other.begin(),_other.end());
        return *this;
    }
    explicit multiset(const linked_multiset<K,Comparer,Allocator> &)=delete;
    multiset(const linked_multiset<K,Comparer,Allocator> &,const allocator_type &)=delete;
    multiset &operator=(const linked_multiset<K,Comparer,Allocator> &)=delete;

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
};

template<typename K,typename Comparer,typename Allocator>
class linked_set:public set<K,Comparer,Allocator,true>{
public:
    typedef set<K,Comparer,Allocator,true> unlinked_type;
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

    const unlinked_type &as_unlinked() const{ return *this; }
    unlinked_type &as_unlinked(){ return *this; }

    using unlinked_type::unlinked_type;
    using rbtree_type::operator=;
    using rbtree_type::sort;

    linked_set(const linked_set &_other,const allocator_type &alloc):unlinked_type(_other.as_unlinked(),alloc){}
    linked_set(linked_set &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other.as_unlinked()),alloc){}
    explicit linked_set(const unlinked_type &_other):unlinked_type(_other){
        rbtree_type::sort();
    }
    linked_set(const unlinked_type &_other,const allocator_type &alloc):unlinked_type(_other,alloc){
        rbtree_type::sort();
    }
    linked_set &operator=(const unlinked_type &_other){
        if(this!=::std::addressof(_other)){
            rbtree_type::template _copy_base<>(_other);
            rbtree_type::template _copy_init<false,false>(_other);
        }
        rbtree_type::sort();
        return *this;
    }
    explicit linked_set(unlinked_type &&_other):unlinked_type(::std::move(_other)){
        rbtree_type::sort();
    }
    linked_set(unlinked_type &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other),alloc){
        rbtree_type::sort();
    }
    linked_set &operator=(unlinked_type &&_other){
        rbtree_type::operator=(::std::move(_other));
        rbtree_type::sort();
        return *this;
    }
    explicit linked_set(const set<K,Comparer,Allocator> &)=delete;
    explicit linked_set(const set<K,Comparer,Allocator> &,const allocator_type &)=delete;

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
    ::std::pair<iterator,bool> emplace(const_iterator pos,Args &&...args){
        return rbtree_type::_lemplace(pos,::std::forward<Args>(args)...);
    }
    //emplace directly before end()
    template<typename ...Args>
    ::std::pair<iterator,bool> emplace(Args &&...args){
        return rbtree_type::_lemplace(cend(),::std::forward<Args>(args)...);
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

    iterator find(const key_type &key){
        return rbtree_type::template _make<iterator>(rbtree_type::_find(key));
    }
    const_iterator find(const key_type &key) const{
        return rbtree_type::template _make<const_iterator>(rbtree_type::_find(key));
    }
};

template<typename K,typename Comparer,typename Allocator>
class linked_multiset:public multiset<K,Comparer,Allocator,true>{
public:
    typedef multiset<K,Comparer,Allocator,true> unlinked_type;
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

    const unlinked_type &as_unlinked() const{ return *this; }
    unlinked_type &as_unlinked(){ return *this; }

    using unlinked_type::unlinked_type;
    using rbtree_type::operator=;
    using rbtree_type::sort;

    linked_multiset(const linked_multiset &_other,const allocator_type &alloc):unlinked_type(_other.as_unlinked(),alloc){}
    linked_multiset(linked_multiset &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other.as_unlinked()),alloc){}
    explicit linked_multiset(const unlinked_type &_other):unlinked_type(_other){
        rbtree_type::sort();
    }
    linked_multiset(const unlinked_type &_other,const allocator_type &alloc):unlinked_type(_other,alloc){
        rbtree_type::sort();
    }
    linked_multiset &operator=(const unlinked_type &_other){
        if(this!=::std::addressof(_other)){
            rbtree_type::template _copy_base<>(_other);
            rbtree_type::template _copy_init<false,false>(_other);
        }
        rbtree_type::sort();
        return *this;
    }
    explicit linked_multiset(unlinked_type &&_other):unlinked_type(::std::move(_other)){
        rbtree_type::sort();
    }
    linked_multiset(unlinked_type &&_other,const allocator_type &alloc):unlinked_type(::std::move(_other),alloc){
        rbtree_type::sort();
    }
    linked_multiset &operator=(unlinked_type &&_other){
        rbtree_type::operator=(::std::move(_other));
        rbtree_type::sort();
        return *this;
    }
    explicit linked_multiset(const multiset<K,Comparer,Allocator> &)=delete;
    explicit linked_multiset(const multiset<K,Comparer,Allocator> &,const allocator_type &)=delete;

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

    iterator find(const key_type &key){
        return rbtree_type::template _make<iterator>(rbtree_type::_find(key));
    }
    const_iterator find(const key_type &key) const{
        return rbtree_type::template _make<const_iterator>(rbtree_type::_find(key));
    }
};

}//namespace htl
