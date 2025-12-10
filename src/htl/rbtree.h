#pragma once
#include"base.h"

namespace htl{

struct _rbnode{
    enum _rbcolor:size_t{
        BLACK       =1,
        RED         =2,
        NODE        =3, // BLACK | RED
        SIZE_UNIT   =4
    };
    _rbnode *_p;// tree::root  | node::parent      | root_node::tree
    _rbnode *_l;// tree::front | node::left child  | leaf_node::nullptr
    _rbnode *_r;// tree::back  | node::right child | leaf_node::nullptr

    //for node:
    // BLACK: 0b000001
    //   RED: 0b000010
    //for tree:
    //  SIZE: 0bsize00
    size_t _c;

    void init(){
        _p=nullptr;
        _l=nullptr;
        _r=nullptr;
        _c=0;
    }
    //assume this!=&_other
    void move_init(_rbnode &_other){
        _p=_other._p;
        _l=_other._l;
        _r=_other._r;
        _c=_other._c;
        _other._p=nullptr;
        _other._l=nullptr;
        _other._r=nullptr;
        _other._c=0;
        HTL_ASSERT(!(_c&NODE));                 // this must be end()
        if(_c)
            _p->_p=this;
    }

    struct _insert_location{
        _rbnode *parent;
        bool on_left;
    };
    struct _bound_location:public _insert_location{
        _rbnode *pbound;
    };
    struct _hint_location:public _insert_location{
        bool duplicate;
    };
    
    //debug
    bool _in_tree() const{ return _c==BLACK||_c==RED; }
    bool _ex_tree(){
        if(!_in_tree())
            return false;
        _p=nullptr;
        _l=nullptr;
        _r=nullptr;
        _c=NODE;
        return true;
    }
private:
    static _rbnode *_min(_rbnode *_ptr){
        _rbnode *_prev;
        while(_prev=_ptr->_l)_ptr=_prev;
        return _ptr;
    }
    static _rbnode *_max(_rbnode *_ptr){
        _rbnode *_next;
        while(_next=_ptr->_r)_ptr=_next;
        return _ptr;
    }
    void _rotate_left(_rbnode *const pnode){
        _rbnode *pr=pnode->_r;
        _rbnode *prl=pr->_l;
        _rbnode *pparent=pnode->_p;
        pr->_p=pparent;
        pnode->_r=prl;
        if(prl)prl->_p=pnode;
        if(pparent==this)_p=pr;
        else if(pnode==pparent->_l)pparent->_l=pr;
        else pparent->_r=pr;
        pnode->_p=pr;
        pr->_l=pnode;
    }
    void _rotate_right(_rbnode *const pnode){
        _rbnode *pl=pnode->_l;
        _rbnode *plr=pl->_r;
        _rbnode *pparent=pnode->_p;
        pl->_p=pparent;
        pnode->_l=plr;
        if(plr)plr->_p=pnode;
        if(pparent==this)_p=pl;
        else if(pnode==pparent->_r)pparent->_r=pl;
        else pparent->_l=pl;
        pnode->_p=pl;
        pl->_r=pnode;
    }
public:
    static _rbnode *next(_rbnode *_ptr){
        HTL_ASSERT(_ptr->_in_tree());           // shall not ++end() or extracted node
        _rbnode *_next=_ptr->_r;
        if(_next)_next=_min(_next);
        else while(((_next=_ptr->_p)->_c&NODE)&&_next->_r==_ptr)
            _ptr=_next;
        return _next;
    }
    static _rbnode *prev(_rbnode *_ptr){
        _rbnode *_prev;
        HTL_ASSERT(_ptr->_c!=NODE);             //shall not --extracted node
        if(!(_ptr->_c&NODE))_prev=_ptr->_r;     // --end() = back
        else if(_prev=_ptr->_l)_prev=_max(_prev);
        else while(((_prev=_ptr->_p)->_c&NODE)&&_prev->_l==_ptr)
            _ptr=_prev;
        HTL_ASSERT(_prev&&_prev->_in_tree());   // shall not --begin()
        return _prev;
    }
    static _rbnode *begin(_rbnode *_pend){
        _rbnode *_begin=_pend->_l;
        return _begin?_begin:_pend;
    }
    void insert_node(_rbnode *const pinsert,_rbnode *pparent,bool on_left){
        HTL_ASSERT(!(_c&NODE));                 // this must be end()
        _c+=SIZE_UNIT;
        pinsert->_p=pparent;
        pinsert->_l=nullptr;
        pinsert->_r=nullptr;
        if(pparent==this){
            pinsert->_c=BLACK;
            _p=pinsert;
            _l=pinsert;
            _r=pinsert;
            return;
        }
        pinsert->_c=RED;
        if(on_left){
            pparent->_l=pinsert;
            if(pparent==_l)_l=pinsert;
        }
        else{
            pparent->_r=pinsert;
            if(pparent==_r)_r=pinsert;
        }
        _rbnode *pnode=pinsert;
        while(pparent->_c==RED){
            _rbnode *puncle=pparent->_p->_l;
            if(pparent==puncle){
                puncle=pparent->_p->_r;
                if(puncle&&puncle->_c==RED){
                    pparent->_c=BLACK;
                    puncle->_c=BLACK;
                    pnode=pparent->_p;
                    pnode->_c=RED;
                    pparent=pnode->_p;
                    continue;
                }
                if(pnode==pparent->_r){
                    _rotate_left(pparent);
                    pparent=pnode;
                }
                pnode=pparent->_p;
                _rotate_right(pnode);
                pparent->_c=BLACK;
                pnode->_c=RED;
            }
            else{
                if(puncle&&puncle->_c==RED){
                    pparent->_c=BLACK;
                    puncle->_c=BLACK;
                    pnode=pparent->_p;
                    pnode->_c=RED;
                    pparent=pnode->_p;
                    continue;
                }
                if(pnode==pparent->_l){
                    _rotate_right(pparent);
                    pparent=pnode;
                }
                pnode=pparent->_p;
                _rotate_left(pnode);
                pparent->_c=BLACK;
                pnode->_c=RED;
            }
            break;
        }
        _p->_c=BLACK;
    }
    void extract_node(_rbnode *const pextract){
        size_t excolor=pextract->_c;
        _rbnode *pnode,*pparent,*pvictim=pextract;
        HTL_ASSERT(pextract->_in_tree());       // shall not extract end() or extracted node
        HTL_ASSERT(!(_c&NODE)&&_c>=SIZE_UNIT);  // this must be end()
        if(_c-=SIZE_UNIT){
            if(!pvictim->_l)pnode=pvictim->_r;
            else if(!pvictim->_r)pnode=pvictim->_l;
            else{
                pvictim=next(pextract);
                pnode=pvictim->_r;
            }

            if(pvictim==pextract){
                pparent=pextract->_p;
                if(pnode)pnode->_p=pparent;
                if(_p==pextract)_p=pnode;
                else (pparent->_l==pextract?pparent->_l:pparent->_r)=pnode;
                if(_l==pextract)_l=pnode?_min(pnode):pparent;
                if(_r==pextract)_r=pnode?_max(pnode):pparent;
            }
            else{
                _rbnode *pexleft=pextract->_l;
                pexleft->_p=pvictim;
                pvictim->_l=pexleft;

                _rbnode *pexright=pextract->_r;
                if(pvictim==pexright)
                    pparent=pvictim;
                else{
                    pparent=pvictim->_p;
                    if(pnode)pnode->_p=pparent;

                    pparent->_l=pnode;
                    pvictim->_r=pexright;
                    pexright->_p=pvictim;
                }

                _rbnode *pexparent=pextract->_p;
                if(_p==pextract)_p=pvictim;
                else (pexparent->_l==pextract?pexparent->_l:pexparent->_r)=pvictim;
                pvictim->_p=pexparent;
                ::std::swap(pvictim->_c,excolor);
            }
        }
        else{
            _p=nullptr;
            _l=nullptr;
            _r=nullptr;
            return;
        }

        if(excolor!=BLACK)
            return;
        
        if(!pnode||pnode->_c==BLACK)do{
            _rbnode *psibling=pparent->_l;
            if(pnode==psibling){
                psibling=pparent->_r;
                if(psibling->_c==RED){
                    psibling->_c=BLACK;
                    pparent->_c=RED;
                    _rotate_left(pparent);
                    psibling=pparent->_r;
                }
                _rbnode *psl=psibling->_l;
                _rbnode *psr=psibling->_r;
                bool slblack=!psl||psl->_c==BLACK;
                bool srblack=!psr||psr->_c==BLACK;
                if(slblack&&srblack){
                    psibling->_c=RED;
                    pnode=pparent;
                    pparent=pnode->_p;
                    continue;
                }
                if(srblack){
                    psl->_c=pparent->_c;
                    _rotate_right(psibling);
                }
                else{
                    psr->_c=BLACK;
                    psibling->_c=pparent->_c;
                }
                pparent->_c=BLACK;
                _rotate_left(pparent);
                break;
            }
            else{
                if(psibling->_c==RED){
                    psibling->_c=BLACK;
                    pparent->_c=RED;
                    _rotate_right(pparent);
                    psibling=pparent->_l;
                }
                _rbnode *psl=psibling->_l;
                _rbnode *psr=psibling->_r;
                bool slblack=!psl||psl->_c==BLACK;
                bool srblack=!psr||psr->_c==BLACK;
                if(slblack&&srblack){
                    psibling->_c=RED;
                    pnode=pparent;
                    pparent=pnode->_p;
                    continue;
                }
                if(slblack){
                    psr->_c=pparent->_c;
                    _rotate_left(psibling);
                }
                else{
                    psl->_c=BLACK;
                    psibling->_c=pparent->_c;
                }
                pparent->_c=BLACK;
                _rotate_right(pparent);
                break;
            }
        } while(pnode!=_p&&pnode->_c==BLACK);

        if(pnode)pnode->_c=BLACK;
    }
    void set_head(_rbnode *const proot,size_t _otherc){
        proot->_p=this;
        _p=proot;
        _l=_min(proot);
        _r=_max(proot);
        _c=_otherc;
        HTL_ASSERT(!(_c&NODE));                 // this must be end()
    }
};

struct _linked_rbnode:public _rbnode{
    _linked_rbnode *_ll;// tree::back  | node::left (prev) | edge_node::tree
    _linked_rbnode *_lr;// tree::front | node::right(next) | edge_node::tree

    void init(){
        _rbnode::init();
        _ll=nullptr;
        _lr=nullptr;
    }
    //assume this!=&_other
    void move_init(_linked_rbnode &_other){
        _rbnode::move_init(_other);
        _ll=_other._ll;
        _lr=_other._lr;
        _other._ll=nullptr;
        _other._lr=nullptr;
        if(_c){
            _ll->_lr=this;
            _lr->_ll=this;
        }
    }
    //debug
    bool _ex_tree(){
        if(!_rbnode::_ex_tree())
            return false;
        _ll=nullptr;
        _lr=nullptr;
        return true;
    }

    static _linked_rbnode *next(_linked_rbnode *_ptr){
        HTL_ASSERT(_ptr->_in_tree());           // shall not ++end() or extracted node
        return _ptr->_lr;
    }
    static _linked_rbnode *prev(_linked_rbnode *_ptr){
        HTL_ASSERT(_ptr->_c!=NODE);             //shall not --extracted node
        _linked_rbnode *_prev=_ptr->_ll;
        HTL_ASSERT(_prev&&_prev->_in_tree());   // shall not --begin()
        return _prev;
    }
    static _linked_rbnode *begin(_linked_rbnode *_pend){
        _linked_rbnode *_begin=_pend->_lr;
        return _begin?_begin:_pend;
    }

    //insert before this==end()
    void insert_node(_linked_rbnode *const pinsert,_linked_rbnode *pparent,bool on_left){
        _rbnode::insert_node(pinsert,pparent,on_left);
        pinsert->_lr=this;
        if(pparent==this){
            pinsert->_ll=this;
            _lr=pinsert;
            _ll=pinsert;
        }
        else{
            pinsert->_ll=_ll;
            _ll->_lr=pinsert;
            _ll=pinsert;
        }
    }
    //insert before pos
    void insert_node(_linked_rbnode *pos,_linked_rbnode *const pinsert,_linked_rbnode *pparent,bool on_left){
        _rbnode::insert_node(pinsert,pparent,on_left);
        if(pparent==this){
            pinsert->_lr=this;
            pinsert->_ll=this;
            _lr=pinsert;
            _ll=pinsert;
        }
        else{
            pinsert->_lr=pos;
            pinsert->_ll=pos->_ll;
            pos->_ll->_lr=pinsert;
            pos->_ll=pinsert;
        }
    }
    void extract_node(_linked_rbnode *const pextract){
        _rbnode::extract_node(pextract);
        _linked_rbnode *pexll=pextract->_ll;
        _linked_rbnode *pexlr=pextract->_lr;
        if(pexll!=pexlr){
            pexll->_lr=pexlr;
            pexlr->_ll=pexll;
        }
        else{
            HTL_ASSERT(pexll==this&&pexlr==this);
            _ll=nullptr;
            _lr=nullptr;
        }
    }
};

template<typename B,typename V>
struct _rbnode_data:public B{
    typedef B node_base;
    typedef V value_type;
    typedef _rbnode_data *node_ptr;
    value_type _d;
    value_type *data(){
        return ::std::addressof(_d);
    }
};

template<typename TreeTraits>
class _rbtree;

template<typename NodeData,typename Allocator,template<typename ...> typename Base,typename ...BaseArgs>
class _rbnode_handle:public Base<_rbnode_handle<NodeData,Allocator,Base,BaseArgs...>,BaseArgs...>{
public:
    typedef Allocator allocator_type;
private:
    template<typename TreeTraits>
    friend class _rbtree;
    friend class Base<_rbnode_handle,BaseArgs...>;

    typedef NodeData node_data;
    typedef typename NodeData::node_base node_base;
    typedef ::std::allocator_traits<allocator_type> alloc_traits;
    typedef typename alloc_traits::template rebind_alloc <node_data> node_allocator;
    typedef typename alloc_traits::template rebind_traits<node_data> node_alloc_traits;
    typedef typename node_alloc_traits::pointer node_ptr;
    typedef _rbnode *rbnode_ptr;
    typedef node_base *base_ptr;

    node_ptr _ptr;
    ::std::optional<node_allocator> _alloc;

    template<typename ...Args>
    static node_ptr _new_node(node_allocator &_al,Args &&...args){
        node_ptr pnode=node_alloc_traits::allocate(_al,1);
        ::new(pnode) node_base;
        node_alloc_traits::construct(_al,pnode->data(),::std::forward<Args>(args)...);
        return pnode;
    }
    static void _delete_node(node_allocator &_al,node_ptr pnode){
        node_alloc_traits::destroy(_al,pnode->data());
        ((base_ptr)pnode)->~node_base();
        node_alloc_traits::deallocate(_al,pnode,1);
    }
    static void _delete_tree(node_allocator &_al,node_ptr pnode){
        while(pnode){
            _delete_tree(_al,(node_ptr)pnode->_r);
            node_ptr pnext=(node_ptr)pnode->_l;
            _delete_node(_al,pnode);
            pnode=pnext;
        }
    }
    void _clear(){
        if(_ptr){
            _delete_node(*_alloc,_ptr);
            _alloc.reset();
            _ptr=nullptr;
        }
    }
    node_ptr _release(node_allocator &_al){
        if constexpr(!node_alloc_traits::is_always_equal::value)
            // according to C++ standard, it is undefined behavior
            // inserting a node to container with unequal allocators.
            HTL_ASSERT(*_alloc==_al);
        node_ptr pnode=_ptr;
        _alloc.reset();
        _ptr=nullptr;
        return pnode;
    }
    _rbnode_handle(node_ptr ptr,const node_allocator &alloc):_ptr(ptr),_alloc(alloc){
        HTL_ASSERT(_ptr&&_ptr->_ex_tree());
    }
public:
    constexpr _rbnode_handle():_ptr(nullptr){}
    _rbnode_handle(_rbnode_handle &&_other) noexcept:_ptr(_other._ptr),_alloc(::std::move(_other._alloc)){
        _other._ptr=nullptr;
        _other._alloc.reset();
    }
    _rbnode_handle &operator=(_rbnode_handle &&_other) noexcept{
        if(!_ptr){
            if(_other._ptr){
                _ptr=_other._ptr;
                _other._ptr=nullptr;
                _alloc=::std::move(*_other._alloc);
                _other._alloc.reset();
            }
        }
        else if(!_other._ptr)
            _clear();
        else if(this!=::std::addressof(_other)){
            _delete_node(*_alloc,_ptr);
            _ptr=_other._ptr;
            _other._ptr=nullptr;
            if constexpr(node_alloc_traits::propagate_on_container_move_assignment::value)
                *_alloc=::std::move(*_other._alloc);
            else
                // according to C++ standard, it is undefined behavior
                // move-assigning between nodes with unequal allocators when propagate_on_container_move_assignment is false.
                HTL_ASSERT(*_alloc==*_other._alloc);
            _other._alloc.reset();
        }
        //else
        //    according to C++ standard, move-assigning an associative_container::node_type to itself
        // effectively clears its content. however, here defined as no effect.
        return *this;
    }
    ~_rbnode_handle(){ _clear(); }

    bool empty() const{ return !_ptr; }
    explicit operator bool() const{ return _ptr; }
    allocator_type get_allocator() const{ return static_cast<allocator_type>(*_alloc); }
    void swap(_rbnode_handle &_other){
        using ::std::swap;
        if constexpr(node_alloc_traits::propagate_on_container_swap::value)
            swap(_alloc,_other._alloc);
        else if(!_alloc||!_other._alloc)
            swap(_alloc,_other._alloc);
        else if constexpr(!node_alloc_traits::is_always_equal::value)
            HTL_ASSERT(*_alloc==*_other._alloc);
        swap(_ptr,_other._ptr);
    }

    friend void swap(_rbnode_handle &_l,_rbnode_handle &_r){ _l.swap(_r); }
};

template<typename I,typename N>
struct _rbinsert_return_type{
    I position;
    bool inserted;
    N node;
};

template<typename TreeTraits>
class _rbtree{
    template<typename TreeTraits2>
    friend class _rbtree;
public:
    typedef typename TreeTraits::key_type     key_type;
    typedef typename TreeTraits::value_type value_type;
    typedef          value_type &reference;
    typedef    const value_type &const_reference;
    typedef typename TreeTraits::allocator_type    allocator_type;
    typedef typename TreeTraits::key_compare       key_compare;
    typedef typename TreeTraits::value_compare     value_compare;
protected:
    static constexpr bool Multi=TreeTraits::Multi;
    static constexpr bool Linked=TreeTraits::Linked;
    static constexpr bool Set=TreeTraits::Set;
    typedef ::std::allocator_traits<allocator_type> alloc_traits;
public:
    typedef typename alloc_traits::pointer         pointer;
    typedef typename alloc_traits::const_pointer   const_pointer;
    typedef typename alloc_traits::size_type       size_type;
    typedef typename alloc_traits::difference_type difference_type;
    typedef typename TreeTraits::node_type node_type;
protected:
    typedef typename TreeTraits::node_base node_base;
    typedef typename TreeTraits::node_data node_data;
    typedef typename alloc_traits::template rebind_alloc <node_data> node_allocator;
    typedef typename alloc_traits::template rebind_traits<node_data> node_alloc_traits;
    typedef typename node_alloc_traits::pointer node_ptr;
    typedef _rbnode *rbnode_ptr;
    typedef node_base *base_ptr;

    template<bool Const,typename N>
    class _iterator{
        friend class _rbtree;
        node_ptr _ptr;
    public:
        typedef _rbtree::difference_type difference_type;
        typedef _rbtree::value_type      value_type;
        typedef ::std::conditional_t<Const,_rbtree::const_pointer,  _rbtree::pointer>   pointer;
        typedef ::std::conditional_t<Const,_rbtree::const_reference,_rbtree::reference> reference;
        typedef ::std::bidirectional_iterator_tag iterator_category;
    public:
        template<bool C2,typename=::std::enable_if_t<C2||!Const>>
        operator _iterator<C2,N>() const{
            _iterator<C2,N> it;
            it._ptr=_ptr;
            return it;
        }

        template<bool C2>
        bool operator==(const _iterator<C2,N> &_other) const{ return _ptr==_other._ptr; }
        template<bool C2>
        bool operator!=(const _iterator<C2,N> &_other) const{ return _ptr!=_other._ptr; }

        _iterator &operator++(){
            _ptr=(node_ptr)N::next(_ptr);
            return *this;
        }
        _iterator &operator--(){
            _ptr=(node_ptr)N::prev(_ptr);
            return *this;
        }
        _iterator operator++(int){
            _iterator _tmp=*this;
            _ptr=(node_ptr)N::next(_ptr);
            return _tmp;
        }
        _iterator operator--(int){
            _iterator _tmp=*this;
            _ptr=(node_ptr)N::prev(_ptr);
            return _tmp;
        }
        reference operator *() const{
            return *_ptr->data();
        }
        pointer operator ->() const{
            return  _ptr->data();
        }
    };
public:
    typedef _iterator<Set, _rbnode>                               iterator;
    typedef _iterator<true,_rbnode>                         const_iterator;
    typedef ::std::reverse_iterator<      iterator>       reverse_iterator;
    typedef ::std::reverse_iterator<const_iterator> const_reverse_iterator;

protected:
    static_assert(size_type(-1)>size_type(0)&&sizeof(size_type)==8,"size should be 64-bit unsigned type.");
    static_assert(::std::is_same_v<value_type,typename allocator_type::value_type>,"value_type of _rbtree and allocator mismatch.");

    compressed_pair<node_allocator,compressed_pair<key_compare,node_base>> _a;

    node_ptr _end() const{ return (node_ptr)::std::addressof(_a.second.second); }
    node_ptr _begin() const{ return (node_ptr)_rbnode::begin(_end()); }
    static const key_type &_node_key(rbnode_ptr _ptr){ return TreeTraits::get_key(((node_ptr)_ptr)->_d); }

    template<bool Lower>
    _rbnode::_bound_location _find_bound(const key_type &k) const{
        const rbnode_ptr pend=_end();
        _rbnode::_bound_location result;
        result.parent=pend;
        result.on_left=false;
        result.pbound=pend;
        rbnode_ptr ptry=pend->_p;
        const auto &_comp=_a.second.get_first();
        while(ptry){
            result.parent=ptry;
            if constexpr(Lower)
                result.on_left=!_comp(_node_key(ptry),k);
            else
                result.on_left=_comp(k,_node_key(ptry));
            if(result.on_left){
                result.pbound=ptry;
                ptry=ptry->_l;
            }
            else
                ptry=ptry->_r;
        }
        return result;
    }
    bool _lower_equal(rbnode_ptr pbound,const key_type &k) const{
        return pbound!=_end()&&!_a.second.get_first()(k,_node_key(pbound));
    }
    _rbnode::_hint_location _find_hint(rbnode_ptr phint,const key_type &k) const{
        const rbnode_ptr pend=_end();
        HTL_ASSERT(phint==pend||phint->_in_tree());
        _rbnode::_hint_location result;
        _rbnode::_insert_location &rebase=result;
        result.parent=pend;
        result.on_left=false;
        result.duplicate=false;
        const auto &_comp=_a.second.get_first();
        if constexpr(Multi){
            if(phint==pend){
                if(_rbnode *pback=pend->_r){
                    if(_comp(k,_node_key(pback)))
                        rebase=_find_bound<false>(k);
                    else
                        rebase.parent=pback;
                }
            }
            else if(_comp(_node_key(phint),k))
                rebase=_find_bound<true>(k);
            else if(phint==pend->_l){
                rebase.parent=phint;
                rebase.on_left=true;
            }
            else{
                _rbnode *const pprev=_rbnode::prev(phint);
                if(_comp(k,_node_key(pprev)))
                    rebase=_find_bound<false>(k);
                else
                    rebase.parent=(rebase.on_left=pprev->_r)?phint:pprev;
            }
        }
        else{
            do{
                if(phint==pend){
                    if(_rbnode *pback=pend->_r){
                        if(!_comp(_node_key(pback),k))break;
                        rebase.parent=pback;
                    }
                }
                else if(phint==pend->_l){
                    if(!_comp(k,_node_key(phint)))break;
                    rebase.parent=phint;
                    rebase.on_left=true;
                }
                else if(_comp(k,_node_key(phint))){
                    _rbnode *const pprev=_rbnode::prev(phint);
                    if(!_comp(_node_key(pprev),k))break;
                    rebase.parent=(rebase.on_left=pprev->_r)?phint:pprev;
                }
                else if(_comp(_node_key(phint),k)){
                    _rbnode *const pnext=_rbnode::next(phint);
                    if(pnext!=pend&&!_comp(k,_node_key(pnext)))break;
                    rebase.parent=(rebase.on_left=phint->_r)?pnext:phint;
                }
                else{
                    rebase.parent=phint;
                    result.duplicate=true;
                }
                return result;
            } while(0);
            const _rbnode::_bound_location _pos=_find_bound<true>(k);
            if(!_lower_equal(_pos.pbound,k))rebase=_pos;
            else{
                rebase.parent=_pos.pbound;
                result.duplicate=true;
            }
        }
        return result;
    }
    rbnode_ptr _find(const key_type &k) const{
        const rbnode_ptr pbound=_find_bound<true>(k).pbound;
        return _lower_equal(pbound,k)?pbound:_end();
    }
    ::std::pair<rbnode_ptr,rbnode_ptr> _equal_range(const key_type &k) const{
        const rbnode_ptr pend=_end();
        rbnode_ptr plower=pend;
        rbnode_ptr pupper=pend;
        rbnode_ptr ptry=pend->_p;
        const auto &_comp=_a.second.get_first();
        while(ptry){
            if(_comp(_node_key(ptry),k))
                ptry=ptry->_r;
            else{
                if(pupper==pend&&_comp(k,_node_key(ptry)))
                    pupper=ptry;
                plower=ptry;
                ptry=ptry->_l;
            }
        }
        ptry=pupper==pend?pend->_p:pupper->_l;
        while(ptry){
            if(_comp(k,_node_key(ptry))){
                pupper=ptry;
                ptry=ptry->_l;
            }
            else
                ptry=ptry->_r;
        }
        return {plower,pupper};
    }
    value_type *_at(const key_type &k) const{
        const rbnode_ptr pbound=_find_bound<true>(k).pbound;
        HTL_RUNTIME_ASSERT(_lower_equal(pbound,k));
        return node_ptr(pbound)->data();
    }
    template<typename N>
    node_ptr _erase(node_ptr pnode){
        node_ptr pnext=(node_ptr)N::next(pnode);
        _end()->extract_node(pnode);
        node_type::_delete_node(_a.get_first(),pnode);
        return pnext;
    }
    template<typename ...Args>
    _rbnode::_bound_location _emplace(Args &&...args){
        _rbnode::_bound_location _pos;
        typedef typename TreeTraits::template get_key_from_args<Args...> key_extractor;
        if constexpr(!Multi&&key_extractor::extractable){
            const key_type &key=key_extractor::get_key(args...);
            _pos=_find_bound<true>(key);
            if(_lower_equal(_pos.pbound,key))
                return _pos.parent=nullptr,_pos;
            _pos.pbound=node_type::_new_node(_a.get_first(),::std::forward<Args>(args)...);
        }
        else{
            node_ptr pnode=node_type::_new_node(_a.get_first(),::std::forward<Args>(args)...);
            const key_type &key=_node_key(pnode);
            _pos=_find_bound<!Multi>(key);
            if constexpr(!Multi){
                if(_lower_equal(_pos.pbound,key)){
                    node_type::_delete_node(_a.get_first(),pnode);
                    return _pos.parent=nullptr,_pos;
                }
            }
            _pos.pbound=pnode;
        }
        return _pos;
    }
    template<typename K,typename ...Args>
    _rbnode::_bound_location _try_emplace(K &&key,Args &&...args){
        static_assert(!Multi&&!Set,"try_emplace works only for map/linked_map");
        _rbnode::_bound_location _pos=_find_bound<!Multi>(key);
        if constexpr(!Multi){
            if(_lower_equal(_pos.pbound,key))
                return _pos.parent=nullptr,_pos;
        }
        _pos.pbound=node_type::_new_node(_a.get_first(),::std::piecewise_construct,
            ::std::forward_as_tuple(::std::forward<K>(key)),
            ::std::forward_as_tuple(::std::forward<Args>(args)...));
        return _pos;
    }
    template<typename K,typename ...Args>
    auto _try_emplace_hint(const_iterator hint,K &&key,Args &&...args){
        static_assert(!Multi&&!Set,"try_emplace works only for map/linked_map");
        _rbnode::_hint_location _pos=_find_hint(hint._ptr,key);
        node_ptr pnode;
        if(_pos.duplicate)pnode=(node_ptr)_pos.parent;
        else{
            pnode=node_type::_new_node(_a.get_first(),::std::piecewise_construct,
                ::std::forward_as_tuple(::std::forward<K>(key)),
                ::std::forward_as_tuple(::std::forward<Args>(args)...));
            _end()->insert_node(pnode,(base_ptr)_pos.parent,_pos.on_left);
        }
        return ::std::pair<iterator,bool>(_make<iterator>(pnode),!_pos.duplicate);
    }
    node_ptr _insert_node(node_type &n){
        if(!n._ptr)return _end();
        const key_type &key=_node_key(n._ptr);
        _rbnode::_bound_location _pos=_find_bound<!Multi>(key);
        if constexpr(!Multi){
            if(_lower_equal(_pos.pbound,key))
                return (node_ptr)_pos.pbound;
        }
        node_ptr pnode=n._release(_a.get_first());
        _end()->insert_node(pnode,(base_ptr)_pos.parent,_pos.on_left);
        return pnode;
    }
    node_ptr _insert_node(node_ptr phint,node_type &n){
        if(!n._ptr)return _end();
        _rbnode::_hint_location _pos=_find_hint(phint,_node_key(n._ptr));
        if constexpr(!Multi){
            if(_pos.duplicate)
                return (node_ptr)_pos.parent;
        }
        node_ptr pnode=n._release(_a.get_first());
        _end()->insert_node(pnode,(base_ptr)_pos.parent,_pos.on_left);
        return pnode;
    }
    node_ptr _linsert_node(node_ptr pos,node_type &n){
        if(!n._ptr)return _end();
        const key_type &key=_node_key(n._ptr);
        _rbnode::_bound_location _pos=_find_bound<!Multi>(key);
        if constexpr(!Multi){
            if(_lower_equal(_pos.pbound,key))
                return (node_ptr)_pos.pbound;
        }
        node_ptr pnode=n._release(_a.get_first());
        _end()->insert_node(pos,pnode,(base_ptr)_pos.parent,_pos.on_left);
        return pnode;
    }
    template<bool Move>
    node_ptr _copy_nodes(const node_ptr psrc){
        const node_ptr psl=(node_ptr)psrc->_l,psr=(node_ptr)psrc->_r;
        node_ptr pl=psl?_copy_nodes<Move>(psl):nullptr;
        node_ptr pnode;
        if constexpr(Move)
            pnode=node_type::_new_node(_a.get_first(),::std::move(psrc->_d));
        else
            pnode=node_type::_new_node(_a.get_first(),psrc->_d);
        pnode->_c=psrc->_c;
        node_ptr pr=psr?_copy_nodes<Move>(psr):nullptr;
        if(pnode->_l=pl)pl->_p=pnode;
        if(pnode->_r=pr)pr->_p=pnode;
        return pnode;
    }
    template<bool Move,bool L=Linked>
    void _copy_init(const _rbtree &_other){
        const node_ptr pend=_end(),potherend=_other._end();
        if constexpr(L){
            pend->init();
            for(node_ptr psrc=_other._lbegin();psrc!=potherend;psrc=(node_ptr)_linked_rbnode::next(psrc)){
                if constexpr(Move)
                    emplace(::std::move(psrc->_d));
                else
                    emplace(psrc->_d);
            }
            HTL_ASSERT(pend->_c==potherend->_c);
        }
        else{
            size_t _otherc=potherend->_c;
            if(!_otherc)
                pend->init();
            else
                pend->set_head(_copy_nodes<Move>((node_ptr)potherend->_p),_otherc);
        }
    }
    template<bool P=node_alloc_traits::propagate_on_container_copy_assignment::value,typename T>
    void _copy_base(const _rbtree<T> &_other){
        clear();
        if constexpr(P)
            _a.get_first()=_other._a.get_first();
        _a.second.get_first()=_other._a.second.get_first();
    }
    template<typename T>
    static const key_compare &_get_compare(const _rbtree<T> &_other){
        return _other._a.second.get_first();
    }
    template<typename I>
    static I _make(rbnode_ptr ptr){
        I it;
        it._ptr=(node_ptr)ptr;
        return it;
    }
protected:// for linked containers
    typedef _iterator<Set, _linked_rbnode>                               linked_iterator;
    typedef _iterator<true,_linked_rbnode>                         const_linked_iterator;
    typedef ::std::reverse_iterator<      linked_iterator>       reverse_linked_iterator;
    typedef ::std::reverse_iterator<const_linked_iterator> const_reverse_linked_iterator;
    node_ptr _lbegin() const{ return (node_ptr)_linked_rbnode::begin(_end()); }
    template<typename ...Args>
    auto _lemplace(const_linked_iterator pos,Args &&...args){
        _rbnode::_bound_location _pos=_emplace(::std::forward<Args>(args)...);
        if(_pos.parent)
            _end()->insert_node(pos._ptr,(node_ptr)_pos.pbound,(base_ptr)_pos.parent,_pos.on_left);
        if constexpr(Multi)
            return _make<linked_iterator>(_pos.pbound);
        else
            return ::std::pair<linked_iterator,bool>(_make<linked_iterator>(_pos.pbound),_pos.parent);
    }
    template<typename K,typename ...Args>
    auto _ltry_emplace(const_linked_iterator pos,K &&key,Args &&...args){
        _rbnode::_bound_location _pos=_try_emplace(::std::forward<K>(key),::std::forward<Args>(args)...);
        if(_pos.parent)
            _end()->insert_node(pos._ptr,(node_ptr)_pos.pbound,(base_ptr)_pos.parent,_pos.on_left);
        return ::std::pair<linked_iterator,bool>(_make<linked_iterator>(_pos.pbound),_pos.parent);
    }
    auto _linsert_node(const_linked_iterator pos,node_type &n){
        node_ptr psrc=n._ptr,pnode=_linsert_node(pos._ptr,n);
        if constexpr(Multi)
            return _make<linked_iterator>(pnode);
        else
            return _rbinsert_return_type<linked_iterator,node_type>{_make<linked_iterator>(pnode),psrc==pnode,::std::move(n)};
    }
public:
    _rbtree(){ _a.second.second.init(); }
    explicit _rbtree(const key_compare &comp):_a(::std::nullopt,comp){
        _a.second.second.init();
    }
    explicit _rbtree(const allocator_type &alloc):_a(alloc){
        _a.second.second.init();
    }
    _rbtree(const key_compare &comp,const allocator_type &alloc)
        :_a(alloc,comp){
        _a.second.second.init();
    }
    _rbtree(const _rbtree &_other)
        :_a(node_alloc_traits::select_on_container_copy_construction(_other._a.get_first()),_other._a.second.get_first()){
        _copy_init<false>(_other);
    }
    _rbtree(const _rbtree &_other,const allocator_type &alloc):_a(alloc,_other._a.second.get_first()){
        _copy_init<false>(_other);
    }
    _rbtree(_rbtree &&_other):_a(::std::move(_other._a.get_first()),_other._a.second.get_first()){
        _a.second.second.move_init(_other._a.second.second);
    }
    _rbtree(_rbtree &&_other,const allocator_type &alloc):_a(alloc,_other._a.second.get_first()){
        if constexpr(!node_alloc_traits::is_always_equal::value){
            if(_other._a.get_first()!=_a.get_first()){
                _copy_init<true>(_other);
                _other.clear();
                return;
            }
        }
        _a.second.second.move_init(_other._a.second.second);
    }
    template<typename I>
    _rbtree(I ibegin,I iend)
        :_rbtree(){
        insert(ibegin,iend);
    }
    template<typename I>
    _rbtree(I ibegin,I iend,const key_compare &comp)
        :_rbtree(comp){
        insert(ibegin,iend);
    }
    template<typename I>
    _rbtree(I ibegin,I iend,const allocator_type &alloc)
        :_rbtree(alloc){
        insert(ibegin,iend);
    }
    template<typename I>
    _rbtree(I ibegin,I iend,const key_compare &comp,const allocator_type &alloc)
        :_rbtree(comp,alloc){
        insert(ibegin,iend);
    }
    _rbtree(::std::initializer_list<value_type> ilist)
        :_rbtree(){
        insert(ilist);
    }
    _rbtree(::std::initializer_list<value_type> ilist,const key_compare &comp)
        :_rbtree(comp){
        insert(ilist);
    }
    _rbtree(::std::initializer_list<value_type> ilist,const allocator_type &alloc)
        :_rbtree(alloc){
        insert(ilist);
    }
    _rbtree(::std::initializer_list<value_type> ilist,const key_compare &comp,const allocator_type &alloc)
        :_rbtree(comp,alloc){
        insert(ilist);
    }
    ~_rbtree(){ clear(); }

    _rbtree &operator=(const _rbtree &_other){
        if(this!=::std::addressof(_other)){
            _copy_base(_other);
            _copy_init<false>(_other);
        }
        return *this;
    }
    _rbtree &operator=(_rbtree &&_other) noexcept(::std::is_nothrow_copy_assignable_v<key_compare>&&
        (node_alloc_traits::propagate_on_container_move_assignment::value||node_alloc_traits::is_always_equal::value)){
        if(this!=::std::addressof(_other)){
            clear();
            _a.second.get_first()=_other._a.second.get_first();
            node_allocator &this_alloc=_a.get_first();
            node_allocator &other_alloc=_other._a.get_first();
            if constexpr(node_alloc_traits::propagate_on_container_move_assignment::value)
                this_alloc=::std::move(other_alloc);
            else if constexpr(!node_alloc_traits::is_always_equal::value){
                if(this_alloc!=other_alloc){
                    _copy_init<true>(_other);
                    _other.clear();
                    return *this;
                }
            }
            _a.second.second.move_init(_other._a.second.second);
        }
        return *this;
    }
    _rbtree &operator=(::std::initializer_list<value_type> ilist){
        clear();
        insert(ilist);
        return *this;
    }
public:
    void swap(_rbtree &_other) noexcept(::std::is_nothrow_swappable_v<key_compare>&&
        (node_alloc_traits::propagate_on_container_swap::value||node_alloc_traits::is_always_equal::value)){
        if(this==::std::addressof(_other))return;

        node_allocator &this_alloc=_a.get_first();
        node_allocator &other_alloc=_other._a.get_first();
        using ::std::swap;
        if constexpr(node_alloc_traits::propagate_on_container_swap::value)
            swap(this_alloc,other_alloc);
        else if constexpr(!node_alloc_traits::is_always_equal::value)
            HTL_ASSERT(this_alloc==other_alloc);
        swap(_a.second.get_first(),_other._a.second.get_first());
        swap(_a.second.second,_other._a.second.second);
    }
    friend void swap(_rbtree &_lhs,_rbtree &_rhs) noexcept(::std::is_nothrow_swappable_v<key_compare>&&
        (node_alloc_traits::propagate_on_container_swap::value||node_alloc_traits::is_always_equal::value)){
        _lhs.swap(_rhs);
    }

    key_compare key_comp() const{ return _a.second.get_first(); }
    value_compare value_comp() const{ return _a.second.get_first(); }
    allocator_type get_allocator() const{ return static_cast<allocator_type>(_a.get_first()); }

    size_type size() const{ return _a.second.second._c>>2; }
    size_type max_size() const{
        const size_type _max=::std::min<size_type>(size_t(-1)>>2,(::std::numeric_limits<difference_type>::max)());
        return ::std::min<size_type>(_max,node_alloc_traits::max_size(_a.get_first()));
    }
    bool empty() const{ return _a.second.second._c==0; }

    void clear(){
        node_ptr pend=_end();
        if constexpr(Linked){
            base_ptr pcur=pend->_ll;
            if(!pcur)return;
            do{
                base_ptr pnext=pcur->_ll;
                node_type::_delete_node(_a.get_first(),(node_ptr)pcur);
                pcur=pnext;
            } while(pcur!=pend);
            pend->_ll=nullptr;
            pend->_lr=nullptr;
        }
        else{
            rbnode_ptr pcur=pend->_p;
            if(!pcur)return;
            node_type::_delete_tree(_a.get_first(),(node_ptr)pcur);
        }
        pend->_p=nullptr;
        pend->_l=nullptr;
        pend->_r=nullptr;
        pend->_c=0;
    }

    iterator begin(){ return _make<iterator>(_begin()); }
    iterator end(){ return _make<iterator>(_end()); }
    const_iterator begin() const{ return _make<const_iterator>(_begin()); }
    const_iterator end() const{ return _make<const_iterator>(_end()); }
    const_iterator cbegin() const{ return _make<const_iterator>(_begin()); }
    const_iterator cend() const{ return _make<const_iterator>(_end()); }
    reverse_iterator rbegin(){ return reverse_iterator(end()); }
    reverse_iterator rend(){ return reverse_iterator(begin()); }
    const_reverse_iterator rbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator rend() const{ return const_reverse_iterator(begin()); }
    const_reverse_iterator crbegin() const{ return const_reverse_iterator(end()); }
    const_reverse_iterator crend() const{ return const_reverse_iterator(begin()); }
protected:
    template<typename ...Args>
    auto try_emplace(const key_type &key,Args &&...args){
        _rbnode::_bound_location _pos=_try_emplace(key,::std::forward<Args>(args)...);
        if(_pos.parent)
            _end()->insert_node((node_ptr)_pos.pbound,(base_ptr)_pos.parent,_pos.on_left);
        return ::std::pair<iterator,bool>(_make<iterator>(_pos.pbound),_pos.parent);
    }
    template<typename ...Args>
    auto try_emplace(key_type &&key,Args &&...args){
        _rbnode::_bound_location _pos=_try_emplace(::std::move(key),::std::forward<Args>(args)...);
        if(_pos.parent)
            _end()->insert_node((node_ptr)_pos.pbound,(base_ptr)_pos.parent,_pos.on_left);
        return ::std::pair<iterator,bool>(_make<iterator>(_pos.pbound),_pos.parent);
    }
    template<typename ...Args>
    iterator try_emplace(const_iterator hint,const key_type &key,Args &&...args){
        return _try_emplace_hint(hint,key,::std::forward<Args>(args)...).first;
    }
    template<typename ...Args>
    iterator try_emplace(const_iterator hint,key_type &&key,Args &&...args){
        return _try_emplace_hint(hint,::std::move(key),::std::forward<Args>(args)...).first;
    }
    void sort(){
        static_assert(Linked,"sort works only for linked types.");
        base_ptr pend=_end(),pprev=pend,pnext=(base_ptr)pend->_l;
        if(pnext)do{
            pprev->_lr=pnext;
            pnext->_ll=pprev;
            if(pnext==pend)break;
            pprev=pnext;
            pnext=(base_ptr)_rbnode::next(pprev);
        } while(1);
    }
public:
    template<typename ...Args>
    auto emplace(Args &&...args){
        _rbnode::_bound_location _pos=_emplace(::std::forward<Args>(args)...);
        if(_pos.parent)
            _end()->insert_node((node_ptr)_pos.pbound,(base_ptr)_pos.parent,_pos.on_left);
        if constexpr(Multi)
            return _make<iterator>(_pos.pbound);
        else
            return ::std::pair<iterator,bool>(_make<iterator>(_pos.pbound),_pos.parent);
    }
    template<typename ...Args>
    iterator emplace_hint(const_iterator hint,Args &&...args){
        typedef typename TreeTraits::template get_key_from_args<Args...> key_extractor;
        node_ptr pnode;
        _rbnode::_hint_location _pos;
        if constexpr(!Multi&&key_extractor::extractable){
            _pos=_find_hint(hint._ptr,key_extractor::get_key(args...));
            if(_pos.duplicate)
                return _make<iterator>(_pos.parent);
            pnode=node_type::_new_node(_a.get_first(),::std::forward<Args>(args)...);
        }
        else{
            pnode=node_type::_new_node(_a.get_first(),::std::forward<Args>(args)...);
            _pos=_find_hint(hint._ptr,_node_key(pnode));
            if constexpr(!Multi){
                if(_pos.duplicate){
                    node_type::_delete_node(_a.get_first(),pnode);
                    return _make<iterator>(_pos.parent);
                }
            }
        }
        _end()->insert_node(pnode,(base_ptr)_pos.parent,_pos.on_left);
        return _make<iterator>(pnode);
    }

    auto insert(const value_type &v){
        return emplace(v);
    }
    auto insert(value_type &&v){
        return emplace(::std::move(v));
    }
    iterator insert(const_iterator hint,const value_type &v){
        return emplace_hint(hint,v);
    }
    iterator insert(const_iterator hint,value_type &&v){
        return emplace_hint(hint,::std::move(v));
    }
    auto insert(node_type &&n){
        node_ptr psrc=n._ptr,pnode=_insert_node(n);
        if constexpr(Multi)
            return _make<iterator>(pnode);
        else
            return _rbinsert_return_type<iterator,node_type>{_make<iterator>(pnode),psrc==pnode,::std::move(n)};
    }
    auto insert(const_iterator hint,node_type &&n){
        node_ptr psrc=n._ptr,pnode=_insert_node(hint._ptr,n);
        if constexpr(Multi)
            return _make<iterator>(pnode);
        else
            return _rbinsert_return_type<iterator,node_type>{_make<iterator>(pnode),psrc==pnode,::std::move(n)};
    }
    template<typename I>
    void insert(I ibegin,I iend){
        const_iterator this_end=cend();
        while(ibegin!=iend){
            emplace_hint(this_end,*ibegin);
            ++ibegin;
        }
    }
    void insert(::std::initializer_list<value_type> ilist){
        insert(ilist.begin(),ilist.end());
    }
    template<typename C,::std::enable_if_t<::std::is_same_v<node_type,typename C::node_type>,int> E=0>
    void merge(C &_other){
        if constexpr(!node_alloc_traits::is_always_equal::value)
            // according to C++ standard, it is undefined behavior
            // merging two container with unequal allocators.
            HTL_ASSERT(_a.get_first()== _other._a.get_first());
        node_ptr pnode,pend=_end(),potherend=_other._end();
        if(pend==potherend)
            return;
        for(typename C::iterator inext=_other.begin();(pnode=inext._ptr)!=potherend;){
            ++inext;
            const key_type &key=_node_key(pnode);
            _rbnode::_bound_location _pos=_find_bound<!Multi>(key);
            if constexpr(!Multi){
                if(_lower_equal(_pos.pbound,key))
                    continue;
            }
            potherend->extract_node(pnode);
            pend->insert_node(pnode,(base_ptr)_pos.parent,_pos.on_left);
        }
    }

    ::std::pair<iterator,iterator> equal_range(const key_type &key){
        auto eqpair=_equal_range(key);
        return ::std::pair<iterator,iterator>(_make<iterator>(eqpair.first),_make<iterator>(eqpair.second));
    }
    ::std::pair<const_iterator,const_iterator> equal_range(const key_type &key) const{
        auto eqpair=_equal_range(key);
        return ::std::pair<const_iterator,const_iterator>(_make<const_iterator>(eqpair.first),_make<const_iterator>(eqpair.second));
    }
    iterator lower_bound(const key_type &key){
        return _make<iterator>(_find_bound<true>(key).pbound);
    }
    const_iterator lower_bound(const key_type &key) const{
        return _make<const_iterator>(_find_bound<true>(key).pbound);
    }
    iterator upper_bound(const key_type &key){
        return _make<iterator>(_find_bound<false>(key).pbound);
    }
    const_iterator upper_bound(const key_type &key) const{
        return _make<const_iterator>(_find_bound<false>(key).pbound);
    }
    iterator find(const key_type &key){
        return _make<iterator>(_find(key));
    }
    const_iterator find(const key_type &key) const{
        return _make<const_iterator>(_find(key));
    }
    size_type count(const key_type &key) const{
        if constexpr(Multi){
            auto eqrange=equal_range(key);
            return ::std::distance(eqrange.first,eqrange.second);
        }
        else return _lower_equal(_find_bound<true>(key).pbound,key);
    }
    bool contains(const key_type &key) const{
        return _lower_equal(_find_bound<true>(key).pbound,key);
    }

    template<bool C,typename N>
    _iterator<false,N> erase(_iterator<C,N> it){
        return _make<_iterator<false,N>>(_erase<N>(it._ptr));
    }
    template<bool C,typename N>
    _iterator<false,N> erase(_iterator<C,N> ibegin,_iterator<C,N> iend){
        node_ptr pbegin=ibegin._ptr,pend=iend._ptr;
        while(pbegin!=pend)
            pbegin=_erase<N>(pbegin);
        return _make<_iterator<false,N>>(pbegin);
    }
    size_type erase(const key_type &key){
        size_type old_size=size();
        auto eqrange=equal_range(key);
        erase(eqrange.first,eqrange.second);
        return old_size-size();
    }
    template<bool C,typename N>
    node_type extract(_iterator<C,N> it){
        node_ptr pend=_end(),pex=it._ptr;
        pend->extract_node(pex);
        return node_type(pex,_a.get_first());
    }
    node_type extract(const key_type &key){
        node_ptr pend=_end(),pex=(node_ptr)_find(key);
        if(pend==pex)return node_type();
        pend->extract_node(pex);
        return node_type(pex,_a.get_first());
    }
};

}//namespace htl
