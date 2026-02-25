#pragma once
#include"base.h"
#include<functional>

namespace htl{

struct unexpect_t{ explicit unexpect_t()=default; };
inline constexpr unexpect_t unexpect{};

template<typename T,typename E,typename=::std::bool_constant<::std::is_void_v<T>>>
class expected;

template<typename E>
class unexpected{
    static_assert( ::std::is_object_v<E>,                   "E must be an object type.");
    static_assert(!::std::is_array_v<E>,                    "E must not be an array type.");
    static_assert(!::std::is_const_v<E>,                    "E must not be const.");
    static_assert(!::std::is_volatile_v<E>,                 "E must not be volatile.");
    static_assert(!is_specialization_v<E,htl::unexpected>,  "E must not be a specialization of unexpected.");
    static constexpr bool _check_unexpected_v=true;

    E _e;
    template<typename G>
    friend class unexpected;
    template<typename U,typename G,typename>
    friend class expected;
public:
    unexpected(const unexpected &)=default;
    unexpected(unexpected &&)=default;
    template<typename G=E,typename=::std::enable_if_t<!::std::is_same_v<unexpected,remove_cvref_t<G>>&&!::std::is_same_v<::std::in_place_t,remove_cvref_t<G>>&&::std::is_constructible_v<E,G>>>
    explicit unexpected(G &&e):_e(::std::forward<G>(e)){}
    template<typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<E,Args...>>>
    explicit unexpected(::std::in_place_t,Args &&...args):_e(::std::forward<Args>(args)...){}
    template<typename G,typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<E,::std::initializer_list<G>&,Args...>>>
    explicit unexpected(::std::in_place_t,::std::initializer_list<G> il,Args &&...args):_e(il,::std::forward<Args>(args)...){}

    const E  &error() const  &noexcept{ return _e; }
          E  &error()        &noexcept{ return _e; }
    const E &&error() const &&noexcept{ return ::std::move(_e); }
          E &&error()       &&noexcept{ return ::std::move(_e); }

    void swap(unexpected &_other) noexcept(::std::is_nothrow_swappable_v<E>){ using ::std::swap; swap(_e,_other._e); }
    friend void swap(unexpected &_l,unexpected &_r) noexcept(::std::is_nothrow_swappable_v<E>){ _l.swap(_r); }

    template<typename G>
    bool operator==(const unexpected<G> &_r) const{ return _e==_r._e; }
    template<typename G>
    bool operator!=(const unexpected<G> &_r) const{ return _e!=_r._e; }
};

template<typename E>
unexpected(E)->unexpected<E>;

template<typename T,typename E>
class expected<T,E,::std::false_type>{
    static_assert(!::std::is_reference_v<T>,                                    "T must not be a reference type.");
    static_assert(!::std::is_function_v<T>,                                     "T must not be a function type.");
    static_assert(!::std::is_array_v<T>,                                        "T must not be an array type.");
    static_assert(!::std::is_same_v<::std::remove_cv_t<T>,::std::in_place_t>,   "T must not be (possibly cv-qualified) in_place_t.");
    static_assert(!::std::is_same_v<::std::remove_cv_t<T>,unexpect_t>,          "T must not be (possibly cv-qualified) unexpect_t.");
    static_assert(!is_specialization_v<::std::remove_cv_t<T>,unexpected>,       "T must not be a (possibly cv-qualified) specialization of unexpected.");
    static_assert(unexpected<E>::_check_unexpected_v);

    union{
        T _t;
        E _e;
    };
    bool _b;

    template<typename ...Args>
    T *_construct_value(Args &&...args){ return ::new (::std::addressof(_t))T(::std::forward<Args>(args)...); }
    template<typename ...Args>
    E *_construct_error(Args &&...args){ return ::new (::std::addressof(_e))E(::std::forward<Args>(args)...); }
    void _destroy_value(){ _t.~T(); }
    void _destroy_error(){ _e.~E(); }

    template<typename N,typename O,typename ...Args>
    static void _assign(N &new_val,O &old_val,Args &&...args){
        old_val.~O();
        ::new (::std::addressof(new_val))N(::std::forward<Args>(args)...);
    }

    const T &_value() const{ HTL_ASSERT( _b); return _t; }
          T &_value()      { HTL_ASSERT( _b); return _t; }
    const E &_error() const{ HTL_ASSERT(!_b); return _e; }
          E &_error()      { HTL_ASSERT(!_b); return _e; }

    template<typename From>
    static constexpr bool may_convertible=
        ::std::disjunction_v<::std::is_same<bool,::std::remove_cv_t<T>>,::std::negation<::std::disjunction<
            ::std::is_constructible<T,From&>,
            ::std::is_constructible<T,From>,
            ::std::is_constructible<T,const From&>,
            ::std::is_constructible<T,const From>,
            ::std::is_convertible<From&,T>,
            ::std::is_convertible<From,T>,
            ::std::is_convertible<const From&,T>,
            ::std::is_convertible<const From,T>>>>
        &&::std::negation_v<::std::disjunction<
            ::std::is_constructible<unexpected<E>,From&>,
            ::std::is_constructible<unexpected<E>,From>,
            ::std::is_constructible<unexpected<E>,const From&>,
            ::std::is_constructible<unexpected<E>,const From>>>;

    template<typename U,typename G>
    struct _convert_from_value{
        static constexpr bool from_value=true;

        static constexpr bool is_copyable_v=::std::is_constructible_v<T,const U&>&&::std::is_constructible_v<E,const G&>&&may_convertible<expected<U,G>>;
        static constexpr bool is_moveable_v=::std::is_constructible_v<T,      U >&&::std::is_constructible_v<E,      G >&&may_convertible<expected<U,G>>;
        static constexpr bool is_copy_explicit_v=!(::std::is_convertible_v<const U&,T>&&::std::is_convertible_v<const G&,E>);
        static constexpr bool is_move_explicit_v=!(::std::is_convertible_v<      U ,T>&&::std::is_convertible_v<      G ,E>);
    };
    template<typename U,typename G>
    struct _convert_from_void{
        static constexpr bool from_value=false;

        static constexpr bool is_copyable_v=::std::is_default_constructible_v<T>&&::std::is_constructible_v<E,const G&>&&may_convertible<expected<U,G>>;
        static constexpr bool is_moveable_v=::std::is_default_constructible_v<T>&&::std::is_constructible_v<E,      G >&&may_convertible<expected<U,G>>;
        static constexpr bool is_copy_explicit_v=!::std::is_convertible_v<const G&,E>;
        static constexpr bool is_move_explicit_v=!::std::is_convertible_v<      G ,E>;
    };

    template<typename U,typename G>
    using convert_from=::std::conditional_t<::std::is_void_v<U>,_convert_from_void<U,G>,_convert_from_value<U,G>>;

    template<typename U>
    struct construct_from{
        static constexpr bool may_constructible_v=::std::negation_v<::std::disjunction<
            ::std::is_same<::std::in_place_t,remove_cvref_t<U>>,
            ::std::is_same<expected,remove_cvref_t<U>>,
            ::std::negation<::std::is_constructible<T,U>>,
            is_specialization<remove_cvref_t<U>,unexpected>,
            ::std::conjunction<::std::is_same<bool,::std::remove_cv_t<T>>,is_specialization<remove_cvref_t<U>,htl::expected>>>>;
        static constexpr bool is_construct_explicit_v=!::std::is_convertible_v<U,T>;
    };

    //here does not require ::std::is_nothrow_move_constructible_v<T>||::std::is_nothrow_move_constructible_v<E>
    // as c++ standard said.
    static constexpr bool is_swappable_v=::std::conjunction_v<
        ::std::is_swappable<T>,
        ::std::is_swappable<E>,
        ::std::is_move_constructible<T>,
        ::std::is_move_constructible<E>>;
    static constexpr bool is_nothrow_swappable_v=::std::conjunction_v<
        ::std::is_nothrow_swappable<T>,
        ::std::is_nothrow_swappable<E>,
        ::std::is_nothrow_move_constructible<T>,
        ::std::is_nothrow_move_constructible<E>>;

    template<typename Ex,typename F,typename R=remove_cvref_t<::std::invoke_result_t<F,decltype(::std::declval<Ex>().value())>>>
    static R _then(Ex &&_this,F &&f){
        //here does not require ::std::is_same_v<E,R::error_type>
        // as c++ standard said.
        static_assert(is_specialization_v<R,htl::expected>,"F must return expect.");
        if(!_this._b)
            return R(unexpect,::std::forward<Ex>(_this)._e);
        else
            return ::std::invoke(::std::forward<F>(f),::std::forward<Ex>(_this)._t);
    }
    template<typename Ex,typename F,typename R=remove_cvref_t<::std::invoke_result_t<F,decltype(::std::declval<Ex>().error())>>>
    static R _else(Ex &&_this,F &&f){
        //here does not require ::std::is_same_v<T,R::value_type>
        // as c++ standard said.
        static_assert(is_specialization_v<R,htl::expected>,"F must return expect.");
        if(!_this._b)
            return ::std::invoke(::std::forward<F>(f),::std::forward<Ex>(_this)._e);
        else
            return R(::std::in_place,::std::forward<Ex>(_this)._t);
    }
    template<typename Ex,typename F,typename U=::std::remove_cv_t<::std::invoke_result_t<F,decltype(::std::declval<Ex>().value())>>>
    static expected<U,E> _fvalue(Ex &&_this,F &&f){
        if(!_this._b)
            return expected<U,E>(unexpect,::std::forward<Ex>(_this)._e);
        else if constexpr(!::std::is_void_v<U>)
            return expected<U,E>(::std::in_place,::std::invoke(::std::forward<F>(f),::std::forward<Ex>(_this)._t));
        else{
            ::std::invoke(::std::forward<F>(f),::std::forward<Ex>(_this)._t);
            return expected<U,E>(::std::in_place);
        }
    }
    template<typename Ex,typename F,typename G=::std::remove_cv_t<::std::invoke_result_t<F,decltype(::std::declval<Ex>().error())>>>
    static expected<T,G> _ferror(Ex &&_this,F &&f){
        if(!_this._b)
            return expected<T,G>(unexpect,::std::invoke(::std::forward<F>(f),::std::forward<Ex>(_this)._e));
        else
            return expected<T,G>(::std::in_place,::std::forward<Ex>(_this)._t);
    }

    template<typename U,typename G,typename>
    friend class expected;
public:
    typedef T value_type;
    typedef E error_type;
    typedef unexpected<E> unexpected_type;
    template<typename U>
    using rebind=expected<U,E>;

    template<typename=::std::enable_if_t<::std::is_default_constructible_v<T>>>
    expected():_b(true){ _construct_value(); }
    expected(const expected &_other):_b(_other._b){
        if(!_b)_construct_error(_other._e);
        else   _construct_value(_other._t);
    }
    template<typename=::std::enable_if_t<::std::is_move_constructible_v<T>&&::std::is_move_constructible_v<E>>>
    expected(expected &&_other) noexcept(::std::is_nothrow_move_constructible_v<T>&&::std::is_nothrow_move_constructible_v<E>)
        :_b(_other._b){
        if(!_b)_construct_error(::std::move(_other._e));
        else   _construct_value(::std::move(_other._t));
    }

    template<typename U,typename G,typename C=convert_from<U,G>,typename=::std::enable_if_t<C::is_copyable_v&&!C::is_copy_explicit_v>>
    expected(const expected<U,G> &_other):_b(_other._b){
        if(!_b)                         _construct_error(_other._e);
        else if constexpr(C::from_value)_construct_value(_other._t);
             else                       _construct_value();
    }
    template<typename U,typename G,typename C=convert_from<U,G>,::std::enable_if_t<C::is_copyable_v&&C::is_copy_explicit_v,int> =0>
    explicit expected(const expected<U,G> &_other):_b(_other._b){
        if(!_b)                         _construct_error(_other._e);
        else if constexpr(C::from_value)_construct_value(_other._t);
             else                       _construct_value();
    }
    template<typename U,typename G,typename C=convert_from<U,G>,typename=::std::enable_if_t<C::is_moveable_v&&!C::is_move_explicit_v>>
    expected(expected<U,G> &&_other):_b(_other._b){
        if(!_b)                         _construct_error(::std::move(_other._e));
        else if constexpr(C::from_value)_construct_value(::std::move(_other._t));
             else                       _construct_value();
    }
    template<typename U,typename G,typename C=convert_from<U,G>,::std::enable_if_t<C::is_moveable_v&&C::is_move_explicit_v,int> =0>
    explicit expected(expected<U,G> &&_other):_b(_other._b){
        if(!_b)                         _construct_error(::std::move(_other._e));
        else if constexpr(C::from_value)_construct_value(::std::move(_other._t));
             else                       _construct_value();
    }

    template<typename U=::std::remove_cv_t<T>,typename C=construct_from<U>,typename=::std::enable_if_t<C::may_constructible_v&&!C::is_construct_explicit_v>>
    expected(U &&v):_b(true){ _construct_value(::std::forward<U>(v)); }
    template<typename U=::std::remove_cv_t<T>,typename C=construct_from<U>,::std::enable_if_t<C::may_constructible_v&&C::is_construct_explicit_v,int> =0>
    explicit expected(U &&v):_b(true){ _construct_value(::std::forward<U>(v)); }

    template<typename G,typename=::std::enable_if_t<::std::is_constructible_v<E,const G&>&&::std::is_convertible_v<const G&,E>>>
    expected(const unexpected<G> &e):_b(false){ _construct_error(e._e); }
    template<typename G,::std::enable_if_t<::std::is_constructible_v<E,const G&>&&!::std::is_convertible_v<const G&,E>,int> =0>
    explicit expected(const unexpected<G> &e):_b(false){ _construct_error(e._e); }
    template<typename G,typename=::std::enable_if_t<::std::is_constructible_v<E,G>&&::std::is_convertible_v<G,E>>>
    expected(unexpected<G> &&e):_b(false){ _construct_error(::std::move(e._e)); }
    template<typename G,::std::enable_if_t<::std::is_constructible_v<E,G>&&!::std::is_convertible_v<G,E>,int> =0>
    explicit expected(unexpected<G> &&e):_b(false){ _construct_error(::std::move(e._e)); }

    template<typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<T,Args...>>>
    explicit expected(::std::in_place_t,Args &&...args):_b(true){ _construct_value(::std::forward<Args>(args)...); }
    template<typename U,typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<T,::std::initializer_list<U>&,Args...>>>
    explicit expected(::std::in_place_t,::std::initializer_list<U> il,Args &&...args):_b(true){ _construct_value(il,::std::forward<Args>(args)...); }
    template<typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<E,Args...>>>
    explicit expected(unexpect_t,Args &&...args):_b(false){ _construct_error(::std::forward<Args>(args)...); }
    template<typename G,typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<E,::std::initializer_list<G>&,Args...>>>
    explicit expected(unexpect_t,::std::initializer_list<G> il,Args &&...args):_b(false){ _construct_error(il,::std::forward<Args>(args)...); }

    ~expected(){
        if(!_b)_destroy_error();
        else   _destroy_value();
    }

    expected &operator=(const expected &_other){
        //here does not require ::std::is_nothrow_move_constructible_v<T>||::std::is_nothrow_move_constructible_v<E>
        // as c++ standard said.
        static_assert(::std::conjunction_v<
            ::std::is_copy_assignable<T>,::std::is_copy_constructible<T>,
            ::std::is_copy_assignable<E>,::std::is_copy_constructible<E>>,
            "This overload is defined as deleted unless all the above conditions are true.");
        if(&_b!=&_other._b){
            if(_b>_other._b){
                _assign(_e,_t,_other._e);
                _b=false;
            }
            else if(_b<_other._b){
                _assign(_t,_e,_other._t);
                _b=true;
            }
            else if(!_b)
                _e=_other._e;
            else
                _t=_other._t;
        }
        return *this;
    }
    //here does not require ::std::is_nothrow_move_constructible_v<T>||::std::is_nothrow_move_constructible_v<E>
    // as c++ standard said.
    template<typename=::std::enable_if_t<::std::conjunction_v<
        ::std::is_move_assignable<T>,::std::is_move_constructible<T>,
        ::std::is_move_assignable<E>,::std::is_move_constructible<E>>>>
    expected &operator=(expected &&_other) noexcept(::std::conjunction_v<
        ::std::is_nothrow_move_constructible<T>,::std::is_nothrow_move_assignable<T>,
        ::std::is_nothrow_move_constructible<E>,::std::is_nothrow_move_assignable<E>>){
        if(&_b!=&_other._b){
            if(_b>_other._b){
                _assign(_e,_t,::std::move(_other._e));
                _b=false;
            }
            else if(_b<_other._b){
                _assign(_t,_e,::std::move(_other._t));
                _b=true;
            }
            else if(!_b)
                _e=::std::move(_other._e);
            else
                _t=::std::move(_other._t);
        }
        return *this;
    }
    //here does not require ::std::is_nothrow_constructible_v<T,U>||::std::is_nothrow_move_constructible_v<T>||::std::is_nothrow_move_constructible_v<E>
    // as c++ standard said.
    template<typename U=::std::remove_cv_t<T>,typename=::std::enable_if_t<::std::conjunction_v<
        ::std::negation<::std::is_same<expected,remove_cvref_t<U>>>,
        ::std::negation<is_specialization<remove_cvref_t<U>,unexpected>>,
        ::std::is_constructible<T,U>,
        ::std::is_assignable<T&,U>>>>
    expected &operator=(U &&v){
        if(!_b){
            _assign(_t,_e,::std::forward<U>(v));
            _b=true;
        }
        else
            _t=::std::forward<U>(v);
        return *this;
    }
    //here does not require ::std::is_nothrow_constructible_v<E,const G&>||::std::is_nothrow_move_constructible_v<T>||::std::is_nothrow_move_constructible_v<E>
    // as c++ standard said.
    template<typename G,typename=::std::enable_if_t<::std::conjunction_v<
        ::std::is_constructible<E,const G&>,
        ::std::is_assignable<E&,const G&>>>>
    expected &operator=(const unexpected<G> &e){
        if(!_b)
            _e=e._e;
        else{
            _assign(_e,_t,e._e);
            _b=false;
        }
        return *this;
    }
    //here does not require ::std::is_nothrow_constructible_v<E,G>||::std::is_nothrow_move_constructible_v<T>||::std::is_nothrow_move_constructible_v<E>
    // as c++ standard said.
    template<typename G,typename=::std::enable_if_t<::std::conjunction_v<
        ::std::is_constructible<E,G>,
        ::std::is_assignable<E&,G>>>>
    expected &operator=(unexpected<G> &&e){
        if(!_b)
            _e=::std::move(e._e);
        else{
            _assign(_e,_t,::std::move(e._e));
            _b=false;
        }
        return *this;
    }

    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &>>>
    auto and_then(F &&f)        &{ return _then(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &>>>
    auto and_then(F &&f) const  &{ return _then(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &&>>>
    auto and_then(F &&f)       &&{ return _then(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &&>>>
    auto and_then(F &&f) const &&{ return _then(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,T &>>>
    auto  or_else(F &&f)        &{ return _else(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,const T &>>>
    auto  or_else(F &&f) const  &{ return _else(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,T &&>>>
    auto  or_else(F &&f)       &&{ return _else(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,const T &&>>>
    auto  or_else(F &&f) const &&{ return _else(::std::move(*this),::std::forward<F>(f)); }
    
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &>>>
    auto transform       (F &&f)        &{ return _fvalue(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &>>>
    auto transform       (F &&f) const  &{ return _fvalue(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &&>>>
    auto transform       (F &&f)       &&{ return _fvalue(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &&>>>
    auto transform       (F &&f) const &&{ return _fvalue(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,T &>>>
    auto  transform_error(F &&f)        &{ return _ferror(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,const T &>>>
    auto  transform_error(F &&f) const  &{ return _ferror(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,T &&>>>
    auto  transform_error(F &&f)       &&{ return _ferror(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<T,const T &&>>>
    auto  transform_error(F &&f) const &&{ return _ferror(::std::move(*this),::std::forward<F>(f)); }

    template<typename ...Args>
    T &emplace(Args &&...args) noexcept(::std::is_nothrow_constructible_v<T,Args...>){
        if(!_b){
            _b=true;
            _destroy_error();
        }
        else    _destroy_value();
        return *_construct_value(::std::forward<Args>(args)...);
    }
    template<typename U,typename ...Args>
    T &emplace(::std::initializer_list<U> il,Args &&...args) noexcept(::std::is_nothrow_constructible_v<T,::std::initializer_list<U>&,Args...>){
        if(!_b){
            _b=true;
            _destroy_error();
        }
        else    _destroy_value();
        return *_construct_value(il,::std::forward<Args>(args)...);
    }

    template<typename=::std::enable_if_t<is_swappable_v>>
    void swap(expected &_other) noexcept(is_nothrow_swappable_v){
        using ::std::swap;
        if(_other._b>_b)_other.swap(*this);
        else if(_other._b)swap(_t,_other._t);
        else if(!_b)swap(_e,_other._e);
        else{
            if constexpr(::std::is_nothrow_move_constructible_v<E>){
                E _tmp(::std::move(_other._e));
                _other._destroy_error();
                _other._construct_value(::std::move(_t));
                _destroy_value();
                _construct_error(::std::move(_tmp));
            }
            else{
                T _tmp(::std::move(_t));
                _destroy_value();
                _construct_error(::std::move(_other._e));
                _other._destroy_error();
                _other._construct_value(::std::move(_tmp));
            }
            _b=false;
            _other._b=true;
        }
    }
    template<typename=::std::enable_if_t<is_swappable_v>>
    friend void swap(expected &_lhs,expected &_rhs) noexcept(is_nothrow_swappable_v){ _lhs.swap(_rhs); }

    const T  *operator->() const   noexcept{ return ::std::addressof(_value()); }
          T  *operator->()         noexcept{ return ::std::addressof(_value()); }
    const T  &operator *() const  &noexcept{ return _value(); }
          T  &operator *()        &noexcept{ return _value(); }
    const T &&operator *() const &&noexcept{ return ::std::move(_value()); }
          T &&operator *()       &&noexcept{ return ::std::move(_value()); }
    const T  &value() const  &{ return _value(); }
          T  &value()        &{ return _value(); }
    const T &&value() const &&{ return ::std::move(_value()); }
          T &&value()       &&{ return ::std::move(_value()); }
    const E  &error() const  &{ return _error(); }
          E  &error()        &{ return _error(); }
    const E &&error() const &&{ return ::std::move(_error()); }
          E &&error()       &&{ return ::std::move(_error()); }
    template<typename U=::std::remove_cv_t<T>>
    T value_or(U &&default_value) const &{ return  _b?            _t :static_cast<T>(::std::forward<U>(default_value)); }
    template<typename U=::std::remove_cv_t<T>>
    T value_or(U &&default_value)      &&{ return  _b?::std::move(_t):static_cast<T>(::std::forward<U>(default_value)); }
    template<typename G=E>
    E error_or(G &&defalut_value) const &{ return !_b?            _e :static_cast<E>(::std::forward<G>(defalut_value)); }
    template<typename G=E>
    E error_or(G &&defalut_value)      &&{ return !_b?::std::move(_e):static_cast<E>(::std::forward<G>(defalut_value)); }

    explicit operator bool() const noexcept{ return _b; }
    bool         has_value() const noexcept{ return _b; }

    template<typename U,typename G,typename=::std::enable_if_t<!::std::is_void_v<U>
        &&::std::is_convertible_v<decltype(::std::declval<T>()==::std::declval<U>()),bool>
        &&::std::is_convertible_v<decltype(::std::declval<E>()==::std::declval<G>()),bool>>>
    bool operator==(const expected<U,G> &_other) const{ return _b==_other._b&&(_b?bool(_t==_other._t):bool(_e==_other._e)); }
    template<typename U,typename=::std::enable_if_t<!is_specialization_v<::std::remove_cv_t<U>,htl::expected>
        &&::std::is_convertible_v<decltype(::std::declval<T>()==::std::declval<U>()),bool>>>
    bool operator==(const U &_other) const{ return _b&&bool(_t==_other); }
    template<typename G,typename=::std::enable_if_t<::std::is_convertible_v<decltype(::std::declval<E>()==::std::declval<G>()),bool>>>
    bool operator==(const unexpected<G> &_other) const{ return !_b&&bool(_e==_other._e); }
    template<typename R,typename=::std::enable_if_t<::std::is_same_v<bool,decltype(::std::declval<expected>().operator==(::std::declval<R>()))>>>
    bool operator!=(const R &_other) const{ return !this->operator==(_other); }
};

template<typename T,typename E>
class expected<T,E,::std::true_type>{
    static_assert(unexpected<E>::_check_unexpected_v);

    union{
        E _e;
    };
    bool _b;

    template<typename ...Args>
    E *_construct_error(Args &&...args){ return ::new (::std::addressof(_e))E(::std::forward<Args>(args)...); }
    void _destroy_error(){ _e.~E(); }

    T _value() const{ HTL_ASSERT(_b); }
    T _value()      { HTL_ASSERT(_b); }
    const E &_error() const{ HTL_ASSERT(!_b); return _e; }
          E &_error()      { HTL_ASSERT(!_b); return _e; }

    template<typename U,typename G>
    struct convert_from{
        static constexpr bool from_value=!::std::is_void_v<U>;

        typedef expected<U,G> From;
        static constexpr bool may_convertible=
            ::std::negation_v<::std::disjunction<
                ::std::is_constructible<unexpected<E>,From&>,
                ::std::is_constructible<unexpected<E>,From>,
                ::std::is_constructible<unexpected<E>,const From&>,
                ::std::is_constructible<unexpected<E>,const From>>>;

        static constexpr bool is_copyable_v=::std::is_constructible_v<E,const G&>&&may_convertible;
        static constexpr bool is_moveable_v=::std::is_constructible_v<E,      G >&&may_convertible;
        static constexpr bool is_copy_explicit_v=!::std::is_convertible_v<const G&,E>;
        static constexpr bool is_move_explicit_v=!::std::is_convertible_v<      G ,E>;
    };

    static constexpr bool is_swappable_v=::std::conjunction_v<
        ::std::is_swappable<E>,
        ::std::is_move_constructible<E>>;
    static constexpr bool is_nothrow_swappable_v=::std::conjunction_v<
        ::std::is_nothrow_swappable<E>,
        ::std::is_nothrow_move_constructible<E>>;

    template<typename Ex,typename F,typename R=remove_cvref_t<::std::invoke_result_t<F>>>
    static R _then(Ex &&_this,F &&f){
        //here does not require ::std::is_same_v<E,R::error_type>
        // as c++ standard said.
        static_assert(is_specialization_v<R,htl::expected>,"F must return expect.");
        if(!_this._b)
            return R(unexpect,::std::forward<Ex>(_this)._e);
        else
            return ::std::invoke(::std::forward<F>(f));
    }
    template<typename Ex,typename F,typename R=remove_cvref_t<::std::invoke_result_t<F,decltype(::std::declval<Ex>().error())>>>
    static R _else(Ex &&_this,F &&f){
        //here does not require ::std::is_same_v<T,R::value_type>
        // as c++ standard said.
        static_assert(is_specialization_v<R,htl::expected>,"F must return expect.");
        if(!_this._b)
            return ::std::invoke(::std::forward<F>(f),::std::forward<Ex>(_this)._e);
        else
            return R();
    }
    template<typename Ex,typename F,typename U=::std::remove_cv_t<::std::invoke_result_t<F>>>
    static expected<U,E> _fvalue(Ex &&_this,F &&f){
        if(!_this._b)
            return expected<U,E>(unexpect,::std::forward<Ex>(_this)._e);
        else if constexpr(!::std::is_void_v<U>)
            return expected<U,E>(::std::in_place,::std::invoke(::std::forward<F>(f)));
        else{
            ::std::invoke(::std::forward<F>(f));
            return expected<U,E>(::std::in_place);
        }
    }
    template<typename Ex,typename F,typename G=::std::remove_cv_t<::std::invoke_result_t<F,decltype(::std::declval<Ex>().error())>>>
    static expected<T,G> _ferror(Ex &&_this,F &&f){
        if(!_this._b)
            return expected<T,G>(unexpect,::std::invoke(::std::forward<F>(f),::std::forward<Ex>(_this)._e));
        else
            return expected<T,G>(::std::in_place);
    }

    template<typename U,typename G,typename>
    friend class expected;
public:
    typedef T value_type;
    typedef E error_type;
    typedef unexpected<E> unexpected_type;
    template<typename U>
    using rebind=expected<U,E>;

    expected():_b(true){}
    expected(const expected &_other):_b(_other._b){
        if(!_b)_construct_error(_other._e);
    }
    template<typename=::std::enable_if_t<::std::is_move_constructible_v<E>>>
    expected(expected &&_other) noexcept(::std::is_nothrow_move_constructible_v<E>)
        :_b(_other._b){
        if(!_b)_construct_error(::std::move(_other._e));
    }

    template<typename U,typename G,typename C=convert_from<U,G>,typename=::std::enable_if_t<C::is_copyable_v&&!C::is_copy_explicit_v>>
    expected(const expected<U,G> &_other):_b(_other._b){
        if(!_b)_construct_error(_other._e);
    }
    template<typename U,typename G,typename C=convert_from<U,G>,::std::enable_if_t<C::is_copyable_v&&C::is_copy_explicit_v,int> =0>
    explicit expected(const expected<U,G> &_other):_b(_other._b){
        if(!_b)_construct_error(_other._e);
    }
    template<typename U,typename G,typename C=convert_from<U,G>,typename=::std::enable_if_t<C::is_moveable_v&&!C::is_move_explicit_v>>
    expected(expected<U,G> &&_other):_b(_other._b){
        if(!_b)_construct_error(::std::move(_other._e));
    }
    template<typename U,typename G,typename C=convert_from<U,G>,::std::enable_if_t<C::is_moveable_v&&C::is_move_explicit_v,int> =0>
    explicit expected(expected<U,G> &&_other):_b(_other._b){
        if(!_b)_construct_error(::std::move(_other._e));
    }

    template<typename G,typename=::std::enable_if_t<::std::is_constructible_v<E,const G&>&&::std::is_convertible_v<const G&,E>>>
    expected(const unexpected<G> &e):_b(false){ _construct_error(e._e); }
    template<typename G,::std::enable_if_t<::std::is_constructible_v<E,const G&>&&!::std::is_convertible_v<const G&,E>,int> =0>
    explicit expected(const unexpected<G> &e):_b(false){ _construct_error(e._e); }
    template<typename G,typename=::std::enable_if_t<::std::is_constructible_v<E,G>&&::std::is_convertible_v<G,E>>>
    expected(unexpected<G> &&e):_b(false){ _construct_error(::std::move(e._e)); }
    template<typename G,::std::enable_if_t<::std::is_constructible_v<E,G>&&!::std::is_convertible_v<G,E>,int> =0>
    explicit expected(unexpected<G> &&e):_b(false){ _construct_error(::std::move(e._e)); }

    template<typename ...Args>
    explicit expected(::std::in_place_t,Args &&...):_b(true){}
    template<typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<E,Args...>>>
    explicit expected(unexpect_t,Args &&...args):_b(false){ _construct_error(::std::forward<Args>(args)...); }
    template<typename G,typename ...Args,typename=::std::enable_if_t<::std::is_constructible_v<E,::std::initializer_list<G>&,Args...>>>
    explicit expected(unexpect_t,::std::initializer_list<G> il,Args &&...args):_b(false){ _construct_error(il,::std::forward<Args>(args)...); }

    ~expected(){
        if(!_b)_destroy_error();
    }

    expected &operator=(const expected &_other){
        static_assert(::std::conjunction_v<
            ::std::is_copy_assignable<E>,::std::is_copy_constructible<E>>,
            "This overload is defined as deleted unless all the above conditions are true.");
        if(&_b!=&_other._b){
            if(_b>_other._b){
                _construct_error(_e,_other._e);
                _b=false;
            }
            else if(_b<_other._b){
                _destroy_error();
                _b=true;
            }
            else if(!_b)
                _e=_other._e;
        }
        return *this;
    }
    template<typename=::std::enable_if_t<::std::conjunction_v<
        ::std::is_move_assignable<E>,::std::is_move_constructible<E>>>>
    expected &operator=(expected &&_other) noexcept(::std::conjunction_v<
            ::std::is_nothrow_move_constructible<E>,::std::is_nothrow_move_assignable<E>>){
        if(&_b!=&_other._b){
            if(_b>_other._b){
                _construct_error(_e,::std::move(_other._e));
                _b=false;
            }
            else if(_b<_other._b){
                _destroy_error();
                _b=true;
            }
            else if(!_b)
                _e=::std::move(_other._e);
        }
        return *this;
    }
    template<typename G,typename=::std::enable_if_t<::std::conjunction_v<
        ::std::is_constructible<E,const G&>,
        ::std::is_assignable<E&,const G&>>>>
    expected &operator=(const unexpected<G> &e){
        if(!_b)
            _e=e._e;
        else{
            _construct_error(e._e);
            _b=false;
        }
        return *this;
    }
    template<typename G,typename=::std::enable_if_t<::std::conjunction_v<
        ::std::is_constructible<E,G>,
        ::std::is_assignable<E&,G>>>>
    expected &operator=(unexpected<G> &&e){
        if(!_b)
            _e=::std::move(e._e);
        else{
            _construct_error(::std::move(e._e));
            _b=false;
        }
        return *this;
    }

    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &>>>
    auto and_then(F &&f)        &{ return _then(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &>>>
    auto and_then(F &&f) const  &{ return _then(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &&>>>
    auto and_then(F &&f)       &&{ return _then(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &&>>>
    auto and_then(F &&f) const &&{ return _then(::std::move(*this),::std::forward<F>(f)); }
    template<typename F>
    auto  or_else(F &&f)        &{ return _else(            *this ,::std::forward<F>(f)); }
    template<typename F>
    auto  or_else(F &&f) const  &{ return _else(            *this ,::std::forward<F>(f)); }
    template<typename F>
    auto  or_else(F &&f)       &&{ return _else(::std::move(*this),::std::forward<F>(f)); }
    template<typename F>
    auto  or_else(F &&f) const &&{ return _else(::std::move(*this),::std::forward<F>(f)); }
    
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &>>>
    auto transform       (F &&f)        &{ return _fvalue(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &>>>
    auto transform       (F &&f) const  &{ return _fvalue(            *this ,::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,E &&>>>
    auto transform       (F &&f)       &&{ return _fvalue(::std::move(*this),::std::forward<F>(f)); }
    template<typename F,typename=::std::enable_if_t<::std::is_constructible_v<E,const E &&>>>
    auto transform       (F &&f) const &&{ return _fvalue(::std::move(*this),::std::forward<F>(f)); }
    template<typename F>
    auto  transform_error(F &&f)        &{ return _ferror(            *this ,::std::forward<F>(f)); }
    template<typename F>
    auto  transform_error(F &&f) const  &{ return _ferror(            *this ,::std::forward<F>(f)); }
    template<typename F>
    auto  transform_error(F &&f)       &&{ return _ferror(::std::move(*this),::std::forward<F>(f)); }
    template<typename F>
    auto  transform_error(F &&f) const &&{ return _ferror(::std::move(*this),::std::forward<F>(f)); }

    T emplace() noexcept{
        if(!_b){
            _b=true;
            _destroy_error();
        }
    }

    template<typename=::std::enable_if_t<is_swappable_v>>
    void swap(expected &_other) noexcept(is_nothrow_swappable_v){
        using ::std::swap;
        if(_other._b>_b)_other.swap(*this);
        else if(_other._b)(void)0;
        else if(!_b)swap(_e,_other._e);
        else{
            _construct_error(::std::move(_other._e));
            _other._destroy_error();
            _b=false;
            _other._b=true;
        }
    }
    template<typename=::std::enable_if_t<is_swappable_v>>
    friend void swap(expected &_lhs,expected &_rhs) noexcept(is_nothrow_swappable_v){ _lhs.swap(_rhs); }

    T operator *() const noexcept{ return _value(); }
    T value() const{ return _value(); }
    const E  &error() const  &{ return _error(); }
          E  &error()        &{ return _error(); }
    const E &&error() const &&{ return ::std::move(_error()); }
          E &&error()       &&{ return ::std::move(_error()); }
    template<typename G=E>
    E error_or(G &&defalut_value) const &{ return !_b?            _e :static_cast<E>(::std::forward<G>(defalut_value)); }
    template<typename G=E>
    E error_or(G &&defalut_value)      &&{ return !_b?::std::move(_e):static_cast<E>(::std::forward<G>(defalut_value)); }

    explicit operator bool() const noexcept{ return _b; }
    bool         has_value() const noexcept{ return _b; }

    template<typename U,typename G,typename=::std::enable_if_t<::std::is_void_v<U>
        &&::std::is_convertible_v<decltype(::std::declval<E>()==::std::declval<G>()),bool>>>
    bool operator==(const expected<U,G> &_other) const{ return _b==_other._b&&(_b||bool(_e==_other._e)); }
    template<typename G,typename=::std::enable_if_t<::std::is_convertible_v<decltype(::std::declval<E>()==::std::declval<G>()),bool>>>
    bool operator==(const unexpected<G> &_other) const{ return !_b&&bool(_e==_other._e); }
    template<typename R,typename=::std::enable_if_t<::std::is_same_v<bool,decltype(::std::declval<expected>().operator==(::std::declval<R>()))>>>
    bool operator!=(const R &_other) const{ return !this->operator==(_other); }
};


}//namespace htl
