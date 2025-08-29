#pragma once

#include<cstdint>
#include<cstddef>
#include<cstdlib>
#include<cstring>
#include<new>
#include<utility>
#include<tuple>
#include<type_traits>
#include<algorithm>
#include<limits>

// for debugging features, change this to #if 0
#if 1
#define HTL_ASSERT(...) ((void)0)
#else
#include<stdexcept>
#define HTL_ASSERT(...) do{if(!(__VA_ARGS__)){fprintf(stderr,"Assertion failed: %s, file %s, line %d\n",#__VA_ARGS__,__FILE__,__LINE__);::std::abort();}}while(0)
#endif

namespace htl{


template<typename E,typename T>
class compressed_pair final:private E{

    template<typename T1,typename T2,size_t ...I1,size_t ...I2>
    constexpr compressed_pair(T1 &t1,T2 &t2,
        ::std::index_sequence<I1...>,
        ::std::index_sequence<I2...>)
            :E(::std::get<I1>(::std::move(t1))...),
        second(::std::get<I2>(::std::move(t2))...){
    }

public:
    T second;

    const E &get_first() const{ return *this; }
    E &get_first(){ return *this; }

    constexpr compressed_pair(){}
    constexpr compressed_pair(const E &e):E(e){}
    constexpr compressed_pair(E &&e):E(::std::move(e)){}

    compressed_pair(const compressed_pair &)=default;
    compressed_pair(compressed_pair &&)=default;

    template<typename ...Args1,typename ...Args2>
    constexpr compressed_pair(::std::piecewise_construct_t,
        ::std::tuple<Args1...> args1,
        ::std::tuple<Args2...> args2)
        :compressed_pair(args1,args2,
            ::std::index_sequence_for<Args1...>{},
            ::std::index_sequence_for<Args2...>{}){
    }
};


}//namespace htl


