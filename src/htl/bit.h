#include"base.h"

namespace htl{

template<typename T,::std::enable_if_t<is_unsigned_integer_v<T>,int> =0>
constexpr int popcount(T x) noexcept{
    constexpr int shr=(sizeof(unsigned long long)-sizeof(T))*CHAR_BIT;
    constexpr T m0(0x5555555555555555ull>>shr);
    constexpr T m1(0x3333333333333333ull>>shr);
    constexpr T m2(0x0f0f0f0f0f0f0f0full>>shr);
    constexpr T f (0x0101010101010101ull>>shr);
    x-=T(x>>1&m0);
    x=(x&m1)+(x>>2&m1);
    x=x+(x>>4)&m2;
    x*=f;
    x>>=(sizeof(T)-1)*CHAR_BIT;
    return int(x);
}

}
