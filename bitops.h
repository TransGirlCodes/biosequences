//
// Created by Ben Ward (EI) on 18/01/2018.
//
#include <limits>
#include <numeric>

#ifndef BIOSEQUENCES_BITOPS_H
#define BIOSEQUENCES_BITOPS_H

namespace bitops{
    template <typename T>
    T repeatbyte(unsigned char byte) {
        static_assert(std::is_unsigned<T>::value,
                      "Type requested of repeatbyte must be unsigned.");
        return (std::numeric_limits<T>::max() / 0xff) * byte;
    }

    template <typename T>
    int count_ones_fallback(T x) {
        static_assert(std::is_unsigned<T>::value,
                      "Value passed to count_ones must be unsigned.");
        x = x - ((x >> 1) & repeatbyte<T>(0x55));
        x = (x & repeatbyte<T>(0x33)) + ((x >> 2) & repeatbyte<T>(0x33));
        return ((x + (x >> 4)) & repeatbyte<T>(0x0F)) * repeatbyte<T>(0x01);
    }

    template<typename T>
    int count_ones(T x) {
        return count_ones_fallback(x);
    }

    template<>
    int count_ones(uint16_t x) {
#if __has_builtin(__builtin_popcount)
        return __builtin_popcount(x);
#else
        return count_ones_fallback(x);
#endif
    }

    template<>
    int count_ones(uint32_t x) {
#if __has_builtin(__builtin_popcountl)
        return __builtin_popcountl(x);
#else
        return count_ones_fallback(x);
#endif
    }

    template<>
    int count_ones(uint64_t x) {
#if __has_builtin(__builtin_popcountll)
        return __builtin_popcountll(x);
#else
        return count_ones_fallback(x);
#endif
    }

    int trailing_zeros(uint8_t x) {
#if __has_builtin(__builtin_ctz)
        return __builtin_ctz(x);
#else
        int c = 8;
        x &= (int8_t)x;
        if(x) c--;
        if(x & 0x0f) c -= 4;
        if(x & 0x33) c -= 2;
        if(x & 0x55) c -= 1;
        return c;
#endif
    }

    int trailing_zeros(uint64_t x) {
#if __has_builtin(__builtin_ctzll)
        return __builtin_ctzll(x);
#else
        int c = 64;
        x &= (int64_t)x;
        if (x) c--;
        if (x & 0x00000000ffffffff) c -= 32;
        if (x & 0x0000ffff0000ffff) c -= 16;
        if (x & 0x00ff00ff00ff00ff) c -= 8;
        if (x & 0x0f0f0f0f0f0f0f0f) c -= 4;
        if (x & 0x3333333333333333) c -= 2;
        if (x & 0x5555555555555555) c -= 1;
        return c;
#endif
    }

    int trailing_zeros(uint32_t x) {
#if __has_builtin(__builtin_ctzl)
        return __builtin_ctzll(x);
#else
        int c = 32;
        x &= (int32_t)x;
        if (x) c--;
        if (x & 0x0000ffff) c -= 16;
        if (x & 0x00ff00ff) c -= 8;
        if (x & 0x0f0f0f0f) c -= 4;
        if (x & 0x33333333) c -= 2;
        if (x & 0x55555555) c -= 1;
        return c;
#endif
    }
}









#endif //BIOSEQUENCES_BITOPS_H
