#pragma once 
#include "uint128.cuh"
#include <stdlib.h>

__device__ unsigned long long modpow128(unsigned long long a, unsigned long long b, unsigned long long mod)
{
    unsigned long long res = 1;

    if (1 & b)
        res = a;

    while (b != 0)
    {
        b = b >> 1;
        uint128_t t128 = host64x2(a, a);
        a = (t128 % mod).low;
        if (b & 1)
        {
            uint128_t r128 = host64x2(res, a);
            res = (r128 % mod).low;
        }

    }
    return res;
}

__device__ unsigned modpow64(unsigned a, unsigned b, unsigned mod)
{
    unsigned res = 1;

    if (1 & b)
        res = a;

    while (b != 0)
    {
        b = b >> 1;
        unsigned long long t64 = (unsigned long long)a * a;
        a = t64 % mod;
        if (b & 1)
        {
            unsigned long long r64 = (unsigned long long)a * res;
            res = r64 % mod;
        }

    }
    return res;
}

__device__ unsigned long long modinv128(unsigned long long a, unsigned long long q)
{
    unsigned long long ainv = modpow128(a, q - 2, q);
    return ainv;
}

__device__ unsigned long long bitReverse(unsigned long long a, int bit_length)
{
    unsigned long long res = 0;

    for (int i = 0; i < bit_length; i++)
    {
        res <<= 1;
        res = (a & 1) | res;
        a >>= 1;
    }

    return res;
}

