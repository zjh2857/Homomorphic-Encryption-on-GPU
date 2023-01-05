#pragma once
#include "uint128.cuh"

// __global__ void polymul(unsigned long long* a,unsigned long long* b,unsigned long long* c,int N,unsigned long long q){
//     int tid = blockDim.x * blockIdx.x +threadIdx.x;
//     a[tid] %= q;
//     b[tid] %= q;
//     c[tid] = (a[tid] * b[tid]) % q;
// }
__global__ void polyadd(unsigned long long* a,unsigned long long* b,unsigned long long* c,int N,unsigned long long q){
    int tid = blockDim.x * blockIdx.x +threadIdx.x;
    c[tid] = (a[tid] + b[tid]) % q;
}

__global__ void polyminus(unsigned long long* a,unsigned long long* b,unsigned long long* c,int N,unsigned long long q){
    int tid = blockDim.x * blockIdx.x +threadIdx.x;
    c[tid] = (a[tid] +  q - b[tid] ) % q;
}

// __global__ void polymuladd(unsigned long long* a,unsigned long long* b,unsigned long long* c,unsigned long long* d,int N,unsigned long long q){
//     int tid = blockDim.x * blockIdx.x +threadIdx.x;
//     unsigned long long t = (a[tid] * b[tid]) % q;
//     d[tid] = (t + c[tid]) % q;
// }
__global__ void polymul(unsigned long long a[], unsigned long long b[], unsigned long long c[],unsigned long long q, unsigned long long mu, int qbit)
{
    register int i = blockIdx.x * blockDim.x + threadIdx.x;

    register unsigned long long ra = a[i];
    register unsigned long long rb = b[i];

    uint128_t rc, rx;

    mul64(ra, rb, rc);

    rx = rc >> (qbit - 2);

    mul64(rx.low, mu, rx);

    uint128_t::shiftr(rx, qbit + 2);

    mul64(rx.low, q, rx);

    sub128(rc, rx);

    if (rc.low < q)
        c[i] = rc.low;
    else
        c[i] = rc.low - q;
}

__global__ void polymuladd(unsigned long long a[], unsigned long long b[], unsigned long long c[],unsigned long long d[],unsigned long long q, unsigned long long mu, int qbit)
{
    register int i = blockIdx.x * blockDim.x + threadIdx.x;

    register unsigned long long ra = a[i];
    register unsigned long long rb = b[i];

    uint128_t rc, rx;

    mul64(ra, rb, rc);

    rx = rc >> (qbit - 2);

    mul64(rx.low, mu, rx);

    uint128_t::shiftr(rx, qbit + 2);

    mul64(rx.low, q, rx);

    sub128(rc, rx);
    d[i] = (rc.low + c[i])%q;

}

__global__ void polymulminus(unsigned long long a[], unsigned long long b[], unsigned long long c[],unsigned long long d[],unsigned long long q, unsigned long long mu, int qbit)
{
    register int i = blockIdx.x * blockDim.x + threadIdx.x;

    register unsigned long long ra = a[i];
    register unsigned long long rb = b[i];

    uint128_t rc, rx;

    mul64(ra, rb, rc);

    rx = rc >> (qbit - 2);

    mul64(rx.low, mu, rx);

    uint128_t::shiftr(rx, qbit + 2);

    mul64(rx.low, q, rx);

    sub128(rc, rx);
    d[i] = (rc.low + q - c[i])%q;

}

__global__ void cudaRescale(unsigned long long *a, unsigned long long *b,unsigned long long q,unsigned long long mu,unsigned long long qbit,unsigned long long qinv){
    register int i = blockIdx.x * blockDim.x + threadIdx.x;
    register unsigned long long remainder = b[i];
    register unsigned long long ra = (a[i] + q - remainder) % q;
    register unsigned long long rb = qinv;
    // register unsigned long long rc = qinv;
    // if(i == 0){
    //     printf("%lld\n",remainder);
    // }
    uint128_t rc, rx;

    mul64(ra, rb, rc);

    rx = rc >> (qbit - 2);

    mul64(rx.low, mu, rx);

    uint128_t::shiftr(rx, qbit + 2);

    mul64(rx.low, q, rx);

    sub128(rc, rx);

    if (rc.low < q)
        a[i] = rc.low;
    else
        a[i] = rc.low - q;
}