#pragma once

#include "polycalc.cuh"
#include "encryptor.cuh"
__global__ void print(unsigned long long* a);


class Evaluator{
    public:
    int N;
    unsigned long long psi;
    unsigned long long psiinv;
    unsigned long long q;
    unsigned long long mu;
    unsigned long long ninv;
    unsigned int q_bit;
    Evaluator(int n){
        N = n * 2;
        getParams(q, psi, psiinv, ninv, q_bit, N);
        uint128_t mu1 = uint128_t::exp2(q_bit * 2);
        mu = (mu1 / q).low;

    }
    void addPlain(cipherText cipter, unsigned long long* plain){
        unsigned long long *a = cipter.a;
        polyadd<<<N/1024,1024>>>(a,plain,a,N,q);
    }
    void addcipter(cipherText cipter1, cipherText cipter2){

        polyadd<<<N/1024,1024>>>(cipter1.a,cipter2.a,cipter1.a,N,q);
        
        polyadd<<<N/1024,1024>>>(cipter1.b,cipter2.b,cipter1.b,N,q);
    }
    void mulPlain(cipherText cipter, unsigned long long* plain){
        print<<<1,1>>>(cipter.a);
        print<<<1,1>>>(plain);
        // printf("%llu,%llu,%u\n",q,mu,q_bit);
        // printf("%p,%p,%p\n",cipter.a,cipter.b,plain);
        polymul<<<N/1024,1024>>>(cipter.a,plain,cipter.a,q,mu,q_bit);
        print<<<1,1>>>(cipter.a);
        polymul<<<N/1024,1024>>>(cipter.b,plain,cipter.b,q,mu,q_bit);
        

    }
    void mulPlain(unsigned long long* plain1, unsigned long long* plain2){
        polymul<<<N/1024,1024>>>(plain1,plain2,plain1,q,mu,q_bit);
        // polymul<<<N/1024,1024>>>(cipter.b,plain,cipter.b,q,mu,q_bit);
    }
    void mulcipter(cipherText cipter1, cipherText cipter2,triplePoly& res){
        unsigned long long *a,*b,*c;

        Check(cudaMalloc((void**)&a, N * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&b, N * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&c, N * sizeof(unsigned long long)));

        polymul<<<N/1024,1024>>>(cipter1.b,cipter2.b,c,q,mu,q_bit);

        polymul<<<N/1024,1024>>>(cipter1.a,cipter2.b,b,q,mu,q_bit);
        polymuladd<<<N/1024,1024>>>(cipter1.b,cipter2.a,b,b,q,mu,q_bit);
        polymul<<<N/1024,1024>>>(cipter1.a,cipter2.a,a,q,mu,q_bit);
        res.set(a,b,c);
    }
};