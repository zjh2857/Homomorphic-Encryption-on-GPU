#pragma once

#include "polycalc.cuh"
#include "encryptor.cuh"
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
        polymul<<<N/1024,1024>>>(cipter.a,plain,cipter.a,q,mu,q_bit);
        polymul<<<N/1024,1024>>>(cipter.b,plain,cipter.b,q,mu,q_bit);
    }
};