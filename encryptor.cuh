#pragma once

#include "random.cuh"
#include "polycalc.cuh"
// #include "evaluator.cuh"
__global__ void print(unsigned long long* a);
class doublePoly{
    public:
    unsigned long long* a;
    unsigned long long* b;
    doublePoly(unsigned long long* a,unsigned long long *b){
        this->a = a;
        this->b = b;
    }
    doublePoly(){
        this->a = nullptr;
        this->b = nullptr;
    }
    void set(unsigned long long* a,unsigned long long *b){
        this->a = a;
        this->b = b;
    }
    void set(unsigned long long *b){
        this->b = b;
    }
    // ~doublePoly(){
    //     if(a)cudaFree(a);
    //     if(b)cudaFree(b);
    // }
};
class triplePoly{
    public:
    unsigned long long* a;
    unsigned long long* b;
    unsigned long long* c;
    triplePoly(unsigned long long* a,unsigned long long *b,unsigned long long *c){
        this->a = a;
        this->b = b;
        this->c = c;
    }
    triplePoly(){
        this->a = nullptr;
        this->b = nullptr;
        this->c = nullptr;
    }
    void set(unsigned long long* a,unsigned long long *b,unsigned long long *c){
        this->a = a;
        this->b = b;
        this->c = c;
    }

    // ~doublePoly(){
    //     if(a)cudaFree(a);
    //     if(b)cudaFree(b);
    // }
};

class cipherText:public doublePoly{

};

class publicKey:public doublePoly{

};

class privateKey:public doublePoly{

};
class keyGen{
    public:
    publicKey pub;
    privateKey pri;
    int N;
    double scale;
    unsigned long long* psiTable;
    unsigned long long* psiinvTable; 
    unsigned long long psi;
    unsigned long long psiinv;
    unsigned long long q;
    unsigned long long mu;
    unsigned long long ninv;
    unsigned int q_bit;
    cudaStream_t ntt;


    keyGen(int n,double scale){
        N = n * 2;
        this->scale = scale;
        getParams(q, psi, psiinv, ninv, q_bit, N);
        cudaMalloc(&psiTable, N * sizeof(unsigned long long));
        cudaMalloc(&psiinvTable, N * sizeof(unsigned long long));
        fillTablePsi128<<<N/1024,1024>>>(psi, q, psiinv, psiTable, psiinvTable, log2(N));
        uint128_t mu1 = uint128_t::exp2(q_bit * 2);
        mu = (mu1 / q).low;
        ntt = 0;

        unsigned long long *pub_a,*pub_b,*pri_b;
        Check(cudaMalloc((void**)&pub_a, N * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&pub_b, N * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&pri_b, N * sizeof(unsigned long long)));
        // Check(cudaMalloc((void**)&e, N * sizeof(unsigned long long)));

        genRandom<<<N/1024,1024>>>(pub_a,1);
        forwardNTT(pub_a,N,ntt,q,mu,q_bit,psiTable);
        
        genRandom<<<N/1024,1024>>>(pub_b,scale);
        forwardNTT(pub_b,N,ntt,q,mu,q_bit,psiTable);

        genRandom<<<N/1024,1024>>>(pri_b,scale);
        forwardNTT(pri_b,N,ntt,q,mu,q_bit,psiTable);
        
        polymuladd<<<N/1024,1024>>>(pri_b,pub_b,pub_a,pub_a,q,mu,q_bit);
        pub.set(pub_a,pub_b);
        pri.set(pri_b);

    };
    // ~keyGen(){
    //     cudaFree(psiTable);
    //     cudaFree(psiinvTable); 
    // }
};

class Encryptor{
    public:
        // Encoder encode;
        int N;
        double scale;
        unsigned long long* psiTable;
        unsigned long long* psiinvTable; 
        unsigned long long psi;
        unsigned long long psiinv;
        unsigned long long q;
        unsigned long long mu;
        unsigned long long ninv;
        unsigned int q_bit;
        cudaStream_t ntt;

        unsigned long long *u,*e;
        Encryptor(int n,double scale){
            // this->n = n;
            N = 2 * n;
            this->scale = scale;
            getParams(q, psi, psiinv, ninv, q_bit, N);
            cudaMalloc(&psiTable, N * sizeof(unsigned long long));
            cudaMalloc(&psiinvTable, N * sizeof(unsigned long long));
            fillTablePsi128<<<N/1024,1024>>>(psi, q, psiinv, psiTable, psiinvTable, log2(N));
            uint128_t mu1 = uint128_t::exp2(q_bit * 2);
            mu = (mu1 / q).low;

            ntt = 0;

            Check(cudaMalloc((void**)&u, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&e, N * sizeof(unsigned long long)));
        };
        cipherText encrypt(unsigned long long* plaintext,publicKey key){
            unsigned long long *a,*b;

            Check(cudaMalloc((void**)&b, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&a, N * sizeof(unsigned long long)));

            
            genRandom<<<N/1024,1024>>>(u,scale);
            forwardNTT(u,N,ntt,q,mu,q_bit,psiTable);
            // print<<<1,1>>>(u);
            genRandom<<<N/1024,1024>>>(e,1);
            forwardNTT(e,N,ntt,q,mu,q_bit,psiTable);
            
            
            polymul<<<N/1024,1024>>>(key.a,u,a,q,mu,q_bit);
            
            // printf("%p,%p,%p,%p\n",key.a,key.b,u,b);
            polymul<<<N/1024,1024>>>(key.b,u,b,q,mu,q_bit);
            // print<<<1,1>>>(b);
            
            polyadd<<<N/1024,1024>>>(a,plaintext,a,N,q);
            polyadd<<<N/1024,1024>>>(a,e,a,N,q);
            polyadd<<<N/1024,1024>>>(b,e,b,N,q);
            // genRandom<<<N/1024,1024>>>(a,scale);

            // forwardNTT(a,N,ntt,q,mu,q_bit,psiTable);

            // Check(cudaMalloc((void**)&s, N * sizeof(unsigned long long)));
            // genRandom<<<N/1024,1024>>>(s,scale);
            // forwardNTT(s,N,ntt,q,mu,q_bit,psiTable);            
            
            // Check(cudaMalloc((void**)&e, N * sizeof(unsigned long long)));
            // genRandom<<<N/1024,1024>>>(e,1);
            // forwardNTT(e,N,ntt,q,mu,q_bit,psiTable);

            // polymul<<<N/1024,1024>>>(s,a,b,N,q);
            // polyadd<<<N/1024,1024>>>(b,e,b,N,q);
            // polymul<<<N/1024,1024>>>(b,plaintext,b,N,q);
            
            // polymul<<<N/1024,1024>>>(a,plaintext,a,N,q);
            // print<<<1,1>>>(a);
            cipherText res;
            res.set(a,b);
            return res;
        }
        unsigned long long* decrypt(cipherText& cipher,privateKey key){
            unsigned long long *a,*b,*s;
            unsigned long long *res;
            // unsigned long long *t;
            // Check(cudaMalloc((void**)&e, N * sizeof(unsigned long long)));
            // Check(cudaMalloc((void**)&t, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&res, N * sizeof(unsigned long long)));

            a = cipher.a;
            b = cipher.b;
            s = key.b;
            
            
            // Check(cudaMalloc((void**)&s, N * sizeof(unsigned long long)));
            // genRandom<<<N/1024,1024>>>(s,scale);
            // forwardNTT(s,N,ntt,q,mu,q_bit,psiTable);

            // polymul<<<N/1024,1024>>>(s,a,t,N,q);
            // genRandom<<<N/1024,1024>>>(e,0);
            
            // forwardNTT(e,N,ntt,q,mu,q_bit,psiTable);
            // print<<<1,1>>>(b);
            polymul<<<N/1024,1024>>>(b,s,res,q,mu,q_bit);
            // print<<<1,1>>>(res);

            polyminus<<<N/1024,1024>>>(a,res,res,N,q);

            return res;
        }
        unsigned long long* decrypt(triplePoly& cipher,privateKey key){
            unsigned long long *a,*b,*c,*s;
            unsigned long long *res;

            a = cipher.a;
            b = cipher.b;
            c = cipher.c;
            s = key.b;
            Check(cudaMalloc((void**)&res, N * sizeof(unsigned long long)));
            polymul<<<N/1024,1024>>>(c,s,res,q,mu,q_bit);
            polymul<<<N/1024,1024>>>(res,s,res,q,mu,q_bit);

            polymulminus<<<N/1024,1024>>>(b,s,res,res,q,mu,q_bit);
            polyminus<<<N/1024,1024>>>(a,res,res,N,q);

            return res;
        }
};