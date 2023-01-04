#pragma once

#include "random.cuh"
#include "polycalc.cuh"
// #include "evaluator.cuh"
__global__ void print(unsigned long long* a);
class doublePoly{
    public:
    unsigned long long* a;
    unsigned long long* b;
    int depth = 0;
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
class relienKey:public doublePoly{

};
class keyGen{
    public:
    publicKey pub;
    privateKey pri;
    relienKey relien;
    int N;
    double scale;
    unsigned long long** psiTable;
    unsigned long long** psiinvTable; 
    unsigned long long* psi;
    unsigned long long* psiinv;
    unsigned long long* q;
    unsigned long long* mu;
    // unsigned long long* ninv;
    unsigned long long* q_bit;
    unsigned long long size;
    cudaStream_t ntt;
    RNS rns;

    keyGen(int n,double scale,int size):rns(size,scale){
        N = n * 2;
        this->scale = scale;
        this->size = size;
        q = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        psi = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        psiinv = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        q_bit = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        getParams(q, psi, psiinv, q_bit, N);
        mu = (unsigned long long*)malloc(size * sizeof(unsigned long long));
        for(int i = 0; i < size; i++){
            uint128_t mu1 = uint128_t::exp2(q_bit[i] * 2);
            mu[i] = (mu1 / q[i]).low;
        }
        psiTable = (unsigned long long**)malloc(size * sizeof(unsigned long long *));
        psiinvTable = (unsigned long long**)malloc(size * sizeof(unsigned long long *));
        for(int i = 0; i < size; i++){
            cudaMalloc(&psiTable[i], N * sizeof(unsigned long long));
            cudaMalloc(&psiinvTable[i], N * sizeof(unsigned long long));
        }
        for(int i = 0; i < size; i++){
            fillTablePsi128<<<N/1024,1024>>>(psi[i], q[i], psiinv[i], psiTable[i], psiinvTable[i], log2(N));
            uint128_t mu1 = uint128_t::exp2(q_bit[i] * 2);
            mu[i] = (mu1 / q[i]).low;
        }
        ntt = 0;

        unsigned long long *pub_a,*pub_b,*pri_b,*relien_a,*relien_b,*e;
        Check(cudaMalloc((void**)&pub_a, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&pub_b, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&relien_a, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&relien_b, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&pri_b, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&e, N * size * sizeof(unsigned long long)));
        // Check(cudaMalloc((void**)&e, N * sizeof(unsigned long long)));

        genRandom<<<N/1024,1024>>>(pub_a,0);
        pub_a = rns.decompose(pub_a,N);
        for(int i = 0; i < size; i++){
            forwardNTT(pub_a + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
        }

        // genRandom<<<N/1024,1024>>>(pub_b,rns.moduleChain_h[0]);
        genRandom<<<N/1024,1024>>>(pub_b,scale);

        pub_b = rns.decompose(pub_b,N);
        print<<<1,1>>>(pub_b);
        // rns.decompose(pub_b);
        for(int  i = 0; i < size; i++){
            forwardNTT(pub_b + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
        }
        // print<<<1,1>>>(pub_b);
        genRandom_s<<<N/1024,1024>>>(pri_b,scale);
        pri_b = rns.decompose(pri_b,N);
        for(int i = 0; i < size; i++){
            forwardNTT(pri_b + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
        }
        
        for(int i = 0; i < size; i++){
            polymuladd<<<N/1024,1024>>>(pri_b + N * i,pub_b + N * i,pub_a + N * i,pub_a + N * i,q[i],mu[i],q_bit[i]);
        }
        print<<<1,1>>>(pub_b);
        pub.set(pub_a,pub_b);
        pri.set(pri_b);

        genRandom<<<N/1024,1024>>>(relien_b,0);
        relien_b = rns.decompose(relien_b,N);
        for(int i = 0; i < size; i++){
            forwardNTT(relien_b + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
        }
        unsigned long long *s_2;
        Check(cudaMalloc((void**)&s_2, N * size * sizeof(unsigned long long)));

        for(int i = 0; i < size; i++){            
            polymul<<<N/1024,1024>>>(pri_b + N * i,pri_b + N * i,s_2 + N * i,q[i],mu[i],q_bit[i]);
        }
        for(int i = 0; i < size; i++){            
            polymul<<<N/1024,1024>>>(pri_b+ N * i,relien_b + N * i,relien_a + N * i,q[i],mu[i],q_bit[i]);
        }
        for(int i = 0; i < size; i++){            
            polyminus<<<N/1024,1024>>>(s_2+ N * i,relien_a + N * i,relien_a + N * i,N,q[i]);
        } 
        
        genRandom<<<N/1024,1024>>>(e,0);
        e = rns.decompose(e,N);

        for(int i = 0; i < size; i++){            
            polyadd<<<N/1024,1024>>>(relien_a+ N * i,e + N * i,relien_a + N * i,N,q[i]);
        }        
        relien.set(relien_a,relien_b);
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
        unsigned long long** psiTable;
        unsigned long long** psiinvTable; 
        unsigned long long* psi;
        unsigned long long* psiinv;
        unsigned long long* q;
        unsigned long long* mu;
        // unsigned long long ninv;
        unsigned long long* q_bit;
        int size;
        cudaStream_t ntt;

        unsigned long long *u,*e;
        RNS rns;
        Encryptor(int n,double scale,int size):rns(size,scale){
            N = 2 * n;
            this->scale = scale;
            this->size = size;
            q = (unsigned long long *)malloc(size * sizeof(unsigned long long));
            psi = (unsigned long long *)malloc(size * sizeof(unsigned long long));
            psiinv = (unsigned long long *)malloc(size * sizeof(unsigned long long));
            q_bit = (unsigned long long *)malloc(size * sizeof(unsigned long long));
            getParams(q, psi, psiinv, q_bit, N);
            mu = (unsigned long long*)malloc(size * sizeof(unsigned long long));
            for(int i = 0; i < size; i++){
                uint128_t mu1 = uint128_t::exp2(q_bit[i] * 2);
                mu[i] = (mu1 / q[i]).low;
            }
            psiTable = (unsigned long long**)malloc(size * sizeof(unsigned long long *));
            psiinvTable = (unsigned long long**)malloc(size * sizeof(unsigned long long *));
            for(int i = 0; i < size; i++){
                cudaMalloc(&psiTable[i], N * sizeof(unsigned long long));
                cudaMalloc(&psiinvTable[i], N * sizeof(unsigned long long));
            }
            for(int i = 0; i < size; i++){
                fillTablePsi128<<<N/1024,1024>>>(psi[i], q[i], psiinv[i], psiTable[i], psiinvTable[i], log2(N));
                uint128_t mu1 = uint128_t::exp2(q_bit[i] * 2);
                mu[i] = (mu1 / q[i]).low;
            }
            ntt = 0;

            Check(cudaMalloc((void**)&u, N * size * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&e, N * size * sizeof(unsigned long long)));
        };
        cipherText encrypt(unsigned long long* plaintext,publicKey key){
            unsigned long long *a,*b;

            Check(cudaMalloc((void**)&b, N * size * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&a, N * size * sizeof(unsigned long long)));

            
            genRandom<<<N/1024,1024>>>(u,scale);
            u = rns.decompose(u,N);
            for(int i = 0; i < size; i++){
                forwardNTT(u + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
            }
            
            // print<<<1,1>>>(u);
            genRandom<<<8 * N/1024,1024>>>(e,0);
            e = rns.decompose(e,N);
            // print<<<1,1>>>(e);
            for(int i = 0; i < size; i++){
                forwardNTT(e + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
            }
            
            // print<<<1,1>>>(key.b);
            // print<<<1,1>>>(u);
            for(int i = 0; i < size; i++){            
                polymul<<<N/1024,1024>>>(key.a + N * i,u + N * i,a + N * i,q[i],mu[i],q_bit[i]);
                polymul<<<N/1024,1024>>>(key.b + N * i,u + N * i,b + N * i,q[i],mu[i],q_bit[i]);
                polyadd<<<N/1024,1024>>>(a + N * i,plaintext + N * i,a + N * i,N,q[i]);
                polyadd<<<N/1024,1024>>>(a + N * i,e + N * i,a + N * i,N,q[i]);
                // polyadd<<<N/1024,1024>>>(b + N * i,e + N * i,b + N * i,N,q[i]);
            }
            // print<<<1,1>>>(key.b);
            cipherText res;
            res.set(a,b);
            // print<<<1,1>>>(b);
            return res;
        }
        unsigned long long* decrypt(cipherText& cipher,privateKey key){
            unsigned long long *a,*b,*s;
            unsigned long long *res;
            // unsigned long long *t;
            // Check(cudaMalloc((void**)&e, N * sizeof(unsigned long long)));
            // Check(cudaMalloc((void**)&t, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&res, N * size * sizeof(unsigned long long)));
            a = cipher.a;
            b = cipher.b;
            s = key.b;
            
            for(int i = 0; i < size; i++){
                polymul<<<N/1024,1024>>>(b + N * i,s + N * i,res + N * i,q[i],mu[i],q_bit[i]);
                polyminus<<<N/1024,1024>>>(a + N * i,res + N * i,res + N * i,N,q[i]);
            }
            return res;
        }
        unsigned long long* decrypt(triplePoly& cipher,privateKey key){
            unsigned long long *a,*b,*c,*s;
            unsigned long long *res;

            a = cipher.a;
            b = cipher.b;
            c = cipher.c;
            s = key.b;
            Check(cudaMalloc((void**)&res, N * size * sizeof(unsigned long long)));
            for(int i = 0; i < size; i++){
                polymul<<<N/1024,1024>>>(c + N * i,s + N * i,res + N * i,q[i],mu[i],q_bit[i]);
                polymul<<<N/1024,1024>>>(res + N * i,s + N * i,res + N * i,q[i],mu[i],q_bit[i]);
                polymulminus<<<N/1024,1024>>>(b + N * i,s + N * i,res + N * i,res + N * i,q[i],mu[i],q_bit[i]);
                polyminus<<<N/1024,1024>>>(a + N * i,res + N * i,res + N * i,N,q[i]);   
            }
            // print<<<1,1>>>(res);
            // printf("aaa\n");
            return res;
        }
};