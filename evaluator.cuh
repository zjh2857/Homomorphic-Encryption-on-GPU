#pragma once

#include "polycalc.cuh"
#include "encryptor.cuh"
#include "helper.cuh"
__global__ void print(unsigned long long* a);


class Evaluator{
    public:
    int N;
    unsigned long long** psiTable;
    unsigned long long** psiinvTable;   
    unsigned long long** psiRotTable;
    unsigned long long* psi;
    unsigned long long* psiinv;
    unsigned long long* q;
    unsigned long long* mu;
    unsigned long long* ninv;
    unsigned long long* q_bit;
    unsigned long long size;
    RNS rns;
    cudaStream_t ntt = 0;

    Evaluator(int n,int size):rns(size,1){
        N = n * 2;
        this->size = size;
        q = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        psi = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        psiinv = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        q_bit = (unsigned long long *)malloc(size * sizeof(unsigned long long));
        psiTable = (unsigned long long**)malloc(size * sizeof(unsigned long long *));
        psiinvTable = (unsigned long long**)malloc(size * sizeof(unsigned long long *));
        psiRotTable = (unsigned long long**)malloc(size * sizeof(unsigned long long *));
        getParams(q, psi, psiinv, q_bit, N);
        mu = (unsigned long long*)malloc(size * sizeof(unsigned long long));
        for(int i = 0; i < size; i++){
            cudaMalloc(&psiTable[i], N * sizeof(unsigned long long));
            cudaMalloc(&psiinvTable[i], N * sizeof(unsigned long long));
            cudaMalloc(&psiRotTable[i], N * sizeof(unsigned long long));
        }
        for(int i = 0; i < size; i++){
            fillTablePsi128<<<N/1024,1024>>>(psi[i], q[i], psiinv[i], psiTable[i], psiinvTable[i], log2(N));
            uint128_t mu1 = uint128_t::exp2(q_bit[i] * 2);
            mu[i] = (mu1 / q[i]).low;
        }
        for(int i = 0; i < size; i++){
            // unsigned long long psiRot = modpow128(psi[i],5,q[i]);
            fillTablePsi128Rot<<<N/1024,1024>>>(q[i],psiTable[i],psiRotTable[i]);
        }
    }
    // void addPlain(cipherText cipter, unsigned long long* plain){
    //     unsigned long long *a = cipter.a;
    //     polyadd<<<N/1024,1024>>>(a,plain,a,N,q);
    // }
    // void addcipter(cipherText cipter1, cipherText cipter2){

    //     polyadd<<<N/1024,1024>>>(cipter1.a,cipter2.a,cipter1.a,N,q);
        
    //     polyadd<<<N/1024,1024>>>(cipter1.b,cipter2.b,cipter1.b,N,q);
    // }
    void mulPlain(cipherText& cipher, unsigned long long* plain){
        for(int i = 0; i < size; i++){
            // print<<<1,1>>>(cipher.a + N * i);
            polymul<<<N/1024,1024>>>(cipher.a + N * i,plain + N * i,cipher.a + N * i,q[i],mu[i],q_bit[i]);
            polymul<<<N/1024,1024>>>(cipher.b + N * i,plain + N * i,cipher.b + N * i,q[i],mu[i],q_bit[i]);
            // print<<<1,1>>>(cipher.a + N * i);
        }
        rescale(cipher);
    }

    void mulPlain(unsigned long long* plain1, unsigned long long* plain2){
        // print<<<1,1>>>(plain1+2 * N + (8448 - 8190));
        for(int i = 0; i < size; i++){
            polymul<<<N/1024,1024>>>(plain1 + N * i,plain2 + N * i,plain1 + N * i,q[i],mu[i],q_bit[i]);
        }
        // for()
        // print<<<1,1>>>(plain1+2 * N + (8448 - 8190));
        
        // <<<N/1024,1024>>>(cipter.b,plain,cipter.b,q,mu,q_bit);
    }
    void rescale(cipherText cipher){
        for(int i = 0; i < size; i++){
            inverseNTT(cipher.a + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
            inverseNTT(cipher.b + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
        }
        // for(int i = 0; i < 8 * N; i += N){
        //     print_d<<<1,1>>>(cipher.a,i);
        // }
        // auto res = rns.compose(cipher.a,N);
        // print<<<1,1>>>(res);

        for(int i = cipher.depth + 1; i < size; i++){
            unsigned long long qinv = modpow128(q[cipher.depth],q[i]-2,q[i]);
            
            cudaRescale<<<N/1024,1024>>>(cipher.a + i * N,cipher.a + cipher.depth * N,q[i],mu[i],q_bit[i],qinv);

            cudaRescale<<<N/1024,1024>>>(cipher.b + i * N,cipher.b + cipher.depth * N,q[i],mu[i],q_bit[i],qinv);
        }
        // for(int i = 0; i < 8 * N; i += N){
        //     print_d<<<1,1>>>(cipher.a,i);
        // }
        // res = rns.compose(cipher.a,N);
        // print<<<1,1>>>(res);
        for(int i = 0; i < size; i++){
            forwardNTT(cipher.a+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
            forwardNTT(cipher.b+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
        }
        cipher.depth++;
    }

    void rescale(unsigned long long* cipher){
        for(int i = 0; i < size; i++){
            inverseNTT(cipher + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
            // inverseNTT(cipher.b + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
        }
        // auto res = rns.compose(cipher.a,N);
        // print<<<1,1>>>(res);
        // for(int i = 0; i < 8 * N; i += N){
        //     print_d<<<1,1>>>(cipher,i);
        // }
        for(int i = 1; i < size; i++){
            unsigned long long qinv = modpow128(q[0],q[i]-2,q[i]);
            
            cudaRescale<<<N/1024,1024>>>(cipher + i * N,cipher,q[i],mu[i],q_bit[i],qinv);

            // cudaRescale<<<N/1024,1024>>>(cipher + i * N,cipher + 1 * N,q[i],mu[i],q_bit[i],qinv);
        }
        // for(int i = 0; i < 8 * N; i += N){
        //     print_d<<<1,1>>>(cipher,i);
        // }

        for(int i = 0; i < size; i++){
            forwardNTT(cipher+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
        }
        // cipher.depth++;
    }
    triplePoly mulcipter(cipherText& cipter1, cipherText& cipter2){
        unsigned long long *a,*b,*c;

        Check(cudaMalloc((void**)&a, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&b, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&c, N * size * sizeof(unsigned long long)));
        // print<<<1,1>>>(cipter2.b);
        for(int i = 0; i < size; i++){
            polymul<<<N/1024,1024>>>(cipter1.b + N * i,cipter2.b + N * i,c + N * i,q[i],mu[i],q_bit[i]);
            polymul<<<N/1024,1024>>>(cipter1.a + N * i,cipter2.b + N * i,b + N * i,q[i],mu[i],q_bit[i]);
            polymuladd<<<N/1024,1024>>>(cipter1.b + N * i,cipter2.a + N * i,b + N * i,b + N * i,q[i],mu[i],q_bit[i]);
            polymul<<<N/1024,1024>>>(cipter1.a + N * i,cipter2.a + N * i,a + N * i,q[i],mu[i],q_bit[i]);
        }

        triplePoly res;
        res.set(a,b,c);
        return res;
    }
    // cipherText mulcipter(cipherText& cipter1, cipherText& cipter2,privateKey key){
    //     unsigned long long *a,*b,*c,*s_2,*res;

    //     Check(cudaMalloc((void**)&a, N * size * sizeof(unsigned long long)));
    //     Check(cudaMalloc((void**)&b, N * size * sizeof(unsigned long long)));
    //     Check(cudaMalloc((void**)&c, N * size * sizeof(unsigned long long)));
    //     Check(cudaMalloc((void**)&s_2, N * size * sizeof(unsigned long long)));
    //     Check(cudaMalloc((void**)&res, N * size * sizeof(unsigned long long)));
    //     // print<<<1,1>>>(cipter2.b);
    //     // for(int i = 0; i < size; i++){
    //     //     polymul<<<N/1024,1024>>>(cipter1.b,cipter2.b,c,q[i],mu[i],q_bit[i]);
    //     //     polymul<<<N/1024,1024>>>(cipter1.a,cipter2.b,b,q[i],mu[i],q_bit[i]);
    //     //     polymuladd<<<N/1024,1024>>>(cipter1.b,cipter2.a,b,b,q[i],mu[i],q_bit[i]);
    //     //     polymul<<<N/1024,1024>>>(cipter1.a,cipter2.a,a,q[i],mu[i],q_bit[i]);
    //     //     polymul<<<N/1024,1024>>>(key.b,key.b,s_2,q[i],mu[i],q_bit[i]);
    //     // }
    //     for(int i = 0; i < size; i++){
    //         polymuladd<<<N/1024,1024>>>(s_2 + N * i,c + N * i,a + N * i,res + N * i,q[i],mu[i],q_bit[i]);
    //     }        
    //     triplePoly res;
    //     res.set(a,b,c);
    //     return res;
    // }
    cipherText relien(triplePoly cipher,relienKey key){
        unsigned long long *a,*b;

        Check(cudaMalloc((void**)&a, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&b, N * size * sizeof(unsigned long long)));
        
        for(int i = 0; i < size; i++){
            polymuladd<<<N/1024,1024>>>(cipher.c + N * i,key.a + N * i,cipher.a + N * i,a + N * i,q[i],mu[i],q_bit[i]);
            polymuladd<<<N/1024,1024>>>(cipher.c + N * i,key.b + N * i,cipher.b + N * i,b + N * i,q[i],mu[i],q_bit[i]);
        }
        cipherText res;
        res.set(a,b);
        return res;
    }
    // void rotation(unsigned long long *plain){
    //     unsigned long long *temp;
    //     Check(cudaMalloc((void**)&temp, N * size * sizeof(unsigned long long)));
    //     cudaMemcpy(temp,plain,N * size * sizeof(unsigned long long),cudaMemcpyDeviceToDevice);
    //     for(int i = 0; i < size; i++){
    //         inverseNTT(plain + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
    //     }
    //     for(int i = 0; i < size; i++){
    //         cudarotation<<<N/1024,1024>>>(plain + N * i,temp + N * i,q[i],1 * 2 + 1,N);
    //     }
    //     cudaMemcpy(plain,temp,N * size * sizeof(unsigned long long),cudaMemcpyDeviceToDevice);
    //     for(int i = 0; i < size; i++){
    //         forwardNTT(plain+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
    //     }
    //     // for(int i = 0; i < size; i++){
    //     //     polyadd<<<N/1024,1024>>>(plain + N * i,temp + N * i,plain + N * i, N,q[i]);
    //     // }
    //     // print<<<1,1>>>(plain);
    // }
    void sum(unsigned long long *plain){
        unsigned long long *temp1,*temp2;
        Check(cudaMalloc((void**)&temp1, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&temp2, N * size * sizeof(unsigned long long)));
        // cudaMemcpy(temp,plain,N * size * sizeof(unsigned long long),cudaMemcpyDeviceToDevice);
        for(int i = 0; i < size; i++){
            inverseNTT(plain + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
        }
        for(int slot_idx = 0; slot_idx < N/2;slot_idx++){
            for(int i = 0; i < size; i++){
                cudarotation<<<N/1024,1024>>>(plain + N * i,plain + N * i,q[i],3,N);
            }
            for(int i = 0; i < size; i++){
                polyadd<<<N/1024,1024>>>(plain + N * i,temp2 + N * i,temp2 + N * i, N,q[i]);
            }
        }

        cudaMemcpy(plain,temp2,N * size * sizeof(unsigned long long),cudaMemcpyDeviceToDevice);
        for(int i = 0; i < size; i++){
            forwardNTT(plain+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
        }
    }
    cipherText sum(cipherText cipher, galoisKey key){
        unsigned long long *temp1,*temp2,*res_a,*res_b;
        cipherText res;
        Check(cudaMalloc((void**)&temp1, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&temp2, N * size * sizeof(unsigned long long)));
        // cudaMemcpy(temp,plain,N * size * sizeof(unsigned long long),cudaMemcpyDeviceToDevice);
        Check(cudaMalloc((void**)&res_a, N * size * sizeof(unsigned long long)));
        Check(cudaMalloc((void**)&res_b, N * size * sizeof(unsigned long long)));
        cudaMemcpy(res_a,cipher.a,N * size * sizeof(unsigned long long),cudaMemcpyDeviceToDevice);
        cudaMemcpy(res_b,cipher.b,N * size * sizeof(unsigned long long),cudaMemcpyDeviceToDevice);

        for(int slot_idx = 0; slot_idx < N/2 - 1;slot_idx++){
            for(int i = 0; i < size; i++){
                inverseNTT(cipher.a + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
                inverseNTT(cipher.b + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
            }
            for(int i = 0; i < size; i++){
                cudarotation<<<N/1024,1024>>>(cipher.a + N * i,cipher.a + N * i,q[i],3,N);
                cudarotation<<<N/1024,1024>>>(cipher.b + N * i,cipher.b + N * i,q[i],3,N);
            }

            for(int i = 0; i < size; i++){
                forwardNTT(cipher.a+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
                forwardNTT(cipher.b+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
            }

            for(int i = 0; i < size; i++){
                polymuladd<<<N/1024,1024>>>(cipher.b + N * i,key.a + N * i,cipher.a + N * i,cipher.a + N * i,q[i],mu[i],q_bit[i]);
                polymul<<<N/1024,1024>>>(cipher.b + N * i,key.b + N * i,cipher.b + N * i,q[i],mu[i],q_bit[i]);
            }
            for(int i = 0; i < size; i++){
                polyadd<<<N/1024,1024>>>(res_a + N * i,cipher.a + N * i,res_a + N * i, N,q[i]);
                polyadd<<<N/1024,1024>>>(res_b + N * i,cipher.b + N * i,res_b + N * i, N,q[i]);
            }
        }
        res.set(res_a,res_b);
        return res;
    }
    cipherText innerProduct(cipherText& cipher1,cipherText& cipher2 ,galoisKey galois,relienKey relien_key){
        auto mul_res = mulcipter(cipher1,cipher2);
        auto relien_res = relien(mul_res,relien_key);

        cipherText res = sum(relien_res,galois);
        return res;
    }

    cipherText innerProduct(cipherText& cipher1,unsigned long long *plain ,galoisKey galois,relienKey relien_key){
        mulPlain(cipher1,plain);
        // auto relien_res = relien(mul_res,relien_key);

        cipherText res = sum(cipher1,galois);
        return res;
    }
};