#pragma once
#include "parameter.cuh"
#include "uint128.cuh"
#include "ntt_60bit.cuh"
#include "cuda_runtime_api.h"
#include "device_launch_parameters.h"
#include "cufft.h"
#include <iostream>
#include <cmath>
#include "cublas_v2.h"
#include "rns.cuh"

// #include "random.cuh"
#include "polycalc.cuh"
__global__ void genRandom(unsigned long long *randomVec,unsigned long long scale);
// typedef unsigned long long long long;
#define Check(call)														\
{																		\
	cudaError_t status = call;											\
	if (status != cudaSuccess)											\
	{																	\
		cout << "行號:" << __LINE__ << endl;							\
		cout << "錯誤:" << cudaGetErrorString(status) << endl;			\
	}																	\
}
using namespace std;
const double pi = 3.14159265359;;
// __global__ void delta(cuDoubleComplex* in,unsigned long long* out,double scale,unsigned long long q,int N){

//     int tid = blockDim.x * blockIdx.x + threadIdx.x;
//     long long t = in[tid].x / N  * scale;
//     t = (t+q)%q;
//     out[tid] = t; 
// }
// __global__ void zerocnt(unsigned long long* a,int N){
//     int cnt = 0;
//     for(int i = 0; i < N; i++){
//         if(a[i] != 0){
//             printf("%d,%llu\t",i,a[i]);
//             cnt++;
//         }
//     }printf("\n");
//     printf("%d\n",cnt);
// }
// __global__ void eqcnt(unsigned long long* a,int N){
//     int cnt = 0;
//     for(int i = 1; i < N; i++){
//         if(a[i] != a[i-1]){
//             cnt++;
//             printf("%d,%llu\n",i,a[i]);
//         }
//     }printf("\n");
//     printf("%d\n",cnt);
// }
// __global__ void deltainv(unsigned long long* in,cuDoubleComplex* out,double scale,unsigned long long q){
//     int tid = blockDim.x * blockIdx.x + threadIdx.x;
//     long long t = in[tid];
//     if(t > q/2){
//         t -= q;
//     }
//     out[tid].x = (double)(t * 1.0 / scale);
//     out[tid].y = 0;
// }
__global__ void print(unsigned long long* a);
__global__ void print_d(unsigned long long* a,int d);

//     printf("%llu\n",a[0]);
// }
__global__ void print(cuDoubleComplex* h_A){
    for(int i = 0;i < 8; i++){
        printf("%lf+%lfi\n",h_A[i].x/(1llu<<40),h_A[i].y);
    }
    printf("\n");
}
__global__ void initfft(cuDoubleComplex* h_A,int N){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int i = tid/N;
    int j = tid%N;
    double real = cos(-pi/N*i*(2*j+1));
    double imag = sin(-pi/N*i*(2*j+1));
    h_A[i*N+j].x = real;
    h_A[i*N+j].y = imag;
}
__global__ void initfftinv(cuDoubleComplex* h_A,int N){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int i = tid/N;
    int j = tid%N;
    double real = cos(pi/N*i*(2*j+1));
    double imag = sin(pi/N*i*(2*j+1));
    h_A[i*N+j].x = real;
    h_A[i*N+j].y = imag;
}
__global__ void product(cuDoubleComplex* h_A){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    h_A[tid].x = h_A[tid].x * h_A[tid].x;
}
class Encoder{
    public:
        int n;
        int N;
        double scale;
        unsigned long long** psiTable;
        unsigned long long** psiinvTable; 
        unsigned long long* psi;
        unsigned long long* psiinv;
        unsigned long long* q;
        unsigned long long* mu;
        unsigned long long* q_bit;
        int BATCH = 1;
        int decomposeSize;
        RNS rns;
        cuDoubleComplex *d_A, *inv_A;

        cublasHandle_t handle;
        cublasHandle_t handleinv;
        // cufftHandle cufftForwrdHandle, cufftInverseHandle;
        Encoder(int n,double scale,int decomposeSize):rns(decomposeSize,scale){
            this->n = n;
            this->decomposeSize = decomposeSize;
            N = n * 2;
            // this->N = N;
            this->scale = scale;
            q = (unsigned long long *)malloc(decomposeSize * sizeof(unsigned long long));
            psi = (unsigned long long *)malloc(decomposeSize * sizeof(unsigned long long));
            psiinv = (unsigned long long *)malloc(decomposeSize * sizeof(unsigned long long));
            q_bit = (unsigned long long *)malloc(decomposeSize * sizeof(unsigned long long));
            
            getParams(q, psi, psiinv, q_bit, N);
            
            psiTable = (unsigned long long**)malloc(decomposeSize * sizeof(unsigned long long *));
            psiinvTable = (unsigned long long**)malloc(decomposeSize * sizeof(unsigned long long *));
            mu = (unsigned long long*)malloc(decomposeSize * sizeof(unsigned long long *));
            for(int i = 0; i < decomposeSize; i++){
                cudaMalloc(&psiTable[i], N * sizeof(unsigned long long));
                cudaMalloc(&psiinvTable[i], N * sizeof(unsigned long long));
                // printf("%p,%p\n",psiTable[i],psiinvTable[i]);
            }
            
            for(int i = 0; i < decomposeSize; i++){
                fillTablePsi128<<<N/1024,1024>>>(psi[i], q[i], psiinv[i], psiTable[i], psiinvTable[i], log2(N));
                uint128_t mu1 = uint128_t::exp2(q_bit[i] * 2);
                // printf("%llu\n",(mu1/q[i]).low);
                mu[i] = (mu1 / q[i]).low;
            }

            cublasCreate(&handle);
            cublasCreate(&handleinv);
            cudaMalloc (
                (void**)&d_A,   
                N*N * sizeof(cuDoubleComplex)    
            );
            cudaMalloc (
                (void**)&inv_A,    
                N*N * sizeof(cuDoubleComplex)    
            );
            initfft<<<N*N/1024,1024>>>(d_A,N);
            initfftinv<<<N*N/1024,1024>>>(inv_A,N);

        };

        unsigned long long* encode(double* plainVec){
            cuDoubleComplex *fft_in,*host_in;
            cuDoubleComplex *fft_out;
            unsigned long long *ntt_in;
            cudaStream_t ntt = 0;
            Check(cudaMallocHost((void**)&host_in, N * sizeof(cuDoubleComplex)));
            Check(cudaMalloc((void**)&ntt_in, N * decomposeSize * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&fft_in, N * sizeof(cuDoubleComplex)));
            Check(cudaMalloc((void**)&fft_out, N * sizeof(cuDoubleComplex)));
            for(int i = 0; i < n; i++){
                host_in[i].x = plainVec[i];
                host_in[i].y = 0;
            }
            for(int i = 0; i < n; i++){
                host_in[N-i-1].x = plainVec[i];
                host_in[N-i-1].y = 0;
            }
            
            Check(cudaMemcpy(fft_in, host_in, N * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
            cuDoubleComplex a, b;
            a.x = 1;a.y = 0;b.x = 0;b.y = 0;

            for(int i = 0; i < N; i++){
                cublasZgemv (
                    handle,    
                    CUBLAS_OP_T,    
                    N,    
                    N,    
                    &a,    
                    d_A,    
                    N,    
                    fft_in,    
                    1,    
                    &b,    
                    fft_out,    
                    1   
                );
            }
            // print<<<1,1>>>(fft_out);
            // product<<<N/1024,1024>>>(fft_out);
            ntt_in = rns.decompose(fft_out,N);
            // print<<<1,1>>>(ntt_in);
            // for(int i = 0; i < N * 8; i += N){
            //     print_d<<<1,1>>>(ntt_in,i);
            // }
            for(int i = 0; i < decomposeSize; i++){
                forwardNTT(ntt_in+N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
            }
            // print<<<1,1>>>(ntt_in);
            // print_d<<<1,1>>>(ntt_in,0);
            // print_d<<<1,1>>>(ntt_in,2048);
            // print_d<<<1,1>>>(ntt_in,4096);

            return ntt_in;
        }
        double* decode(unsigned long long* encodeVec){
            cuDoubleComplex *fft_in,*host_in;
            cuDoubleComplex *fft_out,*fft_out_t;
            unsigned long long *ntt_in;
            cudaStream_t ntt = 0;

            Check(cudaMallocHost((void**)&host_in, n * sizeof(cuDoubleComplex)));
            Check(cudaMalloc((void**)&ntt_in, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&fft_in, N * sizeof(cuDoubleComplex)));
            Check(cudaMalloc((void**)&fft_out, N * sizeof(cuDoubleComplex)));
            Check(cudaMalloc((void**)&fft_out_t, N * sizeof(cuDoubleComplex)));
            
            for(int i = 0; i < decomposeSize; i++){
                inverseNTT(encodeVec + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
            }
            // print<<<1,1>>>(encodeVec);
            // for(int i = 0; i < N * decomposeSize; i += N){
            //     print_d<<<1,1>>>(encodeVec,i);
            // }
            // fft_out_t = rns.compose(encodeVec,N,1);
            fft_out = rns.compose(encodeVec,N,1);
            // print<<<1,1>>>(fft_out);
            // print<<<1,1>>>(fft_out_t);            
            cuDoubleComplex a, b;
            a.x = 1;a.y = 0;b.x = 0;b.y = 0;

            for(int i = 0; i < N; i++){
                cublasZgemv (
                    handleinv,    
                    CUBLAS_OP_N,    
                    N,    
                    N,    
                    &a,    
                    inv_A,   
                    N,    
                    fft_out,   
                    1,    
                    &b,    
                    fft_in,   
                    1   
                );
            }
            
            Check(cudaMemcpy(host_in, fft_in, n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
            double *res;
            res = (double *)malloc(n * sizeof(double));
            for(int i = 0; i < n ; i++){
                res[i] = host_in[i].x / scale;
            }
            return res;
        }
        void test(){
            unsigned long long *temp1;
            unsigned long long *temp2;
            cudaStream_t ntt = 0;
            Check(cudaMallocHost((void**)&temp1, N * sizeof(unsigned long long)));
            // Check(cudaMallocHost((void**)&temp2, N * sizeof(unsigned long long)));
            // for(int i = 0; i < 8; i++){
            //     printf("%llu\n",q[i]);
            // }
            genRandom<<<N/1024,1024>>>(temp1,114514);
            // genRandom<<<N/1024,1024>>>(temp2,2);
            // print<<<1,1>>>(temp1);
            // print<<<1,1>>>(temp2);
            // rns.decompose(temp1,N);
            // // rns.decompose(temp2,N);
            // print<<<1,1>>>(temp1);
            // print<<<1,1>>>(temp2);
            forwardNTT(temp1 ,N,ntt,q[0],mu[0],q_bit[0],psiTable[0]);
            // for(int i = 0; i < decomposeSize; i++){
            //     forwardNTT(temp1 + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
            // }
            print<<<1,1>>>(temp1);
            // for(int i = 0; i < decomposeSize; i++){
            //     // forwardNTT(temp2 + N * i,N,ntt,q[i],mu[i],q_bit[i],psiTable[i]);
            // }
            // print<<<1,1>>>(temp1);
            // print<<<1,1>>>(temp2);
            // for(int i = 0; i < decomposeSize; i++){
            //     // polymul<<<N/1024,1024>>>(temp1 + N * i,temp2 + N * i,temp1 + N * i,q[i],mu[i],q_bit[i]);
            // }
            // print<<<1,1>>>(temp1);
            // for(int i = 0; i < decomposeSize; i++){
            //     inverseNTT(temp1 + N * i,N,ntt,q[i],mu[i],q_bit[i],psiinvTable[i]);
            // }
            // print<<<1,1>>>(temp1);
            // print<<<1,1>>>(temp2);

        }
};