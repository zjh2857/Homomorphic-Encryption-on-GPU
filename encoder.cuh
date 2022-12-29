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
__global__ void delta(cuDoubleComplex* in,unsigned long long* out,double scale,unsigned long long q,int N){

    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    long long t = in[tid].x / N  * scale;
    t = (t+q)%q;
    out[tid] = t; 

    // if(tid < 10 ){
    //     printf("!%d:%lf,%llu\n",tid,in[tid],out[tid]);
    // }

}
__global__ void zerocnt(unsigned long long* a,int N){
    int cnt = 0;
    for(int i = 0; i < N; i++){
        if(a[i] != 0){
            printf("%d,%llu\t",i,a[i]);
            cnt++;
        }
    }printf("\n");
    printf("%d\n",cnt);
}
__global__ void eqcnt(unsigned long long* a,int N){
    int cnt = 0;
    for(int i = 1; i < N; i++){
        if(a[i] != a[i-1]){
            cnt++;
            printf("%d,%llu\n",i,a[i]);
        }
    }printf("\n");
    printf("%d\n",cnt);
}
__global__ void deltainv(unsigned long long* in,cuDoubleComplex* out,double scale,unsigned long long q){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    long long t = in[tid];
    if(t > q/2){
        t -= q;
    }
    out[tid].x = (double)(t * 1.0 / scale);
    out[tid].y = 0;
    // if(out[tid] > q/2){
    //     out[tid] -= q;
    // } 

    // if(tid < 10 ){
    //     printf("$%d:%llu,%lf\n",tid,in[tid],out[tid]);
    // }
    // if(tid == 0){
    //     printf("$%lf",out[tid]);
    // }
}
__global__ void print(unsigned long long* a);
__global__ void print_d(unsigned long long* a,int d);

//     printf("%llu\n",a[0]);
// }
__global__ void print(cuDoubleComplex* h_A){
    for(int i = 0;i < 8; i++){
        printf("%lf+%lfi\n",h_A[i].x,h_A[i].y);
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
class Encoder{
    public:
        int n;
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
            getParams(q, psi, psiinv, ninv, q_bit, N);
            printf("qqqqqqqqqqq%llu\n",q);
            cudaMalloc(&psiTable, N * sizeof(unsigned long long));
            cudaMalloc(&psiinvTable, N * sizeof(unsigned long long));

            fillTablePsi128<<<N/1024,1024>>>(psi, q, psiinv, psiTable, psiinvTable, log2(N));

            uint128_t mu1 = uint128_t::exp2(q_bit * 2);
            mu = (mu1 / q).low;

            cublasCreate(&handle);
            cublasCreate(&handleinv);
            // 在 显存 中为将要计算的矩阵开辟空间
            cudaMalloc (
                (void**)&d_A,    // 指向开辟的空间的指针
                N*N * sizeof(cuDoubleComplex)    //　需要开辟空间的字节数
            );
            cudaMalloc (
                (void**)&inv_A,    // 指向开辟的空间的指针
                N*N * sizeof(cuDoubleComplex)    //　需要开辟空间的字节数
            );
            // 在 显存 中为将要存放运算结果的矩阵开辟空间
            initfft<<<N*N/1024,1024>>>(d_A,N);
            initfftinv<<<N*N/1024,1024>>>(inv_A,N);
            // cufftPlan1d(&cufftForwrdHandle, N, CUFFT_D2Z, BATCH);
	        // cufftPlan1d(&cufftInverseHandle, N, CUFFT_Z2D, BATCH);
        };
        // void encode(unsigned long long* plainVec){
        //     cudaStream_t ntt = 0;
        //     forwardNTT(plainVec,N,ntt,q,mu,q_bit,psiTable);
        //     // return plainVec;
        // }
        // void decode(unsigned long long* encodeVec){
        //     cudaStream_t ntt = 0;
        //     inverseNTT(encodeVec,N,ntt,q,mu,q_bit,psiinvTable);
        //     // return encodeVec;
        // }
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

            // printf("%p,%p\n",fft_in,fft_out);
            for(int i = 0; i < N; i++){
                cublasZgemv (
                    handle,    // blas 库对象
                    CUBLAS_OP_T,    // 矩阵 A 属性参数
                    N,    // A, C 的行数
                    N,    // B, C 的列数
                    &a,    // 运算式的 α 值
                    d_A,    // A 在显存中的地址
                    N,    // lda
                    fft_in,    // B 在显存中的地址
                    1,    // ldb
                    &b,    // 运算式的 β 值
                    fft_out,    // C 在显存中的地址(结果矩阵)
                    1    // ldc
                );
            }
            // print<<<1,1>>>(fft_out);
            delta<<<N/1024,1024>>>(fft_out,ntt_in,scale,q,N);
            print<<<1,1>>>(ntt_in);
            for(int i = 0; i < 16; i++)
            print_d<<<1,1>>>(ntt_in,i * (N/16));
            // print<<<1,1>>>(ntt_in);
            // print_d<<<1,1>>>(ntt_in,1967);
            ntt_in = rns.decompose(ntt_in,N);
            // print<<<1,1>>>(ntt_in);
            // for(int i = 0; i < decomposeSize; i++){
            //     print_d<<<1,1>>>(ntt_in,1967 + i * N);
            // }
            // ntt_in = rns.compose(ntt_in,N);
            // cudaDeviceSynchronize();
            // exit(1);
            // ntt_in = rns.decompose(ntt_in,N);

            // for(int i = 0; i < decomposeSize; i++){
            //     print_d<<<1,1>>>(ntt_in,1967 + i * N);
            // }
// zerocnt<<<1,1>>>(ntt_in,N*decomposeSize);
            for(int i = 0; i < decomposeSize; i++){
                forwardNTT(ntt_in+N * i,N,ntt,q,mu,q_bit,psiTable);
            }
            return ntt_in;
        }
        double* decode(unsigned long long* encodeVec){
            cuDoubleComplex *fft_in,*host_in;
            cuDoubleComplex *fft_out;
            unsigned long long *ntt_in;
            cudaStream_t ntt = 0;

            Check(cudaMallocHost((void**)&host_in, n * sizeof(cuDoubleComplex)));
            Check(cudaMalloc((void**)&ntt_in, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&fft_in, N * sizeof(cuDoubleComplex)));
            Check(cudaMalloc((void**)&fft_out, N * sizeof(cuDoubleComplex)));
            // eqcnt<<<1,1>>>(encodeVec,N*decomposeSize);

            for(int i = 0; i < decomposeSize; i++){
                // printf("%p\n",encodeVec + N * i);
                inverseNTT(encodeVec + N * i,N,ntt,q,mu,q_bit,psiinvTable);
            }
            // print<<<1,1>>>(encodeVec);
            // zerocnt<<<1,1>>>(encodeVec,N*decomposeSize);
            // for(int i = 0; i < decomposeSize; i++){
            //     print_d<<<1,1>>>(encodeVec,1967 + i * N);
            // }
            // print<<<1,1>>>(encodeVec);
            ntt_in = rns.compose(encodeVec,N);
            print<<<1,1>>>(ntt_in);
            // for(int i = 0; i < 16; i++)print_d<<<1,1>>>(ntt_in,i * (N/16));
            deltainv<<<N/1024,1024>>>(ntt_in,fft_out,scale,q);

            cuDoubleComplex a, b;
            a.x = 1;a.y = 0;b.x = 0;b.y = 0;

            for(int i = 0; i < N; i++){
                cublasZgemv (
                    handleinv,    // blas 库对象
                    CUBLAS_OP_N,    // 矩阵 A 属性参数
                    N,    // A, C 的行数
                    N,    // B, C 的列数
                    &a,    // 运算式的 α 值
                    inv_A,    // A 在显存中的地址
                    N,    // lda
                    fft_out,    // B 在显存中的地址
                    1,    // ldb
                    &b,    // 运算式的 β 值
                    fft_in,    // C 在显存中的地址(结果矩阵)
                    1    // ldc
                );
            }
            
            Check(cudaMemcpy(host_in, fft_in, n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
            double *res;
            res = (double *)malloc(n * sizeof(double));
            for(int i = 0; i < n ; i++){
                res[i] = host_in[i].x ;
            }
            return res;
            // return encodeVec;
        }
};