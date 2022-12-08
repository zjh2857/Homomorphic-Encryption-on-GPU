#pragma once
#include "parameter.cuh"
#include "uint128.cuh"
#include "ntt_60bit.cuh"
#include "cuda_runtime_api.h"
#include "device_launch_parameters.h"
#include "cufft.h"
#include <iostream>

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
__global__ void delta(cufftDoubleReal* in,unsigned long long* out,double scale){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    out[tid] = (unsigned long long)(in[tid] * scale); 
}
__global__ void deltainv(unsigned long long* in,cufftDoubleReal* out,double scale){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    out[tid] = (double)(in[tid] * 1.0 / scale); 
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
        cufftHandle cufftForwrdHandle, cufftInverseHandle;
        Encoder(int n,double scale){
            this->n = n;
            N = n * 2;
            // this->N = N;
            this->scale = scale;
            getParams(q, psi, psiinv, ninv, q_bit, N);
            cudaMalloc(&psiTable, N * sizeof(unsigned long long));
            cudaMalloc(&psiinvTable, N * sizeof(unsigned long long));
            fillTablePsi128<<<N/1024,1024>>>(psi, q, psiinv, psiTable, psiinvTable, log2(N));

            uint128_t mu1 = uint128_t::exp2(q_bit * 2);
            mu = (mu1 / q).low;

            cufftPlan1d(&cufftForwrdHandle, N, CUFFT_D2Z, BATCH);
	        cufftPlan1d(&cufftInverseHandle, N, CUFFT_Z2D, BATCH);
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
            cufftDoubleComplex *fft_in,*host_in;
            cufftDoubleReal *fft_out;
            unsigned long long *ntt_in;
            cudaStream_t ntt = 0;

            Check(cudaMallocHost((void**)&host_in, n * sizeof(cufftDoubleComplex)));
            Check(cudaMalloc((void**)&ntt_in, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&fft_in, n * sizeof(cufftDoubleComplex)));
            Check(cudaMalloc((void**)&fft_out, N * sizeof(cufftDoubleReal)));
            for(int i = 0; i < n; i++){
                host_in[i].x = plainVec[i];
                host_in[i].y = 0;
            }

            
            Check(cudaMemcpy(fft_in, host_in, n * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice));
            cufftExecZ2D(cufftInverseHandle, fft_in, fft_out);
            delta<<<N/1024,1024>>>(fft_out,ntt_in,scale);
            forwardNTT(ntt_in,N,ntt,q,mu,q_bit,psiTable);
            return ntt_in;
        }
        double* decode(unsigned long long* encodeVec){
            cufftDoubleComplex *fft_in,*host_in;
            cufftDoubleReal *fft_out;
            unsigned long long *ntt_in;
            cudaStream_t ntt = 0;

            Check(cudaMallocHost((void**)&host_in, n * sizeof(cufftDoubleComplex)));
            Check(cudaMalloc((void**)&ntt_in, N * sizeof(unsigned long long)));
            Check(cudaMalloc((void**)&fft_in, n * sizeof(cufftDoubleComplex)));
            Check(cudaMalloc((void**)&fft_out, N * sizeof(cufftDoubleReal)));


            inverseNTT(encodeVec,N,ntt,q,mu,q_bit,psiinvTable);
            deltainv<<<N/1024,1024>>>(encodeVec,fft_out,scale);

            cufftExecD2Z(cufftInverseHandle, fft_out, fft_in);

            Check(cudaMemcpy(host_in, fft_in, n * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost));
            double *res;
            res = (double *)malloc(n * sizeof(double));
            for(int i = 0; i < n ; i++){
                res[i] = host_in[i].x;
            }
            return res;
            // return encodeVec;
        }
};