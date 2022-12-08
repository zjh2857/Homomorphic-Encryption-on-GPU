#include <cufft.h>
#include <cstdio>
#include <iostream>
#include"cuda_runtime_api.h"
#include"device_launch_parameters.h"
#define NX 1024
#define BATCH 1
using namespace std;
#define Check(status)\
{\
	if (status != cudaSuccess)\
	{\
		cout << "行號:" << __LINE__ << endl;\
		cout << "錯誤:" << cudaGetErrorString(status) << endl;\
	}\
}
int main(){
    cufftComplex *data_h;
    cufftComplex *data;
    Check(cudaMallocHost(&data_h,NX*sizeof(cufftComplex)));
    Check(cudaMalloc(&data,NX * sizeof(cufftComplex)));
    cufftComplex *res = (cufftComplex *)malloc(NX * sizeof(cufftComplex));
    for(int i =0; i < NX/2; i++){
        data_h[i].x = (i * 113 + 71) % 5 * 1.41;
        data_h[i].y = 0;    
    }
    for(int i =0; i < NX/2; i++){
        data_h[NX-i-1].x = (i * 113 + 71) % 5 * 1.41;
        data_h[NX-i-1].y = 0;
    }
    cudaMemcpy(data,data_h,NX * sizeof(cufftComplex),cudaMemcpyHostToDevice);
    cufftHandle plan;
    cufftPlan1d(&plan, NX, CUFFT_C2C,BATCH);

    cufftExecC2C(plan, data, data, CUFFT_FORWARD);

    cufftDestroy(plan);
    Check(cudaMemcpy(res,data,NX * sizeof(cufftComplex),cudaMemcpyDeviceToHost));
    cudaFree(data);
    for(int i = 0;i < 10;i++){
        printf("%lf,%lf\n",res[i].x,res[i].y);
    }
    return 0;
}
