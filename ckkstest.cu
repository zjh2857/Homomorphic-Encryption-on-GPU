#include "encoder.cuh"
#include"cufft.h"
#include "encryptor.cuh"
#include "evaluator.cuh"
__global__ void print(unsigned long long* a){
    for(int i = 0; i < 8;i++){
        printf("%llu ",a[i]);
    }printf("\n");
}
int main(){
    int N = 4096;
    double scale = 100000;
    double a[N];
    for(int i = 0; i < N; i++){
        a[i] = 1000;
    }
    double b[N];
    for(int i = 0; i < N; i++){
        b[i] = 1000;
    }
// a[0] = 2;
    keyGen keygen(N,scale);
    Encoder encoder(N,scale);
    Encryptor encryptor(N,scale);
    Evaluator evaluator(N);

    auto encodeVeca = encoder.encode(a);
    auto encodeVecb = encoder.encode(b);
    // print<<<1,1>>>(encodeVeca);


    auto ciptertexta = encryptor.encrypt(encodeVeca,keygen.pub);
    // print<<<1,1>>>(ciptertexta.a);
    // print<<<1,1>>>(ciptertexta.b);

    auto ciptertextb = encryptor.encrypt(encodeVecb,keygen.pub);
    // printf("%p,%p\n",ciptertexta,ciptertextb);
    //     print<<<1,1>>>(ciptertexta.a);
    // print<<<1,1>>>(ciptertexta.b);
    // print<<<1,1>>>(ciptertextb.a);
    // print<<<1,1>>>(ciptertextb.b);
    evaluator.addcipter(ciptertexta,ciptertextb);
    
    // print<<<1,1>>>(ciptertexta.a);
    // print<<<1,1>>>(ciptertexta.b);    
    // encryptor.decrypt(ciptertextb,keygen.pri);
    unsigned long long* dec = encryptor.decrypt(ciptertexta,keygen.pri);
    auto plain = encoder.decode(dec);
    // print<<<1,1>>>(dec);
    

        for(int i = 0; i < 10; i++){
            printf("%lf\n",plain[i]);
        }
 
        cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
        // Possibly: exit(-1) if program cannot continue....
    }
}