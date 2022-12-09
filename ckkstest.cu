#include "encoder.cuh"
#include"cufft.h"
#include "encryptor.cuh"
__global__ void print(unsigned long long* a){
    for(int i = 0; i < 8;i++){
        printf("%llu ",a[i]);
    }printf("\n");
}
int main(){
    int N = 4096;
    double scale = 10000;
    double a[N];
    for(int i = 0; i < N; i++){
        a[i] = 20000;
    }

// a[0] = 2;
    keyGen keygen(N,scale);
    Encoder encoder(N,scale);
    auto encodeVec = encoder.encode(a);

// print<<<1,1>>>(encodeVec);

    Encryptor encryptor(N,scale);
    auto ciptertext = encryptor.encrypt(encodeVec,keygen.pub);
    unsigned long long* dec = encryptor.decrypt(ciptertext,keygen.pri);
    auto plain = encoder.decode(dec);
    print<<<1,1>>>(dec);


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