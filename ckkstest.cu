#include "encoder.cuh"
#include"cufft.h"
// __global__ void print(unsigned long long* a){
//     printf("%llu\n",a[0]);
// }
int main(){
    int N = 4096;
    double a[N];
    for(int i = 0; i < N; i++){
        a[i] = 1;
    }


    Encoder encoder(N,100);
    auto encodeVec = encoder.encode(a);
    // print<<<1,1>>>(d_a);
    auto plain = encoder.decode(encodeVec);
    // print<<<1,1>>>(d_a);
    for(int i = 0; i < 10; i++){
        printf("%lf\n",plain[i]);
    }
    cudaDeviceSynchronize();
}