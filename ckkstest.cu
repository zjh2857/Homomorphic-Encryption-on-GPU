#include "encoder.cuh"
#include"cufft.h"
#include "encryptor.cuh"
#include "evaluator.cuh"

__global__ void print(unsigned long long* a){
    for(int i = 0; i < 8;i++){
        printf("%llu ",a[i]);
    }printf("\n");
}
__global__ void print_d(unsigned long long* a,int d){
    for(int i = 0; i < 1;i++){
        printf("printf_d%d:%llu ",d,a[d]);
    }printf("\n");
}
int main(){
    int N = 2048;
    double scale = 1 << 20;
    double a[N];
    for(int i = 0; i < N; i++){
        a[i] = i;
    }
    double b[N];
    for(int i = 0; i < N; i++){
        b[i] = 0;
    }
    b[0] = 100;
    keyGen keygen(N,scale);
    Encoder encoder(N,scale,8);
    // Encryptor encryptor(N,scale);
    Evaluator evaluator(N);

    auto encodeVeca = encoder.encode(a);
    auto encodeVecb = encoder.encode(b);
    // print<<<1,1>>>(encodeVeca);


    // auto ciptertexta = encryptor.encrypt(encodeVeca,keygen.pub);
    // // print<<<1,1>>>(ciptertexta.a);
    // // print<<<1,1>>>(ciptertexta.b);

    // auto ciptertextb = encryptor.encrypt(encodeVecb,keygen.pub);
    // printf("%p,%p\n",ciptertexta,ciptertextb);
    //     print<<<1,1>>>(ciptertexta.a);
    // print<<<1,1>>>(ciptertexta.b);
    // print<<<1,1>>>(ciptertextb.a);
    // print<<<1,1>>>(ciptertextb.b);
    // auto plaina = encoder.decode(encodeVeca);
    evaluator.mulPlain(encodeVeca,encodeVecb);
    
    // print<<<1,1>>>(ciptertexta.a);
    // print<<<1,1>>>(ciptertexta.b);    
    // encryptor.decrypt(ciptertextb,keygen.pri);
    // unsigned long long* dec = encryptor.decrypt(ciptertexta,keygen.pri);
    auto plaina = encoder.decode(encodeVeca);
    // print<<<1,1>>>(dec);
    

        for(int i = 0; i < 10; i++){
            printf("%lf\n",plaina[i]);
        }
    //  auto plainb = encoder.decode(encodeVecb);
    // // print<<<1,1>>>(dec);
    

    //     for(int i = 0; i < 10; i++){
    //         printf("%lf\n",plainb[i]);
    //     }
    cudaDeviceSynchronize();
}