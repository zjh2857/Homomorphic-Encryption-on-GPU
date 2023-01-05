#include "encoder.cuh"
#include"cufft.h"
#include "encryptor.cuh"
#include "evaluator.cuh"
#include "freshman.h"
const int testTimes = 100;

__global__ void print(unsigned long long* a){
    for(int i = 0; i < 8;i++){
        printf("%llu,",a[i]);
    }printf("\n");
}
__global__ void print(unsigned long long* a,unsigned long long id){
    printf("%llu\n",id);
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
    printf("===\n");
    int N = 2048;
    double scale = 1llu << 30;
    double a[N];
    for(int i = 0; i < N; i++){
        a[i] = i ;
    }
    // a[N/2-1] = 10000;
    // a[0] = 10;
    double b[N];
    for(int i = 0; i < N; i++){
        b[i] = i;
    }

    keyGen keygen(N,scale,8);
    Encoder encoder(N,scale,8);
    Encryptor encryptor(N,scale,8);
    Evaluator evaluator(N,8);
    // double encodeTimes = cpuSecond();
    // for(int i = 0; i < testTimes; i++){
    //     auto foo = encoder.encode(a);
    // }
    // cudaDeviceSynchronize();
    // printf("encode Times:%lf microseconds",(cpuSecond() - encodeTimes)/testTimes);
    auto encodeVeca = encoder.encode(a);
    auto encodeVecb = encoder.encode(b);
    double encryptoTimes = cpuSecond();
    for(int i = 0; i < testTimes; i++){
        auto foo = encryptor.encrypt(encodeVeca,keygen.pub);
    }
    cudaDeviceSynchronize();
    printf("encrypto Times:%lf seconds",(cpuSecond() - encryptoTimes)/testTimes);
    auto ciptertexta = encryptor.encrypt(encodeVeca,keygen.pub);
    auto ciptertextb = encryptor.encrypt(encodeVecb,keygen.pub);
    auto ciptertextc = evaluator.mulcipter(ciptertexta,ciptertextb);
    auto ciptertextd = evaluator.relien(ciptertextc,keygen.relien);
    evaluator.rescale(ciptertextd);
    unsigned long long* dec = encryptor.decrypt(ciptertextd,keygen.pri);
    auto plaina = encoder.decode(dec);
    for(int i = 0; i < 20; i++){
        printf("%lf\n",plaina[i] * 1179649 / scale);
    }
    //  auto plainb = encoder.decode(encodeVecb);
    // // // print<<<1,1>>>(dec);
    

    // //     for(int i = 0; i < 10; i++){
    // //         printf("%lf\n",plainb[i]);
    // //     }
    cudaDeviceSynchronize();
}