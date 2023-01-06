#include "encoder.cuh"
#include"cufft.h"
#include "encryptor.cuh"
#include "evaluator.cuh"

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
int main(int argc,char *argv[]){
    printf("===\n");
    int N = 2048;
    double scale = 1llu << 30;
    double a[N];
    for(int i = 0; i < N; i++){
        // if(i&1){
        //     a[i] = 1;
        // }
        // else{
        //     a[i] = -1;
        // }
        a[i] = N - i - 1;
    }
    // a[N/2-1] = 10000;
    // a[0] = 10;
    double b[N];
    for(int i = 0; i < N; i++){
        b[i] = i;
    }
    // b[N-1] = a[0];

    keyGen keygen(N,scale,8);
    Encoder encoder(N,scale,8);
    Encryptor encryptor(N,scale,8);
    Evaluator evaluator(N,8);
    // encoder.test();
    auto encodeVeca = encoder.encode(a);
    // auto ciptertexta = encryptor.encrypt(encodeVeca,keygen.pub);

    // int t = atoi(argv[1]);
    // for(int i = 0; i < t;i++){
    //     evaluator.rotation(encodeVeca);
    // }
    // auto ciptertextr = evaluator.sum(ciptertexta,keygen.galois);
    // evaluator.rotation(encodeVeca);
    // for(int i = 0;  i < N; i++){
    //     evaluator.rotation(encodeVeca);
    //     auto plaina = encoder.decode(encodeVeca);
    //     printf("%lf\n",plaina[0]);
    // }
    // evaluator.rotation(encodeVeca);

    // evaluator.rescale(encodeVeca);
    auto encodeVecb = encoder.encode(b);
    // // print<<<1,1>>>(encodeVeca);
    // // print<<<1,1>>>(encodeVecb);
    // // cudaDeviceSynchronize();


    auto ciptertexta = encryptor.encrypt(encodeVeca,keygen.pub);
    // auto ciptertextr = evaluator.sum(ciptertexta,keygen.galois);
    // evaluator.rescale(ciptertexta);
    // // print<<<1,1>>>(keygen.pub.a);
    // // print<<<1,1>>>(keygen.pub.b);
    auto ciptertextb = encryptor.encrypt(encodeVecb,keygen.pub);
    auto ciptertextr = evaluator.innerProduct(ciptertexta,ciptertextb,keygen.galois,keygen.relien);
    // // printf("%p,%p\n",ciptertexta,ciptertextb);
    // //     print<<<1,1>>>(ciptertexta.a);
    // // print<<<1,1>>>(ciptertexta.b);
    // // print<<<1,1>>>(ciptertextb.a);
    // // print<<<1,1>>>(ciptertextb.b);
    // // auto plaina = encoder.decode(encodeVeca);
    // print<<<1,1>>>(encodeVecb);

    // auto ciptertextc = evaluator.mulcipter(ciptertexta,ciptertextb);

    // auto ciptertextd = evaluator.relien(ciptertextc,keygen.relien);
    // evaluator.rescale(ciptertextd);
    // print<<<1,1>>>(ciptertextd.a);
    // print<<<1,1>>>(ciptertextd.b);
    // print<<<1,1>>>(ciptertextd.c);
    // evaluator.mulPlain(ciptertexta,encodeVecb);
    // // evaluator.mulPlain(encodeVeca,encodeVecb);

    // // printf("###\n");
    // // print<<<1,1>>>(ciptertexta.a);
    // // print<<<1,1>>>(ciptertexta.b);    
    // encryptor.decrypt(ciptertextb,keygen.pri);
    unsigned long long* dec = encryptor.decrypt(ciptertextr,keygen.pri);
    // print<<<1,1>>>(dec);
    // // print<<<1,1>>>(dec);
    auto plaina = encoder.decode(dec);
    // auto plaina = encoder.decode(encodeVecd);
    // // // auto plainb =  encoder.decode(encodeVecb);
    // // // print<<<1,1>>>(dec);
    
    // printf("%lf",plaina[0]);
    for(int i = 0; i < 10; i++){
        printf("%d,%lf\n",i,plaina[i]/scale);
    }
    //  auto plainb = encoder.decode(encodeVecb);
    // // // print<<<1,1>>>(dec);
    

    // //     for(int i = 0; i < 10; i++){
    // //         printf("%lf\n",plainb[i]);
    // //     }
    cudaDeviceSynchronize();
}