#include <cstdio>
#include <stdlib.h>
#include "helper.cuh"
using namespace std;

struct uint128{
    unsigned long long low;
    unsigned long long high;
};
// __global__ void print1(unsigned long long* a){
//     for(int i = 0; i < 8; i++){
//         printf("%llu\t",a[i]);
//     }printf("\n");
// }
__constant__ double half_d;
__device__ __host__ __forceinline__ void bigMul(unsigned long long a,unsigned long long b, uint128& res){
    u_int64_t a0 = a & 0xffffffff;
    u_int64_t a1 = a >> 32;
    u_int64_t b0 = b & 0xffffffff;
    u_int64_t b1 = b >> 32;
    u_int64_t low = a0 * b0 ;
    u_int64_t carry = 0;
    if(low > low + ((a0 * b1) << 32llu)){
        carry+=1;
    }
    low = low + ((a0 * b1) << 32llu);
    if(low > low + ((a1 * b0) << 32llu)){
        carry+=1;
    }
    low = low + ((a1 * b0) << 32llu);

    u_int64_t high = a1 * b1 + ((a0 * b1)>> 32llu) + ((a1 * b0)>> 32llu) + carry;
    // high += ((a0 * b0) >> 32) + 
    res.low = low;
    res.high = high;
    // return res;
}
__device__ __host__ __forceinline__ void bigIntegerMul(unsigned long long * a, unsigned long long b,int size){
    uint128 temp; 
    bigMul(a[0],b,temp);
    unsigned long long carry = temp.high;
    a[0] = temp.low;

    for(int i = 1; i < size; i++){
        if(!a[i] && !carry){
            break;
        }
        uint128 temp;
        bigMul(a[i],b,temp);
        // temp = te
        a[i] = temp.low;
        // a[i] = (temp & 0xffffffffffffffff);
        if(a[i] + carry < a[i]){
            a[i] += carry;
            carry = 1 + temp.high;
        }
        else{
            a[i] += carry;
            carry = temp.high;
        }
    }
}

__device__ __host__ __forceinline__ void bigIntegerAdd(unsigned long long * a, unsigned long long *b,int size){
    unsigned long long carry = 0;
    for(int i = 0; i < size; i++){
        if(a[i] + b[i] + carry < a[i]){
            // carry = 1;
            a[i] = a[i] + b[i] + carry;
            carry = 1;
        }
        else{
            a[i] = a[i] + b[i] + carry;
            carry = 0;
        }
    }
}
__device__ __host__ __forceinline__ int isneg(unsigned long long * a, unsigned long long *b ,unsigned long long p,int size){
    // unsigned long long* temp = (unsigned long long*)malloc(size * sizeof(unsigned long long)); 
    unsigned long long temp[64];
    for(int i = 0; i < size; i++){
        temp[i] = b[i];
    }

    bigIntegerMul(temp,p,size);
    unsigned long long borrow = 0;
    for(int i = 0; i < size; i++){
        // if(p == 1742){
        //     printf("&&%d,%llu\n",i,borrow);
        // }
        if(a[i] < borrow){
            borrow = 1;
            temp[i] = a[i] - borrow - temp[i];
            continue;
        }

        if(a[i] - borrow >= temp[i]){
            temp[i] = a[i] - borrow - temp[i];
            borrow = 0;
        }else{
            temp[i] = a[i] - borrow - temp[i];
            borrow = 1;
        }
        
    }
    if(borrow == 1){
        return -1;
    }
    return 1;
}


__device__ __host__ __forceinline__ double bigInteger2udouble(unsigned long long *a,unsigned long long *b,unsigned long long p,int size,int tid){
    // unsigned long long* temp = (unsigned long long*)malloc(size * sizeof(unsigned long long)); 
    unsigned long long temp[64];
    for(int i = 0; i < size; i++){
        temp[i] = b[i];
    }

    bigIntegerMul(temp,p,size);

    unsigned long long borrow = 0;
    for(int i = 0; i < size; i++){

        if(a[i] < borrow){
            borrow = 1;
            temp[i] = a[i] - borrow - temp[i];
            continue;
        }

        if(a[i] - borrow >= temp[i]){
            temp[i] = a[i] - borrow - temp[i];
            borrow = 0;
        }else{
            temp[i] = a[i] - borrow - temp[i];
            borrow = 1;
        }
        
    }

    // if(tid==0){
    //     for(int i = 0; i < size; i++){
    //         printf("%llu * (2 ** %d) + \t",temp[i],i * 64);
    //     }printf("\n");
    //     // printf("%lf,%lf\n",res_n,res_p);
    // }

    double res_p = 0;
    double res = 0;
    double res_n = 0;
    double twopow64 = 18446744073709551616.0;
    double base = 1.0;
    // if(tid==0){
    //     for(int i = 0; i < size; i++){
    //         printf("%llu!\t",temp[i]);
    //     }printf("\n");
    //     // printf("%lf,%lf\n",res_n,res_p);
    // }
    for(int i = 0; i < size; i++){
        res_p += temp[i] * base;
        base*= twopow64;
    }
    for(int i = 0; i < size; i++){
        temp[i] = b[i];
    }

    bigIntegerMul(temp,p+1,size);
    borrow = 0;
    for(int i = 0; i < size; i++){

        if(a[i] < borrow){
            borrow = 1;
            temp[i] = a[i] - borrow - temp[i];
            continue;
        }

        if(a[i] - borrow >= temp[i]){
            temp[i] = a[i] - borrow - temp[i];
            borrow = 0;
        }else{
            temp[i] = a[i] - borrow - temp[i];
            borrow = 1;
        }
        
    }
    for(int i = 0; i < size; i++){
        temp[i] = ~temp[i];
    }
    unsigned long long carry = 1;
    for(int i = 0; i < size; i++){
        temp[i] += carry;
        if(temp[i] == 0 && carry){
            carry = 1;
        }else{
            carry = 0;
        }
    }
    base = 1.0;
    // if(tid==0){
    //     for(int i = 0; i < size; i++){
    //         printf("%llu * (2 ** %d) + \t",temp[i],i * 64);
    //     }printf("\n");
    //     // printf("%lf,%lf\n",res_n,res_p);
    // }
    for(int i = 0; i < size; i++){
        res_n += temp[i] * base;
        base*= twopow64;
    }


    if(res_n < res_p){
        res = -res_n;
    }
    else{
        res = res_p;
    }

    return res;
}
__device__ __host__ __forceinline__ double bigIntegerMod(unsigned long long * a, unsigned long long *b,int size,int tid){
    unsigned long long l = 0;
    unsigned long long r = (1llu << 63);
    unsigned long long guess;


    int cnt = 0;
    unsigned long long p = 0;
    // if(tid==0){
    //     for(int i = 0; i < 8;i++){
    //         printf("%llu\t",a[i]);
    //     }printf("\n");
    //     for(int i = 0; i < 8;i++){
    //         printf("%llu\t",b[i]);
    //     }printf("\n");
    // }
    while(l < r){
        if(cnt++ > 100){
            printf("114514\n");
            return 114514;
        }
        unsigned long long guess = (l + r)/2;
        // if(tid == 114514){
        //     printf("%llu,%llu,%llu\n",l,r,(l+r)/2);
        // }
        int res = isneg(a,b,guess,size);
        if(res == 1){
            p = guess;
            l = guess + 1;

        }
        else if(res == -1){
            r = guess;
        }
    }
    // if(tid == 0){
    //     printf("%llu!!\n",p);
    // }
    return bigInteger2udouble(a,b,p,size,tid);

}


__global__ void cudadecompose(cuDoubleComplex *list,unsigned long long* moduleChain,int listLen,int moduleLen,unsigned long long * decomposeList,unsigned long long scale){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if(list[tid].x < 0){
        unsigned long long temp = -list[tid].x / listLen *scale;
        for(int i = 0; i < moduleLen; i++){
            decomposeList[i * listLen + tid] =moduleChain[i] - (temp % moduleChain[i]);
        }        
    }
    else{
        unsigned long long temp = list[tid].x / listLen *scale;
        for(int i = 0; i < moduleLen; i++){
            decomposeList[i * listLen + tid] = (temp) % moduleChain[i];
        }
    }
}

__global__ void cudadecompose(unsigned long long *list,unsigned long long* moduleChain,int listLen,int moduleLen,unsigned long long * decomposeList){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned long long temp = list[tid];
    for(int i = 0; i < moduleLen; i++){
        decomposeList[i * listLen + tid] = temp % moduleChain[i];
        // if(tid == 0){
        //     printf("%llu,%llu,%llu\n",list[tid],moduleChain[i],decomposeList[i * listLen + tid]);
        // }
    }
}

__global__ void cudacompose(unsigned long long *decomposeList,
                            unsigned long long* moduleChain,
                            int listLen,
                            int moduleLen,
                            unsigned long long* Ni,
                            unsigned long long *bigN,cuDoubleComplex * composeList,unsigned long long* temp1,unsigned long long* temp2,unsigned long long scale){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;

    for(int i = 0; i < moduleLen; i++){
        temp2[tid * moduleLen + i] = 0;
    }
    for(int i = 0; i < moduleLen; i++){
        for(int j = 0; j < moduleLen; j++){
            temp1[tid * moduleLen + j] = Ni[i * moduleLen + j];
        }
        
        bigIntegerMul(&(temp1[tid * moduleLen]),decomposeList[i * listLen + tid]%moduleChain[i],moduleLen);
        bigIntegerAdd(&temp2[tid * moduleLen],&temp1[tid * moduleLen],moduleLen);
    }
    cuDoubleComplex res;
    res.x = bigIntegerMod(&temp2[tid * moduleLen],bigN,moduleLen,tid);
    res.y = 0;
    composeList[tid] = res;
}
__global__ void cudacompose(unsigned long long *decomposeList,
                            unsigned long long* moduleChain,
                            int listLen,
                            int moduleLen,
                            unsigned long long* Ni,
                            unsigned long long *bigN,cuDoubleComplex * composeList,unsigned long long* temp1,unsigned long long* temp2,unsigned long long scale,int rescaleTimes){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;

    for(int i = 0; i < moduleLen; i++){
        temp2[tid * moduleLen + i] = 0;
    }
    for(int i = rescaleTimes; i < moduleLen; i++){
        for(int j = 0; j < moduleLen; j++){
            temp1[tid * moduleLen + j] = Ni[i * moduleLen + j];
        }
        // if(tid == 0){
        //     for(int j = 0; j < moduleLen; j++){
        //         printf("##%llu\n",temp1[tid * moduleLen + j]);
        //     }
        // }
        bigIntegerMul(&(temp1[tid * moduleLen]),decomposeList[i * listLen + tid]%moduleChain[i],moduleLen);
        bigIntegerAdd(&temp2[tid * moduleLen],&temp1[tid * moduleLen],moduleLen);
    }
    cuDoubleComplex res;
    res.x = bigIntegerMod(&temp2[tid * moduleLen],bigN,moduleLen,tid);
    res.y = 0;
    composeList[tid] = res;
    
}


class RNS{
    public:
    int N;
    unsigned long long* moduleChain;
    unsigned long long* Ni;
    unsigned long long *bigN;
    unsigned long long *buff1;
    unsigned long long *buff2;
    unsigned long long* moduleChain_h;
    unsigned long long scale;
    unsigned long long step = 2;
    RNS(int N,unsigned long long scale){
        this->N = N;
        this->scale = scale;
        buff1 = nullptr;
        buff2 = nullptr;
        moduleChain_h = (unsigned long long*)malloc(N * sizeof(unsigned long long));
        unsigned long long temp[8] = {1179649, 1376257, 1769473, 2424833, 2752513, 3604481, 3735553, 5308417};
        double half_h = 1.0;
        for(int i = 0; i < N; i++){
            moduleChain_h[i] = temp[i];
            half_h *= temp[i];
        }
        half_h /= 2.0;
        // cudaMemcpyToSymbol(&half_d, &half_h, sizeof(double));
        // genPrime(moduleChain_h,scale,N);
        unsigned long long** Ni_h = (unsigned long long**)calloc(N,sizeof(unsigned long long**));
        for(int i = 0; i < N;i++){
            Ni_h[i] = (unsigned long long*)calloc(N,sizeof(unsigned long long));
        }
        // unsigned long long bigN[N];
        unsigned long long *bigN_h = (unsigned long long*)calloc(N,sizeof(unsigned long long));
        unsigned long long *ti_h = (unsigned long long*)calloc(N,sizeof(unsigned long long));
        bigN_h[0] = 1;
        for(int i = 0;i < N; i++){
            Ni_h[i][0] = 1;
            ti_h[i] = 1;
        }
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                if(i==j)continue;
                bigIntegerMul(Ni_h[i],moduleChain_h[j],N);
                ti_h[i] = ti_h[i] * modpow128(moduleChain_h[j],moduleChain_h[i]-2,moduleChain_h[i]) % moduleChain_h[i];
            }
        }
        
        for(int i = 0; i < N; i++){
            bigIntegerMul(bigN_h,moduleChain_h[i],N);
        }
        for(int i = 0; i < N; i++){
            bigIntegerMul(Ni_h[i],ti_h[i],N);
        }
        free(ti_h);
        
        cudaMalloc(&moduleChain,N * sizeof(unsigned long long));
        cudaMalloc(&Ni,N * N * sizeof(unsigned long long));
        cudaMalloc(&bigN,N * sizeof(unsigned long long));

        cudaMemcpy(moduleChain,moduleChain_h,N * sizeof(unsigned long long),cudaMemcpyHostToDevice);
        cudaMemcpy(bigN,bigN_h,N * sizeof(unsigned long long),cudaMemcpyHostToDevice);

        for(int i = 0; i < N; i++){
            cudaMemcpy(Ni+(N*i),Ni_h[i],N * sizeof(unsigned long long),cudaMemcpyHostToDevice);
        }
    }
    unsigned long long* decompose(cuDoubleComplex *list,int listLen){
        unsigned long long * decomposeList;
        cudaMalloc(&decomposeList, listLen * N * sizeof(unsigned long long));
        cudadecompose<<<listLen/1024,1024>>>(list,moduleChain,listLen,N,decomposeList,scale);
        return decomposeList;
    }

    unsigned long long* decompose(unsigned long long *list,int listLen){
        unsigned long long * decomposeList;
        cudaMalloc(&decomposeList, listLen * N * sizeof(unsigned long long));
        cudadecompose<<<listLen/1024,1024>>>(list,moduleChain,listLen,N,decomposeList);
        return decomposeList;
    }
    cuDoubleComplex* compose(unsigned long long * decomposeList, int listLen){
        cuDoubleComplex * composeList;
        cudaMalloc(&composeList, listLen * N * sizeof(cuDoubleComplex));

        if(!buff1)cudaMalloc(&buff1, listLen * N * sizeof(unsigned long long));
        if(!buff2)cudaMalloc(&buff2, listLen * N * sizeof(unsigned long long));

        cudacompose<<<listLen/1024,1024>>>(decomposeList,moduleChain,listLen,N,Ni,bigN,composeList,buff1,buff2,scale);
        return composeList;        
    }
    cuDoubleComplex* compose(unsigned long long * decomposeList, int listLen,int rescaleTimes){
        cuDoubleComplex * composeList;
        cudaMalloc(&composeList, listLen * N * sizeof(cuDoubleComplex));

        if(!buff1)cudaMalloc(&buff1, listLen * N * sizeof(unsigned long long));
        if(!buff2)cudaMalloc(&buff2, listLen * N * sizeof(unsigned long long));


        unsigned long long** Ni_h = (unsigned long long**)calloc(N,sizeof(unsigned long long**));
        for(int i = 0; i < N;i++){
            Ni_h[i] = (unsigned long long*)calloc(N,sizeof(unsigned long long));
        }
        // unsigned long long bigN[N];
        unsigned long long *bigN_h = (unsigned long long*)calloc(N,sizeof(unsigned long long));
        unsigned long long *ti_h = (unsigned long long*)calloc(N,sizeof(unsigned long long));
        bigN_h[0] = 1;
        for(int i = 0;i < N; i++){
            Ni_h[i][0] = 1;
            ti_h[i] = 1;
        }
        for(int i = rescaleTimes; i < N; i++){
            for(int j = rescaleTimes; j < N; j++){
                if(i==j)continue;
                bigIntegerMul(Ni_h[i],moduleChain_h[j],N);
                ti_h[i] = ti_h[i] * modpow128(moduleChain_h[j],moduleChain_h[i]-2,moduleChain_h[i]) % moduleChain_h[i];
            }
        }
        
        for(int i = rescaleTimes; i < N; i++){
            bigIntegerMul(bigN_h,moduleChain_h[i],N);
        }
        for(int i = rescaleTimes; i < N; i++){
            bigIntegerMul(Ni_h[i],ti_h[i],N);
        }
        free(ti_h);
        unsigned long long *moduleChain_r,*Ni_r,*bigN_r;
        cudaMalloc(&moduleChain_r,N * sizeof(unsigned long long));
        cudaMalloc(&Ni_r,N * N * sizeof(unsigned long long));
        cudaMalloc(&bigN_r,N * sizeof(unsigned long long));

        cudaMemcpy(moduleChain_r,moduleChain_h,N * sizeof(unsigned long long),cudaMemcpyHostToDevice);
        cudaMemcpy(bigN_r,bigN_h,N * sizeof(unsigned long long),cudaMemcpyHostToDevice);

        for(int i = 1; i < N; i++){
            cudaMemcpy(Ni_r+(N*i),Ni_h[i],N * sizeof(unsigned long long),cudaMemcpyHostToDevice);
        }
        cudacompose<<<listLen/1024,1024>>>(decomposeList,moduleChain,listLen,N,Ni_r,bigN_r,composeList,buff1,buff2,scale,rescaleTimes);
        return composeList;        
    }
    // unsigned long long* compose_u(unsigned long long * decomposeList, int listLen){
    //     unsigned long long* composeList;
    //     cudaMalloc(&composeList, listLen * N * sizeof(unsigned long long));

    //     if(!buff1)cudaMalloc(&buff1, listLen * N * sizeof(unsigned long long));
    //     if(!buff2)cudaMalloc(&buff2, listLen * N * sizeof(unsigned long long));

    //     cudacompose<<<listLen/1024,1024>>>(decomposeList,moduleChain,listLen,N,Ni,bigN,composeList,buff1,buff2);
    //     return composeList;        
    // }
    void getParams(unsigned long long *q, unsigned long long* psi,unsigned long long* psiinv, unsigned long long* q_bit,int polylen){
        //hard
        unsigned long long psi_t[8] = {1034474, 1172569, 1557013, 349058, 785782, 1977521, 627147, 4561077};
        unsigned long long psiinv_t[8] = {441827, 1160428, 233173, 928390, 2113364, 499635, 1712567, 3973130};
        unsigned long long q_bit_t[8] = {21, 21, 21, 22, 22, 22, 22, 23};
        if(polylen == 4096){

            for(int i = 0; i < 8; i++){
                q[i] = moduleChain_h[i];
                psi[i] = psi_t[i];
                psiinv[i] = psiinv_t[i];
                q_bit[i] = q_bit_t[i];
            }
        }else{
            throw "wrong polylen";
        }
    }
    private:
    void genPrime(unsigned long long* moduleChain_h,unsigned long long scale,int N){
        if(scale % 2 == 1){
            scale+=1;
        }
        int cnt = 0;
        while(cnt < N){
            if(MillerRabin(scale+1)){
                moduleChain_h[cnt] = scale + 1;
                scale += step;
                cnt++;
            }
            else{
                scale += step;
            }
        }
    }

    bool MillerRabin(unsigned long long n){

        if(n == 2){
            return true;
        }
        if(n % 2 == 0){
            return false;
        }
        bool res = true;
        for(unsigned long long a = 2; a < 64 && a < n; a++){
            unsigned long long d = n - 1;

            while(d % 2 == 0){
                if(modpow128(a,d,n) == 1){
                    // continue;
                }
                else if(modpow128(a,d,n) == n-1){
                    break;
                }
                else{
                    return false;
                }
                d /=2;
            }
        }
        return true;
    }
};

// __global__ void init(unsigned long long* a){
//     for(int i = 0; i < 4096; i++){
//         a[i] = i ;
//     }
// }

// int main(){
//     RNS rns(8,10000);
//     unsigned long long *ptr;
//     cudaMalloc(&ptr,4096 * sizeof(unsigned long long));
//     init<<<1,1>>>(ptr);
//     unsigned long long *res = rns.decompose(ptr,4096);
//     printf("decompose finish\n");
//     unsigned long long *ori = rns.compose(res,4096);
//     print<<<1,1>>>(ori);
//     cudaDeviceSynchronize();
// }
void bigMulTest(){
    srand((unsigned)time(NULL)); 
    for(int i = 0; i < 1000000; i++){
        unsigned long long r = rand();
        r = (r * 11451419 + 3777);
        unsigned long long s = rand();
        s = (s * 11451419 + 3777);
        uint128 res;
        bigMul(r,s,res);
        unsigned __int128 rr = r;
        unsigned __int128 ss = s;
        if(res.high != (unsigned long long)((rr*ss)>>64)){
            printf("%llu,%llu,%llu,%llu\n",r,s,res.low , (unsigned long long)((rr*ss)>>64));
        }
    }
}
