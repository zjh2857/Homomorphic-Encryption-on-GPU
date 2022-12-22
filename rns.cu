#include <cstdio>
#include <stdlib.h>

using namespace std;
class uint128_t{
    public:

    unsigned long long low;
    unsigned long long high;
    uint128_t(unsigned long long low, unsigned long long high){
        this->low = low;
        this->high = high;
    }
};
__host__ __device__ __forceinline__ uint128_t bigMul(unsigned long long a,unsigned long long b){
    u_int32_t a0 = a;
    u_int32_t a1 = a >> 32;
    u_int32_t b0 = b;
    u_int32_t b1 = b >> 32;
    u_int64_t low = a0 * b0 + ( (a0 * b1 + a1 * b0) & 0xffffffff );
    u_int64_t high = a1 * b1 + ( (a0 * b1 + a1 * b0) >> 32llu );
    uint128_t res(low,high);
    return res;
}
__device__ __host__ void bigIntegerMul(unsigned long long * a, unsigned long long b,int size){
    uint128_t temp = bigMul(a[0],b);
    unsigned long long carry = temp.high;
    
    for(int i = 1; i < size; i++){
        if(!a[i] && !carry){
            break;
        }
        uint128_t temp = bigMul(a[i],b);
        // temp = te
        a[i] = temp.low;
        // a[i] = (temp & 0xffffffffffffffff);
        if(a[i] + carry < a[i]){
            a[i] += carry;
            carry = 1 + temp.high;
        }
        else{
            a[i] += carry;
            carry = 1 + temp.high;
        }
    }
}
__device__ __host__ void bigIntegerAdd(unsigned long long * a, unsigned long long *b,int size){
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
__device__ __host__  int isneg(unsigned long long * a, unsigned long long *b ,unsigned long long p,int size){
    // unsigned long long* temp = (unsigned long long*)malloc(size * sizeof(unsigned long long)); 
    unsigned long long temp[64];
    for(int i = 0; i < size; i++){
        temp[i] = b[i];
    }
    bigIntegerMul(temp,p,size);
    for(int i = size - 1; i > 0; i--){
        if(a[i] > temp[i]){
            // free(temp);
            return 1;
        }
        else if(a[i] < temp[i]){
            // free(temp);
            return -1;
        }
    }
    // free(temp);
    return 0;
}
__device__ __host__ unsigned long long bigIntegerMod(unsigned long long * a, unsigned long long *b,int size){
    unsigned long long l = 0;
    unsigned long long r = (1llu << 63);
    unsigned long long guess;
    while(l <= r){
        printf("%llu,%llu\n",l,r);
        unsigned long long guess = (l + r)/2;
        int res = isneg(a,b,guess,size);
        if(res == 1){
            l = guess + 1;
        }
        else if(res == -1){
            r = guess;
        }
        else{
            return a[0] - (b[0] * guess);
        }
    }
    return 0;
}
__global__ void cudadecompose(unsigned long long *list,unsigned long long* moduleChain,int listLen,int moduleLen,unsigned long long * decomposeList){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    for(int i = 0; i < moduleLen; i++){
        decomposeList[i * listLen + tid] = list[tid] % moduleChain[i];
    }
}
__global__ void cudacompose(unsigned long long *decomposeList,
                            unsigned long long* moduleChain,
                            int listLen,
                            int moduleLen,
                            unsigned long long* Ni,
                            unsigned long long *bigN,unsigned long long * composeList){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned long long temp1[64];
    unsigned long long temp2[64];
    for(int i = 0; i < moduleLen; i++){
        temp1[i] = 0;
        temp2[i] = 0;
    }
    temp1[0] = 1;
    for(int i = 0; i < moduleLen; i++){
        bigIntegerMul(temp1,decomposeList[i * listLen + tid],moduleLen);
        bigIntegerAdd(temp2,temp1,moduleLen);
    }
    unsigned long long res = bigIntegerMod(temp2,bigN,moduleLen);
    composeList[tid] = res;
}


class RNS{
    public:
    int N;
    unsigned long long* moduleChain;
    unsigned long long* Ni;
    unsigned long long *bigN;
    RNS(int N,unsigned long long scale){
        this->N = N;
        unsigned long long* moduleChain_h = (unsigned long long*)malloc(N * sizeof(unsigned long long));
        genPrime(moduleChain_h,scale,N);
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
                bigIntegerMul(Ni_h[i],moduleChain[j],N);
                ti_h[i] = ti_h[i] * qpow(moduleChain[j],moduleChain[i]-2,moduleChain[i]) % moduleChain[i];
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
            cudaMemcpy(Ni+(N*i),bigN_h,N * sizeof(unsigned long long),cudaMemcpyHostToDevice);
        }
    }
    unsigned long long* decompose(unsigned long long *list,int listLen){
        unsigned long long * decomposeList;
        cudaMalloc(&decomposeList, listLen * N * sizeof(unsigned long long));
        cudadecompose<<<listLen/1024,1024>>>(list,moduleChain,listLen,N,decomposeList);
        return decomposeList;
    }
    unsigned long long* compose(unsigned long long * decomposeList, int listLen){
        unsigned long long * composeList;
        cudaMalloc(&composeList, listLen * N * sizeof(unsigned long long));
        cudacompose<<<listLen/1024,1024>>>(decomposeList,moduleChain,listLen,N,Ni,bigN,composeList);
        return composeList;        
    }
    private:
    void genPrime(unsigned long long* moduleChain_h,unsigned long long scale,int N){
        if(scale % 2 == 0){
            scale+=1;
        }
        int cnt = 0;
        while(cnt < N){
            if(MillerRabin(scale)){
                moduleChain_h[cnt++] = scale;
                scale += 2;
            }
            else{
                scale += 2;
            }
        }
    }
    unsigned long long qpow(unsigned long long a,unsigned long long b,unsigned long long q){
        unsigned long long r=1;
        // unsigned long long base = a;
        while(b){
            if(b&1)r = (a * r)%q;
            a = a * a % q;
            b >>= 1;
        }
        return r;
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
            // printf("%llu\n",a);
            while(d % 2 == 0){
                // printf("%llu,%llu,%llu\n",a,d,qpow(a,d,n));
                if(qpow(a,d,n) == 1){
                    // continue;
                }
                else if(qpow(a,d,n) == n-1){
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

__global__ void init(unsigned long long* a){
    for(int i = 0; i < 1024; i++){
        a[i] = i * 10000;
    }
}
__global__ void print(unsigned long long * a){
    printf("1234567\n");
    for(int i = 0; i < 8; i++){
        printf("%llu\t");
    }printf("\n");
}
int main(){
    RNS rns(1024,10000);
    unsigned long long *ptr;
    cudaMalloc(&ptr,1024 * sizeof(unsigned long long));
    init<<<1,1>>>(ptr);
    unsigned long long *res = rns.decompose(ptr,1024);
    unsigned long long *ori = rns.compose(res,1024);
    printf("%p\n",ori);
    print<<<1,1>>>(ori);
    cudaDeviceSynchronize();
}