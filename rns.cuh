#include <cstdio>
#include <stdlib.h>

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
    // printf("%llu,%llu\n",temp.low,temp.high);
    unsigned long long carry = temp.high;
    a[0] = temp.low;
    // printf("%llu,%llu\n",a[0],temp.low);
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
// __device__ __host__ __forceinline__ void bigIntegerMul_d(unsigned long long * a, unsigned long long b,int size){
//     uint128 temp;
//     bigMul(a[0],b,temp);
    
//     // printf("%llu,%llu,%llu,%llu\n",a[0],b,temp.low,temp.high);
//     unsigned long long carry = temp.high;
//     // printf("%llu,%llu\t",a[0],temp.low);
//     a[0] = temp.low;
//     // printf("%llu,%llu\n",a[0],temp.low);
//     for(int i = 1; i < size; i++){
//         if(!a[i] && !carry){
//             break;
//         }
//         uint128 temp;
//         bigMul(a[i],b,temp);
//         // temp = te
//         a[i] = temp.low;
//         // a[i] = (temp & 0xffffffffffffffff);
//         if(a[i] + carry < a[i]){
//             a[i] += carry;
//             carry = 1 + temp.high;
//         }
//         else{
//             a[i] += carry;
//             carry = temp.high;
//         }
//     }
// }
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
    // {for(int i = 0; i < size; i++){
    //     printf("%llu\t",temp[i]);
    // }printf("\n");}
    // {for(int i = 0; i < size; i++){
    //     printf("%llu\t",a[i]);
    // }printf("\n");}
    int res = 0;
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
    // printf("^^%d\n",borrow);
    // for(int i = 0; i < size; i++){
    //     printf("%llu\t",temp[i]);
    // }printf("\n\n\n\n");
    if(borrow == 1){
        return -1;
    }
    return 1;
    // for(int i = size-1; i > 0; i--){
    //     if(temp[i] !=0 ){
    //         return 1;
    //     }
    // }
    // return 0;
    // for(int i = 0; i < size; i++){
    //     if(a[i] < borrow){
    //         return -1;
    //     }
    //     if(a[i] - borrow > temp1[i]){
    //         borrow = 0;
    //     }else{
    //         borrow = 1;
    //     }
    //     // if(a[i] > temp[i]){
    //     //     // free(temp);
    //     //     return 1;
    //     // }
    //     // else if(a[i] < temp[i]){
    //     //     // free(temp);
    //     //     return -1;
    //     // }
    // }
    // free(temp);
    return 0;
}

__device__ __host__ __forceinline__ int isneg_d(unsigned long long * a, unsigned long long *b ,unsigned long long p,int size){
    // unsigned long long* temp = (unsigned long long*)malloc(size * sizeof(unsigned long long)); 
    unsigned long long temp[64];
    for(int i = 0; i < size; i++){
        temp[i] = b[i];
    }

    bigIntegerMul(temp,p,size);
    {for(int i = 0; i < size; i++){
        printf("%llu\t",temp[i]);
    }printf("\n");}
    {for(int i = 0; i < size; i++){
        printf("%llu\t",a[i]);
    }printf("\n");}
    int res = 0;
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
    printf("^^%d\n",borrow);
    for(int i = 0; i < size; i++){
        printf("%llu\t",temp[i]);
    }printf("\n\n\n\n");
    if(borrow == 1){
        return -1;
    }
    return 1;
    // for(int i = size-1; i > 0; i--){
    //     if(temp[i] !=0 ){
    //         return 1;
    //     }
    // }
    // return 0;
}
__device__ __host__ __forceinline__ unsigned long long bigInteger2u64(unsigned long long *a,unsigned long long *b,unsigned long long p,int size,int tid){
    // unsigned long long* temp = (unsigned long long*)malloc(size * sizeof(unsigned long long)); 
    unsigned long long temp[64];
    for(int i = 0; i < size; i++){
        temp[i] = b[i];
    }

    bigIntegerMul(temp,p,size);
    // {for(int i = 0; i < size; i++){
    //     printf("%llu\t",temp[i]);
    // }printf("\n");}
    // {for(int i = 0; i < size; i++){
    //     printf("%llu\t",a[i]);
    // }printf("\n");}
    // int res = 0;
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
    // printf("^^%d\n",borrow);
    // for(int i = 0; i < size; i++){
    //     printf("%llu\t",temp[i]);
    // }printf("\n\n\n\n");
    double res = 0;
    double twopow64 = 18446744073709551616.0;
    double base = 1.0;
    if(tid==0){
        for(int i = 0; i < size; i++){
            printf("@%llu\t",temp[i]);
        }printf("\n");
    }
    if(borrow == 0){
        for(int i = 0; i < size; i++){
            res += temp[i] * base;
            base*= 18446744073709551616.0;
        }
    }
    else{

        for(int i = 0; i < size; i++){
            temp[i] = ~temp[i];
        }
        unsigned long long carry = 1;
        for(int i = 0; i < size; i++){
            temp[i] += carry;
            if(temp[i] == 0){
                carry = 1;
            }else{
                carry = 0;
            }
        }
    if(tid==0){
        for(int i = 0; i < size; i++){
            printf("@%llu\t",temp[i]);
        }printf("\n");
    }        
        for(int i = 0; i < size; i++){
            res += temp[i] * base;
            base*= 18446744073709551616.0;
        }
        res = -res;
    }
    // printf("!!!$$\n");
    if(tid == 0){printf("^^^%lf\n",res);
    if(tid == 0)printf("fff%lf\n",res - (unsigned long long)(res / 288230376135196673) * 288230376135196673.0);
    unsigned long long r = res - (unsigned long long)(res / 288230376135196673) * 288230376135196673.0;
    printf("uuu%llu\n",r);}
    return res - (unsigned long long)(res / 288230376135196673) * 288230376135196673.0;
    // return 0;
}
__device__ __host__ __forceinline__ unsigned long long bigIntegerMod(unsigned long long * a, unsigned long long *b,int size,int tid){
    unsigned long long l = 0;
    unsigned long long r = (1llu << 63);
    unsigned long long guess;
    // printf("%d((((\n",tid);
    // if(a[0] != 0){
    //     for(int i = 0; i < 8; i++){
    //         printf("%llu\t",a[i]);
    //     }printf("\n");
    //     for(int i = 0; i < 8; i++){
    //         printf("%llu\t",b[i]);
    //     }printf("\n");
    // }else{
    //     return 0;
    // }
    // return 0;
    if(tid == 0){
        for(int i = 0; i < 8; i++){
            printf("%llu\t",a[i]);
        }printf("\n");
        for(int i = 0; i < 8; i++){
            printf("%llu\t",b[i]);
        }printf("\n");            
    }
    int cnt = 0;
    unsigned long long p = 0;
    while(l < r){

        if(cnt++ > 100){
            // printf("==============\n");
        //     for(int i = 0; i < 8; i++){
        //         printf("%llu\t",a[i]);
        //     }printf("\n");
        //     for(int i = 0; i < 8; i++){
        //         printf("%llu\t",b[i]);
        //     }printf("\n");            
        //     printf("==============\n");
        // printf("%llu,%llu\n",l,r);
            // printf("%d\n",tid);
            printf("114514\n");
            return 114514;
        }
        if(tid == 0)printf("%llu,%llu,%llu\n",l,r,(l+r)/2);
        unsigned long long guess = (l + r)/2;
        // if(guess == 0){
        //     return 0;
        // }
        int res = isneg(a,b,guess,size);
        if(res == 1){
            p = guess;
            l = guess + 1;

        }
        else if(res == -1){
            r = guess;
        }
        // else{
        //     return a[0] - (b[0] * guess);
        // }
    }
    // if(l != r){
    //     printf("@@@@@@@@@\n");
    //     return 114514;
    // }
    if(tid == 0)printf("%llu&&&&&&&&&&&&&&&\n",p);
    return bigInteger2u64(a,b,p,size,tid);

}


// __device__ __host__ __forceinline__ unsigned long long bigIntegerMod_d(unsigned long long * a, unsigned long long *b,int size){
//     unsigned long long l = 0;
//     unsigned long long r = (1llu << 63);
//     unsigned long long guess;
//     // if(a[0] != 0){
//     //     for(int i = 0; i < 8; i++){
//     //         printf("%llu\t",a[i]);
//     //     }printf("\n");
//     //     for(int i = 0; i < 8; i++){
//     //         printf("%llu\t",b[i]);
//     //     }printf("\n");
//     // }else{
//     //     return 0;
//     // }
//     // return 0;
//     int cnt = 0;
//     while(l <= r){
//         if(cnt++ > 100){
//             return 114514;
//         }
//         printf("%llu,%llu,%llu\n",l,r,(l+r)/2);

//         // printf("%llu,%llu\n",l,r);
//         unsigned long long guess = (l + r)/2;
//         int res = isneg_d(a,b,guess,size);
//         if(res == 1){
//             l = guess + 1;
//         }
//         else if(res == -1){
//             r = guess;
//         }
//         else{
//             return a[0] - (b[0] * guess);
//         }
//     }
//     return 0;
// }
__global__ void cudadecompose(unsigned long long *list,unsigned long long* moduleChain,int listLen,int moduleLen,unsigned long long * decomposeList){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    for(int i = 0; i < moduleLen; i++){
        decomposeList[i * listLen + tid] = list[tid] % moduleChain[i];
    }
    // if(tid == 1){
    //     printf("**%llu\n",list[1]);
    // }
}
__global__ void cudacompose(unsigned long long *decomposeList,
                            unsigned long long* moduleChain,
                            int listLen,
                            int moduleLen,
                            unsigned long long* Ni,
                            unsigned long long *bigN,unsigned long long * composeList,unsigned long long* temp1,unsigned long long* temp2){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    // unsigned long long temp1[64];
    // unsigned long long temp2[64];'
    // printf("!!!!%d\n",tid);
    for(int i = 0; i < moduleLen; i++){
        temp2[tid * moduleLen + i] = 0;
    }
    for(int i = 0; i < moduleLen; i++){
        for(int j = 0; j < moduleLen; j++){
            temp1[tid * moduleLen + j] = Ni[i * moduleLen + j];
        }
        
    // if(i==0)printf("%llu,%llu\n",decomposeList[i * listLen + tid],temp1[tid * moduleLen]);
        // if(tid == 1967){
        //     // printf("%llu,%llu,%llu\n",decomposeList[i * listLen + tid],moduleChain[i],decomposeList[i * listLen + tid]%moduleChain[i]);
        //     for(int i = 0; i < 8; i++){
        //         printf("%llu\t",temp1[tid * moduleLen + i]);
        //     }printf("***\n");

        // }
        bigIntegerMul(&(temp1[tid * moduleLen]),decomposeList[i * listLen + tid]%moduleChain[i],moduleLen);
        // if(i == 0)printf("%llu\n",temp1[tid * moduleLen]);
        // if(tid == 0){
        //     for(int i = 0; i < 8; i++){
        //         printf("%llu\t",temp1[tid * moduleLen + i]);
        //     }printf("***\n");
        // if(tid == 1967){
        //     // printf("%llu,%llu,%llu\n",decomposeList[i * listLen + tid],moduleChain[i],decomposeList[i * listLen + tid]%moduleChain[i]);
        //     for(int i = 0; i < 8; i++){
        //         printf("%llu\t",temp1[tid * moduleLen + i]);
        //     }printf("***\n");

        // }
        // }
        bigIntegerAdd(&temp2[tid * moduleLen],&temp1[tid * moduleLen],moduleLen);
    }
    // if(tid == 114){
    //     printf("%llu\n",decomposeList[i * listLen + tid]);
    //     // printf("%llu\t",Ni[0]);
    //     for(int i = 0; i < 8; i++){
    //         printf("%llu\t",temp1[i]);
    //     }printf("\n");
    //     // for(int i = 0; i < 8; i++){
    //     //     printf("%llu\t",bigN[i]);
    //     // }printf("\n");
    // }
    // }
    // if(tid == 514){
    //     // printf("%llu\t",Ni[0]);
    //     // printf("%p\n",temp1);
    //     for(int i = 0; i < 8; i++){
    //         printf("%llu\t",temp2[tid * moduleLen + i]);
    //     }printf("\n");
    //     for(int i = 0; i < 8; i++){
    //         printf("%llu\t",bigN[i]);
    //     }printf("\n");
    // }
    // if(tid == 1967){
    //     // for(int i = 0; i < 8; i++){
    //     //     printf("%llu\t",temp2[tid * moduleLen + i]);
    //     // }printf("!!!!!!\n");
    //     for(int i = 0; i < 8; i++){
    //         printf("%llu\t",bigN[i]);
    //     }printf("!!!!!!\n");
    // }
    unsigned long long res;
    // if(temp2[tid * moduleLen] == 0 && tid == 3488){
    //     printf("%d\n",tid);
    //     for(int i = 0; i < 8; i++){
    //         printf("(((%llu\t",decomposeList[i * listLen + tid]);
    //     }printf("\n");
    // }

    res = bigIntegerMod(&temp2[tid * moduleLen],bigN,moduleLen,tid);
    // if(res == 114514){
    //     printf("###%d\n",tid);
    //     // for(int i = 0; i < 8; i++){
    //     //     printf("%llu\t",temp2[tid * moduleLen + i]);
    //     // }printf("!!!!!!\n");
    //     // for(int i = 0; i < 8; i++){
    //     //     printf("%llu\t",bigN[tid * moduleLen + i]);
    //     // }printf("!!!!!!\n");
    // }
    // if(tid == 0){
    //     printf("@@@@%llu\n",res);
    // }
    // printf("%llu\n",res);
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
    unsigned long long step = 2;
    RNS(int N,unsigned long long scale){
        this->N = N;
        buff1 = nullptr;
        buff2 = nullptr;
        unsigned long long* moduleChain_h = (unsigned long long*)malloc(N * sizeof(unsigned long long));
        genPrime(moduleChain_h,scale,N);
        // moduleChain_h[N-1] = 288230376135196673;
        // for(int i = 0; i < N; i++){
        //     printf("%llu\t",moduleChain_h[i]);
        // }printf("\n");
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
                ti_h[i] = ti_h[i] * qpow(moduleChain_h[j],moduleChain_h[i]-2,moduleChain_h[i]) % moduleChain_h[i];
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
        // for(int i = 0; i < N; i++){
        //     for(int j = 0; j < N;j++){
        //         printf("%llu\t",Ni_h[i][j]);
        //     }
        //     printf("\n");
        // }
        // exit(1);
        for(int i = 0; i < N; i++){
            cudaMemcpy(Ni+(N*i),Ni_h[i],N * sizeof(unsigned long long),cudaMemcpyHostToDevice);
        }
        // print<<<1,1>>>(Ni);
        // exit(1);
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

        if(!buff1)cudaMalloc(&buff1, listLen * N * sizeof(unsigned long long));
        if(!buff2)cudaMalloc(&buff2, listLen * N * sizeof(unsigned long long));
        // print<<<1,1>>>(decomposeList);

        // cudaDeviceSynchronize();
        printf("##%d\n",listLen);
        cudacompose<<<listLen/1024,1024>>>(decomposeList,moduleChain,listLen,N,Ni,bigN,composeList,buff1,buff2);
        // cudaDeviceSynchronize();

cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess) {
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    // Possibly: exit(-1) if program cannot continue....
}
        // print<<<1,1>>>(composeList);
cudaDeviceSynchronize();
        return composeList;        
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
                // psi[cnt] = GetOrigen(scale+1);
                // psiinv[cnt] = pow(psi[cnt],scale-1,scale+1);
                scale += step;
                cnt++;
            }
            else{
                scale += step;
            }
        }
    }
    // void GetOrigen(unsigned long long prime){
    //     for(unsigned long long p = 2; p < 64){
            
    //     }
        
    // }
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
// int main(){
//     bigMulTest();
// }