#include <iostream>
using namespace std;

const long long M = 998244353;
const long long N = 8;
long long l = 0;
long long Ninv;
long long revphi[N];
long long revphiinv[N];
long long bitrevs[N];
long long bitrevl[N];

void print(long long *a){
    for(int i = 0; i < N; i++){
        cout << a[i] << " ";
    }
    cout << endl;    
}

long long qpow(long long a,long long b,long long n){
    long long ans = 1;
    long long base = a;
    while(b){
        if(b&1){
            ans *= base;
            ans %= n;
        }
        base = base * base % n;
        b >>= 1;
    }
    return ans;
}

long long _bitrevs(long long x){
    long long res = 0;
    for(int i = 0; i < l - 1; i++){
        res = (res << 1) | (x & 1);
        x >>= 1;
    }
    return res;
}
long long _bitrevl(long long x){
    long long res = 0;

    for(int i = 0; i < l; i++){
        res = (res << 1) | (x & 1);
        x >>= 1;
    }
    return res;
}
void init(){
    long long temp = N;
    while(temp != 1){
        temp=temp>>1;
        l++;
    }
    long long g = qpow(3,(M-1)/N,M);
    long long gi = qpow(332748118,(M-1)/N,M);
    for(int i = 0; i < N; i++){
        bitrevs[i] = _bitrevs(i);
        bitrevl[i] = _bitrevl(i);
    }
    for(int i = 0; i < N/2; i++){
        revphi[i] = qpow(g,bitrevs[i],M);
        revphiinv[i] = qpow(gi,bitrevs[i],M);
    }
    Ninv = qpow(N,M-2,M);
    
}
void swapInv(long long* a){
    for(int i = 0; i < N; i++){
        if(i < bitrevl[i])swap(a[i],a[bitrevl[i]]);
    }
}
long long* NTT(long long* a){
    long long t = N;
    long long m = 1;
    while(t > 1){
        t /= 2;
        for(int i = 0; i < m; i++){
            int j1 = 2 * i * t;
            int j2 = j1 + t;
            for(int j = j1; j < j2; j++){
                long long u = a[j] % M;
                long long v = a[j + t] * revphi[i] % M;
                a[j] = (u + v + M) % M;
                a[j + t] = (u - v + M) % M;
            }
        }
        m <<= 1;
    }
    return a;
}

long long* INTT(long long* a){
    swapInv(a);
    long long t = N;
    long long m = 1;
    while(t > 1){
        t /= 2;
        for(int i = 0; i < m; i++){
            int j1 = 2 * i * t;
            int j2 = j1 + t;
            for(int j = j1; j < j2; j++){
                long long u = a[j] % M;
                long long v = a[j + t] * revphiinv[i] % M;
                a[j] = (u + v + M) % M;
                a[j + t] = (u - v + M) % M;
            }
        }

        m <<= 1;
    }
    for(int i = 0; i < N; i++){
        a[i] = a[i] * Ninv % M;
    }
    swapInv(a);
    return a;
}

__global__ void cuNTT(long long* a,long long* revphi,long long l){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for(int i = 0; i < l; i++){
        int t = (N/2) >> i;
        int address = idx / t * t + idx;
        long long u = a[address];
        long long v = a[address+t] * revphi[idx/t] % M;

        a[address] = (u + v + M) % M;
        a[address + t] = (u - v + M) % M;
        __syncthreads();
    }
}


__global__ void cuINTT(long long* a,long long* revphiinv,long long l){
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for(int i = 0; i < l; i++){
        int t = (N/2) >> i;
        int address = idx / t * t + idx;
        long long u = a[address];
        long long v = a[address+t] * revphiinv[idx/t] % M;

        a[address] = (u + v + M) % M;
        a[address + t] = (u - v + M) % M;
        __syncthreads();
    }
}

__global__ void cuSwapInv(long long* a,long long* bitrevl){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < bitrevl[idx]){
        long long t = a[idx];
        a[idx] = a[bitrevl[idx]];
        a[bitrevl[idx]] = t;
    }
}

int main(){
    long long* a = (long long*)malloc(N * sizeof(long long));
    for(int i = 0; i < N;i++){
        a[i] = (i * 141 + 54) % 17;
    }
    init();
    int nByte = N * sizeof(long long);
    long long *a_d,*revphi_d,*revphiinv_d,*bitrevl_d;
    cudaMalloc(&a_d,nByte);
    cudaMemcpy(a_d,a,nByte,cudaMemcpyHostToDevice);
    cudaMalloc(&revphi_d,nByte);
    cudaMemcpy(revphi_d,revphi,nByte,cudaMemcpyHostToDevice);
    cudaMalloc(&revphiinv_d,nByte);
    cudaMemcpy(revphiinv_d,revphiinv,nByte,cudaMemcpyHostToDevice);
    cudaMalloc(&bitrevl_d,nByte);
    cudaMemcpy(bitrevl_d,bitrevl,nByte,cudaMemcpyHostToDevice);
    print(a);
    NTT(a);
    print(a);
    INTT(a);
    print(a);
    
    dim3 block_half(N/2);
    dim3 grid_half(N/2/block_half.x);
    dim3 block(N);
    dim3 grid(N/block.x); 
    long long* res = (long long *)malloc(nByte);
    cuNTT<<<grid_half,block_half>>>(a_d,revphi_d,l);
    cudaMemcpy(res,a_d,nByte,cudaMemcpyDeviceToHost);
    print(res);
    cuSwapInv<<<grid,block>>>(a_d,bitrevl_d);
    cudaMemcpy(res,a_d,nByte,cudaMemcpyDeviceToHost);
    print(res);
    cuINTT<<<grid_half,block_half>>>(a_d,revphiinv_d,l);
    cudaMemcpy(res,a_d,nByte,cudaMemcpyDeviceToHost);
    print(res);
    cuSwapInv<<<grid,block>>>(a_d,bitrevl_d);
    cudaMemcpy(res,a_d,nByte,cudaMemcpyDeviceToHost);
    print(res);
}