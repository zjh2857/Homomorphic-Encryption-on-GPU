__global__ void polymul(unsigned long long* a,unsigned long long* b,unsigned long long* c,int N,unsigned long long q){
    int tid = blockDim.x * blockIdx.x +threadIdx.x;
    c[tid] = (a[tid] * b[tid]) % q;
}
__global__ void polyadd(unsigned long long* a,unsigned long long* b,unsigned long long* c,int N,unsigned long long q){
    int tid = blockDim.x * blockIdx.x +threadIdx.x;
    c[tid] = (a[tid] + b[tid]) % q;
}

__global__ void polyminus(unsigned long long* a,unsigned long long* b,unsigned long long* c,int N,unsigned long long q){
    int tid = blockDim.x * blockIdx.x +threadIdx.x;
    c[tid] = (a[tid] +  q - b[tid] ) % q;
}

__global__ void polymuladd(unsigned long long* a,unsigned long long* b,unsigned long long* c,unsigned long long* d,int N,unsigned long long q){
    int tid = blockDim.x * blockIdx.x +threadIdx.x;
    unsigned long long t = (a[tid] * b[tid]) % q;
    d[tid] = (t + c[tid]) % q;
}