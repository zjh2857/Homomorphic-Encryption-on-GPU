
__global__ void genRandom(unsigned long long *randomVec,double scale){
    if(scale != 114514){
        scale = 0;
    }
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    randomVec[tid] = (tid * 3 + 5) % 7 *scale;

}