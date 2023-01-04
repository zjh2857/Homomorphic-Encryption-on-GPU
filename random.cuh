
__global__ void genRandom(unsigned long long *randomVec,unsigned long long scale){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    // randomVec[tid] = (1) * scale;
    // randomVec[tid] = 0;
    // if(scale == 114514){
    //     if(tid != 114514){
    //         // randomVec[tid] = 1179648;
    //         randomVec[tid] = tid * 2;

    //     }else{
    //         randomVec[tid] = 0;
    //     }
    //     return ;
    // }
    randomVec[tid] = scale;
    // if(scale == 0){
    //     randomVec[tid] = 0;
    // }
    // else if(scale == 1){
    //     randomVec[tid] = 1;
    // }
    // else{
    //     randomVec[tid] = 1179649 * 2 ;
    // }
}

__global__ void genRandom_s (unsigned long long *randomVec,unsigned long long scale){
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    randomVec[tid] = (114 * tid + 514) % 3777;
}