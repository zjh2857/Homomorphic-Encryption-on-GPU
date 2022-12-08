cmake_minimum_required(VERSION 3.21)
project(ckks CUDA )

set(CMAKE_CUDA_STANDARD 14)

# set(CMAKE_CUDA_FLAGS -lcufft)
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lcufft)
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lcufft)

add_executable(ckks ckkstest.cu encoder.cuh helper.cuh parameter.cuh uint128.cuh)
target_link_libraries(ckks -lcufft)
set_target_properties(ckks PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)