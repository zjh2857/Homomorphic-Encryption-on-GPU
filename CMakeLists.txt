cmake_minimum_required(VERSION 3.1)
project(ckks CUDA )

set(CMAKE_CUDA_STANDARD 14)

# set(CMAKE_CUDA_FLAGS -lcufft)
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lcufft)
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -lcufft)

add_executable(ckks ckkstest.cu encoder.cuh helper.cuh parameter.cuh uint128.cuh encryptor.cuh polycalc.cuh random.cuh evaluator.cuh)
target_link_libraries(ckks -lcublas)
set_target_properties(ckks PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)