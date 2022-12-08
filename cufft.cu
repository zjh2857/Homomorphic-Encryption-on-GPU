#include"iostream"
#include"cuda_runtime_api.h"
#include"device_launch_parameters.h"
#include"cufft.h"
using namespace std;
#define Check(call)														\
{																		\
	cudaError_t status = call;											\
	if (status != cudaSuccess)											\
	{																	\
		cout << "行號:" << __LINE__ << endl;							\
		cout << "錯誤:" << cudaGetErrorString(status) << endl;			\
	}																	\
}

//FFT反變換後，用於規範化的函數
__global__ void normalizing(cufftDoubleReal* data, int data_len)
{
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	if (idx<data_len)
	{
		data[idx] /=(data_len);
	}
}

int main()
{
	const int Nt = 4;
	const int BATCH = 1;
	//BATCH用於批量處理一批一維數據，當BATCH=2時
	//則將0-1024，1024-2048作爲兩個一維信號做FFT處理變換
	cufftDoubleReal* host_in,  *device_in;
	cufftDoubleComplex* host_out, *device_out;
	//主機內存申請及初始化--主機鎖頁內存
	Check(cudaMallocHost((void**)&host_in, Nt * sizeof(cufftDoubleReal)));
	//特別要注意：這裏的輸出長度變爲(Nt/2+1)
	Check(cudaMallocHost((void**)&host_out, (Nt / 2 + 1) * sizeof(cufftDoubleComplex)));
	for (int i = 0; i < Nt; i++)
	{
		host_in[i] = 1;
	}
	host_in[1] = -1;
	//設備內存申請
	Check(cudaMalloc((void**)&device_in, Nt * sizeof(cufftDoubleReal)));
	//特別要注意：這裏的輸出長度變爲(Nt/2+1)
	Check(cudaMalloc((void**)&device_out, (Nt / 2 + 1) * sizeof(cufftDoubleComplex)));
	
	//數據傳輸--H2D
	Check(cudaMemcpy(device_in, host_in, Nt * sizeof(cufftDoubleReal), cudaMemcpyHostToDevice));

	//創建cufft句柄
	cufftHandle cufftForwrdHandle, cufftInverseHandle;
	cufftPlan1d(&cufftForwrdHandle, Nt, CUFFT_D2Z, BATCH);
	cufftPlan1d(&cufftInverseHandle, Nt, CUFFT_Z2D, BATCH);

	//執行fft正變換
	cufftExecD2Z(cufftForwrdHandle, device_in, device_out);//由於D2Z的方向是固定的，無需填入參數

	//數據傳輸--D2H
	Check(cudaMemcpy(host_out, device_out, (Nt/2+1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost));

	//設置輸出精度--正變換結果輸出
	cout << "正變換結果:" << endl;
	// cout.setf(20);
	for (int i = 0; i < (Nt / 2 + 1); i++)
	{
		cout << host_out[i].x << "+j*" << host_out[i].y << endl;
	}

	//執行fft反變換
	cufftExecZ2D(cufftInverseHandle, device_out, device_in);//由於Z2D的方向是固定的，無需填入參數
	
	//IFFT結果是真值的N倍，因此要做/N處理
	dim3 grid(ceil((Nt / 2 + 1) / 128.0) + 1);
	dim3 block(128);
	normalizing << <grid, block >> > (device_in, Nt);

	//數據傳輸--D2H
	Check(cudaMemcpy(host_in, device_in, Nt * sizeof(cufftDoubleReal), cudaMemcpyDeviceToHost));

	//設置輸出精度--反變換結果輸出
	cout << "反變換結果:" << endl;
	// cout.setf(20);
	for (int i = 0; i < Nt; i++)
	{
		cout << host_in[i] << endl;
	}
	// cin.get();
	return 0;
}