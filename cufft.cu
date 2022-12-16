#include"iostream"
#include"cuda_runtime_api.h"
#include"device_launch_parameters.h"
#include"cufft.h"
using namespace std;
//FFT反變換後，用於規範化的函數
__global__ void normalizing(cufftDoubleComplex* data,int data_len)
{
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	data[idx].x /= data_len;
	data[idx].y /= data_len;
}
void Check(cudaError_t status)
{
	if (status != cudaSuccess)
	{
		cout << "行號:" << __LINE__ << endl;
		cout << "錯誤:" << cudaGetErrorString(status) << endl;
	}
}
int main()
{
	const int Nt =4;
	const int BATCH = 1;
	//BATCH用於批量處理一批一維數據，當BATCH=2時
	//則將0-1024，1024-2048作爲兩個一維信號做FFT處理變換
	cufftDoubleComplex* host_in, *host_out, *device_in, *device_out;
	//主機內存申請及初始化--主機鎖頁內存
	Check(cudaMallocHost((void**)&host_in, Nt * sizeof(cufftDoubleComplex)));
	Check(cudaMallocHost((void**)&host_out, Nt * sizeof(cufftDoubleComplex)));
	// for (int i = 1; i < Nt; i+=2)
	// {
	// 	host_in[i].x = 1;
	// 	host_in[i].y = 100;
	// }

	// for (int i = 0; i < Nt; i+=1)
	// {
	// 	host_in[i].x = 1;
	// 	// host_in[i].y = -100;
	// }
	host_in[0].x = 1;host_in[0].y = 1;
	host_in[1].x = 3;host_in[1].y = -4;
	host_in[2].x = 3;host_in[2].y=4;
	host_in[3].x = 1;host_in[3].y = -1;

	//設備內存申請
	Check(cudaMalloc((void**)&device_in, Nt * sizeof(cufftDoubleComplex)));
	Check(cudaMalloc((void**)&device_out, Nt * sizeof(cufftDoubleComplex)));
	//數據傳輸--H2D
	Check(cudaMemcpy(device_in, host_in, Nt * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice));

	//創建cufft句柄
	cufftHandle cufftForwrdHandle, cufftInverseHandle;
	cufftPlan1d(&cufftForwrdHandle, Nt, CUFFT_Z2Z, BATCH);
	cufftPlan1d(&cufftInverseHandle, Nt, CUFFT_Z2Z, BATCH);

	//執行fft正變換
	cufftExecZ2Z(cufftForwrdHandle, device_in, device_out, CUFFT_FORWARD);

	//數據傳輸--D2H
	Check(cudaMemcpy(host_out, device_out, Nt * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost));

	//設置輸出精度--正變換結果輸出
	cout << "逆變換結果:" << endl;
	// cout.setf(20);
	for (int i = 0; i < Nt; i++)
	{
		cout << host_out[i].x<< "+j*" << host_out[i].y << endl;
	}

	//執行fft反變換
	// cufftExecZ2Z(cufftInverseHandle,  device_out, device_in, CUFFT_INVERSE);
	
	// //IFFT結果是真值的N倍，因此要做/N處理
	// dim3 grid(Nt/128); 
	// dim3 block(128);
	// normalizing << <grid, block >> > (device_in,Nt);

	// //數據傳輸--D2H
	// Check(cudaMemcpy(host_in, device_in, Nt * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost));

	// //設置輸出精度--反變換結果輸出
	// cout << "反變換結果:" << endl;
	// // cout.setf(20);
	// for (int i = 0; i < Nt; i++)
	// {
	// 	cout << host_in[i].x << "+j*" << host_in[i].y << endl;
	// }
	// // cin.get();
	return 0;
}