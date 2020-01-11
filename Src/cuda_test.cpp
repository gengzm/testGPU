#include <cuda_runtime.h>
#include <stdio.h>
#include <device_launch_parameters.h>
#include <stdlib.h>
#include <iostream>



int dev_msg()
{
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	std::cout << "device count " << deviceCount << std::endl;
	for (int i = 0; i < deviceCount; i++)
	{
		cudaDeviceProp devProp;
		cudaGetDeviceProperties(&devProp, i);
		std::cout << "使用GPU device " << i << ": " << devProp.name << std::endl;
		std::cout << "设备全局内存总量： " << devProp.totalGlobalMem / 1024 / 1024 << "MB" << std::endl;
		std::cout << "SM的数量：" << devProp.multiProcessorCount << std::endl;
		std::cout << "每个线程块的共享内存大小：" << devProp.sharedMemPerBlock / 1024.0 << " KB" << std::endl;
		std::cout << "每个线程块的最大线程数：" << devProp.maxThreadsPerBlock << std::endl;
		std::cout << "设备上一个线程块（Block）种可用的32位寄存器数量： " << devProp.regsPerBlock << std::endl;
		std::cout << "每个EM的最大线程数：" << devProp.maxThreadsPerMultiProcessor << std::endl;
		std::cout << "每个EM的最大线程束数：" << devProp.maxThreadsPerMultiProcessor / 32 << std::endl;
		std::cout << "设备上多处理器的数量： " << devProp.multiProcessorCount << std::endl;
		std::cout << "======================================================" << std::endl;

	}
	return 0;
}

using namespace std;

__global__ void Plus(float A[], float B[], float C[], int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	C[i] = A[i] + B[i];
}

int main()
{

	float*A, *Ad, *B, *Bd, *C, *Cd;
	int n = 1024 * 1024;
	int size = n * sizeof(float);

	// CPU端分配内存
	A = (float*)malloc(size);
	B = (float*)malloc(size);
	C = (float*)malloc(size);

	// 初始化数组
	for (int i = 0; i < n; i++)
	{
		A[i] = 90.0;
		B[i] = 10.0;
	}

	// GPU端分配内存
	cudaMalloc((void**)&Ad, size);
	cudaMalloc((void**)&Bd, size);
	cudaMalloc((void**)&Cd, size);

	// CPU的数据拷贝到GPU端
	cudaMemcpy(Ad, A, size, cudaMemcpyHostToDevice);
	cudaMemcpy(Bd, B, size, cudaMemcpyHostToDevice);
	cudaMemcpy(Bd, B, size, cudaMemcpyHostToDevice);

	// 定义kernel执行配置，（1024*1024/512）个block，每个block里面有512个线程
	dim3 dimBlock(512);
	dim3 dimGrid(n / 512);

	Plus <<< dimGrid, dimBlock >>> (Ad, Bd, Cd, n);

	// 将在GPU端计算好的结果拷贝回CPU端
	cudaMemcpy(C, Cd, size, cudaMemcpyDeviceToHost);

	// 校验误差
	float max_error = 0.0;
	for (int i = 0; i < n; i++)
	{
		max_error += fabs(100.0 - C[i]);
	}

	cout << "max error is " << max_error << endl;

	// 释放CPU端、GPU端的内存
	free(A);
	free(B);
	free(C);
	cudaFree(Ad);
	cudaFree(Bd);
	cudaFree(Cd);

	return 0;
}