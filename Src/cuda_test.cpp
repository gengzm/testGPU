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
		std::cout << "ʹ��GPU device " << i << ": " << devProp.name << std::endl;
		std::cout << "�豸ȫ���ڴ������� " << devProp.totalGlobalMem / 1024 / 1024 << "MB" << std::endl;
		std::cout << "SM��������" << devProp.multiProcessorCount << std::endl;
		std::cout << "ÿ���߳̿�Ĺ����ڴ��С��" << devProp.sharedMemPerBlock / 1024.0 << " KB" << std::endl;
		std::cout << "ÿ���߳̿������߳�����" << devProp.maxThreadsPerBlock << std::endl;
		std::cout << "�豸��һ���߳̿飨Block���ֿ��õ�32λ�Ĵ��������� " << devProp.regsPerBlock << std::endl;
		std::cout << "ÿ��EM������߳�����" << devProp.maxThreadsPerMultiProcessor << std::endl;
		std::cout << "ÿ��EM������߳�������" << devProp.maxThreadsPerMultiProcessor / 32 << std::endl;
		std::cout << "�豸�϶ദ������������ " << devProp.multiProcessorCount << std::endl;
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

	// CPU�˷����ڴ�
	A = (float*)malloc(size);
	B = (float*)malloc(size);
	C = (float*)malloc(size);

	// ��ʼ������
	for (int i = 0; i < n; i++)
	{
		A[i] = 90.0;
		B[i] = 10.0;
	}

	// GPU�˷����ڴ�
	cudaMalloc((void**)&Ad, size);
	cudaMalloc((void**)&Bd, size);
	cudaMalloc((void**)&Cd, size);

	// CPU�����ݿ�����GPU��
	cudaMemcpy(Ad, A, size, cudaMemcpyHostToDevice);
	cudaMemcpy(Bd, B, size, cudaMemcpyHostToDevice);
	cudaMemcpy(Bd, B, size, cudaMemcpyHostToDevice);

	// ����kernelִ�����ã���1024*1024/512����block��ÿ��block������512���߳�
	dim3 dimBlock(512);
	dim3 dimGrid(n / 512);

	Plus <<< dimGrid, dimBlock >>> (Ad, Bd, Cd, n);

	// ����GPU�˼���õĽ��������CPU��
	cudaMemcpy(C, Cd, size, cudaMemcpyDeviceToHost);

	// У�����
	float max_error = 0.0;
	for (int i = 0; i < n; i++)
	{
		max_error += fabs(100.0 - C[i]);
	}

	cout << "max error is " << max_error << endl;

	// �ͷ�CPU�ˡ�GPU�˵��ڴ�
	free(A);
	free(B);
	free(C);
	cudaFree(Ad);
	cudaFree(Bd);
	cudaFree(Cd);

	return 0;
}