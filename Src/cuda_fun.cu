#include "cuda_fun.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#define WIDTH 204
#define HEIGHT 153
#define CELL_W 16
#define CELL_H 16

/*******
********
********  something else
********
********/
extern "C" __global__ void Plus(float A[], float B[], float C[], int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	C[i] = A[i] + B[i];
}

extern "C" __global__ void PictureKernel(float* in, float* out, int n, int m) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (row < m && col < n) {
		out[row*n + col] = 2 * in[row*n + col];
	}
}

extern "C" __global__ void MatrixMulKernel(float *dM, float * dN, float * dP, int width)
{
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	if (row < width && col < width) {
		float value = 0;
		for (int k = 0; k < width; k++) {
			value += dM[row*width + k] * dN[k*width + col];
		}
		dP[row*width + col] = value;
	}
}

/**
***  
***
***
***/
extern "C" __global__ void InitImageData(char* image, int width, int height)
{
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;


	if ((row < height) && (col < width)) {
		double x = (col - width /2 ) / (double)width;
		double y = (row - height /2) / (double)height;
		double x_2 = x * x;
		double y_2 = y * y;

		//double U = 3.0*(1 - 2*x + x_2)*exp(- x_2 - (y_2 + 2*y + 1)) - 10.0*(x / 5.0 - x * x_2 - y * y_2 * y_2)*exp(-x_2 - y_2) - 1.0 / 3.0*exp(-(x_2 + 2*x + 1) - y_2);
		double U = 3.0*(1 - x)*(1 - x)*exp(-x * x - (y + 1)*(y + 1)) - 10.0*(x / 5.0 - x * x*x - y * y*y*y*y)*exp(-x * x - y * y) - 1.0 / 3.0*exp(-(x + 1)*(x + 1) - y * y);

		double t = 127.5*(1.0 + cos(200 * U));

		
		image[col + row * width] = (char)(int)(t + 0.5);
	}

}


extern "C" __global__ void Prefilter_1Dm(double *coefficient, int length, double *pole, double tolerance, double gamma)
{
	//int i, n, k;
	double Lambda;
	Lambda = 6.0 / (6.0 * gamma + 1.0);

	//int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	// Applying the gain to original image
	
	coefficient[col] *= Lambda;

	for (int k = 0; k < 1; k++)	// Testing
	{
		// Compute the first causal coefficient
		*(coefficient) = InitialCausalCoefficient(coefficient, length, pole[k], tolerance);

		// Causal prefilter
		for (int n = 1; n < length; n++)
			coefficient[n] += pole[k] * coefficient[n - 1];

		//Compute the first anticausal coefficient
		*(coefficient + length - 1) = InitialAnticausalCoefficient(coefficient, length, pole[k]);

		//Anticausal prefilter
		for (int n = length - 2; n >= 0; n--)
			coefficient[n] = pole[k] * (coefficient[n + 1] - coefficient[n]);
	}
	
}

// filter one row by one row
extern "C" __global__ void Prefilter_1D(double *coefficient, int length, double *pole, double tolerance, int nPoles)
{
	double Lambda;
	Lambda = 1;
	if (length == 1)
		return;

	/* compute the overall gain */
	for (int k = 0; k < nPoles; k++) {
		Lambda = Lambda * (1.0 - pole[k]) * (1.0 - 1.0 / pole[k]);
	}
	
	// Applying the gain to original image
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	if (col < length) {
		coefficient[col] *= Lambda;
	}
	
	
	
//  顺序执行性很强的语句如何处理？？？？？？？？？？？？？？？？？
	//for (k = 0; k < nPoles; k++)
	//{
	//	// Compute the first causal coefficient
	//	coefficient[0] = InitialCausalCoefficient(coefficient, length, pole[k], tolerance);

	//	// Causal prefilter
	//	for (n = 1; n < length; n++)
	//		coefficient[n] += pole[k] * coefficient[n - 1];

	//	//Compute the first anticausal coefficient
	//	coefficient[ length - 1] = InitialAnticausalCoefficient(coefficient, length, pole[k]);

	//	//Anticausal prefilter
	//	for (n = length - 2; n >= 0; n--)
	//		coefficient[n] = pole[k] * (coefficient[n + 1] - coefficient[n]);
	//}
}

extern "C" __global__ void Prefilter_1D_step2(double *coefficient, int length, double *pole, double tolerance, int nPoles)
{
	int k = 0;
	int n = 0;

	for (k = 0; k < nPoles; k++)
	{
		// Compute the first causal coefficient
		coefficient[0] = InitialCausalCoefficient(coefficient, length, pole[k], tolerance);

		// Causal prefilter
		for (n = 1; n < length; n++)
			coefficient[n] += pole[k] * coefficient[n - 1];

		//Compute the first anticausal coefficient
		coefficient[length - 1] = InitialAnticausalCoefficient(coefficient, length, pole[k]);

		//Anticausal prefilter
		for (n = length - 2; n >= 0; n--)
			coefficient[n] = pole[k] * (coefficient[n + 1] - coefficient[n]);
	}
}


extern "C" __device__ double InitialCausalCoefficient(double *sample, int length, double pole, double tolerance)
{

	double zn, iz, z2n;
	double FirstCausalCoef = 0;
	int n, horizon;
	horizon = (int)(ceil(log(tolerance) / log(fabs(pole))) + 0.01);
	if (horizon < length) {
		/* accelerated loop */
		zn = pole;
		FirstCausalCoef = *(sample);
		for (n = 1; n < horizon; n++) {
			FirstCausalCoef += zn * (*(sample + n));
			zn *= pole;
		}
	}
	else {
		/* full loop */
		zn = pole;
		iz = 1.0 / pole;
		z2n = pow(pole, (double)(length - 1));
		FirstCausalCoef = sample[0] + z2n * sample[length - 1];
		z2n *= z2n * iz;
		for (n = 1; n <= length - 2; n++) {
			FirstCausalCoef += (zn + z2n) * sample[n];
			zn *= pole;
			z2n *= iz;
		}
	}
	return FirstCausalCoef;
	
}


extern "C" __device__ double InitialAnticausalCoefficient(double *CausalCoef, int length, double pole)
{
	return (double)((pole / (pole * pole - 1.0)) * (pole * CausalCoef[length - 2] + CausalCoef[length - 1]));
}


extern "C" __global__ void GetParaFromImage(double * image_param, unsigned char *image, int width, int height)
{

	// 这个应该放在核函数外面生成

	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row < height && col < width) {
		image_param[row*width + col] = (double)image[row*width + col];
		// 原代码是一行一行，或一列一列的算，改成gpu后，应该是全部算完后，再按行列计算
	}

}

extern "C" __global__ void GetRow(double * image_param, double * row, int width, int height, int rownum) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	if (x < width) {
        row[x] = image_param[rownum*width + x];
	}
}

extern "C" __global__ void GetCol(double * image_param, double * col, int width, int height, int colnum) {
	int y = blockIdx.x * blockDim.x + threadIdx.x;

	if (y < height) {
		col[y] = image_param[colnum + y * width];
	}
}

extern "C" __global__ void GetParamFromRow(double * image_param, double * row, int width, int height, int rownum) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	if (x < width) {
		 image_param[rownum*width + x] = row[x];
	}
}


extern "C" __global__ void GetParamFromCol(double * image_param, double * col, int width, int height, int colnum) {
	int y = blockIdx.x * blockDim.x + threadIdx.x;

	if (y < height) {
		image_param[colnum + y * width] = col[y];
	}
}

extern "C" __host__ void GenPoles(int &nPoles, double* pole,double & gamma, double & a, int Interpolation_Algorithm)
{
	if (!pole) {
		return;
	}
	if (Interpolation_Algorithm == 1) // 4-tap
	{
		nPoles = 1;
		pole[0] = sqrt(3.0) - 2.0;
	}
	else if (Interpolation_Algorithm == 2) // 6-tap
	{
		nPoles = 2;
		pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
		pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
	}
	else if (Interpolation_Algorithm == 3) // modified 4-tap
	{
		gamma = 0.0409;
		a = (4.0 - 12.0 * gamma) / (6.0 * gamma + 1.0);
		pole[0] = (-a + sqrt(a * a - 4)) / 2.0;
	}
	else if (Interpolation_Algorithm == 4) // optimized 4-tap
	{
		nPoles = 1;
		pole[0] = (-13.0 + sqrt(105.0)) / 8.0;
	}
	else if (Interpolation_Algorithm == 5) // optimized 6-tap
	{
		nPoles = 2;
		pole[0] = -0.410549185795627524168;
		pole[1] = -0.0316849091024414351363;
	}
	else if (Interpolation_Algorithm == 6) // 8-tap
	{
		nPoles = 2;
		pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0) - 13.0 / 2.0;
		pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0) - 13.0 / 2.0;
	}
}

extern "C" __host__ void GenerateParaSpline(unsigned char* image, double * image_param_out,int width, int height, int Interpolation_Algorithm, double * pole, int nPoles, double a, double gamma) {

	double/* pole[2], a, gamma, */ tolerance = 1e-4;

	double *image_row/*[HEIGHT]*/;
	double *image_col/*[WIDTH]*/;
	double *image_param/*[WIDTH* HEIGHT]*/;
	double *poles;
	cudaMalloc(&image_row, sizeof(double)*WIDTH);
	cudaMalloc(&image_col, sizeof(double)*HEIGHT);
	cudaMalloc(&image_param, sizeof(double)*WIDTH* HEIGHT);
	cudaMalloc(&poles, sizeof(double)*nPoles);
	cudaMemcpy(poles, pole, sizeof(double)*nPoles, cudaMemcpyHostToDevice);

	dim3 Grid(ceil(width / (float)CELL_W), ceil(height / (float)CELL_H));
	dim3 Block(CELL_W, CELL_H);
	// char to double
	GetParaFromImage<<<Grid, Block>>>(image_param, image, width, height);

	for (int y = 0; y < height; y++) {
		dim3 Grid(ceil(width / (float)CELL_W));
		dim3 Block(CELL_W);
		GetRow <<<Grid, Block >> > (image_param, image_row, width, height, y);
		// function GetRow works

		if (Interpolation_Algorithm == 3) {
			
			Prefilter_1Dm<<<Grid, Block>>>(image_row, width, poles, tolerance, gamma);
		}
		else {
			/*if (y == 0) {
				double * temp = (double*)malloc(sizeof(double)*width);
				cudaMemcpy(temp, image_row, sizeof(double)*width, cudaMemcpyDeviceToHost);
				std::ofstream of;
				of.open("D:/Temp/GPU_Prefilter_1D_before.txt");
					for (int j = 0; j < width; j++)
					{
						of << (double)(*(temp + j)) << " ";
					}
				of.close();
				free(temp);
			}*/
			Prefilter_1D <<<Grid, Block >>> (image_row, width, poles, tolerance, nPoles); 
			
			/*if (y == 0) {
				double * temp = (double*)malloc(sizeof(double)*width);
				cudaMemcpy(temp, image_row, sizeof(double)*width, cudaMemcpyDeviceToHost);
				std::ofstream of;
				of.open("D:/Temp/GPU_Prefilter_1D_0.txt");
					for (int j = 0; j < width; j++)
					{
						of << (double)(*(temp + j)) << " ";
					}
				of.close();
				free(temp);
			}*/
			Prefilter_1D_step2 <<<1, 1 >>> (image_row, width, poles, tolerance, nPoles);
		}
        /*if (y == 0) {
			double * temp = (double*)malloc(sizeof(double)*width);
			cudaMemcpy(temp, image_row, sizeof(double)*width, cudaMemcpyDeviceToHost);
			std::ofstream of;
			of.open("D:/Temp/GPU_Prefilter_1D_step2_after_first_cal.txt");
				for (int j = 0; j < width; j++)
				{
					of << (double)(*(temp + j)) << " ";
				}
			of.close();
			free(temp);
		}*/


		GetParamFromRow <<<Grid, Block >>>(image_param_out, image_row, width, height, y);

		/*if (y == 0) {
			double * temp = (double*)malloc(sizeof(double)*width);
			cudaMemcpy(temp, image_row, sizeof(double)*width, cudaMemcpyDeviceToHost);
			std::ofstream of;
			of.open("D:/Temp/GPU_Prefileter.txt");
				for (int j = 0; j < width; j++)
				{
					of << (double)(*(temp + j)) << " ";
				}
			of.close();
			free(temp);
		}*/
	}

	for (int x = 0; x < width; x++) {
		dim3 Grid(ceil(height / (float)CELL_W));
		dim3 Block(CELL_W);
		GetCol <<<Grid, Block >>> (image_param_out, image_col, width, height, x);

		if (Interpolation_Algorithm == 3) {
			Prefilter_1Dm <<<Grid, Block >>> (image_col, height, poles, tolerance, gamma);
		}
		else {
			Prefilter_1D <<<Grid, Block >>> (image_col, height, poles, tolerance, nPoles);
			Prefilter_1D_step2<<<1,1>>> (image_col, height, poles, tolerance, nPoles);
		}

		GetParamFromCol <<<Grid, Block >>> (image_param_out, image_col, width, height, x);
	}

	cudaFree(image_row);
	cudaFree(image_col);
	cudaFree(image_param);

}

extern "C" __device__ void GetValueSpline(double *Para, int width, int height, double X, double Y, double *S, int S_Flag, int Interpolation_Algorithm)
{
	int i, j, width2, height2, xIndex[6], yIndex[6];
	double Para_Value, xWeight[6], yWeight[6], xWeightGradient[6], yWeightGradient[6], w, w2, w4, t, t0, t1, gamma;

	width2 = 2 * width - 2;
	height2 = 2 * height - 2;

	if (Interpolation_Algorithm == 6)
	{
		xIndex[0] = int(X) - 2;
		yIndex[0] = int(Y) - 2;
		for (i = 1; i < 6; i++)
		{
			xIndex[i] = xIndex[i - 1] + 1;
			yIndex[i] = yIndex[i - 1] + 1;
		}
	}
	else if ((Interpolation_Algorithm == 2) || (Interpolation_Algorithm == 5))
	{
		xIndex[0] = int(X + 0.5) - 2;
		yIndex[0] = int(Y + 0.5) - 2;
		for (i = 1; i < 5; i++)
		{
			xIndex[i] = xIndex[i - 1] + 1;
			yIndex[i] = yIndex[i - 1] + 1;
		}
	}
	else
	{
		xIndex[0] = int(X) - 1;
		yIndex[0] = int(Y) - 1;
		for (i = 1; i < 4; i++)
		{
			xIndex[i] = xIndex[i - 1] + 1;
			yIndex[i] = yIndex[i - 1] + 1;
		}
	}

	//Calculate the weights of x,y and their derivatives
	if (Interpolation_Algorithm == 1)
	{
		w = X - (double)xIndex[1];
		xWeight[3] = (1.0 / 6.0) * w * w * w;
		xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
		xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
		xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

		xWeightGradient[3] = w * w / 2.0;
		xWeightGradient[0] = w - 1.0 / 2.0 - xWeightGradient[3];
		xWeightGradient[2] = 1.0 + xWeightGradient[0] - 2.0 * xWeightGradient[3];
		xWeightGradient[1] = -xWeightGradient[0] - xWeightGradient[2] - xWeightGradient[3];

		/* y */
		w = Y - (double)yIndex[1];
		yWeight[3] = (1.0 / 6.0) * w * w * w;
		yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
		yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
		yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

		yWeightGradient[3] = w * w / 2.0;
		yWeightGradient[0] = w - 1.0 / 2.0 - yWeightGradient[3];
		yWeightGradient[2] = 1.0 + yWeightGradient[0] - 2.0 * yWeightGradient[3];
		yWeightGradient[1] = -yWeightGradient[0] - yWeightGradient[2] - yWeightGradient[3];
	}
	else if (Interpolation_Algorithm == 2)
	{
		w = X - (double)xIndex[2];
		w2 = w * w;
		t = (1.0 / 6.0) * w2;
		xWeight[0] = 1.0 / 2.0 - w;
		xWeight[0] *= xWeight[0];
		xWeight[0] *= (1.0 / 24.0) * xWeight[0];
		t0 = w * (t - 11.0 / 24.0);
		t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
		xWeight[1] = t1 + t0;
		xWeight[3] = t1 - t0;
		xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
		xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];

		xWeightGradient[0] = -(1.0 / 2.0 - w) * (1.0 / 2.0 - w) * (1.0 / 2.0 - w) / 6.0;
		xWeightGradient[1] = w * w / 2 - 11.0 / 24.0 + w / 2.0 - 2.0 * w * w * w / 3.0;
		xWeightGradient[3] = -w * w / 2 + 11.0 / 24.0 + w / 2.0 - 2.0 * w * w * w / 3.0;
		xWeightGradient[4] = xWeightGradient[0] + w * w / 2.0 + 1.0 / 24.0;
		xWeightGradient[2] = -xWeightGradient[0] - xWeightGradient[1] - xWeightGradient[3] - xWeightGradient[4];

		/* y */
		w = Y - (double)yIndex[2];
		w2 = w * w;
		t = (1.0 / 6.0) * w2;
		yWeight[0] = 1.0 / 2.0 - w;
		yWeight[0] *= yWeight[0];
		yWeight[0] *= (1.0 / 24.0) * yWeight[0];
		t0 = w * (t - 11.0 / 24.0);
		t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
		yWeight[1] = t1 + t0;
		yWeight[3] = t1 - t0;
		yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
		yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];

		yWeightGradient[0] = -(1.0 / 2.0 - w) * (1.0 / 2.0 - w) * (1.0 / 2.0 - w) / 6.0;
		yWeightGradient[1] = w * w / 2 - 11.0 / 24.0 + w / 2.0 - 2.0 * w * w * w / 3.0;
		yWeightGradient[3] = -w * w / 2 + 11.0 / 24.0 + w / 2.0 - 2.0 * w * w * w / 3.0;
		yWeightGradient[4] = yWeightGradient[0] + w * w / 2.0 + 1.0 / 24.0;
		yWeightGradient[2] = -yWeightGradient[0] - yWeightGradient[1] - yWeightGradient[3] - yWeightGradient[4];
	}
	else if (Interpolation_Algorithm == 3)
	{
		gamma = 0.0409;
		w = X - (double)xIndex[1];
		xWeight[0] = -w * w * w / 6.0 + w * w / 2.0 - (gamma + 0.5) * w + 1.0 / 6.0 + gamma;
		xWeight[1] = w * w * w / 2.0 - w * w + 3 * gamma * w + 2.0 / 3.0 - 2.0 * gamma;
		xWeight[2] = -w * w * w / 2.0 + w * w / 2.0 + (1.0 / 2.0 - 3.0 * gamma) * w + gamma + 1.0 / 6.0;
		xWeight[3] = w * w * w / 6.0 + gamma * w;

		xWeightGradient[0] = -w * w / 2.0 + w - gamma - 0.5;
		xWeightGradient[1] = 3.0 * w * w / 2.0 - 2.0 * w + 3.0 * gamma;
		xWeightGradient[2] = -3.0 * w * w / 2.0 + w + 1.0 / 2.0 - 3.0 * gamma;
		xWeightGradient[3] = w * w / 2.0 + gamma;

		/* y */
		w = Y - (double)yIndex[1];
		yWeight[0] = -w * w * w / 6.0 + w * w / 2.0 - (gamma + 0.5) * w + 1.0 / 6.0 + gamma;
		yWeight[1] = w * w * w / 2.0 - w * w + 3 * gamma * w + 2.0 / 3.0 - 2.0 * gamma;
		yWeight[2] = -w * w * w / 2.0 + w * w / 2.0 + (1.0 / 2.0 - 3.0 * gamma) * w + gamma + 1.0 / 6.0;
		yWeight[3] = w * w * w / 6.0 + gamma * w;

		yWeightGradient[0] = -w * w / 2.0 + w - gamma - 0.5;
		yWeightGradient[1] = 3.0 * w * w / 2.0 - 2.0 * w + 3.0 * gamma;
		yWeightGradient[2] = -3.0 * w * w / 2.0 + w + 1.0 / 2.0 - 3.0 * gamma;
		yWeightGradient[3] = w * w / 2.0 + gamma;
	}
	else if (Interpolation_Algorithm == 4)
	{
		w = X - (double)xIndex[1];
		xWeight[0] = -w * w * w / 6.0 + w * w / 2.0 - 11.0 * w / 21.0 + 4.0 / 21.0;
		xWeight[1] = w * w * w / 2.0 - w * w + 3.0 * w / 42.0 + 13.0 / 21.0;
		xWeight[2] = -w * w * w / 2.0 + w * w / 2.0 + 3.0 * w / 7.0 + 4.0 / 21.0;
		xWeight[3] = w * w * w / 6.0 + w / 42.0;

		xWeightGradient[0] = -w * w / 2.0 + w - 11.0 / 21.0;
		xWeightGradient[1] = 3.0 * w * w / 2.0 - 2.0 * w + 3.0 / 42.0;
		xWeightGradient[2] = -3.0 * w * w / 2.0 + w + 3.0 / 7.0;
		xWeightGradient[3] = w * w / 2.0 + 1.0 / 42.0;

		/* y */
		w = Y - (double)yIndex[1];
		yWeight[0] = -w * w * w / 6.0 + w * w / 2.0 - 11.0 * w / 21.0 + 4.0 / 21.0;
		yWeight[1] = w * w * w / 2.0 - w * w + 3.0 * w / 42.0 + 13.0 / 21.0;
		yWeight[2] = -w * w * w / 2.0 + w * w / 2.0 + 3.0 * w / 7.0 + 4.0 / 21.0;
		yWeight[3] = w * w * w / 6.0 + w / 42.0;

		yWeightGradient[0] = -w * w / 2.0 + w - 11.0 / 21.0;
		yWeightGradient[1] = 3.0 * w * w / 2.0 - 2.0 * w + 3.0 / 42.0;
		yWeightGradient[2] = -3.0 * w * w / 2.0 + w + 3.0 / 7.0;
		yWeightGradient[3] = w * w / 2.0 + 1.0 / 42.0;
	}
	else if (Interpolation_Algorithm == 5)
	{
		w = X - (double)xIndex[2];
		xWeight[0] = w * w * w * w / 24.0 - w * w * w / 12.0 + 11.0 * w * w / 144.0 - 5.0 * w / 144.0 + 743.0 / 120960.0;
		xWeight[1] = -w * w * w * w / 6.0 + w * w * w / 6.0 + 7.0 * w * w / 36.0 - 31.0 * w / 72.0 + 6397.0 / 30240.0;
		xWeight[2] = w * w * w * w / 4.0 - 13.0 * w * w / 24.0 + 11383.0 / 20160.0;
		xWeight[3] = -w * w * w * w / 6.0 - w * w * w / 6.0 + 7.0 * w * w / 36.0 + 31.0 * w / 72.0 + 6397.0 / 30240.0;
		xWeight[4] = w * w * w * w / 24.0 + w * w * w / 12.0 + 11.0 * w * w / 144.0 + 5.0 * w / 144.0 + 743.0 / 120960.0;
		//xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4]; 

		xWeightGradient[0] = w * w * w / 6.0 - w * w / 4.0 + 11.0 * w / 72.0 - 5.0 / 144.0;
		xWeightGradient[1] = -2.0 * w * w * w / 3.0 + w * w / 2.0 + 7.0 * w / 18.0 - 31.0 / 72.0;
		xWeightGradient[2] = w * w * w - 13.0 * w / 12.0;
		xWeightGradient[3] = -2.0 * w * w * w / 3.0 - w * w / 2.0 + 7.0 * w / 18.0 + 31.0 / 72.0;
		xWeightGradient[4] = w * w * w / 6.0 + w * w / 4.0 + 11.0 * w / 72.0 + 5.0 / 144.0;

		/* y */
		w = Y - (double)yIndex[2];
		yWeight[0] = w * w * w * w / 24.0 - w * w * w / 12.0 + 11.0 * w * w / 144.0 - 5.0 * w / 144.0 + 743.0 / 120960.0;
		yWeight[1] = -w * w * w * w / 6.0 + w * w * w / 6.0 + 7.0 * w * w / 36.0 - 31.0 * w / 72.0 + 6397.0 / 30240.0;
		yWeight[2] = w * w * w * w / 4.0 - 13.0 * w * w / 24.0 + 11383.0 / 20160.0;
		yWeight[3] = -w * w * w * w / 6.0 - w * w * w / 6.0 + 7.0 * w * w / 36.0 + 31.0 * w / 72.0 + 6397.0 / 30240.0;
		yWeight[4] = w * w * w * w / 24.0 + w * w * w / 12.0 + 11.0 * w * w / 144.0 + 5.0 * w / 144.0 + 743.0 / 120960.0;
		//yWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4]; 

		yWeightGradient[0] = w * w * w / 6.0 - w * w / 4.0 + 11.0 * w / 72.0 - 5.0 / 144.0;
		yWeightGradient[1] = -2.0 * w * w * w / 3.0 + w * w / 2.0 + 7.0 * w / 18.0 - 31.0 / 72.0;
		yWeightGradient[2] = w * w * w - 13.0 * w / 12.0;
		yWeightGradient[3] = -2.0 * w * w * w / 3.0 - w * w / 2.0 + 7.0 * w / 18.0 + 31.0 / 72.0;
		yWeightGradient[4] = w * w * w / 6.0 + w * w / 4.0 + 11.0 * w / 72.0 + 5.0 / 144.0;
	}
	else if (Interpolation_Algorithm == 6)
	{
		w = X - (double)xIndex[2];
		w2 = w * w;
		xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
		w2 -= w;
		w4 = w2 * w2;
		w -= 1.0 / 2.0;
		t = w2 * (w2 - 3.0);
		xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
		t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
		t1 = (-1.0 / 12.0) * w * (t + 4.0);
		xWeight[2] = t0 + t1;
		xWeight[3] = t0 - t1;
		t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
		t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
		xWeight[1] = t0 + t1;
		xWeight[4] = t0 - t1;

		xWeightGradient[5] = w * w * w * w / 24.0;
		xWeightGradient[0] = (4 * w * w * w - 6 * w * w + 4 * w - 1) / 24.0 - xWeightGradient[5];
		t0 = (4.0 * w * w * w - 6.0 * w * w - 8.0 * w + 5.0) / 24.0;
		t1 = -(5.0 * w * w * w * w - 10.0 * w * w * w - 3.0 * w * w + 8.0 * w + 5.0 / 2.0) / 12.0;
		xWeightGradient[2] = t0 + t1;
		xWeightGradient[3] = t0 - t1;
		t0 = (-4.0 * w * w * w + 6.0 * w * w + 4.0 * w - 3) / 16.0;
		t1 = (5.0 * w * w * w * w - 10.0 * w * w * w + 3.0 * w * w + 2 * w - 11.0 / 2.0) / 24.0;
		xWeightGradient[1] = t0 + t1;
		xWeightGradient[4] = t0 - t1;

		/* y */
		w = Y - (double)yIndex[2];
		w2 = w * w;
		yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
		w2 -= w;
		w4 = w2 * w2;
		w -= 1.0 / 2.0;
		t = w2 * (w2 - 3.0);
		yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
		t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
		t1 = (-1.0 / 12.0) * w * (t + 4.0);
		yWeight[2] = t0 + t1;
		yWeight[3] = t0 - t1;
		t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
		t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
		yWeight[1] = t0 + t1;
		yWeight[4] = t0 - t1;

		yWeightGradient[5] = w * w * w * w / 24.0;
		yWeightGradient[0] = (4 * w * w * w - 6 * w * w + 4 * w - 1) / 24.0 - yWeightGradient[5];
		t0 = (4.0 * w * w * w - 6.0 * w * w - 8.0 * w + 5.0) / 24.0;
		t1 = -(5.0 * w * w * w * w - 10.0 * w * w * w - 3.0 * w * w + 8.0 * w + 5.0 / 2.0) / 12.0;
		yWeightGradient[2] = t0 + t1;
		yWeightGradient[3] = t0 - t1;
		t0 = (-4.0 * w * w * w + 6.0 * w * w + 4.0 * w - 3) / 16.0;
		t1 = (5.0 * w * w * w * w - 10.0 * w * w * w + 3.0 * w * w + 2 * w - 11.0 / 2.0) / 24.0;
		yWeightGradient[1] = t0 + t1;
		yWeightGradient[4] = t0 - t1;
	}
	//***********************************

	/* apply the mirror boundary conditions and calculate the interpolated values */
	S[0] = 0;
	S[1] = 0;
	S[2] = 0;

	if (Interpolation_Algorithm == 6)
	{
		for (i = 0; i < 6; i++)
		{
			xIndex[i] = (width == 1) ? (0) : ((xIndex[i] < 0) ?
				(-xIndex[i] - width2 * ((-xIndex[i]) / width2))
				: (xIndex[i] - width2 * (xIndex[i] / width2)));
			if (width <= xIndex[i]) {
				xIndex[i] = width2 - xIndex[i];
			}
			yIndex[i] = (height == 1) ? (0) : ((yIndex[i] < 0) ?
				(-yIndex[i] - height2 * ((-yIndex[i]) / height2))
				: (yIndex[i] - height2 * (yIndex[i] / height2)));
			if (height <= yIndex[i]) {
				yIndex[i] = height2 - yIndex[i];
			}
		}

		for (i = 0; i < 6; i++)
		{
			for (j = 0; j < 6; j++)
			{
				Para_Value = (*(Para + width * yIndex[i] + xIndex[j]));
				S[0] = S[0] + Para_Value * xWeight[j] * yWeight[i];
				S[1] = S[1] + Para_Value * xWeightGradient[j] * yWeight[i];
				S[2] = S[2] + Para_Value * xWeight[j] * yWeightGradient[i];
			}
		}
	}
	else if ((Interpolation_Algorithm == 2) || (Interpolation_Algorithm == 5))
	{
		for (i = 0; i < 5; i++)
		{
			xIndex[i] = (width == 1) ? (0) : ((xIndex[i] < 0) ?
				(-xIndex[i] - width2 * ((-xIndex[i]) / width2))
				: (xIndex[i] - width2 * (xIndex[i] / width2)));
			if (width <= xIndex[i]) {
				xIndex[i] = width2 - xIndex[i];
			}
			yIndex[i] = (height == 1) ? (0) : ((yIndex[i] < 0) ?
				(-yIndex[i] - height2 * ((-yIndex[i]) / height2))
				: (yIndex[i] - height2 * (yIndex[i] / height2)));
			if (height <= yIndex[i]) {
				yIndex[i] = height2 - yIndex[i];
			}
		}
		for (i = 0; i < 5; i++)
		{
			for (j = 0; j < 5; j++)
			{
				Para_Value = (*(Para + width * yIndex[i] + xIndex[j]));
				S[0] = S[0] + Para_Value * xWeight[j] * yWeight[i];
				S[1] = S[1] + Para_Value * xWeightGradient[j] * yWeight[i];
				S[2] = S[2] + Para_Value * xWeight[j] * yWeightGradient[i];
			}
		}
	}
	else
	{
		for (i = 0; i < 4; i++)
		{
			xIndex[i] = (width == 1) ? (0) : ((xIndex[i] < 0) ?
				(-xIndex[i] - width2 * ((-xIndex[i]) / width2))
				: (xIndex[i] - width2 * (xIndex[i] / width2)));
			if (width <= xIndex[i]) {
				xIndex[i] = width2 - xIndex[i];
			}
			yIndex[i] = (height == 1) ? (0) : ((yIndex[i] < 0) ?
				(-yIndex[i] - height2 * ((-yIndex[i]) / height2))
				: (yIndex[i] - height2 * (yIndex[i] / height2)));
			if (height <= yIndex[i]) {
				yIndex[i] = height2 - yIndex[i];
			}
		}
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 4; j++)
			{
				Para_Value = (*(Para + width * yIndex[i] + xIndex[j]));
				S[0] = S[0] + Para_Value * xWeight[j] * yWeight[i];
				S[1] = S[1] + Para_Value * xWeightGradient[j] * yWeight[i];
				S[2] = S[2] + Para_Value * xWeight[j] * yWeightGradient[i];
			}
		}
	}

	return;
}


extern "C" __global__ void GetSplineData(char *image_data_N, int width_N, int height_N, double *para, int width, int height, int Interpolation_Algorithm) {

	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	float x = col / 2.25f;
	float y = row / 2.25f;
	
	double S[3];
	if (row < height_N && col < width_N)
	{
		GetValueSpline(para, width, height, x, y, S, 0, Interpolation_Algorithm);
	
		image_data_N[row*width_N + col] = (char)((int)(S[0] + 0.5));

	}
	
}

extern "C" __host__ void Process(int Interpolation_Algorithm)
{
	// 分配空间
	char *imagedata_;
	double *imageparam_;
	int width = WIDTH;
	int height = HEIGHT;
	int size = sizeof(char)*width*height;

	imagedata_ = (char*)malloc(size);
	imageparam_ = (double*)malloc(sizeof(double)*width*height);
	memset(imagedata_, 0, size);
	memset(imageparam_, 0, sizeof(double)*width*height);

	char *imagedata;
	double *imageparam;
	cudaMalloc(&imagedata, size);
	cudaMalloc(&imageparam, sizeof(double)*width*height);
	cudaMemset(imagedata, 0, size);
	cudaMemset(imageparam, 0, sizeof(double)*width*height);

	dim3 Grid(ceil(width / (float)CELL_W), ceil(height / (float)CELL_H));
	dim3 Block(CELL_W, CELL_H);
	InitImageData <<<Grid, Block >>> (imagedata, width, height);

	cudaMemcpy(imagedata_, imagedata, size, cudaMemcpyDeviceToHost);

	std::ofstream of;
	/*of.open("D:/Temp/GPU_initim.txt");
	if (of.is_open()) {
		for (int jj = 0; jj < height; jj++)
		{
			for (int ii = 0; ii < width; ii++)
			{
				of << (int)*(imagedata_ + jj * width + ii) << " ";
			}
			of << std::endl;
		}
		of.close();
	}*/

	int nPoles{0};
	double pole[2];
	double gamma{0};
	double a{0};

	GenPoles(nPoles, pole, gamma, a, Interpolation_Algorithm);

	GenerateParaSpline((unsigned char*)imagedata, imageparam, width, height, Interpolation_Algorithm, pole, nPoles, a, gamma);

	cudaMemcpy(imageparam_, imageparam, size, cudaMemcpyDeviceToHost);

	//of.open("D:/Temp/GPU_Para.txt");
	//if (of.is_open()) {
	//	for (int jj = 0; jj < height; jj++)
	//	{
	//		for (int ii = 0; ii < width; ii++)
	//		{
	//			of << (double)*(imageparam_ + jj * width + ii) << " ";
	//		}
	//		of << std::endl;
	//	}
	//	of.close();
	//}


	char *image_data_N;
	int width_N = WIDTH * 255 / 100.0f;
	int height_N = HEIGHT * 255 / 100.0f;
	cudaMalloc(&image_data_N, sizeof(char)*width_N*height_N);
	dim3 Grid_N(ceil(width_N / (float)CELL_W), ceil(height_N / (float)CELL_H));
	dim3 Block_N(CELL_W, CELL_H);
	GetSplineData<<<Grid_N, Block_N>>>(image_data_N, width_N, height_N, imageparam, width, height, Interpolation_Algorithm);
	
	char *image_data_N_c = (char*)malloc(sizeof(char)*width_N*height_N);

	cudaMemcpy(image_data_N_c, image_data_N,sizeof(char)*width_N*height_N, cudaMemcpyDeviceToHost);

	of.open("D:/Temp/testGPU.txt");
	if (of.is_open()) {
		for (int y = 0; y < height_N; y++) {
			for (int x = 0; x < width_N; x++) {
				//	printf("%d ", image_data_N_c[x + y * width_N]);
				of << (int)image_data_N_c[x + y * width_N] << " ";
			}
			of << std::endl;
			//printf("\n");
		}
		of.close();
	}
	// 释放空间

	cudaFree(image_data_N);
	cudaFree(imagedata);
	cudaFree(imageparam);

	free(imagedata_);
	free(imageparam_);
	free(image_data_N_c);
}

using namespace std;
#ifdef __cplusplus
extern "C" 
#endif
void DoAdd()
{

	float*A, *Ad, *B, *Bd, *C, *Cd;
	int n = 1024 * 1024 ;
	int size = n * sizeof(float);

	// CPU端分配内存
	A = (float*)malloc(size);
	B = (float*)malloc(size);
	C = (float*)malloc(size);

	memset(C, 0, size);
	// 初始化数组
	for (int i = 0; i < n; i++)
	{
		A[i] = 90.0;
		B[i] = 10.0;
		C[i] = 1;
	}

	auto time1 = std::chrono::high_resolution_clock::now();

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

	auto time2 = std::chrono::high_resolution_clock::now();
	// 校验误差
	float max_error = 0.0;
	for (int i = 0; i < n; i++)
	{
		max_error += fabs(100.0 - C[i]);
	}
	std::vector<float> v_c;
	v_c.assign(C+0, C+n-1);

	auto time3 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < n; i++)
	{
		C[i] = A[i]+B[i];
	}
	auto time4 = std::chrono::high_resolution_clock::now();

	cout << "GPU time cost = " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << "ms" << std::endl;
	cout << "CPU time cost = " << std::chrono::duration_cast<std::chrono::milliseconds>(time4 - time3).count() << "ms"<<std::endl;
	cout << "max error is " << max_error << endl;

	// 释放CPU端、GPU端的内存
	free(A);
	free(B);
	free(C);
	cudaFree(Ad);
	cudaFree(Bd);
	cudaFree(Cd);

	//return 0;
}

#define BLOCK_SIZE 16

#ifdef __cplusplus
extern "C"
#endif
void PicTest()
{
	const int m = 63; const int n = 76;
	float * image_in = new float[m*n];
	float * image_out = new float[m*n];
	for (int i = 0; i < m*n; i++) {
		image_in[i] = (float)i;
	}
	memset(image_out, 0, sizeof(float)*m*n);

	auto size = sizeof(float)*m*n;
	float* G_in;
	float* G_out;
	// GPU端分配内存
	cudaMalloc((void**)&G_in, size);
	cudaMalloc((void**)&G_out, size);
	cudaMemset(G_out, 0, size);

	cudaMemcpy(G_in, image_in, size, cudaMemcpyHostToDevice);
	//cudaMemcpy(G_out, image_out, size, cudaMemcpyHostToDevice);
	dim3 dimBlock(16,16);
	dim3 dimGrid(ceil(n / 16.0f), ceil(m / 16.0f));

	PictureKernel <<< dimGrid, dimBlock >>> (G_in, G_out, n, m);

	cudaMemcpy(image_out, G_out, size, cudaMemcpyDeviceToHost);

	for (int i = 0; i < m*n; i++) {
		printf(" %10.1f ", image_out[i]);
	}

	cudaFree(G_out);
	cudaFree(G_in);
	delete[] image_in;
	delete[] image_out;

}


#ifdef __cplusplus
extern "C"
#endif
void MatriMul()
{
	const int width = 5;
	float* m1 = new float[width*width];
	float* m2 = new float[width*width];
	float* m3 = new float[width*width];

	auto size = sizeof(float)*width*width;
	float* M1;
	float* M2;
	float* M3;

	cudaMalloc(&M1, size);
	cudaMalloc(&M2, size);
	cudaMalloc(&M3, size);

	for (int i = 0; i < width*width; i++) {
		m1[i] = (float)(i % 100);
		m2[i] = (float)(i % 88) ;
		m3[i] = (float)(0);
	}

	cudaMemcpy(M1, m1, size, cudaMemcpyHostToDevice);
	cudaMemcpy(M2, m2, size, cudaMemcpyHostToDevice);

	dim3 grid(ceil(width / (BLOCK_SIZE * 1.0f)), ceil(width / (BLOCK_SIZE *1.0f)), 1);
	dim3 block(BLOCK_SIZE, BLOCK_SIZE, 1);
	MatrixMulKernel <<< grid, block >>> (M1, M2, M3, width);

	cudaMemcpy(m3, M3, size, cudaMemcpyDeviceToHost);
	printf(" \n\n");
	for (int i = 0; i < width*width; i++) {
		printf(" %10.1f ", m3[i]);
	}

	cudaFree(M1);
	cudaFree(M2);
	cudaFree(M3);
	delete[] m1;
	delete[] m2;
	delete[] m3;


}