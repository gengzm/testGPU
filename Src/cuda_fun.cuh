#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

extern "C" __global__ void Plus(float A[], float B[], float C[], int n);

#ifdef __cplusplus
extern "C"
#endif
void DoAdd();


extern "C" __global__ void PictureKernel(float* in, float* out, int n, int m);

#ifdef __cplusplus
extern "C"
#endif
void PicTest();

extern "C" __global__ void MatrixMulKernel(float *dM, float * dN, float * dP, int width);

#ifdef __cplusplus
extern "C"
#endif
void MatriMul();


extern "C" __host__ void Process(int Interpolation_Algorithm);


extern "C" __global__ void InitImageData(char* image, int width, int height);


extern "C" __global__ void Prefilter_1Dm(double *coefficient, int length, double *pole, double tolerance, double gamma);

extern "C" __global__ void Prefilter_1D(double *coefficient, int length, double *pole, double tolerance, int nPoles);

extern "C" __device__ double InitialCausalCoefficient(double *sample, int length, double pole, double tolerance);


extern "C" __device__ double InitialAnticausalCoefficient(double *CausalCoef, int length, double pole);


extern "C" __global__ void GetParaFromImage(double * image_param, char *image, int width, int height);

extern "C" __global__ void GetRow(double * image_param, double * row, int width, int height, int rownum);

extern "C" __global__ void GetCol(double * image_param, double * col, int width, int height, int colnum);

extern "C" __global__ void GetParamFromRow(double * image_param, double * row, int width, int height, int rownum);

extern "C" __global__ void GetParamFromCol(double * image_param, double * col, int width, int height, int colnum);

extern "C" __host__ void GenPoles(int &nPoles, double* pole, double & gamma, double & a, int Interpolation_Algorithm);

extern "C" __host__ void GenerateParaSpline(char* image, double * image_param_out, int width, int height, int Interpolation_Algorithm, double * pole, double a, double gamma);

extern "C" __device__ void GetValueSpline(double *Para, int width, int height, double X, double Y, double *S, int S_Flag, int Interpolation_Algorithm);

extern "C" __global__ void GetSplineData(char *image_data_N, int width_N, int height_N, double *para, int width, int height, int Interpolation_Algorithm);