//// Spline interpolation functions

#include <stdio.h>
#include <iostream>
#include <fstream>

//// To generate interpolation parametersm, use: 
// template <class m_Type> void Generate_Para_Spline(m_Type *Image, double *Para, int width, int height, int Interpolation_Algorithm);

//// To get the value at point (X, Y), use: 
// extern "C" __declspec(dllexport) void Get_Value_Spline(double *Para, int width, int height, double X, double Y, double *S, int S_Flag, int Interpolation_Algorithm);


extern "C" __declspec(dllexport) double InitialCausalCoefficient(double *sample, int length, double pole, double tolerance)
{
	//AFX_MANAGE_STATE(AfxGetStaticModuleState());

	double zn, iz, z2n;
	double FirstCausalCoef;
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


extern "C" __declspec(dllexport) double InitialAnticausalCoefficient(double *CausalCoef, int length, double pole)
{
	//AFX_MANAGE_STATE(AfxGetStaticModuleState());

	return((pole / (pole * pole - 1.0)) * (pole * CausalCoef[length - 2] + CausalCoef[length - 1]));
}

// prefilter for 4-tap, 6-tap, 8-tap, optimized 4-tap, and optimized 6-tap
extern "C" __declspec(dllexport) void Prefilter_1D(double *coefficient, int length, double *pole, double tolerance, int nPoles)
{
	//AFX_MANAGE_STATE(AfxGetStaticModuleState());

	int i, n, k;
	double Lambda;
	Lambda = 1;
	if (length == 1)
		return;
	/* compute the overall gain */
	for (k = 0; k < nPoles; k++)
		Lambda = Lambda * (1.0 - pole[k]) * (1.0 - 1.0 / pole[k]);

	
	// Applying the gain to original image
	for (i = 0; i < length; i++) {
		*(coefficient + i) = (*(coefficient + i)) * Lambda;
	}


	for (k = 0; k < nPoles; k++)
	{
		// Compute the first causal coefficient
		*(coefficient) = InitialCausalCoefficient(coefficient, length, pole[k], tolerance);

		// Causal prefilter
		for (n = 1; n < length; n++)
			coefficient[n] += pole[k] * coefficient[n - 1];

		//Compute the first anticausal coefficient
		*(coefficient + length - 1) = InitialAnticausalCoefficient(coefficient, length, pole[k]);

		//Anticausal prefilter
		for (n = length - 2; n >= 0; n--)
			coefficient[n] = pole[k] * (coefficient[n + 1] - coefficient[n]);
	}


	static int si = 0;
	if (si == 0) {
		std::ofstream of;
		of.open("D:/Temp/CPU_Prefilter_1D_after_first_cal.txt");
		for (int j = 0; j < length; j++)
		{
			of << (double)(*(coefficient + j)) << " ";
		}
		of.close();

	}
	si++;
}

// Prefilter for modified 4-tap
extern "C" __declspec(dllexport) void Prefilter_1Dm(double *coefficient, int length, double *pole, double tolerance, double gamma)
{
	//AFX_MANAGE_STATE(AfxGetStaticModuleState());

	int i, n, k;
	double Lambda;
	Lambda = 6.0 / (6.0 * gamma + 1.0);
	if (length == 1)
		return;

	// Applying the gain to original image
	for (i = 0; i < length; i++)
		*(coefficient + i) = (*(coefficient + i)) * Lambda;

	for (k = 0; k < 1; k++)	// Testing
	{
		// Compute the first causal coefficient
		*(coefficient) = InitialCausalCoefficient(coefficient, length, pole[k], tolerance);

		// Causal prefilter
		for (n = 1; n < length; n++)
			coefficient[n] += pole[k] * coefficient[n - 1];

		//Compute the first anticausal coefficient
		*(coefficient + length - 1) = InitialAnticausalCoefficient(coefficient, length, pole[k]);

		//Anticausal prefilter
		for (n = length - 2; n >= 0; n--)
			coefficient[n] = pole[k] * (coefficient[n + 1] - coefficient[n]);
	}
}

template <class m_Type> void Generate_Para_Spline(m_Type *Image, double *Para, int width, int height, int Interpolation_Algorithm)
{
	//AFX_MANAGE_STATE(AfxGetStaticModuleState());

	int i, j, nPoles, length = width * height;
	double pole[2], a, gamma, tolerance = 1e-4;

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

	//Perform the 1D prefiltering along the rows
	double *LineWidth = new double[width];
	for (i = 0; i < height; i++)
	{
		//Prefiltering each row
		for (j = 0; j < width; j++)
		{
			*(LineWidth + j) = (double)(*(Image + i * width + j));
		}

		{
			if (i == 0) {
				std::ofstream of;
				of.open("D:/Temp/CPU_Prefilter_1D_before_and_out_image.txt");
				//for (int y = 0; y < HEIGHT; y++)
				{
					for (int x = 0; x < width; x++) {
						of << (double)Image[x] << " ";
					}
					of << std::endl;
				}
				of.close();

				of.open("D:/Temp/CPU_Prefilter_1D_before_and_out.txt");
				//for (int y = 0; y < HEIGHT; y++)
				{
					for (int x = 0; x < width; x++) {
						of << LineWidth[ x] << " ";
					}
					of << std::endl;
				}

				of.close();
			}
		}
		if (Interpolation_Algorithm == 3)
			Prefilter_1Dm(LineWidth, width, pole, tolerance, gamma);
		else
			Prefilter_1D(LineWidth, width, pole, tolerance, nPoles);


		
		// Put the prefiltered coeffiecients into Para array
		for (j = 0; j < width; j++)
		{
			*(Para + i * width + j) = (*(LineWidth + j));
		}
		
		if (i == 0) {
			std::ofstream of;
			of.open("D:/Temp/CPU_Prefileter.txt");
				for (j = 0; j < width; j++)
				{
					of << (*(LineWidth + j)) << " ";
				}
				of.close();
		}
	}
	delete[]LineWidth;

	//Perform the 1D prefiltering along the columns
	double *LineHeight = new double[height];
	for (i = 0; i < width; i++)
	{
		//Prefiltering each comlumn
		for (j = 0; j < height; j++)
		{
			*(LineHeight + j) = (*(Para + j * width + i));
		}
		if (Interpolation_Algorithm == 3)
			Prefilter_1Dm(LineHeight, height, pole, tolerance, gamma);
		else
			Prefilter_1D(LineHeight, height, pole, tolerance, nPoles);

		//Put the prefilterd coefficients into the Para array
		for (j = 0; j < height; j++)
		{
			*(Para + j * width + i) = (*(LineHeight + j));
		}
	}
	delete[]LineHeight;
	return;
}

template __declspec(dllexport) void Generate_Para_Spline(unsigned char *Image, double *Para, int width, int height, int Interpolation_Algorithm);
template __declspec(dllexport) void Generate_Para_Spline(int *Image, double *Para, int width, int height, int Interpolation_Algorithm);
template __declspec(dllexport) void Generate_Para_Spline(double *Image, double *Para, int width, int height, int Interpolation_Algorithm);

template <class m_Type> void Generate_Para_Spline_Parallel(m_Type *Image, double *Para, int width, int height, int Interpolation_Algorithm)
{
	//AFX_MANAGE_STATE(AfxGetStaticModuleState());
#if 0
	int nPoles, length = width * height;
	double pole[2], a, gamma, tolerance = 1e-4;

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

	//Perform the 1D prefiltering along the rows
	concurrency::parallel_for(0, height, [&](int i)
	{
		double *LineWidth = new double[width];
		int j;
		//Prefiltering each row
		for (j = 0; j < width; j++)
		{
			*(LineWidth + j) = (double)(*(Image + i * width + j));
		}
		if (Interpolation_Algorithm == 3)
			Prefilter_1Dm(LineWidth, width, pole, tolerance, gamma);
		else
			Prefilter_1D(LineWidth, width, pole, tolerance, nPoles);

		// Put the prefiltered coeffiecients into Para array
		for (j = 0; j < width; j++)
		{
			*(Para + i * width + j) = (*(LineWidth + j));
		}

		delete[]LineWidth;
	});

	//Perform the 1D prefiltering along the columns
	concurrency::parallel_for(0, width, [&](int i)
	{
		double *LineHeight = new double[height];
		int j;
		//Prefiltering each comlumn
		for (j = 0; j < height; j++)
		{
			*(LineHeight + j) = (*(Para + j * width + i));
		}
		if (Interpolation_Algorithm == 3)
			Prefilter_1Dm(LineHeight, height, pole, tolerance, gamma);
		else
			Prefilter_1D(LineHeight, height, pole, tolerance, nPoles);

		//Put the prefilterd coefficients into the Para array
		for (j = 0; j < height; j++)
		{
			*(Para + j * width + i) = (*(LineHeight + j));
		}
		delete[]LineHeight;
	});
#endif
	return;
}

template __declspec(dllexport) void Generate_Para_Spline_Parallel(unsigned char *Image, double *Para, int width, int height, int Interpolation_Algorithm);
template __declspec(dllexport) void Generate_Para_Spline_Parallel(int *Image, double *Para, int width, int height, int Interpolation_Algorithm);
template __declspec(dllexport) void Generate_Para_Spline_Parallel(double *Image, double *Para, int width, int height, int Interpolation_Algorithm);

extern "C" __declspec(dllexport) void Get_Value_Spline(double *Para, int width, int height, double X, double Y, double *S, int S_Flag, int Interpolation_Algorithm)
{
	//AFX_MANAGE_STATE(AfxGetStaticModuleState());

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
