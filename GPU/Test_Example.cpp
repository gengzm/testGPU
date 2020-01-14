#include "Spline_Interpolation.cpp"

#include <fstream>

void Test()
{
	int ii, jj, width = 204, height = 153, length = width * height;
	int width_N = width * 225 / 100, height_N = height * 225 / 100, length_N = width_N * height_N;
	int Interpolation_Algorithm = 1; // 1~6
	double x, y, U, tt, S[3];
	char *Img_data = new char[length];
	char *Img_data_N = new char[length_N];
	double *Para = new double[length];

	// Simulate an image
	for (jj = 0; jj < height; jj++)
	{
		for (ii = 0; ii < width; ii++)
		{

			double x = ((ii - width /2) / (double)width);
			double y = ((jj - height /2) / (double)height);

			U = 3.0*(1 - x)*(1 - x)*exp(-x * x - (y + 1)*(y + 1)) - 10.0*(x / 5.0 - x * x*x - y * y*y*y*y)*exp(-x * x - y * y) - 1.0 / 3.0*exp(-(x + 1)*(x + 1) - y * y);

			tt = 127.5*(1.0 + cos(200 * U));
			*(Img_data + jj * width + ii) = (char)((int)(tt + 0.5));
		}
	}

	std::ofstream of;
	of.open("D:/Temp/CPU_initim.txt");
	if (of.is_open()) {
		for (jj = 0; jj < height; jj++)
		{
			for (ii = 0; ii < width; ii++)
			{
				of << (int)*(Img_data + jj * width + ii) << " ";
			}
			of << std::endl;
		}
		of.close();
	}


	//	SaveDataToGreyImage(Img_data, width, height, "C:\\Temp\\000.bmp");

	//  强制转换为uchar， -2变成254
	Generate_Para_Spline((unsigned char*)Img_data, Para, width, height, Interpolation_Algorithm);

	of.open("D:/Temp/CPU_Para.txt");
	if (of.is_open()) {
		for (jj = 0; jj < height; jj++)
		{
			for (ii = 0; ii < width; ii++)
			{
				of << *(Para + jj * width + ii) << " ";
			}
			of << std::endl;
		}
		of.close();
	}

	
	of.open("D:/Temp/testCPU.txt");
	// Zoom-in 2.25X
	for (jj = 0; jj < height_N; jj++)
	{
		for (ii = 0; ii < width_N; ii++)
		{
			x = ii / 2.25;
			y = jj / 2.25;

			Get_Value_Spline(Para, width, height, x, y, S, 0, Interpolation_Algorithm); // S_Flag is a reserved para.

			*(Img_data_N + jj * width_N + ii) = (char)((int)(S[0] + 0.5));

			of << (int)Img_data_N[ii + jj * width_N] << " ";

			//printf("%d ", Img_data_N[ii + jj * width_N]);
		}

		of << std::endl;
		//printf("\n");
	}
	of.close();
	
	//	SaveDataToGreyImage(Img_data_N, width_N, height_N, "C:\\Temp\\111.bmp");

	delete[]Para;
	delete[]Img_data_N;
	delete[]Img_data;
	return;
}

int main()
{
	Test();
	return 0;

}