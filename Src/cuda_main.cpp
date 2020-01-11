#include <stdio.h>
#include <iostream>
#include "cuda_fun.cuh"

int main()
{
	std::cout << "Hello NVCC" << std::endl;
	
	Process(2);


	return 0;
}
