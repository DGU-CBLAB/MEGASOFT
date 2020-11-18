#include"MetaGPU.cuh"


__device__ __host__  void gpu_logBeta(double* res, double m, double n)
{
	// Beta Function on CUDA
	// https://stackoverflow.com/questions/15158297/beta-pdf-function-for-cuda/15159945
	*res = log(exp(lgamma(m) + lgamma(n) - lgamma(n + m)));
	//*res = log(boost::math::beta(m, n));
	return;
}
__device__ __host__ void gpu_chiSquareComplemented(double* res, double v, double x) 
{
	// Incomplete Gamma Function
	// https://www.cs.rit.edu/~ark/lectures/cuda01/c2html.php?file=9
	*res = 1 - (gammp(v / 2.0, x / 2.0));
	return;
}