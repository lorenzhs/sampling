#pragma once

#include <cassert>
#include <limits>
#include <random>
#include <curand.h>
#include <cuda_profiler_api.h> //For profiling
#include <chrono> //For profiling
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/system/cuda/execution_policy.h>
#include <omp.h>
#include <map>


// convert a linear index to a row index
struct row_index : public thrust::unary_function<size_t,size_t>
{
	size_t c; // number of columns
	__host__ __device__
	row_index(size_t _c)
      : c(_c) 
	{ }

  __host__ __device__
  size_t operator()(size_t i)
  {
    return i / c;
  }
};

//Sample number with geometric distribution:
template <typename T>
struct geometric_distribution : public thrust::unary_function<T,T>
{
  #define DOUBLES
  #ifdef DOUBLES
  const double p;
  const double maxvalinv;
  #else
  const float p;
  const float maxvalinv;
  #endif

  geometric_distribution(const double _p, const double _maxvalinv)
    : p(_p),
    maxvalinv(_maxvalinv)
  { }

  __device__
  T operator()(T x)
  {
    return __log2f(1.0f-maxvalinv*x)/__log2f(1.0f-p)+1;
    //return log1p(-_maxvalinv*(double)x)/log1p(-_p);
  }
};

//initialize prng only once (seed can be set nevertheless)
template <typename T>
class cuda_generator
{
  private:
    static cuda_generator<T> *_instance;
	//static std::map<int, cuda_generator<T>*> _instance;
    thrust::device_vector<T> device_data;
	//T *output;
    size_t N;
    curandGenerator_t prng;
	cudaStream_t s;
  public:
    cuda_generator(size_t size)
    {
  fprintf(stderr, "creating instance\n");
  fflush(stderr);
	  cudaSetDevice(omp_get_thread_num());
      N = 0;
      device_data = thrust::device_vector<T>(size);
      N = size;
      //curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);
      //curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_MTGP32);
      curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_MT19937);
	  //output = 0;
	  cudaStreamCreate(&s);
	  curandSetStream(prng, s);
    }

    void call_curand_generate()
    {
      curandGenerate(prng, thrust::raw_pointer_cast(device_data.data()), N);
    }

  //public:
  static cuda_generator<T> *instance()
  {
	  int tid = omp_get_thread_num();
  fprintf(stderr, "returning instance %d\n", tid);
  fflush(stderr);
  //return new cuda_generator<T>();
  //#pragma omp critical
    if(_instance)
      _instance = new cuda_generator<T>();
    return _instance;
  }

  template <typename It>
  void generate_block(It dest, size_t size, double p,
                      unsigned long long seed = 0)
  {
    using value_type = typename std::iterator_traits<It>::value_type;			
    assert(p > 0 && p < 1);

    if (seed == 0) {
      seed = std::random_device{}();
    }

    //allocate memory on GPU if blocksize changed
    if(N != size)
    {
	  fprintf(stderr, "reallocating (%d --> %d)\n", N, size);
      device_data = thrust::device_vector<T>(size);
      N = size;
	  /*
	  if(output)
		  cudaFree(output);
	  cudaMalloc(&output, N*sizeof(T));
	  */
    }

    //set seed of prng
    curandSetPseudoRandomGeneratorSeed(prng, seed);

    //fill array with random numbers (should depend on size of the data type)
    call_curand_generate();

    //fprintf(stderr, ".");

    //geometric distribution and prefix sum
    auto tbegin = thrust::make_transform_iterator(device_data.begin(), geometric_distribution<T>(p, 1.0/(std::numeric_limits<T>::max())));
    auto tend = thrust::make_transform_iterator(device_data.end(), geometric_distribution<T>(p, 1.0/(std::numeric_limits<T>::max())));
	thrust::transform(device_data.begin(), device_data.end(), device_data.begin(), geometric_distribution<T>(p, 1.0/(std::numeric_limits<T>::max())));
	#if 0
	//thrust::inclusive_scan(thrust::system::cuda::par, tbegin, tend, device_data.begin());// Determine temporary device storage requirements for inclusive prefix sum
	void     *d_temp_storage = NULL;
	size_t   temp_storage_bytes = 0;
	cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, &device_data[0], output, N,0,1);
	// Allocate temporary storage for inclusive prefix sum
	cudaMalloc(&d_temp_storage, temp_storage_bytes);
	// Run inclusive prefix sum
	cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, &device_data[0], output, N,0,1);
	cudaFree(d_temp_storage);
	#else
	#if 0
    {
	  //thrust::inclusive_scan(thrust::system::cuda::par, device_data.begin(), device_data.end(), output);
      thrust::counting_iterator<size_t> indices(0);
      //auto tbegin = thrust::make_transform_iterator(indices, row_index<size_t>(N/2));
      //auto tend = thrust::make_transform_iterator(indices, row_index<size_t>(N/2))+N*sizeof(T);
      //thrust::inclusive_scan_by_key(tbegin, tend, device_data.begin(), device_data.begin());
	  for(int n = N; n >= 1; n >>= 1)
	  {
		  /*
		  thrust::inclusive_scan_by_key(thrust::make_transform_iterator(indices, row_index(N/2)),
			thrust::make_transform_iterator(indices, row_index(N/2)) + N,
			device_data.begin(),
			device_data.begin());
		  */
		  auto row_iterator = thrust::make_transform_iterator(indices, row_index(n));
		  std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
		  for(int i = 0; i < 100; i++)
			  thrust::inclusive_scan_by_key(row_iterator,
				row_iterator + N,
				device_data.begin(),
				device_data.begin());
		  float dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count();
		  fprintf(stderr, "%d %f\n", n, dt/100);
	  }
    }
	#else
    {
	  thrust::inclusive_scan(thrust::system::cuda::par.on(s), device_data.begin(), device_data.end(), device_data.begin());
    }
	#endif
	#endif
	
    //read back data
    //cudaMemcpy(&dest[0], device_data.data().get(), N*sizeof(T), cudaMemcpyDeviceToHost);
	//cudaMemcpy(&dest[0], output, N*sizeof(T), cudaMemcpyDeviceToHost);
	//cudaDeviceSynchronize();
  }
};
//template <typename T> std::map<int, cuda_generator<T>*> cuda_generator<T>::_instance;
template <typename T> cuda_generator<T>* cuda_generator<T>::_instance;
/*
//TODO: not working!!! curand complains about not using a 64-bit prng.
template <>
void cuda_generator<unsigned long long>::call_curand_generate()
{
  curandGenerateLongLong(prng, thrust::raw_pointer_cast(device_data.data()), N);
}
*/

struct cuda_gen {
  template <typename It>
  static void generate_block(It dest, size_t size, double p,
                             unsigned int seed = 0)
  {
    using value_type = typename std::iterator_traits<It>::value_type;
    cuda_generator<value_type> *generator = cuda_generator<value_type>::instance();
    generator->generate_block(dest, size, p, seed);
  }

  template <typename It>
  static void generate_block(It begin, It end, double p,
                             unsigned int seed = 0)
  {
    generate_block(begin, end-begin, p, seed);
  }
};

#if 1
//For profiling
int main(int argc, char **argv)
{
  unsigned int N = 1024*1024*512;
  std::vector<unsigned int> vec(N);
  fprintf(stderr, "\n");
  //cudaSetDevice(1);
  do {
	  int num_threads = 1;
	  std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
	  #pragma omp parallel num_threads(1)
	  {
		  cuda_generator<unsigned int> *generator = new cuda_generator<unsigned int>(N/num_threads);
		  cudaDeviceSynchronize();
  cudaProfilerStart();
	  for(int i = 0; i < 10; i++)
		generator->generate_block(vec.begin()+omp_get_thread_num()*N/num_threads, N/num_threads, 0.5);
  cudaProfilerStop();
	  }
	  //N /= 2;
	  float dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count();
	  fprintf(stderr, "%d %f\n", N, dt/10);
  } while(N > N);
  //fprintf(stderr, "done\n");
  double avg = 0;
  avg = vec[N-1]/(double)N;
  fprintf(stderr, "avg: %f\n", avg);
  return 0;
}
#endif
