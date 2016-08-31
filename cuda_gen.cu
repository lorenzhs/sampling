#include <cassert>
#include <limits>
#include <random>
#include <cuda.h>
#include <curand_kernel.h>
#include <curand.h>
#include <cuda_profiler_api.h>
#include <chrono>

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

//initialize curand states
__global__ void initCurand(curandState *state, unsigned long long *seed)
{
  unsigned int threadId = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(seed[threadId], 0, 0, &state[threadId]);
}

//sample random numbers with geometric distribution
template <typename T>
__global__ void sample(curandState *state, T *data, const T num, const double oneminusp)
{
  unsigned int threadId = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int idx = threadId * num;
  curandState mState = state[threadId];
  //each thread generates 'num' random numbers with geometric distribution
  for(unsigned int i = 0; i < num; i++, idx++)
    data[idx] = __log2f(1.0f-curand_uniform(&mState))/__log2f(oneminusp)+1;
  state[threadId] = mState;
}

//naive way of removing excessive indices
template <typename T>
__global__ void naive(curandState *state, T *data, const T *split, const T *to_remove)
{
  unsigned int threadId = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int idx = split[threadId];
  unsigned int num = split[threadId+1]-split[threadId];
  curandState mState = state[threadId];

  for(unsigned int i = 0; i < to_remove[threadId]; i++)
  {
    //pick a random element and remove it if it wasn't removed before
    unsigned int r = curand_uniform(&mState)*num;
    if(data[idx+r] == (T)-1)
      --i;
    data[idx+r] = (T)-1;
  }
  state[threadId] = mState;
}

//Vitter's Algorithm A (base case of Algorithm D)
template <typename T>
__global__ void algorithmA(curandState *state, T *data, const T *split, const T *to_remove)
{
  unsigned int threadId = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int idx = split[threadId];
  unsigned int num = split[threadId+1]-split[threadId];
  curandState mState = state[threadId];

  float top = num - to_remove[threadId];
  float Nreal = num;
  float r = curand_uniform(&mState);  
  float quot = top/Nreal;

  for(unsigned int i = 0; i < num; i++, idx++)
  {
    if(quot <= r)
    {
      r = curand_uniform(&mState);
      data[idx+1] = (T)-1;
      quot = 1;
    } else
      top--;
    Nreal--;
    quot *= top/Nreal;
  }
  state[threadId] = mState;
}

//functor for filtering elements less or equal a constant
template <typename T>
struct less_equal_functor
{

  T max_idx;
  less_equal_functor(const T _max_idx)
    : max_idx(_max_idx)
  { }

  __host__ __device__
    bool operator()(T x)
    {
      return x <= max_idx;
    }
};

//functor for filtering elements greater than a constant
template <typename T>
struct greater_than_functor
{
  T max_idx;
  greater_than_functor(const T _max_idx)
    : max_idx(_max_idx)
  { }

  __host__ __device__
    bool operator()(T x)
    {
      return x > max_idx;
    }
};

//initialize prng only once (seed can be set nevertheless)
template <typename T>
class cuda_generator
{
  private:
    static cuda_generator<T> *_instance;

    T *device_data;

    //(GPU-)states for random number generators
    curandState *state;
    //seeds on CPU and GPU for each thread of the random number generator
    unsigned long long *seeds, *seeds_device;
    //
    size_t N;
    const T num_threads = 1<<10;
    const T num_threads_per_group = 1<<4;
    const T num_groups = num_threads / num_threads_per_group;

    cuda_generator()
    {
      N = 0;
      cudaMalloc((void**)&state, num_threads * sizeof(curandState));
      cudaMalloc((void**)&seeds_device, num_threads * sizeof(unsigned long long));
      seeds = new unsigned long long[num_threads];
      device_data = 0;
    }


  public:
    static cuda_generator<T> *instance()
    {
      if(!_instance)
        _instance = new cuda_generator<T>();
      return _instance;
    }
	
	//first call: num_samples, N, 0, num_threads, ...
	static void sampleR(T n, T N, T j, T k, T *startIdx, T *numSamples)
	{
		if(k == j+1)
		{
			numSamples[j] = n;
			startIdx[j] = N;
			return;
		}
		T N2 = N/2;
		T x = n/2;//StochasticLib1::Hypergeometric(n, N2, N);
		sampleR(x, N2, j, (j+k)/2, startIdx, numSamples);
		sampleR(n-x, N-N2, (j+k)/2, k, startIdx, numSamples);
		for(int i = (j+k)/2; i < k; i++)
		  startIdx[i] += N2;
	}

    template <typename It>
      void generate_block(It dest, size_t size, double p,
          unsigned long long seed = 0)
      {
        using value_type = typename std::iterator_traits<It>::value_type;
        assert(p > 0 && p < 1);

        if (seed == 0)
        {
          seed = std::random_device{}();
        }

        //allocate memory on GPU if blocksize changed
        if(N != size)
        {
          if(device_data)
            cudaFree(device_data);
          cudaMalloc((void**)&device_data, sizeof(T)*size);
          N = size;
        }
        for(int i = 0; i < num_threads; i++)
          seeds[i] = seed+i;
        cudaMemcpy(seeds_device, seeds, num_threads*sizeof(unsigned long long), cudaMemcpyHostToDevice);
        //initialize prng-states
        initCurand<<<num_groups, num_threads_per_group>>>(state, seeds_device);
        //fill array with geometrically distributed random numbers
        sample<<<num_groups, num_threads_per_group>>>(state, device_data, (T)(N/num_threads), 1.0f-p);

        //prefix-sum for calculating array-indices
        thrust::inclusive_scan(thrust::device, device_data, device_data+N, device_data);

        //TODO: get universe size from command line parameter
        T universe = 1.5*size;

        //count number of indices that are small enough
        T num_indices = thrust::count_if(thrust::device, device_data, device_data+N, less_equal_functor<T>(universe));

        //divide array using hypergeometric deviates (one subarray per processor)
        T n = 1<<28;
		#if 0
        //TODO!!! - for now just take subarrays of equal size
        //intervals for each thread:
        T *split = new T[num_threads+1];
        T *to_remove = new T[num_threads];
        T remove_remaining = num_indices-n;
        T remove_total = 0;
        for(int i = 0; i < num_threads+1; i++)
        {
          split[i] = (double)i*n/(double)num_threads+0.5;
          unsigned int thread_removes = (remove_remaining)/(double)(num_threads-i);
          remove_total += thread_removes;
          to_remove[i] = thread_removes;
          remove_remaining -= thread_removes;
        }
		#else
        T *split = new T[num_threads+1];
        T *to_remove = new T[num_threads];
		
		sampleR(n, N, 0, num_threads, split+1, to_remove);
		split[0] = 0;
		/*
		for(int i = 0; i <= num_threads; i++)
			fprintf(stderr, "%d ", split[i]);
		for(int i = 0; i < num_threads; i++)
			fprintf(stderr, "%d ", to_remove[i]);
		*/
		#endif

        //remove excessive elements per processor
        T *device_split, *device_to_remove;
        //copy start- and end-indices of subarrays as well as the number
        //of elements to be sampled from these subarrays to GPU memory
        cudaMalloc((void**)&device_split, sizeof(T)*(num_threads+1));
        cudaMalloc((void**)&device_to_remove, sizeof(T)*num_threads);
        cudaMemcpy(device_split, split, sizeof(T)*(num_threads+1), cudaMemcpyHostToDevice);
        cudaMemcpy(device_to_remove, to_remove, sizeof(T)*num_threads, cudaMemcpyHostToDevice);

        //sample elements (elements are removed by setting the index to some value larger
        //than the universe size
        //naive<<<num_groups, num_threads_per_group>>>(state, device_data, device_split, device_to_remove);
        algorithmA<<<num_groups, num_threads_per_group>>>(state, device_data, device_split, device_to_remove);

        //compaction
        T final_count = thrust::copy_if(thrust::device, device_data, device_data+num_indices, device_data, less_equal_functor<T>(universe)) - device_data;
		//T final_count = thrust::remove_copy(thrust::device, device_data, device_data+num_indices, device_data, (T)-1) - device_data;
		fprintf(stderr, "result: %lu / %lu elements\n", final_count, n);
        //read back data
        //cudaMemcpy(&dest[0], device_data, N*sizeof(T), cudaMemcpyDeviceToHost);
      }
};
template <typename T> cuda_generator<T> *cuda_generator<T>::_instance = 0;

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

int main(int argc, char **argv)
{
  unsigned int N = 1024*1024*512;
  std::vector<unsigned int> vec(N);
  fprintf(stderr, "\n");
  cudaSetDevice(3);
  do 
  {
    //unsigned int num_threads = 1;
        cuda_gen::generate_block(vec.begin(), N, 0.5);
    std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
    //#pragma omp parallel num_threads(1)
    {
      //cuda_generator<unsigned int> *generator = new cuda_generator<unsigned int>(N/num_threads);
      cudaDeviceSynchronize();
      cudaProfilerStart();
      for(int i = 0; i < 10; i++)
        //generator->generate_block(vec.begin()/*+omp_get_thread_num()*N/num_threads*/, N/num_threads, 0.5);
        cuda_gen::generate_block(vec.begin(), N, 0.5);
      cudaProfilerStop();
      cudaDeviceSynchronize();
    }
    //N /= 2;
    float dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t).count();
    fprintf(stderr, "%d %f\n", N, dt/10);
  } while(N > N);
  //fprintf(stderr, "done\n");
  double avg = 0;
  avg = vec[N-1]/(double)N;
  fprintf(stderr, "avg: %f\n", avg);
  for(int i = 0; i < 100; i++)
  {
	fprintf(stderr, "%lu ", vec[i]);
  }
  return 0;
}
