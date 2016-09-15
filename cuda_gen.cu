#include <cassert>
#include <limits>
#include <random>
#include <cuda.h>
#include <curand_kernel.h>
#include <curand.h>
#ifdef PROFILING
#include <cuda_profiler_api.h>
#endif
#include "util.h"
#include "timer.h"
#include <stocc/stocc.h>
#include <randomc/randomc.h>

#include <thrust/copy.h>
#include <thrust/find.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

#include <include/arg_parser.h>

#include <algorithm>
#include <cmath>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323814
#endif

//initialize curand states
__global__ void initCurand(curandState *state, unsigned long long *seed)
{
  size_t threadId = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(seed[threadId], 0, 0, &state[threadId]);
}

//sample random numbers with geometric distribution
  template <typename T>
__global__ void sample(curandState *state, T *data, const size_t count, const T num, const double oneminusp)
{
  size_t threadId = blockIdx.x * blockDim.x + threadIdx.x;
  size_t idx = threadId * num;
  curandState mState = state[threadId];
  //each thread generates 'num' random numbers with geometric distribution
  float lpinv = 1/log2(oneminusp);
  for(unsigned int i = 0; i < num && idx < count; i+=4, idx+=4)
  {
    float4 r = make_float4(__log2f(1-curand_uniform(&mState)),__log2f(1-curand_uniform(&mState)),__log2f(1-curand_uniform(&mState)),__log2f(1-curand_uniform(&mState)));
    if(i < num && idx < count)
      data[idx] = r.x*lpinv+1;
    if(i+1 < num && idx+1 < count)
      data[idx+1] = r.y*lpinv+1;
    if(i+2 < num && idx+2 < count)
      data[idx+2] = r.z*lpinv+1;
    if(i+3 < num && idx+3 < count)
      data[idx+3] = r.w*lpinv+1;
  }
  state[threadId] = mState;
}

#if 0
//naive way of removing excessive indices
  template <typename T>
__global__ void naive(curandState *state, T *data, const T *split, const T *to_remove)
{
  size_t threadId = blockIdx.x * blockDim.x + threadIdx.x;
  T idx = split[threadId];
  T num = split[threadId+1]-split[threadId];
  curandState mState = state[threadId];

  for(unsigned int i = 0; i < to_remove[threadId]; i++)
  {
    //pick a random element and remove it if it wasn't removed before
    T r = curand_uniform(&mState)*num;
    if(data[idx+r] == (T)-1)
      --i;
    data[idx+r] = (T)-1;
  }
  state[threadId] = mState;
}
#endif

//Vitter's Algorithm A (base case of Algorithm D)
  template <typename T>
__global__ void algorithmA(curandState *state, T *data, const T *split, const T *to_remove)
{
  size_t threadId = blockIdx.x * blockDim.x + threadIdx.x;
  T idx = split[threadId];
  T num = split[threadId+1]-split[threadId];
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

#if 0
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
#endif

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
    size_t N;
    const T num_threads = 1<<16;
    const T num_threads_per_group = 1<<3;
    const T num_groups = num_threads / num_threads_per_group;
    T *device_split, *device_to_remove, *split, *to_remove;

    cuda_generator()
    {
      N = 0;
      cudaMalloc((void**)&state, num_threads * sizeof(curandState));
      cudaMalloc((void**)&seeds_device, num_threads * sizeof(unsigned long long));
      seeds = new unsigned long long[num_threads];
      device_data = 0;
      split = new T[num_threads+1];
      to_remove = new T[num_threads];

      //remove excessive elements per processor
      cudaMalloc((void**)&device_split, sizeof(T)*(num_threads+1));
      cudaMalloc((void**)&device_to_remove, sizeof(T)*num_threads);
    }

    ~cuda_generator()
    {
      if(split)
        delete [] split;
      if(to_remove)
        delete [] to_remove;
      cudaFree((void**)&device_split);
      cudaFree((void**)&device_to_remove);
      cudaFree((void**)&state);
      cudaFree((void**)&seeds_device);
    }


  public:
    static cuda_generator<T> *instance()
    {
      if(!_instance)
        _instance = new cuda_generator<T>();
      return _instance;
    }

    static void sampleR(T n, T N, T j, T k, T *startIdx, T *numSamples, StochasticLib1 stoc = StochasticLib1(0), T off = 0)
    {
      if(k == j+1)
      {
        numSamples[j] = n;
        startIdx[j] = N+off;
        return;
      }
      T N2 = N/2;
      T x = stoc.Hypergeometric(n, N2, N);
      sampleR(x, N2, j, (j+k)/2, startIdx, numSamples, stoc, off);
      sampleR(n-x, N-N2, (j+k)/2, k, startIdx, numSamples, stoc, off+N2);
    }

    template <typename It>
      void generate_block(It dest, size_t size, size_t k, size_t universe, double p,
          unsigned long long seed = 0)
      {
        //fprintf(stderr, "%d: universe=%Iu\n", __LINE__, universe);
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
        T samples_per_thread = max((T)ceil(N/(double)num_threads), (T)1);
        T num_indices = 0;
        timer t;
        do {
          sample<<<num_groups, num_threads_per_group>>>(state, device_data, N, samples_per_thread, 1.0-p);
          //prefix-sum for calculating array-indices
          thrust::inclusive_scan(thrust::device, device_data, device_data+N, device_data);
          //count_if is faster than the other two.
          num_indices = thrust::count_if(thrust::device, device_data, device_data+N, less_equal_functor<T>((T)universe));
          //num_indices = thrust::find_if(thrust::device, device_data, device_data+N, less_equal_functor<T>((T)universe)) - device_data;
          //num_indices = thrust::copy_if(thrust::device, device_data, device_data+N, device_data, less_equal_functor<T>((T)universe)) - device_data;
        } while(num_indices < k);
        cudaDeviceSynchronize();

        //divide array using hypergeometric deviates (one subarray per processor)
        sampleR(num_indices-k, num_indices, 0, num_threads, split+1, to_remove);
        split[0] = 0;
        //copy start- and end-indices of subarrays as well as the number
        //of elements to be sampled from these subarrays to GPU memory
        cudaMemcpy(device_split, split, sizeof(T)*(num_threads+1), cudaMemcpyHostToDevice);
        cudaMemcpy(device_to_remove, to_remove, sizeof(T)*num_threads, cudaMemcpyHostToDevice);

        //sample elements (elements are removed by setting the index to some value larger
        //than the universe size
        //naive<<<num_groups, num_threads_per_group>>>(state, device_data, device_split, device_to_remove);
        algorithmA<<<num_groups/4, num_threads_per_group*4>>>(state, device_data, device_split, device_to_remove);

        //cudaDeviceSynchronize();

        //compaction
        thrust::copy_if(thrust::device, device_data, device_data+num_indices, device_data, less_equal_functor<T>(universe));
        //read back data
        //cudaMemcpy(&dest[0], device_data, n*sizeof(T), cudaMemcpyDeviceToHost);
      }
};
template <typename T> cuda_generator<T> *cuda_generator<T>::_instance = 0;

struct cuda_gen {
  template <typename It>
    static void generate_block(It dest, size_t size, size_t k, size_t universe, double p,
        unsigned long long seed = 0)
    {
      using value_type = typename std::iterator_traits<It>::value_type;
      cuda_generator<value_type> *generator = cuda_generator<value_type>::instance();

      generator->generate_block(dest, size, k, universe, p, seed);
    }

  template <typename It>
    static void generate_block(It begin, It end, double p,
        unsigned int long long = 0)
    {
      generate_block(begin, end-begin, p, seed);
    }
};

// copied from include/sampler.h
// Formulas from "Sequential Random Sampling" by Ahrens and Dieter, 1985
static std::pair<double, size_t> calc_params(size_t universe, size_t k /* samples */) {
  double r = sqrt(k);
  double a = sqrt(log(1+k/(2*M_PI)));
  a = a + a*a/(3.0 * r);
  size_t b = k + size_t(4 * a * r);
  double p = (k + a * r) / universe;
  return std::make_pair(p, b);
}

int main(int argc, char **argv)
{
  cudaSetDevice(1);
  size_t num_threads = 1;
  arg_parser args(argc, argv);
  size_t universe = args.get<size_t>("n", 1<<30);
  size_t k = args.get<size_t>("k", 1<<20); // sample size

  size_t iterations = args.get<size_t>("i", (1<<30)/k);
  const bool verbose = args.is_set("v") || args.is_set("vv");
  const bool very_verbose = args.is_set("vv");
  const bool quiet = args.is_set("q");

  double p; size_t ssize;
  std::tie(p, ssize) = calc_params(universe, k);
  #ifdef USE64BIT
  std::vector<unsigned long long> vec(ssize);
  #else
  std::vector<unsigned int> vec(ssize);
  #endif
  int version = 0;
  curandGetVersion(&version);
  fprintf(stderr, "%d", version);
  statistics mt_stats;

  // warmup
  // k = number of indices (output), p = probability, universe = maximum index size
  // ssize = size of index array
  cuda_gen::generate_block(vec.begin(), ssize, k, universe, p);
  //std::cout << "Running warmup (" << k << " samples)" << std::endl;

  std::stringstream extra_stream;
  extra_stream << " k=" << k << " b=" << ssize
    << " p=" << p << " N=" << universe;
  auto extra = extra_stream.str();

  //std::cout << "Running measurements..." << std::endl;

  // Measure
  cudaDeviceSynchronize();

  // stats for multi-threaded version including sync / load imbalance
  timer t;
  for(int i = 0; i < iterations; i++)
  {
    cuda_gen::generate_block(vec.begin(), ssize, k, universe, p);
    cudaDeviceSynchronize();
    double time = t.get_and_reset();
    mt_stats.push(time);
  }
  std::cout << " mt_time=" << mt_stats.avg()
    << " mt_dev=" << mt_stats.stddev()
    << " numthreads=" << num_threads
    << " iterations=" << iterations
    << extra << std::endl;
  #ifdef PROFILING
  cudaProfilerStop();
  #endif
}
