/*******************************************************************************
 * cuda_gen.cu
 *
 * Copyright (C) 2016 Emanuel Schrade <schrade@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
******************************************************************************/

#include <cassert>
#include <limits>
#include <random>
#include <cuda.h>
#include <curand_kernel.h>
#include <curand.h>
//#define PROFILING
#ifdef PROFILING
#include <cuda_profiler_api.h>
#endif
#include "util.h"
#include "timer.h"
#include <stocc/stocc.h>
#include <stocc/randomc.h>

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

#ifdef USE_ALGORITHM_H
#include <set>
#else
#define _DEFINITIONS_H_
#include <sampler.h>
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
  float lpinv = 1.0/log2(oneminusp);
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

//set elements to MAX to mark for removal
template <typename T>
__global__ void markRemove(T *data, const size_t count, T *remove, const size_t remove_count)
{
  size_t threadId = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadId < remove_count)
    data[remove[threadId]] = T(-1);
}

//XXX Not working in all cases
#if 0
//replace k elements with k last ones
template <typename T>
__global__ void replace(T *data, const size_t count, T *remove, const size_t remove_count)
{
  size_t threadId = blockIdx.x * blockDim.x + threadIdx.x;
  T idx = (T)-1;
  if(threadId < remove_count)
    idx = remove[threadId];

  if(idx != (T)-1 && data[idx] == (T)-1 && data[count-threadId-1] != (T)-1)
  {
    data[idx] = data[count-threadId-1];
    remove[threadId] = (T)-1;
  }
  if(idx == count-threadId-1)
    remove[threadId] = (T)-1;
  //otherwise replace in next iteration
}
#endif

//functor for filtering elements less than a constant
template <typename T>
struct less_functor
{
  T max_idx;
  less_functor(const T _max_idx)
    : max_idx(_max_idx)
  { }

  __host__ __device__
    bool operator()(T x)
    {
      return x < max_idx;
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
    size_t N;
    const T num_threads = 1<<16;
    const T num_threads_per_group = 1<<3;
    const T num_groups = num_threads / num_threads_per_group;

    cuda_generator()
    {
      N = 0;
      cudaMalloc((void**)&state, num_threads * sizeof(curandState));
      cudaMalloc((void**)&seeds_device, num_threads * sizeof(unsigned long long));
      seeds = new unsigned long long[num_threads];
      device_data = 0;
    }

    ~cuda_generator()
    {
      cudaFree((void**)&state);
      cudaFree((void**)&seeds_device);
      delete [] seeds;
    }


  public:
    static cuda_generator<T> *instance()
    {
      if(!_instance)
        _instance = new cuda_generator<T>();
      return _instance;
    }

    template <typename It>
    void generate_block(It dest, size_t size, size_t k, size_t universe, double p,
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
      T samples_per_thread = max((T)ceil(N/(double)num_threads), (T)1);
      T num_indices = 0;
      timer t;
      do {
        sample<<<num_groups, num_threads_per_group>>>(state, device_data, N, samples_per_thread, 1.0-p);
        //prefix-sum for calculating array-indices
        thrust::inclusive_scan(thrust::device, device_data, device_data+N, device_data);
        //count_if is faster than the other two.
        num_indices = thrust::count_if(thrust::device, device_data, device_data+N, less_functor<T>((T)universe));
        //num_indices = thrust::find_if(thrust::device, device_data, device_data+N, less_functor<T>((T)universe)) - device_data;
        //num_indices = thrust::copy_if(thrust::device, device_data, device_data+N, device_data, less_functor<T>((T)universe)) - device_data;
      } while(num_indices < k);

      //sample num_indices - k elements to remove
      T *removeArray = new T[num_indices - k];
      #ifdef USE_ALGORITHM_H
      std::mt19937 gen(seed);
      std::uniform_int_distribution<T> dist(0, num_indices-1);
      std::set<T> remove;

      while(remove.size() < num_indices - k)
        remove.insert(dist(gen));

      int idx = 0;
      for(auto it = remove.begin(); it != remove.end(); it++)
        removeArray[idx++] = *it;
      #else
	  size_t hole_idx = 0;
	  const size_t basecase = 1024;
	  size_t to_remove = num_indices - k;
	  HashSampling<> hs((ULONG)seed, to_remove);
	  SeqDivideSampling<> s(hs, basecase, (ULONG)seed);
	  // end - begin - 1 because the range is inclusive
	  s.sample(num_indices - 1, to_remove, [&](ULONG pos) {
		  removeArray[hole_idx++] = pos;
	  });
	  assert(hole_idx == to_remove);
      #endif
      //copy indices to GPU
      T *device_removeArray;
      cudaMalloc((void**)&device_removeArray, sizeof(T) * (num_indices-k));
      cudaMemcpyAsync(device_removeArray, removeArray, sizeof(T) * (num_indices-k), cudaMemcpyHostToDevice);

      markRemove<<<(num_indices-k)/num_threads_per_group*4, num_threads_per_group*4>>>(device_data, num_indices, device_removeArray, num_indices-k);
      #if 1
      //compaction
      int removed = thrust::copy_if(thrust::device, device_data, device_data+num_indices, device_data, less_functor<T>(universe)) - device_data;
      #else
      //copy elements and count collisions, restart until no collisions left
      remove_remaining = num_indices-k;
      while(remove_remaining > 0)
      {
        replace<<<remove_remaining, num_threads_per_group>>>(device_data, k+remove_remaining, device_removeArray, num_indices-k);
        remove_remaining = thrust::copy_if(thrust::device, device_removeArray, device_removeArray + remove_remaining, device_removeArray, less_functor<T>((T)universe))-device_removeArray;
      }

      cudaFree(&device_removeArray);
      #endif
      cudaFree(&device_removeArray);
      delete [] removeArray;
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
  cudaSetDevice(3);
  size_t num_threads = 1;
  arg_parser args(argc, argv);
  #ifdef USE64BIT
  size_t universe = args.get<size_t>("n", 1ULL << 50);
  #else
  size_t universe = args.get<size_t>("n", 1<<30);
  #endif
  size_t k = args.get<size_t>("k", 1<<28); // sample size

  size_t iterations = args.get<size_t>("i", 3*(1ULL<<30)/k);
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
  statistics mt_stats;

  // warmup
  // k = number of indices (output), p = probability, universe = maximum index size
  // ssize = size of index array
  cuda_gen::generate_block(vec.begin(), ssize, k, universe, p);

  std::stringstream extra_stream;
  extra_stream << " k=" << k << " b=" << ssize
    << " p=" << p << " N=" << universe;
  auto extra = extra_stream.str();

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
