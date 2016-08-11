#pragma once

#include <cassert>
#include <limits>
#include <random>
#include <curand.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>

//Sample number with geometric distribution:
template <typename T>
struct geometric_distribution : public thrust::unary_function<T,T>
{
  const double p;
  const double maxvalinv;

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
    thrust::device_vector<T> device_data;
    size_t N;
    curandGenerator_t prng;

    cuda_generator()
    {
      N = 0;
      //curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);
      curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_MTGP32);
      //curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_MT19937);
    }

    void call_curand_generate()
    {
      curandGenerate(prng, thrust::raw_pointer_cast(device_data.data()), N);
    }

  public:
  static cuda_generator<T> *instance()
  {
    if(!_instance)
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
      device_data = thrust::device_vector<T>(size);
      N = size;
    }

    //set seed of prng
    curandSetPseudoRandomGeneratorSeed(prng, seed);

    //fill array with random numbers (should depend on size of the data type)
    call_curand_generate();

    fprintf(stderr, ".");

    //geometric distribution and prefix sum
    auto tbegin = thrust::make_transform_iterator(device_data.begin(), geometric_distribution<T>(p, 1.0/(std::numeric_limits<T>::max())));
    auto tend = thrust::make_transform_iterator(device_data.end(), geometric_distribution<T>(p, 1.0/(std::numeric_limits<T>::max())));
    thrust::inclusive_scan(tbegin, tend, device_data.begin());

    //read back data
    cudaMemcpy(&dest[0], device_data.data().get(), N*sizeof(T), cudaMemcpyDeviceToHost);
  }
};
template <typename T> cuda_generator<T> *cuda_generator<T>::_instance = 0;
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
