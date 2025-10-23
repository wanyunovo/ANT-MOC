#include "antmoc/pairwise_sum.h"

#include <cmath>

namespace antmoc {

template <typename T>
double pairwise_sum(T* vector, size_t length) {

  double sum = 0;

  /* Base case: if length is less than 16, perform summation */
  if (length < 16) {

    #pragma omp simd reduction(+:sum)
    for (size_t i=0; i < length; i++)
      sum += vector[i];
  }

  else {  //将数组分成两半，分别计算每一半的和，然后将两部分的结果相加,为什么不直接加呢 因为递归分治策略，减少了计算过程中浮点数加法的误差
    size_t offset = length % 2;
    length = std::floor(length / 2);
    sum = pairwise_sum<T>(&vector[0], length) +
          pairwise_sum<T>(&vector[length], length+offset);
  }

  return sum;
}

// Instantiation
template double pairwise_sum(float*, size_t);
template double pairwise_sum(double*, size_t);

} // namespace antmoc
