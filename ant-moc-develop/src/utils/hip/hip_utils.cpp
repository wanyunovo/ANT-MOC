#include "antmoc/hip/hip_utils.h"
#include "antmoc/log.h"
#include "antmoc/mpi_utils.h"

#include <thrust/device_vector.h>


namespace antmoc {
namespace hip {

///< GPU id for the current MPI rank.
static int _gpu_id = -1;


void hipAssertImpl(hipError_t err, const char* expression, const char *file, int line) {
   if (err != hipSuccess) {
      log::warn("File {}, line {}: {}, failed: {}", file, line, expression,
                hipGetErrorString(err));
      abort();
   }
}


int getGPUId() {

  // FIXME: remove this branch
  // Compute once
  if (_gpu_id < 0) {

      // Number of GPUs
      auto n_gpus = getNumGPUs();

      if (n_gpus < 0) {
        log::error("There is no GPU on this machine");
      }

      _gpu_id = mpi::getRankNodeLocal();

      if (_gpu_id < 0 || _gpu_id >= n_gpus) {
          log::error("Rank {} was assigned an invalid GPU {}, number of GPUs is {}",
                      mpi::getMPIRank(), _gpu_id, n_gpus);
      }
  }

  return _gpu_id;
}


int getNumGPUs() {
  int n_gpus;
  HIP_ASSERT( hipGetDeviceCount(&n_gpus) );

  return n_gpus;
}


void warmUp() {
  thrust::device_vector<float> X(1);
}


void setDevice(int id) {

  HIP_ASSERT( hipDeviceReset() );
  HIP_ASSERT( hipSetDevice(id) );

  log::verbose("Successfully set device id {} for rank {}", id, mpi::getMPIRank());
}


void synchronizeDevice() {

  HIP_ASSERT( hipDeviceSynchronize() );

  log::debug("Device {} synchronized", getGPUId());
}


// Implementation
template <typename T>
void arrayCopyHtoD(T *&ptr, T *ptr_host, size_t length) {
  void *array_d;
  auto size = length * sizeof(T);

  // Allocate memory on device
  HIP_ASSERT( hipMalloc(&array_d, size) );
  // Copy memory from host to device
  HIP_ASSERT( hipMemcpyHtoD(array_d, ptr_host, size) );
  // Copy host pointer to device
  HIP_ASSERT( hipMemcpyHtoD(&ptr, &array_d, sizeof(T*)) );
}

// Instantiation
template void arrayCopyHtoD(int*&, int*, size_t);
template void arrayCopyHtoD(long*&, long*, size_t);
template void arrayCopyHtoD(float*&, float*, size_t);
template void arrayCopyHtoD(double*&, double*, size_t);
template void arrayCopyHtoD(size_t*&, size_t*, size_t);


// Implementation
template <typename T>
void vectorCopyHtoD(T *&ptr, const std::vector<T> &vec_host) {
  arrayCopyHtoD(ptr, const_cast<T*>(vec_host.data()), vec_host.size());
}

// Instantiation
template void vectorCopyHtoD(int*&, const std::vector<int>&);
template void vectorCopyHtoD(long*&, const std::vector<long>&);
template void vectorCopyHtoD(float*&, const std::vector<float>&);
template void vectorCopyHtoD(double*&, const std::vector<double>&);
template void vectorCopyHtoD(size_t*&, const std::vector<size_t>&);


// Implementation
template <typename T>
void matrixCopyHtoD(T *&ptr, const std::vector<std::vector<T>> &mat_host) {
  // Check the shape of this matrix
  auto length1 = mat_host.size();

  if (length1 == 0) {
    log::error("Could not copy empty matrix (2D vector) from host to device");
  }

  auto length2 = mat_host[0].size();

  for (const auto &v : mat_host) {
    if (v.size() == 0) {
      log::error("Could not copy matrix with 0-length element from host to device");
    }
    else if (v.size() != length2) {
      log::error("Only n-by-n matrix could be copied to device by matrixCopyHtoD");
    }
  }

  // Flatten the matrix
  std::vector<T> vec_flatten;
  vec_flatten.reserve(length1 * length2);
  for (const auto &v : mat_host) {
    vec_flatten.insert(vec_flatten.end(), v.begin(), v.end());
  }

  // Copy memory from host to device
  vectorCopyHtoD(ptr, vec_flatten);
}

// Instantiation
template void matrixCopyHtoD(int*&, const std::vector<std::vector<int>>&);
template void matrixCopyHtoD(long*&, const std::vector<std::vector<long>>&);
template void matrixCopyHtoD(float*&, const std::vector<std::vector<float>>&);
template void matrixCopyHtoD(double*&, const std::vector<std::vector<double>>&);
template void matrixCopyHtoD(size_t*&, const std::vector<std::vector<size_t>>&);


} // namespace hip
} // namespace antmoc
