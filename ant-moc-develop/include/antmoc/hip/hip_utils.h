/// \file hip_utils.h
/// \brief HIP utilities.
/// \date December 19, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef HIP_UTILS_H_
#define HIP_UTILS_H_

#ifdef ENABLE_HIP_

#include "hip/hip_runtime.h"
#include <vector>

/// \brief Macro for checking HIP error.
#define HIP_ASSERT(expression) \
   { ::antmoc::hip::hipAssertImpl((expression), #expression, __FILE__, __LINE__); }

namespace antmoc {
namespace hip {

/// \brief Check HIP error.
void hipAssertImpl(hipError_t, const char*, const char *, int);

/// \brief  Get the ID of the GPU assigned to current rank.
/// \return GPU ID (also the node-local rank)
int getGPUId();

/// \brief  Get the number of available GPUs on this node.
/// \return Number of GPUs
int getNumGPUs();

/// \brief Warm up the GPU.
void warmUp();

/// \brief Set the device id for the current MPI rank.
void setDevice(int id);

/// \brief Synchronize device.
void synchronizeDevice();

/// \brief Copy contiguous memory from host to device.
/// \param ptr      Pointer to device memory (device pointer, output).
/// \param ptr_host Pointer to host memory.
/// \param length   Number of array items.
template <typename T>
extern void arrayCopyHtoD(T *&ptr, T *ptr_host, size_t length);

/// \brief Copy contiguous memory from host to device.
/// \param ptr      Pointer to device memory (device pointer, output).
/// \param vec_host Vector on host.
template <typename T>
extern void vectorCopyHtoD(T *&ptr, const std::vector<T> &vec_host);

/// \brief Copy an n-by-n matrix from host to contiguous device memory.
/// \param ptr      Pointer to device memory (device pointer, output).
/// \param mat_host N-by-n matrix on host.
template <typename T>
extern void matrixCopyHtoD(T *&ptr, const std::vector<std::vector<T>> &mat_host);

/// \brief   Print a device object.
/// \details The object must have implemented printString().
template <typename T>
__global__
extern void printObject(T *ptr);


} // namespace hip
} // namespace antmoc

#endif  // ENABLE_HIP_
#endif  // HIP_UTILS_H_
