/// \file hip_info.h
/// \brief Routines for getting some information of HIP.
/// \date December 19, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef HIP_INFO_H_
#define HIP_INFO_H_

namespace antmoc {
namespace hip {

/// \brief Print device properties
void printDeviceProperties();

/// \brief  Check if the machine contains a GPU.
/// \return True if there is at least one GPU, false otherwise.
bool containsGPU();

/// \brief  Get the number of compute nodes
/// \return Number of CUs
int getNumCUs();

}   // namespace hip
}   // namespace antmoc

#endif  // HIP_INFO_H_
