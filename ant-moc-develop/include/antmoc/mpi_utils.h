/// \file mpi_utils.h
/// \brief MPI utilities.
/// \date March 21, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)


#ifndef MPI_UTILS_H_
#define MPI_UTILS_H_


#ifdef ENABLE_MPI_
#include <mpi.h>
#endif

namespace antmoc {

namespace mpi {

//----------------------------------------------------------------------
// Common operations
//----------------------------------------------------------------------
/// \brief Checks whether I am the root in _comm_world_cart.
/// \return Boolean indicating I am the root (true) or not (false).
bool isMPIRoot();

/// \brief Checks whether I am the root in _comm_unique_domains.
/// \return Boolean indicating I am the root (true) or not (false).
bool isRootUniqueDomains();

/// \brief Checks whether I am the root in _comm_shared_domain.
/// \return Boolean indicating I am the root (true) or not (false).
bool isRootSharedDomain();

/// \brief Returns my rank in _comm_world_cart.
/// \return My rank.
int getMPIRank();

/// \brief Returns my rank in _comm_unique_domains.
/// \return My rank.
int getRankUniqueDomains();

/// \brief Returns my rank in _comm_shared_domain.
/// \return My rank.
int getRankSharedDomain();

/// \brief Returns my rank either in _comm_shared_domain or _comm_unique_domains.
int getRankAuto();

/// \brief  Returns my node-local rank.
/// \return Node-local rank.
int getRankNodeLocal();

/// \brief Returns the number of ranks in _comm_world_cart.
/// \return The number of ranks.
int getNumProcs();

/// \brief Returns the number of unique spatial domains.
/// \details Only 1 unique domain exists if the geometry is not spatial
///          decomposed. If both spatial and periodic track decomposition
///          takes effect, the number of periodic track communicators
///          will be returned.
/// \return The number of ranks.
int getNumUniqueDomains();

/// \brief Returns the number of unique periodic track domains sharing the
///        same spatial domain.
/// \return The number of ranks.
int getNumProcsSharedDomain();

/// \brief Returns the number of processes either in _comm_shared_domain
///        or _comm_unique_domains
int getNumProcsAuto();

//----------------------------------------------------------------------
// Decomposition status
//----------------------------------------------------------------------
/// \brief Checks whether the MPI variables are initialized.
/// \return Boolean indicating initialization is completed (true) or not (false).
bool isMPIInitialized();

/// \brief Checks whether any decomposition is perfomed.
/// \return Boolean indicating any decomposition is done (true) or not (false).
bool isDomainDecomposed();

/// \brief Checks whether the space is decomposed.
/// \return Boolean indicating the domain is decomposed (true) or not (false).
bool isSpatialDecomposed();

/// \brief Checks whether periodic tracks are decomposed.
/// \return Boolean indicating tracks are decomposed (true) or not (false).
bool isPrdTrackDecomposed();

//----------------------------------------------------------------------
// Communicator manipulation
//----------------------------------------------------------------------
#ifdef ENABLE_MPI_
/// \brief Initializes variables to default values (MPI only).
void defaultInitialize();

/// \brief Returns the communicator consisting of all processes (MPI only).
/// \return The communicator ID.
MPI_Comm getMPIComm();

/// \brief Returns the communicator in which processes have unique spatial domains
///        (MPI only).
/// \return The communicator ID.
MPI_Comm getCommUniqueDomains();

/// \brief Returns the communicator in which processes have the same spatial domain
///        (MPI only).
/// \return The communicator ID.
MPI_Comm getCommSharedDomain();

/// \brief Returns either _comm_shared_domain or _comm_unique_domains.
MPI_Comm getCommAuto();

/// \brief Sets a name for the main spatial communicator which has the root in it
///        (MPI only).
/// \details This communicator is a 'slice' of _comm_world_cart while it has
///          the root. It is useful for manipulating files.
void setMainCommUniqueDomains();

/// \brief Initializes decomposition for spatial domains and track domains (MPI only).
/// \param comm The MPI communicator to be used to communicate between domains.
/// \param nx The number of MPI domains in the x-direction (defaults to 1).
/// \param ny The number of MPI domains in the y-direction (defaults to 1).
/// \param nz The number of MPI domains in the z-direction (defaults to 1).
/// \param ns The number of MPI ranks sharing the same domain (defaults to 1).
void setDomainDecomposition(MPI_Comm comm = MPI_COMM_WORLD,
                            int nx = 1, int ny = 1, int nz = 1, int ns = 1);
#endif

/// \brief Checks whether I am in the main communicator.
/// \return Boolean indicating I am in the main communicator (true) or not (false).
bool isInMainCommUniqueDomains();

//----------------------------------------------------------------------
// Domain decomposition operations
//----------------------------------------------------------------------
#ifdef ENABLE_MPI_
/// \brief Returns the rank of the specified neighboring domain (MPI only).
/// \param offset_x The shift in the x-direction (number of domains).
/// \param offset_y The shift in the y-direction (number of domains).
/// \param offset_z The shift in the z-direction (number of domains).
/// \return The rank of the neighboring domain.
int getNeighborDomain(int offset_x, int offset_y, int offset_z);
#endif

/// \brief Checks whether this MPI domain is the root domain.
/// \return Boolean indicating this domain is the root domain (true) or not (false).
bool isRootDomain();

/// \brief Returns the number of domains in the X direction.
int getNumDomainsX();

/// \brief Returns the number of domains in the Y direction.
int getNumDomainsY();

/// \brief Returns the number of domains in the Z direction.
int getNumDomainsZ();

/// \brief Returns the index of domain in the X direction.
int getDomainIndexX();

/// \brief Returns the index of domain in the Y direction
int getDomainIndexY();

/// \brief Returns the index of domain in the Z direction
int getDomainIndexZ();

/// \brief Returns the UID of a domain in the rectilinear grid
int getDomainUid();

/// \brief Returns the domain indexes of the current domain.
/// \param indexes Pointer to the array to be filled with the indexes.
void getDomainIndexes(int* indexes);

/// \brief Returns the number of domains in each direction.
/// \param structure Pointer to the array to be filled with domain numbers.
void getDomainStructure(int* structure);

//----------------------------------------------------------------------
// MPI utility
//----------------------------------------------------------------------
#ifdef ENABLE_MPI_
/// \brief Returns the MPI datatype (MPI only).
template <typename T> MPI_Datatype getDatatype();
extern template MPI_Datatype getDatatype<int>();
extern template MPI_Datatype getDatatype<long>();
extern template MPI_Datatype getDatatype<float>();
extern template MPI_Datatype getDatatype<double>();

/// \brief Shows a report for MPI datatypes (MPI only).
void showMPIDatatypes();

/// \brief Synchronizes across all ranks (MPI only).
void mpiBarrier();

/// \brief Synchronizes across ranks that have unique domains (MPI only).
void mpiBarrierUniqueDomains();

/// \brief Synchronizes across ranks that have the same domain (MPI only).
void mpiBarrierSharedDomain();

/// \brief Calls MPI_Allreduce across all processes (MPI only).
/// \param sendbuf Starting address of send buffer.
/// \param recvbuf Starting address of receive buffer (output)
/// \param count Number of elements in send buffer.
/// \param datatype Data type of elements of send buffer.
/// \param op Operation.
/// \return Error value.
int mpiAllreduce(const void *sendbuf, void *recvbuf, int count,
                 MPI_Datatype datatype, MPI_Op op);

/// \brief Calls MPI_Allreduce across ranks that have unique domains (MPI only).
int mpiAllreduceUD(const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op);

/// \brief Calls MPI_Allreduce across ranks that have the same domain (MPI only).
int mpiAllreduceSD(const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op);
#endif

} // namespace mpi

} // namespace antmoc

#endif // MPI_UTILS_H_
