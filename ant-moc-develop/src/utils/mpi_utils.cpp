#include "antmoc/mpi_utils.h"
#include "antmoc/log.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

namespace antmoc
{

  namespace mpi
  {

    //----------------------------------------------------------------------
    // Static variables
    //----------------------------------------------------------------------

    ///< Rank of this process in _comm_world_cart
    static int _rank_world = 0;
    ///< Rank of this process in _comm_unique_domains
    static int _rank_unique_domains = 0;
    ///< Rank of this process in _comm_shared_domain
    static int _rank_shared_domain = 0;
    ///< Rank of this process local to a compute node
    static int _rank_node_local = -1;

    ///< The number of processes in _comm_world_cart
    static int _np_world = 1;
    ///< The number of processes in _comm_unique_domains
    static int _np_unique_domains = 1;
    ///< The number of processes in _comm_shared_domain
    static int _np_shared_domain = 1;

    ///< Number of domains in the X, Y and Z directions
    /**
     * MPI初始化时设置域数，在 mpi_utils.cpp 的 setDomainDecomposition 函数中（第224行附近）
     */
    static int _num_domains_x = 1;
    static int _num_domains_y = 1;
    static int _num_domains_z = 1;

    ///Index of the domain in the whole geometry in the X, Y and Z directions
    static int _domain_index_x = 1;
    static int _domain_index_y = 1;
    static int _domain_index_z = 1;

    //----------------------------------------------------------------------
    // Common operations
    //----------------------------------------------------------------------

    bool isMPIRoot()
    {
      return (_rank_world == 0) ? true : false;
    }

    bool isRootUniqueDomains()
    {
      return (_rank_unique_domains == 0 ? true : false);
    }

    bool isRootSharedDomain()
    {
      return (_rank_shared_domain == 0 ? true : false);
    }

    int getMPIRank()
    {
      return _rank_world;
    }

    int getRankUniqueDomains()
    {
      return _rank_unique_domains;
    }

    int getRankSharedDomain()
    {
      return _rank_shared_domain;
    }

    int getRankAuto()
    {
#ifdef USTB_
      return _rank_shared_domain;
#else
      return _rank_unique_domains;
#endif
    }

    int getNumProcs()
    {
      return _np_world;
    }

    int getNumUniqueDomains()
    {
      return _np_unique_domains;
    }

    int getNumProcsSharedDomain()
    {
      return _np_shared_domain;
    }

    int getNumProcsAuto()
    {
#ifdef USTB_
      return _np_shared_domain;
#else
      return _np_unique_domains;
#endif
    }

    int getRankNodeLocal()
    {

      // Compute only once
      if (_rank_node_local < 0)
      {

        _rank_node_local = 0;

#ifdef ENABLE_MPI_
        // MPI communicator for node-local ranks
        MPI_Comm comm_node;

        // Creates new communicators based on split types and keys
        MPI_Comm_split_type(
            mpi::getMPIComm(),    // MPI_Comm communicator
            MPI_COMM_TYPE_SHARED, // int      split_type
            0,                    // int      key
            MPI_INFO_NULL,        // MPI_Info info
            &comm_node            // MPI_Comm *new_communicator
        );

        MPI_Comm_rank(comm_node, &_rank_node_local);
#endif
      }

      return _rank_node_local;
    }

    //----------------------------------------------------------------------
    // Decomposition status
    //----------------------------------------------------------------------

    ///< Indicate whether MPI variables are initialized
    static bool _mpi_initialized = false;
    ///< Indicate whether the Geometry is domain decomposed
    static bool _spatial_decomposed = false;
    ///< Indicate whether the periodic tracks are decomposed
    static bool _prd_track_decomposed = false;

    bool isMPIInitialized()
    {
      return _mpi_initialized;
    }

    bool isDomainDecomposed()
    {
      return (_spatial_decomposed || _prd_track_decomposed);
    }

    bool isSpatialDecomposed()
    {
      return _spatial_decomposed;
    }

    bool isPrdTrackDecomposed()
    {
      return _prd_track_decomposed;
    }

//----------------------------------------------------------------------
// Communicator manipulation
//----------------------------------------------------------------------
#ifdef ENABLE_MPI_

    ///< The topology including all of the processes
    static MPI_Comm _comm_world_cart = MPI_COMM_NULL;

    ///< The communicator in which processes have unique spatial domains
    static MPI_Comm _comm_unique_domains = MPI_COMM_NULL;

    ///< The communicator in which processes have the same spatial domain
    static MPI_Comm _comm_shared_domain = MPI_COMM_NULL;

    void defaultInitialize()
    {
      _comm_world_cart = MPI_COMM_WORLD;
      MPI_Comm_rank(_comm_world_cart, &_rank_world);
      MPI_Comm_size(_comm_world_cart, &_np_world);
    }

    MPI_Comm getMPIComm()
    {
      if (_comm_world_cart == MPI_COMM_NULL)
      {
        log::warn("Uninitialized communicator returned at getMPIComm()");
      }
      return _comm_world_cart;
    }

    MPI_Comm getCommUniqueDomains()
    {
      if (_comm_unique_domains == MPI_COMM_NULL)
      {
        log::warn("Uninitialized communicator returned at getCommUniqueDomains()");
      }
      return _comm_unique_domains;
    }

    MPI_Comm getCommSharedDomain()
    {
      if (_comm_shared_domain == MPI_COMM_NULL)
      {
        log::warn("Uninitialized communicator returned at getCommSharedDomain()");
      }
      return _comm_shared_domain;
    }

    MPI_Comm getCommAuto()
    {
#ifdef USTB_
      return getCommSharedDomain();
#else
      return getCommUniqueDomains();
#endif
    }

    void setMainCommUniqueDomains()
    {
      // Compute the index of _comm_unique_domains
      int color = _rank_world % _np_shared_domain;
      std::string comm_name;
      if (color == 0)
      {
        comm_name = "main communicator";
      }
      else
      {
        comm_name = "spatial communicator " + std::to_string(color);
      }

      MPI_Comm_set_name(_comm_unique_domains, comm_name.c_str());
    }

    void setDomainDecomposition(MPI_Comm comm, int nx, int ny, int nz, int ns)
    {
      if (_mpi_initialized)
        return;

      if (nx < 1 || ny < 1 || nz < 1 || ns < 1)
      {
        log::error("Invalid number of domains: {},{},{},{}", nx, ny, nz, ns);
      }

      // Calculate number of domains and get the number of MPI ranks
      MPI_Comm_size(comm, &_np_world);
      _np_unique_domains = nx * ny * nz;
      _np_shared_domain = ns;

      log::debug("Number of processes which have unique domains: {}", _np_unique_domains);
      log::debug("Number of processes which shares the same sub-domain: {}", _np_shared_domain);

      // Check that the number of domains equals the number of ranks
      if (_np_world != _np_unique_domains * _np_shared_domain)
      {
        log::error("Issued ranks ({}) is not equivalent to the required ({}x{}x{}x{})",
                   _np_world, nx, ny, nz, ns);
      }
      else
      {
        log::debug("Issued ranks = {}, required ranks = {}x{}x{}x{}",
                   _np_world, nx, ny, nz, ns);
      }

      _num_domains_x = nx;
      _num_domains_y = ny;
      _num_domains_z = nz;

      // Debug info
      log::debug("Set the number of spatial domains: {}, {}, {}",
                 _num_domains_x, _num_domains_y, _num_domains_z);

      // Create the MPI Communicator for all processes
      int dims[4] = {nx, ny, nz, ns};
      int wrap[4] = {false, false, false, false};
      MPI_Cart_create(comm, 4, dims, wrap, true, &_comm_world_cart);

      // Create sub communicators with nx-ny-nz-1 layout
      // This operation generates several communicators each of which is
      // a shared-spatial-domain communicator.
      int remain_spatial[4] = {true, true, true, false};
      MPI_Cart_sub(_comm_world_cart, remain_spatial, &_comm_unique_domains);

      // Create sub communicators with 1-1-1-_np_shared_domain layout
      int remain_track[4] = {false, false, false, true};
      MPI_Cart_sub(_comm_world_cart, remain_track, &_comm_shared_domain);

      // Determine my ranks
      MPI_Comm_rank(_comm_shared_domain, &_rank_shared_domain);
      MPI_Comm_rank(_comm_world_cart, &_rank_world);
      MPI_Comm_rank(_comm_unique_domains, &_rank_unique_domains);

      // Name the main spatial communicator which has the root in it
      setMainCommUniqueDomains();

      // Determine the spatial domain indexes
      int cart_coords[3];
      MPI_Cart_coords(_comm_unique_domains, _rank_unique_domains, 3, cart_coords);
      _domain_index_x = cart_coords[0];
      _domain_index_y = cart_coords[1];
      _domain_index_z = cart_coords[2];

      // Make note of the domain decomposition
      if (_np_unique_domains > 1)
      {
        _spatial_decomposed = true;

        // Case 1: spatial + periodic track decomposition
        if (_np_shared_domain > 1)
        {
          _prd_track_decomposed = true;
          // Send information to the user
          log::info("Periodic track domain decomposition = enabled");
        }
        // Case 2: purely spatial decomposed
        else
        {
          _prd_track_decomposed = false;
          log::info("Periodic track domain decomposition = disabled");
        }
        // Send information to the user
        log::info("Spatial domain decomposition        = enabled");
      }
      // Case 3: purely periodic track decomposed
      else if (_np_shared_domain > 1)
      {
        _spatial_decomposed = false;
        _prd_track_decomposed = true;

        // Determine my ranks
        MPI_Comm_rank(_comm_shared_domain, &_rank_shared_domain);
        _rank_world = _rank_shared_domain;

        log::info("Periodic track domain decomposition = enabled");
        log::info("Spatial domain decomposition        = disabled");
      }
      // Case 4: no decomposition
      else
      {
        _spatial_decomposed = false;
        _prd_track_decomposed = false;

        log::info("Periodic track domain decomposition = disabled");
        log::info("Spatial domain decomposition        = disabled");
      }

      _mpi_initialized = true;

      // Assertion
      if (_rank_unique_domains != getNeighborDomain(0, 0, 0))
      {
        log::error("rank_unique_domains {} != getNeighborDomain(0,0,0) {}",
                   _rank_unique_domains, getNeighborDomain(0, 0, 0));
      }

      // Debugging
      log::debug("Rank(W,U,S) = {},{},{}, sub-domain indx = {},{},{}",
                 _rank_world, _rank_unique_domains, _rank_shared_domain,
                 _domain_index_x, _domain_index_y, _domain_index_z);

      log::debug("Domain {} indx = {},{},{}, neighbor(x+1,y+1,z+1) = {},{},{}",
                 getNeighborDomain(0, 0, 0),
                 _domain_index_x, _domain_index_y, _domain_index_z,
                 getNeighborDomain(1, 0, 0),
                 getNeighborDomain(0, 1, 0),
                 getNeighborDomain(0, 0, 1));
    }

#endif // ENABLE_MPI_

    bool isInMainCommUniqueDomains()
    {
#ifdef ENABLE_MPI_
      if (isSpatialDecomposed())
      {
        char *comm_name = new char[MPI_MAX_OBJECT_NAME];
        int len_name;
        MPI_Comm_get_name(_comm_unique_domains, comm_name, &len_name);

        std::string name(comm_name, len_name);

        log::profile("Name of the main spatial communicator: {}", name);

        return name == "main communicator";
      }
      else
      {
        return isRootSharedDomain();
      }
#else
      return true;
#endif
    }

    //----------------------------------------------------------------------
    // Domain decomposition operations
    //----------------------------------------------------------------------

#ifdef ENABLE_MPI_
    int getNeighborDomain(int offset_x, int offset_y, int offset_z)
    {
      int neighbor_rank = -1;
      int neighbor_coords[3];
      neighbor_coords[0] = offset_x + _domain_index_x;
      neighbor_coords[1] = offset_y + _domain_index_y;
      neighbor_coords[2] = offset_z + _domain_index_z;

      if (neighbor_coords[0] >= 0 && neighbor_coords[0] < _num_domains_x &&
          neighbor_coords[1] >= 0 && neighbor_coords[1] < _num_domains_y &&
          neighbor_coords[2] >= 0 && neighbor_coords[2] < _num_domains_z)
        MPI_Cart_rank(_comm_unique_domains, neighbor_coords, &neighbor_rank);

      return neighbor_rank;
    }
#endif

    bool isRootDomain()
    {
      bool first_domain = true;
      if (_spatial_decomposed)
        if (_domain_index_x != 0 || _domain_index_y != 0 || _domain_index_z != 0)
          first_domain = false;
      return first_domain;
    }

    int getNumDomainsX() { return _num_domains_x; }
    int getNumDomainsY() { return _num_domains_y; }
    int getNumDomainsZ() { return _num_domains_z; }
    int getDomainIndexX() { return _domain_index_x; }
    int getDomainIndexY() { return _domain_index_y; }
    int getDomainIndexZ() { return _domain_index_z; }

    /// \details FIXME: this is not for grids with hexagons
    int getDomainUid()
    {
      int uid = _domain_index_x +
                _domain_index_y * _num_domains_x +
                _domain_index_z * _num_domains_x * _num_domains_y;
      return uid;
    }

    void getDomainIndexes(int *indexes)
    {
      indexes[0] = _domain_index_x;
      indexes[1] = _domain_index_y;
      indexes[2] = _domain_index_z;
    }

    void getDomainStructure(int *structure)
    {
      structure[0] = _num_domains_x;
      structure[1] = _num_domains_y;
      structure[2] = _num_domains_z;
    }

//----------------------------------------------------------------------
// MPI utility
//----------------------------------------------------------------------
#ifdef ENABLE_MPI_
    template <typename T>
    MPI_Datatype getDatatype()
    {

      if (std::is_same<T, int>::value)
      {
        return MPI_INT32_T;
      }
      else if (std::is_same<T, long>::value)
      {
        return MPI_INT64_T;
      }
      else if (std::is_same<T, float>::value)
      {
        return MPI_FLOAT;
      }
      else if (std::is_same<T, double>::value)
      {
        return MPI_DOUBLE;
      }
      else
      {
        log::error("Unsupported type mapping from {} to MPI datatype", typeid(T).name());
        return MPI_BYTE;
      }
    }

    template MPI_Datatype getDatatype<int>();
    template MPI_Datatype getDatatype<long>();
    template MPI_Datatype getDatatype<float>();
    template MPI_Datatype getDatatype<double>();

    void showMPIDatatypes()
    {

      int size = 0;

      MPI_Type_size(MPI_INT, &size);
      log::verbose_once("sizeof(MPI_INT) = {}, sizeof(int) = {}", size, sizeof(int));

      MPI_Type_size(MPI_LONG, &size);
      log::verbose_once("sizeof(MPI_LONG) = {}, sizeof(long) = {}", size, sizeof(long));

      MPI_Type_size(MPI_FLOAT, &size);
      log::verbose_once("sizeof(MPI_FLOAT) = {}, sizeof(float) = {}", size, sizeof(float));

      MPI_Type_size(MPI_DOUBLE, &size);
      log::verbose_once("sizeof(MPI_DOUBLE) = {}, sizeof(double) = {}", size, sizeof(double));
    }

    void mpiBarrier()
    {
      MPI_Barrier(getMPIComm());
    }

    void mpiBarrierUniqueDomains()
    {
      MPI_Barrier(getCommUniqueDomains());
    }

    void mpiBarrierSharedDomain()
    {
      MPI_Barrier(getCommSharedDomain());
    }

    int mpiAllreduce(const void *sendbuf, void *recvbuf, int count,
                     MPI_Datatype datatype, MPI_Op op)
    {
      return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, getMPIComm());
    }

    int mpiAllreduceUD(const void *sendbuf, void *recvbuf, int count,
                       MPI_Datatype datatype, MPI_Op op)
    {
      return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, getCommUniqueDomains());
    }

    int mpiAllreduceSD(const void *sendbuf, void *recvbuf, int count,
                       MPI_Datatype datatype, MPI_Op op)
    {
      return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, getCommSharedDomain());
    }

#endif // ENABLE_MPI_

  } // namespace mpi

} // namespace antmoc
