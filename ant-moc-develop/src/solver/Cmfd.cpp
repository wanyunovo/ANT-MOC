#include "antmoc/Cmfd.h"
#include "antmoc/Geometry.h"
#include "antmoc/Lattice.h"
#include "antmoc/LocalCoords.h"
#include "antmoc/Material.h"
#include "antmoc/Quadrature.h"
#include "antmoc/Timer.h"
#include "antmoc/Track3D.h"

#include <fstream>
#include <iomanip>

namespace antmoc
{

  /**
   * @brief Constructor initializes boundaries and variables that describe
   *          the CMFD object.
   * @details The construcor initializes the many variables that describe
   *          the CMFD mesh and are used to solve the nonlinear diffusion
   *          acceleration problem.
   */
  Cmfd::Cmfd()
  {

    /* Initialize Geometry and Mesh-related attribute 初始化几何和网格相关的属性*/
    _quadrature = NULL;
    _geometry = NULL;
    _materials = NULL;

    /* Communication structs */
    _convergence_data = NULL;
    _domain_communicator = NULL;

    /* Global variables used in solving CMFD problem */
    _source_convergence_threshold = 1E-5;
    _num_x = 1;
    _num_y = 1;
    _num_z = 1;
    _num_r = 1;
    _local_num_x = 1;
    _local_num_y = 1;
    _local_num_z = 1;
    _width_x = 0.;
    _width_y = 0.;
    _width_z = 0.;
    _width_r = 0.;
    _cell_width_x = 0.;
    _cell_width_y = 0.;
    _cell_width_z = 0.;
    _flux_update_on = true;
    _centroid_update_on = true;
    _use_axial_interpolation = 0;
    _flux_limiting = true;
    _balance_sigma_t = false;
    _k_nearest = 1;
    _SOR_factor = 1.0;
    _num_FSRs = 0;
    _empty_cells_num = 0;
#ifndef THREED
    _SOLVE_3D = false;
#endif
    _total_tally_size = 0;
    _tallies_allocated = false;
    _domain_communicator_allocated = false;
    _linear_source = false;
    _check_neutron_balance = false;
    _old_dif_surf_valid = false;
    _non_uniform = false;
    _widths_adjusted_for_domains = false;
    _hexlattice_enable = false;
    _orientation = Orientation::y;

    /* Energy group and polar angle problem parameters */
    _num_moc_groups = 0;
    _num_cmfd_groups = 0;
    _hex_num_groups = 0;
    _num_backup_groups = 1;
    _num_polar = 0;
    _num_azim = 0;

    /* Set matrices and arrays to NULL */
    _A = NULL;
    _M = NULL;
    _k_eff = 1.0;
    _relaxation_factor = 1.0;
    _old_flux = NULL;
    _new_flux = NULL;
    _old_dif_surf_corr = NULL;
    _old_source = NULL;
    _new_source = NULL;
    _flux_moments = NULL;
    _group_indices = NULL;
    _group_indices_map = NULL;
    _user_group_indices = false;
    _surface_currents = NULL;
    _starting_currents = NULL;
    _net_currents = NULL;
    _full_surface_currents = NULL;
    _cell_locks = NULL;
    _volumes = NULL;
    _lattice = NULL;
    _azim_spacings = NULL;
    _polar_spacings = NULL;
    _temporary_currents = NULL;
    _backup_cmfd = NULL;
    _cmfd_group_to_backup_group = NULL;
    _backup_group_structure.resize(0);

    /* Initialize boundaries to be reflective */
    /* reclattice have 6 surface, hexlattice have 8 surface */
    _boundaries = new boundaryType[8];

    _boundaries[0] = REFLECTIVE;
    _boundaries[1] = REFLECTIVE;
    _boundaries[2] = REFLECTIVE;
    _boundaries[3] = REFLECTIVE;
    _boundaries[4] = REFLECTIVE;
    _boundaries[5] = REFLECTIVE;
    _boundaries[6] = REFLECTIVE;
    _boundaries[7] = REFLECTIVE;

    /*
    _boundaries[0] = VACUUM;
    _boundaries[1] = VACUUM;
    _boundaries[2] = VACUUM;
    _boundaries[3] = VACUUM;
    _boundaries[4] = VACUUM;
    _boundaries[5] = VACUUM;
    _boundaries[6] = VACUUM;
    _boundaries[7] = VACUUM;
    */

    /* Initialize CMFD timer */
    _timer = new Timer();
  }

  /**
   * @brief Destructor deletes arrays of A and M row insertion arrays.
   */
  Cmfd::~Cmfd()
  {

    if (_cell_locks != NULL)
      delete[] _cell_locks;

    if (_boundaries != NULL)
      delete[] _boundaries;

    if (_group_indices != NULL)
      delete[] _group_indices;

    if (_group_indices_map != NULL)
      delete[] _group_indices_map;

    /* Delete the Matrix and Vector objects */
    if (_M != NULL)
      delete _M;

    if (_A != NULL)
      delete _A;

    if (_old_source != NULL)
      delete _old_source;

    if (_new_source != NULL)
      delete _new_source;

    if (_old_flux != NULL)
      delete _old_flux;

    if (_new_flux != NULL)
      delete _new_flux;

    if (_surface_currents != NULL)
      delete _surface_currents;

    if (_starting_currents != NULL)
      delete _starting_currents;

    if (_net_currents != NULL)
      delete _net_currents;

    if (_full_surface_currents != NULL)
      delete _full_surface_currents;

    if (_volumes != NULL)
      delete _volumes;

    if (_azim_spacings != NULL)
      delete[] _azim_spacings;

    if (_polar_spacings != NULL)
    {
      for (int a = 0; a < _num_azim / 4; a++)
        delete[] _polar_spacings[a];
      delete[] _polar_spacings;
    }

    /* Delete CMFD materials array */
    if (_materials != NULL)
    {
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
        delete _materials[i];
      delete[] _materials;
    }

    /* Delete the CMFD lattice */
    if (_lattice != NULL)
      delete _lattice;

    /* Clear the _cell_fsrs vector of vectors */
    std::vector<std::vector<long>>::iterator iter1;
    for (iter1 = _cell_fsrs.begin(); iter1 != _cell_fsrs.end(); ++iter1)
      iter1->clear();
    _cell_fsrs.clear();

    /* Clear the _k_nearest_stencils map of vectors */
    std::map<int, std::vector<std::pair<int, double>>>::iterator iter2;
    for (iter2 = _k_nearest_stencils.begin(); iter2 != _k_nearest_stencils.end();
         ++iter2)
      iter2->second.clear();
    _k_nearest_stencils.clear();

    /* Delete tally information */
    if (_tallies_allocated)
    {

      delete[] _tally_memory;

      delete[] _reaction_tally;
      delete[] _volume_tally;
      delete[] _diffusion_tally;
    }

    int num_threads = omp_get_max_threads();
    if (_temporary_currents != NULL)
    {
      for (int t = 0; t < num_threads; t++)
        delete[] _temporary_currents[t];
      delete[] _temporary_currents;
    }

    /* De-allocate domain communicator */
    // unused
    // int num_cells_local = _local_num_x * _local_num_y * _local_num_z;
    if (_domain_communicator != NULL)
    {
      if (_domain_communicator_allocated)
      {
        for (int rb = 0; rb < 2; rb++)
        {
          for (int f = 0; f < NUM_FACES; f++)
          {
            delete[] _domain_communicator->indexes[rb][f];
            delete[] _domain_communicator->domains[rb][f];
            delete[] _domain_communicator->coupling_coeffs[rb][f];
            delete[] _domain_communicator->fluxes[rb][f];
          }
          delete[] _domain_communicator->num_connections[rb];
          delete[] _domain_communicator->indexes[rb];
          delete[] _domain_communicator->domains[rb];
          delete[] _domain_communicator->coupling_coeffs[rb];
          delete[] _domain_communicator->fluxes[rb];
        }

        delete[] _domain_communicator->num_connections;
        delete[] _domain_communicator->indexes;
        delete[] _domain_communicator->fluxes;
        delete[] _domain_communicator->coupling_coeffs;
        delete _domain_communicator;

        delete[] _inter_domain_data;
        for (int s = 0; s < NUM_FACES; s++)
        {
          delete[] _boundary_volumes[s];
          delete[] _boundary_reaction[s];
          delete[] _boundary_diffusion[s];
          delete[] _old_boundary_flux[s];
          delete[] _boundary_surface_currents[s];
        }

        delete[] _boundary_volumes;
        delete[] _boundary_reaction;
        delete[] _boundary_diffusion;
        delete[] _old_boundary_flux;
        delete[] _boundary_surface_currents;
      }
    }

    for (size_t r = 0; r < _axial_interpolants.size(); r++)
      delete[] _axial_interpolants.at(r);

    if (_backup_cmfd != NULL)
      delete _backup_cmfd;

    delete _timer;
  }

  /**
   * @brief Set the bool _hexlattice_enable.
   * @param bool create HexLattice or RecLattice
   */
  void Cmfd::setHexLatticeEnable(bool hexlattice_enable)
  {

    if (hexlattice_enable)
    { // Hexagonal
      _hexlattice_enable = true;
    }
    else
    { // Rectangular
      _hexlattice_enable = false;
    }
  }

  bool Cmfd::GetHexLatticeEnable()
  {
    return _hexlattice_enable;
  }

  void Cmfd::convertEdgeToSurfaces(int surface, int *partial_surfaces)
  {
    if (surface == HEX_SURFACE_BETA_MIN_GAMMA_MIN)
    {
      partial_surfaces[0] = 0;
      partial_surfaces[1] = 1;
    }
    else if (surface == HEX_SURFACE_BETA_MIN_DELTA_MIN)
    {
      partial_surfaces[0] = 0;
      partial_surfaces[1] = 2;
    }
    else if (surface == HEX_SURFACE_GAMMA_MIN_DELTA_MAX)
    {
      partial_surfaces[0] = 1;
      partial_surfaces[1] = 2;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_GAMMA_MAX)
    {
      partial_surfaces[0] = 4;
      partial_surfaces[1] = 5;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_DELTA_MAX)
    {
      partial_surfaces[0] = 4;
      partial_surfaces[1] = 6;
    }
    else if (surface == HEX_SURFACE_GAMMA_MAX_DELTA_MIN)
    {
      partial_surfaces[0] = 5;
      partial_surfaces[1] = 6;
    }
    else if (surface == HEX_SURFACE_BETA_MIN_Z_MIN)
    {
      partial_surfaces[0] = 0;
      partial_surfaces[1] = 3;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_Z_MIN)
    {
      partial_surfaces[0] = 4;
      partial_surfaces[1] = 3;
    }
    else if (surface == HEX_SURFACE_BETA_MIN_Z_MAX)
    {
      partial_surfaces[0] = 0;
      partial_surfaces[1] = 7;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_Z_MAX)
    {
      partial_surfaces[0] = 4;
      partial_surfaces[1] = 7;
    }
    else if (surface == HEX_SURFACE_GAMMA_MIN_Z_MIN)
    {
      partial_surfaces[0] = 1;
      partial_surfaces[1] = 3;
    }
    else if (surface == HEX_SURFACE_GAMMA_MAX_Z_MIN)
    {
      partial_surfaces[0] = 5;
      partial_surfaces[1] = 3;
    }
    else if (surface == HEX_SURFACE_GAMMA_MIN_Z_MAX)
    {
      partial_surfaces[0] = 1;
      partial_surfaces[1] = 7;
    }
    else if (surface == HEX_SURFACE_GAMMA_MAX_Z_MAX)
    {
      partial_surfaces[0] = 5;
      partial_surfaces[1] = 7;
    }
    else if (surface == HEX_SURFACE_DELTA_MIN_Z_MIN)
    {
      partial_surfaces[0] = 2;
      partial_surfaces[1] = 3;
    }
    else if (surface == HEX_SURFACE_DELTA_MAX_Z_MIN)
    {
      partial_surfaces[0] = 6;
      partial_surfaces[1] = 3;
    }
    else if (surface == HEX_SURFACE_DELTA_MIN_Z_MAX)
    {
      partial_surfaces[0] = 2;
      partial_surfaces[1] = 7;
    }
    else if (surface == HEX_SURFACE_DELTA_MAX_Z_MAX)
    {
      partial_surfaces[0] = 6;
      partial_surfaces[1] = 7;
    }
  }

  void Cmfd::convertVertexToSurfaces(int surface, int *remainder_surfaces, int *partial_surfaces)
  {
    if (surface == HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MIN)
    {
      remainder_surfaces[0] = 8;
      remainder_surfaces[1] = 18;
      remainder_surfaces[2] = 14;
      partial_surfaces[0] = 3;
      partial_surfaces[1] = 0;
      partial_surfaces[2] = 1;
    }
    else if (surface == HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MAX)
    {
      remainder_surfaces[0] = 8;
      remainder_surfaces[1] = 20;
      remainder_surfaces[2] = 16;
      partial_surfaces[0] = 7;
      partial_surfaces[1] = 0;
      partial_surfaces[2] = 1;
    }
    else if (surface == HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MIN)
    {
      remainder_surfaces[0] = 9;
      remainder_surfaces[1] = 22;
      remainder_surfaces[2] = 14;
      partial_surfaces[0] = 3;
      partial_surfaces[1] = 0;
      partial_surfaces[2] = 2;
    }
    else if (surface == HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MAX)
    {
      remainder_surfaces[0] = 9;
      remainder_surfaces[1] = 24;
      remainder_surfaces[2] = 16;
      partial_surfaces[0] = 7;
      partial_surfaces[1] = 0;
      partial_surfaces[2] = 2;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MIN)
    {
      remainder_surfaces[0] = 11;
      remainder_surfaces[1] = 19;
      remainder_surfaces[2] = 15;
      partial_surfaces[0] = 3;
      partial_surfaces[1] = 4;
      partial_surfaces[2] = 5;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MAX)
    {
      remainder_surfaces[0] = 11;
      remainder_surfaces[1] = 21;
      remainder_surfaces[2] = 17;
      partial_surfaces[0] = 7;
      partial_surfaces[1] = 4;
      partial_surfaces[2] = 5;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MIN)
    {
      remainder_surfaces[0] = 12;
      remainder_surfaces[1] = 23;
      remainder_surfaces[2] = 15;
      partial_surfaces[0] = 3;
      partial_surfaces[1] = 4;
      partial_surfaces[2] = 6;
    }
    else if (surface == HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MAX)
    {
      remainder_surfaces[0] = 12;
      remainder_surfaces[1] = 25;
      remainder_surfaces[2] = 17;
      partial_surfaces[0] = 7;
      partial_surfaces[1] = 4;
      partial_surfaces[2] = 6;
    }
    else if (surface == HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MIN)
    {
      remainder_surfaces[0] = 10;
      remainder_surfaces[1] = 23;
      remainder_surfaces[2] = 18;
      partial_surfaces[0] = 3;
      partial_surfaces[1] = 1;
      partial_surfaces[2] = 6;
    }
    else if (surface == HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MAX)
    {
      remainder_surfaces[0] = 10;
      remainder_surfaces[1] = 25;
      remainder_surfaces[2] = 20;
      partial_surfaces[0] = 7;
      partial_surfaces[1] = 1;
      partial_surfaces[2] = 6;
    }
    else if (surface == HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MIN)
    {
      remainder_surfaces[0] = 13;
      remainder_surfaces[1] = 22;
      remainder_surfaces[2] = 19;
      partial_surfaces[0] = 3;
      partial_surfaces[1] = 5;
      partial_surfaces[2] = 2;
    }
    else if (surface == HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MAX)
    {
      remainder_surfaces[0] = 13;
      remainder_surfaces[1] = 24;
      remainder_surfaces[2] = 21;
      partial_surfaces[0] = 7;
      partial_surfaces[1] = 5;
      partial_surfaces[2] = 2;
    }
  }

  /// \brief Set the orientation of the HexLattice.
  void Cmfd::setOrientation(Orientation orientation)
  {
    _orientation = orientation;
  }

  /// \brief Set the orientation of the HexLattice.
  /// \details The orientation is defined according to the direction
  ///          of the flat side.
  void Cmfd::setOrientation(std::string orientation)
  {
    stringutils::trim(orientation);
    stringutils::toUpper(orientation);

    if (orientation == "Y")
    {
      _orientation = Orientation::y;
    }
    else
    {
      _orientation = Orientation::x;
    }
  }

  /**
   * @brief Set the number of Mesh cells in a row.
   * @param number of Mesh cells in a row
   */
  void Cmfd::setNumX(int num_x)
  {

    if (num_x < 1)
      log::ferror("The number of lattice cells in the x direction "
                  "must be > 0. Input value: %i",
                  num_x);

    _num_x = num_x;        // 整个几何上x方向上CMFD网格的数量
    _local_num_x = _num_x; // 本进程（MPI并行中的单个进程）负责的X方向网格数（默认等于全局）

    /*
    _domain_communicator：MPI通信器对象（用于多进程协调）
    如果启用了MPI并行：
    _num_domains_x：X方向被分成几个子域（进程数）
    将全局网格数均分给各个进程
    例如：100个网格分给4个进程 → 每个进程25个
    */
    if (_domain_communicator != NULL)
      _local_num_x = _num_x / _domain_communicator->_num_domains_x;

    // 计算单个网格宽度
    if (_width_x != 0.)
      _cell_width_x = _width_x / _num_x; // _width_x整个几何上X方向的宽度/整个几何上x方向上CMFD网格的数量，就是一个cmfd网格在x方向上的长度
  }

  /**
   * @brief Set the number of Mesh cells in a column
   * @param number of Mesh cells in a column
   */
  void Cmfd::setNumY(int num_y)
  {

    if (num_y < 1)
      log::ferror("The number of lattice cells in the y direction "
                  "must be > 0. Input value: %i",
                  num_y);

    _num_y = num_y;
    _local_num_y = _num_y;
    if (_domain_communicator != NULL)
      _local_num_y = _num_y / _domain_communicator->_num_domains_y;
    if (_width_y != 0.)
      _cell_width_y = _width_y / _num_y;
  }

  /**
   * @brief Set the number of Mesh cells in a column
   * @param number of Mesh cells in a column
   */
  void Cmfd::setNumZ(int num_z)
  {

    if (num_z < 1)
      log::ferror("The number of lattice cells in the z direction "
                  "must be > 0. Input value: %i",
                  num_z);

    _num_z = num_z;
    _local_num_z = _num_z;
    if (_domain_communicator != NULL)
      _local_num_z = _num_z / _domain_communicator->_num_domains_z; // 每个域中的z轴晶格数量
    if (_width_z != 0.)
      _cell_width_z = _width_z / _num_z;
    if (_width_z == std::numeric_limits<double>::infinity()) // 如果整个几何的高度为正无穷大，则每个CMFD网格的高度设为1
      _cell_width_z = 1.0;
  }

  /**
   * @brief Set the number of Mesh cells in a column
   * @param number of Mesh cells in a column
   */
  void Cmfd::setNumR(int num_r)
  {

    if (num_r < 1)
      log::ferror("The number of lattice cells in the r direction "
                  "must be > 0. Input value: %i",
                  num_r);

    _num_r = num_r;
  }

  /**
   * @brief Get the number of Mesh cells in a row.
   * @return number of Mesh cells in a row
   */
  int Cmfd::getNumX()
  {
    return _num_x;
  }

  /**
   * @brief Get the number of Mesh cells in a column
   * @return number of Mesh cells in a column
   */
  int Cmfd::getNumY()
  {
    return _num_y;
  }

  /**
   * @brief Get the number of Mesh cells in the z-direction
   * @return number of Mesh cells in the z-direction
   */
  int Cmfd::getNumZ()
  {
    return _num_z;
  }

  /**
   * @brief Get the Vector of surface currents.
   * @return pointer to a vector containing the surface currents.
   */
  Vector *Cmfd::getLocalCurrents()
  {
    return _surface_currents;
  }

  /**
   * @brief Get the array of surface currents on the boundaries.
   * @return 3D array containing the boundary surface currents.
   */
  CMFD_PRECISION ***Cmfd::getBoundarySurfaceCurrents()
  {
    return _boundary_surface_currents;
  }

  /**
   * @brief Set Mesh width in the x-direction
   * 设置 CMFD 网格在 x 方向的总几何宽度，并在已知单元数时同步更新每个单元的平均宽度。
   * @param width Physical width of Mesh in the x-direction
   */
  void Cmfd::setWidthX(double width)
  {
    _width_x = width;
    if (_num_x != 0)
      _cell_width_x = _width_x / _num_x;
  }

  /**
   * @brief Set Mesh width in the y-direction
   * @param width Physical width of Mesh in the y-direction
   */
  void Cmfd::setWidthY(double width)
  {
    _width_y = width;
    if (_num_y != 0)
      _cell_width_y = _width_y / _num_y;
  }

  /**
   * @brief Set Mesh width in the z-direction
   * @param width Physical width of Mesh in the z-direction
   */
  void Cmfd::setWidthZ(double width)
  {
    _width_z = width;
    if (_num_z != 0)
      _cell_width_z = _width_z / _num_z;
  }

  /**
   * @brief Set Mesh width in the z-direction
   * @param width Physical width of Mesh in the z-direction
   */
  void Cmfd::setWidthR(double width)
  {
    _width_r = width;
  }

  /**
   * @brief Set Mesh width in the z-direction
   * @param width Physical width of Mesh in the z-direction
   */
  void Cmfd::setWidthsZ(DoubleVec widths)
  {
    _non_uniform = true;
    _widths_z = widths;
  }

#ifdef ENABLE_MPI_
  /**
   * @brief Set the number of domains in each direction.
   * @param num_x number of domains in the X direction
   * @param num_y number of domains in the Y direction
   * @param num_z number of domains in the Z direction
   */
  void Cmfd::setNumDomains(int num_x, int num_y, int num_z)
  {

    if (_domain_communicator == NULL)
    {
      _domain_communicator = new DomainCommunicator; // DomainCommunicator：存储MPI通信相关信息的结构体
      _domain_communicator->_MPI_cart = _geometry->getMPICart();
    }

    _domain_communicator->_num_domains_x = num_x;
    _domain_communicator->_num_domains_y = num_y;
    _domain_communicator->_num_domains_z = num_z;

    _local_num_x = _num_x / num_x;
    _local_num_y = _num_y / num_y;
    _local_num_z = _num_z / num_z;
  }

  /**
   * @brief Set the indexes of the domain among the global lattice of domains.
   */
  void Cmfd::setDomainIndexes(int idx_x, int idx_y, int idx_z)
  {

    if (_domain_communicator == NULL)
    {
      _domain_communicator = new DomainCommunicator;
      _domain_communicator->_MPI_cart = _geometry->getMPICart();
    }

    _domain_communicator->_domain_idx_x = idx_x;
    _domain_communicator->_domain_idx_y = idx_y;
    _domain_communicator->_domain_idx_z = idx_z;
  }
#endif

  /**
   * @brief Collapse cross-sections and fluxes for each CMFD cell by
   *        energy condensing and volume averaging cross sections from
   *        the MOC sweep.
   * 通过MOC扫描的能群压缩和体积平均截面，对每个CMFD单元的横截面和通量数据进行压缩 (就是将 MOC 细网格数据压缩到 CMFD 粗网格上，以便加速计算过程)
   * @details This method performs a cell-wise energy condensation and volume
   *          average of the cross sections of the fine, unstructured FSR mesh.
   *          The cross sections are condensed such that all reaction rates and
   *          the neutron production rate from fission are conserved. It is
   *          important to note that the volume averaging is performed before
   *          energy condensation in order to properly collapse the diffusion
   *          coefficients.
   * 该方法按单元进行能量缩并和体积平均，从MOC细网格FSR中获取截面数据。截面数据的压缩确保了所有反应率和来自裂变的中子生产率保持不变
   * 为了正确地压缩扩散系数，体积平均在能量缩并之前执行
   */
  void Cmfd::collapseXS()
  {

    log::fverbose("Collapsing cross-sections onto CMFD mesh...");

    /* Record net currents over cells if neutron balance of sigma-t requested */
    if (_balance_sigma_t)
      recordNetCurrents(); // 统计净中子流数据

    /* Check to see that CMFD tallies have been allocated */
    if (!_tallies_allocated)
      log::ferror("Tallies need to be allocated before collapsing "
                  "cross-sections");

    /* Split vertex and edge currents to side surfaces */
    splitVertexCurrents(); // 顶点拆分
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
      communicateSplits(false);
#endif
    splitEdgeCurrents(); // 边拆分
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
      communicateSplits(true); // 通信域的流拆分,传入true代表从接收的电流数据中提取相关的面电流，并将其累加到 _surface_currents数据中
#endif

#pragma omp parallel
    {

      /* Initialize variables for FSR properties*/
      FP_PRECISION volume, flux;
      FP_PRECISION tot, nu_fis, chi; // 其中nu_fis表示平均裂变中子数*裂变截面,chi表示中子裂变能谱
      FP_PRECISION *scat;

      double *scat_tally = new double[_num_cmfd_groups];
      double *chi_tally = new double[_num_cmfd_groups];

      /* Pointers to material objects */
      Material *fsr_material;
      Material *cell_material;

      /* Loop over CMFD cells */
      /*
       * 遍历每个CMFD网格，进行截面和通量统计
       * 原理对应：限制操作的核心 - 将多个FSR的数据聚合到单个CMFD网格
       */
#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        std::vector<long>::iterator iter;
        cell_material = _materials[i];

        /* Zero group-wise fission terms */
        double neutron_production_tally = 0.0;
        for (int e = 0; e < _num_cmfd_groups; e++)
          chi_tally[e] = 0.0;

        /* Loop over FSRs in CMFD cell */
        for (iter = _cell_fsrs.at(i).begin();
             iter != _cell_fsrs.at(i).end(); ++iter)
        {
          fsr_material = _FSR_materials[*iter];
          volume = _FSR_volumes[*iter];

          /* Calculate total neutron production in the FSR */
          /*
           * 计算FSR的中子产生率
           * 用于后续计算裂变能谱的体积加权平均
           * 公式对应：文档中的体积加权平均公式
           */
          double neutron_production = 0.0; // 每个FSR的中子产生
          for (int h = 0; h < _num_moc_groups; h++)
          { // 遍历该FSR的moc能群
            neutron_production += fsr_material->getNuSigmaFByGroup(h + 1) *
                                  _FSR_fluxes[(*iter) * _num_moc_groups + h] * volume; // 产生中子=FSR材料的平均裂变中子*裂变截面*FSR通量*FSR体积
            /*
            if(fsr_material->getNuSigmaFByGroup(h+1) > 0 || neutron_production > 0) {
              log::fdebug("cell %d fsr %d group %d, flux is %.2f, volume is %.2f, nusigmaf is %.2f, neutron_production is %.6f",
              i, *iter, h, _FSR_fluxes[(*iter)*_num_moc_groups+h], volume, fsr_material->getNuSigmaFByGroup(h+1), neutron_production);
            }
            */
          }

          /* Calculate contribution to all CMFD groups */
          /*
           * 计算chi（裂变能谱）的体积加权平均
           * 将MOC多群的chi压缩到CMFD少群
           * 文档B.10 - CMFD群裂变能谱的计算
           */
          for (int e = 0; e < _num_cmfd_groups; e++)
          {
            chi = 0;
            for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++) // 遍历对应CMFD能群下的MOC能群
              chi += fsr_material->getChiByGroup(h + 1);                    // 获取FSR材料的中子裂变能谱
            chi_tally[e] += chi * neutron_production;                       // 每个CMFD能群的中子裂变能谱统计=(对应MOC能群的chi总和)*产生中子
          }

          /* Add to total neutron production within the CMFD cell */
          neutron_production_tally += neutron_production; // 每个CMFD的产生中子统计=该网格下所有FSR的产生中子总和
        }

        /*
        if(neutron_production_tally > 0) {
          log::fdebug("Cell id is %d, neutron production tally is %.2f", i, neutron_production_tally);
        }
        */

        /* Set chi 设置CMFD材料的 chi 值*/
        if (fabs(neutron_production_tally) > FLT_EPSILON)
        {

          /* Calculate group-wise fission contriubtions */
          for (int e = 0; e < _num_cmfd_groups; e++)
            cell_material->setChiByGroup(chi_tally[e] / neutron_production_tally,
                                         e + 1); // 对应公式B.10
        }
        else
        {
          /* Calculate group-wise chi to zero */
          for (int e = 0; e < _num_cmfd_groups; e++)
            cell_material->setChiByGroup(0.0, e + 1);
        }

        /* ========================================
         * 计算总截面、散射截面、裂变截面的通量体积加权平均
         * ========================================
         * 原理对应：B.11-B.13公式
         *
         * 关键点：
         * - 通量体积加权：Σ = (Σ_flux*volume) / (flux*volume)
         * - 这保证了反应率守恒：Σ_new * φ_CMFD = Σ_old * φ_FSR（在网格内）
         */
        /* Loop over CMFD coarse energy groups 遍历CMFD能组，计算各项统计数据 ,后缀为tally皆为从FSR数据累加到CMFD粗网格的平均值*/
        for (int e = 0; e < _num_cmfd_groups; e++)
        {

          /* Zero tallies for this group */
          double nu_fission_tally = 0.0;
          double total_tally = 0.0;

          _diffusion_tally[i][e] = 0.0; // CMFD
          _reaction_tally[i][e] = 0.0;
          _volume_tally[i][e] = 0.0;

          /* Zero each group-to-group scattering tally */
          for (int g = 0; g < _num_cmfd_groups; g++)
            scat_tally[g] = 0.0;

          /* Loop over MOC energy groups within this CMFD coarse group 遍历CMFD能群对应的MOC能群*/
          for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
          {

            /* Reset volume tally for this MOC group */
            _volume_tally[i][e] = 0.0;
            double rxn_tally_group = 0.0;
            double trans_tally_group = 0.0;

            /* Loop over FSRs in CMFD cell */
            for (iter = _cell_fsrs.at(i).begin();
                 iter != _cell_fsrs.at(i).end(); ++iter)
            {

              /* Gets FSR volume, material, and cross sections */
              fsr_material = _FSR_materials[*iter];
              volume = _FSR_volumes[*iter];
              scat = fsr_material->getSigmaS(); // 散射截面
              flux = _FSR_fluxes[(*iter) * _num_moc_groups + h];
              tot = fsr_material->getSigmaTByGroup(h + 1);      // FSR对应能群总截面
              nu_fis = fsr_material->getNuSigmaFByGroup(h + 1); // FSR对应能群裂变截面*平均裂变中子数

              /* Increment tallies for this group */
              total_tally += tot * flux * volume; // FSR总截面*通量*体积
              nu_fission_tally += nu_fis * flux * volume;
              _reaction_tally[i][e] += flux * volume;
              _volume_tally[i][e] += volume;
              /*log::fdebug("cell is %d, group is %d, flux is %.2f, volume is %.2f, reaction tally is %.6f",
              i, e, flux, volume, _reaction_tally[i][e]);*/

              /* Increment diffusion MOC group-wise tallies */
              rxn_tally_group += flux * volume;         // FSR通量*体积累加
              trans_tally_group += tot * flux * volume; // FSR总截面*通量*体积累加

              /* Scattering tallies */
              for (int g = 0; g < _num_moc_groups; g++)
              {
                scat_tally[getCmfdGroup(g)] +=
                    scat[g * _num_moc_groups + h] * flux * volume;
              }
            }
            if (fabs(rxn_tally_group) > FLT_EPSILON &&
                fabs(trans_tally_group) > FLT_EPSILON)
            {
              CMFD_PRECISION flux_avg_sigma_t = trans_tally_group /
                                                rxn_tally_group; // B.13,求CMFD总截面
              _diffusion_tally[i][e] += rxn_tally_group /
                                        (3.0 * flux_avg_sigma_t); // B.24 好像只对应分子
            }
          }

          /* Save cross-sections to material */
          double rxn_tally = _reaction_tally[i][e]; //

          if (fabs(rxn_tally) < FLT_EPSILON)
          {
            log::fwarn("Zero reaction tally calculated in CMFD cell %d "
                       "in CMFD group %d",
                       i, e);
            rxn_tally = ZERO_SIGMA_T;
            _reaction_tally[i][e] = ZERO_SIGMA_T;
            _diffusion_tally[i][e] = ZERO_SIGMA_T;
          }

          cell_material->setSigmaTByGroup(total_tally / rxn_tally, e + 1);        // 设置CMFD材料对应能群总截面
          cell_material->setNuSigmaFByGroup(nu_fission_tally / rxn_tally, e + 1); // 设置CMFD材料对应能群裂变截面

          /* Set scattering xs 设置散射截面*/
          for (int g = 0; g < _num_cmfd_groups; g++)
          {
            cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, e + 1,
                                            g + 1); //_sigma_s[_num_groups*(destination-1) + (origin-1)] = xs;
          }
        }
      }
      delete[] scat_tally;
      delete[] chi_tally;
    }

#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
    {
      if (_domain_communicator != NULL)
      {

        /* Start recording MPI communication time */
        _timer->startTimer();

        /* Do the Ghost cell exchange */
        ghostCellExchange();

        /* Tally the MPI communication time */
        _timer->stopTimer("CMFD MPI communication time");
      }
    }
#endif

    /* Calculate (local) old fluxes and set volumes */
#pragma omp parallel for
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {

      /* Loop over CMFD coarse energy groups */
      for (int e = 0; e < _num_cmfd_groups; e++)
      {

        /* Load tallies at this cell and energy group */
        double vol_tally = _volume_tally[i][e];
        double rxn_tally = _reaction_tally[i][e];
        _old_flux->setValue(i, e, rxn_tally / vol_tally);

        /* Set the Mesh cell properties with the tallies */
        _volumes->setValue(i, 0, vol_tally);
      }
    }

    /* Loop over boundary CMFD cells and set cross sections */
    if (mpi::isSpatialDecomposed())
    {
#pragma omp parallel for
      for (int s = 0; s < NUM_FACES; s++)
      {

        /* Loop over all CMFD cells on the current surface */
        std::map<int, int>::iterator it;
        for (it = _boundary_index_map.at(s).begin();
             it != _boundary_index_map.at(s).end(); ++it)
        {

          int idx = it->second;

          /* Loop over CMFD coarse energy groups */
          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            /* Load tallies at this cell and energy group */
            double vol_tally = _boundary_volumes[s][idx][0];
            double rxn_tally = _boundary_reaction[s][idx][e];
            _old_boundary_flux[s][idx][e] = rxn_tally / vol_tally;
          }
        }
      }
    }
  }

  /**
   * @brief Collapse cross-sections and fluxes for each CMFD cell by
   *        energy condensing and volume averaging cross sections from
   *        the MOC sweep.
   * @details This method performs a cell-wise energy condensation and volume
   *          average of the cross sections of the fine, unstructured FSR mesh.
   *          The cross sections are condensed such that all reaction rates and
   *          the neutron production rate from fission are conserved. It is
   *          important to note that the volume averaging is performed before
   *          energy condensation in order to properly collapse the diffusion
   *          coefficients.
   */
  void Cmfd::hexCollapseXS()
  {

    log::fverbose("Collapsing cross-sections onto CMFD mesh...");

    /* Record net currents over cells if neutron balance of sigma-t requested */
    if (_balance_sigma_t)
      recordNetCurrents();

    /* Check to see that CMFD tallies have been allocated */
    if (!_tallies_allocated)
      log::ferror("Tallies need to be allocated before collapsing "
                  "cross-sections");

    /* Split vertex and edge currents to side surfaces */
    if (_hexlattice_enable)
      splitVertexCurrentsHex();
    else
      splitVertexCurrents();

#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
      communicateSplits(false);
#endif
    if (_hexlattice_enable)
      splitEdgeCurrentsHex();
    else
      splitEdgeCurrents();
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
      communicateSplits(true);
#endif

#pragma omp parallel
    {

      /* Initialize variables for FSR properties*/
      FP_PRECISION volume, flux;
      FP_PRECISION tot, nu_fis, chi;
      FP_PRECISION *scat;

      double *scat_tally = new double[_num_cmfd_groups];
      double *chi_tally = new double[_num_cmfd_groups];

      /* Pointers to material objects */
      Material *fsr_material;
      Material *cell_material;

      /* Loop over CMFD cells */
#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {
        if (_empty_fsrs_cells[i])
        {
          log::fdebug("CMFD Cell %d is empty, dont have any fsrs", i);
          continue;
        }
        else
        {
          int j = _logical_actual_map[i];
          std::vector<long>::iterator iter;
          cell_material = _materials[j];

          /* Zero group-wise fission terms */
          double neutron_production_tally = 0.0;
          for (int e = 0; e < _num_cmfd_groups; e++)
            chi_tally[e] = 0.0;

          /* Loop over FSRs in CMFD cell */
          for (iter = _cell_fsrs.at(i).begin();
               iter != _cell_fsrs.at(i).end(); ++iter)
          {
            fsr_material = _FSR_materials[*iter];
            volume = _FSR_volumes[*iter];

            /* Calculate total neutron production in the FSR */
            double neutron_production = 0.0;
            for (int h = 0; h < _num_moc_groups; h++)
            {
              neutron_production += fsr_material->getNuSigmaFByGroup(h + 1) *
                                    _FSR_fluxes[(*iter) * _num_moc_groups + h] * volume;
              /*
              if(fsr_material->getNuSigmaFByGroup(h+1) > 0 || neutron_production > 0) {
                log::fdebug("cell %d fsr %d group %d, flux is %.2f, volume is %.2f, nusigmaf is %.2f, neutron_production is %.6f",
                i, *iter, h, _FSR_fluxes[(*iter)*_num_moc_groups+h], volume, fsr_material->getNuSigmaFByGroup(h+1), neutron_production);
              }
              */
            }

            /* Calculate contribution to all CMFD groups */
            for (int e = 0; e < _num_cmfd_groups; e++)
            {
              chi = 0;
              for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
                chi += fsr_material->getChiByGroup(h + 1);

              chi_tally[e] += chi * neutron_production;
              // log::fdebug("CMFD Cell %d chi_tally[%d] is %.6f", j, e, chi_tally[e]);
            }

            /* Add to total neutron production within the CMFD cell */
            neutron_production_tally += neutron_production;
          }

          /*
          if(neutron_production_tally > 0) {
            log::fdebug("Cell id is %d, neutron production tally is %.2f", i, neutron_production_tally);
          }
          */

          /* Set chi */
          if (fabs(neutron_production_tally) > FLT_EPSILON)
          {

            /* Calculate group-wise fission contriubtions */
            for (int e = 0; e < _num_cmfd_groups; e++)
              cell_material->setChiByGroup(chi_tally[e] / neutron_production_tally,
                                           e + 1);
          }
          else
          {
            /* Calculate group-wise chi to zero */
            for (int e = 0; e < _num_cmfd_groups; e++)
              cell_material->setChiByGroup(0.0, e + 1);
          }

          /* Loop over CMFD coarse energy groups */
          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            /* Zero tallies for this group */
            double nu_fission_tally = 0.0;
            double total_tally = 0.0;

            _diffusion_tally[j][e] = 0.0;
            _reaction_tally[j][e] = 0.0;
            _volume_tally[j][e] = 0.0;

            /* Zero each group-to-group scattering tally */
            for (int g = 0; g < _num_cmfd_groups; g++)
              scat_tally[g] = 0.0;

            /* Loop over MOC energy groups within this CMFD coarse group */
            for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
            {

              /* Reset volume tally for this MOC group */
              _volume_tally[j][e] = 0.0;
              double rxn_tally_group = 0.0;
              double trans_tally_group = 0.0;

              /* Loop over FSRs in CMFD cell */
              for (iter = _cell_fsrs.at(i).begin();
                   iter != _cell_fsrs.at(i).end(); ++iter)
              {

                /* Gets FSR volume, material, and cross sections */
                fsr_material = _FSR_materials[*iter];
                volume = _FSR_volumes[*iter];
                scat = fsr_material->getSigmaS();
                flux = _FSR_fluxes[(*iter) * _num_moc_groups + h];
                tot = fsr_material->getSigmaTByGroup(h + 1);
                nu_fis = fsr_material->getNuSigmaFByGroup(h + 1);

                /* Increment tallies for this group */
                total_tally += tot * flux * volume;
                nu_fission_tally += nu_fis * flux * volume;
                _reaction_tally[j][e] += flux * volume;
                _volume_tally[j][e] += volume;
                /*log::fdebug("cell is %d, group is %d, flux is %.2f, volume is %.2f, reaction tally is %.6f",
                i, e, flux, volume, _reaction_tally[i][e]);*/

                /* Increment diffusion MOC group-wise tallies */
                rxn_tally_group += flux * volume;
                trans_tally_group += tot * flux * volume;

                /* Scattering tallies */
                for (int g = 0; g < _num_moc_groups; g++)
                {
                  scat_tally[getCmfdGroup(g)] +=
                      scat[g * _num_moc_groups + h] * flux * volume;
                }
              }
              if (fabs(rxn_tally_group) > FLT_EPSILON &&
                  fabs(trans_tally_group) > FLT_EPSILON)
              {
                CMFD_PRECISION flux_avg_sigma_t = trans_tally_group /
                                                  rxn_tally_group;
                _diffusion_tally[j][e] += rxn_tally_group /
                                          (3.0 * flux_avg_sigma_t);
              }
            }

            /* Save cross-sections to material */
            double rxn_tally = _reaction_tally[j][e];

            if (fabs(rxn_tally) < FLT_EPSILON)
            {
              log::fwarn("Zero reaction tally calculated in CMFD cell %d "
                         "in CMFD group %d",
                         i, e);
              rxn_tally = ZERO_SIGMA_T;
              _reaction_tally[j][e] = ZERO_SIGMA_T;
              _diffusion_tally[j][e] = ZERO_SIGMA_T;
            }

            cell_material->setSigmaTByGroup(total_tally / rxn_tally, e + 1);
            cell_material->setNuSigmaFByGroup(nu_fission_tally / rxn_tally, e + 1);

            /* Set scattering xs */
            for (int g = 0; g < _num_cmfd_groups; g++)
            {
              cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, e + 1,
                                              g + 1);
            }
          }
        }
      }
      delete[] scat_tally;
      delete[] chi_tally;
    }

#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
    {
      if (_domain_communicator != NULL)
      {

        /* Start recording MPI communication time */
        _timer->startTimer();

        /* Do the Ghost cell exchange */
        ghostCellExchange();

        /* Tally the MPI communication time */
        _timer->stopTimer("CMFD MPI communication time");
      }
    }
#endif

    /* Calculate (local) old fluxes and set volumes */
#pragma omp parallel for
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {
      if (_empty_fsrs_cells[i])
        continue;
      else
      {
        int j = _logical_actual_map[i];
        /* Loop over CMFD coarse energy groups */
        for (int e = 0; e < _num_cmfd_groups; e++)
        {

          /* Load tallies at this cell and energy group */
          double vol_tally = _volume_tally[j][e];
          double rxn_tally = _reaction_tally[j][e];
          _old_flux->setValue(j, e, rxn_tally / vol_tally);

          /* Set the Mesh cell properties with the tallies */
          _volumes->setValue(j, 0, vol_tally);
          log::fdebug("CMFD Cell %d energy %d Old flux is %.2f, volume is %.2f",
                      j, e, rxn_tally / vol_tally, vol_tally);
        }
      }
    }

    /* Loop over boundary CMFD cells and set cross sections */
    if (mpi::isSpatialDecomposed())
    {
#pragma omp parallel for
      for (int s = 0; s < HEX_NUM_FACES; s++)
      {

        /* Loop over all CMFD cells on the current surface */
        std::map<int, int>::iterator it;
        for (it = _boundary_index_map.at(s).begin();
             it != _boundary_index_map.at(s).end(); ++it)
        {

          int idx = it->second;

          /* Loop over CMFD coarse energy groups */
          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            /* Load tallies at this cell and energy group */
            double vol_tally = _boundary_volumes[s][idx][0];
            double rxn_tally = _boundary_reaction[s][idx][e];
            _old_boundary_flux[s][idx][e] = rxn_tally / vol_tally;
          }
        }
      }
    }
  }

  /**
   * @brief Computes the diffusion coefficient for a given CMFD cell and CMFD
   *        energy group.
   * @details This method computes the diffusion coefficient for a CMFD cell and
   *          CMFD energy group by spatially collapsing the total/transport xs
   *          in each FSR contained within the CMFD cell and then energy
   *          collapsing the diffusion coefficient (\f$1 / (3 * \Sigma_t)\f$) for
   *          all MOC groups in the given CMFD energy group.
   * @param cmfd_cell A CMFD cell
   * @param group A CMFD energy group
   * @return The diffusion coefficient
   */
  CMFD_PRECISION Cmfd::getDiffusionCoefficient(int cmfd_cell, int group)
  {
    return _diffusion_tally[cmfd_cell][group] /
           _reaction_tally[cmfd_cell][group];
  }

  /**
   * @brief Compute the surface diffusion coefficient for a given CMFD cell,
   *        cell surface, and group.
   * @details This method uses finite differencing to compute the surface
   *          diffusion coefficient (\f$ \hat{D} \f$) or surface diffusion
   *          coefficient correction (\f$ \tilde{D} \f$) for a given CMFD cell,
   *          cell surface, and CMFD energy group. If the MOC iteration is zero,
   *          (\f$ \tilde{D} \f$) is returned as zero. Since (\f$ \hat{D} \f$) and
   *          (\f$ \tilde{D} \f$) are dependent on each other, they must be
   *          computed together; therefore, the boolean correction is used to
   *          indicate which value is to be returned.
   * @param cmfd_cell A CMFD cell
   * @param surface A surface of the CMFD cell
   * @param group A CMFD energy group
   * @param moc_iteration MOC iteration number
   * @param correction Boolean indicating whether (\f$ \hat{D} \f$) or
   *                   (\f$ \tilde{D} \f$) is to be returned
   * @return The surface diffusion coefficient, (\f$ \hat{D} \f$) or
   *         (\f$ \tilde{D} \f$)
   */
  CMFD_PRECISION Cmfd::getSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                                      int group, int moc_iteration,
                                                      bool correction)
  {

    CMFD_PRECISION dif_surf, dif_surf_corr;
    FP_PRECISION current, current_out, current_in; // 存储当前、流出的电流和流入的电流
    CMFD_PRECISION flux_next;                      // 相邻单元的通量

    /* Get diffusivity and flux for Mesh cell */
    CMFD_PRECISION dif_coef = getDiffusionCoefficient(cmfd_cell, group); // 当前单元的扩散系数
    int global_cmfd_cell = getGlobalCMFDCell(cmfd_cell);
    int global_cmfd_cell_next = getCellNext(global_cmfd_cell, surface);
    CMFD_PRECISION flux = _old_flux->getValue(cmfd_cell, group);
    CMFD_PRECISION delta_interface = getSurfaceWidth(surface);    // 得到该面的面积
    CMFD_PRECISION delta = getPerpendicularSurfaceWidth(surface); // 垂直于该面的CMFD边长长度
    int sense = getSense(surface);

    /* Correct the diffusion coefficient with Larsen's effective diffusion
     * coefficient correction factor 如果没有线性源 (_linear_source 为 false)，则根据 Larsen 有效扩散系数修正因子调整扩散系数。*/
    if (!_linear_source) // 是否使用线性源近似
      dif_coef *= computeLarsensEDCFactor(dif_coef, delta);

    /* If surface is on a boundary with REFLECTIVE or VACUUM BCs, choose
     * approipriate BC */
    if (global_cmfd_cell_next == -1)
    { // 检查是否到达边界，即相邻单元是否不存在

      /* REFLECTIVE BC 表面是反射边界 (REFLECTIVE)，表面扩散系数和表面修正系数都设为 0*/
      if (_boundaries[surface] == REFLECTIVE)
      {
        dif_surf = 0.0;
        dif_surf_corr = 0.0;
      }

      /* VACUUM BC 如果是真空边界 (VACUUM)，*/
      else if (_boundaries[surface] == VACUUM)
      {

        /* Compute the surface-averaged current leaving the cell计算离开单元的表面平均电流*/
        current_out = sense * _surface_currents->getValue // sense 用于确定电流的方向
                              (cmfd_cell, surface * _num_cmfd_groups + group) /
                      delta_interface;

        /* Set the surface diffusion coefficient and MOC correction 计算真空边界的表面扩散系数和表面修正扩散系数*/
        dif_surf = 2 * dif_coef / delta / (1 + 4 * dif_coef / delta);
        dif_surf_corr = (sense * dif_surf * flux - current_out) / flux;

        /* Weight the old and new corrected diffusion coefficients by the
           relaxation factor 如果 dif_surf_corr 已存在旧的修正值，则通过松弛因子 _relaxation_factor 混合旧的和新的修正值*/
        if (_old_dif_surf_valid)
        {
          CMFD_PRECISION old_dif_surf_corr = _old_dif_surf_corr->getValue(cmfd_cell, surface * _num_cmfd_groups + group);
          dif_surf_corr = _relaxation_factor * dif_surf_corr +
                          (1.0 - _relaxation_factor) * old_dif_surf_corr;
        }
      }
    }

    /* If surface is an interface or PERIODIC BC, use finite differencing 如果表面是界面或周期性边界条件，则使用有限差分*/
    else
    {

      /* Get the surface index for the surface in the neighboring cell 获取相邻单元的表面索引*/
      int surface_next = (surface + NUM_FACES / 2) % NUM_FACES; // 比如0->3;3->0

      /* Get the outward current on surface 获取流出表面的电流*/
      current_out = _surface_currents->getValue(cmfd_cell, surface * _num_cmfd_groups + group);

      /* Set diffusion coefficient and flux for the neighboring cell */
      int cmfd_cell_next = getLocalCMFDCell(global_cmfd_cell_next); // 获取相邻单元的域局部一维索引
      CMFD_PRECISION dif_coef_next;
      if (cmfd_cell_next == -1)
      { // 如果 cmfd_cell_next == -1，说明相邻单元不在当前域中

        /* Get the currents in cells touching this boundary 获取边界表面电流*/
        CMFD_PRECISION **boundary_currents = _boundary_surface_currents[surface];

        int idx = _boundary_index_map.at(surface)[global_cmfd_cell_next]; // 获取相邻单元在边界上的索引
        dif_coef_next = _boundary_diffusion[surface][idx][group] /
                        _boundary_reaction[surface][idx][group]; // 计算相邻单元的扩散系数
        flux_next = _old_boundary_flux[surface][idx][group];     // 获取相邻单元的通量

        /* Get the inward current on the surface 相邻单元流入当前表面的电流*/
        current_in = boundary_currents[idx][surface_next * _num_cmfd_groups + group];
      }
      else
      { // 如果相邻单元在当前域中
        // 获取相邻单元的扩散系数、通量和电流
        dif_coef_next = getDiffusionCoefficient(cmfd_cell_next, group);
        flux_next = _old_flux->getValue(cmfd_cell_next, group);

        /* Get the inward current on the surface */
        current_in = _surface_currents->getValue(cmfd_cell_next, surface_next * _num_cmfd_groups + group);
      }

      /* Correct the diffusion coefficient with Larsen's effective diffusion
       * coefficient correction factor */
      if (!_linear_source)
        dif_coef_next *= computeLarsensEDCFactor(dif_coef_next, delta); // 修正相邻单元的扩散系数（如果没有线性源近似）

      /* Compute the surface diffusion coefficient 计算表面扩散系数*/
      dif_surf = 2.0 * dif_coef * dif_coef_next / (delta * dif_coef + delta * dif_coef_next); // 与B.23公式对应,其中delta为δrh的两倍,所以分子多乘了个2

      /* Compute the surface-averaged net current across the surface 计算流经表面的净电流*/
      current = sense * (current_out - current_in) / delta_interface; // 除以delta_interface得的是平均值

      /* Compute the surface diffusion coefficient correction 计算修正的表面扩散系数*/
      dif_surf_corr = -(sense * dif_surf * (flux_next - flux) + current) / (flux_next + flux); // 对应公式B.26,这个current好像多乘了个方向sense(不一定,有可能净流项的计算本身需要乘以sense方向)

      /* Flux limiting condition 通量限制条件*/
      if (_flux_limiting && moc_iteration > 0)
      { // 如果启用了通量限制 (_flux_limiting) 并且迭代次数大于零，检查修正系数与原始系数的比率
        double ratio = dif_surf_corr / dif_surf;
        if (std::abs(ratio) > 1.0)
        { // 如果比率大于 1.0，表示修正过大，需要限制,根据电流和通量重新计算 dif_surf 和 dif_surf_corr。
          // 调整扩散系数
          if (sense * current > 0.0)
            dif_surf = std::abs(current / (2.0 * flux));
          else
            dif_surf = std::abs(current / (2.0 * flux_next));
          // 调整修正扩散系数
          dif_surf_corr = -(sense * dif_surf * (flux_next - flux) + current) / (flux_next + flux);
        }
      }

      /* Weight the old and new corrected diffusion coefficients by the
         relaxation factor 如果已有旧的修正值，使用松弛因子 _relaxation_factor 结合新旧值*/
      if (_old_dif_surf_valid)
      {
        CMFD_PRECISION old_dif_surf_corr = _old_dif_surf_corr->getValue(cmfd_cell, surface * _num_cmfd_groups + group);
        dif_surf_corr = _relaxation_factor * dif_surf_corr +
                        (1.0 - _relaxation_factor) * old_dif_surf_corr;
      }
    }

    /* If it is the first MOC iteration, solve the straight diffusion problem
     * with no MOC correction 如果是第一次MOC迭代 (moc_iteration == 0)，直接将修正系数设为 0.0
     在第一轮迭代中使用一个未修正的扩散系数 dif_surf 可以避免在早期迭代中引入不准确的修正
     在第一轮迭代中不使用修正项，使得初始迭代仅基于基础的扩散系数进行计算。
     这提供了一个干净的基准，从而在后续迭代中可以逐步引入和调整修正系数，使其更加精确.
     */
    if (moc_iteration == 0)
      dif_surf_corr = 0.0;

    /* Determine which surface diffusion coefficient is corrected */
    if (correction) // 返回最终的表面扩散系数或修正扩散系数
      return dif_surf_corr;
    else
      return dif_surf;
  }

  /**
   * @brief Compute the surface diffusion coefficient for a given CMFD cell,
   *        cell surface, and group.
   * @details This method uses finite differencing to compute the surface
   *          diffusion coefficient (\f$ \hat{D} \f$) or surface diffusion
   *          coefficient correction (\f$ \tilde{D} \f$) for a given CMFD cell,
   *          cell surface, and CMFD energy group. If the MOC iteration is zero,
   *          (\f$ \tilde{D} \f$) is returned as zero. Since (\f$ \hat{D} \f$) and
   *          (\f$ \tilde{D} \f$) are dependent on each other, they must be
   *          computed together; therefore, the boolean correction is used to
   *          indicate which value is to be returned.
   * @param cmfd_cell A CMFD cell
   * @param surface A surface of the CMFD cell
   * @param group A CMFD energy group
   * @param moc_iteration MOC iteration number
   * @param correction Boolean indicating whether (\f$ \hat{D} \f$) or
   *                   (\f$ \tilde{D} \f$) is to be returned
   * @return The surface diffusion coefficient, (\f$ \hat{D} \f$) or
   *         (\f$ \tilde{D} \f$)
   */
  CMFD_PRECISION Cmfd::getHexSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                                         int group, int moc_iteration,
                                                         bool correction)
  {
    int logical_cmfd_cell = _logical_actual_map[cmfd_cell];

    CMFD_PRECISION dif_surf, dif_surf_corr;
    FP_PRECISION current, current_out, current_in;
    CMFD_PRECISION flux_next;

    /* Get diffusivity and flux for Mesh cell */
    CMFD_PRECISION dif_coef = getDiffusionCoefficient(logical_cmfd_cell, group); // Cell内扩散系数计算
    int global_cmfd_cell = getGlobalCMFDCell(cmfd_cell);
    int global_cmfd_cell_next = getCellNext(global_cmfd_cell, surface);
    CMFD_PRECISION flux = _old_flux->getValue(logical_cmfd_cell, group);
    CMFD_PRECISION delta_interface = getHexSurfaceWidth(surface);    // surface对应面的面积
    CMFD_PRECISION delta = getHexPerpendicularSurfaceWidth(surface); // 垂直于surface边的长度
    int sense = getSense(surface);                                   // current的朝向

    /* Correct the diffusion coefficient with Larsen's effective diffusion
     * coefficient correction factor */
    if (!_linear_source)
      dif_coef *= computeLarsensEDCFactor(dif_coef, delta);

    /* If surface is on a boundary with REFLECTIVE or VACUUM BCs, choose
     * approipriate BC */
    /* 对于边界处Cell的Current处理方式与位于内部处Cell的Current处理方式有所不同 */
    if (CellinXYZBoundary(global_cmfd_cell, surface) || global_cmfd_cell_next == -1)
    {

      /* REFLECTIVE BC */
      if (_boundaries[surface] == REFLECTIVE)
      {
        dif_surf = 0.0;
        dif_surf_corr = 0.0;
      }

      /* VACUUM BC */
      else if (_boundaries[surface] == VACUUM)
      {

        /* Compute the surface-averaged current leaving the cell */
        current_out = sense * _surface_currents->getValue(logical_cmfd_cell, surface * _num_cmfd_groups + group) / delta_interface;

        /* Set the surface diffusion coefficient and MOC correction */
        dif_surf = 2 * dif_coef / delta / (1 + 4 * dif_coef / delta);
        dif_surf_corr = (sense * dif_surf * flux - current_out) / flux;

        /* Weight the old and new corrected diffusion coefficients by the
           relaxation factor */
        if (_old_dif_surf_valid)
        {
          CMFD_PRECISION old_dif_surf_corr = _old_dif_surf_corr->getValue(logical_cmfd_cell, surface * _num_cmfd_groups + group);
          dif_surf_corr = _relaxation_factor * dif_surf_corr +
                          (1.0 - _relaxation_factor) * old_dif_surf_corr;
        }
      }
    }

    /* If surface is an interface or PERIODIC BC, use finite differencing */
    else
    {

      /* Get the surface index for the surface in the neighboring cell */
      int surface_next = (surface + HEX_NUM_FACES / 2) % HEX_NUM_FACES;

      /* Get the outward current on surface */
      current_out = _surface_currents->getValue(logical_cmfd_cell, surface * _num_cmfd_groups + group);

      /* Set diffusion coefficient and flux for the neighboring cell */
      int cmfd_cell_next = getLocalCMFDCell(global_cmfd_cell_next);
      CMFD_PRECISION dif_coef_next;

      /* Find local CMFD next cell in domain
         Current Hex does not support Spatial Decomposed*/
      if (cmfd_cell_next == -1)
      {
        log::finfo("Use Domain Decomposed");
        /* Get the currents in cells touching this boundary */
        CMFD_PRECISION **boundary_currents = _boundary_surface_currents[surface];

        int idx = _boundary_index_map.at(surface)[global_cmfd_cell_next];
        dif_coef_next = _boundary_diffusion[surface][idx][group] /
                        _boundary_reaction[surface][idx][group];
        flux_next = _old_boundary_flux[surface][idx][group];

        /* Get the inward current on the surface */
        current_in = boundary_currents[idx][surface_next * _num_cmfd_groups + group];
      }
      else
      {
        if (!_empty_fsrs_cells[cmfd_cell_next])
        {
          dif_coef_next = getDiffusionCoefficient(_logical_actual_map[cmfd_cell_next], group);
          flux_next = _old_flux->getValue(_logical_actual_map[cmfd_cell_next], group);

          /* Get the inward current on the surface */
          current_in = _surface_currents->getValue(_logical_actual_map[cmfd_cell_next], surface_next * _num_cmfd_groups + group);
        }
      }

      /* Correct the diffusion coefficient with Larsen's effective diffusion
       * coefficient correction factor */
      if (!_linear_source)
        dif_coef_next *= computeLarsensEDCFactor(dif_coef_next, delta);

      /* Compute the surface diffusion coefficient */
      dif_surf = 2.0 * dif_coef * dif_coef_next / (delta * dif_coef + delta * dif_coef_next);

      /* Compute the surface-averaged net current across the surface */
      current = sense * (current_out - current_in) / delta_interface;

      /* Compute the surface diffusion coefficient correction */
      dif_surf_corr = -(sense * dif_surf * (flux_next - flux) + current) / (flux_next + flux);

      /* Flux limiting condition */
      if (_flux_limiting && moc_iteration > 0)
      {
        double ratio = dif_surf_corr / dif_surf;
        if (std::abs(ratio) > 1.0)
        {

          if (sense * current > 0.0)
            dif_surf = std::abs(current / (2.0 * flux));
          else
            dif_surf = std::abs(current / (2.0 * flux_next));

          dif_surf_corr = -(sense * dif_surf * (flux_next - flux) + current) / (flux_next + flux);
        }
      }

      /* Weight the old and new corrected diffusion coefficients by the
         relaxation factor */
      if (_old_dif_surf_valid)
      {
        CMFD_PRECISION old_dif_surf_corr = _old_dif_surf_corr->getValue(logical_cmfd_cell, surface * _num_cmfd_groups + group);
        dif_surf_corr = _relaxation_factor * dif_surf_corr +
                        (1.0 - _relaxation_factor) * old_dif_surf_corr;
      }
    }

    /* If it is the first MOC iteration, solve the straight diffusion problem
     * with no MOC correction */
    if (moc_iteration == 0)
      dif_surf_corr = 0.0;

    /* Determine which surface diffusion coefficient is corrected */
    if (correction)
      return dif_surf_corr;
    else
      return dif_surf;
  }

  /**
   * @brief Solve the nonlinear diffusion acceleration problem to accelerate the
   *        convergence of the MOC problem.求解非线性扩散加速问题，以加速 MOC（特征线法）问题的收敛
   * @details This method uses the information from the last MOC transport sweep
   *          and solves a simplified nonlinear diffusion problem. The diffusion
   *          problem is tightly converged and the solution is used to update the
   *          the solution of the MOC problem.
   * 该方法利用上一次 MOC 传输扫描的信息来求解简化的非线性扩散问题
   * 扩散问题一定要收敛，其解将用于更新 MOC 问题的解
   *  @param moc_iteration MOC iteration number
   *  @return The dominant eigenvalue of the nonlinear diffusion problem
   */
  double Cmfd::computeKeff(int moc_iteration)
  {

    log::fverbose("Running diffusion solver...");

    /* Start recording total CMFD time */
    _timer->startTimer();

    /* Create matrix and vector objects */
    if (_A == NULL)
    {
      log::ferror("Unable to compute k-eff in CMFD since the CMFD "
                  "linear algebra matrices and arrays have not been created.");
    }

    /* Start recording XS collapse time */
    _timer->startTimer();

    // printHex();

    /* Copy surface currents if neutron balance check requested 暂时不看*/
    if (_check_neutron_balance)
    {
      if (_hexlattice_enable)
        copyHexFullSurfaceCurrents();
      else
        copyFullSurfaceCurrents();
    }

    /*
     * 限制操作（Restriction）
     * 作用：将MOC细网格的FSR通量、截面等数据"压缩"（体积加权平均）到CMFD粗网格
     */
    /* Collapse the cross sections onto the CMFD mesh */
    if (_hexlattice_enable)
      hexCollapseXS();
    else
      collapseXS(); // 就是将 MOC 细网格数据压缩到 CMFD 粗网格上

    /* Tally the XS collpase time */
    _timer->stopTimer("Total collapse time");

    /* ========================================
     * 构建CMFD矩阵
     * 作用：根据压缩后的截面数据，构建CMFD类扩散方程的矩阵 A和M
     * 矩阵A：包含损失项（吸收+散射损失）和流动项（扩散项）
     * 矩阵M：包含裂变源项
     */
    /* Construct matrices */
    if (_hexlattice_enable)
      hexConstructMatrices(moc_iteration);
    else
      constructMatrices(moc_iteration); // 构建CMFD矩阵

    /* Check neturon balance if requested 暂时不看*/
    if (_check_neutron_balance)
      if (_hexlattice_enable)
        checkNeutronBalanceHex();
      else
        checkNeutronBalance();

    /* Copy old flux to new flux 把_old_flux复制到_new_flux中*/
    _old_flux->copyTo(_new_flux);

    //_old_flux->printString();

    /* Start recording CMFD solve time */
    _timer->startTimer();

    //_A->printString();
    //_M->printString();
    double k_eff;
    /* Solve the eigenvalue problem 求解CMFD特征值问题*/
    if (_hexlattice_enable)
    {
      k_eff = hexeigenvalueSolve(_A, _M, _new_flux, _k_eff,
                                 _source_convergence_threshold, _SOR_factor,
                                 _convergence_data, _domain_communicator, _empty_cells_num,
                                 _empty_fsrs_cells, _logical_actual_map); // 00000
    }
    else
    {
      k_eff = eigenvalueSolve(_A, _M, _new_flux, _k_eff,
                              _source_convergence_threshold, _SOR_factor,
                              _convergence_data, _domain_communicator); // 使用幂迭代方法（Power Iteration Method）来求解广义特征值问题的算法
    }

    /* Try to use a few-group solver to remedy convergence issues 使用少群求解器解决收敛问题*/
    bool reduced_group_solution = false;
    // 如果 k_eff 的计算结果不收敛，并且群数多于备份群数，尝试使用少群求解器
    if (fabs(k_eff + 1) < FLT_EPSILON && _num_cmfd_groups > _num_backup_groups)
    {

      log::finfo("Switching to a %d group CMFD solver on this iteration",
                 _num_backup_groups);

      if (_backup_cmfd == NULL) // 备份 CMFD 求解器为空，初始化它
        initializeBackupCmfdSolver();

      copyCurrentsToBackup();                           // 当前的中字流复制到备份求解器
      k_eff = _backup_cmfd->computeKeff(moc_iteration); // 使用少群求解器重新计算 k_eff
      reduced_group_solution = true;                    // 标记使用了少群求解器
    }

    /* Tally the CMFD solver time */
    _timer->stopTimer("Total solver time");

    /* Check for a legitimate solve */
    if (fabs(k_eff + 1) > FLT_EPSILON)
      _k_eff = k_eff;
    else
      return _k_eff; // 如果无效，直接返回之前的 k_eff

    /* Do not prolong again if the few-group solution was used */
    if (reduced_group_solution)
    { // 如果使用少群求解器
      log::finfo("The %d group CMFD solver was successful",
                 _num_backup_groups);
      return _k_eff;
    }

    /* Rescale the old and new flux 重新标定旧的和新的通量*/
    rescaleFlux(); // 归一化操作,将初始和收敛后的通量进行重缩放，以确保每个网格单元和每个能群中的平均裂变源项都为标准化的值1.0

    /*
     * 延拓操作（Prolongation）
     * 作用：将收敛后的CMFD粗网格标量通量插值回MOC细网格，修正FSR通量
     *
     * 关键点：
     * - 计算更新比率：ratio = φ_CMFD_new / φ_CMFD_old
     * - 更新FSR通量：φ_FSR_new = φ_FSR_old * ratio
     * - 这样可以为下一次MOC迭代提供更好的源项
     */
    /* Update the MOC flux */
    if (_hexlattice_enable)
      updateHexMOCFlux(); // 更新六边形 MOC 通量
    else
      updateMOCFlux();

    /* Tally the total CMFD time */
    _timer->stopTimer("Total CMFD time");

    return _k_eff;
  }

  /**
   * @
   * brief Rescale the initial and converged flux arrays.
   * @details The diffusion problem is a generalized eigenvalue problem and
   *          therefore the solution is independent of flux level. This method
   *          rescales the input flux and converged flux to both have an average
   *          fission source of 1.0 in each group in each cell.
   */
  void Cmfd::rescaleFlux()
  {

    /* Rescale the new and old flux to have an avg source of 1.0 */
    matrixMultiplication(_M, _new_flux, _new_source);
    matrixMultiplication(_M, _old_flux, _old_source);

    double new_source_sum = _new_source->getSum();
    double old_source_sum = _old_source->getSum();
#ifdef ENABLE_MPI_
    if (_domain_communicator != NULL)
    {
      double temp_sum_new = new_source_sum;
      MPI_Allreduce(&temp_sum_new, &new_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                    _domain_communicator->_MPI_cart);
      double temp_sum_old = old_source_sum;
      MPI_Allreduce(&temp_sum_old, &old_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                    _domain_communicator->_MPI_cart);
    }
#endif
    // 归一化操作,将初始和收敛的通量重缩放，使得在每个能群和每个网格单元的平均源项为标准化的值1.0.
    _new_flux->scaleByValue(1.0 / new_source_sum);
    _old_flux->scaleByValue(1.0 / old_source_sum);
  }

  /**
   * @brief Construct the loss + streaming matrix (A) and the fission gain
   *         matrix (M) in preparation for solving the eigenvalue problem.
   * 构建两个矩阵：损失+流动矩阵（_A 矩阵）和 裂变增益矩阵（_M 矩阵），以准备求解特征值问题
   * @details This method loops over all mesh cells and energy groups and
   *          accumulates the iteraction and streaming terms into their
   *          approipriate positions in the loss + streaming matrix and
   *          fission gain matrix.
   */
  void Cmfd::constructMatrices(int moc_iteration)
  {

    log::fverbose("Constructing matrices...");

    /* Zero _A and _M matrices 清空 _A 和 _M 矩阵，确保矩阵在开始时没有任何累积的值*/
    _A->clear();
    _M->clear();

    /* Zero the number of connections */
    if (_domain_communicator != NULL)
    {
      int num_local_cells = _local_num_x * _local_num_y * _local_num_z;
      for (int c = 0; c < 2; c++)
      {
        for (int ncg = 0; ncg < num_local_cells * _num_cmfd_groups; ncg++)
        {
          _domain_communicator->num_connections[c][ncg] = 0;
        }
      }
    }
#pragma omp parallel
    {

      FP_PRECISION value, volume, delta;
      CMFD_PRECISION dif_surf, dif_surf_corr;
      int sense;
      Material *material;

      // Unused, FIXME
      // int x_start = 0;
      // int y_start = 0;
      // int z_start = 0;
      // int x_end = _num_x;
      // int y_end = _num_y;
      // int z_end = _num_z;
      // if (mpi::isSpatialDecomposed()) {
      //   if (_domain_communicator != NULL) {
      //     x_start = _domain_communicator->_domain_idx_x * _local_num_x;
      //     x_end = x_start + _local_num_x;
      //     y_start = _domain_communicator->_domain_idx_y * _local_num_y;
      //     y_end = y_start + _local_num_y;
      //     z_start = _domain_communicator->_domain_idx_z * _local_num_z;
      //     z_end = z_start + _local_num_z;
      //   }
      // }

      /* Loop over cells */
#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        int global_ind = getGlobalCMFDCell(i); // 得到域全局一维索引
        /*结果为 0，则颜色为黑色；如果结果为 1，则颜色为红色。
        这个方法确保在三维空间中相邻单元格的颜色是不同的，有利于红/黑 SOR 算法的执行
         */
        int color = getCellColor(global_ind); // color = (ix + iy + iz) % 2
        // 获取其材料属性和体积
        material = _materials[i];
        volume = _volumes->getValue(i, 0);
        log::fdebug("CMFD Cell %d volume is %.6f", i, volume);

        /* Loop over groups */
        for (int e = 0; e < _num_cmfd_groups; e++)
        {

          /* ========================================
           * 构建矩阵A的对角线元素 - 净去除项
           * ========================================
           * 原理对应：文档B.28公式的最左侧项
           * 物理意义：该CMFD网格内的中子损失（吸收+散射出群）
           */
          /* Net removal term 净去除项*/
          value = material->getSigmaTByGroup(e + 1) * volume; // 总截面（SigmaT），乘以单元体积  B.28最左侧
          _A->incrementValue(i, e, i, e, value);              // 行和列都一样,说明这个值加到了对角线元素上  将此值累加到 _A 矩阵的对角线元素中
          log::fdebug("Get the Material's total cross section for energy group %d. value is %.6f", e + 1, value);

          /* Re-compute diagonal if neutron re-balance requested 如果要求重新平衡中子，则重新计算对角线元素*/
          if (_balance_sigma_t)
          {
            enforceBalanceOnDiagonal(i, e); // 注:不是很明白这个函数
          }

          /* ========================================
           * 构建矩阵A的非对角线元素 - 散射增益项
           * ========================================
           * 原理对应：B.28公式左侧的散射项
           * 物理意义：从其他能群g散射到当前能群e的中子增益
           */
          /* Scattering gain from all groups  散射增益项*/
          for (int g = 0; g < _num_cmfd_groups; g++)
          {                                                             // 这表示从 g 能群向当前 e 能群的散射贡献
            value = -material->getSigmaSByGroup(g + 1, e + 1) * volume; // 对应B.28的散射截面*体积,左侧最后一项
            _A->incrementValue(i, g, i, e, value);                      // 反映在矩阵中相应的非对角元素上
            log::fdebug("Get the Material's scattering cross section for energy group from %d to %d. value is %.6f, ", g + 1, e + 1, value);
          }

          /* ========================================
           * 构建矩阵A的流动项 - 扩散项
           * ========================================
           * 原理对应：B.28公式的流动项（含扩散系数修正）
           * 物理意义：中子从当前网格流向邻居网格的扩散流动
           *
           * 关键点：
           * - dif_surf：表面扩散系数（基于体积平均截面）
           * - dif_surf_corr：表面修正扩散系数（保证表面净电流守恒）
           * - sense：方向因子，决定流动的方向性
           */
          /* Streaming to neighboring cells
          遍历所有的网格单元面，计算网格单元的流动项，包括表面扩散系数、修正的扩散系数以及它们在对角线和非对角线元素中的贡献*/
          for (int s = 0; s < NUM_FACES; s++)
          {

            sense = getSense(s);        // 获取表面的方向，决定该表面与邻居单元的关系 对应u(j,h)
            delta = getSurfaceWidth(s); // 获取该面的面积  对应A

            /* Set transport term on diagonal */
            dif_surf = getSurfaceDiffusionCoefficient(
                i, s, e, moc_iteration, false); // 表面扩散系数获取
            dif_surf_corr = getSurfaceDiffusionCoefficient(
                i, s, e, moc_iteration, true); // 表面修正扩散系数获取
            log::fdebug("Set transport term on diagonal by cell:%d surface:%d energy:%d, dif_surf is %.6f, dif_surf_corr is %.6f", i, s, e, dif_surf, dif_surf_corr);

            /* Record the corrected diffusion coefficient 记录表面修正扩散系数*/
            _old_dif_surf_corr->setValue(i, s * _num_cmfd_groups + e, dif_surf_corr);
            _old_dif_surf_valid = true; // 表明修正的扩散系数已经被计算并存储

            /* Set the diagonal term 设置对角线值*/
            value = (dif_surf - sense * dif_surf_corr) * delta; // 对应B.28第一项的第二项吗?
            _A->incrementValue(i, e, i, e, value);
            log::fdebug("Set the diagonal term for energy group %d, value is %.6f", e + 1, value);

            /* Set the off diagonal term 设置非对角线值*/
            if (getCellNext(i, s, false, false) != -1)
            {                                                      // 检查当前表面 s 是否有邻居单元格
              value = -(dif_surf + sense * dif_surf_corr) * delta; // 对应B.28第二项?这个sense怎么感觉有点问题
              _A->incrementValue(getCellNext(i, s, false, false), e, i, e, value);
            }

            /* Check for cell in neighboring domain if applicable */
            else if (mpi::isSpatialDecomposed())
            {
              if (_domain_communicator != NULL)
              {
                if (getCellNext(i, s, false, true) != -1)
                {
                  int neighbor_cell = getCellNext(i, s, false, true);
                  int row = i * _num_cmfd_groups + e;
                  int idx = _domain_communicator->num_connections[color][row];
                  value = -(dif_surf + sense * dif_surf_corr) * delta;
                  _domain_communicator->indexes[color][row][idx] = neighbor_cell;
                  _domain_communicator->domains[color][row][idx] = s;
                  _domain_communicator->coupling_coeffs[color][row][idx] = value;
                  _domain_communicator->num_connections[color][row]++;
                }
              }
            }
          }

          /* Fission source term 构建矩阵M - 裂变源项*/
          for (int g = 0; g < _num_cmfd_groups; g++)
          {
            value = material->getChiByGroup(e + 1) * material->getNuSigmaFByGroup(g + 1) * volume; // 裂变能谱*裂变截面*平均裂变中子数*体积 B.28的右侧项
            _M->incrementValue(i, g, i, e, value);
            log::fdebug("Fission source term cell from %d to %d value is %.6f", g + 1, e + 1, value);
          }
        }
      }
    }

    log::fverbose("Done constructing matrices...");
  }

  /**
   * @brief Construct the loss + streaming matrix (A) and the fission gain
   *         matrix (M) in preparation for solving the eigenvalue problem.
   * @details This method loops over all mesh cells and energy groups and
   *          accumulates the iteraction and streaming terms into their
   *          approipriate positions in the loss + streaming matrix and
   *          fission gain matrix.
   */
  void Cmfd::hexConstructMatrices(int moc_iteration)
  {

    log::fverbose("Constructing matrices...");

    /* Zero _A and _M matrices */
    _A->clear();
    _M->clear();

    /* Zero the number of connections */
    if (_domain_communicator != NULL)
    {
      int num_local_cells = _local_num_x * _local_num_y * _local_num_z;
      for (int c = 0; c < 2; c++)
      {
        for (int ncg = 0; ncg < num_local_cells * _num_cmfd_groups; ncg++)
        {
          _domain_communicator->num_connections[c][ncg] = 0;
        }
      }
    }
#pragma omp parallel
    {

      FP_PRECISION value, volume, delta;
      CMFD_PRECISION dif_surf, dif_surf_corr;
      int sense;
      Material *material;

      // Unused, FIXME
      // int x_start = 0;
      // int y_start = 0;
      // int z_start = 0;
      // int x_end = _num_x;
      // int y_end = _num_y;
      // int z_end = _num_z;
      // if (mpi::isSpatialDecomposed()) {
      //   if (_domain_communicator != NULL) {
      //     x_start = _domain_communicator->_domain_idx_x * _local_num_x;
      //     x_end = x_start + _local_num_x;
      //     y_start = _domain_communicator->_domain_idx_y * _local_num_y;
      //     y_end = y_start + _local_num_y;
      //     z_start = _domain_communicator->_domain_idx_z * _local_num_z;
      //     z_end = z_start + _local_num_z;
      //   }
      // }

      /* Loop over cells */
#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        if (_empty_fsrs_cells[i])
          continue;
        else
        {
          int j = _logical_actual_map[i];
          log::fdebug("Loop over with cell %d", j);
          // int global_ind = getGlobalCMFDCell(i);
          int color = getCellColor(i);

          material = _materials[j];
          volume = _volumes->getValue(j, 0);

          /* Loop over groups */
          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            /* Net removal term */
            value = material->getSigmaTByGroup(e + 1) * volume;
            _A->incrementValue(j, e, j, e, value);
            log::fdebug("Get the Material's total cross section for energy group %d. value is %.6f", e + 1, value);

            /* Re-compute diagonal if neutron re-balance requested */
            if (_balance_sigma_t)
            {
              enforceBalanceOnDiagonal(j, e);
            }

            /* Scattering gain from all groups */
            for (int g = 0; g < _num_cmfd_groups; g++)
            {
              value = -material->getSigmaSByGroup(g + 1, e + 1) * volume;
              _A->incrementValue(j, g, j, e, value);
              log::fdebug("Get the Material's scattering cross section for energy group from %d to %d. value is %.6f, ", g + 1, e + 1, value);
            }

            /* Streaming to neighboring cells */
            for (int s = 0; s < HEX_NUM_FACES; s++)
            {

              sense = getSense(s);
              delta = getHexSurfaceWidth(s);

              /* Set transport term on diagonal */
              dif_surf = getHexSurfaceDiffusionCoefficient(
                  i, s, e, moc_iteration, false);
              dif_surf_corr = getHexSurfaceDiffusionCoefficient(
                  i, s, e, moc_iteration, true);
              log::fdebug("Set transport term on diagonal by cell:%d surface:%d energy:%d, dif_surf is %.6f, dif_surf_corr is %.6f", j, s, e, dif_surf, dif_surf_corr);

              /* Record the corrected diffusion coefficient */
              _old_dif_surf_corr->setValue(j, s * _num_cmfd_groups + e, dif_surf_corr);
              _old_dif_surf_valid = true;

              /* Set the diagonal term */
              value = (dif_surf - sense * dif_surf_corr) * delta;
              _A->incrementValue(j, e, j, e, value);
              log::fdebug("Set the diagonal term for energy group %d, value is %.6f", e + 1, value);

              /* Set the off diagonal term */
              if ((getCellNext(i, s, false, false) != -1) && !CellinXYZBoundary(i, s))
              {
                value = -(dif_surf + sense * dif_surf_corr) * delta;
                if (!_empty_fsrs_cells[getCellNext(i, s, false, false)])
                  _A->incrementValue(_logical_actual_map[getCellNext(i, s, false, false)], e, j, e, value);
              }

              /* Check for cell in neighboring domain if applicable */
              else if (mpi::isSpatialDecomposed())
              {
                if (_domain_communicator != NULL)
                {
                  if (getCellNext(i, s, false, true) != -1)
                  {
                    int neighbor_cell = getCellNext(i, s, false, true);
                    int row = i * _num_cmfd_groups + e;
                    int idx = _domain_communicator->num_connections[color][row];
                    value = -(dif_surf + sense * dif_surf_corr) * delta;
                    _domain_communicator->indexes[color][row][idx] = neighbor_cell;
                    _domain_communicator->domains[color][row][idx] = s;
                    _domain_communicator->coupling_coeffs[color][row][idx] = value;
                    _domain_communicator->num_connections[color][row]++;
                  }
                }
              }
            }

            /* Fission source term */
            for (int g = 0; g < _num_cmfd_groups; g++)
            {
              value = material->getChiByGroup(e + 1) * material->getNuSigmaFByGroup(g + 1) * volume;
              _M->incrementValue(j, g, j, e, value);
              log::fdebug("Fission source term cell from %d to %d value is %.6f", g + 1, e + 1, value);
            }
          }
        }
      }
    }

    log::fverbose("Done constructing matrices...");
  }

  /**
   * @brief Update the MOC flux in each FSR.更新每个 FSR的 MOC 通量
   * @details This method uses the condensed flux from the last MOC transport
   *          sweep and the converged flux from the eigenvalue problem to
   *          update the MOC flux in each FSR.
   *
   *  * 原理对应关系：
   * - 对应文档1.1节"延拓操作"：用收敛后的CMFD粗网格标量通量插值修正MOC细网格通量
   * - 这是多重网格方法的第二步，提供更好的源项给下一次MOC迭代
   *
   * 核心思想：
   * 1. 计算更新比率：ratio = φ_CMFD_new / φ_CMFD_old
   * 2. 更新FSR通量：φ_FSR_new = φ_FSR_old * ratio
   * 3. 这样修正后的通量分布更接近收敛解，加速MOC迭代
   */
  void Cmfd::updateMOCFlux()
  {

    log::fverbose("Updating MOC flux...");

    /* Set max prolongation factor 如果 _convergence_data 不为空，则将最大延续因子 pf 设置为 1.0,用于初始化*/
    if (_convergence_data != NULL)
      _convergence_data->pf = 1.0;

    /* ========================================
     * 遍历所有CMFD网格及其包含的FSR
     * ========================================
     */
    /* Loop over mesh cells */
#pragma omp parallel for
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    { // 遍历CMFD网格

      std::vector<long>::iterator iter;

      /* Loop over CMFD groups */
      for (int e = 0; e < _num_cmfd_groups; e++)
      {

        /* Loop over FRSs in mesh cell */
        for (iter = _cell_fsrs.at(i).begin();
             iter != _cell_fsrs.at(i).end(); ++iter)
        {

          /* ========================================
           * 获取更新比率
           * ========================================
           * 原理对应：延拓操作的核心
           * 作用：计算CMFD粗网格新旧通量的比率，用于修正FSR通量
           */
          /* Get the update ratio */
          CMFD_PRECISION update_ratio = getUpdateRatio(i, e, *iter);

          /* Limit the update ratio */
          if (update_ratio > 20.0)
            update_ratio = 20.0;
          if (update_ratio < 0.05)
            update_ratio = 0.05;

          if (_convergence_data != NULL)
            if (std::abs(std::log(update_ratio)) > std::abs(std::log(_convergence_data->pf)))
              _convergence_data->pf = update_ratio;

          /* ========================================
           * 更新FSR通量（延拓操作）
           * ========================================
           * 原理对应：文档1.1节"延拓：用粗网格解插值修正细网格解"
           *
           * 物理意义：
           * - φ_FSR_old：MOC上一次迭代得到的FSR通量
           * - update_ratio：CMFD修正后的通量比率
           * - φ_FSR_new：修正后的FSR通量，用于下一次MOC迭代
           *
           * 关键点：
           * - 这个修正保证了粗网格和细网格解的一致性
           * - 为下一次MOC迭代提供了更好的初始猜测
           */
          for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
          {

            /* Update FSR flux using ratio of old and new CMFD flux */
            _FSR_fluxes[*iter * _num_moc_groups + h] *= update_ratio;

            /* Update flux moments if they were set 如果设置了通量矩，则更新通量矩*/
            if (_linear_source)
            {
              _flux_moments[(*iter) * 3 * _num_moc_groups + h * 3] *= update_ratio;
              _flux_moments[(*iter) * 3 * _num_moc_groups + h * 3 + 1] *= update_ratio;
              _flux_moments[(*iter) * 3 * _num_moc_groups + h * 3 + 2] *= update_ratio;
            }

            log::fdebug("Updating flux in FSR: %d, cell: %d, MOC group: "
                        "%d, CMFD group: %d, ratio: %f",
                        *iter, i, h, e, update_ratio);
          }
        }
      }
    }
#ifdef ENABLE_MPI_
    if (_domain_communicator != NULL && _convergence_data != NULL)
    {
      double max_pf = _convergence_data->pf;
      MPI_Allreduce(&max_pf, &_convergence_data->pf, 1, MPI_DOUBLE, MPI_MAX,
                    _domain_communicator->_MPI_cart);
    }
#endif
  }

  /**
   * @brief Update the Hex Lattice MOC flux in each FSR.
   * @details This method uses the condensed flux from the last MOC transport
   *          sweep and the converged flux from the eigenvalue problem to
   *          update the MOC flux in each FSR.
   */
  void Cmfd::updateHexMOCFlux()
  {

    log::fverbose("Updating MOC flux...");

    /* Set max prolongation factor */
    if (_convergence_data != NULL)
      _convergence_data->pf = 1.0;

    /* Loop over mesh cells */
#pragma omp parallel for
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {
      if (!_empty_fsrs_cells[i])
      {
        std::vector<long>::iterator iter;

        /* Loop over CMFD groups */
        for (int e = 0; e < _num_cmfd_groups; e++)
        {

          /* Loop over FRSs in mesh cell */
          for (iter = _cell_fsrs.at(i).begin();
               iter != _cell_fsrs.at(i).end(); ++iter)
          {

            /* Get the update ratio */
            CMFD_PRECISION update_ratio = getHexUpdateRatio(i, e, *iter);
            // log::finfo("Before cell: %d, energy group: %d, ratio: %.2f", _logical_actual_map[i], e, update_ratio);

            /* Limit the update ratio */
            if (update_ratio > 20.0)
              update_ratio = 20.0;
            if (update_ratio < 0.05)
              update_ratio = 0.05;

            // log::finfo("After cell: %d, energy group: %d, ratio: %.2f", _logical_actual_map[i], e, update_ratio);

            if (_convergence_data != NULL)
              if (std::abs(std::log(update_ratio)) > std::abs(std::log(_convergence_data->pf)))
                _convergence_data->pf = update_ratio;

            for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
            {

              /* Update FSR flux using ratio of old and new CMFD flux */
              _FSR_fluxes[*iter * _num_moc_groups + h] *= update_ratio;

              /* Update flux moments if they were set */
              if (_linear_source)
              {
                _flux_moments[(*iter) * 3 * _num_moc_groups + h * 3] *= update_ratio;
                _flux_moments[(*iter) * 3 * _num_moc_groups + h * 3 + 1] *= update_ratio;
                _flux_moments[(*iter) * 3 * _num_moc_groups + h * 3 + 2] *= update_ratio;
              }
            }
          }
        }
      }
    }
#ifdef ENABLE_MPI_
    if (_domain_communicator != NULL && _convergence_data != NULL)
    {
      double max_pf = _convergence_data->pf;
      MPI_Allreduce(&max_pf, &_convergence_data->pf, 1, MPI_DOUBLE, MPI_MAX,
                    _domain_communicator->_MPI_cart);
    }
#endif
  }

  /**
   * @brief Compute Larsen's effective diffusion coefficient correction factor.
   * @details By conserving reaction and leakage rates within cells, CMFD
   *          guarantees preservation of area-averaged scalar fluxes and net
   *          surface currents from the MOC fixed source iteration if the CMFD
   *          equations can be converged. However, when the MOC mesh cell size
   *          becomes significantly larger than the neutron mean free path in that
   *          cell, the step characteristics no longer preserve the linear
   *          infinite medium solution to the transport equation. While the
   *          surface diffusion coefficient correction term in CMFD is guaranteed
   *          to preserve reaction rates and surface net currents for any choice
   *          of diffusion coefficient, convergence (and convergence rate) of the
   *          nonlinear iteration acceleration of CMFD is affected by the choice
   *          of diffusion coefficient. All flat source methods, when applied for
   *          thick optical meshes, artificially distribute neutrons in space.
   *          This is the reason that Larsen's effective diffusion coefficient is
   *          useful in ensuring that the CMFD acceleration equations have a
   *          diffusion coefficient (on the flux gradient term) that is
   *          consistent, not with the physical transport problem, but with the
   *          transport problem that is being accelerated by the CMFD equations.
   *          Larsen's effective diffusion coefficient is precisely this term in
   *          the one-dimensional limit. The following publications provide
   *          further background on how this term is derived and used:
   *
   *            [1] E. Larsen, "Infinite Medium Solutions to the transport
   *                equation, Sn discretization schemes, and the diffusion
   *                approximation", M&C 2001.
   *            [2] S. Shaner, "Transient Method of Characteristics via the
   *                Adiabatic, Theta, and Multigrid Amplitude Function Methods",
   *                Masters Thesis, MIT 2014.
   * @param dif_coef Diffusion coefficient before applying correction factor
   * @param delta Width of the cell in the direction of interest
   * @return The diffusion coefficient correction factor
   */
  CMFD_PRECISION Cmfd::computeLarsensEDCFactor(CMFD_PRECISION dif_coef,
                                               CMFD_PRECISION delta)
  {

    /* Initialize variables */
    CMFD_PRECISION alpha, mu, expon;
    double rho = 0.0;

    /* Loop over polar angles */
    for (int p = 0; p < _num_polar / 2; p++)
    {
      mu = cos(asin(_quadrature->getSinTheta(0, p)));
      expon = exp(-delta / (3 * dif_coef * mu));
      alpha = (1 + expon) / (1 - expon) - 2 * (3 * dif_coef * mu) / delta;
      rho += 2.0 * mu * _quadrature->getPolarWeight(0, p) * alpha;
    }

    /* Compute the correction factor */
    CMFD_PRECISION correction = 1.0 + delta * rho / (2 * dif_coef);

    return correction;
  }

  /**
   * @brief Set the FSR materials array pointer.
   * @param pointer to FSR_materials array
   */
  void Cmfd::setFSRMaterials(Material **FSR_materials)
  {
    _FSR_materials = FSR_materials;
  }

  /**
   * @brief Set the pointer to the array of FSR_volumes.
   * @param array of FSR volumes
   */
  void Cmfd::setFSRVolumes(FP_PRECISION *FSR_volumes)
  {
    _FSR_volumes = FSR_volumes;
  }

  /**
   * @brief Set pointer to FSR flux array.
   * @param pointer to FSR flux array
   */
  void Cmfd::setFSRFluxes(FP_PRECISION *scalar_flux)
  {
    _FSR_fluxes = scalar_flux;
  }

  /**
   * @brief Set pointer to FSR source array.
   * @param pointer to FSR source array
   */
  void Cmfd::setFSRSources(FP_PRECISION *sources)
  {
    _FSR_sources = sources;
  }

  /**
   * @brief Set pointer to source region flux moments array
   * @param pointer to source region flux moments array
   */
  void Cmfd::setFluxMoments(FP_PRECISION *flux_moments)
  {
    _flux_moments = flux_moments;
    _linear_source = true;
  }

  /**
   * @brief Set successive over-relaxation relaxation factor.
   * @param over-relaxation factor
   */
  void Cmfd::setSORRelaxationFactor(double SOR_factor)
  {

    if (SOR_factor <= 0.0 || SOR_factor >= 2.0)
      log::ferror("The successive over-relaxation relaxation factor "
                  "must be > 0 and < 2. Input value: %i",
                  SOR_factor);

    _SOR_factor = SOR_factor;
  }

  /**
   * @brief Set the CMFD relaxation factor applied to diffusion coefficients
   * @param CMFD relaxation factor
   */
  void Cmfd::setCMFDRelaxationFactor(double relaxation_factor)
  {

    if (relaxation_factor <= 0.0 || relaxation_factor > 1.0)
      log::ferror("The successive over-relaxation relaxation factor "
                  "must be greater than 0 and less than or equal to 1. "
                  "Input value: %i",
                  relaxation_factor);

    _relaxation_factor = relaxation_factor;
    if (fabs(relaxation_factor - 1.0) > FLT_EPSILON)
      log::finfo("CMFD relaxation factor: %6.4f", _relaxation_factor);
  }

  /**
   * @brief Forces CMFD to check neutron balance on every solve
   */
  void Cmfd::checkBalance()
  {
    _check_neutron_balance = true;
  }

  /**
   * @brief Get the number of coarse CMFD energy groups.
   * @return the number of CMFD energy groups
   */
  int Cmfd::getNumCmfdGroups()
  {
    return _num_cmfd_groups;
  }

  /**
   * @brief Set a coarse energy group structure for CMFD.
   * @details CMFD does not necessarily need to have the same energy group
   *          structure as the MOC problem. This function can be used to set
   *          a sparse energy group structure to speed up the CMFD solve. An
   *          example of how this may be called from Python to use a coarse
   *          2-group CMFD structure atop a fine 7-group MOC structure is
   *          illustrated below:
   *
   * @code
   *          cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
   * @endcode
   *
   * @param group_indices A nested vector of MOC-to-CMFD group mapping
   */
  void Cmfd::setGroupStructure(std::vector<std::vector<int>> group_indices)
  {

    _user_group_indices = true;

    /* Delete old group indices array if it exists */
    if (_group_indices != NULL)
      delete[] _group_indices;

    /* Allocate memory for new group indices */
    _num_cmfd_groups = group_indices.size();
    _group_indices = new int[_num_cmfd_groups + 1]; // 数组大小 = CMFD能群数 + 1

    /* Initialize first group index to 0 */
    int last_moc_group = 0; // 设置这个变量的作用是保证CMFD对应的MOC能群的值是单调线性递增的

    /* Set MOC group bounds for rest of CMFD energy groups */
    for (int i = 0; i < _num_cmfd_groups; i++)
    {
      for (size_t j = 0; j < group_indices[i].size(); j++)
      {
        if (group_indices[i][j] <= last_moc_group) // 每个数字必须比前一个大
          log::ferror("The CMFD coarse group indices are not "
                      "monotonically increasing");
        last_moc_group = group_indices[i][j];
      }
      _group_indices[i] = group_indices[i][0] - 1; //_group_indices[0]  = 0   _group_indices[1] = 3
      log::fdebug("CMFD group indices %d: %d", i, _group_indices[i]);
    }

    /*
    _group_indices[0] = 1 - 1 = 0  // CMFD第0群从MOC第0群开始
    _group_indices[1] = 4 - 1 = 3  // CMFD第1群从MOC第3群开始
    _group_indices[2] = 7          // 总共7个MOC能群
       */
    _group_indices[_num_cmfd_groups] =
        group_indices[_num_cmfd_groups - 1].back(); // back()获取最后一个元素
    log::fdebug("CMFD group indices %d: %d",
                _num_cmfd_groups, _group_indices[_num_cmfd_groups]);
  }

  /**
   * @brief Initialize the CMFD materials.
   */
  void Cmfd::initializeMaterials()
  {

    Material *material;

    /* Delete old CMFD materials if it exists */
    if (_materials != NULL)
    {
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
        delete _materials[i];
      delete[] _materials;
    }

    try
    {
      _materials = new Material *[_local_num_x * _local_num_y * _local_num_z];
      for (int z = 0; z < _local_num_z; z++)
      {
        for (int y = 0; y < _local_num_y; y++)
        {
          for (int x = 0; x < _local_num_x; x++)
          {
            int ind = z * _local_num_x * _local_num_y + y * _local_num_x + x;
            material = new Material(ind);
            material->setNumEnergyGroups(_num_cmfd_groups);
            _materials[ind] = material;
          }
        }
      }
    }
    catch (std::exception &e)
    {
      log::ferror("Could not allocate memory for the Mesh cell materials. "
                  "Backtrace:%s",
                  e.what());
    }
  }

  /**
   * @brief Initialize the CMFD materials.
   */
  void Cmfd::initializeHexMaterials()
  {

    Material *material;

    /* Delete old CMFD surface currents vector if it exists */
    if (_materials != NULL)
    {
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z - _empty_cells_num; i++)
        delete _materials[i];
      delete[] _materials;
    }

    try
    {
      _materials = new Material *[_local_num_x * _local_num_y * _local_num_z - _empty_cells_num];
      for (int z = 0; z < _local_num_z; z++)
      {
        for (int y = 0; y < _local_num_y; y++)
        {
          for (int x = 0; x < _local_num_x; x++)
          {
            int ind = z * _local_num_x * _local_num_y + y * _local_num_x + x;
            if (!_empty_fsrs_cells[ind])
            {
              ind = _logical_actual_map[ind];
              material = new Material(ind);
              material->setNumEnergyGroups(_num_cmfd_groups);
              _materials[ind] = material;
            }
          }
        }
      }
    }
    catch (std::exception &e)
    {
      log::ferror("Could not allocate memory for the Mesh cell materials. "
                  "Backtrace:%s",
                  e.what());
    }
  }

  /**
   * @brief Initializes CMFD surface currents Vector prior to first MOC iteration.
   * 在第一次MOC迭代之前初始化CMFD表面中子流向量
   */
  void Cmfd::initializeCurrents()
  {

    /* Delete old CMFD surface currents vector if it exists */
    if (_surface_currents != NULL)
      delete _surface_currents;

    /* Allocate memory for the CMFD Mesh surface and corner currents Vectors  存储每个 CMFD 单元的表面电流*/
    _surface_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                   _local_num_z, _num_cmfd_groups * NUM_FACES);

    if (_balance_sigma_t)
    { // 在每次扫描计算时，根据 MOC 的解重新平衡和调整总截面（sigma-t）的值，以确保计算的一致性和准确性*/
      /* Allocate memory for the actual starting currents on boundary CMFD cells
      存储由起始边界通量产生的总电流*/
      _starting_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                      _local_num_z, _num_cmfd_groups);

      /* Allocate memory for the net currents of all CMFD cells 净电流是指单元内中子流入和流出之间的差值，为出-入*/
      _net_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                 _local_num_z, _num_cmfd_groups);
    }
  }

  void Cmfd::initializeCurrentsHex()
  {

    /* Delete old CMFD surface currents vector if it exists */
    if (_surface_currents != NULL)
      delete _surface_currents;

    /* Allocate memory for the Hex CMFD Mesh surface and corner currents Vectors */
    _surface_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                   _local_num_z, _num_cmfd_groups * HEX_NUM_FACES, _empty_cells_num);

    if (_balance_sigma_t)
    {
      /* Allocate memory for the actual starting currents on boundary CMFD cells */
      _starting_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                      _local_num_z, _num_cmfd_groups, _empty_cells_num);

      /* Allocate memory for the net currents of all CMFD cells */
      _net_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                 _local_num_z, _num_cmfd_groups, _empty_cells_num);
    }
  }

  /**
   * @brief Initializes the vector of vectors that links CMFD cells with FSRs
   * @details This method is called by the geometry once the CMFD mesh has been
   *          initialized by the geometry. This method allocates a vector for
   *          each CMFD cell that is used to store the FSR ids contained within
   *          that cell.
   * 建立 CMFD粗网格 -> FSR的索引关系
   * _cell_fsrs[cmfd_cell_id] 返回该粗网格内包含的所有FSR ID列表
   * 执行前：_cell_fsrs = []
   * 执行后：_cell_fsrs = [[], [], [], ..., []]  // 共 _local_num_x × _local_num_y × _local_num_z 个空向量
   *_local_num_x：当前进程负责的X方向CMFD单元数
   _local_* 前缀表示并行计算中的区域分解：整个计算域被分割给多个进程，每个进程只处理自己的局部区域
   */
  void Cmfd::initializeCellMap()
  {

    /* Allocate memory for mesh cell FSR vectors 只是给一个域中的所有CMFD单元对应的FSR数组分配内存（每一个CMFD所包含的FSR数组）*/
    for (int z = 0; z < _local_num_z; z++)
    {
      for (int y = 0; y < _local_num_y; y++)
      {
        // 一维index = z * (_local_num_y * _local_num_x) + y * _local_num_x + x
        for (int x = 0; x < _local_num_x; x++)
          _cell_fsrs.push_back(std::vector<long>()); // 创建一个空的long类型动态数组（初始大小为0）
      }
    }
  }

  /**
   * @brief Allocates memory for the CMFD tallies.
   * @details This method is called by the CMFD initialization routine, and
   *          allocates memory for the diffusion, reaction and volume tallies for
   *          every CMFD cells.
   */
  void Cmfd::allocateTallies()
  {

    if (_num_x * _num_y * _num_z == 0)
      log::ferror("Zero cells in CMFD mesh. Please set CMFD mesh before "
                  "initializing CMFD tallies.");

    if (_num_cmfd_groups == 0)
      log::ferror("Zero CMFD gropus. Please set CMFD group structure "
                  "before initializing CMFD tallies.");

    /* Determine tally sizes */
    // unused
    // int num_cells = _num_x * _num_y * _num_z;
    int local_num_cells = _local_num_x * _local_num_y * _local_num_z;
    int tally_size = local_num_cells * _num_cmfd_groups;
    _total_tally_size = 3 * tally_size; // 扩散记录，反应记录和体积记录
    _tally_memory = new CMFD_PRECISION[_total_tally_size];
    CMFD_PRECISION **all_tallies[3];
    for (int t = 0; t < 3; t++)
    { // t = 0 代表是扩散统计（计数器），1代表反应统计（计数器），2代表体积统计（计数器）这里的tally翻译为计数器或统计数值
      all_tallies[t] = new CMFD_PRECISION *[local_num_cells];
      for (int i = 0; i < local_num_cells; i++)
      {                                                  // i为CMFD网格数量的下标
        int idx = i * _num_cmfd_groups + t * tally_size; // idx为_total_tally_size总数量中的全局下标
        all_tallies[t][i] = &_tally_memory[idx];
      }
    }

    /* Assign tallies to allocated data */
    _diffusion_tally = all_tallies[0];
    _reaction_tally = all_tallies[1];
    _volume_tally = all_tallies[2];
    _tallies_allocated = true;
  }

  /**
   * @brief Allocates memory for the Hex CMFD tallies.
   * @details This method is called by the CMFD initialization routine, and
   *          allocates memory for the diffusion, reaction and volume tallies for
   *          every CMFD cells.
   */
  void Cmfd::allocateHexTallies()
  {

    if (_num_x * _num_y * _num_z == 0)
      log::ferror("Zero cells in CMFD mesh. Please set CMFD mesh before "
                  "initializing CMFD tallies.");

    if (_num_cmfd_groups == 0)
      log::ferror("Zero CMFD gropus. Please set CMFD group structure "
                  "before initializing CMFD tallies.");

    /* Determine tally sizes */
    // unused
    // int num_cells = _num_x * _num_y * _num_z;
    int local_num_cells = _local_num_x * _local_num_y * _local_num_z - _empty_cells_num;
    int tally_size = local_num_cells * _num_cmfd_groups;
    _total_tally_size = 3 * tally_size;
    _tally_memory = new CMFD_PRECISION[_total_tally_size];
    CMFD_PRECISION **all_tallies[3];
    for (int t = 0; t < 3; t++)
    {
      all_tallies[t] = new CMFD_PRECISION *[local_num_cells];
      for (int i = 0; i < local_num_cells; i++)
      {
        int idx = i * _num_cmfd_groups + t * tally_size;
        all_tallies[t][i] = &_tally_memory[idx];
      }
    }

    /* Assign tallies to allocated data */
    _diffusion_tally = all_tallies[0];
    _reaction_tally = all_tallies[1];
    _volume_tally = all_tallies[2];
    _tallies_allocated = true;
  }

  /**
   * @brief Initialize and set array that links the MOC energy groups to the
   *        CMFD energy groups.
   * @details This method initializes the _group_indices_map, which is a 1D array
   *           of length _num_moc_groups that maps the MOC energy groups to CMFD
   *           energy groups. The indices into _group_indices_map are the MOC
   *           energy groups and the values are the CMFD energy groups.
   */
  void Cmfd::initializeGroupMap()
  {

    /* Setup one-to-one fine-to-coarse group map if not specified by user */
    if (!_user_group_indices)
    {
      _num_cmfd_groups = _num_moc_groups;

      /* Delete old group indices array if it exists */
      if (_group_indices != NULL)
        delete[] _group_indices;

      /* Allocate memory for new group indices */
      _group_indices = new int[_num_cmfd_groups + 1];

      /* Populate a 1-to-1 mapping from MOC to CMFD groups */
      for (int i = 0; i <= _num_cmfd_groups; i++)
      {
        _group_indices[i] = i;
      }
    }
    else
    {
      if (_num_moc_groups != _group_indices[_num_cmfd_groups])
        log::ferror("The CMFD coarse group mapping is specified for "
                    "%d groups, but the MOC problem contains %d groups",
                    _group_indices[_num_cmfd_groups], _num_moc_groups);
    }

    /* Delete old group indices map if it exists */
    if (_group_indices_map != NULL)
      delete[] _group_indices_map;

    /* Allocate memory for new group indices map */
    _group_indices_map = new int[_num_moc_groups];

    /* Create group indices map */
    for (int e = 0; e < _num_cmfd_groups; e++)
    {
      for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
      {
        _group_indices_map[h] = e;
      }
    }
  }

  /**
   * @brief Find the CMFD surface that a LocalCoords object lies on.
   * @details If the coords is not on a surface, -1 is returned. Otherwise,
   *          the surface ID is returned.
   * @param The CMFD cell ID that the local coords is in.
   * @param The coords being evaluated.
   * @return The surface ID.
   *
   * 查找 LocalCoords 对象所在的 CMFD 表面索引
   */
  int Cmfd::findCmfdSurface(int cell, LocalCoords *coords, double azim, double polar)
  {
    double surface_polar = polar; // 把输入 polar 角复制到局部变量，防止后续修改原参数
    Point *point = coords->getHighestLevel()->getPoint();
    cell = getGlobalCMFDCell(cell); // 获取一维全局ID
    return _lattice->getLatticeSurface(cell, point, azim, surface_polar);
  }

  /*
   * @brief Quickly finds a 3D CMFD surface given a cell, global coordinate, and
   *        2D CMFD surface. Intended for use in axial on-the-fly ray tracing.
   * @details If the coords is not on a surface, -1 is returned. If there is
   *          no 2D CMFD surface intersection, -1 should be input for the 2D CMFD
   *          surface.
   * @param cell_id The CMFD cell ID that the local coords is in.
   * @param z the axial height in the root universe of the point being evaluated.
   * @param surface_2D The ID of the 2D CMFD surface that the LocalCoords object
   *        intersects. If there is no 2D intersection, -1 should be input.
   */
  int Cmfd::findCmfdSurfaceOTF(int cell_id, double z, int surface_2D)
  {
    if (_hexlattice_enable) // 0000000000
      return _lattice->getLatticeSurfaceOTF(cell_id, z, surface_2D);
    else
    {
      int global_cell_id = getGlobalCMFDCell(cell_id);
      return _lattice->getLatticeSurfaceOTF(global_cell_id, z, surface_2D);
    }
  }

  /**
   * @brief Find the CMFD cell that a LocalCoords object is in.
   * @param The coords being evaluated.
   * @return The CMFD cell ID.
   */
  int Cmfd::findCmfdCell(LocalCoords *coords)
  {
    Point *point = coords->getHighestLevel()->getPoint();
    int global_cmfd_cell = _lattice->getLatticeCell(point); // 全局一维索引
    if (_hexlattice_enable)
      return global_cmfd_cell;
    else
    {
      int local_cmfd_cell = getLocalCMFDCell(global_cmfd_cell);
      return local_cmfd_cell;
    }
  }

  /**
   * @brief The structure of the Lattice to be used as the CMFD mesh.
   * @param num_x The number of cells in the x direction.
   * @param num_y The number of cells in the y direction.
   * @param num_z The number of cells in the z direction.
   */
  void Cmfd::setLatticeStructure(int num_x, int num_y, int num_z)
  {
    setNumX(num_x);
    setNumY(num_y);
    setNumZ(num_z);
  }

  /**
   * @brief Returns the Lattice object used as the CMFD mesh.
   * @return A pointer to a Lattice object.
   */
  Lattice *Cmfd::getLattice()
  {
    return _lattice;
  }

  /**
   * @brief Add an FSR ID to a vector that contains all the FSR IDs
   *        contained within a CMFD mesh cell.
   * @param The CMFD cell ID.
   * @param The FSR ID.
   */
  void Cmfd::addFSRToCell(int cmfd_cell, long fsr_id)
  {
    if (_hexlattice_enable)
    { // 0000000000000
      if (CellinHexLattice(cmfd_cell))
      {
        _cell_fsrs.at(cmfd_cell).push_back(fsr_id);
      }
    }
    else
    {
      /** Vector of vectors of FSRs containing in each cell 这个是二维vector数组，第一维确定是具体的哪个CMFD单元，第二个是单个CMFD中所关联的所有FSR
       * std::vector<std::vector<long>> _cell_fsrs;
       * */
      _cell_fsrs.at(cmfd_cell).push_back(fsr_id); // 把fsr_id数据添加到具体的cmfd_cell单元相关联
    }
  }

  /**
   * Print the _cell_fsrs
   */
  void Cmfd::printCellFSRs()
  {
    for (int i = 0; i < _cell_fsrs.size(); i++)
    {
      std::stringstream string; // 用于构建每个单元格信息的stringstream
      string << "CMFD Cell " << i << ": ";
      // 避免在内部循环中重复检查_size，提升性能
      if (_cell_fsrs[i].size() == 0)
      {
        string << "Not within the hexagonal grid";
        string << std::endl;
        log::finfo(string.str().c_str());
        continue;
      }
      else
      {
        for (size_t j = 0; j < _cell_fsrs[i].size(); ++j)
        {
          string << _cell_fsrs[i][j] << ",";
        }
      }

      // 移除最后一个逗号并添加换行符
      if (!string.str().empty())
      {
        string.seekp(-1, std::ios_base::end);
        string << std::endl;
      }

      log::finfo(string.str().c_str());
    }
  }

  /**
   * @brief Set the number of MOC energy groups.
   * @param number of MOC energy groups
   */
  void Cmfd::setNumMOCGroups(int num_groups)
  {
    _num_moc_groups = num_groups;
  }

  /**
   * @brief Get the number of MOC energy groups.
   * @return the number of MOC energy groups
   */
  int Cmfd::getNumMOCGroups()
  {
    return _num_moc_groups;
  }

  /**
   * @brief Get the number of CMFD cells.
   * @return the number of CMFD cells
   */
  int Cmfd::getNumCells()
  {
    return _num_x * _num_y * _num_z;
  }

  /**
   * @brief Get the CMFD group given an MOC group.
   * @param group the MOC energy group
   * @return the CMFD energy group
   */
  int Cmfd::getCmfdGroup(int group)
  {
    return _group_indices_map[group];
  }

  /**
   * @brief set the number of FSRs.
   * @param the number of FSRs
   */
  void Cmfd::setNumFSRs(long num_fsrs)
  {
    _num_FSRs = num_fsrs;
  }

  /**
   * @brief Split the currents of the Mesh cell vertices to the adjacent faces and
   *        edges.
   * 顶点流分割的原理就是,首先找到每个粗网格所属的顶点surface,判断该顶点surface的idx(全局表面索引*能群数)
   * 在_edge_corner_currents存的有没有值,如果没有中子流值,自然也就不用划分了,
   * 如果有则将该值通过getVertexSplitSurfaces函数找到与当前粗网格该顶点相关联的三个面,以及该顶点下一个方向的next粗网格编号的三条边,
   * 并将这些surfaceID存入vector<int> surfaces中,最后取出该顶点在_edge_corner_currents数组中存的值,并除以三,在当前网格的该顶点相关联的
   * 三个面_surface_currents各加1/3值,再通过三条边传递到下一个网格的_edge_corner_currents存入1/3中子流值
   * @details This method takes the currents tallied across the vertices of a CMFD
   *          cell and splits them evenly across the adjacent faces and edges. In
   *          order to transport the current through to the diagonal cell, the
   *          current is also tallied on the edges of the adjacent cells.
   *          Essentially, the tracks that cross through vertices are split into
   *          three one-third-weight tracks as shown in the illustration below.
   *          Face crossings are denoted as "o" and edge crossings are denoted as
   *          "e". As shown, each partial-weight track first crosses a face and
   *          then crosses through an edge to get into an adjacent cell. After all
   *          partial-weight tracks reach the diagonal cell, they are recombined
   *          into one full-weight track. Note tracks are not physically split
   *          into partial-weight tracks for ray tracing; rather tracks cross
   *          through the vertices and the current through each vertex is tallied.
   *
   *                                                       . .      .
   *                                                    .   \|  .
   *                                        |    /   .   .   .
   *                                        |   / .   .       \
   *                                        |  e   .           .
   *                                        | / .           .
   *                                     .  e/           .
   *                 x -------------------.-+---------e-----------
   *                               o   o   /|       .
   *                            .   .     / |    .
   *                         .   .       /  | .
   *                          \ |       /  o|
   *                           .       /.   |
   *                        .   \    ./     |
   *                     .        .   y     z
   *                  .
   *
   */
  void Cmfd::splitVertexCurrents()
  {

    log::fverbose("Splitting CMFD vertex currents...");

    int ncg = _num_cmfd_groups;
    // unused
    // int nf = NUM_FACES;
    // int ne = NUM_EDGES;
    int ns = NUM_SURFACES;
    if (_hexlattice_enable)
      ns = HEX_NUM_SURFACES;

#pragma omp parallel
    {

      FP_PRECISION current;
      std::vector<int> surfaces;
      std::vector<int>::iterator iter;
      std::map<int, CMFD_PRECISION>::iterator it;
      int cell, surface;

#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        int global_id = getGlobalCMFDCell(i); // 域在全局的一维索引

        for (int v = NUM_FACES + NUM_EDGES; v < NUM_SURFACES; v++)
        { // 划分的是所有的顶点

          /* Check if this surface is contained in the map */
          int ind = i * NUM_SURFACES * ncg + v * ncg;
          it = _edge_corner_currents.find(ind);
          if (it == _edge_corner_currents.end()) // 如果该顶点不在map中,则不划分
            continue;

          getVertexSplitSurfaces(global_id, v, &surfaces); // 将该顶点划分,其中将对应的surfaceID存入surfaces数组中,方法包含如何划分以及流的归向

          for (int g = 0; g < ncg; g++)
          {
            /* Divide vertex current by 3 since we will split to 3 surfaces,
             * which propagate through 3 edges 将顶点电流除以3，因为我们将分裂为3个surface，这些surface通过3条边传播到下一个cell粗网格*/
            current = _edge_corner_currents[ind + g] / 3;

            /* Increment current for faces and edges adjacent to this vertex 增加与此顶点相邻的面和边的电流*/
            for (iter = surfaces.begin(); iter != surfaces.end(); ++iter)
            {                         // 循环6次
              cell = (*iter) / ns;    // 取出cell值
              surface = (*iter) % ns; // 取出surface值

              /* Look for the CMFD cell on-domain */
              int local_cell = getLocalCMFDCell(cell); // 域一维局部索引

              /* Add face contributions 添加1/3面的贡献值*/
              if (local_cell != -1)
              {
                if (surface < NUM_FACES)
                { // 其中有三个值是面的
                  _surface_currents->incrementValue(local_cell,
                                                    surface * ncg + g, current);
                }
                else
                { // 添加边的电流(其余三个是边的,也是下一个CMFD网格的电流)

                  /* Check for new index in map */
                  int new_ind = (local_cell * NUM_SURFACES + surface) * ncg + g;
                  std::map<int, CMFD_PRECISION>::iterator it =
                      _edge_corner_currents.find(new_ind);

                  /* If it doesn't exist, initialize to zero */
                  if (it == _edge_corner_currents.end())
                    _edge_corner_currents[new_ind] = 0.0;

                  /* Add the contribution */
                  _edge_corner_currents[new_ind] += current;
                }
              }

              /* Look for the CMFD cell off-domain 如果该网格不在域中*/
              else
              {

                /* Look for the boundary containing the cell */
                for (int s = 0; s < NUM_FACES; s++)
                {

                  std::map<int, int>::iterator it =
                      _boundary_index_map.at(s).find(cell);

                  if (it != _boundary_index_map.at(s).end())
                  {

                    int idx = it->second;

                    /* Add the current to the off-domain split currents cell */
                    _off_domain_split_currents[s][idx][surface * ncg + g] +=
                        current;
                    break;
                  }
                }
              }
            }
            _edge_corner_currents[ind + g] = 0.0; // 将已经处理过的顶点电流值清零
          }
        }
      }
    }
  }

  void Cmfd::splitVertexCurrentsHex()
  {

    log::fverbose("Splitting CMFD vertex currents...");

    int ncg = _num_cmfd_groups;
    // unused
    // int nf = NUM_FACES;
    // int ne = NUM_EDGES;
    int ns = HEX_NUM_SURFACES;

#pragma omp parallel
    {

      FP_PRECISION current;
      std::vector<int> surfaces;
      std::vector<int>::iterator iter;
      std::map<int, CMFD_PRECISION>::iterator it;
      int cell, surface;

#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        int global_id = getGlobalCMFDCell(i);
        if (!_empty_fsrs_cells[global_id])
        {
          int logical_id = _logical_actual_map[global_id];
          for (int v = HEX_NUM_FACES + HEX_NUM_EDGES; v < HEX_NUM_SURFACES; v++)
          {

            /* Check if this surface is contained in the map */
            int ind = logical_id * HEX_NUM_SURFACES * ncg + v * ncg;
            it = _edge_corner_currents.find(ind);
            if (it == _edge_corner_currents.end())
              continue;

            getVertexSplitSurfacesHex(global_id, v, &surfaces);

            for (int g = 0; g < ncg; g++)
            {
              /* Divide vertex current by 3 since we will split to 3 surfaces,
               * which propagate through 3 edges */
              current = _edge_corner_currents[ind + g] / 3;

              /* Increment current for faces and edges adjacent to this vertex */
              for (iter = surfaces.begin(); iter != surfaces.end(); ++iter)
              {
                cell = (*iter) / ns;
                surface = (*iter) % ns;

                /* Look for the CMFD cell on-domain */
                int local_cell = getLocalCMFDCell(cell);

                /* Add face contributions */
                if (local_cell != -1)
                {

                  if (!_empty_fsrs_cells[local_cell])
                  {

                    local_cell = _logical_actual_map[local_cell];

                    if (surface < HEX_NUM_FACES)
                    {
                      _surface_currents->incrementValue(local_cell,
                                                        surface * ncg + g, current);
                    }
                    else
                    {

                      /* Check for new index in map */
                      int new_ind = (local_cell * HEX_NUM_SURFACES + surface) * ncg + g;
                      std::map<int, CMFD_PRECISION>::iterator it =
                          _edge_corner_currents.find(new_ind);

                      /* If it doesn't exist, initialize to zero */
                      if (it == _edge_corner_currents.end())
                        _edge_corner_currents[new_ind] = 0.0;

                      /* Add the contribution */
                      _edge_corner_currents[new_ind] += current;
                    }
                  }
                }

                /* Look for the CMFD cell off-domain */
                else
                {

                  /* Look for the boundary containing the cell */
                  for (int s = 0; s < HEX_NUM_FACES; s++)
                  {

                    std::map<int, int>::iterator it =
                        _boundary_index_map.at(s).find(cell);

                    if (it != _boundary_index_map.at(s).end())
                    {

                      int idx = it->second;

                      /* Add the current to the off-domain split currents cell */
                      _off_domain_split_currents[s][idx][surface * ncg + g] +=
                          current;
                      break;
                    }
                  }
                }
              }
              _edge_corner_currents[ind + g] = 0.0;
            }
          }
        }
      }
    }
  }

  /**
   * @brief Split the currents of the Mesh cell edges to the adjacent faces.
   * @details This method takes the currents tallied across the edges (or corners)
   *          of a CMFD cell and splits them evenly across the adjacent surfaces
   *          (locations 1 and 2). In order to transport the current through to
   *          the diagonal cell, the current is also tallied on the surfaces of
   *          the adjacent cells (locations 3 and 4). Essentially, the tracks that
   *          cross through edges are split into two half-weight tracks as shown
   *          in the illustration below:
   *
   *                                       |    /
   *                                       | __/_
   *                                       |/   /
   *                                     3 /   /
   *                                      /|  / 4
   *                   ------------------/-+-/------------------
   *                                  1 /  |/
   *                                   /   / 2
   *                                  /___/|
   *                                   /   |
   *                                  /    |
   *
   */
  void Cmfd::splitEdgeCurrents()
  {

    log::fverbose("Splitting CMFD edge currents...");

    int ncg = _num_cmfd_groups;
    // unused
    // int nf = NUM_FACES;
    // int ne = NUM_EDGES;
    int ns = NUM_SURFACES;

#pragma omp parallel
    {

      FP_PRECISION current;
      std::vector<int> surfaces;
      std::vector<int>::iterator iter;
      std::map<int, CMFD_PRECISION>::iterator it;
      int cell, surface;

#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        int global_id = getGlobalCMFDCell(i); // 域全局一维索引

        for (int e = NUM_FACES; e < NUM_FACES + NUM_EDGES; e++)
        { // 划分边

          /* Check if this surface is contained in the map */
          int ind = i * NUM_SURFACES * ncg + e * ncg;
          it = _edge_corner_currents.find(ind);
          if (it == _edge_corner_currents.end())
            continue;

          getEdgeSplitSurfaces(global_id, e, &surfaces); // 拆分边为两个面,将surfaceID放入surfaces数组中

          for (int g = 0; g < ncg; g++)
          {
            /* Divide edge current by 2 since we will split to 2 surfaces,
             * which propagate through 2 surfaces 将边电流除以2，因为我们将分裂为2个表面，这些表面通过2个表面传播*/
            current = _edge_corner_currents[ind + g] / 2;

            /* Increment current for faces and edges adjacent to this vertex */
            for (iter = surfaces.begin(); iter != surfaces.end(); ++iter)
            { // 循环4次
              cell = (*iter) / ns;
              surface = (*iter) % ns;

              /* Look for the CMFD cell on-domain */
              int local_cell = getLocalCMFDCell(cell);
              if (local_cell != -1)
              { // 2次加当前网格,2次加方向对应的下一个网格,他们的local_cell不一样
                _surface_currents->incrementValue(local_cell, surface * ncg + g,
                                                  current);
              }

              /* Look for the CMFD cell off-domain */
              else
              {

                /* Look for the boundary containing the cell */
                for (int s = 0; s < NUM_FACES; s++)
                {

                  std::map<int, int>::iterator it =
                      _boundary_index_map.at(s).find(cell);

                  if (it != _boundary_index_map.at(s).end())
                  {

                    int idx = it->second;

                    /* Add the current to the off-domain split currents cell */
                    _off_domain_split_currents[s][idx][surface * ncg + g] +=
                        current;
                    break;
                  }
                }
              }
            }
            _edge_corner_currents[ind + g] = 0.0;
          }
        }
      }
    }
  }

  void Cmfd::splitEdgeCurrentsHex()
  {

    log::fverbose("Splitting CMFD edge currents...");

    int ncg = _num_cmfd_groups;
    // unused
    // int nf = NUM_FACES;
    // int ne = NUM_EDGES;
    int ns = HEX_NUM_SURFACES;

#pragma omp parallel
    {

      FP_PRECISION current;
      std::vector<int> surfaces;
      std::vector<int>::iterator iter;
      std::map<int, CMFD_PRECISION>::iterator it;
      int cell, surface;

#pragma omp for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        int global_id = getGlobalCMFDCell(i);
        if (!_empty_fsrs_cells[global_id])
        {

          for (int e = HEX_NUM_FACES; e < HEX_NUM_FACES + HEX_NUM_EDGES; e++)
          {
            int logical_id = _logical_actual_map[global_id];
            /* Check if this surface is contained in the map */
            int ind = logical_id * HEX_NUM_SURFACES * ncg + e * ncg;
            it = _edge_corner_currents.find(ind);
            if (it == _edge_corner_currents.end())
              continue;

            getEdgeSplitSurfacesHex(global_id, e, &surfaces);

            for (int g = 0; g < ncg; g++)
            {
              /* Divide edge current by 2 since we will split to 2 surfaces,
               * which propagate through 2 surfaces */
              current = _edge_corner_currents[ind + g] / 2;

              /* Increment current for faces and edges adjacent to this vertex */
              for (iter = surfaces.begin(); iter != surfaces.end(); ++iter)
              {
                cell = (*iter) / ns;
                surface = (*iter) % ns;

                /* Look for the CMFD cell on-domain */
                int local_cell = getLocalCMFDCell(cell);
                if (local_cell != -1)
                {
                  if (!_empty_fsrs_cells[local_cell])
                  {
                    local_cell = _logical_actual_map[local_cell];
                    _surface_currents->incrementValue(local_cell, surface * ncg + g,
                                                      current);
                  }
                }

                /* Look for the CMFD cell off-domain */
                else
                {

                  /* Look for the boundary containing the cell */
                  for (int s = 0; s < HEX_NUM_FACES; s++)
                  {

                    std::map<int, int>::iterator it =
                        _boundary_index_map.at(s).find(cell);

                    if (it != _boundary_index_map.at(s).end())
                    {

                      int idx = it->second;

                      /* Add the current to the off-domain split currents cell */
                      _off_domain_split_currents[s][idx][surface * ncg + g] +=
                          current;
                      break;
                    }
                  }
                }
              }
              _edge_corner_currents[ind + g] = 0.0;
            }
          }
        }
      }
    }
  }

  /**
   * @brief Get the faces and edges to split the currents of the Mesh cell
   *        vertices.
   * 获取某个顶点处经过的面和边的ID，用于在网格单元的顶点分割电流时使用
   * @details The process by which the current of tracks passing through vertices
   *          is split is described in the comment for
   *          Cmfd::splitVertexCurrents(). This method takes in the cell and
   *          vertex that is being split as well as a std::vector used to store
   *          the IDs of surfaces that are crossed by the partial-weight tracks.
   *          This method properly accounts for crossings on the geometry
   *          boundaries by applying the corresponding boundary conditions to
   *          split the currents.
   * @param cell The CMFD cell ID that the vertex is in.
   * @param vertex The vertex that the track crosses through.
   * @param surfaces A std::vector that is populated with the IDs of surfaces that
   *        are crossed by partial-weight tracks.
   */
  void Cmfd::getVertexSplitSurfaces(int cell, int vertex,
                                    std::vector<int> *surfaces)
  {

    surfaces->clear();
    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y); // 将一维索引转为三维索引
    int ns = NUM_SURFACES;

    int cell_indexes[3] = {x, y, z};
    int cell_limits[3] = {_num_x, _num_y, _num_z};

    int direction[3];
    convertSurfaceToDirection(vertex, direction);

    /* Get the partial surfaces composing the edge split */
    int remainder_surfaces[3];
    int partial_surfaces[3];
    // 下面两层循环的作用是把顶点对应的三个方向分别传给partial_direction数组,比如direction={1,1,1}就转为了partial_direction[0]={1,0,0}
    // partial_direction[1]={0,1,0},partial_direction[2]={0,0,1},其中remainder_direction是用来传递给下一个CMFD的
    // 也就是说partial_direction用于当前网格(也就是1/3),其余的(2/3)传递给下一个网格
    for (int i = 0; i < 3; i++)
    { // 遍历三个维度x,y,z
      int remainder_direction[3];
      int partial_direction[3];

      for (int j = 0; j < 3; j++)
      { // 遍历三个维度x,y,z,给remainder_direction和partial_direction赋值
        if (i == j)
        {                                      // 对应上就把direction相应的值传给partial_direction,其余为0
          remainder_direction[j] = 0;          //{0}
          partial_direction[j] = direction[j]; //{1}
        }
        else
        {                                        // 对应不上就把direction相应的值传给remainder_direction,其余为0
          remainder_direction[j] = direction[j]; //{0,1}->{0,1,1}
          partial_direction[j] = 0;              //{1,0}->{1,0,0}
        }
      }
      remainder_surfaces[i] = convertDirectionToSurface(remainder_direction); // 将方向值转为surface编号值,其中remainder_direction中有两个非0,转为了对应边的序号
      partial_surfaces[i] = convertDirectionToSurface(partial_direction);     // partial_direction有一个非零,转为了对应面的surface序号
    }

    /* Treat all partial surfaces */
    for (int i = 0; i < 3; i++)
    {

      int remainder_surface = remainder_surfaces[i]; // 取出该surface值
      int partial_surface = partial_surfaces[i];

      surfaces->push_back(cell * ns + partial_surface); // 转为全局surface编号

      /* Tally current on neighboring cell or appropriate boundary */
      int cell_next = getCellNext(cell, partial_surface);
      if ((cell_indexes[i] == 0 && direction[i] == -1) ||
          (cell_indexes[i] == cell_limits[i] - 1 && direction[i] == +1))
      {
        if (_boundaries[partial_surface] == REFLECTIVE)
        {                                                     // 处于边界情况
          surfaces->push_back(cell * ns + remainder_surface); // 当前部分表面是反射边界，则将电流反射回当前单元
        }
        else if (_boundaries[partial_surface] == PERIODIC)
        { // 当前部分表面是周期性边界，则将电流传播到下一个单元的相应位置
          surfaces->push_back(cell_next * ns + remainder_surface);
        }
      }
      else
      {
        surfaces->push_back(cell_next * ns + remainder_surface); // 其余情况将电流传播到下一个单元的相应位置
      }
    }
  }

  void Cmfd::getVertexSplitSurfacesHex(int cell, int vertex,
                                       std::vector<int> *surfaces)
  {

    surfaces->clear();
    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y);
    int ns = HEX_NUM_SURFACES;

    int cell_indexes[3] = {x, y, z};
    int cell_limits[3] = {_num_x, _num_y, _num_z};

    int direction[3];

    convertSurfaceToDirectionHex(vertex, direction);

    /* Get the partial surfaces composing the edge split */
    int remainder_surfaces[3];
    int partial_surfaces[3];

    convertVertexToSurfaces(vertex, remainder_surfaces, partial_surfaces);

    /* Treat all partial surfaces */
    for (int i = 0; i < 3; i++)
    {

      int remainder_surface = remainder_surfaces[i];
      int partial_surface = partial_surfaces[i];

      surfaces->push_back(cell * ns + partial_surface);

      /* Tally current on neighboring cell or appropriate boundary */
      int cell_next = getCellNext(cell, partial_surface);

      if (CellinXYZBoundary(x, y, z, partial_surface))
      {
        if (_boundaries[partial_surface] == REFLECTIVE)
        {
          surfaces->push_back(cell * ns + remainder_surface);
        }
        else if (_boundaries[partial_surface] == PERIODIC)
        {
          if (cell_next != -1 && !_empty_fsrs_cells[cell_next])
            surfaces->push_back(cell_next * ns + remainder_surface);
        }
      }
      else
      {
        if (cell_next == -1 || _empty_fsrs_cells[cell_next])
          continue;
        else
        {
          surfaces->push_back(cell_next * ns + remainder_surface);
        }
      }
    }
  }

  /**
   * @brief Hex Lattice uses an XYZ coordinate system to store the layout,
   *        Detect whether the Cell is located at the radial/axial boundary.
   * @details
   *  y
   *  |———————————————————————————
   *  | *  *  *  *  *  -  -  -  - |
   *  | *  *  *  *  *  *  -  -  - |
   *  | *  *  *  *  *  *  *  -  - |
   *  | *  *  *  *  *  *  *  *  - |
   *  | *  *  *  *  *  *  *  *  * |
   *  | -  *  *  *  *  *  *  *  * |
   *  | -  -  *  *  *  *  *  *  * |
   *  | -  -  -  *  *  *  *  *  * |
   *  | -  -  -  -  *  *  *  *  * |
   *   ——————————————————————————————x
   *
   */
  bool Cmfd::CellinXYZBoundary(int x, int y, int z, int surface)
  {
    bool inboundary = false;
    int r = _num_r;
    int length = 2 * (r - 1);
    int mid = r - 1;
    int high = (3 * r) - 3;

    if (x == 0 && y >= mid)
    { // 左边界,包含上下两端点
      if (y == length && (surface == 0 || surface == 1 || surface == 6))
        inboundary = true;
      else if (y == mid && (surface == 0 || surface == 1 || surface == 2))
        inboundary = true;
      else if (y > mid && y < length && (surface == 0 || surface == 1))
        inboundary = true;
    }
    else if (x == length && y <= mid)
    { // 右边界,包含上下两端点
      if (y == 0 && (surface == 2 || surface == 4 || surface == 5))
        inboundary = true;
      else if (y == mid && (surface == 4 || surface == 5 || surface == 6))
        inboundary = true;
      else if (y > 0 && y < mid && (surface == 4 || surface == 5))
        inboundary = true;
    }
    else if (y == 0 && x >= mid && x < length)
    { // 下边界，只包含左边界端点
      if (x == mid && (surface == 0 || surface == 2 || surface == 5))
        inboundary = true;
      else if (x > mid && (surface == 2 || surface == 5))
        inboundary = true;
    }
    else if (y == length && x > 0 && x <= mid)
    { // 上边界，只包含右边界端点
      if (x == mid && (surface == 1 || surface == 4 || surface == 6))
        inboundary = true;
      else if (x < mid && (surface == 1 || surface == 6))
        inboundary = true;
    }
    else if (x < mid && y < mid && (x + y == mid))
    { // 左下边界，不含端点
      if (surface == 0 || surface == 3)
        inboundary = true;
    }
    else if (x > mid && y > mid && (x + y == high))
    { // 右上边界，不含端点
      if (surface == 4 || surface == 6)
        inboundary = true;
    }

    if (z == 0 && surface == 3)
      inboundary = true;
    else if (z == (_num_z - 1) && surface == 7)
      inboundary = true;

    return inboundary;
  }

  /**
   * @brief Hex Lattice uses an XYZ coordinate system to store the layout,
   *        Detect whether the Cell is located at the radial/axial boundary.
   */
  bool Cmfd::CellinXYZBoundary(int cell, int surface)
  {
    bool inboundary = false;

    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y);

    int r = _num_r;
    int length = 2 * (r - 1);
    int mid = r - 1;
    int high = (3 * r) - 3;

    if (x == 0 && y >= mid)
    { // 左边界,包含上下两端点
      if (y == length && (surface == 0 || surface == 1 || surface == 6))
        inboundary = true;
      else if (y == mid && (surface == 0 || surface == 1 || surface == 2))
        inboundary = true;
      else if (y > mid && y < length && (surface == 0 || surface == 1))
        inboundary = true;
    }
    else if (x == length && y <= mid)
    { // 右边界,包含上下两端点
      if (y == 0 && (surface == 2 || surface == 4 || surface == 5))
        inboundary = true;
      else if (y == mid && (surface == 4 || surface == 5 || surface == 6))
        inboundary = true;
      else if (y > 0 && y < mid && (surface == 4 || surface == 5))
        inboundary = true;
    }
    else if (y == 0 && x >= mid && x < length)
    { // 下边界，只包含左边界端点
      if (x == mid && (surface == 0 || surface == 2 || surface == 5))
        inboundary = true;
      else if (x > mid && (surface == 2 || surface == 5))
        inboundary = true;
    }
    else if (y == length && x > 0 && x <= mid)
    { // 上边界，只包含右边界端点
      if (x == mid && (surface == 1 || surface == 4 || surface == 6))
        inboundary = true;
      else if (x < mid && (surface == 1 || surface == 6))
        inboundary = true;
    }
    else if (x < mid && y < mid && (x + y == mid))
    { // 左下边界，不含端点
      if (surface == 0 || surface == 3)
        inboundary = true;
    }
    else if (x > mid && y > mid && (x + y == high))
    { // 右上边界，不含端点
      if (surface == 4 || surface == 6)
        inboundary = true;
    }

    if (z == 0 && surface == 3)
      inboundary = true;
    else if (z == (_num_z - 1) && surface == 7)
      inboundary = true;

    return inboundary;
  }

  /**
   * @brief Hex Lattice uses an XY coordinate system to store the layout,
   *        Detect whether the Cell is located at the radial/axial boundary.
   */
  bool Cmfd::CellinXYBoundary(int cell)
  {
    bool inboundary = false;

    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y);

    int r = _num_r;
    int length = 2 * (r - 1);
    int mid = r - 1;
    int high = (3 * r) - 3;

    if (x == 0 && y >= mid)
    { // 左边界
      inboundary = true;
    }
    else if (x == length && y <= mid)
    { // 右边界,包含上下两端点
      inboundary = true;
    }
    else if (y == 0 && x >= mid && x < length)
    { // 下边界，只包含左边界端点
      inboundary = true;
    }
    else if (y == length && x > 0 && x <= mid)
    { // 上边界，只包含右边界端点
      inboundary = true;
    }
    else if (x < mid && y < mid && (x + y == mid))
    { // 左下边界，不含端点
      inboundary = true;
    }
    else if (x > mid && y > mid && (x + y == high))
    { // 右上边界，不含端点
      inboundary = true;
    }

    return inboundary;
  }

  /**
   * @brief Hex Lattice uses an XY coordinate system to store the layout,
   *        Detect whether the Cell is located at the radial/axial boundary.
   */
  bool Cmfd::CellinXYBoundaryWithStencil(int cell, int stencil_index)
  {
    bool stencilnext = false;

    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y);

    int r = _num_r;
    int length = 2 * (r - 1);
    int mid = r - 1;
    int high = (3 * r) - 3;

    if (x == 0 && y >= mid)
    { // 左边界,包含上下两端点
      if (y == length && stencil_index != 2 && stencil_index != 5 && stencil_index != 6)
        stencilnext = true;
      else if (y == mid && stencil_index != 0 && stencil_index != 2 && stencil_index != 5)
        stencilnext = true;
      else if (y > mid && y < length && stencil_index != 2 && stencil_index != 5)
        stencilnext = true;
    }
    else if (x == length && y <= mid)
    { // 右边界,包含上下两端点
      if (y == 0 && stencil_index != 0 && stencil_index != 1 && stencil_index != 4)
        stencilnext = true;
      else if (y == mid && stencil_index != 1 && stencil_index != 4 && stencil_index != 6)
        stencilnext = true;
      else if (y > 0 && y < mid && stencil_index != 1 && stencil_index != 4)
        stencilnext = true;
    }
    else if (y == 0 && x >= mid && x < length)
    { // 下边界，只包含左边界端点
      if (x == mid && stencil_index != 0 && stencil_index != 1 && stencil_index != 2)
        stencilnext = true;
      else if (x > mid && stencil_index != 0 && stencil_index != 1)
        stencilnext = true;
    }
    else if (y == length && x > 0 && x <= mid)
    { // 上边界，只包含右边界端点
      if (x == mid && stencil_index != 4 && stencil_index != 5 && stencil_index != 6)
        stencilnext = true;
      else if (x < mid && stencil_index != 5 && stencil_index != 6)
        stencilnext = true;
    }
    else if (x < mid && y < mid && (x + y == mid))
    { // 左下边界，不含端点
      if (stencil_index != 0 && stencil_index != 2)
        stencilnext = true;
    }
    else if (x > mid && y > mid && (x + y == high))
    { // 右上边界，不含端点
      if (stencil_index != 4 && stencil_index != 6)
        stencilnext = true;
    }

    return stencilnext;
  }

  /**
   * @brief Get the faces to split the currents of the Mesh cell edges.
   * @details The process by which the current of tracks passing through edges
   *          is split is described in the comment for Cmfd::splitEdgeCurrents().
   *          This method takes in the cell and edge that is being split as well
   *          as a std::vector used to store the IDs of surfaces that are crossed
   *          by the partial-weight tracks. This method properly accounts for
   *          crossings on the geometry boundaries by applying the corresponding
   *          boundary conditions to split the currents.
   * @param cell The CMFD cell ID that the edge is in.
   * @param edge The edge that the track crosses through.
   * @param surfaces A std::vector that is populated with the IDs of surfaces that
   *        are crossed by partial-weight tracks.
   */
  void Cmfd::getEdgeSplitSurfaces(int cell, int edge,
                                  std::vector<int> *surfaces)
  {

    surfaces->clear();
    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y);
    int ns = NUM_SURFACES;
    if (_hexlattice_enable)
      ns = HEX_NUM_SURFACES;

    int cell_indexes[3] = {x, y, z};
    int cell_limits[3] = {_num_x, _num_y, _num_z};

    int direction[3];
    convertSurfaceToDirection(edge, direction); // 对于边来说,direction中有一个非0

    /* Get the partial surfaces composing the edge split */
    int partial_surfaces[2];
    // unused
    // int opposite_surfaces[2];
    int ind = 0;

    for (int i = 0; i < 3; i++)
    { // 比如边为{1,1,0},则被拆分为partial_direction[0]={1,0,0},partial_direction[0]={0,1,0}
      if (direction[i] != 0)
      {
        int partial_direction[3] = {0, 0, 0};
        partial_direction[i] = direction[i];
        partial_surfaces[ind] = convertDirectionToSurface(partial_direction); // 获取对应面的surface
        ind++;
      }
    }

    /* Treat all partial surfaces */
    ind = 0;
    for (int i = 0; i < 3; i++)
    {
      if (direction[i] != 0)
      {

        int partial_surface = partial_surfaces[ind];   // 取出面surface
        int other_surface = partial_surfaces[1 - ind]; // 取出另一个surface

        surfaces->push_back(cell * ns + partial_surface);

        /* Tally current on neighboring cell or appropriate boundary */
        int cell_next = getCellNext(cell, partial_surface); // 找到对应的下一个cell
        if (cell_next == -1)
          continue;
        else
        {
          if ((cell_indexes[i] == 0 && direction[i] == -1) ||
              (cell_indexes[i] == cell_limits[i] - 1 && direction[i] == +1))
          {
            if (_boundaries[partial_surface] == REFLECTIVE)
            {
              surfaces->push_back(cell * ns + other_surface); // 当前部分表面是反射边界，则将电流反射回当前单元
            }
            else if (_boundaries[partial_surface] == PERIODIC)
            {
              surfaces->push_back(cell_next * ns + other_surface); // 当前部分表面是周期性边界，则将电流传播到下一个单元的相应位置
            }
          }
          else
          {
            surfaces->push_back(cell_next * ns + other_surface); // 其余情况将电流传播到下一个单元的相应位置
          }
          ind++;
        }
      }
    }
  }

  void Cmfd::getEdgeSplitSurfacesHex(int cell, int edge,
                                     std::vector<int> *surfaces)
  {

    surfaces->clear();
    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y);
    int ns = HEX_NUM_SURFACES;

    int cell_indexes[3] = {x, y, z};
    int cell_limits[3] = {_num_x, _num_y, _num_z};

    int direction[3];
    convertSurfaceToDirectionHex(edge, direction);

    /* Get the partial surfaces composing the edge split */
    int partial_surfaces[2];
    // unused
    // int opposite_surfaces[2];
    int ind = 0;

    convertEdgeToSurfaces(edge, partial_surfaces);

    /* Treat all partial surfaces */
    ind = 0;

    for (int i = 0; i < 2; i++)
    {
      int partial_surface = partial_surfaces[ind];
      int other_surface = partial_surfaces[1 - ind];

      surfaces->push_back(cell * ns + partial_surface);

      /* Tally current on neighboring cell or appropriate boundary */
      int cell_next = getCellNext(cell, partial_surface);

      if (CellinXYZBoundary(x, y, z, partial_surface))
      {
        if (_boundaries[partial_surface] == REFLECTIVE)
        {
          surfaces->push_back(cell * ns + other_surface);
        }
        else if (_boundaries[partial_surface] == PERIODIC)
        {
          if (cell_next != -1 && !_empty_fsrs_cells[cell_next])
            surfaces->push_back(cell_next * ns + other_surface);
        }
      }
      else
      {
        if (cell_next == -1 || _empty_fsrs_cells[cell_next])
          continue;
        else
          surfaces->push_back(cell_next * ns + other_surface);
      }
      ind++;
    }
  }

  /**
   * @brief Get the ID of the Mesh cell next to given Mesh cell.
   * @param cell index of the current CMFD cell
   * @param surface_id id of the surface between the current cell and the next
   * @param global work at the global (all domains together) level
   * @param neighbor give cell in neighboring domain
   * @return neighboring CMFD cell ID
   */
  int Cmfd::getCellNext(int cell, int surface_id, bool global, bool neighbor)
  {

    int cell_next = -1;

    int x, y, z;
    int nx, ny, nz;
    int x_global, y_global, z_global;
    if (global || _domain_communicator == NULL)
    { // 如果没有域分解
      x_global = (cell % (_num_x * _num_y)) % _num_x;
      y_global = (cell % (_num_x * _num_y)) / _num_x;
      z_global = cell / (_num_x * _num_y);
      x = x_global;
      y = y_global;
      z = z_global;
      nx = _num_x;
      ny = _num_y;
      nz = _num_z;
    }
    else
    {
      x = (cell % (_local_num_x * _local_num_y)) % _local_num_x;
      y = (cell % (_local_num_x * _local_num_y)) / _local_num_x;
      z = cell / (_local_num_x * _local_num_y);
      x_global = x + _domain_communicator->_domain_idx_x * _local_num_x;
      y_global = y + _domain_communicator->_domain_idx_y * _local_num_y;
      z_global = z + _domain_communicator->_domain_idx_z * _local_num_z;
      nx = _local_num_x;
      ny = _local_num_y;
      nz = _local_num_z;
    }

    if (_hexlattice_enable)
    {
      /* Find the cell on the other side of the surface */
      if (surface_id == HEX_SURFACE_BETA_MIN)
      {
        if (x != 0 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell - 1;
        else if (neighbor && !global && x_global != 0)
          cell_next = z * _local_num_y + y;
        else if (_boundaries[HEX_SURFACE_BETA_MIN] == PERIODIC && CellinXYZBoundary(cell, surface_id))
        {
          // cell + (_num_x-1);
          if (y <= ((_num_y - 1) / 2))
            cell_next = (z * (_num_x * _num_y)) + (y * _num_x) + (_num_x - 1);
          else
            cell_next = (z * (_num_x * _num_y)) + (y * _num_x) + (_num_x - (y % ((_num_y - 1) / 2)) - 1);
        }
      }

      else if (surface_id == HEX_SURFACE_GAMMA_MIN)
      {
        if (x != 0 && y != ny - 1 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell - 1 + nx;
        else if (neighbor && !global && y_global != 0)
          cell_next = z * _local_num_x + x;
        else if (_boundaries[HEX_SURFACE_GAMMA_MIN] == PERIODIC && CellinXYZBoundary(cell, surface_id))
        {
          // cell_next = y + x * _num_x;
          cell_next = (z * (_num_x * _num_y)) + (x * _num_x) + y;
        }
      }

      else if (surface_id == HEX_SURFACE_DELTA_MIN)
      {
        if (y != 0 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell - nx;
        else if (neighbor && !global && y_global != 0)
          cell_next = z * _local_num_x + x;
        else if (_boundaries[HEX_SURFACE_DELTA_MIN] == PERIODIC && CellinXYZBoundary(cell, surface_id))
        {
          // cell_next = cell + _num_x*(_num_y-1);
          if (x <= ((_num_x - 1) / 2))
            cell_next = (z * (_num_x * _num_y)) + (_num_x * (_num_y - 1)) + x;
          else
            cell_next = (z * (_num_x * _num_y)) + (_num_x * (_num_y - (x % ((_num_x - 1) / 2)) - 1)) + x;
        }
      }

      else if (surface_id == HEX_SURFACE_Z_MIN)
      {
        if (z != 0 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell - nx * ny;
        else if (neighbor && !global && z_global != 0)
          cell_next = y * _local_num_x + x;
        else if (_boundaries[HEX_SURFACE_Z_MIN] == PERIODIC)
          cell_next = cell + _num_x * _num_y * (_num_z - 1);
      }

      else if (surface_id == HEX_SURFACE_BETA_MAX)
      {
        if (x != nx - 1 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell + 1;
        else if (neighbor && !global && x_global != _num_x - 1)
          cell_next = z * _local_num_y + y;
        else if (_boundaries[HEX_SURFACE_BETA_MAX] == PERIODIC && CellinXYZBoundary(cell, surface_id))
        {
          // cell_next = cell - (_num_x-1);
          if (y >= ((_num_y - 1) / 2))
            cell_next = (z * (_num_x * _num_y)) + (y * _num_x);
          else
            cell_next = (z * (_num_x * _num_y)) + (y * _num_x) + (((_num_x - 1) / 2) - y);
        }
      }

      else if (surface_id == HEX_SURFACE_GAMMA_MAX)
      {
        if (x != nx - 1 && y != 0 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell + 1 - nx;
        else if (neighbor && !global && y_global != _num_y - 1)
          cell_next = z * _local_num_x + x;
        else if (_boundaries[HEX_SURFACE_GAMMA_MAX] == PERIODIC && CellinXYZBoundary(cell, surface_id))
          // cell_next = y + x * _num_x;
          cell_next = (z * (_num_x * _num_y)) + (x * _num_x) + y;
      }

      else if (surface_id == HEX_SURFACE_DELTA_MAX)
      {
        if (y != ny - 1 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell + nx;
        else if (neighbor && !global && y_global != _num_y - 1)
          cell_next = z * _local_num_x + x;
        else if (_boundaries[HEX_SURFACE_DELTA_MAX] == PERIODIC && CellinXYZBoundary(cell, surface_id))
        {
          // cell_next = cell - _num_x*(_num_y-1);
          if (x >= ((_num_x - 1) / 2))
            cell_next = (z * (_num_x * _num_y)) + x;
          else
            cell_next = (z * (_num_x * _num_y)) + (_num_x * ((_num_x - 1) / 2) - x) + x;
        }
      }

      else if (surface_id == HEX_SURFACE_Z_MAX)
      {
        if (z != nz - 1 && !CellinXYZBoundary(cell, surface_id))
          cell_next = cell + nx * ny;
        else if (neighbor && !global && z_global != _num_z - 1)
          cell_next = y * _local_num_x + x;
        else if (_boundaries[HEX_SURFACE_Z_MAX] == PERIODIC)
          cell_next = cell - _num_x * _num_y * (_num_z - 1);
      }
    }
    else
    {
      /* Find the cell on the other side of the surface */
      if (surface_id == SURFACE_X_MIN)
      {
        if (x != 0)
          cell_next = cell - 1; // 一般x负方向CMFD编号直接减1
        else if (neighbor && !global && x_global != 0)
          cell_next = z * _local_num_y + y; // 默认global=true,neighbor=false
        else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
          cell_next = cell + (_num_x - 1); // 如果是周期边界,则最小的x变为最大的cell + (_num_x-1)
      }

      else if (surface_id == SURFACE_Y_MIN)
      {
        if (y != 0)
          cell_next = cell - nx; // 顺序是先x后y后z,所以这里是-nx
        else if (neighbor && !global && y_global != 0)
          cell_next = z * _local_num_x + x;
        else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
          cell_next = cell + _num_x * (_num_y - 1);
      }

      else if (surface_id == SURFACE_Z_MIN)
      {
        if (z != 0)
          cell_next = cell - nx * ny; // 同上
        else if (neighbor && !global && z_global != 0)
          cell_next = y * _local_num_x + x;
        else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
          cell_next = cell + _num_x * _num_y * (_num_z - 1);
      }

      else if (surface_id == SURFACE_X_MAX)
      {
        if (x != nx - 1)
          cell_next = cell + 1;
        else if (neighbor && !global && x_global != _num_x - 1)
          cell_next = z * _local_num_y + y;
        else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
          cell_next = cell - (_num_x - 1);
      }

      else if (surface_id == SURFACE_Y_MAX)
      {
        if (y != ny - 1)
          cell_next = cell + nx;
        else if (neighbor && !global && y_global != _num_y - 1)
          cell_next = z * _local_num_x + x;
        else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
          cell_next = cell - _num_x * (_num_y - 1);
      }

      else if (surface_id == SURFACE_Z_MAX)
      {
        if (z != nz - 1)
          cell_next = cell + nx * ny;
        else if (neighbor && !global && z_global != _num_z - 1)
          cell_next = y * _local_num_x + x;
        else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
          cell_next = cell - _num_x * _num_y * (_num_z - 1);
      }
    }
    // log::fdebug("Cell is %d, surface id is %d, cell next is %d, x:%d,y:%d,z:%d", cell, surface_id, cell_next, x, y, z);
    return cell_next;
  }

  /**
   * @brief Set the CMFD boundary type for a given surface.
   * @details The CMFD boundary is assumed to be rectangular with the
   *          surfaces identified by constants in the constants.h file.
   * @param side The CMFD surface UID.
   * @param boundary The boundaryType of the surface.
   */
  void Cmfd::setBoundary(int side, boundaryType boundary)
  {
    _boundaries[side] = boundary;
  }

  /**
   * @brief Get the boundaryType for one side of the CMFD mesh.
   * @param the CMFD mesh surface ID.
   * @return the boundaryType for the surface.
   */
  int Cmfd::getBoundary(int side)
  {
    return _boundaries[side];
  }

  /**
   * @brief Return the CMFD cell ID that a FSR lies in.
   * @details Note that a CMFD cell is not an actual Cell object; rather, a CMFD
   *          cell is just a way of describing each of the rectangular regions
   *          that make up a CMFD lattice. CMFD cells are numbered with 0 in the
   *          lower left corner and monotonically increasing from left to right,
   *          from bottom to top. For example, the indices for a 4 x 4 lattice
   *          are:
   *                  12  13  14  15
   *                  8    9  10  11
   *                  4    5   6   7
   *                  0    1   2   3
   * @param fsr_id the FSR ID.
   * @return The CMFD cell ID. Return -1 if cell is not found.
   */
  int Cmfd::convertFSRIdToCmfdCell(long fsr_id)
  {

    std::vector<long>::iterator iter;
    for (int cell_id = 0; cell_id < _local_num_x * _local_num_y * _local_num_z;
         cell_id++)
    {
      for (iter = _cell_fsrs.at(cell_id).begin();
           iter != _cell_fsrs.at(cell_id).end(); ++iter)
      {
        if (*iter == fsr_id)
          return cell_id;
      }
    }

    return -1;
  }

  void Cmfd::findEmptyCmfdCells()
  {
    int count = 0;
    // 遍历 _cell_fsrs 中的每个CMFD cell
    for (int i = 0; i < _cell_fsrs.size(); ++i)
    {
      if (_cell_fsrs[i].empty())
      {
        _empty_fsrs_cells.push_back(true);
      }
      else
      {
        _empty_fsrs_cells.push_back(false);
      }
    }

    for (int i = 0; i < _empty_fsrs_cells.size(); ++i)
    {
      if (_empty_fsrs_cells[i])
      {
        log::fdebug("CMFD Cell %d is empty", i);
        count += 1;
      }
    }
    _empty_cells_num = count;

    _logical_actual_map.resize(_cell_fsrs.size(), -1);
    int j = 0;
    for (int i = 0; i < _empty_fsrs_cells.size(); ++i)
    {
      if (!_empty_fsrs_cells[i])
      {
        _logical_actual_map[i] = j;
        log::fdebug("Logical cmfd cell is %d, actual cmfd cell is %d", j, i);
        j += 1;
      }
    }

    log::fdebug("Total CMFD Cell is %d", _cell_fsrs.size());
    log::fdebug("Total empty CMFD Cell is %d", count);
  }

  /**
   * @brief Return the CMFD cell ID that a FSR lies in.
   * @param global_fsr_id The global FSR ID.
   * @return The CMFD cell ID.
   */
  int Cmfd::convertGlobalFSRIdToCmfdCell(long global_fsr_id)
  {

    /* Determine the domain and local FSR ID */
    int cmfd_cell = -1;
    if (!mpi::isSpatialDecomposed())
    {
      cmfd_cell = convertFSRIdToCmfdCell(global_fsr_id);
    }
#ifdef ENABLE_MPI_
    else
    {

      long fsr_id;
      int domain;
      _geometry->getLocalFSRId(global_fsr_id, fsr_id, domain);

      /* Get the FSR centroid in the correct domain */
      int rank;
      MPI_Comm comm = _geometry->getMPICart();
      MPI_Comm_rank(comm, &rank);
      int temp_cmfd_cell = 0;
      if (rank == domain)
        temp_cmfd_cell = convertFSRIdToCmfdCell(fsr_id);

      /* Broadcast the temp_cmfd_cell */
      MPI_Allreduce(&temp_cmfd_cell, &cmfd_cell, 1, MPI_INT, MPI_SUM, comm);
    }
#endif
    return cmfd_cell;
  }

  /**
   * @brief Return a pointer to the vector of vectors that contains
   *        the FSRs that lie in each cell.
   * @return Vector of vectors containing FSR IDs in each cell.
   */
  std::vector<std::vector<long>> *Cmfd::getCellFSRs()
  {
    return &_cell_fsrs;
  }

  /**
   * @brief Set the vector of vectors that contains.
   *        the FSRs that lie in each cell.
   * @param Vector of vectors containing FSR IDs in each cell.
   */
  void Cmfd::setCellFSRs(std::vector<std::vector<long>> *cell_fsrs)
  {

    if (!_cell_fsrs.empty())
    {
      std::vector<std::vector<long>>::iterator iter;
      for (iter = _cell_fsrs.begin(); iter != _cell_fsrs.end(); ++iter)
        iter->clear();
      _cell_fsrs.clear();
    }

    _cell_fsrs = *cell_fsrs;
  }

  /**
   * @brief Set flag indicating whether to update the MOC flux.
   * @param flux_update_on Flag saying whether to update MOC flux.
   */
  void Cmfd::setFluxUpdateOn(bool flux_update_on)
  {
    _flux_update_on = flux_update_on;
  }

  /**
   * @brief Sets the a ConvergenceData object to record diagnostics
   * @details The ConvergenceData object records the number of fission source
   *          and flux iterations for the CMFD solver as well as the maximum
   *          magnitude prolongation factor
   * @param convergence_data The convergence data object
   */
  void Cmfd::setConvergenceData(ConvergenceData *convergence_data)
  {
    _convergence_data = convergence_data;
  }

  /**
   * @brief Set the flag indicating whether to use quadratic axial interpolation
   *        for update ratios.
   * @param interpolate flag meaning No interpolation(0), FSR axially averaged
            value(1) or centroid z-coordinate evaluated value(2)
   */
  void Cmfd::useAxialInterpolation(int interpolate)
  {

    if (interpolate < 0 || interpolate > 2)
      log::ferror("interpolate can only has value 0, 1, or 2, respectively"
                  " meaning No interpolation, FSR axially averaged value or"
                  " centroid z-coordinate evaluated value");
    if (interpolate == 1 || interpolate == 2)
      log::fwarn_once("Axial interpolation CMFD prolongation may only"
                      " be effective when all the FSRs are axially homogeneous");
    _use_axial_interpolation = interpolate;
  }

  /**
   * @brief Turns on the flux limiting condition
   * @details If the CMFD correction diffusion coefficient is larger than the
   *          diffusion coefficient, recompute the diffusion coefficient as the
   *          ratio of current to twice the flux, and re-compute a correction
   *          diffusion coefficient.
   * @param flux_limiting whether to turn on the flux limiting condition
   */
  void Cmfd::useFluxLimiting(bool flux_limiting)
  {
    _flux_limiting = flux_limiting;
  }

  /**
   * @brief Modifies the diagonal element to be consistent with the MOC solve
   * @details This function re-computes a new total cross-section x volume that
   *          maintains consistency with the MOC solution. Generall, this will
   *          not change the diagonal element at all since CMFD should be
   *          consistent with MOC. However, if negative fluxes are corrected to
   *          zero after the MOC transport sweep, there will be an inconsistency.
   *          This function modifies sigma-t so that there is consistency with
   *          the altered solution.
   * @param cmfd_cell The cmfd cell of the element to adjust
   * @param group The cmfd group of the element to adjust
   */
  void Cmfd::enforceBalanceOnDiagonal(int cmfd_cell, int group)
  {

    /* Initialize tallies */
    // unused
    // Material* material = _materials[cmfd_cell];
    // double cmfd_volume = _volumes->getValue(cmfd_cell, 0);

    /* Loop over FSRs in CMFD cell to tally the total neutron source */
    double moc_source = 0.0;
    for (size_t j = 0; j < _cell_fsrs.at(cmfd_cell).size(); j++)
    {

      long fsr_id = _cell_fsrs.at(cmfd_cell).at(j);
      FP_PRECISION volume = _FSR_volumes[fsr_id];

      /* Loop over MOC energy groups within this CMFD coarse group */
      for (int h = _group_indices[group]; h < _group_indices[group + 1]; h++)
        moc_source += 4 * M_PI * volume *
                      _FSR_sources[fsr_id * _num_moc_groups + h]; // MOC源项=4π*FSR体积*FSR源项
    }

    if (fabs(moc_source) < FLT_EPSILON)
      moc_source = 1e-20;

    /* Compute updated value */
    double flux = _old_flux->getValue(cmfd_cell, group);
    CMFD_PRECISION net_current = _net_currents->getValue(cmfd_cell, group);
    CMFD_PRECISION updated_value = (moc_source - net_current) / flux; //(MOC源项-净流项)/老通量

    if (updated_value < 0.0)
      log::ferror("Negative Total XS of %6.4f computed in CMFD rebalance",
                  updated_value);

    /* Update the diagonal element */
    _A->setValue(cmfd_cell, group, cmfd_cell, group, updated_value); // 对角线设置为这个值,而不是在之前的基础上增加
  }

  /**
   * @brief Rebalances the total cross section to be consistent with the MOC
   *        solution on every sweep
   * @param balance_sigma_t Wheter to compute the rebalanced total cross-section
   */
  void Cmfd::rebalanceSigmaT(bool balance_sigma_t)
  {
    _balance_sigma_t = balance_sigma_t;
  }

  /**
   * @brief Returns a flag indicating whether the sigma-t rebalance is on
   * @return A flag indicating whether the rebalance is on
   */
  bool Cmfd::isSigmaTRebalanceOn()
  {
    return _balance_sigma_t;
  }

  /**
   * @brief Get flag indicating whether to update the MOC flux.
   * @return Flag saying whether to update MOC flux.
   */
  bool Cmfd::isFluxUpdateOn()
  {
    return _flux_update_on;
  }

  /**
   * @brief Set flag indicating whether to use FSR centroids to update
   *        the MOC flux.
   * @param centroid_update_on Flag saying whether to use centroids to
   *        update MOC flux.
   */
  void Cmfd::setCentroidUpdateOn(bool centroid_update_on)
  {
    _centroid_update_on = centroid_update_on;
  }

  /**
   * @brief Get flag indicating whether to use FSR centroids to update
   *        the MOC flux.
   * @return Flag saying whether to use centroids to update MOC flux.
   */
  bool Cmfd::isCentroidUpdateOn()
  {
    return _centroid_update_on;
  }

  /**
   * @brief Sets the threshold for CMFD source convergence (>0)
   * @param the threshold for source convergence
   */
  void Cmfd::setSourceConvergenceThreshold(double source_thresh)
  {

    if (source_thresh <= 0.0)
      log::ferror("Unable to set the CMFD source convergence threshold to"
                  " %f since the threshold must be positive.",
                  source_thresh);

    _source_convergence_threshold = source_thresh;
  }

  /**
   * @brief Sets the PolarQuad object in use by the MOC Solver.
   * @param quadrature a PolarQuad object pointer from the Solver
   */
  void Cmfd::setQuadrature(QuadraturePtr quadrature)
  {
    _quadrature = quadrature;
    _num_polar = quadrature->getNumPolarAngles();
    _num_azim = quadrature->getNumAzimAngles();
  }

  /**
   * @brief Generate the k-nearest neighbor CMFD cell stencil for each FSR.为每个FSR生成k个最近邻CMFD单元模板。
   * @details This method finds the k-nearest CMFD cell stencil for each FSR
   *          and saves the stencil, ordered from the closest-to-furthest
   *          CMFD cell, in the _k_nearest_stencils map. The stencil of cells
   *          surrounding the current cell is defined as:
   *
   *          RecLattice         6 7 8
   *                             3 4 5
   *                             0 1 2
   *
   *          where 4 is the given CMFD cell. If the cell is on the edge or corner
   *          of the geometry and there are less than k nearest neighbor cells,
   *          k is reduced to the number of neighbor cells for that instance.
   */
  void Cmfd::generateKNearestStencils()
  {
    std::vector<std::pair<int, double>>::iterator stencil_iter;
    std::vector<long>::iterator fsr_iter;
    Point *centroid;
    long fsr_id;

    if (_centroid_update_on)
    {
      /* Number of cells in stencil 一个模具有九个CMFD，分别代表左右上下角、左右上下边以及当前的CMFD，每个数字对应哪一个在getDistanceToCentroid函数中给出*/
      int num_cells_in_stencil = 9;

      /* Loop over mesh cells  遍历当前domain 中的所有CMFD网格*/
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        int global_ind = getGlobalCMFDCell(i);

        /* Loop over FRSs in mesh cell 遍历一个CMFD网格的所有FSR */
        for (fsr_iter = _cell_fsrs.at(i).begin();
             fsr_iter != _cell_fsrs.at(i).end(); ++fsr_iter)
        {

          fsr_id = *fsr_iter; // 对迭代器进行解引用,相当于赋值操作

          /* Get centroid */
          centroid = _geometry->getFSRCentroid(fsr_id);

          /* Create new stencil 创建一个模具，它本身是一个Map，大小是一个CMFD中所有的FSR数量，用来存储每个fsr的k个最近模板的映射，映射值为一个二维Vector容器*/
          _k_nearest_stencils[fsr_id] =
              std::vector<std::pair<int, double>>();

          /* Get distance to all cells that touch current cell 遍历9个（如上图的九个CMFD的位置）来计算每个FSR质心到该FSR所在的CMFD网格的距离*/
          for (int j = 0; j < num_cells_in_stencil; j++)
          {
            _k_nearest_stencils[fsr_id]
                .push_back(std::make_pair<int, double>(int(j), getDistanceToCentroid(centroid, global_ind, j)));
          }

          /* Sort the distances 按距离长度排个序从小到大，以便于后续取出最小的几个*/
          std::sort(_k_nearest_stencils[fsr_id].begin(),
                    _k_nearest_stencils[fsr_id].end(), stencilCompare);

          /* Remove ghost cells that are outside the geometry boundaries 移除不在几何边界的CMFD，如果不在距离会被设置为max，下面就是将为距离为max的移除出_k_nearest_stencils*/
          stencil_iter = _k_nearest_stencils[fsr_id].begin();
          while (stencil_iter != _k_nearest_stencils[fsr_id].end())
          {
            if (stencil_iter->second == std::numeric_limits<double>::max())
              stencil_iter = _k_nearest_stencils[fsr_id].erase(stencil_iter++);
            else
              ++stencil_iter;
          }

          /* Resize stencil to be of size <= _k_nearest 截取模具（给定下标了，它就是一个二维Vector容器）的前_k_nearest
          （它的作用其实就是用于更新MOC通量的CMFD网格数量，初始化为3）本来有9个，截取为前3个，也是距离最近的前3个*/
          _k_nearest_stencils[fsr_id].resize(std::min(_k_nearest, int(_k_nearest_stencils[fsr_id].size())));
        }
      }

      /* Precompute (1.0 - cell distance / total distance) of each FSR centroid to
       * its k-nearest CMFD cells 根据距离计算更新moc时的权重（占比）*/
      double total_distance;
      for (long i = 0; i < _num_FSRs; i++)
      {
        total_distance = 1.e-10;

        /* Compute the total distance of each FSR centroid to its k-nearest CMFD
         * cells   计算每个FSR质心到各个CMFD前k个最短距离的总长度，一个求和操作*/
        for (stencil_iter = _k_nearest_stencils[i].begin();
             stencil_iter < _k_nearest_stencils[i].end(); ++stencil_iter)
          total_distance += stencil_iter->second;

        /* Reset the second stencil value to
        * (1.0 - cell_distance / total_distance) 把这个前K个模具的的值改为：1-（本来的长度/刚刚求的总长度）
        这个操作其实就是算每个相邻CMFD的对这个FSR通量更新的贡献，因为（本来的长度/刚刚求的总长度）
        这个数据为长度在总长度中的占比，1-这个数据，表示距离越小（距离FSR更近的CMFD），贡献越大*/
        for (stencil_iter = _k_nearest_stencils[i].begin();
             stencil_iter < _k_nearest_stencils[i].end(); ++stencil_iter)
          stencil_iter->second = 1.0 - stencil_iter->second / total_distance;
      }
    }

    /* Compute axial quadratic interpolation values if requested
    上述代码的作用是FSR的通量更新的二维CMFD模具，这部分计算Z轴的轴向二次插值，FSR通量更新比率也使用轴向插值 对应CMFD文档B.33部分*/
    if (_use_axial_interpolation && _local_num_z >= 3)
    { // 选择插值器来保持三个节点中的每一个节点(当前单元和两个轴向邻居)中的平均通量,Z轴向CMFD的数量不少于3

      /* Initialize axial quadratic interpolant values 初始化轴向二次插值*/
      _axial_interpolants.resize(_num_FSRs);
      for (long r = 0; r < _num_FSRs; r++)
      {
        _axial_interpolants.at(r) = new double[3](); // 分配一个包含 3 个 double 类型元素的数组，并将其初始化为零
      }

      /* Calculate common factors 计算常用系数*/
      double dz = _cell_width_z;
      // unused
      // double dz_2 = _cell_width_z * _cell_width_z;

      /* Loop over mesh cells 循环所有的CMFD网格*/
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        /* Calculate the CMFD cell z-coordinate */
        int z_ind = i / (_local_num_x * _local_num_y);            // 局部Z坐标索引
        double z_cmfd = (z_ind + 0.5) * dz + _lattice->getMinZ(); // 求的是该CMFD网格在Z轴的中心坐标，也就是质心高度
        if (_domain_communicator != NULL)
          z_cmfd += _domain_communicator->_domain_idx_z * dz * _local_num_z; // 获取全局高度

        /* Loop over FRSs in mesh cell 循环一个CMFD中的所有FSR*/
        int num_fissionable_FSRs = 0;
        for (fsr_iter = _cell_fsrs.at(i).begin();
             fsr_iter != _cell_fsrs.at(i).end(); ++fsr_iter)
        {

          /* Get centroid and calculate relative z-coordinate */
          fsr_id = *fsr_iter;
          Point *centroid = _geometry->getFSRCentroid(fsr_id); // 求得该FSR的质心
          double zc = (centroid->getZ() - z_cmfd) / dz;        // ZC对应文档中的（z-zjc），多除了个高度dz,这里暂时对应不上！！！！！！
          if (std::abs(zc) > 0.5)
            log::ferror("Found FSR %d with z-centroid offset in z "
                        "from CMFD cell %d by %6.4f, whereas the CMFD z-spacing"
                        " is %6.4f. Coordinates: (%6.4f, %6.4f, %6.4f), cmfd z: "
                        "%6.4f",
                        fsr_id, i, zc * dz, dz, centroid->getX(),
                        centroid->getY(), centroid->getZ(), z_cmfd);

          /* Check that the CMFD cell is not an end cell 如果该CMFD网格在z轴上为第一个或最后一个，对FSR质心和它所在的CMFD网格质心
          的Z轴距离差所占CMFD网格高度的比重（也就是ZC）进行-或+1*/
          if (z_ind == 0)
            zc -= 1.0; // 是为了将 zc 向下调整一个单位，使得在底部边界处进行插值时，能够正确地考虑底部邻接的虚拟单元
          else if (z_ind == _local_num_z - 1)
            zc += 1.0; // 是为了将 zc 向上调整一个单位，使得在顶部边界处进行插值时，能够正确地考虑顶部邻接的虚拟单元

          /* Calculate components for quadratic interpolation  计算二次插值的分量 */
          _axial_interpolants.at(fsr_id)[0] = zc * zc / 2.0 - zc / 2.0 - 1.0 / 24.0; // 这个公式也没对应上，第二项多除了个2！！！！！！！
          _axial_interpolants.at(fsr_id)[1] = -zc * zc + 26.0 / 24.0;
          _axial_interpolants.at(fsr_id)[2] = zc * zc / 2.0 + zc / 2.0 - 1.0 / 24.0; // 这里也是没对应上

          /* Set zero axial prolongation for cells with no fissionalbe material 没有裂变物质的FSR设置零轴向延长*/
          if (_FSR_materials[fsr_id]->isFissionable())
            num_fissionable_FSRs++; // 统计了FSR中材料为裂变材料的FSR数量，但这个数据好像没有用到，是不是如果FSR材料不是裂变材料就不进行轴向二次插值（个人猜测）
        }
      }
    }
  }

  /**
   * @brief Generate the k-nearest neighbor CMFD cell stencil for each FSR.
   * @details This method finds the k-nearest CMFD cell stencil for each FSR
   *          and saves the stencil, ordered from the closest-to-furthest
   *          CMFD cell, in the _k_nearest_stencils map. The stencil of cells
   *          surrounding the current cell is defined as:
   *
   *          HexLattice         5 6
   *                             2 3 4
   *                               0 1
   *
   *          where 3 is the given CMFD cell. If the cell is on the edge or corner
   *          of the geometry and there are less than k nearest neighbor cells,
   *          k is reduced to the number of neighbor cells for that instance.
   */
  void Cmfd::generateKNearestStencilsHex()
  {
    std::vector<std::pair<int, double>>::iterator stencil_iter;
    std::vector<long>::iterator fsr_iter;
    Point *centroid;
    long fsr_id;

    if (_centroid_update_on)
    {
      /* Number of cells in stencil */
      int num_cells_in_stencil = 7;

      /* Loop over mesh cells */
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        int global_ind = getGlobalCMFDCell(i);

        if (!_empty_fsrs_cells[global_ind])
        {

          /* Loop over FRSs in mesh cell */
          for (fsr_iter = _cell_fsrs.at(i).begin();
               fsr_iter != _cell_fsrs.at(i).end(); ++fsr_iter)
          {

            fsr_id = *fsr_iter;

            /* Get centroid */
            centroid = _geometry->getFSRCentroid(fsr_id);

            /* Create new stencil */
            _k_nearest_stencils[fsr_id] =
                std::vector<std::pair<int, double>>();

            /* Get distance to all cells that touch current cell */
            for (int j = 0; j < num_cells_in_stencil; j++)
            {
              _k_nearest_stencils[fsr_id]
                  .push_back(std::make_pair<int, double>(int(j), getDistanceToCentroidHex(centroid, global_ind, j)));
            }

            /* Sort the distances */
            std::sort(_k_nearest_stencils[fsr_id].begin(),
                      _k_nearest_stencils[fsr_id].end(), stencilCompare);

            /* Remove ghost cells that are outside the geometry boundaries */
            stencil_iter = _k_nearest_stencils[fsr_id].begin();
            while (stencil_iter != _k_nearest_stencils[fsr_id].end())
            {
              if (stencil_iter->second == std::numeric_limits<double>::max())
                stencil_iter = _k_nearest_stencils[fsr_id].erase(stencil_iter++);
              else
                ++stencil_iter;
            }

            /* Resize stencil to be of size <= _k_nearest */
            _k_nearest_stencils[fsr_id].resize(std::min(_k_nearest, int(_k_nearest_stencils[fsr_id].size())));
          }
        }
      }

      /* Precompute (1.0 - cell distance / total distance) of each FSR centroid to
       * its k-nearest CMFD cells */
      double total_distance;
      for (long i = 0; i < _num_FSRs; i++)
      {
        total_distance = 1.e-10;

        /* Compute the total distance of each FSR centroid to its k-nearest CMFD
         * cells */
        for (stencil_iter = _k_nearest_stencils[i].begin();
             stencil_iter < _k_nearest_stencils[i].end(); ++stencil_iter)
          total_distance += stencil_iter->second;

        /* Reset the second stencil value to
         * (1.0 - cell_distance / total_distance) */
        for (stencil_iter = _k_nearest_stencils[i].begin();
             stencil_iter < _k_nearest_stencils[i].end(); ++stencil_iter)
          stencil_iter->second = 1.0 - stencil_iter->second / total_distance;
      }
    }

    /* Compute axial quadratic interpolation values if requested */
    if (_use_axial_interpolation && _local_num_z >= 3)
    {

      /* Initialize axial quadratic interpolant values */
      _axial_interpolants.resize(_num_FSRs);
      for (long r = 0; r < _num_FSRs; r++)
      {
        _axial_interpolants.at(r) = new double[3]();
      }

      /* Calculate common factors */
      double dz = _cell_width_z;
      // unused
      // double dz_2 = _cell_width_z * _cell_width_z;

      /* Loop over mesh cells */
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        if (!_empty_fsrs_cells[i])
        {
          /* Calculate the CMFD cell z-coordinate */
          int z_ind = i / (_local_num_x * _local_num_y);
          double z_cmfd = (z_ind + 0.5) * dz + _lattice->getMinZ();
          if (_domain_communicator != NULL)
            z_cmfd += _domain_communicator->_domain_idx_z * dz * _local_num_z;

          /* Loop over FRSs in mesh cell */
          int num_fissionable_FSRs = 0;
          for (fsr_iter = _cell_fsrs.at(i).begin();
               fsr_iter != _cell_fsrs.at(i).end(); ++fsr_iter)
          {

            /* Get centroid and calculate relative z-coordinate */
            fsr_id = *fsr_iter;
            Point *centroid = _geometry->getFSRCentroid(fsr_id);
            double zc = (centroid->getZ() - z_cmfd) / dz;
            if (std::abs(zc) > 0.5)
              log::ferror("Found FSR %d with z-centroid offset in z "
                          "from CMFD cell %d by %6.4f, whereas the CMFD z-spacing"
                          " is %6.4f. Coordinates: (%6.4f, %6.4f, %6.4f), cmfd z: "
                          "%6.4f",
                          fsr_id, i, zc * dz, dz, centroid->getX(),
                          centroid->getY(), centroid->getZ(), z_cmfd);

            /* Check that the CMFD cell is not an end cell */
            if (z_ind == 0)
              zc -= 1.0;
            else if (z_ind == _local_num_z - 1)
              zc += 1.0;

            /* Calculate components for quadratic interpolation */
            _axial_interpolants.at(fsr_id)[0] = zc * zc / 2.0 - zc / 2.0 - 1.0 / 24.0;
            _axial_interpolants.at(fsr_id)[1] = -zc * zc + 26.0 / 24.0;
            _axial_interpolants.at(fsr_id)[2] = zc * zc / 2.0 + zc / 2.0 - 1.0 / 24.0;

            /* Set zero axial prolongation for cells with no fissionalbe material */
            if (_FSR_materials[fsr_id]->isFissionable())
              num_fissionable_FSRs++;
          }
        }
      }
    }
  }

  /**
   * @brief Get the ID of the Mesh cell given a stencil ID and Mesh cell ID.
   * @details The stencil of cells surrounding the current cell is defined as:
   *
   *                             6 7 8
   *                             3 4 5
   *                             0 1 2
   *
   * @param cell_id Current Mesh cell ID
   * @param stencil_id CMFD cell stencil ID
   * @return Neighboring CMFD cell ID
   */
  int Cmfd::getCellByStencil(int cell_id, int stencil_id)
  {

    int cell_next_id = -1;
    int x = (cell_id % (_local_num_x * _local_num_y)) % _local_num_x;
    int y = (cell_id % (_local_num_x * _local_num_y)) / _local_num_x;

    if (stencil_id == 0)
    {
      if (x != 0 && y != 0)
        cell_next_id = cell_id - _local_num_x - 1;
    }
    else if (stencil_id == 1)
    {
      if (y != 0)
        cell_next_id = cell_id - _local_num_x;
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        cell_next_id = cell_id + _local_num_x * (_local_num_y - 1);
    }
    else if (stencil_id == 2)
    {
      if (x != _local_num_x - 1 && y != 0)
        cell_next_id = cell_id - _local_num_x + 1;
    }
    else if (stencil_id == 3)
    {
      if (x != 0)
        cell_next_id = cell_id - 1;
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        cell_next_id = cell_id + (_local_num_x - 1);
    }
    else if (stencil_id == 4)
    {
      cell_next_id = cell_id;
    }
    else if (stencil_id == 5)
    {
      if (x != _local_num_x - 1)
        cell_next_id = cell_id + 1;
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        cell_next_id = cell_id - (_local_num_x - 1);
    }
    else if (stencil_id == 6)
    {
      if (x != 0 && y != _local_num_y - 1)
        cell_next_id = cell_id + _local_num_x - 1;
    }
    else if (stencil_id == 7)
    {
      if (y != _local_num_y - 1)
        cell_next_id = cell_id + _local_num_x;
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        cell_next_id = cell_id - _local_num_x * (_local_num_y - 1);
    }
    else if (stencil_id == 8)
    {
      if (x != _local_num_x - 1 && y != _local_num_y - 1)
        cell_next_id = cell_id + _local_num_x + 1;
    }

    return cell_next_id;
  }

  /**
   * @brief Get the ID of the Mesh cell given a stencil ID and Mesh cell ID.
   * @details The stencil of cells surrounding the current cell is defined as:
   *
   *                             5 6
   *                             2 3 4
   *                               0 1
   *
   * @param cell_id Current Mesh cell ID
   * @param stencil_id CMFD cell stencil ID
   * @return Neighboring CMFD cell ID
   */
  int Cmfd::getHexCellByStencil(int cell_id, int stencil_id)
  {

    int cell_next_id = -1;

    int x = (cell_id % (_local_num_x * _local_num_y)) % _local_num_x;
    int y = (cell_id % (_local_num_x * _local_num_y)) / _local_num_x;

    if (stencil_id == 0)
    {
      if (CellinXYBoundaryWithStencil(cell_id, stencil_id))
        cell_next_id = cell_id - _local_num_x;
      else if (_boundaries[HEX_SURFACE_DELTA_MIN] == PERIODIC)
        cell_next_id = cell_id + _local_num_x * (_local_num_y - 1);
    }
    else if (stencil_id == 1)
    {
      if (CellinXYBoundaryWithStencil(cell_id, stencil_id))
        cell_next_id = cell_id - _local_num_x + 1;
      else if (_boundaries[HEX_SURFACE_GAMMA_MAX] == PERIODIC)
        cell_next_id = cell_id + _local_num_x * (_local_num_y - 2) + 1;
    }
    else if (stencil_id == 2)
    {
      if (CellinXYBoundaryWithStencil(cell_id, stencil_id))
        cell_next_id = cell_id - 1;
      else if (_boundaries[HEX_SURFACE_BETA_MIN] == PERIODIC)
        cell_next_id = cell_id + (_local_num_x - 1);
    }
    else if (stencil_id == 3)
    {
      cell_next_id = cell_id;
    }
    else if (stencil_id == 4)
    {
      if (CellinXYBoundaryWithStencil(cell_id, stencil_id))
        cell_next_id = cell_id + 1;
      else if (_boundaries[HEX_SURFACE_BETA_MAX] == PERIODIC)
        cell_next_id = cell_id - (_local_num_x - 1);
    }
    else if (stencil_id == 5)
    {
      if (CellinXYBoundaryWithStencil(cell_id, stencil_id))
        cell_next_id = cell_id + _local_num_x - 1;
      else if (_boundaries[HEX_SURFACE_GAMMA_MIN] == PERIODIC)
        cell_next_id = cell_id - _local_num_x * (_local_num_y - 2) - 1;
    }
    else if (stencil_id == 6)
    {
      if (CellinXYBoundaryWithStencil(cell_id, stencil_id))
        cell_next_id = cell_id + _local_num_x;
      else if (_boundaries[HEX_SURFACE_DELTA_MAX] == PERIODIC)
        cell_next_id = cell_id - _local_num_x * (_local_num_y - 1);
    }

    return cell_next_id;
  }

  int Cmfd::getHexCellByStencilPre(int cell_id, int stencil_id)
  {

    int cell_next_id = -1;

    int x = (cell_id % (_local_num_x * _local_num_y)) % _local_num_x;
    int y = (cell_id % (_local_num_x * _local_num_y)) / _local_num_x;

    /* Function from HexLattice "isValidIndex(int index)"
    int nx = 2*_num_r - 1;
    int nxy = nx*nx;  // nx == ny
    int z = cell_id / nxy;
    int y = (cell_id - nxy * z) / nx;
    int x = cell_id - nxy * z - nx * y;
    */

    if (stencil_id == 0)
    {
      if (y != 0)
        cell_next_id = cell_id - _local_num_x;
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        cell_next_id = cell_id + _local_num_x * (_local_num_y - 1);
    }
    else if (stencil_id == 1)
    {
      if (x != _local_num_x - 1 && y != 0)
        cell_next_id = cell_id - _local_num_x + 1;
    }
    else if (stencil_id == 2)
    {
      if (x != 0)
        cell_next_id = cell_id - 1;
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        cell_next_id = cell_id + (_local_num_x - 1);
    }
    else if (stencil_id == 3)
    {
      cell_next_id = cell_id;
    }
    else if (stencil_id == 4)
    {
      if (x != _local_num_x - 1)
        cell_next_id = cell_id + 1;
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        cell_next_id = cell_id - (_local_num_x - 1);
    }
    else if (stencil_id == 5)
    {
      if (x != 0 && y != _local_num_y - 1)
        cell_next_id = cell_id + _local_num_x - 1;
    }
    else if (stencil_id == 6)
    {
      if (y != _local_num_y - 1)
        cell_next_id = cell_id + _local_num_x;
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        cell_next_id = cell_id - _local_num_x * (_local_num_y - 1);
    }

    return cell_next_id;
  }

  /**
   * @brief Get the ratio used to update the FSR flux after converging CMFD.
   * @details This method takes in a CMFD cell, a MOC energy group, and a FSR
   *          and returns the ratio used to update the FSR flux. There are two
   *          methods that can be used to update the flux, conventional and
   *          k-nearest centroid updating. The k-nearest centroid updating uses
   *          the k-nearest cells (with k between 1 and 9) of the current CMFD
   *          cell and the 8 neighboring CMFD cells. The stencil of cells
   *          surrounding the current cell is defined as:
   *
   *                             6 7 8
   *                             3 4 5
   *                             0 1 2
   *
   *          where 4 is the given CMFD cell. If the cell is on the edge or corner
   *          of the geometry and there are less than k nearest neighbor cells,
   *          k is reduced to the number of neighbor cells for that instance.
   * @param cell_id The CMFD cell ID containing the FSR.
   * @param group The CMFD energy group being updated.
   * @param fsr The fsr being updated.
   * @return the ratio used to update the FSR flux.
   * 考虑FSR在CMFD网格内的位置
   */
  CMFD_PRECISION Cmfd::getUpdateRatio(int cell_id, int group, int fsr)
  {

    CMFD_PRECISION ratio = 0.0;
    std::vector<std::pair<int, double>>::iterator iter;
    int cell_next_id;

    if (_centroid_update_on)
    {
      /* ========================================
       * 质心更新方法：使用k近邻CMFD网格
       * ========================================
       * 原理：考虑FSR在CMFD网格内的实际位置（质心）
       * 作用：基于距离加权，使用多个CMFD网格的比率
       */

      /* Compute the ratio for all the surrounding cells */
      for (iter = _k_nearest_stencils[fsr].begin();
           iter != _k_nearest_stencils[fsr].end(); ++iter)
      {
        if (iter->first != 4)
        { // 对于每个质心（不包括中心单元，即 iter->first != 4），计算与这些质心相关的 CMFD 单元的 ID
          cell_next_id = getCellByStencil(getLocalCMFDCell(cell_id), iter->first);

          ratio += iter->second * getFluxRatio(cell_next_id, group, fsr);
        }
      }

      /* INTERNAL */
      if (_k_nearest_stencils[fsr].size() == 1)
        ratio += getFluxRatio(cell_id, group, fsr);
      else
      {
        ratio += _k_nearest_stencils[fsr][0].second *
                 getFluxRatio(cell_id, group, fsr);
        ratio /= (_k_nearest_stencils[fsr].size() - 1);
      }
    }
    else
      /* ========================================
       * 基础方法：直接使用CMFD网格的比率
       * ========================================
       */
      ratio = getFluxRatio(cell_id, group, fsr);

    return ratio;
  }

  /**
   * @brief Get the Hex Lattice ratio used to update the FSR flux after converging CMFD.
   * @details This method takes in a CMFD cell, a MOC energy group, and a FSR
   *          and returns the ratio used to update the FSR flux. There are two
   *          methods that can be used to update the flux, conventional and
   *          k-nearest centroid updating. The k-nearest centroid updating uses
   *          the k-nearest cells (with k between 1 and 7) of the current CMFD
   *          cell and the 6 neighboring CMFD cells. The stencil of cells
   *          surrounding the current cell is defined as:
   *
   *                             5 6
   *                             2 3 4
   *                               0 1
   *
   *          where 3 is the given CMFD cell. If the cell is on the edge or corner
   *          of the geometry and there are less than k nearest neighbor cells,
   *          k is reduced to the number of neighbor cells for that instance.
   * @param cell_id The CMFD cell ID containing the FSR.
   * @param group The CMFD energy group being updated.
   * @param fsr The fsr being updated.
   * @return the ratio used to update the FSR flux.
   */
  CMFD_PRECISION Cmfd::getHexUpdateRatio(int cell_id, int group, int fsr)
  {

    CMFD_PRECISION ratio = 0.0;
    std::vector<std::pair<int, double>>::iterator iter;
    int cell_next_id;

    if (_centroid_update_on)
    {

      /* Compute the ratio for all the surrounding cells */
      for (iter = _k_nearest_stencils[fsr].begin();
           iter != _k_nearest_stencils[fsr].end(); ++iter)
      {

        if (iter->first != 3)
        {
          // getHexCellByStencilPre() need modified
          cell_next_id = getHexCellByStencil(cell_id, iter->first);
          if (cell_next_id != -1)
          {
            if (!_empty_fsrs_cells[cell_next_id])
            {
              ratio += iter->second * getHexFluxRatio(cell_next_id, group, fsr);
            }
          }
        }
      }

      /* INTERNAL */
      if (_k_nearest_stencils[fsr].size() == 1)
      {
        ratio += getHexFluxRatio(cell_id, group, fsr);
      }
      else
      {
        ratio += _k_nearest_stencils[fsr][0].second *
                 getHexFluxRatio(cell_id, group, fsr);
        ratio /= (_k_nearest_stencils[fsr].size() - 1);
      }
    }
    else
      ratio = getHexFluxRatio(cell_id, group, fsr);

    return ratio;
  }

  /**
   * @brief Retreives the ratio of pre- and post- CMFD solve fluxes
   * @details The CMFD flux ratio is returned for the given FSR. A quadratic
   *          axial interpolant is used to estimate the value at the FSR.
   * @param cell_id The CMFD cell ID containing the FSR.
   * @param group The CMFD energy group being updated.
   * @param fsr The fsr being updated.
   * @return the ratio of CMFD fluxes
   *  * 关键点：
   * - 可以支持轴向插值（3D问题）
   * - 如果旧通量为零，比率设为1.0（避免除零）
   */
  CMFD_PRECISION Cmfd::getFluxRatio(int cell_id, int group, int fsr)
  {

    double *interpolants;
    double ratio = 1.0;
    /* 如果启用轴向插值（3D问题） */
    if (_use_axial_interpolation)
      interpolants = _axial_interpolants.at(fsr);
    if (_use_axial_interpolation && _local_num_z >= 3 &&
        (fabs(interpolants[0]) > FLT_EPSILON ||
         fabs(interpolants[2]) > FLT_EPSILON))
    {
      int z_ind = cell_id / (_local_num_x * _local_num_y);
      int cell_mid = cell_id;
      if (z_ind == 0)
        cell_mid += _local_num_x * _local_num_y;
      else if (z_ind == _local_num_z - 1)
        cell_mid -= _local_num_x * _local_num_y;

      int cell_prev = cell_mid - _local_num_x * _local_num_y;
      int cell_next = cell_mid + _local_num_x * _local_num_y;

      double old_flux_prev = _old_flux->getValue(cell_prev, group);
      double new_flux_prev = _new_flux->getValue(cell_prev, group);

      double old_flux_next = _old_flux->getValue(cell_next, group);
      double new_flux_next = _new_flux->getValue(cell_next, group);

      double old_flux_mid = _old_flux->getValue(cell_mid, group);
      double new_flux_mid = _new_flux->getValue(cell_mid, group);

      double old_flux = interpolants[0] * old_flux_prev +
                        interpolants[1] * old_flux_mid +
                        interpolants[2] * old_flux_next;
      double new_flux = interpolants[0] * new_flux_prev +
                        interpolants[1] * new_flux_mid +
                        interpolants[2] * new_flux_next;

      if (fabs(old_flux) > FLT_EPSILON)
        ratio = new_flux / old_flux;

      if (ratio < 0)
      {

        /* Try a linear interpolation */
        double zc_2 = 26.0 / 24.0 - interpolants[1];
        if (zc_2 < 0.0)
          zc_2 = 0.0;
        double zc = sqrt(zc_2);
        if (z_ind < _num_z / 2)
        {
          old_flux = zc * (old_flux_mid - old_flux_prev) + old_flux_mid;
          new_flux = zc * (new_flux_mid - new_flux_prev) + new_flux_mid;
          if (fabs(old_flux) > FLT_EPSILON)
            ratio = new_flux / old_flux;
          else
            ratio = 0;
        }
        else
        {
          old_flux = zc * (old_flux_next - old_flux_mid) + old_flux_mid;
          new_flux = zc * (new_flux_next - new_flux_mid) + new_flux_mid;
          if (fabs(old_flux) > FLT_EPSILON)
            ratio = new_flux / old_flux;
          else
            ratio = 0;
        }

        /* Fallback: using the cell average flux ratio */
        if (ratio < 0)
        {
          if (fabs(_old_flux->getValue(cell_id, group)) > FLT_EPSILON)
            ratio = _new_flux->getValue(cell_id, group) /
                    _old_flux->getValue(cell_id, group);
          else
            ratio = 0.0;
        }
      }

      return ratio;
    }
    else
    {
      /* 直接计算比率 */
      if (fabs(_old_flux->getValue(cell_id, group)) > FLT_EPSILON)
        return _new_flux->getValue(cell_id, group) /
               _old_flux->getValue(cell_id, group);
      else
        return 0.0;
    }
  }

  CMFD_PRECISION Cmfd::getHexFluxRatio(int cell_id, int group, int fsr)
  {

    double *interpolants;
    double ratio = 1.0;
    if (_use_axial_interpolation)
      interpolants = _axial_interpolants.at(fsr);
    if (_use_axial_interpolation && _local_num_z >= 3 &&
        (fabs(interpolants[0]) > FLT_EPSILON ||
         fabs(interpolants[2]) > FLT_EPSILON))
    {
      int z_ind = cell_id / (_local_num_x * _local_num_y);
      int cell_mid = cell_id;
      if (z_ind == 0)
        cell_mid += _local_num_x * _local_num_y;
      else if (z_ind == _local_num_z - 1)
        cell_mid -= _local_num_x * _local_num_y;

      int cell_prev = cell_mid - _local_num_x * _local_num_y;
      int cell_next = cell_mid + _local_num_x * _local_num_y;

      double old_flux_prev;
      double new_flux_prev;
      double old_flux_next;
      double new_flux_next;
      double old_flux_mid;
      double new_flux_mid;

      if (!_empty_fsrs_cells[cell_prev])
      {
        cell_prev = _logical_actual_map[cell_prev];
        double old_flux_prev = _old_flux->getValue(cell_prev, group);
        double new_flux_prev = _new_flux->getValue(cell_prev, group);
      }

      if (!_empty_fsrs_cells[cell_next])
      {
        cell_next = _logical_actual_map[cell_next];
        double old_flux_next = _old_flux->getValue(cell_next, group);
        double new_flux_next = _new_flux->getValue(cell_next, group);
      }

      if (!_empty_fsrs_cells[cell_mid])
      {
        cell_mid = _logical_actual_map[cell_mid];
        double old_flux_mid = _old_flux->getValue(cell_mid, group);
        double new_flux_mid = _new_flux->getValue(cell_mid, group);
      }

      double old_flux = interpolants[0] * old_flux_prev +
                        interpolants[1] * old_flux_mid +
                        interpolants[2] * old_flux_next;
      double new_flux = interpolants[0] * new_flux_prev +
                        interpolants[1] * new_flux_mid +
                        interpolants[2] * new_flux_next;

      if (fabs(old_flux) > FLT_EPSILON)
        ratio = new_flux / old_flux;

      if (ratio < 0)
      {

        /* Try a linear interpolation */
        double zc_2 = 26.0 / 24.0 - interpolants[1];
        if (zc_2 < 0.0)
          zc_2 = 0.0;
        double zc = sqrt(zc_2);
        if (z_ind < _num_z / 2)
        {
          old_flux = zc * (old_flux_mid - old_flux_prev) + old_flux_mid;
          new_flux = zc * (new_flux_mid - new_flux_prev) + new_flux_mid;
          if (fabs(old_flux) > FLT_EPSILON)
            ratio = new_flux / old_flux;
          else
            ratio = 0;
        }
        else
        {
          old_flux = zc * (old_flux_next - old_flux_mid) + old_flux_mid;
          new_flux = zc * (new_flux_next - new_flux_mid) + new_flux_mid;
          if (fabs(old_flux) > FLT_EPSILON)
            ratio = new_flux / old_flux;
          else
            ratio = 0;
        }

        /* Fallback: using the cell average flux ratio */
        if (ratio < 0)
        {
          cell_id = _logical_actual_map[cell_id];
          if (fabs(_old_flux->getValue(cell_id, group)) > FLT_EPSILON)
            ratio = _new_flux->getValue(cell_id, group) /
                    _old_flux->getValue(cell_id, group);
          else
            ratio = 0.0;
        }
      }

      return ratio;
    }
    else
    {
      cell_id = _logical_actual_map[cell_id];
      if (fabs(_old_flux->getValue(cell_id, group)) > FLT_EPSILON)
        return _new_flux->getValue(cell_id, group) /
               _old_flux->getValue(cell_id, group);
      else
        return 0.0;
    }
  }

  /**
   * @brief Get the distances from an FSR centroid to a given CMFD cell.
   * @details This method takes in a FSR centroid, a CMFD cell, and a stencil index
   *          to a cell located in the 9-point stencil encompassing the CMFD
   *          cell an all its possible neighbors. The CMFD cell stencil is:
   *
   *                             6 7 8
   *                             3 4 5
   *                             0 1 2
   *
   *          where 4 is the given CMFD cell. If a CMFD edge or corner cells is
   *          given and the stencil indexed cell lies outside the geometry, the
   *          maximum allowable double value is returned.
   * @param centroid The numerical centroid an FSR in the cell.
   * @param cell_id The CMFD cell containing the FSR.
   * @param stencil_index The index of the cell in the stencil that we want to
   *        get the distance from.
   * @return the distance from the CMFD cell centroid to the FSR centroid.
   */
  double Cmfd::getDistanceToCentroid(Point *centroid, int cell_id,
                                     int stencil_index)
  {

    int x = (cell_id % (_num_x * _num_y)) % _num_x; // 全局 x 方向的索引
    int y = (cell_id % (_num_x * _num_y)) / _num_x;
    double dist_x, dist_y;
    bool found = false;
    double centroid_x = centroid->getX();
    double centroid_y = centroid->getY();

    /* LOWER LEFT CORNER 左下 0*/
    if (x > 0 && y > 0 && stencil_index == 0)
    {                                                                                // x > 0 && y > 0  保证存在左下的这个CMFD 以下同
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x - 0.5) * _cell_width_x), 2.0); // x方向上CMFD质心与FSR质心距离
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y - 0.5) * _cell_width_y), 2.0); // y方向上CMFD质心与FSR质心距离
      found = true;
    }

    /* BOTTOM SIDE 下 1*/
    else if (y > 0 && stencil_index == 1)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x + 0.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y - 0.5) * _cell_width_y), 2.0);
      found = true;
    }

    /* LOWER RIGHT CORNER 右下 2*/
    else if (x < _num_x - 1 && y > 0 && stencil_index == 2)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x + 1.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y - 0.5) * _cell_width_y), 2.0);
      found = true;
    }

    /* LEFT SIDE  左 3*/
    else if (x > 0 && stencil_index == 3)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x - 0.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y + 0.5) * _cell_width_y), 2.0);
      found = true;
    }

    /* CURRENT  当前CMFD（中间） 4*/
    else if (stencil_index == 4)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x + 0.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y + 0.5) * _cell_width_y), 2.0);
      found = true;
    }

    /* RIGHT SIDE  右 5*/
    else if (x < _num_x - 1 && stencil_index == 5)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x + 1.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y + 0.5) * _cell_width_y), 2.0);
      found = true;
    }

    /* UPPER LEFT CORNER 左上 6*/
    else if (x > 0 && y < _num_y - 1 && stencil_index == 6)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x - 0.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y + 1.5) * _cell_width_y), 2.0);
      found = true;
    }

    /* TOP SIDE  上 7*/
    else if (y < _num_y - 1 && stencil_index == 7)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x + 0.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y + 1.5) * _cell_width_y), 2.0);
      found = true;
    }

    /* UPPER RIGHT CORNER 右上 8*/
    else if (x < _num_x - 1 && y < _num_y - 1 && stencil_index == 8)
    {
      dist_x = pow(centroid_x - (-_width_x / 2.0 + (x + 1.5) * _cell_width_x), 2.0);
      dist_y = pow(centroid_y - (-_width_y / 2.0 + (y + 1.5) * _cell_width_y), 2.0);
      found = true;
    }

    if (found)
      return pow(dist_x + dist_y, 0.5);
    else
      return std::numeric_limits<double>::max();
  }

  /**
   * @brief Get the distances from an FSR centroid to a given CMFD cell.
   * @details This method takes in a FSR centroid, a CMFD cell, and a stencil index
   *          to a cell located in the 9-point stencil encompassing the CMFD
   *          cell an all its possible neighbors. The CMFD cell stencil is:
   *
   *                             5 6
   *                             2 3 4
   *                               0 1
   *
   *          where 3 is the given CMFD cell. If a CMFD edge or corner cells is
   *          given and the stencil indexed cell lies outside the geometry, the
   *          maximum allowable double value is returned.
   * @param centroid The numerical centroid an FSR in the cell.
   * @param cell_id The CMFD cell containing the FSR.
   * @param stencil_index The index of the cell in the stencil that we want to
   *        get the distance from.
   * @return the distance from the CMFD cell centroid to the FSR centroid.
   */
  double Cmfd::getDistanceToCentroidHexPre(Point *centroid, int cell_id,
                                           int stencil_index)
  {

    int x = (cell_id % (_num_x * _num_y)) % _num_x;
    int y = (cell_id % (_num_x * _num_y)) / _num_x;
    int z = (cell_id / (_num_x * _num_y));
    double dist_x, dist_y;
    bool found = false;
    double lattice_centroid_x, lattice_centroid_y;
    double centroid_x = centroid->getX();
    double centroid_y = centroid->getY();
    Point lattice_centroid;

    /* BOTTOM SIDE */
    if (y > 0 && stencil_index == 0)
    {
      lattice_centroid = _lattice->getLatticeCellCenter(x, y - 1, z);
      lattice_centroid_x = lattice_centroid.getX();
      lattice_centroid_y = lattice_centroid.getY();
      dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
      dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
      found = true;
    }

    /* LOWER RIGHT CORNER */
    else if (x < _num_x - 1 && y > 0 && stencil_index == 1)
    {
      lattice_centroid = _lattice->getLatticeCellCenter(x + 1, y - 1, z);
      lattice_centroid_x = lattice_centroid.getX();
      lattice_centroid_y = lattice_centroid.getY();
      dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
      dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
      found = true;
    }

    /* LEFT SIDE */
    else if (x > 0 && stencil_index == 2)
    {
      lattice_centroid = _lattice->getLatticeCellCenter(x - 1, y, z);
      lattice_centroid_x = lattice_centroid.getX();
      lattice_centroid_y = lattice_centroid.getY();
      dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
      dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
      found = true;
    }

    /* CURRENT */
    else if (stencil_index == 3)
    {
      lattice_centroid = _lattice->getLatticeCellCenter(x, y, z);
      lattice_centroid_x = lattice_centroid.getX();
      lattice_centroid_y = lattice_centroid.getY();
      dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
      dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
      found = true;
    }

    /* RIGHT SIDE */
    else if (x < _num_x - 1 && stencil_index == 4)
    {
      lattice_centroid = _lattice->getLatticeCellCenter(x + 1, y, z);
      lattice_centroid_x = lattice_centroid.getX();
      lattice_centroid_y = lattice_centroid.getY();
      dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
      dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
      found = true;
    }

    /* UPPER LEFT CORNER */
    else if (x > 0 && y < _num_y - 1 && stencil_index == 5)
    {
      lattice_centroid = _lattice->getLatticeCellCenter(x - 1, y + 1, z);
      lattice_centroid_x = lattice_centroid.getX();
      lattice_centroid_y = lattice_centroid.getY();
      dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
      dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
      found = true;
    }

    /* TOP SIDE */
    else if (y < _num_y - 1 && stencil_index == 6)
    {
      lattice_centroid = _lattice->getLatticeCellCenter(x, y + 1, z);
      lattice_centroid_x = lattice_centroid.getX();
      lattice_centroid_y = lattice_centroid.getY();
      dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
      dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
      found = true;
    }

    if (found)
      return pow(dist_x + dist_y, 0.5);
    else
      return std::numeric_limits<double>::max();
  }

  /**
   * @brief Get the distances from an FSR centroid to a given CMFD cell.
   * @details This method takes in a FSR centroid, a CMFD cell, and a stencil index
   *          to a cell located in the 9-point stencil encompassing the CMFD
   *          cell an all its possible neighbors. The CMFD cell stencil is:
   *
   *                             5 6
   *                             2 3 4
   *                               0 1
   *
   *          where 3 is the given CMFD cell. If a CMFD edge or corner cells is
   *          given and the stencil indexed cell lies outside the geometry, the
   *          maximum allowable double value is returned.
   * @param centroid The numerical centroid an FSR in the cell.
   * @param cell_id The CMFD cell containing the FSR.
   * @param stencil_index The index of the cell in the stencil that we want to
   *        get the distance from.
   * @return the distance from the CMFD cell centroid to the FSR centroid.
   */
  double Cmfd::getDistanceToCentroidHex(Point *centroid, int cell_id,
                                        int stencil_index)
  {

    int x = (cell_id % (_num_x * _num_y)) % _num_x;
    int y = (cell_id % (_num_x * _num_y)) / _num_x;
    int z = (cell_id / (_num_x * _num_y));
    double dist_x, dist_y;
    bool found = false;
    double lattice_centroid_x, lattice_centroid_y;
    double centroid_x = centroid->getX();
    double centroid_y = centroid->getY();
    Point lattice_centroid;

    if (CellinXYBoundary(cell_id))
    {
      /* BOTTOM SIDE */
      if (CellinXYBoundaryWithStencil(cell_id, stencil_index) && stencil_index == 0)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x, y - 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* LOWER RIGHT CORNER */
      else if (CellinXYBoundaryWithStencil(cell_id, stencil_index) && stencil_index == 1)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x + 1, y - 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* LEFT SIDE */
      else if (CellinXYBoundaryWithStencil(cell_id, stencil_index) && stencil_index == 2)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x - 1, y, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* CURRENT */
      else if (CellinXYBoundaryWithStencil(cell_id, stencil_index) && stencil_index == 3)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x, y, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* RIGHT SIDE */
      else if (CellinXYBoundaryWithStencil(cell_id, stencil_index) && stencil_index == 4)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x + 1, y, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* UPPER LEFT CORNER */
      else if (CellinXYBoundaryWithStencil(cell_id, stencil_index) && stencil_index == 5)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x - 1, y + 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* TOP SIDE */
      else if (CellinXYBoundaryWithStencil(cell_id, stencil_index) && stencil_index == 6)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x, y + 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }
    }
    else
    {
      /* BOTTOM SIDE */
      if (stencil_index == 0)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x, y - 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* LOWER RIGHT CORNER */
      else if (stencil_index == 1)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x + 1, y - 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* LEFT SIDE */
      else if (stencil_index == 2)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x - 1, y, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* CURRENT */
      else if (stencil_index == 3)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x, y, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* RIGHT SIDE */
      else if (stencil_index == 4)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x + 1, y, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* UPPER LEFT CORNER */
      else if (stencil_index == 5)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x - 1, y + 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }

      /* TOP SIDE */
      else if (stencil_index == 6)
      {
        lattice_centroid = _lattice->getLatticeCellCenter(x, y + 1, z);
        lattice_centroid_x = lattice_centroid.getX();
        lattice_centroid_y = lattice_centroid.getY();
        dist_x = pow(centroid_x - lattice_centroid_x, 2.0);
        dist_y = pow(centroid_y - lattice_centroid_y, 2.0);
        found = true;
      }
    }

    if (found)
      return pow(dist_x + dist_y, 0.5);
    else
      return std::numeric_limits<double>::max();
  }

  /** @brief Set a pointer to the Geometry.
   * @param goemetry A pointer to a Geometry object.
   */
  void Cmfd::setGeometry(Geometry *geometry)
  {
    _geometry = geometry;
  }

  /** @brief Set a number of k-nearest neighbor cells to use in updating
   *         the FSR flux.设置用于更新 FSR 通量时使用的 “k 近邻” CMFD 单元数量。
   * @param k_nearest The number of nearest neighbor CMFD cells. 选取的最近邻 CMFD 单元数
   */
  void Cmfd::setKNearest(int k_nearest)
  {

    if (k_nearest < 1 || k_nearest > 9)
      log::ferror("Unable to set CMFD k-nearest to %i. k-nearest "
                  "must be between 1 and 9.",
                  k_nearest);
    else
      _k_nearest = k_nearest; // 用于更新MOC通量的CMFD数
  }

  /**
   * @brief Zero the surface currents for each mesh cell and energy group.
   */
  void Cmfd::zeroCurrents()
  {

    _surface_currents->clear();

    if (_balance_sigma_t)
    {
      _starting_currents->clear();
      _net_currents->clear();
    }

    // Clear boundary currents
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
    {
#pragma omp parallel for
      for (int s = 0; s < NUM_FACES; s++)
      {

        // Loop over all CMFD cells on the current surface
        std::map<int, int>::iterator it;
        for (it = _boundary_index_map.at(s).begin(); //_boundary_index_map.at(s).begin(),应该是边界域的中的CMFD单元的全局索引
             it != _boundary_index_map.at(s).end(); ++it)
        {

          int idx = it->second; // 它的第一个是域的左侧第一个，不是本域第一个，它的下一个才是第一个CMFD全局索引

          // Loop over CMFD coarse energy groups
          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            // Loop over cell faces
            for (int f = 0; f < NUM_FACES; f++)
              _boundary_surface_currents[s][idx][f * _num_cmfd_groups + e] = 0.0; // 第三维大小为：comm_data_size + internal

            // Loop over all cell faces and edges
            for (int fe = 0; fe < NUM_FACES + NUM_EDGES; fe++)
            {
              _off_domain_split_currents[s][idx][fe * _num_cmfd_groups + e] = 0.0;
              _received_split_currents[s][idx][fe * _num_cmfd_groups + e] = 0.0;
            }
          }
        }
      }
    }
#endif
  }

  /**
   * @brief Zero the Hex surface currents for each mesh cell and energy group.
   */
  void Cmfd::zeroHexCurrents()
  {

    _surface_currents->clear();

    if (_balance_sigma_t)
    {
      _starting_currents->clear();
      _net_currents->clear();
    }

    // Clear boundary currents
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
    {
#pragma omp parallel for
      for (int s = 0; s < HEX_NUM_FACES; s++)
      {
        // Loop over all CMFD cells on the current surface
        std::map<int, int>::iterator it;
        for (it = _boundary_index_map.at(s).begin();
             it != _boundary_index_map.at(s).end(); ++it)
        {

          int idx = it->second;

          // Loop over CMFD coarse energy groups
          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            // Loop over cell faces
            for (int f = 0; f < HEX_NUM_FACES; f++)
              _boundary_surface_currents[s][idx][f * _num_cmfd_groups + e] = 0.0;

            // Loop over all cell faces and edges
            for (int fe = 0; fe < HEX_NUM_FACES + HEX_NUM_EDGES; fe++)
            {
              _off_domain_split_currents[s][idx][fe * _num_cmfd_groups + e] = 0.0;
              _received_split_currents[s][idx][fe * _num_cmfd_groups + e] = 0.0;
            }
          }
        }
      }
    }
#endif
  }

  /**
   * @brief Initialize the Matrix and Vector objects, k-nearest stencils, the
   *        CMFD cell currents and MOC materials.
   * 初始化CMFD相关的矩阵和矢量对象、k-nearest模板、CMFD中子流和MOC材料
   */
  void Cmfd::initialize()
  {

    /* Delete old Matrix and Vector objects if they exist */
    if (_A != NULL) // 对应CMFD文档中的B.29公式
      delete _A;
    if (_M != NULL)
      delete _M;
    if (_old_source != NULL) // The old source vector
      delete _old_source;
    if (_new_source != NULL)
      delete _new_source;
    if (_old_flux != NULL) // B.32，CMFD求解前的通量
      delete _old_flux;
    if (_new_flux != NULL)
      delete _new_flux;
    if (_old_dif_surf_corr != NULL) // 上一次迭代的表面修正扩散系数
      delete _old_dif_surf_corr;
    if (_volumes != NULL)
      delete _volumes;
    if (_cell_locks != NULL) // CMFD互斥锁
      delete[] _cell_locks;

    /* Calculate the number of elements */
    int num_cells = _local_num_x * _local_num_y * _local_num_z;
    int ncg = _num_cmfd_groups;

    try
    {

      /* Allocate temporary tally vectors for surface currents by thread 按线程为表面中子流分配临时记录向量 */
      int num_threads = omp_get_max_threads();
      _temporary_currents = new CMFD_PRECISION *[num_threads];
      for (int t = 0; t < num_threads; t++)
        _temporary_currents[t] = new CMFD_PRECISION[ncg]; //_temporary_currents:用于记录每个线程计算出的临时表面中子流

      /* Allocate array of OpenMP locks for each CMFD cell 分配锁空间*/
      _cell_locks = new omp_lock_t[num_cells];

      /* Loop over all cells to initialize OpenMP locks 初始化锁*/
#pragma omp parallel for schedule(guided)
      for (int r = 0; r < num_cells; r++)
        omp_init_lock(&_cell_locks[r]);
      omp_init_lock(&_edge_corner_lock);

      // 为各个变量分配空间
      _M = new Matrix(_cell_locks, _local_num_x, _local_num_y, _local_num_z,
                      ncg);
      _A = new Matrix(_cell_locks, _local_num_x, _local_num_y, _local_num_z,
                      ncg);
      _old_source = new Vector(_cell_locks, _local_num_x, _local_num_y,
                               _local_num_z, ncg);
      _new_source = new Vector(_cell_locks, _local_num_x, _local_num_y,
                               _local_num_z, ncg);
      _old_flux = new Vector(_cell_locks, _local_num_x, _local_num_y,
                             _local_num_z, ncg);
      _new_flux = new Vector(_cell_locks, _local_num_x, _local_num_y,
                             _local_num_z, ncg);
      _old_dif_surf_corr = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                      _local_num_z, NUM_FACES * ncg);
      _old_dif_surf_corr->setAll(0.0);
      _volumes = new Vector(_cell_locks, _local_num_x, _local_num_y, _local_num_z, 1);

      /* Initialize k-nearest stencils, currents, flux, materials and tallies */
      generateKNearestStencils(); // 设置周围的K个CMFD网格对FSR通量更新的贡献（百分比），数据保存在了_k_nearest_stencils这个模具中
      initializeCurrents();       // 在第一次MOC迭代之前初始化CMFD表面中子流向量,有表面电流，净电流，起始边界电流
      initializeMaterials();      // 初始化CMFD材料内存
      allocateTallies();          // 为每个CMFD单元分配扩散、反应和体积统计的内存all_tallies

      /* Initialize domain communicator */
      if (_domain_communicator != NULL) // 这里是区域分解的初始化代码
      {

        /* Size of domain in each direction */
        _local_num_x = _num_x / _domain_communicator->_num_domains_x;
        _local_num_y = _num_y / _domain_communicator->_num_domains_y;
        _local_num_z = _num_z / _domain_communicator->_num_domains_z;
        _domain_communicator->stop = false;
        int offset = _domain_communicator->_domain_idx_x * _local_num_x + //_domain_idx_x:当前域在x方向上的索引
                     _domain_communicator->_domain_idx_y * _local_num_y +
                     _domain_communicator->_domain_idx_z * _local_num_z;
        _domain_communicator->_offset = offset; // 计算当前进程（或当前子域）在全局CMFD网格中的起始位置
        _domain_communicator->_local_num_x = _local_num_x;
        _domain_communicator->_local_num_y = _local_num_y;
        _domain_communicator->_local_num_z = _local_num_z;
        _domain_communicator->num_groups = ncg;

        int dir_sizes[3] = {num_cells / _local_num_x, num_cells / _local_num_y,
                            num_cells / _local_num_z}; // 当前子域在x，y，z方向上每层包含的CMFD单元格数量

        /* Allocate arrays to contain information about the domain's neighbors 存储每个表面单元格连接邻居的信息 */
        _domain_communicator->num_connections = new int *[2]; // 连接邻居的数量
        _domain_communicator->indexes = new int **[2];        // 邻居的索引
        _domain_communicator->domains = new int **[2];        // 邻居的表面编号

        /* Arrays to contain data to communicate to/receive from other domains 包含要与其他域通信/从其他域接收的数据的阵列 */
        _domain_communicator->fluxes = new CMFD_PRECISION **[2];          // 每个表面单元的连接邻居的通量
        _domain_communicator->coupling_coeffs = new CMFD_PRECISION **[2]; // 每个表面单元的连接邻居与其自身之间的耦合系数
        _domain_communicator->buffer = new CMFD_PRECISION *[NUM_FACES];   // 用于向连接邻居发送通量/从连接邻居接收通量的缓冲区
        for (int rb = 0; rb < 2; rb++)
        {
          _domain_communicator->num_connections[rb] = new int[num_cells * ncg];
          _domain_communicator->indexes[rb] = new int *[num_cells * ncg];
          _domain_communicator->domains[rb] = new int *[num_cells * ncg];
          _domain_communicator->fluxes[rb] = new CMFD_PRECISION *[NUM_FACES];
          _domain_communicator->coupling_coeffs[rb] =
              new CMFD_PRECISION *[num_cells * ncg];

          for (int coord = 0; coord < 3; coord++)
          {
            for (int d = 0; d < 2; d++)
            {
              int surf = coord + 3 * d; // 0,3,1,4,2,5
              _domain_communicator->fluxes[rb][surf] =
                  new CMFD_PRECISION[dir_sizes[coord] * ncg];
              _domain_communicator->buffer[surf] =
                  new CMFD_PRECISION[2 * dir_sizes[coord] * ncg];
            }
          }
          for (int nsc = 0; nsc < num_cells * ncg; nsc++)
          {
            _domain_communicator->num_connections[rb][nsc] = 0;
            _domain_communicator->indexes[rb][nsc] = new int[NUM_FACES];
            _domain_communicator->domains[rb][nsc] = new int[NUM_FACES];
            _domain_communicator->coupling_coeffs[rb][nsc] =
                new CMFD_PRECISION[NUM_FACES];
          }
          _domain_communicator_allocated = true;
        }

        int storage_per_cell = ((2 + NUM_FACES) * ncg + 1); // 每个单元格中要存储的基本信息的数量,这个2和1代表什么？
        int num_per_side[3] = {_local_num_y * _local_num_z,
                               _local_num_x * _local_num_z,
                               _local_num_x * _local_num_y};

        /* Count total number of cells at all faces of the domain 计算域所有面上的单元格总数 */
        int num_boundary_cells = 0;
        for (int s = 0; s < NUM_FACES; s++)
          num_boundary_cells += num_per_side[s % 3];

        int internal = ncg * num_boundary_cells;
        int comm_data_size = storage_per_cell * num_boundary_cells; // 每个单元格中要存储的基本信息的数量*所有面上的单元格总数

        _inter_domain_data = new CMFD_PRECISION[comm_data_size + internal]; // 包含接收数据的缓冲区
        _send_domain_data = new CMFD_PRECISION[comm_data_size];             // 包含发送数据的缓冲区

        /* Allocate memory for communication of off-domain quantities 为域外数量的通信分配内存*/
        _domain_data_by_surface = new CMFD_PRECISION *[NUM_FACES]; // 对于每个面（数组的第一个维度），将包含接收到的数据
        _send_data_by_surface = new CMFD_PRECISION *[NUM_FACES];   // 对于每个面（数组的第一个维度），将包含要发送的数据
        _boundary_volumes = new CMFD_PRECISION **[NUM_FACES];      // 域边界通信的体积缓冲区 以下4个变量同
        _boundary_reaction = new CMFD_PRECISION **[NUM_FACES];
        _boundary_diffusion = new CMFD_PRECISION **[NUM_FACES];
        _boundary_surface_currents = new CMFD_PRECISION **[NUM_FACES];

        _old_boundary_flux = new CMFD_PRECISION **[NUM_FACES];

        int start = 0;
        int ext = 0;
        for (int s = 0; s < NUM_FACES; s++)
        {

          _domain_data_by_surface[s] = &_inter_domain_data[start];
          _send_data_by_surface[s] = &_send_domain_data[start];
          _boundary_volumes[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _boundary_reaction[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _boundary_diffusion[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _boundary_surface_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];

          _old_boundary_flux[s] = new CMFD_PRECISION *[num_per_side[s % 3]];

          for (int idx = 0; idx < num_per_side[s % 3]; idx++)
          {

            _boundary_volumes[s][idx] = &_inter_domain_data[start];
            _boundary_reaction[s][idx] = &_inter_domain_data[start + 1];
            _boundary_diffusion[s][idx] = &_inter_domain_data[start + ncg + 1];
            _boundary_surface_currents[s][idx] =
                &_inter_domain_data[start + 2 * ncg + 1];

            _old_boundary_flux[s][idx] = &_inter_domain_data[comm_data_size + ext];

            ext += ncg;
            start += storage_per_cell; // ((2 + NUM_FACES) * ncg + 1)
          }
        }

        /* Allocate memory for split current communication 为分流通信分配内存*/
        int ns = NUM_FACES + NUM_EDGES; // 面+边 = 18
        int vec_size = ns * ncg * sizeof(CMFD_PRECISION);
        int split_current_size = ncg * ns * num_boundary_cells; // num_boundary_cells：域中所有面（6个面）的CMFD数量
        _send_split_current_data = new CMFD_PRECISION[split_current_size];
        _receive_split_current_data = new CMFD_PRECISION[split_current_size];

        _send_split_currents_array = new CMFD_PRECISION *[NUM_FACES];
        _receive_split_currents_array = new CMFD_PRECISION *[NUM_FACES];
        _off_domain_split_currents = new CMFD_PRECISION **[NUM_FACES];
        _received_split_currents = new CMFD_PRECISION **[NUM_FACES];

        start = 0;
        for (int s = 0; s < NUM_FACES; s++)
        {

          _send_split_currents_array[s] =
              &_send_split_current_data[start];
          _receive_split_currents_array[s] =
              &_receive_split_current_data[start];
          _off_domain_split_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _received_split_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];

          for (int idx = 0; idx < num_per_side[s % 3]; idx++)
          {

            _off_domain_split_currents[s][idx] =
                &_send_split_current_data[start];
            memset(_off_domain_split_currents[s][idx], 0.0, vec_size);
            _received_split_currents[s][idx] =
                &_receive_split_current_data[start];
            memset(_received_split_currents[s][idx], 0.0, vec_size);

            start += ns * ncg;
          }
        }

        /* Allocate memory for communication of on-domain quantities 为域内数量的通信分配内存 */
        _send_volumes = new CMFD_PRECISION **[NUM_FACES];
        _send_reaction = new CMFD_PRECISION **[NUM_FACES];
        _send_diffusion = new CMFD_PRECISION **[NUM_FACES];
        _send_currents = new CMFD_PRECISION **[NUM_FACES];

        start = 0;
        for (int s = 0; s < NUM_FACES; s++)
        {
          _send_volumes[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _send_reaction[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _send_diffusion[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _send_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          for (int idx = 0; idx < num_per_side[s % 3]; idx++)
          {
            _send_volumes[s][idx] = &_send_domain_data[start];
            _send_reaction[s][idx] = &_send_domain_data[start + 1]; //_send_domain_data = new CMFD_PRECISION[comm_data_size];  //包含发送数据的缓冲区
            _send_diffusion[s][idx] = &_send_domain_data[start + ncg + 1];
            _send_currents[s][idx] = &_send_domain_data[start + 2 * ncg + 1];
            start += storage_per_cell;
          }
        }

        /* Calculate the starting and ending indexes of on-domain CMFD cells 计算域上CMFD单元的起始索引和结束索引 */
        int x_start = _domain_communicator->_domain_idx_x * _local_num_x;
        int x_end = x_start + _local_num_x;
        int y_start = _domain_communicator->_domain_idx_y * _local_num_y;
        int y_end = y_start + _local_num_y;
        int z_start = _domain_communicator->_domain_idx_z * _local_num_z;
        int z_end = z_start + _local_num_z;

        _boundary_index_map.resize(NUM_FACES); // 界面的索引映射信息

        /* Map connecting cells on x-surfaces 映射x表面上的连接单元*/
        int global_ind;
        for (int y = 0; y < _local_num_y; y++)
        {
          for (int z = 0; z < _local_num_z; z++)
          {
            if (x_start - 1 >= 0)
            { // 用于检查当前域左侧是否有相邻的域单元。如果 x_start - 1 >= 0，说明在当前域的左侧还有域单元存在，左侧边界不是整个几何的边界。
              global_ind = ((z_start + z) * _num_y + y + y_start) *
                               _num_x +
                           x_start - 1; // 左侧域CMFD单元的全局索引
              /*
              计算得到全局索引 global_ind，并将其映射到 _boundary_index_map 中，表示 SURFACE_X_MIN 处的连接单元。
              这个映射表明当前域左侧边界的单元和相邻域的哪个单元相连接。
               */
              _boundary_index_map.at(SURFACE_X_MIN)[global_ind] = z * _local_num_y + y;
            }
            if (x_end < _num_x)
            { // 用于检查当前域右侧是否到达几何边界的右侧
              global_ind = ((z_start + z) * _num_y + y + y_start) *
                               _num_x +
                           x_end;

              _boundary_index_map.at(SURFACE_X_MAX)[global_ind] = z * _local_num_y + y;
            }
          }
        }

        /* Map connecting cells on y-surfaces 同上*/
        for (int x = 0; x < _local_num_x; x++)
        {
          for (int z = 0; z < _local_num_z; z++)
          {
            if (y_start - 1 >= 0)
            {
              global_ind = ((z_start + z) * _num_y + y_start - 1) *
                               _num_x +
                           x + x_start;
              _boundary_index_map.at(SURFACE_Y_MIN)[global_ind] = z * _local_num_x + x;
            }
            if (y_end < _num_y)
            {
              global_ind = ((z_start + z) * _num_y + y_end) * _num_x + x + x_start;
              _boundary_index_map.at(SURFACE_Y_MAX)[global_ind] = z * _local_num_x + x;
            }
          }
        }

        /* Map connecting cells on z-surfaces 同上*/
        for (int x = 0; x < _local_num_x; x++)
        {
          for (int y = 0; y < _local_num_y; y++)
          {
            if (z_start - 1 >= 0)
            {
              global_ind = ((z_start - 1) * _num_y + y + y_start) *
                               _num_x +
                           x + x_start;
              _boundary_index_map.at(SURFACE_Z_MIN)[global_ind] = y * _local_num_x + x;
            }
            if (z_end < _num_z)
            {
              global_ind = ((z_end)*_num_y + y + y_start) *
                               _num_x +
                           x + x_start;
              _boundary_index_map.at(SURFACE_Z_MAX)[global_ind] = y * _local_num_x + x;
            }
          }
        }
      }
    }
    catch (std::exception &e)
    {
      log::ferror("Could not allocate memory for the CMFD mesh objects. "
                  "Backtrace:%s",
                  e.what());
    }
  }

  /**
   * @brief Initialize the Hex Lattice Matrix and Vector objects, k-nearest stencils, the
   *        CMFD cell currents and MOC materials.
   */
  void Cmfd::initializeHex()
  {

    /* Delete old Matrix and Vector objects if they exist */
    if (_A != NULL)
      delete _A;
    if (_M != NULL)
      delete _M;
    if (_old_source != NULL)
      delete _old_source;
    if (_new_source != NULL)
      delete _new_source;
    if (_old_flux != NULL)
      delete _old_flux;
    if (_new_flux != NULL)
      delete _new_flux;
    if (_old_dif_surf_corr != NULL)
      delete _old_dif_surf_corr;
    if (_volumes != NULL)
      delete _volumes;
    if (_cell_locks != NULL)
      delete[] _cell_locks;

    /* Calculate the number of elements */
    int num_cells = _local_num_x * _local_num_y * _local_num_z - _empty_cells_num;
    int ncg = _num_cmfd_groups;

    try
    {

      /* Allocate temporary tally vectors for surface currents by thread */
      int num_threads = omp_get_max_threads();
      _temporary_currents = new CMFD_PRECISION *[num_threads];
      for (int t = 0; t < num_threads; t++)
        _temporary_currents[t] = new CMFD_PRECISION[ncg];

      /* Allocate array of OpenMP locks for each CMFD cell */
      _cell_locks = new omp_lock_t[num_cells];

      /* Loop over all cells to initialize OpenMP locks */
#pragma omp parallel for schedule(guided)
      for (int r = 0; r < num_cells; r++)
        omp_init_lock(&_cell_locks[r]);
      omp_init_lock(&_edge_corner_lock);

      /* Allocate memory for matrix and vector objects */
      _M = new Matrix(_cell_locks, _local_num_x, _local_num_y, _local_num_z,
                      ncg, _empty_cells_num);
      _A = new Matrix(_cell_locks, _local_num_x, _local_num_y, _local_num_z,
                      ncg, _empty_cells_num);
      _old_source = new Vector(_cell_locks, _local_num_x, _local_num_y,
                               _local_num_z, ncg, _empty_cells_num);
      _new_source = new Vector(_cell_locks, _local_num_x, _local_num_y,
                               _local_num_z, ncg, _empty_cells_num);
      _old_flux = new Vector(_cell_locks, _local_num_x, _local_num_y,
                             _local_num_z, ncg, _empty_cells_num);
      _new_flux = new Vector(_cell_locks, _local_num_x, _local_num_y,
                             _local_num_z, ncg, _empty_cells_num);
      _old_dif_surf_corr = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                      _local_num_z, HEX_NUM_SURFACES * ncg, _empty_cells_num);
      _old_dif_surf_corr->setAll(0.0);
      _volumes = new Vector(_cell_locks, _local_num_x, _local_num_y, _local_num_z, 1, _empty_cells_num);

      /* Initialize k-nearest stencils, currents, flux, materials and tallies */
      generateKNearestStencilsHex();
      initializeCurrentsHex();
      initializeHexMaterials();
      allocateHexTallies();

      /* need fixed */
      /* Initialize domain communicator */
      if (_domain_communicator != NULL)
      {
        /* Size of domain in each direction */
        _local_num_x = _num_x / _domain_communicator->_num_domains_x;
        _local_num_y = _num_y / _domain_communicator->_num_domains_y;
        _local_num_z = _num_z / _domain_communicator->_num_domains_z;
        _domain_communicator->stop = false;
        int offset = _domain_communicator->_domain_idx_x * _local_num_x +
                     _domain_communicator->_domain_idx_y * _local_num_y +
                     _domain_communicator->_domain_idx_z * _local_num_z;
        _domain_communicator->_offset = offset;
        _domain_communicator->_local_num_x = _local_num_x;
        _domain_communicator->_local_num_y = _local_num_y;
        _domain_communicator->_local_num_z = _local_num_z;
        _domain_communicator->num_groups = ncg;

        int dir_sizes[3] = {num_cells / _local_num_x, num_cells / _local_num_y,
                            num_cells / _local_num_z};

        /* Allocate arrays to contain information about the domain's neighbors */
        _domain_communicator->num_connections = new int *[2];
        _domain_communicator->indexes = new int **[2];
        _domain_communicator->domains = new int **[2];

        /* Arrays to contain data to communicate to/receive from other domains */
        _domain_communicator->fluxes = new CMFD_PRECISION **[2];
        _domain_communicator->coupling_coeffs = new CMFD_PRECISION **[2];
        _domain_communicator->buffer = new CMFD_PRECISION *[HEX_NUM_FACES];
        for (int rb = 0; rb < 2; rb++)
        {
          _domain_communicator->num_connections[rb] = new int[num_cells * ncg];
          _domain_communicator->indexes[rb] = new int *[num_cells * ncg];
          _domain_communicator->domains[rb] = new int *[num_cells * ncg];
          _domain_communicator->fluxes[rb] = new CMFD_PRECISION *[HEX_NUM_FACES];
          _domain_communicator->coupling_coeffs[rb] =
              new CMFD_PRECISION *[num_cells * ncg];

          for (int coord = 0; coord < 3; coord++)
          {
            for (int d = 0; d < 2; d++)
            {
              int surf = coord + 3 * d;
              _domain_communicator->fluxes[rb][surf] =
                  new CMFD_PRECISION[dir_sizes[coord] * ncg];
              _domain_communicator->buffer[surf] =
                  new CMFD_PRECISION[2 * dir_sizes[coord] * ncg];
            }
          }
          for (int nsc = 0; nsc < num_cells * ncg; nsc++)
          {
            _domain_communicator->num_connections[rb][nsc] = 0;
            _domain_communicator->indexes[rb][nsc] = new int[HEX_NUM_FACES];
            _domain_communicator->domains[rb][nsc] = new int[HEX_NUM_FACES];
            _domain_communicator->coupling_coeffs[rb][nsc] =
                new CMFD_PRECISION[HEX_NUM_FACES];
          }
          _domain_communicator_allocated = true;
        }

        int storage_per_cell = ((2 + HEX_NUM_FACES) * ncg + 1);
        int num_per_side[3] = {_local_num_y * _local_num_z,
                               _local_num_x * _local_num_z,
                               _local_num_x * _local_num_y};

        /* Count total number of cells at all faces of the domain */
        int num_boundary_cells = 0;
        for (int s = 0; s < HEX_NUM_FACES; s++)
          num_boundary_cells += num_per_side[s % 3];

        int internal = ncg * num_boundary_cells;
        int comm_data_size = storage_per_cell * num_boundary_cells;

        _inter_domain_data = new CMFD_PRECISION[comm_data_size + internal];
        _send_domain_data = new CMFD_PRECISION[comm_data_size];

        /* Allocate memory for communication of off-domain quantities */
        _domain_data_by_surface = new CMFD_PRECISION *[HEX_NUM_FACES];
        _send_data_by_surface = new CMFD_PRECISION *[HEX_NUM_FACES];
        _boundary_volumes = new CMFD_PRECISION **[HEX_NUM_FACES];
        _boundary_reaction = new CMFD_PRECISION **[HEX_NUM_FACES];
        _boundary_diffusion = new CMFD_PRECISION **[HEX_NUM_FACES];
        _boundary_surface_currents = new CMFD_PRECISION **[HEX_NUM_FACES];

        _old_boundary_flux = new CMFD_PRECISION **[HEX_NUM_FACES];

        int start = 0;
        int ext = 0;
        for (int s = 0; s < HEX_NUM_FACES; s++)
        {

          _domain_data_by_surface[s] = &_inter_domain_data[start];
          _send_data_by_surface[s] = &_send_domain_data[start];
          _boundary_volumes[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _boundary_reaction[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _boundary_diffusion[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _boundary_surface_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];

          _old_boundary_flux[s] = new CMFD_PRECISION *[num_per_side[s % 3]];

          for (int idx = 0; idx < num_per_side[s % 3]; idx++)
          {

            _boundary_volumes[s][idx] = &_inter_domain_data[start];
            _boundary_reaction[s][idx] = &_inter_domain_data[start + 1];
            _boundary_diffusion[s][idx] = &_inter_domain_data[start + ncg + 1];
            _boundary_surface_currents[s][idx] =
                &_inter_domain_data[start + 2 * ncg + 1];

            _old_boundary_flux[s][idx] = &_inter_domain_data[comm_data_size + ext];

            ext += ncg;
            start += storage_per_cell;
          }
        }

        /* Allocate memory for split current communication */
        int ns = HEX_NUM_FACES + HEX_NUM_EDGES;
        int vec_size = ns * ncg * sizeof(CMFD_PRECISION);
        int split_current_size = ncg * ns * num_boundary_cells;
        _send_split_current_data = new CMFD_PRECISION[split_current_size];
        _receive_split_current_data = new CMFD_PRECISION[split_current_size];

        _send_split_currents_array = new CMFD_PRECISION *[HEX_NUM_FACES];
        _receive_split_currents_array = new CMFD_PRECISION *[HEX_NUM_FACES];
        _off_domain_split_currents = new CMFD_PRECISION **[HEX_NUM_FACES];
        _received_split_currents = new CMFD_PRECISION **[HEX_NUM_FACES];

        start = 0;
        for (int s = 0; s < HEX_NUM_FACES; s++)
        {

          _send_split_currents_array[s] =
              &_send_split_current_data[start];
          _receive_split_currents_array[s] =
              &_receive_split_current_data[start];
          _off_domain_split_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _received_split_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];

          for (int idx = 0; idx < num_per_side[s % 3]; idx++)
          {

            _off_domain_split_currents[s][idx] =
                &_send_split_current_data[start];
            memset(_off_domain_split_currents[s][idx], 0.0, vec_size);
            _received_split_currents[s][idx] =
                &_receive_split_current_data[start];
            memset(_received_split_currents[s][idx], 0.0, vec_size);

            start += ns * ncg;
          }
        }

        /* Allocate memory for communication of on-domain quantities */
        _send_volumes = new CMFD_PRECISION **[HEX_NUM_FACES];
        _send_reaction = new CMFD_PRECISION **[HEX_NUM_FACES];
        _send_diffusion = new CMFD_PRECISION **[HEX_NUM_FACES];
        _send_currents = new CMFD_PRECISION **[HEX_NUM_FACES];

        start = 0;
        for (int s = 0; s < HEX_NUM_FACES; s++)
        {
          _send_volumes[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _send_reaction[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _send_diffusion[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          _send_currents[s] = new CMFD_PRECISION *[num_per_side[s % 3]];
          for (int idx = 0; idx < num_per_side[s % 3]; idx++)
          {
            _send_volumes[s][idx] = &_send_domain_data[start];
            _send_reaction[s][idx] = &_send_domain_data[start + 1];
            _send_diffusion[s][idx] = &_send_domain_data[start + ncg + 1];
            _send_currents[s][idx] = &_send_domain_data[start + 2 * ncg + 1];
            start += storage_per_cell;
          }
        }

        /* Calculate the starting and ending indexes of on-domain CMFD cells */
        int x_start = _domain_communicator->_domain_idx_x * _local_num_x;
        int x_end = x_start + _local_num_x;
        int y_start = _domain_communicator->_domain_idx_y * _local_num_y;
        int y_end = y_start + _local_num_y;
        int z_start = _domain_communicator->_domain_idx_z * _local_num_z;
        int z_end = z_start + _local_num_z;

        _boundary_index_map.resize(HEX_NUM_FACES);

        /* Map connecting cells on x-surfaces */
        int global_ind;
        for (int y = 0; y < _local_num_y; y++)
        {
          for (int z = 0; z < _local_num_z; z++)
          {
            if (x_start - 1 >= 0)
            {
              global_ind = ((z_start + z) * _num_y + y + y_start) *
                               _num_x +
                           x_start - 1;
              _boundary_index_map.at(HEX_SURFACE_BETA_MIN)[global_ind] = z * _local_num_y + y;
            }
            if (x_end < _num_x)
            {
              global_ind = ((z_start + z) * _num_y + y + y_start) *
                               _num_x +
                           x_end;

              _boundary_index_map.at(HEX_SURFACE_BETA_MAX)[global_ind] = z * _local_num_y + y;
            }
          }
        }

        /* Map connecting cells on y-surfaces */
        for (int x = 0; x < _local_num_x; x++)
        {
          for (int z = 0; z < _local_num_z; z++)
          {
            if (y_start - 1 >= 0)
            {
              global_ind = ((z_start + z) * _num_y + y_start - 1) *
                               _num_x +
                           x + x_start;
              _boundary_index_map.at(HEX_SURFACE_DELTA_MIN)[global_ind] = z * _local_num_x + x;
            }
            if (y_end < _num_y)
            {
              global_ind = ((z_start + z) * _num_y + y_end) * _num_x + x + x_start;
              _boundary_index_map.at(HEX_SURFACE_DELTA_MAX)[global_ind] = z * _local_num_x + x;
            }
          }
        }

        /* Map connecting cells on z-surfaces */
        for (int x = 0; x < _local_num_x; x++)
        {
          for (int y = 0; y < _local_num_y; y++)
          {
            if (z_start - 1 >= 0)
            {
              global_ind = ((z_start - 1) * _num_y + y + y_start) *
                               _num_x +
                           x + x_start;
              _boundary_index_map.at(HEX_SURFACE_Z_MIN)[global_ind] = y * _local_num_x + x;
            }
            if (z_end < _num_z)
            {
              global_ind = ((z_end)*_num_y + y + y_start) *
                               _num_x +
                           x + x_start;
              _boundary_index_map.at(HEX_SURFACE_Z_MAX)[global_ind] = y * _local_num_x + x;
            }
          }
        }
      }
    }
    catch (std::exception &e)
    {
      log::ferror("Could not allocate memory for the CMFD mesh objects. "
                  "Backtrace:%s",
                  e.what());
    }
  }

  /**
   * @brief Initialize the CMFD lattice.
   */
  void Cmfd::initializeLattice(Point *offset)
  {

    if (_hexlattice_enable)
    { // Create HexLattice

      _cell_width_z = _width_z / _num_z;
      _cell_widths_z.resize(_num_z, _cell_width_z);
      _accumulate_z.resize(_num_z + 1, 0.0);

      for (int i = 0; i < _num_z; i++)
        _accumulate_z[i + 1] = _accumulate_z[i] + _cell_widths_z[i];

      /* Delete old lattice if it exists */
      if (_lattice != NULL)
        delete _lattice;

      /* Initialize the lattice */
      _lattice = new HexLattice();
      _lattice->setNumR(_num_r);
      _lattice->setNumZ(_num_z);

      _lattice->setWidths(_width_r, _cell_widths_z);
      _lattice->setOrientation(_orientation);

      _lattice->setOffset(offset->getX(), offset->getY(), offset->getZ());
      _lattice->computeSizes();

      setNumX(_lattice->getNumX());
      setNumY(_lattice->getNumY());
    }
    else
    { // Create RecLattice - 创建四边形lattice
      if (_non_uniform)
      { // _non_uniform：布尔标志，表示网格是否非均匀；非均匀：每个网格宽度不同，例如 [2, 1, 3]
        // 从已设置的宽度数组获取网格数量
        setNumX(_cell_widths_x.size());
        setNumY(_cell_widths_y.size());
        setNumZ(_cell_widths_z.size());
      }
      else
      {
        // 均匀：每个方向的单元宽度 = 总宽度 / 单元数
        _cell_width_x = _width_x / _num_x;
        _cell_width_y = _width_y / _num_y;
        _cell_width_z = _width_z / _num_z;

        _cell_widths_x.resize(_num_x, _cell_width_x); // 如果是均匀的CMFD，则将_cell_widths_x数组的所有宽度设为一样的大小
        _cell_widths_y.resize(_num_y, _cell_width_y);
        _cell_widths_z.resize(_num_z, _cell_width_z);
      }

      // 计算三个方向的累计宽度（前缀和），用于定位单元边界与做宽度一致性校验
      _accumulate_x.resize(_num_x + 1, 0.0);
      // 将 _accumulate_x 调整为长度 _num_x + 1 的向量，并把所有元素初始化为 0.0 。num_x 个单元对应 num_x + 1 个边界（含起点和终点）。例如3个格子需要4个边界位置。
      _accumulate_y.resize(_num_y + 1, 0.0);
      _accumulate_z.resize(_num_z + 1, 0.0);

      //_accumulate_x = [0.0, 2.0, 4.0, 6.0...]
      for (int i = 0; i < _num_x; i++)
        _accumulate_x[i + 1] = _accumulate_x[i] + _cell_widths_x[i]; // X 方向累计宽度

      for (int i = 0; i < _num_y; i++)
        _accumulate_y[i + 1] = _accumulate_y[i] + _cell_widths_y[i]; // Y 方向累计宽度

      for (int i = 0; i < _num_z; i++)
        _accumulate_z[i + 1] = _accumulate_z[i] + _cell_widths_z[i]; // Z 方向累计宽度

      // 一致性校验：检查累积宽度是否等于总宽度（允许微小误差）
      if (fabs(_width_x - _accumulate_x[_num_x]) > FLT_EPSILON ||
          fabs(_width_y - _accumulate_y[_num_y]) > FLT_EPSILON ||
          fabs(_width_z - _accumulate_z[_num_z]) > FLT_EPSILON)
      {
        log::ferror("The sum of non-uniform mesh widths are not consistent "
                    "with geometry dimensions. width_x = %20.17E, width_y = %20.17E"
                    ", width_z = %20.17E, sum_x = %20.17E, sum_y = %20.17E, sum_z ="
                    " %20.17E, diff_x = %20.17E, diff_y = %20.17E, diff_z = %20.17E"
                    ", FLT_EPSILON = %20.17E",
                    _width_x, _width_y, _width_z,
                    _accumulate_x[_num_x], _accumulate_y[_num_y],
                    _accumulate_z[_num_z], fabs(_width_x - _accumulate_x[_num_x]),
                    fabs(_width_y - _accumulate_y[_num_y]),
                    fabs(_width_z - _accumulate_z[_num_z]), FLT_EPSILON);
      }

      /* Delete old lattice if it exists - 若存在旧 lattice 则删除 */
      if (_lattice != NULL)
        delete _lattice;

      /* Initialize the lattice */
      _lattice = new RecLattice();
      _lattice->setNumX(_num_x); // 设置 x 方向单元数
      _lattice->setNumY(_num_y); // 设置 y 方向单元数
      _lattice->setNumZ(_num_z); // 设置 z 方向单元数

      if (_non_uniform)
        _lattice->setWidths(_cell_widths_x, _cell_widths_y, _cell_widths_z); // 非均匀：传入每个方向的单元宽度数组
      else
        _lattice->setWidth(_cell_width_x, _cell_width_y, _cell_width_z); // 均匀：传入每个方向的单元宽度即可

      _lattice->setOffset(offset->getX(), offset->getY(), offset->getZ()); // 设置lattice偏移量，offset是整个几何的中心点
      _lattice->computeSizes();                                            // 计算lattice三个方向的累计宽度（前缀和），用于定位单元边界，与前面的_accumulate_x的逻辑一样
    }
  }

  /**
   * @brief Set width of non-uniform meshes in x y z directions.
   * @details An example of how this may be called from Python illustrated below:
   *
   * @code
   *          cmfd.setWidths([[1,2,3], [4,5,6,7], [3.3,2.4]])
   * @endcode
   *
   * @param widths A vector of 3 vectors for the x y z sizes of non-uniform meshes
   */
  void Cmfd::setWidths(std::vector<std::vector<double>> widths)
  {

    if (widths.size() == 3)
      _cell_widths_z = widths[2];
    else if (widths.size() == 2)
      _cell_widths_z.push_back(1.0);
    else
      log::ferror("CMFD lattice widths must have dimension 2 or 3.");

    _non_uniform = true;
    _cell_widths_x = widths[0];
    _cell_widths_y = widths[1];
  }

  /**
   * @brief Initializes a backup CMFD solver.
   * @details This backup solver is not necessary to run simulations, but may be
   *          used if the regular solver fails and the user wants to try another
   *          group structure without restarting the simulation.
   * 用户希望在不重新启动模拟的情况下尝试另一种少能群结构
   */
  void Cmfd::initializeBackupCmfdSolver()
  {

    /* Initialize new CMFD object 配置备份 CMFD 的参数*/
    _backup_cmfd = new Cmfd();
    _backup_cmfd->useAxialInterpolation(_use_axial_interpolation);
    _backup_cmfd->setLatticeStructure(_num_x, _num_y, _num_z);
    _backup_cmfd->setKNearest(_k_nearest);
    _backup_cmfd->setSORRelaxationFactor(_SOR_factor);
    _backup_cmfd->setCMFDRelaxationFactor(_relaxation_factor);
    _backup_cmfd->useFluxLimiting(_flux_limiting);
    _backup_cmfd->setHexLatticeEnable(_hexlattice_enable);

    /* Set one-group group structure */
    if (_backup_group_structure.size() == 0)
    { // 如果没有指定备份群结构，则自动将所有群设置为一群

      std::vector<int> all_groups;
      for (int e = 0; e < _num_moc_groups; e++)
        all_groups.push_back(e + 1);
      _backup_group_structure.push_back(all_groups);

      _cmfd_group_to_backup_group = new int[_num_cmfd_groups];
      for (int e = 0; e < _num_cmfd_groups; e++)
        _cmfd_group_to_backup_group[e] = 0;
    }
    _backup_cmfd->setGroupStructure(_backup_group_structure);

    /* Set CMFD mesh boundary conditions */
    if (_hexlattice_enable)
    {
      for (int i = 0; i < 8; i++)
        _backup_cmfd->setBoundary(i, _boundaries[i]);
    }
    else
    {
      for (int i = 0; i < 6; i++)
        _backup_cmfd->setBoundary(i, _boundaries[i]);
    }

    /* Set CMFD mesh dimensions 设置备份 CMFD 网格在 X、Y、Z 方向的尺寸*/
    _backup_cmfd->setWidthX(_width_x);
    _backup_cmfd->setWidthY(_width_y);
    _backup_cmfd->setWidthZ(_width_z);

    /* Intialize CMFD Maps 初始化 CMFD 映射*/
    _backup_cmfd->initializeCellMap();

    /* Initialize the CMFD lattice 初始化 CMFD 晶格*/
    _backup_cmfd->initializeLattice(_lattice->getOffset());
    _backup_cmfd->setGeometry(_geometry);

#ifdef ENABLE_MPI_
    if (_domain_communicator != NULL)
    {

      _backup_cmfd->setNumDomains(_domain_communicator->_num_domains_x,
                                  _domain_communicator->_num_domains_y,
                                  _domain_communicator->_num_domains_z);
      _backup_cmfd->setDomainIndexes(_domain_communicator->_domain_idx_x,
                                     _domain_communicator->_domain_idx_y,
                                     _domain_communicator->_domain_idx_z);
    }
#endif

    /* Initialize the backup CMFD solver 初始化备份 CMFD 求解器*/
    if (_hexlattice_enable)
      _backup_cmfd->initializeHex();
    else
      _backup_cmfd->initialize();

    /* Intialize the CMFD energy group structure 设置能群结构*/
    _backup_cmfd->setSourceConvergenceThreshold(_source_convergence_threshold);
    _backup_cmfd->setNumMOCGroups(_num_moc_groups);
    _backup_cmfd->initializeGroupMap();

    /* Give CMFD number of FSRs and FSR property arrays CMFD的FSR参数*/
    _backup_cmfd->setSolve3D(_SOLVE_3D);
    _backup_cmfd->setNumFSRs(_num_FSRs);
    _backup_cmfd->setFSRVolumes(_FSR_volumes);
    _backup_cmfd->setFSRMaterials(_FSR_materials);
    _backup_cmfd->setFSRFluxes(_FSR_fluxes);
    _backup_cmfd->setFSRSources(_FSR_sources);
    _backup_cmfd->setQuadrature(_quadrature);
    if (_flux_moments != NULL)
      _backup_cmfd->setFluxMoments(_flux_moments);

    /* Add FSRs to cells CMFD所关联的FSR*/
    _backup_cmfd->setCellFSRs(&_cell_fsrs);

    /* Initialize the backup CMFD solver */
    if (_hexlattice_enable)
      _backup_cmfd->initializeHex();
    else
      _backup_cmfd->initialize();
    _backup_cmfd->setConvergenceData(_convergence_data);
  }

  /**
   * @brief Copies the current from the regular to the backup CMFD solver.
   * @details The currents are condensed to the backup solver's energy structure
   *          when transfered as well.
   */
  void Cmfd::copyCurrentsToBackup()
  {

    /* Clear currents */
    _backup_cmfd->zeroCurrents();

    /* Get the number of backup groups */
    int nbg = _backup_group_structure.size();

    /* Get the local current array */
    Vector *backup_currents = _backup_cmfd->getLocalCurrents();

    if (_hexlattice_enable)
    {
      /* Copy on-node surface currents */
#pragma omp parallel for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        for (int f = 0; f < HEX_NUM_FACES; f++)
        {

          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            /* Sum group contributions and add to currents */
            int bg = _cmfd_group_to_backup_group[e];
            CMFD_PRECISION val =
                _surface_currents->getValue(i, f * _num_cmfd_groups + e);
            backup_currents->incrementValue(i, f * nbg + bg, val);
          }
        }
      }

#ifdef ENABLE_MPI_
      /* Copy off-node surface currents */
      if (_domain_communicator != NULL)
      {

        CMFD_PRECISION ***off_node_currents =
            _backup_cmfd->getBoundarySurfaceCurrents();

        for (int surface = 0; surface < HEX_NUM_FACES; surface++)
        {

          /* Extract arrays on surface */
          CMFD_PRECISION **boundary_currents = _boundary_surface_currents[surface];
          CMFD_PRECISION **backup_currents = off_node_currents[surface];

          /* Loop over all CMFD cells on the current surface */
          std::map<int, int>::iterator it;
          for (it = _boundary_index_map.at(surface).begin();
               it != _boundary_index_map.at(surface).end(); ++it)
          {

            int idx = it->second;

            /* Loop over cell faces */
            for (int f = 0; f < HEX_NUM_FACES; f++)
            {

              /* Loop over CMFD coarse energy groups */
              for (int e = 0; e < _num_cmfd_groups; e++)
              {
                int bg = _cmfd_group_to_backup_group[e];
                backup_currents[idx][f * nbg + bg] +=
                    boundary_currents[idx][f * _num_cmfd_groups + e];
              }
            }
          }
        }
      }
#endif
    }
    else
    {
      /* Copy on-node surface currents */
#pragma omp parallel for
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        for (int f = 0; f < NUM_FACES; f++)
        {

          for (int e = 0; e < _num_cmfd_groups; e++)
          {

            /* Sum group contributions and add to currents */
            int bg = _cmfd_group_to_backup_group[e];
            CMFD_PRECISION val =
                _surface_currents->getValue(i, f * _num_cmfd_groups + e);
            backup_currents->incrementValue(i, f * nbg + bg, val);
          }
        }
      }

#ifdef ENABLE_MPI_
      /* Copy off-node surface currents */
      if (_domain_communicator != NULL)
      {

        CMFD_PRECISION ***off_node_currents =
            _backup_cmfd->getBoundarySurfaceCurrents();

        for (int surface = 0; surface < NUM_FACES; surface++)
        {

          /* Extract arrays on surface */
          CMFD_PRECISION **boundary_currents = _boundary_surface_currents[surface];
          CMFD_PRECISION **backup_currents = off_node_currents[surface];

          /* Loop over all CMFD cells on the current surface */
          std::map<int, int>::iterator it;
          for (it = _boundary_index_map.at(surface).begin();
               it != _boundary_index_map.at(surface).end(); ++it)
          {

            int idx = it->second;

            /* Loop over cell faces */
            for (int f = 0; f < NUM_FACES; f++)
            {

              /* Loop over CMFD coarse energy groups */
              for (int e = 0; e < _num_cmfd_groups; e++)
              {
                int bg = _cmfd_group_to_backup_group[e];
                backup_currents[idx][f * nbg + bg] +=
                    boundary_currents[idx][f * _num_cmfd_groups + e];
              }
            }
          }
        }
      }
#endif
    }
  }

  /**
   * @brief Returns the width of a given surface
   * @param surface A surface index, from 0 to NUM_FACES - 1
   * @return The surface width
   */
  CMFD_PRECISION Cmfd::getSurfaceWidth(int surface)
  {

    CMFD_PRECISION width;

    if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
      width = _cell_width_y * _cell_width_z;
    else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX)
      width = _cell_width_x * _cell_width_z;
    else
      width = _cell_width_x * _cell_width_y;

    return width;
  }

  /**
   * @brief Returns the width of a given surface
   * @param surface A surface index, from 0 to HEX_NUM_FACES - 1
   * @return The surface width
   */
  CMFD_PRECISION Cmfd::getHexSurfaceWidth(int surface)
  {

    /* width实际代指的是surface对应的面的面积 */
    CMFD_PRECISION width;
    double r = _width_r;
    double cos30 = std::cos(M_PI / 6);
    double sin30 = 0.5;
    double length = r / (2 * cos30);

    if (surface == HEX_SURFACE_BETA_MIN || surface == HEX_SURFACE_BETA_MAX ||
        surface == HEX_SURFACE_GAMMA_MIN || surface == HEX_SURFACE_GAMMA_MAX ||
        surface == HEX_SURFACE_DELTA_MIN || surface == HEX_SURFACE_DELTA_MAX)
      width = length * _cell_width_z;
    else
      width = 1.5 * length * r;

    return width;
  }

  /**
   * @brief Returns the width of the surface perpendicular to a given surface
   * @param surface A surface index, from 0 to NUM_FACES - 1
   * @return The perpendicular surface width
   */
  CMFD_PRECISION Cmfd::getPerpendicularSurfaceWidth(int surface)
  {

    if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
      return _cell_width_x;
    else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX)
      return _cell_width_y;
    else
      return _cell_width_z;
  }

  /**
   * @brief Returns the width of the surface perpendicular to a given surface
   * @param surface A surface index, from 0 to NUM_FACES - 1
   * @return The perpendicular surface width
   */
  CMFD_PRECISION Cmfd::getHexPerpendicularSurfaceWidth(int surface)
  {

    if (surface == HEX_SURFACE_BETA_MIN || surface == HEX_SURFACE_BETA_MAX ||
        surface == HEX_SURFACE_GAMMA_MIN || surface == HEX_SURFACE_GAMMA_MAX ||
        surface == HEX_SURFACE_DELTA_MIN || surface == HEX_SURFACE_DELTA_MAX)
      return _width_r;
    else
      return _cell_width_z;
  }

  /**
   * @brief Returns the sense of a given surface
   * @details The sense of minimum surfaces (e.g. SURFACE_X_MIN) is defined to be
   *          -1 while maximum surfaces (e.g. SURFACE_X_MAX) are defined to have a
   *          sense of +1. This is based on the current exiting a cell from a
   *          minimum surface being in the direction of negative net current and
   *          the current leaving a cell from a maximum surface being in the
   *          direction of positive net current.
   * @param surface A surface index, from 0 to NUM_FACES - 1
   * @return The sense of the surface
   */
  int Cmfd::getSense(int surface)
  {
    if (_hexlattice_enable)
    {
      if (surface == HEX_SURFACE_BETA_MIN || surface == HEX_SURFACE_GAMMA_MIN ||
          surface == HEX_SURFACE_DELTA_MIN || surface == HEX_SURFACE_Z_MIN)
        return -1;
      else
        return 1;
    }
    else
    {
      if (surface == SURFACE_X_MIN || surface == SURFACE_Y_MIN ||
          surface == SURFACE_Z_MIN)
        return -1;
      else
        return 1;
    }
  }

  /**
   * @brief Sets a flag to indicate whether a 2D or 3D problem is being solved.
   * @param solve_3D A boolean indicate whether a 2D or 3D problem is being
   *        solved.
   */
  void Cmfd::setSolve3D(bool solve_3D)
  {
    _SOLVE_3D = solve_3D;
  }

  /**
   * @brief Sets the azimuthal spacings.
   * @param azim_spacings An array of azimuthal spacings for each azimuthal angle.
   * @param num_azim the number of azimuthal angles.
   */
  void Cmfd::setAzimSpacings(const std::vector<double> &azim_spacings,
                             int num_azim)
  {

    if (_azim_spacings != NULL)
      delete[] _azim_spacings;

    _azim_spacings = new double[num_azim / 4];

    for (int a = 0; a < num_azim / 4; a++)
      _azim_spacings[a] = double(azim_spacings[a]);
  }

  /**
   * @brief Sets the polar spacings.
   * @param polar_spacings A 2D array of polar spacings for each azimuthal and
   *        polar angle combination. 每个方位角和极角的组合
   * @param num_azim the number of azimuthal angles.
   * @param num_polar the number of polar angles.
   */
  void Cmfd::setPolarSpacings(const std::vector<std::vector<double>> &
                                  polar_spacings,
                              int num_azim, int num_polar)
  {

    if (_polar_spacings != NULL)
    {
      for (int a = 0; a < num_azim / 4; a++)
        delete[] _polar_spacings[a];
      delete[] _polar_spacings;
    }

    _polar_spacings = new double *[num_azim / 4];
    for (int a = 0; a < num_azim / 4; a++)
      _polar_spacings[a] = new double[num_polar / 2];

    for (int a = 0; a < num_azim / 4; a++)
    {
      for (int p = 0; p < num_polar / 2; p++)
        _polar_spacings[a][p] = double(polar_spacings[a][p]);
    }
  }

  /**
   * @brief Set the value of the k effective for the CMFD solver. This is meant
   *        for research / debugging purposes.
   * @param k_eff the k_eff value to set.
   */
  void Cmfd::setKeff(double k_eff)
  {
    _k_eff = k_eff;
  }

  /**
   * @brief Set the backup CMFD solver's group structure. It is necessarily
   *        coarser than and must align with the regular CMFD group structure.
   * @param group_indices the indices of the CMFD groups in the MOC groups
   */
  void Cmfd::setBackupGroupStructure(std::vector<std::vector<int>>
                                         group_indices)
  {

    /* Assign the number of backup energy groups */
    _num_backup_groups = group_indices.size();

    /* Initialize mappings */
    _backup_group_structure = group_indices;
    _cmfd_group_to_backup_group = new int[_num_cmfd_groups];
    for (int e = 0; e < _num_cmfd_groups; e++)
      _cmfd_group_to_backup_group[e] = -1;

    /* Check that the mapping is valid and assign CMFD groups to backup groups */
    int cmfd_group = -1;
    int moc_group = 0;
    for (long i = 0; static_cast<size_t>(i) < group_indices.size(); i++)
    {
      for (long j = 0; static_cast<size_t>(j) < group_indices.at(i).size(); j++)
      {
        if (group_indices.at(i).at(j) != moc_group + 1)
        {
          log::ferror("Invalid backup group structure: indices must be "
                      "monotonic and include all MOC groups.");
        }
        if (moc_group >= _group_indices[cmfd_group + 1])
        {
          cmfd_group++;
          _cmfd_group_to_backup_group[cmfd_group] = i;
        }

        if (i != _cmfd_group_to_backup_group[cmfd_group])
          log::ferror("Invalid backup group structure: indices of backup "
                      "group structure must align with boundaries of CMFD group "
                      "structure.");

        moc_group++;
      }
    }

    /* Ensure that every CMFD group has a backup group */
    for (int e = 0; e < _num_cmfd_groups; e++)
    {
      if (_cmfd_group_to_backup_group[e] == -1)
        log::ferror("Invalid backup group structure: failed to find "
                    "matching index for CMFD group %d",
                    e);
    }
  }

  /**
   * @brief A function that prints a summary of the CMFD input parameters.
   */
  void Cmfd::printInputParamsSummary()
  {

    if (_flux_update_on)
      log::finfo("CMFD acceleration: ON");
    else
      log::finfo("CMFD acceleration: OFF (no MOC flux update)");

    // Print CMFD relaxation information
    if (std::abs(_SOR_factor - 1) > FLT_EPSILON)
      log::finfo("CMFD inner linear solver SOR factor: %f", _SOR_factor);
    log::finfo("CMFD corrected diffusion coef. relaxation factor: %f",
               _relaxation_factor);

    // Print CMFD interpolation techniques
    if (_centroid_update_on)
      log::finfo("CMFD K-nearest scheme: %d neighbors", _k_nearest);
    if (_use_axial_interpolation == 1)
      log::finfo("CMFD axial interpolation with axially averaged update "
                 "ratios");
    else if (_use_axial_interpolation == 2)
      log::finfo("CMFD axial interpolation with update ratios evaluated "
                 "at centroid Z-coordinate");

    // Print other CMFD modifications
    if (_flux_limiting)
      log::fverbose_once("CMFD corrected diffusion coef. bounded by "
                         "regular diffusion coef.");
    if (_balance_sigma_t)
      log::fverbose_once("CMFD total cross sections adjusted for matching MOC "
                         "reaction rates");

    // Print CMFD space and energy mesh information
    log::finfo("CMFD Mesh: %d x %d x %d", _num_x, _num_y, _num_z);
    if (!_hexlattice_enable)
    {
      if (_num_cmfd_groups != _num_moc_groups)
      {
        log::finfo("CMFD Group Structure:");
        log::finfo("\t MOC Group \t CMFD Group");
        for (int g = 0; g < _num_moc_groups; g++)
          log::finfo("\t %d \t\t %d", g + 1, getCmfdGroup(g) + 1);
      }
      else
        log::finfo("CMFD and MOC group structures match");
    }
  }

  /**
   * @brief Report the physical time use by major components of the CMFD solver.
   */
  void Cmfd::printTimerReport()
  {

    std::string msg_string;

    /* Get the total CMFD time */
    double tot_time = _timer->getSplit("Total CMFD time");
    msg_string = "Total CMFD computation time";
    msg_string.resize(53, '.');
    log::fresult("%s%1.4E sec", msg_string.c_str(), tot_time);

    /* Get the total XS collapse time */
    double xs_collapse_time = _timer->getSplit("Total collapse time");
    msg_string = "Total XS collapse time";
    msg_string.resize(53, '.');
    log::fresult("%s%1.4E sec", msg_string.c_str(), xs_collapse_time);

    /* Get the MPI communication time */
    double comm_time = _timer->getSplit("CMFD MPI communication time");
    msg_string = "CMFD MPI communication time";
    msg_string.resize(53, '.');
    log::fresult("%s%1.4E sec", msg_string.c_str(), comm_time);

    /* Get the total solver time */
    double solver_time = _timer->getSplit("Total solver time");
    msg_string = "Total CMFD solver time";
    msg_string.resize(53, '.');
    log::fresult("%s%1.4E sec", msg_string.c_str(), solver_time);
  }

  /**
   * @brief Forms a full copy of the surface currents on every surface
   * @details The copy contains all surface currents including edge and corner
   *          currents explicitly. It is stored in the _full_surface_currents
   *          vector for use in debugging and diagnostics.
   */
  void Cmfd::copyFullSurfaceCurrents()
  {
    log::finfo("Forms a full copy of the surface currents on every surface");

    /* Allocate full surface currents if necessary */
    if (_full_surface_currents == NULL)
      _full_surface_currents = new Vector(_cell_locks, _local_num_x,
                                          _local_num_y, _local_num_z,
                                          _num_cmfd_groups * NUM_SURFACES);

    /* Clear the currently saved surface currents */
    _full_surface_currents->clear();

    /* Copy surface currents from surface faces */
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {
      for (int s = 0; s < NUM_FACES; s++)
      {
        for (int g = 0; g < _num_cmfd_groups; g++)
        {
          FP_PRECISION current =
              _surface_currents->getValue(i, s * _num_cmfd_groups + g);
          _full_surface_currents->incrementValue(i, s * _num_cmfd_groups + g,
                                                 current);
        }
      }
    }

    /* Copy surface currents from edges and corners */
    std::map<int, CMFD_PRECISION>::iterator it;
    for (it = _edge_corner_currents.begin();
         it != _edge_corner_currents.end(); ++it)
    {
      int key = it->first;
      int cell = key / (_num_cmfd_groups * NUM_SURFACES);
      int surf_group = key % (_num_cmfd_groups * NUM_SURFACES);
      _full_surface_currents->incrementValue(cell, surf_group, it->second);
    }
  }

  /**
   * @brief Forms a full copy of the hex surface currents on every hex surface
   * @details The copy contains all surface currents including edge and corner
   *          currents explicitly. It is stored in the _full_surface_currents
   *          vector for use in debugging and diagnostics.
   */
  void Cmfd::copyHexFullSurfaceCurrents()
  {
    log::finfo("Forms a full copy of the surface currents on every surface");

    /* Allocate full surface currents if necessary */
    if (_full_surface_currents == NULL)
      _full_surface_currents = new Vector(_cell_locks, _local_num_x,
                                          _local_num_y, _local_num_z,
                                          _num_cmfd_groups * HEX_NUM_SURFACES,
                                          _empty_cells_num);

    /* Clear the currently saved surface currents */
    _full_surface_currents->clear();

    /* Copy surface currents from surface faces */
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {
      if (_empty_fsrs_cells[i])
        continue;
      else
      {
        int j = _logical_actual_map[i];
        for (int s = 0; s < HEX_NUM_FACES; s++)
        {
          for (int g = 0; g < _num_cmfd_groups; g++)
          {
            FP_PRECISION current =
                _surface_currents->getValue(j, s * _num_cmfd_groups + g);
            _full_surface_currents->incrementValue(j, s * _num_cmfd_groups + g,
                                                   current);
          }
        }
      }
    }

    /* Copy surface currents from edges and corners */
    std::map<int, CMFD_PRECISION>::iterator it;
    for (it = _edge_corner_currents.begin();
         it != _edge_corner_currents.end(); ++it)
    {
      int key = it->first;
      int cell = key / (_num_cmfd_groups * HEX_NUM_SURFACES);
      if (!_empty_fsrs_cells[cell])
      {
        cell = _logical_actual_map[cell];
        int surf_group = key % (_num_cmfd_groups * HEX_NUM_SURFACES);
        _full_surface_currents->incrementValue(cell, surf_group, it->second);
      }
    }
  }

  /**
   * @brief Computes the neutron balance over each CMFD cell for both MOC and CMFD
   * @details This routine can be used once the CMFD matrices have been formed to
   *          compute the neutron balance in the CMFD cell. With regards to MOC,
   *          it loops over all fsrs in the cell to compute all reaction rates and
   *          currents.
   */
  void Cmfd::checkNeutronBalance(bool pre_split)
  {

    /* Initialize variables */
    omp_lock_t *cell_locks = _old_flux->getCellLocks();
    // unused
    // int num_rows = _old_flux->getNumRows();
    int num_x = _old_flux->getNumX();
    int num_y = _old_flux->getNumY();
    int num_z = _old_flux->getNumZ();
    int num_groups = _old_flux->getNumGroups();
    Vector m_phi(cell_locks, num_x, num_y, num_z, num_groups);
    Vector a_phi(cell_locks, num_x, num_y, num_z, num_groups);

    /* Compute CMFD balance */

    /* Compute neutron production */
    matrixMultiplication(_M, _old_flux, &m_phi);

    /* Compute neutron transfer and loss */
    matrixMultiplication(_A, _old_flux, &a_phi);
    CMFD_PRECISION *a_phi_array = a_phi.getArray();
    int *coupling_sizes = NULL;
    int **coupling_indexes = NULL;
    CMFD_PRECISION **coupling_coeffs = NULL;
    CMFD_PRECISION **coupling_fluxes = NULL;
    int offset = 0;
    for (int color = 0; color < 2; color++)
    {

#ifdef ENABLE_MPI_
      getCouplingTerms(_domain_communicator, color, coupling_sizes,
                       coupling_indexes, coupling_coeffs, coupling_fluxes,
                       _old_flux->getArray(), offset);
#endif

#pragma omp parallel for collapse(2)
      for (int iz = 0; iz < _local_num_z; iz++)
      {
        for (int iy = 0; iy < _local_num_y; iy++)
        {
          for (int ix = (iy + iz + color + offset) % 2; ix < _local_num_x; ix += 2)
          {
            int cell = (iz * _local_num_y + iy) * _local_num_x + ix;
            for (int g = 0; g < _num_cmfd_groups; g++)
            {

              int row = cell * _num_cmfd_groups + g;

              for (int i = 0; i < coupling_sizes[row]; i++)
              {
                int idx = coupling_indexes[row][i] * _num_cmfd_groups + g;
                int domain = _domain_communicator->domains[color][row][i];
                CMFD_PRECISION flux = coupling_fluxes[domain][idx];
                a_phi_array[row] += coupling_coeffs[row][i] * flux;
              }
            }
          }
        }
      }
    }

    double max_imbalance = 0.0;
    int max_imbalance_cell = -1;
    int max_imbalance_grp = -1;

    /* Loop over CMFD cells */
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {

      int x = (i % (_local_num_x * _local_num_y)) % _local_num_x;
      int y = (i % (_local_num_x * _local_num_y)) / _local_num_x;
      int z = i / (_local_num_x * _local_num_y);

      // unused
      // Material* cell_material = _materials[i];

      /* Loop over CMFD coarse energy groups */
      for (int e = 0; e < _num_cmfd_groups; e++)
      {

        /* Initialize tallies */
        double total = 0.0;
        double in_scattering = 0.0;
        double fission = 0.0;

        /* Loop over FSRs in CMFD cell */
        for (size_t j = 0; j < _cell_fsrs.at(i).size(); j++)
        {

          long fsr_id = _cell_fsrs.at(i).at(j);
          Material *fsr_material = _FSR_materials[fsr_id];
          FP_PRECISION volume = _FSR_volumes[fsr_id];
          FP_PRECISION *scat = fsr_material->getSigmaS();
          FP_PRECISION *flux = &_FSR_fluxes[fsr_id * _num_moc_groups];

          /* Loop over MOC energy groups within this CMFD coarse group */
          double chi = 0.0;
          for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
            chi += fsr_material->getChiByGroup(h + 1);

          /* Calculate total fission and in-scattering in the FSR */
          double tot_fission = 0.0;
          for (int g = 0; g < _num_moc_groups; g++)
          {

            /* Tally total fission */
            double nu_fis = fsr_material->getNuSigmaFByGroup(g + 1);
            tot_fission += nu_fis * flux[g] * volume;

            /* Loop over MOC energy groups within this CMFD coarse group */
            for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
              in_scattering += scat[h * _num_moc_groups + g] * flux[g] * volume;
          }

          /* Calculate fission contribution to this CMFD coarse group */
          fission += chi * tot_fission / _k_eff;

          /* Calcualte total reaction rate in this CMFD coarse group */
          for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
          {
            double tot = fsr_material->getSigmaTByGroup(h + 1);
            total += tot * flux[h] * volume;
          }
        }

        /* Calculate net current out of the cell */
        double net_current = 0.0;

        /* Use currents before splitting edges/corners if requested */
        if (pre_split)
        {

          /* Create arrays of cell indexes and bounds */
          int cell_limits[3] = {_local_num_x, _local_num_y, _local_num_z};
          int cell_ind[3];
          cell_ind[0] = i % _local_num_x;
          cell_ind[1] = (i / _local_num_x) % _local_num_y;
          cell_ind[2] = i / (_local_num_x * _local_num_y);

          /* Tally current from all surfaces including edges and corners */
          for (int s = 0; s < NUM_SURFACES; s++)
          {

            /* Compute index and vector direction */
            int idx = s * _num_cmfd_groups + e;
            int direction[3];
            convertSurfaceToDirection(s, direction);

            /* Copute the next CMFD cell from the cell indexes and direction */
            int cmfd_cell_next = 0;
            int cell_next_ind[3];
            for (int d = 0; d < 3; d++)
              cell_next_ind[d] = cell_ind[d] + direction[d];

            cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x + cell_next_ind[2] * (_local_num_x * _local_num_y);

            /* Compute the opposite direction vector */
            int op_direction[3];
            for (int d = 0; d < 3; d++)
              op_direction[d] = -1 * direction[d];

            /* Determine if the next CMFD cell is within the bounds */
            for (int d = 0; d < 3; d++)
              if (cell_next_ind[d] < 0 || cell_next_ind[d] >= cell_limits[d])
                cmfd_cell_next = -1;

            /* If the cell is outside the bounds, handle boundaries */
            if (cmfd_cell_next == -1)
            {

              /* Booleans for determining surface boundary type */
              bool vacuum = false;
              bool reflective = false;
              bool transmit_avail = false;

              int transmit_direction[3] = {0, 0, 0};

              /* Loop over all directions to handle boundaries */
              for (int d = 0; d < 3; d++)
              {
                if (cell_next_ind[d] < 0 || cell_next_ind[d] >= cell_limits[d])
                {

                  /* Form the surface for each direction */
                  int partial_direction[3] = {0, 0, 0};
                  partial_direction[d] = direction[d];
                  int partial_surface =
                      convertDirectionToSurface(partial_direction);

                  /* Look at the boundary type in this direction */
                  if (_boundaries[partial_surface] == VACUUM)
                  {
                    vacuum = true;
                  }
                  else if (_boundaries[partial_surface] == REFLECTIVE)
                  {
                    reflective = true;
                    op_direction[d] *= -1;
                  }
                }

                /* For non-boundary surfaces, save the direction */
                else if (direction[d] != 0)
                {
                  transmit_avail = true;
                  transmit_direction[d] = direction[d];
                }
              }

              /* For vacuum boundaries, tally the leakage */
              if (vacuum)
              {
                net_current += _full_surface_currents->getValue(i, idx);
              }

              /* For reflective boundaries, find the appropriate cell to
                 deliver current if available */
              else if (reflective && transmit_avail)
              {
                for (int d = 0; d < 3; d++)
                  cell_next_ind[d] = cell_ind[d] + transmit_direction[d];
                cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x + cell_next_ind[2] * (_local_num_x * _local_num_y);
              }
            }

            /* Transmit current to available cells */
            if (cmfd_cell_next != -1)
            {
              int surface_next = convertDirectionToSurface(op_direction);
              int idx_next = surface_next * _num_cmfd_groups + e;
              net_current += _full_surface_currents->getValue(i, idx) -
                             _full_surface_currents->getValue(cmfd_cell_next, idx_next);
            }
          }
        }

        /* Use post-split edges/corner currents if requested */
        else
        {

          /* Tally current only over surface faces */
          for (int s = 0; s < NUM_FACES; s++)
          {
            int idx = s * _num_cmfd_groups + e;
            int cmfd_cell_next = getCellNext(i, s);
            int surface_next = (s + NUM_FACES / 2) % NUM_FACES;
            int idx_next = surface_next * _num_cmfd_groups + e;
            if (cmfd_cell_next == -1)
            {
              if (_boundaries[s] == VACUUM)
              {
                net_current += _surface_currents->getValue(i, idx);
                int direction[3];
                convertSurfaceToDirection(s, direction);
              }
            }
            else
            {
              net_current += _surface_currents->getValue(i, idx) -
                             _surface_currents->getValue(cmfd_cell_next, idx_next);
              int direction[3];
              int op_direction[3];
              convertSurfaceToDirection(s, direction);
              convertSurfaceToDirection(surface_next, op_direction);
            }
          }
        }

        double moc_balance = in_scattering + fission - total - net_current;

        double cmfd_balance = m_phi.getValue(i, e) / _k_eff -
                              a_phi.getValue(i, e);

        double tmp_imbalance = std::max(std::abs(moc_balance),
                                        std::abs(cmfd_balance));
        if (tmp_imbalance > max_imbalance)
        {
          max_imbalance = tmp_imbalance;
          max_imbalance_cell = i;
          max_imbalance_grp = e;
        }

        if (std::abs(moc_balance - cmfd_balance) > 1e-14)
        {
          log::finfo("MOC neutron balance in cell (%d, %d, %d) for CMFD "
                     "group %d = %g",
                     x, y, z, e, moc_balance);

          log::finfo("CMFD neutron balance in cell (%d, %d, %d) for CMFD "
                     "group %d = %g",
                     x, y, z, e, cmfd_balance);

          log::finfo("Net neutron balance in cell (%d, %d, %d) for CMFD "
                     "group %d = %g",
                     x, y, z, e, moc_balance - cmfd_balance);
        }
      }
    }
    log::finfo("Maximum neutron imbalance of %g at cell %i and group "
               "%d.",
               max_imbalance, max_imbalance_cell, max_imbalance_grp);
  }

  void Cmfd::checkNeutronBalanceHex(bool pre_split)
  {

    /* Initialize variables */
    omp_lock_t *cell_locks = _old_flux->getCellLocks();
    // unused
    // int num_rows = _old_flux->getNumRows();
    int num_x = _old_flux->getNumX();
    int num_y = _old_flux->getNumY();
    int num_z = _old_flux->getNumZ();
    int num_groups = _old_flux->getNumGroups();
    Vector m_phi(cell_locks, num_x, num_y, num_z, num_groups, _empty_cells_num);
    Vector a_phi(cell_locks, num_x, num_y, num_z, num_groups, _empty_cells_num);

    /* Compute CMFD balance */

    /* Compute neutron production */
    matrixMultiplication(_M, _old_flux, &m_phi);

    /* Compute neutron transfer and loss */
    matrixMultiplication(_A, _old_flux, &a_phi);
    CMFD_PRECISION *a_phi_array = a_phi.getArray();
    int *coupling_sizes = NULL;
    int **coupling_indexes = NULL;
    CMFD_PRECISION **coupling_coeffs = NULL;
    CMFD_PRECISION **coupling_fluxes = NULL;
    int offset = 0;
    for (int color = 0; color < 2; color++)
    {

#ifdef ENABLE_MPI_
      getCouplingTerms(_domain_communicator, color, coupling_sizes,
                       coupling_indexes, coupling_coeffs, coupling_fluxes,
                       _old_flux->getArray(), offset);
#endif

#pragma omp parallel for collapse(2)
      for (int iz = 0; iz < _local_num_z; iz++)
      {
        for (int iy = 0; iy < _local_num_y; iy++)
        {
          for (int ix = (iy + iz + color + offset) % 2; ix < _local_num_x; ix += 2)
          {
            int cell = (iz * _local_num_y + iy) * _local_num_x + ix;
            if (!_empty_fsrs_cells[cell])
            {
              for (int g = 0; g < _num_cmfd_groups; g++)
              {

                int row = cell * _num_cmfd_groups + g;

                for (int i = 0; i < coupling_sizes[row]; i++)
                {
                  int idx = coupling_indexes[row][i] * _num_cmfd_groups + g;
                  int domain = _domain_communicator->domains[color][row][i];
                  CMFD_PRECISION flux = coupling_fluxes[domain][idx];
                  a_phi_array[row] += coupling_coeffs[row][i] * flux;
                }
              }
            }
          }
        }
      }
    }

    double max_imbalance = 0.0;
    int max_imbalance_cell = -1;
    int max_imbalance_grp = -1;

    /* Loop over CMFD cells */
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {

      int x = (i % (_local_num_x * _local_num_y)) % _local_num_x;
      int y = (i % (_local_num_x * _local_num_y)) / _local_num_x;
      int z = i / (_local_num_x * _local_num_y);

      // unused
      // Material* cell_material = _materials[i];

      /* Loop over CMFD coarse energy groups */
      for (int e = 0; e < _num_cmfd_groups; e++)
      {

        /* Initialize tallies */
        double total = 0.0;
        double in_scattering = 0.0;
        double fission = 0.0;

        /* Loop over FSRs in CMFD cell */
        for (size_t j = 0; j < _cell_fsrs.at(i).size(); j++)
        {

          long fsr_id = _cell_fsrs.at(i).at(j);
          Material *fsr_material = _FSR_materials[fsr_id];
          FP_PRECISION volume = _FSR_volumes[fsr_id];
          FP_PRECISION *scat = fsr_material->getSigmaS();
          FP_PRECISION *flux = &_FSR_fluxes[fsr_id * _num_moc_groups];

          /* Loop over MOC energy groups within this CMFD coarse group */
          double chi = 0.0;
          for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
            chi += fsr_material->getChiByGroup(h + 1);

          /* Calculate total fission and in-scattering in the FSR */
          double tot_fission = 0.0;
          for (int g = 0; g < _num_moc_groups; g++)
          {

            /* Tally total fission */
            double nu_fis = fsr_material->getNuSigmaFByGroup(g + 1);
            tot_fission += nu_fis * flux[g] * volume;

            /* Loop over MOC energy groups within this CMFD coarse group */
            for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
              in_scattering += scat[h * _num_moc_groups + g] * flux[g] * volume;
          }

          /* Calculate fission contribution to this CMFD coarse group */
          fission += chi * tot_fission / _k_eff;

          /* Calcualte total reaction rate in this CMFD coarse group */
          for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
          {
            double tot = fsr_material->getSigmaTByGroup(h + 1);
            total += tot * flux[h] * volume;
          }
        }

        /* Calculate net current out of the cell */
        double net_current = 0.0;

        /* Use currents before splitting edges/corners if requested */
        if (pre_split)
        {

          /* Create arrays of cell indexes and bounds */
          int cell_limits[3] = {_local_num_x, _local_num_y, _local_num_z};
          int cell_ind[3];
          cell_ind[0] = i % _local_num_x;
          cell_ind[1] = (i / _local_num_x) % _local_num_y;
          cell_ind[2] = i / (_local_num_x * _local_num_y);

          /* Tally current from all surfaces including edges and corners */
          for (int s = 0; s < HEX_NUM_SURFACES; s++)
          {

            /* Compute index and vector direction */
            int idx = s * _num_cmfd_groups + e;
            int direction[3];
            convertSurfaceToDirection(s, direction);

            /* Copute the next CMFD cell from the cell indexes and direction */
            int cmfd_cell_next = 0;
            int cell_next_ind[3];
            for (int d = 0; d < 3; d++)
              cell_next_ind[d] = cell_ind[d] + direction[d];

            cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x + cell_next_ind[2] * (_local_num_x * _local_num_y);

            /* Compute the opposite direction vector */
            int op_direction[3];
            for (int d = 0; d < 3; d++)
              op_direction[d] = -1 * direction[d];

            /* Determine if the next CMFD cell is within the bounds */
            for (int d = 0; d < 3; d++)
              if (cell_next_ind[d] < 0 || cell_next_ind[d] >= cell_limits[d])
                cmfd_cell_next = -1;

            /* If the cell is outside the bounds, handle boundaries */
            if (cmfd_cell_next == -1)
            {

              /* Booleans for determining surface boundary type */
              bool vacuum = false;
              bool reflective = false;
              bool transmit_avail = false;

              int transmit_direction[3] = {0, 0, 0};

              /* Loop over all directions to handle boundaries */
              for (int d = 0; d < 3; d++)
              {
                if (cell_next_ind[d] < 0 || cell_next_ind[d] >= cell_limits[d])
                {

                  /* Form the surface for each direction */
                  int partial_direction[3] = {0, 0, 0};
                  partial_direction[d] = direction[d];
                  int partial_surface =
                      convertDirectionToSurface(partial_direction);

                  /* Look at the boundary type in this direction */
                  if (_boundaries[partial_surface] == VACUUM)
                  {
                    vacuum = true;
                  }
                  else if (_boundaries[partial_surface] == REFLECTIVE)
                  {
                    reflective = true;
                    op_direction[d] *= -1;
                  }
                }

                /* For non-boundary surfaces, save the direction */
                else if (direction[d] != 0)
                {
                  transmit_avail = true;
                  transmit_direction[d] = direction[d];
                }
              }

              /* For vacuum boundaries, tally the leakage */
              if (vacuum)
              {
                net_current += _full_surface_currents->getValue(i, idx);
              }

              /* For reflective boundaries, find the appropriate cell to
                 deliver current if available */
              else if (reflective && transmit_avail)
              {
                for (int d = 0; d < 3; d++)
                  cell_next_ind[d] = cell_ind[d] + transmit_direction[d];
                cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x + cell_next_ind[2] * (_local_num_x * _local_num_y);
              }
            }

            /* Transmit current to available cells */
            if (cmfd_cell_next != -1)
            {
              int surface_next = convertDirectionToSurface(op_direction);
              int idx_next = surface_next * _num_cmfd_groups + e;
              net_current += _full_surface_currents->getValue(i, idx) -
                             _full_surface_currents->getValue(cmfd_cell_next, idx_next);
            }
          }
        }

        /* Use post-split edges/corner currents if requested */
        else
        {

          /* Tally current only over surface faces */
          for (int s = 0; s < HEX_NUM_FACES; s++)
          {
            int idx = s * _num_cmfd_groups + e;
            int cmfd_cell_next = getCellNext(i, s);
            int surface_next = (s + HEX_NUM_FACES / 2) % HEX_NUM_FACES;
            int idx_next = surface_next * _num_cmfd_groups + e;
            if (cmfd_cell_next == -1)
            {
              if (_boundaries[s] == VACUUM)
              {
                net_current += _surface_currents->getValue(i, idx);
                int direction[3];
                convertSurfaceToDirection(s, direction);
              }
            }
            else
            {
              net_current += _surface_currents->getValue(i, idx) -
                             _surface_currents->getValue(cmfd_cell_next, idx_next);
              int direction[3];
              int op_direction[3];
              convertSurfaceToDirection(s, direction);
              convertSurfaceToDirection(surface_next, op_direction);
            }
          }
        }

        double moc_balance = in_scattering + fission - total - net_current;

        double cmfd_balance = m_phi.getValue(i, e) / _k_eff -
                              a_phi.getValue(i, e);

        double tmp_imbalance = std::max(std::abs(moc_balance),
                                        std::abs(cmfd_balance));
        if (tmp_imbalance > max_imbalance)
        {
          max_imbalance = tmp_imbalance;
          max_imbalance_cell = i;
          max_imbalance_grp = e;
        }

        if (std::abs(moc_balance - cmfd_balance) > 1e-14)
        {
          log::finfo("MOC neutron balance in cell (%d, %d, %d) for CMFD "
                     "group %d = %g",
                     x, y, z, e, moc_balance);

          log::finfo("CMFD neutron balance in cell (%d, %d, %d) for CMFD "
                     "group %d = %g",
                     x, y, z, e, cmfd_balance);

          log::finfo("Net neutron balance in cell (%d, %d, %d) for CMFD "
                     "group %d = %g",
                     x, y, z, e, moc_balance - cmfd_balance);
        }
      }
    }
    log::finfo("Maximum neutron imbalance of %g at cell %i and group "
               "%d.",
               max_imbalance, max_imbalance_cell, max_imbalance_grp);
  }

  /**
   * @brief Returns the color of a CMFD cell in the red/black SOR solver
   * @param cmfd_cell The cmfd cell's global ID
   */
  int Cmfd::getCellColor(int cmfd_cell)
  {
    int ix = cmfd_cell % _num_x;
    int iy = (cmfd_cell % (_num_x * _num_y)) / _num_x;
    int iz = cmfd_cell / (_num_x * _num_y);
    int color = (ix + iy + iz) % 2;
    return color;
  }

  /**
   * @brief Packs reaction rates and currents into buffers for communication
   * 将每个网格单元的体积数据、反应率、扩散率和电流数据打包到相应的发送缓冲区中
   */
  void Cmfd::packBuffers()
  {

    int current_idx[6] = {0, 0, 0, 0, 0, 0};
    bool found_surfaces[NUM_FACES];

    for (int z = 0; z < _local_num_z; z++)
    {
      for (int y = 0; y < _local_num_y; y++)
      {
        for (int x = 0; x < _local_num_x; x++)
        {
          for (int s = 0; s < NUM_FACES; s++)
            found_surfaces[s] = false;
          if (x == 0)
            found_surfaces[SURFACE_X_MIN] = true;
          if (x == _local_num_x - 1)
            found_surfaces[SURFACE_X_MAX] = true;
          if (y == 0)
            found_surfaces[SURFACE_Y_MIN] = true;
          if (y == _local_num_y - 1)
            found_surfaces[SURFACE_Y_MAX] = true;
          if (z == 0)
            found_surfaces[SURFACE_Z_MIN] = true;
          if (z == _local_num_z - 1)
            found_surfaces[SURFACE_Z_MAX] = true;
          for (int s = 0; s < NUM_FACES; s++)
          {
            if (found_surfaces[s])
            {
              int idx = current_idx[s];
              int cell_id = ((z * _local_num_y) + y) * _local_num_x + x;
              _send_volumes[s][idx][0] = _volume_tally[cell_id][0];
              for (int e = 0; e < _num_cmfd_groups; e++)
              {
                _send_reaction[s][idx][e] = _reaction_tally[cell_id][e];
                _send_diffusion[s][idx][e] = _diffusion_tally[cell_id][e];
                for (int f = 0; f < NUM_FACES; f++)
                {
                  _send_currents[s][idx][f * _num_cmfd_groups + e] =
                      _surface_currents->getValue(cell_id, f * _num_cmfd_groups + e);
                }
              }
              current_idx[s]++;
            }
          }
        }
      }
    }
  }

/**
 * @brief Exchanges ghost cell buffers in 3D cartesian (i.e., 6 directions)
 * 三维笛卡尔网格计算中进行 ghost cell（幽灵单元）数据的交换。
 * 该函数使用 MPI 实现数据的并行传输，以确保计算域的边界和邻接区域数据的一致性
 * @param comm The cartesian MPI domain communicator object that is configured
 *        for the CMFD exchange
 * @param send_buffers A 2D array of floating point data. The outer dimension
 *        corresponds to each face of the domain,
 *        while the inner dimension is the serialized buffer corresponding to
 *        the number of 2D cells to exchange times the number of energy groups.
 * @param recv_buffers A 2D array of floating point data. The outer dimension
 *        corresponds to each face of the domain,
 *        while the inner dimension is the serialized buffer corresponding to
 *        the number of 2D cells to exchange times the number of energy groups.
 */
#ifdef ENABLE_MPI_
  void Cmfd::ghostCellExchange()
  {

    packBuffers();

    MPI_Request requests[2 * NUM_FACES]; // 为每个面准备两个 MPI 请求（一个发送请求和一个接收请求）

    MPI_Datatype precision;
    if (sizeof(CMFD_PRECISION) == 4)
      precision = MPI_FLOAT;
    else
      precision = MPI_DOUBLE;

    int storage_per_cell = ((2 + NUM_FACES) * _num_cmfd_groups + 1);

    int sizes[NUM_FACES]; // 每个面数据的大小
    for (int coord = 0; coord < 3; coord++)
    { // 遍历三个维度
      for (int d = 0; d < 2; d++)
      { // 最大面和最小面两个方向

        int dir = 2 * d - 1;          //-1/+1
        int surf = coord + 3 * d;     // 当前面索引
        int op_surf = surf - 3 * dir; // 对面索引
        int source, dest;

        // Figure out serialized buffer length for this face
        int size = 0;
        if (surf == SURFACE_X_MIN)
        {
          size = _local_num_y * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_X_MAX)
        {
          size = _local_num_y * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_Y_MIN)
        {
          size = _local_num_x * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_Y_MAX)
        {
          size = _local_num_x * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_Z_MIN)
        {
          size = _local_num_x * _local_num_y * storage_per_cell;
        }
        else if (surf == SURFACE_Z_MAX)
        {
          size = _local_num_x * _local_num_y * storage_per_cell;
        }

        sizes[surf] = size;

        MPI_Cart_shift(_domain_communicator->_MPI_cart, coord, dir, &source, &dest); // 确定数据交换的源和目标进程

        // Post send
        MPI_Isend(_send_data_by_surface[surf], size, precision,
                  dest, 0, _domain_communicator->_MPI_cart, &requests[2 * surf]);

        // Post receive
        MPI_Irecv(_domain_data_by_surface[op_surf], size, precision,
                  source, 0, _domain_communicator->_MPI_cart, &requests[2 * surf + 1]);
      }
    }

    // Block for communication round to complete  等待发送和接收完成
    bool round_complete = false;
    while (!round_complete)
    {

      round_complete = true;
      int flag;
      MPI_Status send_stat;
      MPI_Status recv_stat;

      for (int coord = 0; coord < 3; coord++)
      {
        for (int d = 0; d < 2; d++)
        {
          int surf = coord + 3 * d;
          // 使用 MPI_Test() 检查每个请求的完成状态。只有当所有的发送和接收操作都完成时，才会退出循环
          MPI_Test(&requests[2 * surf], &flag, &send_stat);
          if (flag == 0)
            round_complete = false;

          MPI_Test(&requests[2 * surf + 1], &flag, &recv_stat);
          if (flag == 0)
            round_complete = false;
        }
      }
    }
  }

  /**
   * @brief Communicate split (at corners and edges) currents (respectively edge
   *        and face currents) to other domains.
   * 用于在不同的域之间通信分裂的电流（即在角和边上的电流，分别是边电流和面电流）
   * @param faces whether the currents are for edges or faces, for unpacking
   * 当前通信的是面电流（face currents）还是边电流（edge currents）。用于解包时的处理
   */
  void Cmfd::communicateSplits(bool faces)
  {

    // TODO: Form into communicateEdgeCurrents and communicateFaceCurrents
    //  1. communicate edge currents use array of length NUM_EDGES, called after
    //     vertex splits
    //  2. communicate face currents use array of length NUM_FACES, called after
    //     edge splits
    //  NOTE: only communicate currents that are saved OFF DOMAIN at each step
    //  NOTE: communicateFaceCurrents will use currents formed from vertex splits

    MPI_Request requests[2 * NUM_FACES];

    MPI_Datatype precision;
    if (sizeof(CMFD_PRECISION) == 4)
      precision = MPI_FLOAT;
    else
      precision = MPI_DOUBLE;

    int storage_per_cell = (NUM_FACES + NUM_EDGES) * _num_cmfd_groups;

    int sizes[NUM_FACES]; // 用于存储每个面的大小
    /*
    根据每个面的尺寸（SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURFACE_Y_MAX, SURFACE_Z_MIN, SURFACE_Z_MAX）
    计算序列化缓冲区的大小。
    */
    for (int coord = 0; coord < 3; coord++)
    { // 遍历 X, Y, Z 三个方向的坐标
      for (int d = 0; d < 2; d++)
      { // 遍历每个方向的两个面（最小面和最大面）

        int dir = 2 * d - 1;          // d=0 时 dir=-1，d=1 时 dir=1
        int surf = coord + 3 * d;     // 当前面的索引
        int op_surf = surf - 3 * dir; // 对面的索引
        int source, dest;

        // Figure out serialized buffer length for this face  确定每个面的序列化缓冲区大小
        int size = 0;
        if (surf == SURFACE_X_MIN)
        {
          size = _local_num_y * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_X_MAX)
        {
          size = _local_num_y * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_Y_MIN)
        {
          size = _local_num_x * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_Y_MAX)
        {
          size = _local_num_x * _local_num_z * storage_per_cell;
        }
        else if (surf == SURFACE_Z_MIN)
        {
          size = _local_num_x * _local_num_y * storage_per_cell;
        }
        else if (surf == SURFACE_Z_MAX)
        {
          size = _local_num_x * _local_num_y * storage_per_cell;
        }

        sizes[surf] = size;

        MPI_Cart_shift(_domain_communicator->_MPI_cart, coord, dir, &source, &dest); // 使用 MPI_Cart_shift 函数计算通信的源和目标进程

        // Post send  使用 MPI_Isend 和 MPI_Irecv 非阻塞地发送和接收数据
        MPI_Isend(_send_split_currents_array[surf], size, precision,
                  dest, 0, _domain_communicator->_MPI_cart, &requests[2 * surf]);

        // Post receive
        MPI_Irecv(_receive_split_currents_array[op_surf], size, precision,
                  source, 0, _domain_communicator->_MPI_cart, &requests[2 * surf + 1]);
      }
    }

    // Block for communication round to complete  等待所有通信完成
    bool round_complete = false;
    while (!round_complete)
    {

      round_complete = true;
      int flag;
      MPI_Status send_stat;
      MPI_Status recv_stat;

      for (int coord = 0; coord < 3; coord++)
      {
        for (int d = 0; d < 2; d++)
        {
          int surf = coord + 3 * d;

          MPI_Test(&requests[2 * surf], &flag, &send_stat); // 检查发送请求 requests[2 * surf] 是否完成
          if (flag == 0)
            round_complete = false;

          MPI_Test(&requests[2 * surf + 1], &flag, &recv_stat); // 检查接收请求 requests[2 * surf + 1] 是否完成
          if (flag == 0)
            round_complete = false;
        }
      }
    }

    unpackSplitCurrents(faces);
  }
#endif

  /**
   * @brief Unpacks communicated split current data
   * @param faces Whether to split the currents onto surface faces
   * 解包已经通信完成的分割电流数据，并将其应用到对应的 CMFD单元中
   */
  void Cmfd::unpackSplitCurrents(bool faces)
  {

    int current_idx[6] = {0, 0, 0, 0, 0, 0}; // 跟踪每个面的当前索引
    bool found_surfaces[NUM_FACES];          // 用于标记当前单元是否在每个面上

    /* Loop over all CMFD cells */
    for (int z = 0; z < _local_num_z; z++)
    {
      for (int y = 0; y < _local_num_y; y++)
      {
        for (int x = 0; x < _local_num_x; x++)
        {

          /* Look for boundaries touching the CMFD cell 标记当前单元所属的面*/
          for (int s = 0; s < NUM_FACES; s++)
            found_surfaces[s] = false;
          if (x == 0)
            found_surfaces[SURFACE_X_MIN] = true;
          if (x == _local_num_x - 1)
            found_surfaces[SURFACE_X_MAX] = true;
          if (y == 0)
            found_surfaces[SURFACE_Y_MIN] = true;
          if (y == _local_num_y - 1)
            found_surfaces[SURFACE_Y_MAX] = true;
          if (z == 0)
            found_surfaces[SURFACE_Z_MIN] = true;
          if (z == _local_num_z - 1)
            found_surfaces[SURFACE_Z_MAX] = true;

          /* Handle all boundaries 处理每个边界*/
          for (int s = 0; s < NUM_FACES; s++)
          {
            if (found_surfaces[s])
            {

              /* Convert the (x,y,z) indexes to a cell ID and boundary index */
              int cell_id = ((z * _local_num_y) + y) * _local_num_x + x; // 域局部一维索引
              int idx = current_idx[s];

              /* Copy the appropriate face or edge information */
              if (faces)
              { // 若 faces 为 true，解包面电流数据并将非零值累加到 cell_id 对应的 _surface_currents 中

                /* Treat CMFD cell face currents */
                for (int f = 0; f < NUM_FACES; f++)
                {
                  for (int g = 0; g < _num_cmfd_groups; g++)
                  {

                    /* Get the face current value */
                    CMFD_PRECISION value =
                        _received_split_currents[s][idx][f * _num_cmfd_groups + g];

                    /* Treat nonzero values */
                    if (fabs(value) > FLT_EPSILON)
                      _surface_currents->incrementValue(cell_id,
                                                        f * _num_cmfd_groups + g,
                                                        value); // 添加接收到的值
                  }
                }
              }
              else
              { // 若 faces 为 false，则解包边缘电流数据

                /* Treat CMFD cell edge currents */
                for (int e = NUM_FACES; e < NUM_EDGES + NUM_FACES; e++)
                {

                  int surf_idx = cell_id * NUM_SURFACES * _num_cmfd_groups + e *
                                                                                 _num_cmfd_groups;

                  for (int g = 0; g < _num_cmfd_groups; g++)
                  {

                    /* Get the edge current value */
                    CMFD_PRECISION value =
                        _received_split_currents[s][idx][e * _num_cmfd_groups + g];

                    /* Treat nonzero values 对于非零值，
                    首先检查是否在 _edge_corner_currents 中存在相应索引。如果不存在则初始化为零，然后将电流值累加到对应位置*/
                    if (fabs(value) > FLT_EPSILON)
                    {

                      /* Check for new index in map */
                      int new_ind = surf_idx + g;
                      std::map<int, CMFD_PRECISION>::iterator it =
                          _edge_corner_currents.find(new_ind);

                      /* If it doesn't exist, initialize to zero */
                      if (it == _edge_corner_currents.end())
                        _edge_corner_currents[new_ind] = 0.0;

                      /* Add the contribution */
                      _edge_corner_currents[new_ind] += value;
                    }
                  }
                }
              }

              /* Increment the boundary index */
              current_idx[s]++;
            }
          }
        }
      }
    }
  }

  /**
   * @brief Converts a global CMFD cell ID into its local ID
   * @details Marked for deletion, but still used thoroughly.
   * @param cmfd_cell The global CMFD cell ID
   * @return The local CMFD cell ID, -1 if not in the domain.
   */
  int Cmfd::getLocalCMFDCell(int cmfd_cell)
  {

    int x_start = 0;
    int y_start = 0;
    int z_start = 0;
    int x_end = _num_x;
    int y_end = _num_y;
    int z_end = _num_z;
    if (mpi::isSpatialDecomposed())
    {
      if (_domain_communicator != NULL)
      {
        x_start = _domain_communicator->_domain_idx_x * _local_num_x;
        x_end = x_start + _local_num_x;
        y_start = _domain_communicator->_domain_idx_y * _local_num_y;
        y_end = y_start + _local_num_y;
        z_start = _domain_communicator->_domain_idx_z * _local_num_z;
        z_end = z_start + _local_num_z;
      }
    }

    int ix = (cmfd_cell % (_num_x * _num_y)) % _num_x; // 获取全局三维索引
    int iy = (cmfd_cell % (_num_x * _num_y)) / _num_x;
    int iz = cmfd_cell / (_num_x * _num_y);

    int local_cmfd_cell;
    if (ix < x_start || ix >= x_end || iy < y_start || iy >= y_end ||
        iz < z_start || iz >= z_end)
    {
      local_cmfd_cell = -1;
    }
    else // 转为局部一维索引
      local_cmfd_cell = ((iz - z_start) * _local_num_y + iy - y_start) * _local_num_x + ix - x_start;
    return local_cmfd_cell;
  }

  /**
   * @brief Converts a local CMFD cell ID into its global ID
   * @param cmfd_cell The local CMFD cell ID
   * @return The global CMFD cell ID
   *
   * 把一个“ MPI 子域里的 CMFD 单元索引”换算成全局单元索引。
   * 这样得到的整数就是该单元在整个 CMFD 网格中的唯一编号，和 MPI 分解无关。
   */
  int Cmfd::getGlobalCMFDCell(int cmfd_cell)
  {

    int x_start = 0;
    int y_start = 0;
    int z_start = 0;
    /*
    初始为 0；如果 mpi::isSpatialDecomposed()（说明 CMFD 网格被分块到不同 MPI 进程），并且 _domain_communicator 存在，就用该域在 X/Y/Z 方向的索引号乘以每个域本地单元数 _local_num_x 等，得到本域在全局坐标上的起始偏移。
    */
    if (mpi::isSpatialDecomposed())
    {
      if (_domain_communicator != NULL)
      {
        x_start = _domain_communicator->_domain_idx_x * _local_num_x;
        y_start = _domain_communicator->_domain_idx_y * _local_num_y;
        z_start = _domain_communicator->_domain_idx_z * _local_num_z;
      }
    }

    int ix = cmfd_cell % _local_num_x;                                   // 计算 cmfd_cell 在 x ,y,z方向上的局部索引
    int iy = (cmfd_cell % (_local_num_x * _local_num_y)) / _local_num_x; // 先取 XY 平面内的偏移，再除以 _local_num_x 得到 Y 位置。
    int iz = cmfd_cell / (_local_num_x * _local_num_y);

    return ((iz + z_start) * _num_y + iy + y_start) * _num_x + ix + x_start;
  }

  /**
   * @brief Converts a 3 integer vector direction into a surface
   * @details The direction is a tuplet with each value taking either
   *          +1 (positive directed), 0 (neutral, or -1 (negative directed)
   * @param direction The integer vector describing the direction
   * @return The surface associated with traveling the provided direction from
   *         the origin of the cell
   */
  int Cmfd::convertDirectionToSurface(int *direction)
  {
    int surface = 0;
    int num_crossings = std::abs(direction[0]) + std::abs(direction[1]) +
                        std::abs(direction[2]); // 交叉数量,1一个就是面,2个是边,3个是顶点
    if (num_crossings == 1)
    { // 即方向向量中只有一个非零值，其余两个值为0,也就是一个面
      for (int i = 0; i < 3; i++)
      {
        int present = std::abs(direction[i]); // direction如果为0,surface仍为0,即直到找到那个非零的方向
        int fwd = (direction[i] + 1) / 2;     // 计算当前方向是正方向（1）还是负方向（0）
        surface += present * (3 * fwd + i);   // 如果是负方向(-1)，编号增加 i,正方向编号增加 3 + i(比如SURFACE_Y_MIN =1,SURFACE_Y_MAX =4)
      }
    }
    else if (num_crossings == 2)
    {                       // 处理边的情况,边只有1个为0,其余两个非零
      surface += NUM_FACES; // 首先+6,是因为边最小值就是6,从6开始的
      int ind1 = 0;
      int ind2 = 0; // ind1,ind2确定该边是由哪两个方向的面交叉而来的
      if (direction[0] == 0)
      { // 如果x方向为零,说明是y和z方向交叉的
        ind1 = direction[1];
        ind2 = direction[2];
        surface += 8; // 因为涉及y和z的边是从14开始的,所以要先+8;
      }
      else if (direction[1] == 0)
      { // 如果y方向为零,说明是x和z方向交叉的
        ind1 = direction[0];
        ind2 = direction[2];
        surface += 4; // 因为涉及x和z的边是从10开始的,所以要先+4;
      }
      else if (direction[2] == 0)
      {
        ind1 = direction[0];
        ind2 = direction[1]; // 因为涉及x和y的边是从6开始的,所以不需要加了;
      }
      ind1 = (ind1 + 1) / 2; // 将方向值从 -1、-1 转换为 0 和 1
      ind2 = (ind2 + 1) / 2;
      surface += 2 * ind2 + ind1; // 前面的ind1(边的第一个,可能为x或y)跨一个单位,后面的ind2(边的第二个,可能为y,z)跨两个单位,参考编号值
    }
    else if (num_crossings == 3)
    {                                   // 处理顶点的情况,三个全非零
      surface += NUM_FACES + NUM_EDGES; // 同上
      int fwd[3];
      for (int i = 0; i < 3; i++)
        fwd[i] = (direction[i] + 1) / 2;           // 将方向值从 -1、-1 转换为 0 和 1
      surface += 4 * fwd[0] + 2 * fwd[1] + fwd[2]; // x跨四个,y跨2,z跨1,参考编号值
    }
    else
    {
      log::ferror("Invalid number of surface crossings");
    }
    return surface;
  }

  /**
   * @brief Converts a surface into a 3 integer vector direction 将一个表面转换为一个3整数的方向向量
   * @details The direction is a tuplet with each value taking either
   *          +1 (positive directed), 0 (neutral, or -1 (negative directed)
   * @param surface The surface of interest
   * @param direction The integer vector describing the direction
   */
  void Cmfd::convertSurfaceToDirection(int surface, int *direction)
  {
    direction[0] = 0;
    direction[1] = 0;
    direction[2] = 0;

    if (surface < NUM_FACES)
    {
      int ind = surface % 3;           // 确定哪个轴方向
      int dir = 2 * (surface / 3) - 1; // 确定方向（+1或-1）
      direction[ind] = dir;            // 表面ID 0、1、2 对应 x、y、z 轴的负方向（即 -1），表面ID 3、4、5 对应 x、y、z 轴的正方向（即 +1）
    }
    else if (surface < NUM_FACES + NUM_EDGES)
    {
      surface -= NUM_FACES;    // 将表面ID转换为边ID
      int group = surface / 4; // 分组,每个组包含四条边
      int skipped = 2 - group; // 计算需要跳过的轴。由于有三个轴 (x, y, z)，对于每个 group，一个轴会被跳过，另外两个轴用于确定边的方向
      surface = surface % 4;   // 将surface其值调整到当前组内的范围，即0到3
      int ind[2];              // 该数组用于表示边在两个方向上的方向。ind[0] 表示在第一方向上的方向，ind[1] 表示在第二方向上的方向。
      ind[0] = surface % 2;
      ind[1] = (surface - ind[0]) / 2;
      int n = 0;
      for (int i = 0; i < 3; i++)
      { // 比如位于XY最小边 SURFACE_X_MIN_Y_MIN 6 第一行减去6,surface=0;
        if (i != skipped)
        {                                // 在三个轴中，跳过 skipped 轴，其他两个轴的方向由 ind 数组决定
          direction[i] = 2 * ind[n] - 1; // x轴方向direction[0] = -1,位于负方向,y轴方向direction[0] = -1 也位于负方向,Z轴跳过
          n++;
        }
      }
    }
    else if (surface < NUM_SURFACES)
    {                                       // 表面ID在NUM_FACES + NUM_EDGES和NUM_SURFACES之间（即它是一个顶点），计算该顶点对应的方向
      surface -= NUM_FACES + NUM_EDGES;     // 将表面ID转换为顶点ID
      direction[0] = 2 * (surface / 4) - 1; // 比如SURFACE_X_MIN_Y_MIN_Z_MIN = 18-18 = 0,推论正确
      direction[1] = 2 * ((surface / 2) % 2) - 1;
      direction[2] = 2 * (surface % 2) - 1;
    }
    else
    {
      log::ferror("Invalid surface ID %d", surface);
    }
  }

  void Cmfd::convertSurfaceToDirectionHex(int surface, int *direction)
  {
    direction[0] = 0;
    direction[1] = 0;
    direction[2] = 0;

    if (surface < HEX_NUM_FACES)
    {
      if (surface == HEX_SURFACE_BETA_MIN)
        direction[0] = -1;
      else if (surface == HEX_SURFACE_BETA_MAX)
        direction[0] = 1;
      else if (surface == HEX_SURFACE_GAMMA_MIN)
      {
        direction[0] = -1;
        direction[1] = 1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MAX)
      {
        direction[0] = 1;
        direction[1] = -1;
      }
      else if (surface == HEX_SURFACE_DELTA_MIN)
        direction[1] = -1;
      else if (surface == HEX_SURFACE_DELTA_MAX)
        direction[1] = 1;
      else if (surface == HEX_SURFACE_Z_MIN)
        direction[2] = -1;
      else if (surface == HEX_SURFACE_Z_MAX)
        direction[2] = 1;
    }
    else if (surface < HEX_NUM_FACES + HEX_NUM_EDGES)
    {
      if (surface == HEX_SURFACE_BETA_MIN_GAMMA_MIN)
      {
        direction[0] = -2;
        direction[1] = 1;
      }
      else if (surface == HEX_SURFACE_BETA_MIN_DELTA_MIN)
      {
        direction[0] = -1;
        direction[1] = -1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MIN_DELTA_MAX)
      {
        direction[0] = -1;
        direction[1] = 2;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_GAMMA_MAX)
      {
        direction[0] = 2;
        direction[1] = -1;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_DELTA_MAX)
      {
        direction[0] = 1;
        direction[1] = 1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MAX_DELTA_MIN)
      {
        direction[0] = 1;
        direction[1] = -2;
      }
      else if (surface == HEX_SURFACE_BETA_MIN_Z_MIN)
      {
        direction[0] = -1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_Z_MIN)
      {
        direction[0] = 1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_BETA_MIN_Z_MAX)
      {
        direction[0] = -1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_Z_MAX)
      {
        direction[0] = 1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MIN_Z_MIN)
      {
        direction[0] = -1;
        direction[1] = 1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MAX_Z_MIN)
      {
        direction[0] = 1;
        direction[1] = -1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MIN_Z_MAX)
      {
        direction[0] = -1;
        direction[1] = 1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MAX_Z_MAX)
      {
        direction[0] = 1;
        direction[1] = -1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_DELTA_MIN_Z_MIN)
      {
        direction[1] = -1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_DELTA_MAX_Z_MIN)
      {
        direction[1] = 1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_DELTA_MIN_Z_MAX)
      {
        direction[1] = -1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_DELTA_MAX_Z_MAX)
      {
        direction[1] = 1;
        direction[2] = 1;
      }
    }
    else if (surface < HEX_NUM_SURFACES)
    {
      if (surface == HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MIN)
      {
        direction[0] = -2;
        direction[1] = 1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MAX)
      {
        direction[0] = -2;
        direction[1] = 1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MIN)
      {
        direction[0] = -1;
        direction[1] = -1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MAX)
      {
        direction[0] = -1;
        direction[1] = -1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MIN)
      {
        direction[0] = 2;
        direction[1] = -1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MAX)
      {
        direction[0] = 2;
        direction[1] = -1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MIN)
      {
        direction[0] = 1;
        direction[1] = 1;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MAX)
      {
        direction[0] = 1;
        direction[1] = 1;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MIN)
      {
        direction[0] = -1;
        direction[1] = 2;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MAX)
      {
        direction[0] = -1;
        direction[1] = 2;
        direction[2] = 1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MIN)
      {
        direction[0] = 1;
        direction[1] = -2;
        direction[2] = -1;
      }
      else if (surface == HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MAX)
      {
        direction[0] = 1;
        direction[1] = -2;
        direction[2] = 1;
      }
    }
    else
    {
      log::ferror("Invalid surface ID %d", surface);
    }
  }

  /**
   * @brief Returns the surface name associated with the 3 integer vector
   *        direction
   * @details The direction is a tuplet with each value taking either
   *          +1 (positive directed), 0 (neutral, or -1 (negative directed)
   * @param direction The integer vector describing the direction
   * @return A string containing the surface name
   */
  std::string Cmfd::getSurfaceNameFromDirection(int *direction)
  {
    std::string str = "SURFACE";
    std::string variables = "XYZ";
    for (int i = 0; i < 3; i++)
    {
      if (direction[i] != 0)
      {
        str += "_";
        str += variables.at(i);
        if (direction[i] < 0)
          str += "_MIN";
        else
          str += "_MAX";
      }
    }
    return str;
  }

  /**
   * @brief Returns the surface name associated with a surface
   * @param surface The surface of interest
   * @return A string containing the surface name
   */
  std::string Cmfd::getSurfaceNameFromSurface(int surface)
  {
    int direction[3];
    convertSurfaceToDirection(surface, direction);
    return getSurfaceNameFromDirection(direction);
  }

  /**
   * @brief A debugging tool that prints all prolongation facotrs to file
   */
  void Cmfd::printProlongationFactors(int iteration)
  {

    /* Loop over CMFD groups */
    for (int e = 0; e < _num_cmfd_groups; e++)
    {

      /* Create arrays for spatial data */
      double *log_ratios = new double[_num_x * _num_y * _num_z];
      for (int i = 0; i < _num_x * _num_y * _num_z; i++)
        log_ratios[i] = 0.0;
      for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      {

        double old_flux = _old_flux->getValue(i, e);
        double new_flux = _new_flux->getValue(i, e);
        int cell_id = getGlobalCMFDCell(i);
        log_ratios[cell_id] = std::log(new_flux / old_flux);
      }

#ifdef ENABLE_MPI_
      if (mpi::isSpatialDecomposed())
      {
        double temp_log_ratios[_num_x * _num_y * _num_z];
        for (int i = 0; i < _num_x * _num_y * _num_z; i++)
          temp_log_ratios[i] = log_ratios[i];
        MPI_Allreduce(temp_log_ratios, log_ratios, _num_x * _num_y * _num_z,
                      MPI_DOUBLE, MPI_SUM, _geometry->getMPICart());
      }
#endif

      /* Print prolongation factors distribution to file */
      if (_geometry->isRootDomain())
      {
        long long iter = iteration;
        long long group = e;
        std::string fname = "pf_group_";
        std::string group_num = std::to_string(group);
        std::string iter_num = std::to_string(iter);
        fname += group_num;
        fname += "_iter_";
        fname += iter_num;
        std::ofstream out(fname);

        out << "[NORMAL]  Spatial distribution of prolongation factors:"
            << std::endl;
        for (int z = 0; z < _num_z; z++)
        {
          out << " -------- z = " << z << " ----------" << std::endl;
          for (int y = 0; y < _num_y; y++)
          {
            for (int x = 0; x < _num_x; x++)
            {
              int ind = (z * _num_y + y) * _num_x + x;
              out << log_ratios[ind] << " ";
            }
            out << std::endl;
          }
        }
        out.close();
      }
      delete[] log_ratios;
    }
  }

  /**
   * @brief Tallies the current contribution from this segment across the
   *        the appropriate CMFD mesh cell surface.
   * 用于在 CMFD网格单元表面上累积当前段（segment）的中子流贡献。
   * 这个函数通过处理段的流出角通量，计算并累积在每个网格单元表面的中子流
   * @param curr_segment the current Track segment
   * @param track_flux the outgoing angular flux for this segment 指向当前段的出射角通量的指针
   * @param polar_weights array of polar weights for some azimuthal angle
   * @param fwd boolean indicating direction of integration along segment
   */
  void Cmfd::tallyCurrent(segment *curr_segment, float *track_flux,
                          int azim_index, int polar_index, bool fwd)
  {

    int surf_id, cell_id, cmfd_group;
    int ncg = _num_cmfd_groups;
    int tid = omp_get_thread_num();
    CMFD_PRECISION *currents = _temporary_currents[tid];
    memset(currents, 0.0, sizeof(CMFD_PRECISION) * _num_cmfd_groups);
    std::map<int, CMFD_PRECISION>::iterator it;

    /* Check if the current needs to be tallied */
    bool tally_current = false;

    if (curr_segment->_cmfd_surface_fwd != -1 && fwd)
    { // 如果当前轨迹处于CMFD网格surface上且方向为正向
      surf_id = curr_segment->_cmfd_surface_fwd % NUM_SURFACES;
      cell_id = curr_segment->_cmfd_surface_fwd / NUM_SURFACES;
      tally_current = true;
    }
    else if (curr_segment->_cmfd_surface_bwd != -1 && !fwd)
    { // 如果当前轨迹处于CMFD网格surface上且方向为反向
      surf_id = curr_segment->_cmfd_surface_bwd % NUM_SURFACES;
      cell_id = curr_segment->_cmfd_surface_bwd / NUM_SURFACES;
      tally_current = true;
    }

    /* Tally current if necessary */
    if (tally_current)
    { // 如果tally_current仍为false,则不需要统计中子流因为该线段不处于CMFD的surface上

      int local_cell_id = getLocalCMFDCell(cell_id); // 取出域中的一维索引

      if (_SOLVE_3D)
      { // 三维求解
        double wgt = _quadrature->getWeightInline(azim_index, polar_index);
        for (int e = 0; e < _num_moc_groups; e++)
        {

          /* Get the CMFD group */
          cmfd_group = getCmfdGroup(e); // moc能群对应的cmfd能群

          /* Increment the surface group */
          currents[cmfd_group] += track_flux[e] * wgt;
        }

        /* Increment currents NUM_FACES = 6 */
        if (surf_id < NUM_FACES)
        {                                                                                                      // 统计面上的中子流
          _surface_currents->incrementValues(local_cell_id, surf_id * ncg, (surf_id + 1) * ncg - 1, currents); // 将结果存放在_array数组中
        }
        else
        { // 统计CMFD粗网网格边和顶点上的中子流

          omp_set_lock(&_edge_corner_lock);

          int first_ind = (local_cell_id * NUM_SURFACES + surf_id) * ncg; // 计算在 _edge_corner_currents 中的起始索引* ncg
          it = _edge_corner_currents.find(first_ind);                     // 查找起始索引 first_ind 是否已经存在于 _edge_corner_currents 中
          if (it == _edge_corner_currents.end())                          // 如果 it 等于 end()，则表示 first_ind 处没有找到元素,说明之前没有计算过,赋值为0.0
            for (int g = 0; g < ncg; g++)
              _edge_corner_currents[first_ind + g] = 0.0; // 初始化为 0

          for (int g = 0; g < ncg; g++)
            _edge_corner_currents[first_ind + g] += currents[g]; // 累积当前段的中子流贡献到 _edge_corner_currents

          omp_unset_lock(&_edge_corner_lock);
        }
      }
      else
      { // 如果是二维求解
        int pe = 0;
        for (int e = 0; e < _num_moc_groups; e++)
        {

          /* Get the CMFD group */
          cmfd_group = getCmfdGroup(e);

          for (int p = 0; p < _num_polar / 2; p++)
          { // 如果是二维,则把极角和能群的角通量全部累加起来
            currents[cmfd_group] += track_flux[pe] * _quadrature->getWeightInline(azim_index, p);
            pe++;
          }
        }

        /* Increment currents */
        if (surf_id < NUM_FACES)
        {
          _surface_currents->incrementValues(local_cell_id, surf_id * ncg, (surf_id + 1) * ncg - 1, currents);
        }
        else
        {
          omp_set_lock(&_edge_corner_lock);

          int first_ind = (local_cell_id * NUM_SURFACES + surf_id) * ncg;
          it = _edge_corner_currents.find(first_ind);
          if (it == _edge_corner_currents.end())
            for (int g = 0; g < ncg; g++)
              _edge_corner_currents[first_ind + g] = 0.0;

          for (int g = 0; g < ncg; g++)
            _edge_corner_currents[first_ind + g] += currents[g];

          omp_unset_lock(&_edge_corner_lock);
        }
      }
    }
  }

  void Cmfd::tallyHexCurrent(segment *curr_segment, float *track_flux,
                             int azim_index, int polar_index, bool fwd)
  {
    int surf_id, cell_id, cmfd_group;
    int ncg = _num_cmfd_groups;
    int tid = omp_get_thread_num();
    CMFD_PRECISION *currents = _temporary_currents[tid];
    memset(currents, 0.0, sizeof(CMFD_PRECISION) * _num_cmfd_groups);
    std::map<int, CMFD_PRECISION>::iterator it;

    /* Check if the current needs to be tallied */
    bool tally_current = false;
    if (curr_segment->_cmfd_surface_fwd != -1 && fwd)
    {
      surf_id = curr_segment->_cmfd_surface_fwd % HEX_NUM_SURFACES;
      cell_id = curr_segment->_cmfd_surface_fwd / HEX_NUM_SURFACES;
      tally_current = true;
    }
    else if (curr_segment->_cmfd_surface_bwd != -1 && !fwd)
    {
      surf_id = curr_segment->_cmfd_surface_bwd % HEX_NUM_SURFACES;
      cell_id = curr_segment->_cmfd_surface_bwd / HEX_NUM_SURFACES;
      tally_current = true;
    }

    /* Tally current if necessary */
    if (tally_current)
    {

      int local_cell_id = getLocalCMFDCell(cell_id);
      if (!_empty_fsrs_cells[local_cell_id])
      {
        // log::finfo("CMFD Cell:%d, FSR:%ld, Fwd surface:(%d, %d)", cell_id, curr_segment->_region_id, curr_segment->_cmfd_surface_fwd, surf_id);
        int logical_cell_id = _logical_actual_map[local_cell_id];

        if (_SOLVE_3D)
        {
          double wgt = _quadrature->getWeightInline(azim_index, polar_index);
          for (int e = 0; e < _num_moc_groups; e++)
          {

            /* Get the CMFD group */
            cmfd_group = getCmfdGroup(e);

            /* Increment the surface group */
            currents[cmfd_group] += track_flux[e] * wgt;
          }

          /* Increment currents */
          if (surf_id < HEX_NUM_FACES)
          {
            _surface_currents->incrementValues(logical_cell_id, surf_id * ncg, (surf_id + 1) * ncg - 1, currents);
          }
          else
          {

            omp_set_lock(&_edge_corner_lock);

            int first_ind = (logical_cell_id * HEX_NUM_SURFACES + surf_id) * ncg;
            it = _edge_corner_currents.find(first_ind);
            if (it == _edge_corner_currents.end())
              for (int g = 0; g < ncg; g++)
                _edge_corner_currents[first_ind + g] = 0.0;

            for (int g = 0; g < ncg; g++)
            {
              _edge_corner_currents[first_ind + g] += currents[g];
            }

            omp_unset_lock(&_edge_corner_lock);
          }
        }
        else
        {
          int pe = 0;
          for (int e = 0; e < _num_moc_groups; e++)
          {

            /* Get the CMFD group */
            cmfd_group = getCmfdGroup(e);

            for (int p = 0; p < _num_polar / 2; p++)
            {
              currents[cmfd_group] += track_flux[pe] * _quadrature->getWeightInline(azim_index, p);
              pe++;
            }
          }

          /* Increment currents */
          if (surf_id < HEX_NUM_FACES)
          {
            _surface_currents->incrementValues(logical_cell_id, surf_id * ncg, (surf_id + 1) * ncg - 1, currents);
          }
          else
          {
            omp_set_lock(&_edge_corner_lock);

            int first_ind = (logical_cell_id * HEX_NUM_SURFACES + surf_id) * ncg;
            it = _edge_corner_currents.find(first_ind);
            if (it == _edge_corner_currents.end())
              for (int g = 0; g < ncg; g++)
                _edge_corner_currents[first_ind + g] = 0.0;

            for (int g = 0; g < ncg; g++)
              _edge_corner_currents[first_ind + g] += currents[g];

            omp_unset_lock(&_edge_corner_lock);
          }
        }
      }
    }
  }

  /**
   * @brief This function tallies the current impinging on the domain from
   *        starting fluxes
   * @details Incoming currents are tallied for use in diagnostics, debugging,
   *          and adjusting sigma-t to enforce consistency with the MOC solution,
   *          if requested
   * @param point The point where the fluxes enter the geometry
   * @param delta_x The a small x-nudge in the direction of travel
   * @param delta_y The a small y-nudge in the direction of travel
   * @param delta_z The a small z-nudge in the direction of travel
   * @param track_flux The angular fluxes impinging on the domain
   * @param weight The weight of the Track
   */
  void Cmfd::tallyStartingCurrent(Point *point, double delta_x, double delta_y,
                                  double delta_z, float *track_flux,
                                  double weight)
  {

    /* Check for non-zero current */
    bool non_zero = false;
    for (int e = 0; e < _num_moc_groups; e++)
    {
      if (fabs(track_flux[e]) > FLT_EPSILON)
      {
        non_zero = true;
        break;
      }
    }
    if (!non_zero)
      return;

    /* Create local coordinate */
    LocalCoords coords;
    coords.setUniverse(_geometry->getRootUniverse());
    coords.setX(point->getX());
    coords.setY(point->getY());
    coords.setZ(point->getZ());

    /* Find the CMFD cell */
    coords.adjustCoords(delta_x, delta_y, delta_z);
    int cell = findCmfdCell(&coords);
    coords.adjustCoords(-delta_x, -delta_y, -delta_z);

    /* Check the CMFD cell */
    if (cell == -1)
      log::ferror("Failed to find starting CMFD cell for track start "
                  "point");
    int cell_x = cell % _local_num_x;
    int cell_y = (cell % (_local_num_x * _local_num_y)) / _local_num_x;
    int cell_z = cell / (_local_num_x * _local_num_y);
    int bounds[3];
    bool singular[3] = {_local_num_x == 1, _local_num_y == 1, _local_num_z == 1};
    bounds[0] = -1 * (cell_x == 0) + (cell_x == _local_num_x - 1);
    bounds[1] = -1 * (cell_y == 0) + (cell_y == _local_num_y - 1);
    bounds[2] = -1 * (cell_z == 0) + (cell_z == _local_num_z - 1);
    if ((bounds[0] == 0 && !singular[0]) && (bounds[1] == 0 && !singular[1]) &&
        (bounds[2] == 0 && !singular[2]))
      log::ferror("Track start point not on a boundary CMFD cell. "
                  "Cell = %d (%d, %d, %d) from Track: (%3.2f, %3.2f, %3.2f) "
                  "adjusted (%3.2e, %3.2e, %3.2e)",
                  cell, cell_x, cell_y, cell_z,
                  point->getX(), point->getY(), point->getZ(), delta_x, delta_y,
                  delta_z);

    int tid = omp_get_thread_num();
    CMFD_PRECISION *currents = _temporary_currents[tid];
    memset(currents, 0.0, sizeof(CMFD_PRECISION) * _num_cmfd_groups);

    /* Tally currents to each CMFD group locally */
    for (int e = 0; e < _num_moc_groups; e++)
    {

      /* Get the CMFD group */
      int cmfd_group = getCmfdGroup(e);

      /* Increment the surface group */
      currents[cmfd_group] += track_flux[e] * weight;
    }

    /* Tally starting currents to cell */
    _starting_currents->incrementValues(cell, 0, _num_cmfd_groups - 1, currents);
  }

  /**
   * @brief Records net currents (leakage) on every CMFD cell for every group
   *  * 对于CMFD网格来说,净电流首先减去该网格所有能群对应的起始电流
   * 再对于每个网格的所有表面(面,边,点),净电流 = 净电流 - 起始电流 + 该网格的表面(边或顶点)电流 - 对应方向的下一个网格的表(边或顶点)电流
   */
  void Cmfd::recordNetCurrents()
  {

#pragma omp parallel for
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++)
    {

      for (int e = 0; e < _num_cmfd_groups; e++)
        _net_currents->incrementValue(i, e,
                                      -1 * _starting_currents->getValue(i, e)); // 净流量首先减去起始流量(因为val=-1 * _starting_currents->getValue(i,e))

      /* Compute cell indexes */
      int cell_ind[3];
      cell_ind[0] = i % _local_num_x;
      cell_ind[1] = (i / _local_num_x) % _local_num_y;
      cell_ind[2] = i / (_local_num_x * _local_num_y);

      if (_hexlattice_enable)
      {
        /* Tally current from all surfaces including edges and corners */
        for (int s = 0; s < HEX_NUM_SURFACES; s++)
        {

          /* Check if edge/corner exists */
          if (s >= HEX_NUM_FACES)
          {
            int idx = i * HEX_NUM_SURFACES * _num_cmfd_groups + s * _num_cmfd_groups;
            std::map<int, CMFD_PRECISION>::iterator it =
                _edge_corner_currents.find(idx);
            if (it == _edge_corner_currents.end())
              continue;
          }

          /* Compute index and vector direction */
          int direction[3];
          convertSurfaceToDirection(s, direction);

          /* Copute the next CMFD cell from the cell indexes and direction */
          int cmfd_cell_next = 0;
          int cell_next_ind[3];
          for (int d = 0; d < 3; d++)
            cell_next_ind[d] = cell_ind[d] + direction[d]; // 下一个CMFD索引为当前索引+当前索引的方向

          cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x + cell_next_ind[2] * (_local_num_x * _local_num_y); // 域一维索引
          if (cell_next_ind[0] < 0 || cell_next_ind[0] >= _local_num_x ||
              cell_next_ind[1] < 0 || cell_next_ind[1] >= _local_num_y ||
              cell_next_ind[2] < 0 || cell_next_ind[2] >= _local_num_z)
            cmfd_cell_next = -1;

          /* Tally net currents */
          if (s < HEX_NUM_FACES)
          {
            int idx = s * _num_cmfd_groups;
            for (int e = 0; e < _num_cmfd_groups; e++)
            {
              double current = 1 * _surface_currents->getValue(i, idx + e);
              _net_currents->incrementValue(i, e, current);
            }

            if (cmfd_cell_next != -1)
            {
              for (int e = 0; e < _num_cmfd_groups; e++)
              {
                double current = -1 * _surface_currents->getValue(i, idx + e);
                _net_currents->incrementValue(cmfd_cell_next, e, current);
              }
            }
          }
          else
          {
            int idx = i * HEX_NUM_SURFACES * _num_cmfd_groups + s * _num_cmfd_groups;
            for (int e = 0; e < _num_cmfd_groups; e++)
            {
              double current = _edge_corner_currents.at(idx + e);
              _net_currents->incrementValue(i, e, current);
            }

            if (cmfd_cell_next != -1)
            {
              for (int e = 0; e < _num_cmfd_groups; e++)
              {
                double current = -1 * _edge_corner_currents.at(idx + e);
                _net_currents->incrementValue(cmfd_cell_next, e, current);
              }
            }
          }
        }
      }
      else
      {
        /* Tally current from all surfaces including edges and corners */
        for (int s = 0; s < NUM_SURFACES; s++)
        {

          /* Check if edge/corner exists */
          if (s >= NUM_FACES)
          {                                                                       // 检查位于边和角的净电流
            int idx = i * NUM_SURFACES * _num_cmfd_groups + s * _num_cmfd_groups; // idx=(表面的全局索引)*CMFD能群数量
            std::map<int, CMFD_PRECISION>::iterator it =
                _edge_corner_currents.find(idx);
            if (it == _edge_corner_currents.end()) // 如果在_edge_corner_currents没有找到idx索引,则直接进入下一次循环
              continue;
          }

          /* Compute index and vector direction */
          int direction[3];
          convertSurfaceToDirection(s, direction);

          /* Copute the next CMFD cell from the cell indexes and direction */
          int cmfd_cell_next = 0;
          int cell_next_ind[3];
          for (int d = 0; d < 3; d++)
            cell_next_ind[d] = cell_ind[d] + direction[d]; // 下一个CMFD索引为当前索引+当前索引的方向

          cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x + cell_next_ind[2] * (_local_num_x * _local_num_y); // 域一维索引
          if (cell_next_ind[0] < 0 || cell_next_ind[0] >= _local_num_x ||
              cell_next_ind[1] < 0 || cell_next_ind[1] >= _local_num_y ||
              cell_next_ind[2] < 0 || cell_next_ind[2] >= _local_num_z)
            cmfd_cell_next = -1;

          /* Tally net currents */
          if (s < NUM_FACES)
          { // 面的统计
            int idx = s * _num_cmfd_groups;
            for (int e = 0; e < _num_cmfd_groups; e++)
            {
              double current = 1 * _surface_currents->getValue(i, idx + e);
              _net_currents->incrementValue(i, e, current); // 净电流+=当前CMFD的表面电流
            }

            if (cmfd_cell_next != -1)
            {
              for (int e = 0; e < _num_cmfd_groups; e++)
              {
                double current = -1 * _surface_currents->getValue(i, idx + e);
                _net_currents->incrementValue(cmfd_cell_next, e, current); // 净电流-=下一个CMFD的表面电流
              }
            }
          }
          else
          {
            int idx = i * NUM_SURFACES * _num_cmfd_groups + s * _num_cmfd_groups;
            for (int e = 0; e < _num_cmfd_groups; e++)
            {
              double current = _edge_corner_currents.at(idx + e);
              _net_currents->incrementValue(i, e, current); // 净电流+=当前CMFD的(边和顶点)电流
            }

            if (cmfd_cell_next != -1)
            {
              for (int e = 0; e < _num_cmfd_groups; e++)
              {
                double current = -1 * _edge_corner_currents.at(idx + e);
                _net_currents->incrementValue(cmfd_cell_next, e, current); // 净电流-=下一个CMFD的(边和顶点)电流
              }
            }
          }
        }
      }
    }
  }

  bool Cmfd::CellinHexLattice(int cell)
  {
    bool inhexlattice = true;

    int x = (cell % (_num_x * _num_y)) % _num_x;
    int y = (cell % (_num_x * _num_y)) / _num_x;
    int z = cell / (_num_x * _num_y);

    int r = _num_r;
    int length = 2 * (r - 1);
    int mid = r - 1;
    int high = (3 * r) - 3;

    if (x < mid || y < mid)
    {
      if ((x + y) < mid)
        inhexlattice = false;
    }
    else if (x > mid || y > mid)
    {
      if ((x + y) > (3 * mid))
        inhexlattice = false;
    }

    return inhexlattice;
  }

  void Cmfd::printHex()
  {

    std::stringstream string;

    for (int i = _local_num_z - 1; i >= 0; i--)
    {
      string << "Z: " << i << std::endl;
      for (int j = _local_num_y - 1; j >= 0; j--)
      {
        for (int k = 0; k < _local_num_x; k++)
        {
          int cell = i * (_local_num_y * _local_num_x) + j * _local_num_x + k;
          if (!_empty_fsrs_cells[cell])
            string << " " << std::setw(4) << std::right << cell << " ";
          else
            string << " " << std::setw(4) << std::setfill(' ') << "-" << " ";
        }
        string << std::endl;
      }
      string << std::endl;
    }

    std::cout << string.str() << std::endl;
    log::finfo(string.str().c_str());
  }

  /* Ensure surfaces are x-y surfaces (no z-crossings) */
  /* Note: this code takes advantage of the numeric representation of
    surfaces to find a mapping that removes z-surfaces   删除z表面的映射*/
  void Cmfd::GetXYSurfaces(int *cmfd_surfaces)
  {
    if (GetHexLatticeEnable())
    {
      00000000000000000000 for (int d = 0; d < 2; d++)
      {
        if (cmfd_surfaces[d] < 0)
          continue;
        int local_surface = cmfd_surfaces[d] % HEX_NUM_SURFACES;

        if (local_surface == 3 || local_surface == 7)
        {
          cmfd_surfaces[d] = -1;
        }
        else if (local_surface > 13)
        {
          int cell = cmfd_surfaces[d] / HEX_NUM_SURFACES;
          int half_surf = local_surface / 2;
          if (local_surface > 25)
          {
            int quart_surf = half_surf / 2;
            local_surface = 2 + quart_surf + (half_surf == 2 * quart_surf);
            cmfd_surfaces[d] = cell * HEX_NUM_SURFACES + local_surface;
          }
          else
          {
            if (local_surface > 13 && local_surface < 18)
            {
              if ((local_surface % 2) == 0)
                local_surface = 0;
              else
                local_surface = 4;
            }
            else if (local_surface > 17 && local_surface < 22)
            {
              if ((local_surface % 2) == 0)
                local_surface = 1;
              else
                local_surface = 5;
            }
            else
            {
              if ((local_surface % 2) == 0)
                local_surface = 2;
              else
                local_surface = 6;
            }
            cmfd_surfaces[d] = cell * HEX_NUM_SURFACES + local_surface;
          }
        }
      }
    }
    else
    {
      for (int d = 0; d < 2; d++)
      {
        int local_surface = cmfd_surfaces[d] % NUM_SURFACES;
        if (local_surface == 2 || local_surface == 5)
        { // 2和5都是这个点只在Z轴的表面，不关联x和y轴
          cmfd_surfaces[d] = -1;
        }
        else if (local_surface > 9)
        {                                             // 大于9说明这个点一定有和Z轴相关的边和角（因为如果只有XY轴，只有4边，4点，2个上述已经排除的）
          int cell = cmfd_surfaces[d] / NUM_SURFACES; // 计算出CMFD单元格索引
          int half_surf = local_surface / 2;
          if (local_surface > 17)
          { // 说明有和Z轴相关的顶角  可能为：18,19,20,21,22,23,24,25（8个）
            int quart_surf = half_surf / 2;
            local_surface = 2 + quart_surf + (half_surf == 2 * quart_surf); // 通过这个计算公式把local_surface固定在【6-9】，把三维顶点转为二维xy顶点
            cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface;
          }
          else
          { // 说明有和Z轴相关的边，可能为：10,11,12,13,14,15,16,17（8个）
            local_surface = (half_surf > 6) + 3 *
                                                  (local_surface != 2 * half_surf);
            cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface; // 通过这个计算公式把local_surface固定在【0-1,3-4】，把涉及z轴的二维边转为一维的面
          } // 通过上式的计算，成功的把cmfd_surfaces的local_surface全部转为【0-1,3-4】，【6-9】区间，也就是映射到XY平面上，没有Z轴
        }
      }
    }
  }

} /* namespace antmoc */
