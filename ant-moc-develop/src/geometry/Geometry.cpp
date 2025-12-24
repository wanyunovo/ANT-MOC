#include "antmoc/Geometry.h"
#include "antmoc/Cell.h"
#include "antmoc/Cmfd.h"
#include "antmoc/Lattice.h"
#include "antmoc/LocalCoords.h"
#include "antmoc/log.h"
#include "antmoc/Material.h"
#include "antmoc/Progress.h"
#include "antmoc/string_utils.h"
#include "antmoc/Surface.h"
#include "antmoc/Track.h"
#include "antmoc/Track3D.h"

#include <sys/stat.h>
#include <sys/types.h>

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <omp.h>
#include <set>
#include <sstream>

namespace antmoc
{

  /**
   * @brief Resets the auto-generated unique IDs for Materials, Surfaces,
   *        Cells and Universes/Lattices to 10000.
   */
  void reset_auto_ids()
  {
    reset_material_id();
    reset_surface_id();
    reset_cell_id();
    reset_universe_id();
  }

  /**
   * @brief Constructor initializes an empty Geometry.
   */
  Geometry::Geometry()
  {

    /* Initialize CMFD object to nullptr */
    _cmfd = nullptr;
    _overlaid_mesh = nullptr;
    _root_universe = nullptr;
    _num_modules_x = 1;
    _num_modules_y = 1;
    _num_modules_z = 1;
    _domain_FSRs_counted = false;
    _contains_FSR_centroids = false;
    _twiddle = false;
    _loaded_from_file = false;
    _energy_group_checked = false;
    _extruded_FSRs_reduced = false;
  }

  /**
   * @brief Destructor clears FSR to Cells and Materials maps.
   */
  Geometry::~Geometry()
  {

    if (_loaded_from_file)
    {
      /* Free all materials */
      for (auto &m : _root_universe->getAllMaterials())
      {
        delete m.second;
      }

      /* Free all surfaces */
      for (auto &s : getAllSurfaces())
      {
        delete s.second;
      }

      /* Free all cells */
      for (auto &c : getAllCells())
      {
        delete c.second;
      }

      /* Free all universes */
      for (auto &u : getAllUniverses())
      {
        delete u.second;
      }
    }

    // Free the material map
    _all_materials.clear();

    /* Free FSR maps if they were initialized */
    if (_FSR_keys_map.size() != 0)
    {

      auto values = _FSR_keys_map.values();
      for (size_t i = 0; i < _FSR_keys_map.size(); i++)
        delete values[i];
      delete[] values;

      // ExtrudedFSR will release memory for its members
      auto extruded_fsrs = _extruded_FSR_keys_map.values();
      for (size_t i = 0; i < _extruded_FSR_keys_map.size(); i++)
        delete extruded_fsrs[i];
      delete[] extruded_fsrs;

      _FSR_keys_map.clear();
      _extruded_FSR_keys_map.clear();
      _FSRs_to_keys.clear();
      _FSRs_to_material_IDs.clear();
    }

    delete _overlaid_mesh;
  }

  /**
   * @brief Returns the total width in the x-direction of the Geometry in cm.
   * @return the total width of the Geometry in the x-direction (cm)
   */
  double Geometry::getWidthX()
  {
    return (getMaxX() - getMinX());
  }

  /**
   * @brief Returns the total width in the y-direction of the Geometry in cm.
   * @return the total width of the Geometry in the y-direction (cm)
   */
  double Geometry::getWidthY()
  {
    return (getMaxY() - getMinY());
  }

  /**
   * @brief Returns the total width in the z-direction of the Geometry in cm.
   * @return the total width of the Geometry in the z-direction (cm)
   */
  double Geometry::getWidthZ()
  {
    return (getMaxZ() - getMinZ());
  }

  /**
   * @brief Return the minimum x-coordinate contained by the Geometry.
   * @return the minimum x-coordinate (cm)
   */
  double Geometry::getMinX()
  {
    if (isDomainDecomposed())
    {
      double geometry_min_x = _root_universe->getMinX();
      double geometry_max_x = _root_universe->getMaxX();
      double domain_width_x = (geometry_max_x - geometry_min_x) / mpi::getNumDomainsX();
      return geometry_min_x + mpi::getDomainIndexX() * domain_width_x;
    }
    else
    {
      return _root_universe->getMinX();
    }
  }

  /**
   * @brief Return the maximum x-coordinate contained by the Geometry.
   * @return the maximum x-coordinate (cm)
   */
  double Geometry::getMaxX()
  {
    if (isDomainDecomposed())
    {
      double geometry_min_x = _root_universe->getMinX();
      double geometry_max_x = _root_universe->getMaxX();
      double domain_width_x = (geometry_max_x - geometry_min_x) / mpi::getNumDomainsX();
      int reverse_index_x = mpi::getNumDomainsX() - mpi::getDomainIndexX() - 1;
      return geometry_max_x - reverse_index_x * domain_width_x;
    }
    else
    {
      return _root_universe->getMaxX();
    }
  }

  /**
   * @brief Return the minimum y-coordinate contained by the Geometry.
   * @return the minimum y-coordinate (cm)
   */
  double Geometry::getMinY()
  {
    if (isDomainDecomposed())
    {
      double geometry_min_y = _root_universe->getMinY();
      double geometry_max_y = _root_universe->getMaxY();
      double domain_width_y = (geometry_max_y - geometry_min_y) / mpi::getNumDomainsY();
      return geometry_min_y + mpi::getDomainIndexY() * domain_width_y;
    }
    else
    {
      return _root_universe->getMinY();
    }
  }

  /**
   * @brief Return the maximum y-coordinate contained by the Geometry.
   * @return the maximum y-coordinate (cm)
   */
  double Geometry::getMaxY()
  {
    if (isDomainDecomposed())
    {
      double geometry_min_y = _root_universe->getMinY();
      double geometry_max_y = _root_universe->getMaxY();
      double domain_width_y = (geometry_max_y - geometry_min_y) / mpi::getNumDomainsY();
      int reverse_index_y = mpi::getNumDomainsY() - mpi::getDomainIndexY() - 1;
      return geometry_max_y - reverse_index_y * domain_width_y;
    }
    else
    {
      return _root_universe->getMaxY();
    }
  }

  /**
   * @brief Return the minimum z-coordinate contained by the Geometry.
   * @return the minimum z-coordinate (cm)
   */
  double Geometry::getMinZ()
  {
    if (isDomainDecomposed())
    {
      double geometry_min_z = _root_universe->getMinZ();
      double geometry_max_z = _root_universe->getMaxZ();
      double domain_width_z = (geometry_max_z - geometry_min_z) / mpi::getNumDomainsZ();
      return geometry_min_z + mpi::getDomainIndexZ() * domain_width_z;
    }
    else
    {
      return _root_universe->getMinZ();
    }
  }

  /**
   * @brief Return the maximum z-coordinate contained by the Geometry.
   * @return the maximum z-coordinate (cm)
   */
  double Geometry::getMaxZ()
  {
    if (isDomainDecomposed())
    {
      double geometry_min_z = _root_universe->getMinZ();
      double geometry_max_z = _root_universe->getMaxZ();
      double domain_width_z = (geometry_max_z - geometry_min_z) / mpi::getNumDomainsZ();
      int reverse_index_z = mpi::getNumDomainsZ() - mpi::getDomainIndexZ() - 1;
      return geometry_max_z - reverse_index_z * domain_width_z;
    }
    else
    {
      return _root_universe->getMaxZ();
    }
  }

  /**
   * @brief Return the global minimum x-coordinate contained by the Geometry.
   * @return the minimum x-coordinate (cm)
   */
  double Geometry::getGlobalMinX()
  {
    return _root_universe->getMinX();
  }

  /**
   * @brief Return the global maximum x-coordinate contained by the Geometry.
   * @return the maximum x-coordinate (cm)
   */
  double Geometry::getGlobalMaxX()
  {
    return _root_universe->getMaxX();
  }

  /**
   * @brief Return the global minimum y-coordinate contained by the Geometry.
   * @return the minimum y-coordinate (cm)
   */
  double Geometry::getGlobalMinY()
  {
    return _root_universe->getMinY();
  }

  /**
   * @brief Return the global maximum y-coordinate contained by the Geometry.
   * @return the maximum y-coordinate (cm)
   */
  double Geometry::getGlobalMaxY()
  {
    return _root_universe->getMaxY();
  }

  /**
   * @brief Return the global minimum z-coordinate contained by the Geometry.
   * @return the minimum z-coordinate (cm)
   */
  double Geometry::getGlobalMinZ()
  {
    return _root_universe->getMinZ();
  }

  /**
   * @brief Return the global maximum z-coordinate contained by the Geometry.
   * @return the maximum z-coordinate (cm)
   */
  double Geometry::getGlobalMaxZ()
  {
    return _root_universe->getMaxZ();
  }

  /**
   * @brief Returns the boundary conditions at the minimum x-coordinate in the
   *        Geometry / domain if the geometry is domain-decomposed.
   * @return the boundary conditions for the minimum x-coordinate in the domain
   */
  boundaryType Geometry::getMinXBoundaryType()
  {
    if (isDomainDecomposed() && mpi::getDomainIndexX() > 0)
      return INTERFACE;
    else // 几何结构没有进行域分解或者当前进程在 x 方向上是第一个域
      return _root_universe->getMinXBoundaryType();
  }

  /**
   * @brief Returns the boundary conditions at the maximum x-coordinate in the
   *        Geometry / domain if the geometry is domain-decomposed.
   * @return the boundary conditions for the maximum x-coordinate in the domain
   */
  boundaryType Geometry::getMaxXBoundaryType()
  {
    if (isDomainDecomposed() && mpi::getDomainIndexX() < mpi::getNumDomainsX() - 1)
      return INTERFACE;
    else
      return _root_universe->getMaxXBoundaryType();
  }

  /**
   * @brief Returns the boundary conditions at the minimum y-coordinate in the
   *        Geometry / domain if the geometry is domain-decomposed.
   * @return the boundary conditions for the minimum y-coordinate in the domain
   */
  boundaryType Geometry::getMinYBoundaryType()
  {
    if (isDomainDecomposed() && mpi::getDomainIndexY() > 0)
      return INTERFACE;
    else
      return _root_universe->getMinYBoundaryType();
  }

  /**
   * @brief Returns the boundary conditions at the maximum y-coordinate in the
   *        Geometry / domain if the geometry is domain-decomposed.
   * @return the boundary conditions for the maximum y-coordinate in the domain
   */
  boundaryType Geometry::getMaxYBoundaryType()
  {
    if (isDomainDecomposed() && mpi::getDomainIndexY() < mpi::getNumDomainsY() - 1)
      return INTERFACE;
    else
      return _root_universe->getMaxYBoundaryType();
  }

  /**
   * @brief Returns the boundary conditions at the minimum z-coordinate in the
   *        Geometry / domain if the geometry is domain-decomposed.
   * @return the boundary conditions for the minimum z-coordinate in the domain
   */
  boundaryType Geometry::getMinZBoundaryType()
  {
    if (isDomainDecomposed() && mpi::getDomainIndexZ() > 0)
      return INTERFACE;
    else
      return _root_universe->getMinZBoundaryType();
  }

  /**
   * @brief Returns the boundary conditions at the maximum z-coordinate in the
   *        Geometry / domain if the geometry is domain-decomposed.
   * @return the boundary conditions for the maximum z-coordinate in the domain
   */
  boundaryType Geometry::getMaxZBoundaryType()
  {
    if (isDomainDecomposed() && mpi::getDomainIndexZ() < mpi::getNumDomainsZ() - 1)
      return INTERFACE;
    else
      return _root_universe->getMaxZBoundaryType();
  }

  /**
   * @brief Returns the number of source regions in the Geometry domain.
   * @return number of FSRs
   */
  long Geometry::getNumFSRs()
  {
    return _FSR_keys_map.size();
  }

  /**
   * @brief Returns the number of extruded FSRs in the Geometry domain.
   * @return number of extruded FSRs
   */
  long Geometry::getNumExtrudedFSRs()
  {
    return _extruded_FSR_keys_map.size();
  }

  /// \brief Returns the number of source regions in the entire Geometry.
  /// \return The total number of FSRs.
  long Geometry::getNumTotalFSRs()
  {
    long domain_fsrs = _FSRs_to_keys.size();
    long total_fsrs = domain_fsrs;
#ifdef ENABLE_MPI_
    if (isDomainDecomposed())
      MPI_Allreduce(&domain_fsrs, &total_fsrs, 1, mpi::getDatatype<long>(), MPI_SUM, getMPICart());
#endif
    return total_fsrs;
  }

  /// \brief Returns the number of energy groups for each Material's nuclear data.
  /// \details This method suffers from performance penalty due to getAllMaterials
  ///          because all the materials will be gathered by iterating the whole
  ///          geometry. A solution provided by OpenMOC is to define a macro named
  ///          NGROUPS to fix the number of energy group. Another solution is to
  ///          return the cached map everytime this method is invoked. The second
  ///          solution must be treated carefully to make sure that materials won't
  ///          be changed at runtime.
  /// \return the number of energy groups
  int Geometry::getNumEnergyGroups()
  {

    auto &materials = getAllMaterials();

    if (materials.size() == 0)
      log::ferror("Unable to return the number of energy groups from "
                  "the Geometry since it does not contain any Materials");

    int num_groups = materials.begin()->second->getNumEnergyGroups();

    // Check energy group consistent only once. If there is a need to re-check
    // energy groups, the flag must be set to false.
    if (!_energy_group_checked)
    {
      for (const auto &m : materials)
      {
        if (m.second->getNumEnergyGroups() != num_groups)
          log::ferror("Unable to return the number of energy groups from "
                      "the Geometry since it contains different numbers of groups: "
                      "%d and %d",
                      num_groups, m.second->getNumEnergyGroups());
      }

      _energy_group_checked = true;
    }

    return num_groups;
  }

  /**
   * @brief Returns the number of Materials in the Geometry.
   * @return the number of Materials
   */
  int Geometry::getNumMaterials()
  {

    auto &all_materials = getAllMaterials();
    int num_materials = all_materials.size();

    return num_materials;
  }

  /**
   * @brief Returns the number of Cells in the Geometry.
   * @return the number of Cells
   */
  int Geometry::getNumCells()
  {

    int num_cells = 0;

    if (_root_universe != nullptr)
    {
      std::map<int, Cell *> all_cells = _root_universe->getAllCells();
      num_cells = all_cells.size();
    }

    return num_cells;
  }

  /**
   * @brief Return a std::map container of Surface IDs (keys) with Surfaces
   *        pointers (values).
   * @return a std::map of Surfaces indexed by Surface ID in the geometry
   */
  std::map<int, Surface *> Geometry::getAllSurfaces()
  {

    std::map<int, Surface *> all_surfaces;

    if (_root_universe != nullptr)
    {
      std::map<int, Cell *> all_cells = getAllCells();

      for (const auto &c : all_cells)
      {
        Cell *cell = c.second;
        auto surfaces = cell->getSurfaces();

        for (const auto &s : surfaces)
        {
          Surface *surf = s.second;
          all_surfaces[surf->getId()] = surf;
        }
      }
    }

    return all_surfaces;
  }

  /// \brief Retrieve all materials from the root universe to this geometry.
  ///        This method is supposed to be invoked before any further usage
  ///        of getAllMaterials().
  ///        拿到所有cell里用到的材料，存放在_all_materials里，键是材料id，值是材料指针。 该方法应在任何需要调用getAllMaterials()之前执行。
  void Geometry::retrieveAllMaterials()
  {
/*
宏定义：这是 OpenMP 的编译器指令（pragma directive）
作用：只让主线程执行紧随其后的代码块，其他线程跳过
为什么需要：在并行环境中，收集材料信息只需一个线程做一次即可，避免重复工作和数据竞争
*/
#pragma omp master
    {
      /* Gather all materials from all cells */
      if (_root_universe != nullptr)
      {
        std::map<int, Cell *> all_cells = _root_universe->getAllCells();

        for (const auto &c : all_cells)
        {
          if (c.second->getType() == MATERIAL) // 只处理材料类型的 Cell
          {
            Material *material = c.second->getFillMaterial();

            if (material != nullptr)
              _all_materials[material->getId()] = material;
          }
        }
      }
    }

/*
宏定义：OpenMP 的同步障碍指令
作用：所有线程必须等待，直到主线程完成材料收集后，所有线程都到达这个点，才能继续执行，保证数据一致性
*/
#pragma omp barrier
  }

  /// \brief Return a std::map container of Material IDs (keys) with Materials
  ///        pointers (values). This method should only be invoked after
  ///        retrieveAllMaterials().
  std::map<int, Material *> &Geometry::getAllMaterials()
  {
    if (_all_materials.empty())
      log::fwarn("Empty _all_materials is returned, which may be buggy");

    return _all_materials;
  }

  /**
   * @brief Modify scattering and total cross sections to study MOC stability
   */
  void Geometry::manipulateXS()
  {

    auto &all_materials = getAllMaterials();

    for (auto m : all_materials)
    {

      Material *mat = m.second;
      int ng = mat->getNumEnergyGroups();
      for (int g = 0; g < ng; g++)
      {
        double sigma_t = mat->getSigmaTByGroup(g + 1);
        for (int gp = 0; gp < ng; gp++)
          if (mat->getSigmaSByGroup(g + 1, gp + 1) < 0)
            mat->setSigmaSByGroup(0.0, g + 1, gp + 1);
        double sigma_s = 0.0;
        for (int gp = 0; gp < ng; gp++)
          sigma_s += mat->getSigmaSByGroup(g + 1, gp + 1);
        sigma_s *= 1.1;
        if (sigma_s > sigma_t)
          mat->setSigmaTByGroup(sigma_s, g + 1);
      }
    }
  }

  /**
   * @brief Return a std::map container of Cell IDs (keys) with Cells
   *        pointers (values).
   * @return a std::map of Cells indexed by Cell ID in the geometry
   */
  std::map<int, Cell *> Geometry::getAllCells()
  {

    std::map<int, Cell *> all_cells;

    if (_root_universe != nullptr)
      all_cells = _root_universe->getAllCells();

    return all_cells;
  }

  /**
   * @brief Return a std::map container of Cell IDs (keys) with Cells
   *        pointers (values).
   * @return a std::map of Cells indexed by Cell ID in the geometry
   */
  std::map<int, Cell *> Geometry::getAllMaterialCells()
  {

    std::map<int, Cell *> all_material_cells;

    if (_root_universe != nullptr)
    {
      std::map<int, Cell *> all_cells = _root_universe->getAllCells();

      for (const auto &c : all_cells)
      {
        Cell *cell = c.second;

        if (cell->getType() == MATERIAL)
          all_material_cells[cell->getId()] = cell;
      }
    }

    return all_material_cells;
  }

  /**
   * @brief Return a std::map container of Universe IDs (keys) with Unierses
   *        pointers (values).
   * @return a std::map of Universes indexed by Universe ID in the geometry
   */
  std::map<int, Universe *> Geometry::getAllUniverses()
  {

    std::map<int, Universe *> all_universes;

    if (_root_universe != nullptr)
      all_universes = _root_universe->getAllUniverses();

    return all_universes;
  }

  /**
   * @brief Returns the Universe at the root node in the CSG tree.
   * @return the root Universe
   */
  Universe *Geometry::getRootUniverse()
  {
    return _root_universe;
  }

  /**
   * @brief Returns a pointer to the CMFD object.
   * @return A pointer to the CMFD object
   */
  Cmfd *Geometry::getCmfd()
  {
    return _cmfd;
  }

  /// \brief See mpi::isSpatialDecomposed
  bool Geometry::isDomainDecomposed()
  {
    return mpi::isSpatialDecomposed();
  }

  /// \brief See mpi::isRootDomain
  bool Geometry::isRootDomain()
  {
    return mpi::isRootDomain();
  }

  /// \brief See mpi::getDomainIndexes
  void Geometry::getDomainIndexes(int *indexes)
  {
    mpi::getDomainIndexes(indexes);
  }

  /// \brief See mpi::getDomainStructure
  void Geometry::getDomainStructure(int *structure)
  {
    mpi::getDomainStructure(structure);
  }

  /**
   * @brief Sets the root Universe for the CSG tree.
   * @param root_universe the root Universe of the CSG tree.
   */
  void Geometry::setRootUniverse(Universe *root_universe)
  {
    _root_universe = root_universe;
  }

  // FIXME: add setDefaultDomainDecomposition() function

  /**
   * @brief Sets how many modular track laydown domains are in each MPI domain
   * @param num_x The number of modular domains in the x-direction per MPI domain
   * @param num_y The number of modular domains in the y-direction per MPI domain
   * @param num_z The number of modular domains in the z-direction per MPI domain
   */
  void Geometry::setNumDomainModules(int num_x, int num_y, int num_z)
  {
    _num_modules_x = num_x;
    _num_modules_y = num_y;
    _num_modules_z = num_z;
  }

  /**
   * @brief Get the number of modular domains in the x-direction per MPI domain
   * @return _num_modules_x number of modular domains in the x-direction in domain
   */
  int Geometry::getNumXModules()
  {
    return _num_modules_x;
  }

  /**
   * @brief Get the number of modular domains in the y-direction per MPI domain
   * @return _num_modules_y number of modular domains in the y-direction in domain
   */
  int Geometry::getNumYModules()
  {
    return _num_modules_y;
  }

  /**
   * @brief Get the number of modular domains in the z-direction per MPI domain
   * @return _num_modules_z number of modular domains in the z-direction in domain
   */
  int Geometry::getNumZModules()
  {
    return _num_modules_z;
  }

  /**
   * @brief Sets the pointer to a CMFD object used for acceleration.
   * @param cmfd a pointer to the CMFD object
   */
  void Geometry::setCmfd(Cmfd *cmfd)
  {
    _cmfd = cmfd;
  }

  /**
   * @brief Sets a global overlaid mesh with the given mesh height
   * @details The global overlaid mesh is overlaid across the entire Geometry
   * @param axial_mesh_height The desired height of axial mesh cells
   * @param num_x number of divisions in the X direction
   * @param num_y number of divisions in the Y direction
   * @param num_radial_domains number of radial domains
   * @param radial_domains array with the indexes of each domain in X and Y
   */
  void Geometry::setOverlaidMesh(double axial_mesh_height, int num_x, int num_y,
                                 int num_radial_domains, int *radial_domains)
  {

    /* Get the global Geometry boundaries */
    double min_x = _root_universe->getMinX();
    double max_x = _root_universe->getMaxX();
    double min_y = _root_universe->getMinY();
    double max_y = _root_universe->getMaxY();
    double min_z = _root_universe->getMinZ();
    double max_z = _root_universe->getMaxZ();

    /* Create the lattice */
    _overlaid_mesh = new RecLattice();
    int real_num_x = 1;
    int real_num_y = 1;
    if (num_x > 0 && num_y > 0)
    {
      if (radial_domains == nullptr)
      {
        real_num_x = num_x;
        real_num_y = num_y;
      }
      else
      {
        for (int i = 0; i < num_radial_domains; i++)
        {
          if (radial_domains[2 * i] == mpi::getDomainIndexX() &&
              radial_domains[2 * i + 1] == mpi::getDomainIndexY())
          {
            real_num_x = num_x;
            real_num_y = num_y;
          }
        }
      }
    }
    num_x = real_num_x;
    num_y = real_num_y;
    _overlaid_mesh->setNumX(num_x);
    _overlaid_mesh->setNumY(num_y);

    /* Determine actual axial mesh spacing from desired spacing */
    double total_width_z = max_z - min_z;
    int num_cells_z = total_width_z / axial_mesh_height;
    axial_mesh_height = total_width_z / num_cells_z;
    double mesh_width_x = (max_x - min_x) / num_x;
    double mesh_width_y = (max_y - min_y) / num_y;
    _overlaid_mesh->setNumZ(num_cells_z);
    _overlaid_mesh->setWidth(mesh_width_x, mesh_width_y, axial_mesh_height);

    /* Create the center point */
    Point offset;
    offset.setX(min_x + (max_x - min_x) / 2.0);
    offset.setY(min_y + (max_y - min_y) / 2.0);
    offset.setZ(min_z + (max_z - min_z) / 2.0);
    _overlaid_mesh->setOffset(offset.getX(), offset.getY(), offset.getZ());

    _overlaid_mesh->computeSizes();

    log::finfo("Set global axial mesh of width %6.4f cm",
               axial_mesh_height);
  }

  /**
   * @brief Clears all boundary conditions from the Geometry
   */
  void Geometry::clearBoundaries()
  {

    /* Extract all surfaces from the Geometry */
    std::map<int, Surface *> all_surfaces = getAllSurfaces();
    std::map<int, Surface *>::iterator it;

    /* Iterate over all surfaces */
    for (it = all_surfaces.begin(); it != all_surfaces.end(); ++it)
    {

      /* Remove boundary conditions */
      Surface *surface = it->second;
      if (surface->getBoundaryType() != BOUNDARY_NONE)
        surface->setBoundaryType(BOUNDARY_NONE);
    }
  }

  /**
   * @brief Find the Cell that this LocalCoords object is in at the lowest level
   *        of the nested Universe hierarchy.
   * @details This method assumes that the LocalCoords has been initialized
   *          with coordinates and a Universe ID. The method will recursively
   *          find the Cell on the lowest level of the nested Universe hierarchy
   *          by building a linked list of LocalCoords from the LocalCoord
   *          passed in as an argument down to the lowest level Cell found. In
   *          the process it will set the coordinates at each level of the
   *          hierarchy for each LocalCoord in the linked list for the Lattice
   *          or Universe that it is in. If the LocalCoords is outside the bounds
   *          of the Geometry or on the boundaries this method will return nullptr;
   *          otherwise it will return a pointer to the Cell that is found by the
   *          recursive Geometry::findCell(...) method.
   * @param coords pointer to a LocalCoords object
   * @return returns a pointer to a Cell if found, nullptr if no Cell found
   */
  Cell *Geometry::findCellContainingCoords(LocalCoords *coords)
  {

    Universe *univ = coords->getUniverse();
    Cell *cell;

    /* Check if the coords are inside the geometry bounds */
    if (univ->getId() == _root_universe->getId())
    {
      if (!withinBounds(coords))
        return nullptr;
    }

    if (univ->getType() == SIMPLE)
      cell = univ->findCell(coords);
    else
      cell = static_cast<Lattice *>(univ)->findCell(coords);

    return cell;
  }

  /**
   * @brief Find the first Cell of a Track segment with a starting Point that is
   *        represented by the LocalCoords method parameter.
   * @details This method assumes that the LocalCoords has been initialized
   *          with coordinates and a Universe ID. This method will move the
   *          initial starting point by a small amount along the direction of
   *          the Track in order to ensure that the track starts inside of a
   *          distinct FSR rather than on the boundary between two of them.
   *          The method will recursively find the LocalCoords by building a
   *          linked list of LocalCoords from the LocalCoords passed in as an
   *          argument down to the Cell found in the lowest level of the nested
   *          Universe hierarchy. In the process, the method will set the
   *          coordinates at each level in the nested Universe hierarchy for
   *          each LocalCoord in the linked list for the Lattice or Universe
   *          that it is in.
   * @param coords pointer to a LocalCoords object
   * @param azim the azimuthal angle for a trajectory projected from the LocalCoords
   * @param polar the polar angle for a trajectory projected from the LocalCoords
   * @return returns a pointer to a cell if found, nullptr if no cell found
   */
  Cell *Geometry::findFirstCell(LocalCoords *coords, double azim, double polar)
  {
    double delta_x = cos(azim) * sin(polar) * TINY_MOVE;
    double delta_y = sin(azim) * sin(polar) * TINY_MOVE;
    double delta_z = cos(polar) * TINY_MOVE;
    coords->adjustCoords(delta_x, delta_y, delta_z);
    return findCellContainingCoords(coords);
  }

  /**
   * @brief Find the Material for a flat source region ID.
   * @details  This method finds the fsr_id within the
   *           _FSR_to_material_IDs map and returns the corresponding
   *           pointer to the Material object.
   * @param fsr_id a FSR id
   * @return a pointer to the Material that this FSR is in
   */
  Material *Geometry::findFSRMaterial(long fsr_id)
  {

    auto &all_materials = getAllMaterials();

    int mat_id = _FSRs_to_material_IDs.at(fsr_id);
    if (all_materials.find(mat_id) == all_materials.end())
      log::ferror("Failed to find FSR Material for FSR %ld with "
                  "Material ID %d",
                  fsr_id, mat_id);

    return all_materials[mat_id];
  }

  /**
   * @brief Finds the next Cell for a LocalCoords object along a trajectory
   *        defined by two angles : azimuthal from 0 to 2Pi, polar from 0 to Pi.
   * @details The method will update the LocalCoords passed in as an argument
   *          to be the one at the boundary of the next Cell crossed along the
   *          given trajectory. It will do this by finding the minimum distance
   *          to the surfaces at all levels of the coords hierarchy.
   *          If the LocalCoords is outside the bounds of the Geometry or on
   *          the boundaries this method will return nullptr; otherwise it will
   *          return a pointer to the Cell that the LocalCoords will reach
   *          next along its trajectory.
   * 该方法将更新作为参数传入的LocalCoords，使其成为沿给定轨迹交叉的下一个Cell的边界处的LocalCoords。
   * 它将通过沿给定方向找到轨迹与最近面的交点来实现这一点。如果LocalCoords在Geometry的边界之外或在边界上，
   * 则此方法将返回nullptr；否则，它将返回一个指向LocalCoords将沿着其轨迹下一个到达的Cell的指针。
   * @param coords pointer to a LocalCoords object  传入时是起点，函数返回时它已经是终点
   * @param azim the azimuthal angle of the trajectory
   * @param polar the polar angle of the trajectory
   * @return a pointer to a Cell if found, nullptr if no Cell found  线段终点所在的Cell
   */
  Cell *Geometry::findNextCell(LocalCoords *coords, double azim, double polar)
  {

    double dist;
    double min_dist = std::numeric_limits<double>::infinity(); // 记录这条射线到各级边界的最小距离。初值 ∞

    /* Save coords and angles in case of a translated/rotated cell  保存原来的坐标与角度，防止在多级 Universe 中做平移/旋转时破坏上层信息，后面会恢复。*/
    double old_position[3] = {coords->getX(), coords->getY(), coords->getZ()};
    double old_azim = azim;
    double old_polar = polar;

    /* Get highest level coords */
    coords = coords->getHighestLevel();

    /* If the current coords is outside the domain, return nullptr */
    /*
    如果几何启用了区域分解 (isDomainDecomposed())，先看当前点是否还在本域边界 _domain_bounds 里；若不在直接返回 nullptr。若 coords->getLowestLevel()->getCell() == NULL，说明还没定位到任何 Cell（可能在空白区域），返回 nullptr。
     */
    if (isDomainDecomposed())
    {
      Point *point = coords->getPoint();
      if (!_domain_bounds->containsPoint(point))
        return nullptr;
    }

    /* If the current coords is not in any Cell, return nullptr */
    if (coords->getLowestLevel()->getCell() == NULL)
      return nullptr;

    /* If the current coords is inside a Cell, look for next Cell */
    else
    {

      /* Descend universes/coord until at the lowest level.
       * At each universe/lattice level get distance to next
       * universe or lattice cell. Compare to get new min_dist. 遍历这个点所在链中的每一级CSG，比较该点与面的最小距离，找到该点所有级别中到最近面的最小dist*/
      while (coords != nullptr)
      {

        /* If we reach a LocalCoord in a Lattice, find the distance to the
         * nearest lattice cell boundary */
        if (coords->getType() == LAT)
        { // Lattice 边界计算最小距离
          auto lattice = coords->getLattice();
          dist = lattice->minSurfaceDist(coords->getPoint(), azim, polar);
        }
        /* If we reach a LocalCoord in a Universe, find the distance to the
         * nearest cell surface */
        else
        { // cell 边界计算最小距离
          Cell *cell = coords->getCell();
          dist = cell->minSurfaceDist(coords->getPoint(), azim, polar);

          // 如果 Cell 有平移或旋转属性，需要把 coords 和方向向量 (sinθ·cosφ, sinθ·sinφ, cosθ) 变换到 Cell 的局部坐标系再做距离计算，这就有那段平移/旋转矩阵的代码。
          /* Apply translation to position */
          if (cell->isTranslated())
          {
            double *translation = cell->getTranslation();
            double new_x = coords->getX() - translation[0];
            double new_y = coords->getY() - translation[1];
            double new_z = coords->getZ() - translation[2];
            coords->setX(new_x);
            coords->setY(new_y);
            coords->setZ(new_z);
          }

          /* Apply rotation to position and direction */
          if (cell->isRotated())
          {
            double x = coords->getX();
            double y = coords->getY();
            double z = coords->getZ();
            double *matrix = cell->getRotationMatrix();
            double new_x = matrix[0] * x + matrix[1] * y + matrix[2] * z;
            double new_y = matrix[3] * x + matrix[4] * y + matrix[5] * z;
            double new_z = matrix[6] * x + matrix[7] * y + matrix[8] * z;
            coords->setX(new_x);
            coords->setY(new_y);
            coords->setZ(new_z);

            double uvw[3] = {sin(polar) * cos(azim), sin(polar) * sin(azim),
                             cos(polar)};
            double rot_uvw[3] = {matrix[0] * uvw[0] + matrix[1] * uvw[1] +
                                     matrix[2] * uvw[2],
                                 matrix[3] * uvw[0] +
                                     matrix[4] * uvw[1] + matrix[5] * uvw[2],
                                 matrix[6] * uvw[0] + matrix[7] * uvw[1] +
                                     matrix[8] * uvw[2]};
            polar = acos(rot_uvw[2]);
            double sin_p = sqrt(1 - rot_uvw[2] * rot_uvw[2]);
            int sgn = (asin(rot_uvw[1] / sin_p) > 0) - (asin(rot_uvw[1] / sin_p) < 0);
            azim = acos(rot_uvw[0] / sin_p) * sgn;
          }
        }

        /* Recheck min distance */
        min_dist = std::min(dist, min_dist); // 始终记录所有层中最小的那一个。

        /* Descend one level, exit if at the lowest level of coordinates */
        if (coords->getNext() == nullptr)
          break;
        else
          coords = coords->getNext();
      }

      /* Reset coords position in case there was a translated or rotated cell */
      // 在 Cell 有平移/旋转的过程中我们改动了 coords；现在通过 prune() 清掉临时子节点，还原到原来的位置和角度。
      coords = coords->getHighestLevel();
      coords->prune();
      coords->setX(old_position[0]);
      coords->setY(old_position[1]);
      coords->setZ(old_position[2]);
      azim = old_azim;
      polar = old_polar;

      /* Check for distance to an overlaid mesh 检查到（_overlaid_mesh）覆盖网格的最小距离*/
      if (_overlaid_mesh != nullptr)
      {
        dist = _overlaid_mesh->minSurfaceDist(coords->getPoint(), azim, polar);
        min_dist = std::min(dist, min_dist);
      }

      /* Check for distance to nearest CMFD mesh cell boundary  检查到CMFD网格的最小距离*/
      if (_cmfd != nullptr)
      {
        auto lattice = _cmfd->getLattice();
        dist = lattice->minSurfaceDist(coords->getPoint(), azim, polar); // 找到在CMFD Lattice 中距离下一下X，Y ，Z面中的最小距离
        min_dist = std::min(dist, min_dist);
      }

      /* Check for distance to nearest domain boundary */
      bool domain_boundary = false;
      if (isDomainDecomposed())
      {
        dist = _domain_bounds->minSurfaceDist(coords->getPoint(), azim, polar);
        if (dist - min_dist < ON_SURFACE_THRESH)
        {
          min_dist = dist;
          domain_boundary = true;
        }
      }

      /* Move point and get next cell 将该点沿给定方向移动最小距离*/
      double delta_x = cos(azim) * sin(polar) * (min_dist + TINY_MOVE);
      double delta_y = sin(azim) * sin(polar) * (min_dist + TINY_MOVE);
      double delta_z = cos(polar) * (min_dist + TINY_MOVE);
      coords->adjustCoords(delta_x, delta_y, delta_z);

      if (domain_boundary)
        return nullptr;
      else
        return findCellContainingCoords(coords); // 返回这个点所在的cell
    }
  }

  /**
   * @brief Get the symmetries used to restrict the domain
   * @return a boolean indicating if the symmetry along this axis is used
   */
  bool Geometry::getSymmetry(int axis)
  {
    return _symmetries[axis];
  }

  /**
   * @brief Find and return the ID of the flat source region that a given
   *        LocalCoords object resides within.
   * @param coords a LocalCoords object pointer
   * @return the FSR ID for a given LocalCoords object
   * findFSRId 输入一个 LocalCoords*，返回这个坐标所属的 FSR（Flat Source Region）编号。
   * _FSR_keys_map建立fsr_id与FSRData* 的映射
   *
   * 1.定位并生成“唯一标签”：沿着 coords 的几何层级（CMFD→外加网格→lattice/universe→cell）拼出一串字符串，这串 key 描述了“这个点是通过哪几层几何定位到的”。同一 FSR 的点一定会得到同一 key。
   * _FSR_keys_map 是线程安全的 字符串-->FSRData* 映射。若 map 中还没有这条 key，就调用 insert_and_get_count 分配一个新的 fsr_id。这里存在多线程竞争，所以：
   * 如果返回 -1，说明别的线程抢先插入了同一 key，本线程就自旋等待 map 里对应条目不为空后读取它；
   * 否则当前线程负责填充新的 FSRData：创建记录点坐标的 Point，调用 findCellContainingCoords(curr) 拿到所属 Cell，再把该 Cell 的 Material ID 记录下来；若启用了 CMFD，还会查它的 CMFD cell ID。最后把 FSRData 填入 map，供其他线程复用。
   *
   */
  long Geometry::findFSRId(LocalCoords *coords)
  {

    long fsr_id;
    LocalCoords *curr = coords;
    curr = coords->getLowestLevel();

    /* Generate unique FSR key */
    int thread_id = omp_get_thread_num();
    std::string &fsr_key = _fsr_keys[thread_id];

    /*getFSRKeyFast会先统计字符串的长度再 resize 字符串。由于要输出的信息既有字符又有数字，
    这个函数使用了 replace 和自定义的 printToString 来处理。getFSRKey 则用流ostringstream
    来处理*/
    getFSRKeyFast(coords, fsr_key);

    /* If FSR has not been encountered, update FSR maps and vectors */
    if (!_FSR_keys_map.contains(fsr_key))
    { // 若map里面还没有这个标识串

      /* Try to get a clean copy of the fsr_id, adding the FSR data
         if necessary where -1 indicates the key was already added 插入串并返回id*/
      fsr_id = _FSR_keys_map.insert_and_get_count(fsr_key, nullptr);

      if (fsr_id == -1)
      { // 插入失败，表明另一个线程已经插入了fsr_key，等待
        FSRData volatile *fsr;
        do
        {
          fsr = _FSR_keys_map.at(fsr_key);
        } while (fsr == nullptr);
        fsr_id = fsr->_fsr_id;
      }
      else
      {

        /* Add FSR information to FSR key map and FSR_to vectors */
        auto fsr = new FSRData(); // 新的data实例
        fsr->_fsr_id = fsr_id;
        _FSR_keys_map.update(fsr_key, fsr); // 插入到map中

        // 设置fsr的特征点
        Point *point = new Point();
        point->setCoords(coords->getHighestLevel()->getX(),
                         coords->getHighestLevel()->getY(),
                         coords->getHighestLevel()->getZ());

        /* Get the cell that contains coords */
        /* 找到fsr对应的cell，获取其材料的id */
        Cell *cell = findCellContainingCoords(curr);
        fsr->_point = point;
        fsr->_mat_id = cell->getFillMaterial()->getId();

        /* If CMFD acceleration is on, add FSR CMFD cell to FSR data 如果启用了CMFD加速，则将FSR的CMFD单元添加到FSR数据中*/
        if (_cmfd != nullptr)
        {
          if (_cmfd->findCmfdCell(coords->getHighestLevel()) != -1)
          {
            fsr->_cmfd_cell = _cmfd->findCmfdCell(coords->getHighestLevel());
          }
        }
      }
    }

    /* If FSR has already been encountered, get the fsr id from map */
    else
    { // 若map里面已经有这个标识串，获取它
      FSRData volatile *fsr;
      do
      {
        fsr = _FSR_keys_map.at(fsr_key);
      } while (fsr == nullptr);

      fsr_id = fsr->_fsr_id;
    }

    return fsr_id;
  }

  /**
   * @brief Finds and returns a pointer to the axially extruded flat source
   *        region that a given LocalCoords object resides within.
   * @param coords a LocalCoords object pointer
   * @return the ID of an extruded FSR for a given LocalCoords object
   * 根据给定的局部坐标点返回相应的 FSR id。执行函数的过程也就是给 FSR 编号的过程。
   */
  int Geometry::findExtrudedFSR(LocalCoords *coords)
  {

    int fsr_id;
    // unused
    // LocalCoords* curr = coords;
    // curr = coords->getLowestLevel();

    /* Generate unique FSR key */
    std::string fsr_key = getFSRKey(coords);

    /* If FSR has not been encountered, update FSR maps and vectors */
    if (!_extruded_FSR_keys_map.contains(fsr_key))
    {

      /* Try to get a clean copy of the fsr_id, adding the FSR data
         if necessary where -1 indicates the key was already added */
      fsr_id = _extruded_FSR_keys_map.insert_and_get_count(fsr_key, nullptr);

      if (fsr_id == -1)
      {
        ExtrudedFSR volatile *fsr;
        int count = 0;
        do
        {
          fsr = _extruded_FSR_keys_map.at(fsr_key);
          count++;
          if (count > 1e8)
            log::ferror("Application stuck in finding extruded FSR");
        } while (fsr == nullptr);
        fsr_id = fsr->_fsr_id;
      }
      else
      {

        /* Add FSR information to FSR key map and FSR_to vectors */
        auto fsr = new ExtrudedFSR();
        fsr->_fsr_id = fsr_id;
        fsr->_num_fsrs = 0;
        fsr->_coords = new LocalCoords(0, 0, 0, true);
        coords->copyCoords(fsr->_coords);
        _extruded_FSR_keys_map.update(fsr_key, fsr);
      }
    }

    /* If FSR has already been encountered, get the fsr id from map */
    else
    {
      ExtrudedFSR volatile *fsr;
      int count = 0;
      do
      {
        fsr = _extruded_FSR_keys_map.at(fsr_key);
        count++;
        if (count > 1e8)
          log::ferror("Application stuck in finding extruded FSR");
      } while (fsr == nullptr);

      fsr_id = fsr->_fsr_id;
    }

    /*
    log::finfo("x:%.2f, y:%.2f, z:%.2f", coords->getX(), coords->getY(), coords->getZ());
    log::finfo("FSR Key :%s", fsr_key);
    log::finfo("FSR id: %d", fsr_id);
    */

    return fsr_id;
  }

  /**
   * @brief Return the ID of the flat source region that a given
   *        LocalCoords object resides within.
   * @param coords a LocalCoords object pointer
   * @param err_check whether to fail instead of returning -1 if not found
   * @return the FSR ID for a given LocalCoords object
   */
  long Geometry::getFSRId(LocalCoords *coords, bool err_check)
  {

    long fsr_id = 0;
    int thread_id = omp_get_thread_num();

    try
    {
      /* Generate unique FSR key */
      std::string &fsr_key = _fsr_keys[thread_id];
      getFSRKeyFast(coords, fsr_key);
      if (!_FSR_keys_map.contains(fsr_key) && !err_check)
        return -1;
      fsr_id = _FSR_keys_map.at(fsr_key)->_fsr_id;
    }
    catch (std::exception &e)
    {
      if (err_check)
      {
        log::ferror("Could not find FSR ID with key: %s. Try creating "
                    "geometry with finer track spacing",
                    _fsr_keys[thread_id].c_str());
      }
      else
      {
        return -1;
      }
    }

    return fsr_id;
  }

  /**
   * @brief Returns the rank of the domain containing the given coordinates
   * @param coords The coordinates used to search the domain rank
   * @return The rank of the domain containing the coordinates
   */
  int Geometry::getDomainByCoords(LocalCoords *coords)
  {
    int domain = 0;
#ifdef ENABLE_MPI_
    if (isDomainDecomposed())
    {
      int domain_idx[3];
      double width_x = _root_universe->getMaxX() - _root_universe->getMinX();
      domain_idx[0] = (coords->getPoint()->getX() - _root_universe->getMinX()) * mpi::getNumDomainsX() / width_x;
      double width_y = _root_universe->getMaxY() - _root_universe->getMinY();
      domain_idx[1] = (coords->getPoint()->getY() - _root_universe->getMinY()) * mpi::getNumDomainsY() / width_y;
      double width_z = _root_universe->getMaxZ() - _root_universe->getMinZ();
      domain_idx[2] = (coords->getPoint()->getZ() - _root_universe->getMinZ()) * mpi::getNumDomainsZ() / width_z;

      MPI_Cart_rank(getMPICart(), domain_idx, &domain);
    }
#endif
    return domain;
  }

  /**
   * @brief Returns a map from cells to FSRs.
   * @return a map from cells to FSRs contained in those cells
   */
  std::map<Cell *, std::vector<long>> Geometry::getCellsToFSRs()
  {

    std::map<int, Cell *> all_cells = _root_universe->getAllCells();
    Cell *fsr_cell;
    std::map<Cell *, std::vector<long>> cells_to_fsrs;
    size_t num_FSRs = _FSR_keys_map.size();

    for (auto &c : all_cells)
    {

      /* Search for this Cell in all FSRs */
      for (size_t r = 0; r < num_FSRs; r++)
      {
        fsr_cell = findCellContainingFSR(r);
        if (c.first == fsr_cell->getId())
          cells_to_fsrs[c.second].push_back(r);
      }
    }

    return cells_to_fsrs;
  }

  /**
   * @brief Return the global ID of the flat source region that a given
   *        LocalCoords object resides within.
   * @param coords a LocalCoords object pointer
   * @param err_check whether to fail instead of returning -1 if not found
   * @return the FSR ID for a given LocalCoords object
   */
  long Geometry::getGlobalFSRId(LocalCoords *coords, bool err_check)
  {

    /* Check if the Geometry is domain decomposed */
    if (!isDomainDecomposed())
    {
      return getFSRId(coords, err_check);
    }
    else
    {
      long global_fsr_id = 0;
#ifdef ENABLE_MPI_
      // FIXME
      long temp_fsr_id = 0;
      int domain = getDomainByCoords(coords);
      int rank = mpi::getRankUniqueDomains();
      if (domain == rank)
        temp_fsr_id = getFSRId(coords);
      MPI_Allreduce(&temp_fsr_id, &global_fsr_id, 1, mpi::getDatatype<long>(), MPI_SUM,
                    getMPICart());

      /* Count FSRs on each domain if not already counted */
      if (!_domain_FSRs_counted)
        countDomainFSRs();

      /* Add FSR count from prior domains */
      for (long i = 0; i < domain; i++)
        global_fsr_id += _num_domain_FSRs.at(i);

#endif
      return global_fsr_id;
    }
  }

  /**
   * @brief Return the characteristic point for a given FSR ID
   * @param fsr_id the FSR ID
   * @return the FSR's characteristic point
   */
  Point *Geometry::getFSRPoint(long fsr_id)
  {

    Point *point = nullptr;

    try
    {
      std::string &key = _FSRs_to_keys[fsr_id];
      point = _FSR_keys_map.at(key)->_point;
    }
    catch (std::exception &e)
    {
      log::ferror("Could not find characteristic point in FSR: %d", fsr_id);
    }

    return point;
  }

  /**
   * @brief Return the centroid for a given FSR ID
   * @param fsr_id the FSR ID
   * @return the FSR's centroid
   */
  Point *Geometry::getFSRCentroid(long fsr_id)
  {

    if (fsr_id >= 0 && (size_t)fsr_id < _FSR_keys_map.size())
      return _FSRs_to_centroids[fsr_id];
    else
      log::ferror("Could not find centroid in FSR: %d.", fsr_id);
    return nullptr;
  }

  /**
   * @brief Returns whether any FSR centroids have been set
   */
  bool Geometry::containsFSRCentroids()
  {
    return _contains_FSR_centroids;
  }

  /**
   * @brief Return the CMFD cell for a given FSR ID
   * @param fsr_id the FSR ID
   * @return the CMFD cell
   */
  int Geometry::getCmfdCell(long fsr_id)
  {
    return _FSRs_to_CMFD_cells[fsr_id];
  }

  /**
   * @brief reserves the FSR key strings
   * @param num_threads the number of threads
   */
  void Geometry::reserveKeyStrings(int num_threads)
  {
    int string_size = 255;
    _fsr_keys.resize(num_threads);
    for (int i = 0; i < num_threads; ++i)
    {
      _fsr_keys[i].reserve(string_size);
    }
  }

  /**
   * @brief Generate a string FSR "key" that identifies an FSR by its
   *        unique hierarchical lattice/universe/cell structure.
   * @details Since not all FSRs will reside on the absolute lowest universe
   *          level and Cells might overlap other cells, it is important to
   *          have a method for uniquely identifying FSRs. This method
   *          creates a unique FSR key by constructing a structured string
   *          that describes the hierarchy of lattices/universes/cells.
   * @param coords a LocalCoords object pointer
   * @param key a reference to the FSR key
   * 给定局部坐标点链表,把链表对应的 CSG 对象信息构造  成字符串再返回
   * TODO: 用":"做分割符不太行
   */
  void Geometry::getFSRKeyFast(LocalCoords *coords, std::string &key)
  {

    // 先把 coords 升到最高层，读取当前点并计算最终 key 的长度，以便一次性 resize 字符串避免反复扩容。
    LocalCoords *curr = coords->getHighestLevel();

    /* Assess the string size of the key 统计字符串的长度*/
    Point *point = curr->getPoint();
    int total_size = 0;
    if (_cmfd != nullptr)
    {
      if (_cmfd->getLattice()->areValidIndices(_cmfd->getLattice()->getLatX(point),
                                               _cmfd->getLattice()->getLatY(point), _cmfd->getLattice()->getLatZ(point)))
      {
        total_size += getNumDigits(_cmfd->getLattice()->getLatX(point));
        total_size += getNumDigits(_cmfd->getLattice()->getLatY(point));
        total_size += getNumDigits(_cmfd->getLattice()->getLatZ(point));
        total_size += 9; // "CMFD(" + three digits包含的two commas + "):" (total punctuation)
      }
    }
    if (_overlaid_mesh != nullptr)
    {
      total_size += getNumDigits(_overlaid_mesh->getLatX(point));
      total_size += getNumDigits(_overlaid_mesh->getLatY(point));
      total_size += getNumDigits(_overlaid_mesh->getLatZ(point));
      total_size += 4; // "A", "R", ',' and ':' separators
    }
    while (curr != nullptr)
    {
      if (curr->getType() == LAT)
      {
        total_size += getNumDigits(curr->getLattice()->getId());
        total_size += getNumDigits(curr->getLatticeX());
        total_size += getNumDigits(curr->getLatticeY());
        total_size += getNumDigits(curr->getLatticeZ());
        total_size += 6; // "L", parentheses, commas, trailing ":"
      }
      else
      {
        total_size += getNumDigits(curr->getUniverse()->getId());
        total_size += 2; // 'U' prefix and ':' suffix
      }
      if (curr->getNext() == nullptr)
        break;
      else
        curr = curr->getNext();
    }
    total_size += getNumDigits(curr->getCell()->getId());
    total_size += 1; // 'C' prefix before cell id
    int version_num = coords->getVersionNum();
    if (version_num != 0)
    {
      total_size += getNumDigits(version_num);
      total_size += 2; // ":V"
    }

    /* Resize key resize 字符串*/
    if (total_size > 255)
      log::ferror("Found key exceeding the 255 character threshold");
    key.resize(total_size);
    curr = curr->getHighestLevel();

    // 长度确认后，从头逐段写入：可选的 CMFD 坐标、可选的全局覆盖网格、各层 lattice/universe 信息、最底层 Cell ID，再附加非零版本号
    /* If CMFD is on, get CMFD latice cell and write to key */
    int ind = 0;
    if (_cmfd != nullptr)
    {
      if (_cmfd->getLattice()->areValidIndices(_cmfd->getLattice()->getLatX(point),
                                               _cmfd->getLattice()->getLatY(point), _cmfd->getLattice()->getLatZ(point)))
      {
        key.replace(ind, 5, "CMFD(");
        ind += 5;

        printToString(key, ind, _cmfd->getLattice()->getLatX(point));
        key.at(ind) = ',';
        ind++;

        printToString(key, ind, _cmfd->getLattice()->getLatY(point));
        key.at(ind) = ',';
        ind++;

        printToString(key, ind, _cmfd->getLattice()->getLatZ(point));
        key.replace(ind, 2, "):");
        ind += 2;
      }
    }

    /* If a global overlaid mesh is present, get the axial mesh cell */
    if (_overlaid_mesh != nullptr)
    {
      key.at(ind) = 'A';
      ind++;
      printToString(key, ind, _overlaid_mesh->getLatZ(point));
      key.at(ind) = 'R';
      ind++;
      printToString(key, ind, _overlaid_mesh->getLatX(point));
      key.at(ind) = ',';
      ind++;
      printToString(key, ind, _overlaid_mesh->getLatY(point));
      key.at(ind) = ':';
      ind++;
    }

    /* Descend the linked list hierarchy until the lowest level has
     * been reached */
    while (curr != nullptr)
    {

      /* write lattice to key */
      if (curr->getType() == LAT)
      {
        key.at(ind) = 'L';
        ind++;
        printToString(key, ind, curr->getLattice()->getId());
        key.at(ind) = '(';
        ind++;
        printToString(key, ind, curr->getLatticeX());
        key.at(ind) = ',';
        ind++;
        printToString(key, ind, curr->getLatticeY());
        key.at(ind) = ',';
        ind++;
        printToString(key, ind, curr->getLatticeZ());
        key.replace(ind, 2, "):");
        ind += 2;
      }
      else
      {

        /* write universe ID to key */
        key.at(ind) = 'U';
        ind++;
        printToString(key, ind, curr->getUniverse()->getId());
        key.at(ind) = ':';
        ind++;
      }

      /* If lowest coords reached break; otherwise get next coords */
      if (curr->getNext() == nullptr)
        break;
      else
        curr = curr->getNext();
    }

    /* write cell id to key */
    key.at(ind) = 'C';
    ind++;
    printToString(key, ind, curr->getCell()->getId());

    /* write version number to key */
    if (version_num != 0)
    {
      key.replace(ind, 2, ":V");
      ind += 2;
      printToString(key, ind, version_num);
    }
  }

  // FIXME Find a better way to do this, without a function call
  /**Using std::stringstream would be more clear.
   * @brief Get the number of digits in base 10 of a number
   * @param number the number of interest
   * @return the number of digits in base 10 of a number
   *
   * 求数字的位数
   */
  int Geometry::getNumDigits(int number)
  {
    if (number < 0)
      log::ferror("Trying to get the digits of negative number %d", number);
    int ref = 10;
    int num_digits = 1;
    while (number >= ref)
    {
      num_digits++;
      ref *= 10;
    }
    return num_digits;
  }

  // FIXME Find a better way to do this, without a function call
  /**Using std::stringstream would be more clear.
   * @brief Print a number to a given String.
   * @param str the string to print to
   * @param index the last index in that string
   * @param value the number to print
   * @return the number of digits in base 10 of a number
   *
   * 把整数按十进制形式拷贝进目标字符串
   *
   * 读取 value 的十进制各位，逐位写入已经分配好的字符串 str。
   * 写入位置从 index 开始，写完以后把 index 增加写入的位数，方便后续继续往后写。
   * 为了避免依赖 std::stringstream，它先算出 value 有几位 (getNumDigits(value))，再通过 base = 10^(剩余位数-1) 的方式依次取出最高位到最低位，并用
   * 本地的 digits 查表映射成字符 '0'…'9'。
   */
  void Geometry::printToString(std::string &str, int &index, int value)
  {

    char digits[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

    int num_digits = getNumDigits(value);

    for (int i = 0; i < num_digits; i++)
    {
      int base = 1;
      for (int j = 0; j < num_digits - i - 1; j++)
        base *= 10;
      char digit = digits[value / base];
      str.at(index + i) = digit;
      value -= (value / base) * base;
    }
    index += num_digits;
  }

  /**
   * @brief Generate a string FSR "key" for the FSR where the point reside in. A
            string FSR "key" identifies an FSR by its unique hierarchical
   *        lattice/universe/cell structure.
   * @details Since not all FSRs will reside on the absolute lowest universe
   *          level and Cells might overlap other cells, it is important to
   *          have a method for uniquely identifying FSRs. This method
   *          creates a unique FSR key by constructing a structured string
   *          that describes the hierarchy of lattices/universes/cells.
   * @param coords a LocalCoords object pointer
   * @return the FSR key
   * 先给定局部坐标点链表,把链表对应的 CSG 对象信息构造成字符串再返回。
   * getFSRKey 用流ostringstream来处理
   * LocalCoords 的成员变量里有一个为version_num(版本号)。在创建 FSR 标识串时,这个值会加入到标识串中。
   * 也就是说,当两个局部坐标表示三维空间中同一个点 (x, y, z) 时,它们的版本号如果不同也会使他们对应不同的轴向挤出 FSR
   */
  std::string Geometry::getFSRKey(LocalCoords *coords)
  {

    std::stringstream key;
    LocalCoords *curr = coords->getHighestLevel();
    std::ostringstream curr_level_key;

    /* If CMFD is on, get CMFD latice cell and write to key */
    if (_cmfd != nullptr)
    {
      /*
      if (_cmfd->getLattice()->getLatX(curr->getPoint()) > 0 &&
          _cmfd->getLattice()->getLatY(curr->getPoint()) > 0 &&
          _cmfd->getLattice()->getLatZ(curr->getPoint()) > 0) {
        curr_level_key << _cmfd->getLattice()->getLatX(curr->getPoint());
        key << "CMFD(" << curr_level_key.str() << ",";
        curr_level_key.str(std::string());
        curr_level_key << _cmfd->getLattice()->getLatY(curr->getPoint());
        key << curr_level_key.str() << ",";
        curr_level_key.str(std::string());
        curr_level_key << _cmfd->getLattice()->getLatZ(curr->getPoint());
        key << curr_level_key.str() << "):";
      }*/
      if (_cmfd->getLattice()->areValidIndices(_cmfd->getLattice()->getLatX(curr->getPoint()),
                                               _cmfd->getLattice()->getLatY(curr->getPoint()), _cmfd->getLattice()->getLatZ(curr->getPoint())))
      {
        curr_level_key << _cmfd->getLattice()->getLatX(curr->getPoint());
        key << "CMFD(" << curr_level_key.str() << ",";
        curr_level_key.str(std::string());
        curr_level_key << _cmfd->getLattice()->getLatY(curr->getPoint());
        key << curr_level_key.str() << ",";
        curr_level_key.str(std::string());
        curr_level_key << _cmfd->getLattice()->getLatZ(curr->getPoint());
        key << curr_level_key.str() << "):";
      }
    }

    /* If a global overlaid mesh is present, record mesh cells */
    if (_overlaid_mesh != nullptr)
    {
      curr_level_key.str(std::string());
      curr_level_key << _overlaid_mesh->getLatZ(curr->getPoint());
      key << "A" << curr_level_key.str() << ":";
      curr_level_key.str(std::string());
      curr_level_key << _overlaid_mesh->getLatX(curr->getPoint());
      key << "R" << curr_level_key.str() << ",";
      curr_level_key.str(std::string());
      curr_level_key << _overlaid_mesh->getLatY(curr->getPoint());
      key << curr_level_key.str() << ":";
    }

    /* Descend the linked list hierarchy until the lowest level has
     * been reached */
    while (curr != nullptr)
    {

      /* Clear string stream */
      curr_level_key.str(std::string());

      if (curr->getType() == LAT)
      {

        /* Write lattice ID and lattice cell to key */
        curr_level_key << curr->getLattice()->getId();
        key << "L" << curr_level_key.str() << "(";
        curr_level_key.str(std::string());
        curr_level_key << curr->getLatticeX();
        key << curr_level_key.str() << ",";
        curr_level_key.str(std::string());
        curr_level_key << curr->getLatticeY();
        key << curr_level_key.str() << ",";
        curr_level_key.str(std::string());
        curr_level_key << curr->getLatticeZ();
        key << curr_level_key.str() << "):";
      }
      else
      {
        /* write universe ID to key */
        curr_level_key << curr->getUniverse()->getId();
        key << "U" << curr_level_key.str() << " : ";
      }

      /* If lowest coords reached break; otherwise get next coords */
      if (curr->getNext() == nullptr)
        break;
      else
        curr = curr->getNext();
    }

    /* clear string stream */
    curr_level_key.str(std::string());

    /* write cell id to key */
    curr_level_key << curr->getCell()->getId();
    key << "C" << curr_level_key.str();

    /* write version number to key */
    int version_num = coords->getVersionNum();
    if (version_num != 0)
    {
      curr_level_key.str(std::string());
      curr_level_key << version_num;
      key << ":V" << curr_level_key.str();
    }

    return key.str();
  }

  /**
   * @brief Return a pointer to an ExtrudedFSR by its extruded FSR ID
   * @param extruded_fsr_id the extruded FSR ID
   * @return a pointer to the ExtrudedFSR
   */
  ExtrudedFSR *Geometry::getExtrudedFSR(int extruded_fsr_id)
  {
    return _extruded_FSR_lookup[extruded_fsr_id];
  }

  /**
   * @brief 将几何中所有单元格（Cell）细分为与z轴对齐的环（rings）和扇区（angular sectors）。
   * @details 此方法默认由Geometry::initializeFlatSourceRegions()方法调用，
   *          若有需要，用户也可在Python中直接调用，示例如下：
   *
   * @code
   *          geometry.subdivideCells()  // Python中的调用示例
   * @endcode
   */
  void Geometry::subdivideCells()
  {

    /* 计算与当前几何模型面积等效的“等效半径”
    该半径将作为所有需“环形细分”的单元格的最大半径限制
    */
    // width_x：_root_universe在x方向的宽度（最大x坐标 - 最小x坐标）
    double width_x = _root_universe->getMaxX() - _root_universe->getMinX();
    // width_y：_root_universe在y方向的宽度（最大y坐标 - 最小y坐标）
    double width_y = _root_universe->getMaxY() - _root_universe->getMinY();
    // max_radius：等效半径，计算逻辑为“对角线长度的一半”（勾股定理）
    double max_radius = sqrt(width_x * width_x + width_y * width_y) / 2;

    /* 递归地将所有单元格细分为环和扇区 */
    _root_universe->subdivideCells(max_radius);
  }

  /**
   * @brief Compute the number of flat source regions in the Geometry and
   *        initialize CMFD. - 计算几何中的平面源区域数量并初始化 CMFD
   * @details This method is intended to be called by the user before initiating
   *          source iteration. This method first subdivides all Cells by calling
   *          the Geometry::subdivideCells() method. Then it initializes the CMFD
   *          object. - 该方法应在开始源迭代前调用：先通过 Geometry::subdivideCells() 细分所有单元，再初始化 CMFD 对象
   */
  void Geometry::initializeFlatSourceRegions()
  {

    log::finfo("Initializing flat source regions...");

    /* Subdivide Cells into sectors and rings - 将单元细分为扇区与环形 */
    subdivideCells();

    /* Build collections of neighbor Cells for optimized ray tracing - 为优化光线追踪构建邻居单元集合 */
    // FIXME
    //_root_universe->buildNeighbors();

    /* Create map of Material IDs to Material pointers - 建立材料 ID 到材料指针的映射 */
    retrieveAllMaterials();
    auto &all_materials = getAllMaterials();

    /* Initialize absorption XS if CMFD absent - 若CMFD模块为空，则初始化材料的吸收截面（Σₐ） */
    if (_cmfd == nullptr)
    {
      /*
      const auto &pair：
      const：不修改元素
      auto：自动推导类型（这里 pair 是键值对）
      &：引用，避免复制
      */
      for (const auto &pair : all_materials) // 遍历：遍历材料映射表中的每个键值对
        pair.second->getSigmaA();            // 初始化时预先计算每个材料在所有能群下的吸收截面
    }

    if (_cmfd != nullptr) /* 若CMFD模块已存在，则初始化CMFD */
      initializeCmfd();
  }
  /*
   */

  /**
   * @brief This method performs ray tracing to create Track segments within each
   *        flat source region in the Geometry.
   * @details This method starts at the beginning of a Track and finds successive
   *          intersection points with FSRs as the Track crosses through the
   *          Geometry and creates segment structs and adds them to the Track.
   * @param track a pointer to a track to segmentize
   * @param z_coord the axial height at which the 2D plane of the geometry is
   *        formed
   */
  void Geometry::segmentize2D(Track *track, double z_coord, bool FSRonly)
  {

    /* Track starting Point coordinates and azimuthal angle */
    double x0 = track->getStart()->getX();
    double y0 = track->getStart()->getY();
    double z0 = z_coord;
    double phi = track->getPhi();
    double delta_x, delta_y;
    // unused
    // double delta_z;

    /* ID of each FSR containing segments */
    int fsr_id;

    /* Use a LocalCoords for the start and end of each segment */
    LocalCoords start(x0, y0, z0, true);
    LocalCoords end(x0, y0, z0, true);
    start.setUniverse(_root_universe);
    end.setUniverse(_root_universe);

    /* Find the Cell containing the Track starting Point */
    Cell *curr = findFirstCell(&end, phi);
    Cell *prev;

    /* If starting Point was outside the bounds of the Geometry */
    if (curr == nullptr)
      log::ferror("Could not find a Cell containing the start Point "
                  "of this Track: %s",
                  track->toString().c_str());

    /* While the end of the segment's LocalCoords is still within the Geometry,
     * move it to the next Cell, create a new segment, and add it to the
     * Geometry */
    while (curr != nullptr)
    {

      end.copyCoords(&start);

      /* Find the next Cell along the Track's trajectory */
      // 找下一个Cell，主要算法所在
      prev = curr;
      curr = findNextCell(&end, phi);

      /* Checks that segment does not have the same start and end Points */
      if (fabs(start.getX() - end.getX()) < FLT_EPSILON && fabs(start.getY() - end.getY()) < FLT_EPSILON)
        log::error("Created 2D segment with same start and end point: "
                   "x = {}, y = {}, z = {}",
                   start.getX(), start.getY(), start.getZ());

      // Find FSR, create a new FSR if it doesn't exist.
      fsr_id = findFSRId(&start);

      // Create a new Track segment as needed.
      segment new_segment;
      if (!FSRonly)
      {
        // Set the segment length, Material and FSR ID.
        new_segment._length = double(end.getPoint()->distanceToPoint(start.getPoint()));
        new_segment._material = prev->getFillMaterial();
        new_segment._region_id = fsr_id;
      }

      log::debug("segment start x = {}, y = {}; end x = {}, y = {}",
                 start.getX(), start.getY(), end.getX(), end.getY());

      /* Save indices of CMFD Mesh surfaces that the Track segment crosses */
      if (!FSRonly && _cmfd != nullptr)
      {

        /* Find cmfd cell that segment lies in */
        int cmfd_cell = _cmfd->findCmfdCell(&start); // 找到起点所在的 CMFD cell。

        if (cmfd_cell < 0)
          break;

        /* Reverse nudge from surface to determine whether segment start or end
         * points lie on a CMFD surface. */
        delta_x = cos(phi) * TINY_MOVE;
        delta_y = sin(phi) * TINY_MOVE;
        start.adjustCoords(-delta_x, -delta_y);
        end.adjustCoords(-delta_x, -delta_y); // 为防判断在边界上产生歧义，沿着反方向微调 start/end（相当于“从表面往里面推一点”）。

        /* Calculate CMFD surfaces */
        int cmfd_surfaces[2];
        cmfd_surfaces[0] = _cmfd->findCmfdSurface(cmfd_cell, &end, phi, M_PI / 2.0);
        cmfd_surfaces[1] = _cmfd->findCmfdSurface(cmfd_cell, &start, phi, M_PI / 2.0);

        /* Ensure surfaces are x-y surfaces (no z-crossings) */
        /* Note: this code takes advantage of the numeric representation of
           surfaces to find a mapping that removes z-surfaces */
        // 接下来的 for 循环清理/映射 surface index：去掉 z 方向面的编码，只保留 x/y 方向，或把复杂编号归约。这部分利用了 NUM_SURFACES 的编号规则。
        for (int d = 0; d < 2; d++)
        {
          int local_surface = cmfd_surfaces[d] % NUM_SURFACES;
          if (local_surface == 2 || local_surface == 5)
          {
            cmfd_surfaces[d] = -1;
          }
          else if (local_surface > 9)
          {
            int cell = cmfd_surfaces[d] / NUM_SURFACES;
            int half_surf = local_surface / 2;
            if (local_surface > 17)
            {
              int quart_surf = half_surf / 2;
              local_surface = 2 + quart_surf + (half_surf == 2 * quart_surf);
              cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface;
            }
            else
            {
              local_surface = (half_surf > 6) + 3 *
                                                    (local_surface != 2 * half_surf);
              cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface;
            }
          }
        }

        /* Save CMFD surfaces */
        new_segment._cmfd_surface_fwd = cmfd_surfaces[0];
        new_segment._cmfd_surface_bwd = cmfd_surfaces[1];

        /* Re-nudge segments from surface. */
        // 最后把 start/end 再调回原位（加回 delta_x/delta_y）。
        start.adjustCoords(delta_x, delta_y);
        end.adjustCoords(delta_x, delta_y);
      }

      if (!FSRonly)
      {
        /* Calculate the local centroid of the segment if available */
        // 如果FSR中心可用,计算线段起点相对于FSR中心的局部坐标
        //  FIXME Consider reversing nudge
        Point *starting_point = start.getHighestLevel()->getPoint();
        new_segment._starting_position[0] = starting_point->getX();
        new_segment._starting_position[1] = starting_point->getY();
        // if (_contains_FSR_centroids)/ 一个布尔值，指示是否已设置任何质心
        {
          Point *centroid = getFSRCentroid(fsr_id); // FSR 的质心
          double x_start = starting_point->getX() - centroid->getX();
          double y_start = starting_point->getY() - centroid->getY();
          new_segment._starting_position[0] = x_start;
          new_segment._starting_position[1] = y_start;
        }

        /* Add the segment to the Track */
        track->addSegment(new_segment);
      }

    } // end of while

    log::fdebug("Created %d segments for Track: %s",
                track->getNumSegments(), track->toString().c_str());

    /* Truncate the linked list for the LocalCoords */
    start.prune();
    end.prune();
  }

  /**
   * @brief This method performs ray tracing to create Track segments within each
   *        flat source region in the Geometry.执行射线追踪，以在几何体中的每个平面源区域内创建轨迹段
   * @details This method starts at the beginning of a Track and finds successive
   *          intersection points with FSRs as the Track crosses through the
   *          Geometry and creates segment structs and adds them to the Track.
   * 从轨迹开始，在轨迹穿过几何体时找到与FSR的连续交点，并创建分段结构并将其添加到轨迹中
   * @param track a pointer to a track to segmentize
   * @param OTF_setup whether this routine is called during OTF ray tracing setup  传入的OTF_setup = true
   */
  void Geometry::segmentize3D(Track3D *track, bool OTF_setup)
  {

    /* Track starting Point coordinates and azimuthal angle */
    double x0 = track->getStart()->getX();
    double y0 = track->getStart()->getY();
    double z0 = track->getStart()->getZ();
    double phi = track->getPhi(); // 方位角和极角都是0
    double theta = track->getTheta();

    /* Length of each segment */
    double length;

    /* Use a LocalCoords for the start and end of each segment */
    LocalCoords start(x0, y0, z0, true);
    LocalCoords end(x0, y0, z0, true);
    start.setUniverse(_root_universe);
    end.setUniverse(_root_universe);

    if (x0 != x0 || y0 != y0)
    {
      log::ferror("Nan is found in starting points at the beginning of "
                  "Geometry::segmentize3D: x = %f, y = %f",
                  x0, y0);
    }

    /* Find the Cell containing the Track starting Point 找到包含该起点的cell*/
    Cell *curr = findFirstCell(&end, phi, theta);
    Cell *prev;

    /* Vector to fill coordinates if necessary */
    std::vector<LocalCoords *> fsr_coords; // 目前没用上，不管

    /* If starting Point was outside the bounds of the Geometry */
    if (curr == nullptr)
      log::ferror("Could not find a Cell containing the start Point "
                  "of this Track3D: %s",
                  track->toString().c_str());

    /* While the end of the segment's LocalCoords is still within the Geometry,
     * move it to the next Cell, create a new segment, and add it to the
     * Geometry 当线段的LocalCoords末端仍在几何体中时，将其移动到下一个单元，创建一个新线段，并将其添加到几何体*/
    while (curr != nullptr)
    {

      end.copyCoords(&start);

      /* Find the next Cell along the Track's trajectory */
      prev = curr;
      curr = findNextCell(&end, phi, theta); // 与径向找最小线段类似，这个是在垂直方向上找最小线段

      /* Checks to make sure that new Segment does not have the same start
       * and end Points 检查以确保新线段的起点和终点不相同*/
      if (fabs(start.getX() - end.getX()) < FLT_EPSILON &&
          fabs(start.getY() - end.getY()) < FLT_EPSILON &&
          fabs(start.getZ() - end.getZ()) < FLT_EPSILON)
      {
        log::ferror("Created a Track3D segment with the same start and end "
                    "point: x = %f, y = %f, z = %f",
                    start.getX(),
                    start.getY(), start.getZ());
      }

      /* Find the segment length between the segment's start and end points 计算线段长度*/
      length = double(end.getPoint()->distanceToPoint(start.getPoint()));
      long fsr_id = findFSRId(&start); // 找到该起点处的FSR_ID

      /* Create a new Track segment */
      segment new_segment;
      new_segment._material = prev->getFillMaterial();
      new_segment._length = length;
      new_segment._region_id = fsr_id;

      log::fdebug("segment start x = %f, y = %f, z = %f; "
                  "end x = %f, y = %f, z = %f",
                  start.getX(), start.getY(), start.getZ(),
                  end.getX(), end.getY(), end.getZ());

      /* Save indices of CMFD Mesh surfaces that the Track segment crosses */
      if (_cmfd != nullptr && !OTF_setup)
      { // 由于OTF_setup = true 这个if语句不执行
        /* Find cmfd cell that segment lies in */
        int cmfd_cell = _cmfd->findCmfdCell(&start);

        if (cmfd_cell >= 0)
        {
          /* Reverse nudge from surface to determine whether segment start or end
           * points lie on a cmfd surface. */
          double delta_x = cos(phi) * sin(theta) * TINY_MOVE;
          double delta_y = sin(phi) * sin(theta) * TINY_MOVE;
          double delta_z = cos(theta) * TINY_MOVE;
          start.adjustCoords(-delta_x, -delta_y, -delta_z);
          end.adjustCoords(-delta_x, -delta_y, -delta_z);

          new_segment._cmfd_surface_fwd =
              _cmfd->findCmfdSurface(cmfd_cell, &end, phi, theta);
          new_segment._cmfd_surface_bwd =
              _cmfd->findCmfdSurface(cmfd_cell, &start, phi, theta);

          /* Re-nudge segments from surface. */
          start.adjustCoords(delta_x, delta_y, delta_z);
          end.adjustCoords(delta_x, delta_y, delta_z);
        }
      }

      /* For regular 3D tracks, get starting position relative to FSR centroid */

      if (!OTF_setup)
      {
        Point *starting_point = start.getHighestLevel()->getPoint();
        new_segment._starting_position[0] = starting_point->getX();
        new_segment._starting_position[1] = starting_point->getY();
        new_segment._starting_position[2] = starting_point->getZ();
        if (_contains_FSR_centroids)
        {
          Point *centroid = getFSRCentroid(fsr_id);
          double x_start = starting_point->getX() - centroid->getX();
          double y_start = starting_point->getY() - centroid->getY();
          double z_start = starting_point->getZ() - centroid->getZ();
          new_segment._starting_position[0] = x_start;
          new_segment._starting_position[1] = y_start;
          new_segment._starting_position[2] = z_start;
        }
      }

      /* Add the segment to the Track */
      track->addSegment(new_segment); // 添加轨迹段
    }

    log::fdebug("Created %d segments for Track3D: %s",
                track->getNumSegments(), track->toString().c_str());

    /* Truncate the linked list for the LocalCoords */
    start.prune();
    end.prune();
  }

  /**
   * @brief This method performs ray tracing to create extruded track segments
   *        within each flat source region in the implicit 2D superposition plane
   *        of the Geometry by 2D ray tracing across input heights that encompass
   *        all radial geometric details.
   * @details This method starts at the beginning of an extruded track and finds
   *          successive intersection points with FSRs as the extruded track
   *          crosses radially through the Geometry at defined z-coords. The
   *          minimum distance to intersection of all z-coords is chosen leading
   *          to implicitly capturing all geometric radial detail at the defined
   *          z-heights, saving the lengths and region IDs to the extruded track
   *          and initializing ExtrudedFSR structs in the traversed FSRs.
   * @param flattened_track a pointer to a 2D track to segmentize into regions of
   *        extruded FSRs
   * @param z_coords a vector of axial heights in the root geometry at which
   *        the Geometry is segmentized radially
   * flattened_track：要分段的这条 2D 轨迹
   * z_coords：所有轴向层的 z 高度。
   * FSRonly：如果为 true，只建 FSR 不真正往轨迹上加 segment（一般为 false）。
   */
  void Geometry::segmentizeExtruded(Track *flattened_track,
                                    std::vector<double> z_coords, bool FSRonly)
  {

    /* Track starting Point coordinates and azimuthal angle */
    double x0 = flattened_track->getStart()->getX();
    double y0 = flattened_track->getStart()->getY();
    double z0 = z_coords[0];
    double phi = flattened_track->getPhi();
    double delta_x, delta_y; // 后面处理 CMFD 时需要做一个微小平移（nudging）
    // unused
    // double delta_z;

    log::fdebug("segmentizeExtruded starts at (%f, %f, %f) along "
                "angle %f",
                x0, y0, z0, phi);

    /* Length of each segment */
    double length;
    int min_z_ind; // 对于当前段来说，哪个 z 层给出的交点最近;最短距离对应的z_coords索引
    int region_id; // 找到的 ExtrudedFSR 的全局 id（也就是 FSR id）

    /* Use a LocalCoords for the start and end of each segment */
    LocalCoords start(x0, y0, z0, true);
    LocalCoords end(x0, y0, z0, true);
    start.setUniverse(_root_universe);
    end.setUniverse(_root_universe);

    /* Find the Cell containing the Track starting Point*/
    Cell *curr = findFirstCell(&end, phi); // 起点所在的最底层 Cell，这里的 &end 不是“按引用传参”，而是取对象 end 的地址

    /* If starting Point was outside the bounds of the Geometry */
    if (curr == nullptr)
    {
      int dom = mpi::getDomainUid();
      log::ferror("Could not find a Cell containing the start Point "
                  "of this Track on domain %d with bounds [%f, %f] x [%f, %f]"
                  " x [%f, %f]. Track: %s. Current coords: %s",
                  dom, getMinX(),
                  getMaxX(), getMinY(), getMaxY(), getMinZ(), getMaxZ(),
                  flattened_track->toString().c_str(),
                  end.getPoint()->toString().c_str());
    }

    /* While the end of the segment's LocalCoords is still within the Geometry,
     * move it to the next Cell, create a new segment, and add it to the
     * Geometry */
    int find_cell_count = 0;
    while (curr != nullptr) // 大循环：沿着轨迹一步步往前走,只要终点 end 仍然落在几何里的某个 Cell 中，就继续分下一段
    {

      /* Check if stuck in loop */
      find_cell_count++; // find_cell_count 防止无限循环
      if (find_cell_count > 1e6)
        log::ferror("Caught in infinite loop finding next cell");

      /* Records the minimum length to a 2D intersection */
      double min_length = std::numeric_limits<double>::infinity();
      region_id = -1;
      min_z_ind = -1; //

      end.copyCoords(&start); // 会把 end 的坐标和整条 LocalCoords 链复制给 start（把终点作为起点），从而开启另一段最小线段的寻找

      /* Loop over all z-heights to find shortest 2D intersection */
      for (size_t i = 0; i < z_coords.size(); i++) // 内层循环：对每一个 z 层试一遍交点，选最近的那个
      {

        /* Change z-height and copy starting coordinates to end */
        // 不同z对应的点所在层不同，得到的Cell也可能不同
        start.setZ(z_coords[i]);
        start.prune(); // 只保留这一级 LocalCoords，把下面的链都删掉，避免上一次循环的 CSG 信息干扰现在这一次；

        findCellContainingCoords(&start);
        start.copyCoords(&end); // 把start的值赋给end，重新寻找最小的线段长度

        /* Find the next Cell along the Track's trajectory */
        curr = findNextCell(&end, phi); // 找从当前点出发，沿方向 phi，最近会碰到哪一个几何边界，并把 end 移动到那个交点；同时返回交点所在的 Cell。

        /* Checks that segment does not have the same start and end Points 起点和终点不能相同*/
        if (fabs(start.getX() - end.getX()) < FLT_EPSILON &&
            fabs(start.getY() - end.getY()) < FLT_EPSILON)
          log::ferror("Created segment with same start and end "
                      "point: x = %f, y = %f",
                      start.getX(), start.getY());
        if (end.getX() != end.getX() ||
            end.getY() != end.getY())
        {
          log::ferror("Nan is found in points when loop over z-heights "
                      "during segmenting: x = %f, y = %f, z_index = %d",
                      end.getX(), end.getY(), i);
        }

        /* Find the segment length and extruded FSR */
        length = double(end.getPoint()->distanceToPoint(start.getPoint())); // 计算特定Z轴高度下的最小的长度

        /* Check if the segment length is the smallest found */
        if (length < min_length)
        { // 把最短线段的长度和所在轴向层高度的索引记录下来
          min_length = length;
          min_z_ind = i;
        }
      }
      // 到这里为止： 我们已经知道了，从当前起点出发，在所有轴向层里，哪个 z 层的面最先被撞到，以及对应的最短距离 min_length 和层索引 min_z_ind。

      // 把 start 的 CSG 链清掉；把 z 设置成刚才找到的那一层；再根据这个坐标重新找 Cell，建好 LocalCoords 链。
      start.prune();
      start.setZ(z_coords[min_z_ind]);
      findCellContainingCoords(&start);

      // Create a unique ExtrudedFSR by comparing ExtrudedFSRs with the same CSG
      region_id = createUniqueExtrudedFSR(start, z_coords); // 根据最小线段的Z高度来创建轴向挤出FSR

      /* Move the coordinates to the next intersection */
      start.copyCoords(&end);
      curr = findNextCell(&end, phi, M_PI_2); // 在选定的层上重新走一遍 findNextCell，移动end到终点，并返回其所在的Cell

#ifdef ENABLE_DEBUG_
      log::fdebug("segment start x = %f, y = %f; end x = %f, y = %f",
                  start.getX(), start.getY(), end.getX(), end.getY());
#endif

      /* Check that the region ID is valid */
      if (region_id == -1)
        log::ferror("Failed to find a valid FSR during axial extruded "
                    "segmentation");

      /* Create a new 2D Track segment with extruded region ID  给 segment 填数据*/
      segment new_segment;
      if (!FSRonly)
      { // FSRonly 默认为false
        new_segment._length = min_length;
        new_segment._region_id = region_id; // 这一段所在的 ExtrudedFSR 编号
      }

      /* Save indices of CMFD Mesh surfaces that the Track segment crosses 保存轨迹段所穿过的CMFD网格surface面的索引*/
      if (!FSRonly && _cmfd != nullptr)
      {

        /* Find cmfd cell that segment lies in */
        int cmfd_cell = _cmfd->findCmfdCell(&start); // 计算线段起点所在的CMFD cell，如果不在当前域上（返回 -1），就不处理 CMFD，直接跳到 add_segment 标签处。

        if (cmfd_cell == -1)
        {
          log::fdebug("segment not in CMFD Cell start(%lf,%lf), end(%lf,%lf)", // 该线段不在当前域中
                      start.getX(), start.getY(), end.getX(), end.getY());
          goto add_segment;
        }

        if (_cmfd->GetHexLatticeEnable())
        { // 000000000000
          if (!_cmfd->CellinHexLattice(cmfd_cell))
          {
            log::fdebug("segment not in Hex CMFD Cell start(%lf,%lf), end(%lf,%lf)",
                        start.getX(), start.getY(), end.getX(), end.getY());
            goto add_segment;
          }
        }

        /* Reverse nudge from surface to determine whether segment start or end
         * points lie on a CMFD surface. */
        delta_x = cos(phi) * TINY_MOVE;
        delta_y = sin(phi) * TINY_MOVE;
        start.adjustCoords(-delta_x, -delta_y, 0); // 为了避免起点/终点正好落在 CMFD 网格的顶点或面上，会先往反方向退一点点（TINY_MOVE）； 这有助于判断“是穿过哪一个面”
        end.adjustCoords(-delta_x, -delta_y, 0);

        /* Calculate CMFD surfaces */
        int cmfd_surfaces[2];

        // cmfd_surfaces[0] = _cmfd->findCmfdSurface(cmfd_cell, &end, phi, M_PI_2);
        // cmfd_surfaces[1] = _cmfd->findCmfdSurface(cmfd_cell, &start, phi, M_PI_2);

        cmfd_surfaces[0] = _cmfd->findCmfdSurface(cmfd_cell, &end, phi, M_PI_2);          //  从终点沿正方向 phi 延伸，看它是从哪个 CMFD 面出界；
        cmfd_surfaces[1] = _cmfd->findCmfdSurface(cmfd_cell, &start, phi + M_PI, M_PI_2); // 从起点沿反方向 phi+π 看，是从哪个 CMFD 面进来的。

        /* Ensure surfaces are x-y surfaces (no z-crossings) */
        /* Note: this code takes advantage of the numeric representation of
           surfaces to find a mapping that removes z-surfaces   删除z表面的映射，转为XY平面映射*/
        _cmfd->GetXYSurfaces(cmfd_surfaces);

        /* Save CMFD surfaces */
        new_segment._cmfd_surface_fwd = cmfd_surfaces[0];
        new_segment._cmfd_surface_bwd = cmfd_surfaces[1];

        /* Re-nudge segments from surface. */
        start.adjustCoords(delta_x, delta_y, 0);
        end.adjustCoords(delta_x, delta_y, 0);
      }

    add_segment: // goto标签
      if (!FSRonly)
      {
        /* Add the segment to the 2D track */
        flattened_track->addSegment(new_segment); // 把这段加入到 2D 轨迹的 segments list
        // log::fdebug("segment start(%lf,%lf), end(%lf,%lf), length:%lf, fwd:(%d,%d), bwd:(%d,%d)",
        //             start.getX(), start.getY(), end.getX(), end.getY(),
        //             new_segment._length,new_segment._cmfd_surface_fwd,new_segment._cmfd_surface_fwd % HEX_NUM_SURFACES,
        //             new_segment._cmfd_surface_bwd,new_segment._cmfd_surface_bwd % HEX_NUM_SURFACES);
      }
    }

    /* Truncate the linked list for the LocalCoords  把 LocalCoords 链表后面的节点都删掉，释放内存。*/
    start.prune();
    end.prune();
  }

  /// \brief Create a unique ExtrudedFSR on the superposition plane.
  /// \details ExtrudedFSRs are different from each other in two ways:
  ///          1. which have different CSG information are different,
  ///             just like 3-D FSRs.
  ///          2. which have different axial information are different,
  ///             which is distinguished by a version number.
  ///          Thus, unique ExtrudedFSRs which have the same CSG information
  ///          must have different version numbers.
  /// \param start A point in the ExtrudedFSR whose version number will be set.
  /// \param z_coords Z-coordinates of the superposition plane.
  /// \return ID of the newly created ExtrudedFSR.
  /*
  start 是当前轨迹段中、落在这个挤出区里的一个局部坐标点；
  这个函数会修改它的 version 号，让它指向对应的 ExtrudedFSR。
  z_coords：重叠平面上的所有 z 层高度列表；返回：这个 ExtrudedFSR 的 id（也就是 segment 的 _region_id）
  */
  int Geometry::createUniqueExtrudedFSR(LocalCoords &start,
                                        const std::vector<double> &z_coords)
  {
    // Default region id
    int region_id = -1;

    // Create two localCoords to check results
    LocalCoords test_ext_coords(0, 0, 0, true);
    LocalCoords test_start_coords(0, 0, 0, true);

    // unused
    // bool found_coordinate = false;

    int next_version = 0; // 统计已经尝试了几个 version 号
    // 外层循环：每次用一个新的 version 号来尝试创建 一个 ExtrudedFSR；
    // 理论上大部分情况 v=0 一次就成功，只有在比较极端的几何情况才会走到 v=1,2,...。
    for (int v = 0; v < MAX_VERSION_NUM; v++)
    {

      // Find FSR ID using starting coordinate, create a new FSR if there is not an FSR associated with the LocalCoords.
      start.setVersionNum(v);
      region_id = findExtrudedFSR(&start);
      std::string fsr_key = getFSRKey(&start); // start所在FSR的标识串，在2D情况下用的是getFSRKeyFast，这里再取一遍，是为了后面从 _extruded_FSR_keys_map 里取出刚刚关联的那块挤出区。

      // 获取extruded FSR的坐标
      // 有两种可能的结果： (1) 检索到的坐标是我们刚刚处理的坐标，也就是说，start==retrieved_coords
      // (2) 检索到的坐标是在另一个线程中创建的。
      // 对于 (1)，没有什么可比较的，因为它们是相同的。
      // 对于 (2)，我们必须检查检查当前局部点和表中存储的 FSR 特征点是否位于重叠平面上的同一 FSR 中，否则我们必须为 LocalCoords 创建一个新的挤出 FSR。新创建的 FSR 与现有 FSR 之间的区别就是所谓的版本号。
      LocalCoords *volatile retrieved_coords = nullptr;
      // volatile 是 C/C++ 中的一个关键字，用于修饰可能被意外改变的变量。其主要作用是告诉编译器不要对该变量进行优化，每次都直接从内存中读取其值
      // 多线程使用：在多线程环境中，当多个线程访问并修改同一变量时，使用 volatile 可以防止编译器将该变量的值缓存到寄存器中，从而确保每次都能获取最新的值。
      // 误区：volatile 不能保证原子性，也不是线程锁，不能防止竞争条件。
      do
      {
        retrieved_coords = _extruded_FSR_keys_map.at(fsr_key)->_coords;
      } while (retrieved_coords == nullptr);

      LocalCoords *ext_coords = retrieved_coords;

      // 拷贝一份坐标信息，避免直接动原对象
      // 接下来我们要反复修改它们的 z 坐标、调用 prune()、重新跑 findCellContainingCoords；
      ext_coords->copyCoords(&test_ext_coords);
      start.copyCoords(&test_start_coords);

      bool coords_contained = true; // 初始设为 true，表示“默认认为这两个点在所有 z 层上都属于同一个挤出区”
      // 遍历轴向层，看看两个点在所有 z 层的CSG层次信息是否都一致
      for (size_t i = 0; i < z_coords.size(); i++)
      {

        /* Check the FSR key at this level */
        test_start_coords.setZ(z_coords[i]);
        test_start_coords.prune();
        test_start_coords.setVersionNum(0); // 注意：这里特地把 version 号设回 0； 因为我们现在是在“比较几何结构是否一致”，不希望 version 号影响 FSR key
        // 取得该点的CSG层次信息，在这个新高度下，根据 (x,y,z) 重新在整个几何里找它落在哪个 Cell；同时重建 LocalCoords 的 CSG 链（Universe → Cell → Lattice…）
        findCellContainingCoords(&test_start_coords);
        fsr_key = getFSRKey(&test_start_coords);

        test_ext_coords.setZ(z_coords[i]);
        test_ext_coords.prune();
        test_ext_coords.setVersionNum(0);
        // 取得该点的CSG层次信息
        findCellContainingCoords(&test_ext_coords);
        std::string ext_fsr_key = getFSRKey(&test_ext_coords);

        /* Check that FSR keys match */
        if (fsr_key != ext_fsr_key)
        {
          coords_contained = false;
          break;
        }
      }

      /* Check if we found a valid coordinate */
      // 前面比较的两个点如果信息都一致，这里就会跳出循环
      if (coords_contained)
      {
        // unused, FIXME
        // found_coordinate = true;
        break;
      }

      // A unique ExtrudedFSR was found, increase the version number
      next_version++;
    }

    if (next_version >= MAX_VERSION_NUM)
    {
      std::string fsr_key = getFSRKey(&start);
      log::error("Exceeded the maximum version number of 2D extruded FSRs, FSR:\n{}\n{}",
                 fsr_key.c_str(), start.toString().c_str());
    }

    return region_id;
  }

  /**
   * @brief Fixes the FSR map size so that the map is static and fast
   */
  void Geometry::fixFSRMaps()
  {
    _FSR_keys_map.setFixedSize();
    _extruded_FSR_keys_map.setFixedSize();
  }

  /**
   * @brief Rays are shot vertically through each ExtrudedFSR struct to calculate
   *        the axial mesh and initialize 3D FSRs
   * @details From a 2D point within each FSR, a temporary 3D track is created
   *          starting at the bottom of the geometry and extending vertically to
   *          the top of the geometry. These tracks are segmented using the
   *          segmentize3D routine to calculate the distances between axial
   *          intersections forming the axial meshes if necessary and
   *          initializing the 3D FSRs as new regions are traversed.
   * 现在，我们已经有了一堆 ExtrudedFSR（每个代表 x-y 平面里的一块“柱底”区域），但还没真正生成 3D FSR。
   * 这个函数会对 每一个 ExtrudedFSR 做：从它的 2D 特征点出发，在 z 方向打出一条竖直 3D 轨迹，从整个几何底部 min_z 到顶部 max_z，
   * 用 segmentize3D 分段（或者直接用统一的 global_z_mesh 切分），得到一串 3D FSR 的材料、fsr_id、以及轴向网格 _mesh。
   * 最后统计总共有多少 3D FSR，并且重新整理一下 FSR 的编号。
   * 重点求：3D FSR的_FSR_keys_map的key，和value（FSRData）
   * extruded_FSR->_num_fsrs、extruded_FSR->_materials、extruded_FSR->_fsr_ids和extruded_FSR->_mesh（如果不使用全局轴向网）
   * @param global_z_mesh A global z mesh used for ray tracing. If the vector's
   *        length is zero, z meshes are local and need to be created for every
   *        ExtrudedFSR.
   */
  void Geometry::initializeAxialFSRs(std::vector<double> global_z_mesh)
  {

    log::finfo("Initializing 3D FSRs in axially extruded regions...");

    /* Determine the extent of the axial geometry 确定轴向几何形状的范围 */
    double min_z = getMinZ();
    double max_z = getMaxZ();

    auto extruded_FSRs = _extruded_FSR_keys_map.values(); // 返回所有 ExtrudedFSR* 组成的数组（new 出来的，要记得 delete）

    std::string msg = "initializing 3D FSRs";
    Progress progress(_extruded_FSR_keys_map.size(), msg, 0.1, this, true);

    int anticipated_size = 2 * _extruded_FSR_keys_map.size(); // 提前给 _FSR_keys_map（3D FSR 的 key → FSRData* 映射）预留更大的容量
    if (_overlaid_mesh != nullptr)                            // 估计：3D FSR 的数量大约是 ExtrudedFSR 数量的 2 倍（若定义了_overlaid_mesh， ×_overlaid_mesh->getNumZ()）
      anticipated_size *= _overlaid_mesh->getNumZ();

    _FSR_keys_map.realloc(anticipated_size); // 扩容 3D FSR 的哈希表。
    long total_number_fsrs_in_stack = 0;     // 累加所有挤出 FSR 内含的 3D FSR 数量，用于统计内存需求。

    /* Loop over extruded FSRs 遍历轴向挤出FSR*/
#pragma omp parallel for // 多线程并行
    for (size_t i = 0; i < _extruded_FSR_keys_map.size(); i++)
    {

      progress.incrementCounter();

      /* Extract coordinates of extruded FSR */
      auto extruded_FSR = extruded_FSRs[i];
      double x0 = extruded_FSR->_coords->getX(); // 之前在 segmentizeExtruded 时保存的“特征点”，start点
      double y0 = extruded_FSR->_coords->getY();

      /* Determine if there is a global mesh or local meshes should be created 确定是否应创建全局轴向网或局部轴向网*/
      if (global_z_mesh.size() > 0)
      { // 提供了全局轴向网（global_z_mesh）

        size_t num_regions = global_z_mesh.size() - 1; // 分割的区域数量
        extruded_FSR->_num_fsrs = num_regions;
        extruded_FSR->_materials = new Material *[num_regions]; //_materials[n]：第 n 个 3D FSR 的材料指针
        extruded_FSR->_fsr_ids = new long[num_regions];

        /* 遍历 num_regions 个 3D FSR,确定其材料和 id*/
#pragma omp critical
        {
          for (size_t n = 0; n < num_regions; n++)
          {

            /* Set the axial coordinate at the midpoint of mesh boundaries 设置网格边界中点的轴向坐标*/
            double midpt = (global_z_mesh[n] + global_z_mesh[n + 1]) / 2;

            LocalCoords coord(x0, y0, midpt); // 实例化一个局部点,以 (x0, y0, zmid) 为坐标,定在 root universe 上
            coord.setUniverse(_root_universe);

            /* Get the FSR ID and material */
            Cell *cell = findCellContainingCoords(&coord);
            long fsr_id = findFSRId(&coord); // 找到点所在 FSR 的 id(findFSRId),此时会给轴向 FSR 逐个编号
            Material *material = cell->getFillMaterial();

            /* Set the FSR ID and material */
            extruded_FSR->_fsr_ids[n] = fsr_id;
            extruded_FSR->_materials[n] = material;
          }
        }
      }
      else
      { // 没有提供全局轴向网

        // 实例化一条起点为 (x0, y0, zmin),终点为 (x0, y0, zmax) 的垂直轨迹，从几何体的最低点（min_z）到最高点（max_z）
        Track3D track;
        track.setValues(x0, y0, min_z, x0, y0, max_z, 0, 0);

        segmentize3D(&track, true); // 沿垂直轨迹的起点分段，并保存每个线段的长度及FSR_ID等数据

        /* Extract segments from track */
        int num_segments = track.getNumSegments();
        segment *segments = track.getSegments();

        /* Allocate materials and mesh in extruded FSR */
        extruded_FSR->_num_fsrs = (size_t)num_segments;
        extruded_FSR->_materials = new Material *[num_segments];
        extruded_FSR->_fsr_ids = new long[num_segments];
        extruded_FSR->_mesh = new double[num_segments + 1];

        /* Initialize values in extruded FSR */
        for (int s = 0; s < num_segments; s++)
        {
          extruded_FSR->_materials[s] = segments[s]._material;
          extruded_FSR->_fsr_ids[s] = segments[s]._region_id;
        }

        /* Initialize z mesh  初始化Z网格的高度/初始化轴向网的值
         * The accumulation performed below may cause a non-negligible
         * rounding error which makes the program never take the conditional
         * branch. 下面执行的累加可能会导致不可忽略的舍入误差，这使得程序永远不会采用条件分支（待）*/
        int s;
        double level = min_z;
        extruded_FSR->_mesh[0] = level;
        for (s = 0; s < num_segments - 1; s++)
        {
          level += segments[s]._length;
          // if (std::abs(level - max_z) < FLT_EPSILON)
          //   level = max_z;
          extruded_FSR->_mesh[s + 1] = level;
        }
        extruded_FSR->_mesh[s + 1] = max_z;
      }
/* Keep track of the number of FSRs in extruded FSRs  记录所有的轴向挤出FSR里的FSR数量（累加）*/
#pragma omp atomic update
      total_number_fsrs_in_stack += extruded_FSR->_num_fsrs;
    }

    delete[] extruded_FSRs;

    // Print the number of ExtrudedFSRs
    printNumExtrudedFSRs();

    // Print the storage requirement of ExtrudedFSRs
    printMemUsageExtrudedFSRs(total_number_fsrs_in_stack);

    /* Re-order FSR IDs so they are sequential in the axial direction */
    reorderFSRIDs(); // 调用reorderFSRIDs 把 3D FSR 重新编号，保证每个轴向挤出区的 3D FSR的 id 都是连续的,使其在轴向方向上连续

    log::info("Finished initialization of axial 3D FSRs");
  }

  /**
   * @brief Reorders FSRs so that they are contiguous in the axial direction
   * @details Only FSR ids changed. For example,
   *           13                 3
   *            9    12 11        2     8 11
   *            8 10  4  7  6  => 1  5  7 10 13
   *            0  3  1  2  5     0  4  6  9 12
   *           --------------- => --------------
   *            0  1  2  3  4     0  1  2  3  4
   *
   * _extruded_FSR_keys_map：保存了每个挤出 FSR 的数据结构（含 _fsr_ids 等）。
   * _FSR_keys_map：所有 FSR 的大哈希表，value 是 FSRData*，其中 _fsr_id 可能需要同步更新。
   */
  void Geometry::reorderFSRIDs()
  {

    log::fprofile("Reordering 3D FSR IDs...");

    /* Extract list of extruded FSRs */
    auto extruded_FSRs = _extruded_FSR_keys_map.values();

    std::string msg = "reordering FSR IDs";
    Progress progress(_extruded_FSR_keys_map.size(), msg, 0.1, this, true);

    /* Get the FSR data objects */
    // unused
    // long curr_id = 0;
    auto value_list = _FSR_keys_map.values(); // 获取 _FSR_keys_map 中所有 FSRData 指针（包含 _fsr_id、体积等）。
    // 用于修改结构体的字段
    auto fsr_data_objects = new FSRData *[_FSR_keys_map.size()]; // 准备一个数组，按照 ID 直接索引到对应的 FSRData*，后面好改 _fsr_id。

    /* Create a mapping of old to new IDs */
    // 分配两个辅助数组：id_mapping 存“老 ID -> 新 ID”的映射，id_remapped 标记某个老 ID 是否已经被分配过新号。
    long *id_mapping = new long[_FSR_keys_map.size()];
    bool *id_remapped = new bool[_FSR_keys_map.size()];
#pragma omp parallel for
    for (size_t i = 0; i < _FSR_keys_map.size(); i++)
    {
      long id = value_list[i]->_fsr_id;
      fsr_data_objects[id] = value_list[i]; // 这样后面只要凭旧 ID 就能拿到对象。
      id_mapping[i] = i;
      id_remapped[i] = false;
    }

    /* Loop over extruded FSRs */
    long count = 0;
    for (size_t i = 0; i < _extruded_FSR_keys_map.size(); i++) // 遍历每个挤出 FSR，逐层处理
    {

      progress.incrementCounter();

      /* Extract coordinates of extruded FSR */
      auto extruded_FSR = extruded_FSRs[i]; // 取得该柱体数据。

      /* Get the number of FSRs in this axial region */
      size_t num_local_fsrs = extruded_FSR->_num_fsrs; // 它在轴向上有多少段。

      /* Re-assign the IDs of all axial FSRs */
      for (size_t j = 0; j < num_local_fsrs; j++) // 对这根柱体里的每个 3D FSR：
      {

        long previous_id = extruded_FSR->_fsr_ids[j]; // 旧的 ID
        if (!id_remapped[previous_id])                // 如果这个 ID 还没分配过新号：
        {
          id_mapping[previous_id] = count; // 把它映射到当前 count，并标记为已重新编号，然后 count 自增。
          id_remapped[previous_id] = true;
          count++;
        }
        long new_id = id_mapping[previous_id];

        fsr_data_objects[previous_id]->_fsr_id = new_id; // 把 FSRData 对象里保存的 _fsr_id 更新成新号，保证全局数据一致。
        extruded_FSR->_fsr_ids[j] = new_id;
      }
    }

    delete[] extruded_FSRs;
    delete[] value_list;
    delete[] fsr_data_objects;
    delete[] id_mapping;
    delete[] id_remapped;
  }

  /**
   * @brief Initialize key and material ID vectors for lookup by FSR ID
   * @detail This function initializes and sets reverse lookup vectors by FSR ID.
   *      This is called after the FSRs have all been identified and allocated
   *      during segmentation. This function must be called after
   *      Geometry::segmentize() has completed. It should not be called if tracks
   *      are loaded from a file.
   * 在前面的分段过程中，FSR 是这样被记录的：
   * _FSR_keys_map：key → FSRData*
   * key：一个字符串，唯一描述一个 FSR 在 CSG 结构中的位置（类似“几何指纹”）；
   * FSRData 结构体里有 _fsr_id、_mat_id、_point（质心）等。
   *
   * 但是后续求解时，我们通常是通过“fsr_id”去查：
   * 这个编号对应的 key 是啥？（有时要 debug 用）
   * 他的特征点是哪个？
   * 它的材料 id 是多少？
   * 它在 CMFD 的哪个粗网格 cell 里？
   * initializeFSRVectors 就是把FSRData 结构体的数据，整理成几个“数组/向量”，让你可以通过 fsr_id O(1) 访问到对应信息。
   */
  void Geometry::initializeFSRVectors()
  {

    /* Get keys and values from map */
    log::info("Initializing FSR lookup vectors");
    auto key_list = _FSR_keys_map.keys();
    auto value_list = _FSR_keys_map.values();

    /* Allocate vectors */
    /* 分配向量 */
    size_t num_FSRs = _FSR_keys_map.size();
    _FSRs_to_keys = std::vector<std::string>(num_FSRs);           // FSR_ID与key的对应关系
    _FSRs_to_centroids = std::vector<Point *>(num_FSRs, nullptr); // FSR_ID与centroids的对应关系,不在这里赋值
    _FSRs_to_material_IDs = std::vector<int>(num_FSRs);           // FSR_ID与material_IDs的对应关系
    _FSRs_to_CMFD_cells = std::vector<int>(num_FSRs);             // FSR_ID与CMFD_cells的对应关系 它在 CMFD 的哪个粗网格 cell 里？
                                                                  // 创建四个向量，分别存储FSR ID到key、质心点、材料ID和CMFD单元的映射关系

#pragma omp parallel for
    // 是一个 OpenMP 指令，告诉编译器“把接下来的 for 循环拆成多个线程同时执行”。编译时如果启用了 OpenMP，这条指令会让循环体（通常是紧跟着的 for (...) { ... }）自动分配给多个线程，提高并行性能；如果没启用 OpenMP，它会被当作普通注释忽略，不影响串行执行。
    for (size_t i = 0; i < num_FSRs; i++) // 填充_FSRs_to_keys, _FSRs_to_material_IDs
    {
      auto key = key_list[i];
      auto fsr = value_list[i];
      long fsr_id = fsr->_fsr_id;
      _FSRs_to_keys.at(fsr_id) = key;
      _FSRs_to_material_IDs.at(fsr_id) = fsr->_mat_id;
    }

    /* Add cmfd information serially */
    /* 串行添加CMFD信息 */
    if (_cmfd != nullptr)
    {
      for (size_t i = 0; i < num_FSRs; i++)
      {
        auto fsr = value_list[i];
        auto fsr_id = fsr->_fsr_id;
        // 每个 FSRData 里有 _cmfd_cell（所在的粗网格单元）和 _point（特征点）；
        auto point = fsr->_point;
        /*
        下面两行对CMFD方法非常重要：建立FSR与CMFD单元之间的双向映射关系
        1. 将FSR添加到对应的CMFD单元中
        2. 记录每个FSR所属的CMFD单元
        std::vector<std::vector<long>> _cell_fsrs;
        这个是二维vector数组，第一维确定是具体的哪个CMFD单元，第二个是单个CMFD中所关联的所有FSR
        */
        _cmfd->addFSRToCell(fsr->_cmfd_cell, fsr_id);     // 让 CMFD 知道这个 cell 里有哪些 FSR
        _FSRs_to_CMFD_cells.at(fsr_id) = fsr->_cmfd_cell; // fsr_id与CMFD cell的对应关系
        log::fdebug("FSR %ld Point x = %.2f, y = %.2f, z = %.2f", fsr_id, point->getX(), point->getY(), point->getZ());
        log::fdebug("cmfd cell is %d, fsr id is %d", fsr->_cmfd_cell, fsr_id);
      }
      if (_cmfd->GetHexLatticeEnable()) // 00000000
        _cmfd->findEmptyCmfdCells();
    }

    //_cmfd->printCellFSRs();

    /* Check if extruded FSRs are present 检查是否存在轴向挤出FSR*/
    size_t num_extruded_FSRs = _extruded_FSR_keys_map.size();
    if (num_extruded_FSRs > 0)
    {

      /* _extruded_FSR_keys_map 是 key → ExtrudedFSR* 的 map；
      同样，我们希望通过 extruded_fsr_id 直接拿到对应的 ExtrudedFSR*/
      _extruded_FSR_lookup = std::vector<ExtrudedFSR *>(num_extruded_FSRs);
      auto extruded_value_list = _extruded_FSR_keys_map.values();

#pragma omp parallel for
      for (size_t i = 0; i < num_extruded_FSRs; i++)
      {
        long fsr_id = extruded_value_list[i]->_fsr_id;
        _extruded_FSR_lookup[fsr_id] = extruded_value_list[i];
      }

      delete[] extruded_value_list;
    }

    /* Delete key and value lists */
    delete[] key_list;
    delete[] value_list;

    // Output storage requirement of FSRs
    // 输出FSR的存储需求
    printNumFSRs();
    printMemUsageFSRs();

    log::verbose_once("Finished initialization of FSR lookup arrays");
  }

  /**
   * @brief Determines the fissionability of each Universe within this Geometry.
   * @details A Universe is determined fissionable if it contains a Cell
   *          filled by a Material with a non-zero fission cross-section. Note
   *          that this method recurses through all Universes at each level in
   *          the nested Universe hierarchy. Users should only call this method
   *          without a parameter (the default) from Python as follows to ensure
   *          that the recursion starts from the uppermost Universe level:
   *
   * @code
   *          geometry.computeFissionability()
   * @endcode
   *
   * @param univ the Universe of interest (default is nullptr)
   */
  void Geometry::computeFissionability(Universe *univ)
  {

    bool fissionable = false;

    std::map<int, Material *> materials;
    std::map<int, Universe *> universes;

    /* If no Universe was passed in as an argument, then this is the first
     * recursive call from a user via Python, so get the base Universe */
    if (univ == nullptr)
      univ = _root_universe;

    /* If a Universe was passed in as an argument, then this is a recursive
     * call with a Universe at a lower level in the nested Universe hierarchy */
    if (univ->getType() == SIMPLE)
    {
      materials = univ->getAllMaterials();
      universes = univ->getAllUniverses();
    }

    else
      universes = static_cast<Lattice *>(univ)->getAllUniverses();

    /* Loop over the nested Universes first to ensure that fissionability
     * is set at each nested Universe level */
    for (auto &u : universes)
    {
      /* Recursively check whether this nested Universe is fissionable */
      computeFissionability(u.second);

      if (u.second->isFissionable())
        fissionable = true;
    }

    /* Loop over the Materials in this Universe at this level */
    for (auto &m : materials)
    {
      /* Check whether this Material is fissionable or not */
      if (m.second->isFissionable())
        fissionable = true;
    }

    /* Set this Universe's fissionability based on the nested Universes
     * and Materials within it */
    univ->setFissionability(fissionable);
  }

  /**
   * @brief Get the material, cell or FSR IDs on a 2D spatial grid.
   * @details This is a helper method for the openmoc.plotter module.
   *          This method may also be called by the user in Python if needed.
   *          A user must initialize NumPy arrays with the x and y grid
   *          coordinates input to this function. This function then fills
   *          a NumPy array with the domain IDs for each coordinate. An example
   *          of how this function might be called in Python is as follows:
   *
   * @code
   *          grid_x = numpy.arange(-2., +2., 100)
   *          grid_y = numpy.arange(-2., +2., 100)
   *          domain_ids = geometry.getSpatialDataOnGrid(
   *              grid_x, grid_y, 20., 'xy', 'material')
   * @endcode
   *
   * @param dim1 a numpy array of the first dimension's coordinates
   * @param dim2 a numpy array of the second dimension's coordinates
   * @param offset The coordinate at which the plane is located
   * @param plane The plane for which data is gathered ('xy', 'xz', 'yz')
   * @param domain_type the type of domain ('fsr', 'material', 'cell')
   * @return a NumPy array or list of the domain IDs
   */
  std::vector<long> Geometry::getSpatialDataOnGrid(std::vector<double> dim1,
                                                   std::vector<double> dim2,
                                                   double offset,
                                                   const char *plane,
                                                   const char *domain_type)
  {

    /* Instantiate a vector to hold the domain IDs */
    std::vector<long> domains(dim1.size() * dim2.size());

    /* Extract the source region IDs */
    // #pragma omp parallel for
    for (size_t i = 0; i < dim1.size(); i++)
    {
      for (size_t j = 0; j < dim2.size(); j++)
      {

        Cell *cell;
        LocalCoords *point = nullptr;

        /* Find the Cell containing this point */
        if (strcmp(plane, "xy") == 0)
          point = new LocalCoords(dim1[i], dim2[j], offset, true);
        else if (strcmp(plane, "xz") == 0)
          point = new LocalCoords(dim1[i], offset, dim2[j], true);
        else if (strcmp(plane, "yz") == 0)
          point = new LocalCoords(offset, dim1[i], dim2[j], true);
        else
          log::ferror("Unable to extract spatial data for "
                      "unsupported plane %s",
                      plane);

        point->setUniverse(_root_universe);
        cell = _root_universe->findCell(point);
        domains[i + j * dim1.size()] = -1;

        /* Extract the ID of the domain of interest */
        if (withinGlobalBounds(point) && cell != NULL)
        {
          if (strcmp(domain_type, "fsr") == 0)
            domains[i + j * dim1.size()] = getGlobalFSRId(point, false);
          else if (strcmp(domain_type, "material") == 0)
            domains[i + j * dim1.size()] = cell->getFillMaterial()->getId();
          else if (strcmp(domain_type, "cell") == 0)
            domains[i + j * dim1.size()] = cell->getId();
          else
            log::ferror("Unable to extract spatial data for "
                        "unsupported domain type %s",
                        domain_type);
        }

        /* Deallocate memory for LocalCoords */
        point = point->getHighestLevel();
        delete point;
      }
    }

    /* Return the domain IDs */
    return domains;
  }

  /**
   * @brief Converts this Geometry's attributes to a character array.
   * @details This method calls the toString() method for all Materials,
   *          Surfaces, Cell, Universes and Lattices contained by the Geometry.
   * @return a character array of this Geometry's class attributes
   */
  std::string Geometry::toString()
  {

    std::stringstream string;

    std::map<int, Cell *> all_cells = _root_universe->getAllCells();
    std::map<int, Universe *> all_universes = _root_universe->getAllUniverses();

    string << "\n\tCells:\n\t\t";
    for (const auto &c : all_cells)
      string << c.second->toString() << "\n\t\t";

    string << "\n\tUniverses:\n\t\t";
    for (const auto &u : all_universes)
      string << u.second->toString() << "\n\t\t";

    std::string formatted_string = string.str();
    formatted_string.erase(formatted_string.end() - 3);

    return formatted_string;
  }

  /**
   * @brief Prints a string representation of all of the Geometry's attributes to
   *        the console.
   * @details This method calls the printString() method for all Materials,
   *          Surfaces, Cell, Universes and Lattices contained by the Geometry.
   */
  void Geometry::printString()
  {
    log::fresult("%s", toString().c_str());
  }

  /**
   * @brief Prints FSR layout to file
   * @details This provides a way to get the functionality of the
   *              plot_flat_source_regions Python function without Python
   * @param plane The "xy", "xz", or "yz" plane in which to extract flat source
   *        regions
   * @param gridsize The number of points to plot in each direction
   * @param offset The offset of the plane in the third dimension
   * @param bounds_x a two valued array for the plotted x-limits
   * @param bounds_y a two valued array for the plotted y-limits
   * @param bounds_z a two valued array for the plotted z-limits
   */
  void Geometry::printFSRsToFile(const char *plane, int gridsize, double offset,
                                 double *bounds_x, double *bounds_y,
                                 double *bounds_z)
  {

    /* Get geometry min and max */
    double min_x = _root_universe->getMinX();
    double max_x = _root_universe->getMaxX();
    double min_y = _root_universe->getMinY();
    double max_y = _root_universe->getMaxY();
    double min_z = _root_universe->getMinZ();
    double max_z = _root_universe->getMaxZ();

    if (bounds_x != nullptr)
    {
      min_x = bounds_x[0];
      max_x = bounds_x[1];
    }
    if (bounds_y != nullptr)
    {
      min_y = bounds_y[0];
      max_y = bounds_y[1];
    }
    if (bounds_z != nullptr)
    {
      min_z = bounds_z[0];
      max_z = bounds_z[1];
    }

    /* Create coordinate vectors */
    std::vector<double> dim1(gridsize);
    std::vector<double> dim2(gridsize);

    /* Determine minimum and maximum values */
    double dim1_min = -1;
    double dim1_max = -1;
    double dim2_min = -1;
    double dim2_max = -1;

    /* x-y plane */
    if (strcmp(plane, "xy") == 0)
    {
      dim1_min = min_x;
      dim1_max = max_x;
      dim2_min = min_y;
      dim2_max = max_y;
    }

    /* x-z plane */
    else if (strcmp(plane, "xz") == 0)
    {
      dim1_min = min_x;
      dim1_max = max_x;
      dim2_min = min_z;
      dim2_max = max_z;
    }

    /* y-z plane */
    else if (strcmp(plane, "yz") == 0)
    {
      dim1_min = min_y;
      dim1_max = max_y;
      dim2_min = min_z;
      dim2_max = max_z;
    }

    else
    {
      log::ferror("Plane type %s unrecognized", plane);
    }

    /* Create grid */
    double width1 = (dim1_max - dim1_min) / (gridsize + 1);
    double width2 = (dim2_max - dim2_min) / (gridsize + 1);
    for (int i = 0; i < gridsize; i++)
    {
      dim1.at(i) = dim1_min + (i + 1) * width1;
      dim2.at(i) = dim2_min + (i + 1) * width2;
    }

    /* Retrieve data */
    log::finfo("Getting FSR layout on domains");
    std::vector<long> domain_data = getSpatialDataOnGrid(dim1, dim2, offset,
                                                         plane, "fsr");

    long *fsr_array = new long[domain_data.size()];
#pragma omp parallel for
    for (size_t i = 0; i < domain_data.size(); i++)
    {
      fsr_array[i] = domain_data.at(i) + 1;
    }

#ifdef ENABLE_MPI_
    if (isDomainDecomposed())
    {
      log::finfo("Communicating FSR layout accross domains");
      long *reduced_fsr_array = new long[domain_data.size()];
      MPI_Allreduce(fsr_array, reduced_fsr_array, domain_data.size(),
                    mpi::getDatatype<long>(), MPI_SUM, getMPICart());
#pragma omp parallel for
      for (size_t i = 0; i < domain_data.size(); i++)
        fsr_array[i] = reduced_fsr_array[i];
      delete[] reduced_fsr_array;
    }
#endif

    /* Print to file */
    log::finfo("Printing FSRs to file");
    if (isRootDomain())
    {
      std::ofstream out("fsr-printout.txt");
      out << "[HEADER] FSR printout" << std::endl;
      out << "[HEADER] Plane = " << plane << std::endl;
      out << "[HEADER] Offset = " << offset << std::endl;
      out << "[HEADER] Bounds = (" << dim1_min << ", " << dim1_max << ") x ("
          << dim2_min << ", " << dim2_max << ")" << std::endl;
      out << "[HEADER] Gridsize = " << gridsize << std::endl;
      for (size_t i = 0; i < domain_data.size(); i++)
      {
        out << fsr_array[i] << " ";
      }
    }

    delete[] fsr_array;
  }

  /**
   * @brief This is a method that initializes the CMFD Lattice and sets
   *          CMFD parameters.
   * 初始化CMFD网格并设置相关参数
   */
  void Geometry::initializeCmfd()
  {

    /* Get the global Geometry boundary conditions
    获取几何模型六个面的边界条件
    */
    boundaryType min_x_bound = _root_universe->getMinXBoundaryType(); // 计算方法主要是calculateBoundaries();
    boundaryType max_x_bound = _root_universe->getMaxXBoundaryType();
    boundaryType min_y_bound = _root_universe->getMinYBoundaryType();
    boundaryType max_y_bound = _root_universe->getMaxYBoundaryType();
    boundaryType min_z_bound = _root_universe->getMinZBoundaryType();
    boundaryType max_z_bound = _root_universe->getMaxZBoundaryType();

    /* Get the global Geometry boundaries
    获取几何模型六个边界的坐标值
    */
    double min_x = _root_universe->getMinX();
    double max_x = _root_universe->getMaxX();
    double min_y = _root_universe->getMinY();
    double max_y = _root_universe->getMaxY();
    double min_z = _root_universe->getMinZ();
    double max_z = _root_universe->getMaxZ();

    if (_cmfd->GetHexLatticeEnable()) // 这部分是HexLattice的CMFD初始化代码
    {                                 // HexLattice   000000000000
      /* Set CMFD mesh boundary conditions */

      _cmfd->setBoundary(HEX_SURFACE_BETA_MIN, min_x_bound);
      _cmfd->setBoundary(HEX_SURFACE_DELTA_MIN, min_y_bound);
      _cmfd->setBoundary(HEX_SURFACE_GAMMA_MIN, max_x_bound);
      _cmfd->setBoundary(HEX_SURFACE_Z_MIN, min_z_bound);
      _cmfd->setBoundary(HEX_SURFACE_BETA_MAX, max_x_bound);
      _cmfd->setBoundary(HEX_SURFACE_DELTA_MAX, max_y_bound);
      _cmfd->setBoundary(HEX_SURFACE_GAMMA_MAX, max_x_bound);
      _cmfd->setBoundary(HEX_SURFACE_Z_MAX, max_z_bound);

      /* Set CMFD mesh dimensions */
      _cmfd->setWidthX(max_x - min_x);
      _cmfd->setWidthY(max_y - min_y);
      _cmfd->setWidthZ(max_z - min_z);

      /* Initialize the CMFD lattice */
      Point offset;
      offset.setX(min_x + (max_x - min_x) / 2.0);
      offset.setY(min_y + (max_y - min_y) / 2.0);

      if (std::abs(min_z + (max_z - min_z) / 2.0) < FLT_INFINITY)
        offset.setZ(min_z + (max_z - min_z) / 2.0);
      else
        offset.setZ(0.);

      log::fdebug("CMFD Min x = %.2f, Max x = %.2f, Min y = %.2f, Max y = %.2f, Min z = %.2f, Max z = %.2f",
                  min_x, max_x, min_y, max_y, min_z, max_z);
      log::fdebug("CMFD Offset x = %.2f, y = %.2f ,z = %.2f", offset.getX(), offset.getY(), offset.getZ());

      _cmfd->setGeometry(this);
      _cmfd->initializeLattice(&offset);

#ifdef ENABLE_MPI_
      if (isDomainDecomposed())
      {
        /* Check that CMFD mesh is compatible with domain decomposition */
        _cmfd->setNumDomains(mpi::getNumDomainsX(),
                             mpi::getNumDomainsY(),
                             mpi::getNumDomainsZ());
        _cmfd->setDomainIndexes(mpi::getDomainIndexX(),
                                mpi::getDomainIndexY(),
                                mpi::getDomainIndexZ());
      }
#endif
      /* Initialize CMFD Maps -. 初始化每个 CMFD 单元与 FSR 关联的二维向量*/
      _cmfd->initializeCellMap();
    }
    else
    { // RecLattice - 四边形
      // Set CMFD mesh boundary conditions
      _cmfd->setBoundary(SURFACE_X_MIN, min_x_bound);
      _cmfd->setBoundary(SURFACE_Y_MIN, min_y_bound);
      _cmfd->setBoundary(SURFACE_Z_MIN, min_z_bound);
      _cmfd->setBoundary(SURFACE_X_MAX, max_x_bound);
      _cmfd->setBoundary(SURFACE_Y_MAX, max_y_bound);
      _cmfd->setBoundary(SURFACE_Z_MAX, max_z_bound);

      /* Set CMFD mesh dimensions 设置CMFD网格在三个方向的总宽度，把universe几何体的物理尺寸传递给CMFD */
      _cmfd->setWidthX(max_x - min_x);
      _cmfd->setWidthY(max_y - min_y);
      _cmfd->setWidthZ(max_z - min_z);

      /* Initialize the CMFD lattice 初始化 CMFD lattice */
      Point offset; // 选择几何中心作为参考点
      offset.setX(min_x + (max_x - min_x) / 2.0);
      offset.setY(min_y + (max_y - min_y) / 2.0);
      if (std::abs(min_z + (max_z - min_z) / 2.0) < FLT_INFINITY) // 如果Z中心坐标有效（不是无穷大），使用计算值；否则设为0
        offset.setZ(min_z + (max_z - min_z) / 2.0);
      else
        offset.setZ(0.);

      _cmfd->setGeometry(this);

      _cmfd->initializeLattice(&offset); //  设置表示整个CMFD的lattice的沿各方向的单元宽度以及累积宽度，并设置此Lattice的几何中心在父Universe坐标系中的位置

#ifdef ENABLE_MPI_
      if (isDomainDecomposed()) // 检查是否进行了域分解，域分解: 将计算区域分割给多个进程并行计算
      {
        /*
         mpi::getNumDomainsX(): 获取X方向的域数量（MPI进程数）
         mpi::getDomainIndexX():  获取当前进程负责的子区域在X方向的索引
         */
        _cmfd->setNumDomains(mpi::getNumDomainsX(),
                             mpi::getNumDomainsY(),
                             mpi::getNumDomainsZ());

        _cmfd->setDomainIndexes(mpi::getDomainIndexX(),
                                mpi::getDomainIndexY(),
                                mpi::getDomainIndexZ());
      }
#endif
      /* Initialize CMFD Maps  只是给一个域中的所有CMFD单元对应的FSR数组分配内存*/
      _cmfd->initializeCellMap();
    }
  }

  /**
   * @brief This is a method that initializes the initial spectrum calculator
   */
  void Geometry::initializeSpectrumCalculator(Cmfd *spectrum_calculator)
  {

    /* Setup the CMFD lattice with the domain dimensions */
    spectrum_calculator->setLatticeStructure(mpi::getNumDomainsX(),
                                             mpi::getNumDomainsY(),
                                             mpi::getNumDomainsZ());

    /* Get the global Geometry boundary conditions */
    boundaryType min_x_bound = _root_universe->getMinXBoundaryType();
    boundaryType max_x_bound = _root_universe->getMaxXBoundaryType();
    boundaryType min_y_bound = _root_universe->getMinYBoundaryType();
    boundaryType max_y_bound = _root_universe->getMaxYBoundaryType();
    boundaryType min_z_bound = _root_universe->getMinZBoundaryType();
    boundaryType max_z_bound = _root_universe->getMaxZBoundaryType();

    /* Get the global Geometry boundaries */
    double min_x = _root_universe->getMinX();
    double max_x = _root_universe->getMaxX();
    double min_y = _root_universe->getMinY();
    double max_y = _root_universe->getMaxY();
    double min_z = _root_universe->getMinZ();
    double max_z = _root_universe->getMaxZ();

    /* Set spectrum caclulator boundary conditions */
    spectrum_calculator->setBoundary(SURFACE_X_MIN, min_x_bound);
    spectrum_calculator->setBoundary(SURFACE_Y_MIN, min_y_bound);
    spectrum_calculator->setBoundary(SURFACE_Z_MIN, min_z_bound);
    spectrum_calculator->setBoundary(SURFACE_X_MAX, max_x_bound);
    spectrum_calculator->setBoundary(SURFACE_Y_MAX, max_y_bound);
    spectrum_calculator->setBoundary(SURFACE_Z_MAX, max_z_bound);

    /* Set spectrum calculator dimensions */
    spectrum_calculator->setWidthX(max_x - min_x);
    spectrum_calculator->setWidthY(max_y - min_y);
    spectrum_calculator->setWidthZ(max_z - min_z);

    /* Initialize CMFD Maps */
    spectrum_calculator->initializeCellMap();

    /* Initialize the CMFD lattice */
    Point offset;
    offset.setX(min_x + (max_x - min_x) / 2.0);
    offset.setY(min_y + (max_y - min_y) / 2.0);
    offset.setZ(min_z + (max_z - min_z) / 2.0);

    spectrum_calculator->initializeLattice(&offset);
    spectrum_calculator->setGeometry(this);

#ifdef ENABLE_MPI_
    if (isDomainDecomposed())
    {
      spectrum_calculator->setNumDomains(mpi::getNumDomainsX(),
                                         mpi::getNumDomainsY(),
                                         mpi::getNumDomainsZ());
      spectrum_calculator->setDomainIndexes(mpi::getDomainIndexX(),
                                            mpi::getDomainIndexY(),
                                            mpi::getDomainIndexZ());
    }
#endif

    /* Add FSRs to domain cell */
    for (long r = 0; r < getNumFSRs(); r++)
      spectrum_calculator->addFSRToCell(0, r);
  }

  /**
   * @brief Returns a pointer to the map that maps FSR keys to FSR IDs
   * @return pointer to _FSR_keys_map map of FSR keys to FSR IDs
   */
  ParallelHashMap<std::string, FSRData *> &Geometry::getFSRKeysMap()
  {
    return _FSR_keys_map;
  }

  /**
   * @brief Returns a pointer to the map that maps FSR keys to extruded FSRs
   * @return pointer to _FSR_keys_map map of FSR keys to extruded FSRs
   */
  ParallelHashMap<std::string, ExtrudedFSR *> &Geometry::getExtrudedFSRKeysMap()
  {
    return _extruded_FSR_keys_map;
  }

  /**
   * @brief Returns the vector that maps FSR IDs to FSR key hashes
   * @return _FSR_keys_map map of FSR keys to FSR IDs
   */
  std::vector<std::string> &Geometry::getFSRsToKeys()
  {
    return _FSRs_to_keys;
  }

  /**
   * @brief Returns the vector that maps FSR IDs to extruded FSRs
   * @return _extruded_FSR_lookup map of FSR keys to extruded FSRs
   */
  std::vector<ExtrudedFSR *> &Geometry::getExtrudedFSRLookup()
  {
    return _extruded_FSR_lookup;
  }

  /**
   * @brief Return a vector indexed by flat source region IDs which contains
   *        the corresponding Material IDs.
   * @return an integer vector of FSR-to-Material IDs indexed by FSR ID
   */
  std::vector<int> &Geometry::getFSRsToMaterialIDs()
  {
    return _FSRs_to_material_IDs;
  }

  /**
   * @brief Return a vector indexed by flat source region IDs which contains
   *        pointers to the corresponding Centroid Information.
   * @return an array of centroid pointers indexed by FSR ID
   */
  std::vector<Point *> &Geometry::getFSRsToCentroids()
  {
    return _FSRs_to_centroids;
  }

  /**
   * @brief Return a vector indexed by flat source region IDs which contains
   *        the corresponding CMFD cell.
   * @return an integer vector of FSR to CMFD cell IDs indexed by FSR ID
   */
  std::vector<int> &Geometry::getFSRsToCMFDCells()
  {
    return _FSRs_to_CMFD_cells;
  }

  /// \brief   Returns a map of Material IDs to indices in the array of materials
  /// \details This method is basically used by GPU-related routines. The indices
  ///          is starting from 0 and is contiguous. With this map, tracks on GPU
  ///          could reference contiguous IDs of materials rather than hold
  ///          pointers to materials.
  /// \return A map of Material IDs to indices.
  std::map<int, int> Geometry::getMaterialIDsToIndices()
  {
    auto &materials = getAllMaterials();

    if (materials.empty())
    {
      log::error("Map of materials is empty, which should not be used to "
                 "generate the ID-to-indice map");
    }

    int material_index = 0;
    std::map<int, int> material_IDs_to_indices;

    for (const auto &m : materials)
      material_IDs_to_indices[m.second->getId()] = material_index++;

    return material_IDs_to_indices;
  }

  /**
   * @brief Determines whether a point is within the bounding box of the domain.
   * @param coords a populated LocalCoords linked list
   * @return boolean indicating whether the coords is within the domain
   */
  bool Geometry::withinBounds(LocalCoords *coords)
  {

    double x = coords->getX();
    double y = coords->getY();
    double z = coords->getZ();

    if (x <= getMinX() - FLT_EPSILON || x >= getMaxX() + FLT_EPSILON ||
        y <= getMinY() - FLT_EPSILON || y >= getMaxY() + FLT_EPSILON ||
        z <= getMinZ() - FLT_EPSILON || z >= getMaxZ() + FLT_EPSILON)
      return false;
    else
      return true;
  }

  /**
   * @brief Determines whether a point is within the bounding box of the Geometry.
   * @param coords a populated LocalCoords linked list
   * @return boolean indicating whether the coords is within the geometry
   */
  bool Geometry::withinGlobalBounds(LocalCoords *coords)
  {

    double x = coords->getX();
    double y = coords->getY();
    double z = coords->getZ();

    if (x <= _root_universe->getMinX() || x >= _root_universe->getMaxX() ||
        y <= _root_universe->getMinY() || y >= _root_universe->getMaxY() ||
        z <= _root_universe->getMinZ() || z >= _root_universe->getMaxZ())
      return false;
    else
      return true;
  }

  /**
   * @brief Finds the Cell containing a given fsr ID.
   * @param fsr_id an FSR ID.
   */
  Cell *Geometry::findCellContainingFSR(long fsr_id)
  {

    std::string &key = _FSRs_to_keys[fsr_id];
    Point *point = _FSR_keys_map.at(key)->_point;
    LocalCoords *coords = new LocalCoords(point->getX(), point->getY(),
                                          point->getZ(), true);
    coords->setUniverse(_root_universe);
    Cell *cell = findCellContainingCoords(coords);

    delete coords;

    return cell;
  }

  /**
   * @brief Sets the centroid for an FSR
   * @details The _FSR_keys_map stores a hash of a std::string representing
   *          the Lattice/Cell/Universe hierarchy for a unique region
   *          and the associated FSR data. _centroid is a point that represents
   *          the numerical centroid of an FSR computed using all segments
   *          contained in the FSR. This method is used by the TrackGenerator
   *          to set the centroid after segments have been created. It is
   *          important to note that this method is a helper function for the
   *          TrackGenerator and should not be explicitly called by the user.
   * @param fsr a FSR id
   * @param centroid a Point representing the FSR centroid
   */
  void Geometry::setFSRCentroid(long fsr, Point *centroid)
  {
    _contains_FSR_centroids = true;
    std::string &key = _FSRs_to_keys[fsr];
    _FSR_keys_map.at(key)->_centroid = centroid;
    _FSRs_to_centroids[fsr] = centroid;
  }

  /**
   * @brief Returns a vector of z-coords defining a superposition of all axial
   *        boundaries in the Geometry.
   * @details The Geometry is traversed to retrieve all Z-planes and implicit
   *          z-boundaries, such as lattice boundaries. The levels of all these
   *          z-boundaries are rounded and added to a set containing no
   *          duplicates, creating a mesh.
   * @param include_overlaid_mesh whether to include an overlaid mesh in the
   *        set of unique z-coords
   * @return a vector of z-coords
   */
  std::vector<double> Geometry::getUniqueZHeights(bool include_overlaid_mesh)
  {

    /* Get the bounds of the geometry */
    double min_z = getMinZ();
    double max_z = getMaxZ();

    /* Initialize set for axial mesh */
    std::set<double> unique_mesh;
    unique_mesh.insert(min_z);
    unique_mesh.insert(max_z);

    /* Initialize vector of unvisited universes and add the root universe */
    std::vector<Universe *> universes;
    universes.push_back(_root_universe);

    /* Initialize vector of offsets */
    std::vector<double> offsets;
    offsets.push_back(0.0);

    /* Cycle through known universes */
    while (!universes.empty())
    {

      /* Get the last universe and explore it */
      Universe *curr_universe = universes.back();
      universes.pop_back();

      /* Get the z-offset of the universe */
      // 最外层中心到当前层中心的z偏移量
      double z_offset = offsets.back();
      offsets.pop_back();

      /* Store a vector of the z_heights before rounding */
      // 保存这一层的z高度
      std::vector<double> z_heights;

      /* Check if universe is actually a lattice */
      universeType type = curr_universe->getType();
      if (type == LATTICE)
      {

        /* Get lattice dimensions */
        RecLattice *lattice = static_cast<RecLattice *>(curr_universe);

        /* Calculate z-intersections */
        const std::vector<double> &accumulatez = lattice->getAccumulateZ();
        // 最外层中心到内层底边的z偏移量
        double offset = z_offset + lattice->getMinZ();
        for (int k = 0; k < lattice->getNumZ() + 1; k++)
        {
          // 每层lattice cell的高度都算上
          double z_height = accumulatez[k] + offset;
          z_heights.push_back(z_height);
        }

        /* Add universes to unvisited universes vector */
        for (auto &u : *lattice)
        {
          universes.push_back(u);
          offsets.push_back(z_offset);
        }
      }

      /* Otherwise check if universe is simple, contains cells */
      else if (type == SIMPLE)
      {

        /* Get all cells in the universe */
        std::map<int, Cell *> cells = curr_universe->getCells();

        /* Cycle through all cells */
        for (const auto &c : cells)
        {

          /* Get surfaces bounding the cell */
          SurfaceMap surfaces =
              c.second->getSurfaces();

          /* Cycle through all surfaces and add them to the set */
          for (const auto &s : surfaces)
          {

            /* Extract surface type */
            Surface *surface = s.second;
            surfaceType surf_type = surface->getSurfaceType();

            /* Treat surface types */
            if (surf_type == PLANE)
            {

              /* Extract plane paramters */
              Plane *plane = static_cast<Plane *>(surface);
              double A = plane->getA();
              double B = plane->getB();
              double C = plane->getC();
              double D = plane->getD();

              /* Check if there is a z-component */
              if (fabs(C) > FLT_EPSILON)
              {

                /* Check if plane has a continuous varying slope */
                if (fabs(A) > FLT_EPSILON || fabs(B) > FLT_EPSILON)
                  log::ferror("Continuous axial variation found in the "
                              "Geometry during axial on-the-fly ray tracing. "
                              "Axial on-the-fly ray tracing only supports "
                              "geometries that are capable of having an axially "
                              "extruded representation");

                /* Otherwise, surface is a z-plane */
                // 此时是一个表达式为Cz+D=0的z平面
                else
                  z_heights.push_back(-D / C + z_offset);
              }
            }

            /* Treat explicit z-planes */
            else if (surf_type == ZPLANE)
            {
              ZPlane *zplane = static_cast<ZPlane *>(surface);
              z_heights.push_back(zplane->getZ() + z_offset);
            }
          }

          /* Add min and max z-height to the cell */
          double z_limits[2];
          z_limits[0] = c.second->getMinZ();
          z_limits[1] = c.second->getMaxZ();
          for (int i = 0; i < 2; i++)
          {
            if (std::abs(z_limits[i]) != std::numeric_limits<double>::infinity())
              z_heights.push_back(z_limits[i] + z_offset);
          }

          /* See if cell is filled with universes or lattices */
          cellType cell_type = c.second->getType();
          if (cell_type == FILL)
          {
            Universe *new_universe = c.second->getFillUniverse();
            universes.push_back(new_universe);
            offsets.push_back(z_offset);
          }
        }
      }

      /* Add rounded z-heights to the set of unique z-heights */
      for (size_t i = 0; i < z_heights.size(); i++)
      {

        /* Round z-height */
        // 如果z_heights不是整数，会变小
        z_heights[i] = floor(z_heights[i] / ON_SURFACE_THRESH) * ON_SURFACE_THRESH;

        /* Add the rounded z-height to the set */
        if (z_heights[i] > min_z && z_heights[i] < max_z)
          unique_mesh.insert(z_heights[i]);
      }
    }

    /* Include overlaid mesh heights if requested */
    if (include_overlaid_mesh && _overlaid_mesh != nullptr)
    {
      int num_z = _overlaid_mesh->getNumZ();
      double dz = (max_z - min_z) / num_z;
      for (int i = 1; i < num_z; i++)
        unique_mesh.insert(min_z + i * dz);
    }

    /* Get a vector of the unique z-heights in the Geometry */
    std::vector<double> unique_heights;
    std::set<double>::iterator iter;
    for (iter = unique_mesh.begin(); iter != unique_mesh.end(); ++iter)
      unique_heights.push_back(static_cast<double>(*iter));

    std::sort(unique_heights.begin(), unique_heights.end());

    /* Output the unique Z heights for debugging */
    log::verbose("Found {} unique Z heights with bounds (cm): {}",
                 unique_heights.size(), stringutils::join(unique_heights, ", "));

    return unique_heights;
  }

  /**
   * @brief Returns a vector of z-coords defining potential unique radial planes
   *        in the Geometry
   * @details The Geometry is traversed to retrieve all Z-planes and implicit
   *          z-boundaries, such as lattice boundaries. The mid points of this
   *          mesh are then used to construct a vector of all potential unique
   *          radial planes and returned to the user.
   * @return a vector of z-coords
   *
   * 函数先拿到几何里所有不同的 Z 坐标集合，然后返回相邻节点之间的中点位置，避免与真实边界重合，选取代表性的“径向平面”
   */
  std::vector<double> Geometry::getUniqueZPlanes()
  {

    /* Get a vector of all unique z-heights in the Geometry */
    std::vector<double> unique_heights = getUniqueZHeights(); // 先调用另一个工具函数，把几何里所有不同的 Z 坐标集合取出来。这个向量应该是已排序且去重的。

    /* Use the midpoints to construct all possible unique radial planes */
    std::vector<double> unique_z_planes;
    for (size_t i = 1; i < unique_heights.size(); i++)
    {
      double mid = (unique_heights[i - 1] + unique_heights[i]) / 2;
      unique_z_planes.push_back(mid);
    }

    /* Output the unique Z planes for debugging */
    log::verbose("Found {} unique Z planes (cm): {}",
                 unique_z_planes.size(), stringutils::join(unique_z_planes, ", "));

    return unique_z_planes;
  }

  /**
   * @brief Prints all Geometry and Material details to a Geometry restart file
   * @param filename The name of the file where the data is printed
   */
  void Geometry::dumpToFile(std::string filename)
  {

    FILE *out;
    out = fopen(filename.c_str(), "w");

    /* Print number of energy groups */
    int num_groups = getNumEnergyGroups();
    fwrite(&num_groups, sizeof(int), 1, out);

    /* Print all material information */
    auto &all_materials = getAllMaterials();
    int num_materials = all_materials.size();
    fwrite(&num_materials, sizeof(int), 1, out);
    for (const auto &m : all_materials)
    {
      int key = m.first;
      Material *mat = m.second;
      int id = mat->getId();
      char *name = mat->getName();

      /* Print key and general material information */
      fwrite(&key, sizeof(int), 1, out);
      fwrite(&id, sizeof(int), 1, out);
      if (strcmp(name, "") == 0)
      {
        int length = 0;
        fwrite(&length, sizeof(int), 1, out);
      }
      else
      {
        int length = std::char_traits<char>::length(name);
        fwrite(&length, sizeof(int), 1, out);
        fwrite(name, sizeof(char), length, out);
      }

      /* Print total cross-section */
      for (int g = 0; g < num_groups; g++)
      {
        double value = mat->getSigmaTByGroup(g + 1);
        fwrite(&value, sizeof(double), 1, out);
      }

      /* Print fission cross-section */
      for (int g = 0; g < num_groups; g++)
      {
        double value = mat->getSigmaFByGroup(g + 1);
        fwrite(&value, sizeof(double), 1, out);
      }

      /* Print nu * fisison cross-section */
      for (int g = 0; g < num_groups; g++)
      {
        double value = mat->getNuSigmaFByGroup(g + 1);
        fwrite(&value, sizeof(double), 1, out);
      }

      /* Print neutron emission spectrum (chi) */
      for (int g = 0; g < num_groups; g++)
      {
        double value = mat->getChiByGroup(g + 1);
        fwrite(&value, sizeof(double), 1, out);
      }

      /* Print absorption cross section */
      // NOTE This can be used to transfer another XS, like U238 absorption
      FP_PRECISION *xs = mat->getSigmaA();
      for (int g = 0; g < num_groups; g++)
      {
        double value = xs[g];
        fwrite(&value, sizeof(double), 1, out);
      }

      /* Print scattering cross-section */
      for (int g = 0; g < num_groups; g++)
      {
        for (int gp = 0; gp < num_groups; gp++)
        {
          double value = mat->getSigmaSByGroup(g + 1, gp + 1);
          fwrite(&value, sizeof(double), 1, out);
        }
      }
    }

    /* Print root universe ID */
    int root_id = _root_universe->getId();
    fwrite(&root_id, sizeof(int), 1, out);

    /* Retrieve all surfaces, cells, and universes */
    std::map<int, Surface *> all_surfaces = getAllSurfaces();
    std::map<int, Cell *> all_cells = _root_universe->getAllCells();
    std::map<int, Universe *> all_universes = _root_universe->getAllUniverses();

    /* Print all surface information */
    int num_surfaces = all_surfaces.size();
    fwrite(&num_surfaces, sizeof(int), 1, out);
    for (const auto &s : all_surfaces)
    {

      /* Get key, value pair and general surface information */
      int key = s.first;
      Surface *value = s.second;
      int id = value->getId();
      char *name = value->getName();
      surfaceType st = value->getSurfaceType();
      boundaryType bt = value->getBoundaryType();

      /* Print key and general surface information */
      fwrite(&key, sizeof(int), 1, out);
      fwrite(&id, sizeof(int), 1, out);
      if (strcmp(name, "") == 0)
      {
        int length = 0;
        fwrite(&length, sizeof(int), 1, out);
      }
      else
      {
        int length = std::char_traits<char>::length(name);
        fwrite(&length, sizeof(int), 1, out);
        fwrite(name, sizeof(char), length, out);
      }
      fwrite(&st, sizeof(surfaceType), 1, out);
      fwrite(&bt, sizeof(boundaryType), 1, out);

      /* Treat specific surface types */
      if (st == PLANE)
      {
        Plane *pl = static_cast<Plane *>(value);
        double a = pl->getA();
        double b = pl->getB();
        double c = pl->getC();
        double d = pl->getD();
        fwrite(&a, sizeof(double), 1, out);
        fwrite(&b, sizeof(double), 1, out);
        fwrite(&c, sizeof(double), 1, out);
        fwrite(&d, sizeof(double), 1, out);
      }
      else if (st == ZCYLINDER)
      {
        ZCylinder *zcyl = static_cast<ZCylinder *>(value);
        double x = zcyl->getX0();
        double y = zcyl->getY0();
        double radius = zcyl->getRadius();
        fwrite(&x, sizeof(double), 1, out);
        fwrite(&y, sizeof(double), 1, out);
        fwrite(&radius, sizeof(double), 1, out);
      }
      else if (st == XPLANE)
      {
        XPlane *xpl = static_cast<XPlane *>(value);
        double x = xpl->getX();
        fwrite(&x, sizeof(double), 1, out);
      }
      else if (st == YPLANE)
      {
        YPlane *ypl = static_cast<YPlane *>(value);
        double y = ypl->getY();
        fwrite(&y, sizeof(double), 1, out);
      }
      else if (st == ZPLANE)
      {
        ZPlane *zpl = static_cast<ZPlane *>(value);
        double z = zpl->getZ();
        fwrite(&z, sizeof(double), 1, out);
      }
      else
      {
        log::ferror("Unsupported surface type for surface ID: %d, name:",
                    " %s", id, name);
      }
    }

    /* Print all cell information */
    int num_cells = all_cells.size();
    fwrite(&num_cells, sizeof(int), 1, out);
    for (const auto &c : all_cells)
    {

      /* Get key, value pair and general cell information */
      int key = c.first;
      Cell *cell = c.second;
      int id = cell->getId();
      char *name = cell->getName();
      cellType ct = cell->getType();

      /* Print key and general cell information */
      fwrite(&key, sizeof(int), 1, out);
      fwrite(&id, sizeof(int), 1, out);
      if (strcmp(name, "") == 0)
      {
        int length = 0;
        fwrite(&length, sizeof(int), 1, out);
      }
      else
      {
        int length = std::char_traits<char>::length(name);
        fwrite(&length, sizeof(int), 1, out);
        fwrite(name, sizeof(char), length, out);
      }
      fwrite(&ct, sizeof(cellType), 1, out);

      /* Print cell fill information */
      if (ct == MATERIAL)
      {
        Material *mat = cell->getFillMaterial();
        int mat_id = mat->getId();
        fwrite(&mat_id, sizeof(int), 1, out);
      }
      else if (ct == FILL)
      {
        Universe *univ = cell->getFillUniverse();
        int univ_id = univ->getId();
        fwrite(&univ_id, sizeof(int), 1, out);
      }

      /* Print cell rotations */
      bool rot = cell->isRotated();
      fwrite(&rot, sizeof(bool), 1, out);
      if (rot)
      {
        double rotation[3];
        rotation[0] = cell->getPhi("radians");
        rotation[1] = cell->getTheta("radians");
        rotation[2] = cell->getPsi("radians");
        fwrite(rotation, sizeof(double), 3, out);
      }

      /* Print cell translations */
      bool trans = cell->isTranslated();
      fwrite(&trans, sizeof(bool), 1, out);
      if (trans)
      {
        double *translation = cell->getTranslation();
        fwrite(translation, sizeof(double), 3, out);
      }

      /* Print ring / sector information */
      int num_rings = cell->getNumRings();
      int num_sectors = cell->getNumSectors();
      fwrite(&num_rings, sizeof(int), 1, out);
      fwrite(&num_sectors, sizeof(int), 1, out);

      /* Print parent cell */
      bool has_parent = cell->hasParent();
      fwrite(&has_parent, sizeof(bool), 1, out);
      if (has_parent)
      {
        int parent_id = cell->getParent()->getId();
        fwrite(&parent_id, sizeof(int), 1, out);
      }

      /* Print bounding surfaces */
      // FIXME, print regions
      int num_cell_surfaces = cell->getNumSurfaces();
      fwrite(&num_cell_surfaces, sizeof(int), 1, out);
      for (const auto &hs : *cell)
      { // cell iterator
        int surface_id = hs.surface->getId();
        int halfspace = hs.halfspace;
        fwrite(&surface_id, sizeof(int), 1, out);
        fwrite(&halfspace, sizeof(int), 1, out);
      }

      // FIXME Print neighbors or decide to re-compute them
    }

    /* Print all universe information */
    int num_universes = all_universes.size();
    fwrite(&num_universes, sizeof(int), 1, out);
    for (const auto &u : all_universes)
    {

      /* Get key, value pair and general universe information */
      int key = u.first;
      Universe *universe = u.second;
      int id = universe->getId();
      char *name = universe->getName();
      universeType ut = universe->getType();

      /* Print key and general universe information */
      fwrite(&key, sizeof(int), 1, out);
      fwrite(&id, sizeof(int), 1, out);
      if (strcmp(name, "") == 0)
      {
        int length = 0;
        fwrite(&length, sizeof(int), 1, out);
      }
      else
      {
        int length = std::char_traits<char>::length(name);
        fwrite(&length, sizeof(int), 1, out);
        fwrite(name, sizeof(char), length, out);
      }
      fwrite(&ut, sizeof(universeType), 1, out);

      if (ut == SIMPLE)
      {
        /* Print all cells in the universe */
        std::map<int, Cell *> cells = universe->getCells();
        int num_universe_cells = cells.size();
        fwrite(&num_universe_cells, sizeof(int), 1, out);
        for (const auto &c : cells)
        {
          int cell_id = c.first;
          fwrite(&cell_id, sizeof(int), 1, out);
        }
      }
      else if (ut == LATTICE)
      {
        /* Print lattice information */
        RecLattice *lattice = static_cast<RecLattice *>(universe);
        int num_x = lattice->getNumX();
        int num_y = lattice->getNumY();
        int num_z = lattice->getNumZ();
        double width_x = lattice->getWidthX();
        double width_y = lattice->getWidthY();
        double width_z = lattice->getWidthZ();
        double *offset = lattice->getOffset()->getXYZ();
        fwrite(&num_x, sizeof(int), 1, out);
        fwrite(&num_y, sizeof(int), 1, out);
        fwrite(&num_z, sizeof(int), 1, out);
        fwrite(&width_x, sizeof(double), 1, out);
        fwrite(&width_y, sizeof(double), 1, out);
        fwrite(&width_z, sizeof(double), 1, out);
        fwrite(offset, sizeof(double), 3, out);
        bool non_uniform = lattice->getNonUniform();
        fwrite(&non_uniform, sizeof(bool), 1, out);
        if (non_uniform)
        {
          const std::vector<double> widths_x = lattice->getWidthsX();
          const std::vector<double> widths_y = lattice->getWidthsY();
          const std::vector<double> widths_z = lattice->getWidthsZ();
          fwrite(&widths_x[0], sizeof(double), num_x, out);
          fwrite(&widths_y[0], sizeof(double), num_y, out);
          fwrite(&widths_z[0], sizeof(double), num_z, out);
        }

        /* Get universes */
        for (const auto &lat_universe : *lattice)
        {
          int universe_id = lat_universe->getId();
          fwrite(&universe_id, sizeof(int), 1, out);
        }
      }
    }

    /* Close the output file */
    fclose(out);
  }

  /**
   * @brief Loads all Geometry and Material details from a Geometry restart file
   * @param filename The name of the file where the data is loaded
   * @param twiddle Whether the bytes are inverted (BGQ) or not
   */
  void Geometry::loadFromFile(std::string filename, bool twiddle)
  {

    _twiddle = twiddle;
    _loaded_from_file = true;

    FILE *in;
    in = fopen(filename.c_str(), "r");

    delete _root_universe;

    log::finfo("Reading Geometry from %s", filename.c_str());

    std::map<int, Surface *> all_surfaces;
    std::map<int, Cell *> all_cells;
    std::map<int, Universe *> all_universes;
    std::map<int, Material *> all_materials;

    std::map<int, int> fill_cell_universes;
    std::map<int, int> cell_parent;
    std::map<int, int *> lattice_universes;

    /* Read number of energy groups */
    int num_groups;
    int ret = twiddleRead(&num_groups, sizeof(int), 1, in);

    /* Read all material infromation */
    int num_materials;
    ret = twiddleRead(&num_materials, sizeof(int), 1, in);
    for (int i = 0; i < num_materials; i++)
    {

      /* Get key, value pair and cross section information */
      int key, id;
      int length;
      char *str;
      const char *name;
      ret = twiddleRead(&key, sizeof(int), 1, in);
      ret = twiddleRead(&id, sizeof(int), 1, in);
      ret = twiddleRead(&length, sizeof(int), 1, in);
      if (length > 0)
      {
        str = new char[length + 1];
        ret = twiddleRead(str, sizeof(char), length, in);
        str[length] = '\0';
        name = str;
      }
      else
      {
        name = "";
      }

      /* Create Material */
      all_materials[key] = new Material(id, name);
      if (strcmp(name, "") != 0)
        delete[] name;
      Material *mat = all_materials[key];
      mat->setNumEnergyGroups(num_groups);

      /* Set total cross-section */
      double value;
      for (int g = 0; g < num_groups; g++)
      {
        ret = twiddleRead(&value, sizeof(double), 1, in);
        mat->setSigmaTByGroup(value, g + 1);
      }

      /* Set fission cross-section */
      for (int g = 0; g < num_groups; g++)
      {
        ret = twiddleRead(&value, sizeof(double), 1, in);
        mat->setSigmaFByGroup(value, g + 1);
      }

      /* Set nu * fisison cross-section */
      for (int g = 0; g < num_groups; g++)
      {
        ret = twiddleRead(&value, sizeof(double), 1, in);
        mat->setNuSigmaFByGroup(value, g + 1);
      }

      /* Set neutron emission spectrum (chi) */
      for (int g = 0; g < num_groups; g++)
      {
        ret = twiddleRead(&value, sizeof(double), 1, in);
        mat->setChiByGroup(value, g + 1);
      }

      /* Set absorption cross section */
      for (int g = 0; g < num_groups; g++)
      {
        ret = twiddleRead(&value, sizeof(double), 1, in);
        mat->setSigmaAByGroup(value, g + 1);
      }

      /* Set scattering cross-section */
      for (int g = 0; g < num_groups; g++)
      {
        for (int gp = 0; gp < num_groups; gp++)
        {
          ret = twiddleRead(&value, sizeof(double), 1, in);
          mat->setSigmaSByGroup(value, g + 1, gp + 1);
        }
      }
    }

    /* Read root universe ID */
    int root_id;
    ret = twiddleRead(&root_id, sizeof(int), 1, in);

    /* Read all surface information */
    int num_surfaces;
    ret = twiddleRead(&num_surfaces, sizeof(int), 1, in);
    for (int i = 0; i < num_surfaces; i++)
    {

      /* Get key, value pair and general surface information */
      int key, id;
      int length;
      const char *name;
      ret = twiddleRead(&key, sizeof(int), 1, in);
      ret = twiddleRead(&id, sizeof(int), 1, in);
      ret = twiddleRead(&length, sizeof(int), 1, in);
      if (length > 0)
      {
        char *str = new char[length + 1];
        ret = twiddleRead(str, sizeof(char), length, in);
        str[length] = '\0';
        name = str;
      }
      else
      {
        name = "";
      }
      surfaceType st;
      boundaryType bt;
      ret = twiddleRead(&st, sizeof(surfaceType), 1, in);
      ret = twiddleRead(&bt, sizeof(boundaryType), 1, in);

      /* Treat specific surface types */
      if (st == PLANE)
      {
        double a, b, c, d;
        ret = twiddleRead(&a, sizeof(double), 1, in);
        ret = twiddleRead(&b, sizeof(double), 1, in);
        ret = twiddleRead(&c, sizeof(double), 1, in);
        ret = twiddleRead(&d, sizeof(double), 1, in);
        all_surfaces[key] = new Plane(a, b, c, d, id, name);
      }
      else if (st == ZCYLINDER)
      {
        double x, y, radius;
        ret = twiddleRead(&x, sizeof(double), 1, in);
        ret = twiddleRead(&y, sizeof(double), 1, in);
        ret = twiddleRead(&radius, sizeof(double), 1, in);
        all_surfaces[key] = new ZCylinder(x, y, radius, id, name);
      }
      else if (st == XPLANE)
      {
        double x;
        ret = twiddleRead(&x, sizeof(double), 1, in);
        all_surfaces[key] = new XPlane(x, id, name);
      }
      else if (st == YPLANE)
      {
        double y;
        ret = twiddleRead(&y, sizeof(double), 1, in);
        all_surfaces[key] = new YPlane(y, id, name);
      }
      else if (st == ZPLANE)
      {
        double z;
        ret = twiddleRead(&z, sizeof(double), 1, in);
        all_surfaces[key] = new ZPlane(z, id, name);
      }
      else
      {
        log::ferror("Unsupported surface type %s", name);
      }
      if (strcmp(name, "") != 0)
        delete[] name;

      /* Check that the key and ID match */
      if (key != id)
      {
        std::string str = all_surfaces[key]->toString();
        log::ferror("Surface key %d does not match it's corresponding ID "
                    "%d for surface:\n%s",
                    key, id, str.c_str());
      }

      /* Set boundary */
      all_surfaces[key]->setBoundaryType(bt);
    }

    /* Read all cell information */
    int num_cells;
    ret = twiddleRead(&num_cells, sizeof(int), 1, in);
    for (int i = 0; i < num_cells; i++)
    {

      /* Get key, value pair and general cell information */
      int key, id;
      int length;
      const char *name;
      ret = twiddleRead(&key, sizeof(int), 1, in);
      ret = twiddleRead(&id, sizeof(int), 1, in);
      ret = twiddleRead(&length, sizeof(int), 1, in);
      if (length > 0)
      {
        char *str = new char[length + 1];
        ret = twiddleRead(str, sizeof(char), length, in);
        str[length] = '\0';
        name = str;
      }
      else
      {
        name = "";
      }
      cellType ct;
      ret = twiddleRead(&ct, sizeof(cellType), 1, in);

      /* Create the cell */
      all_cells[key] = new Cell(id, name);
      if (strcmp(name, "") != 0)
        delete[] name;

      /* Fill the cell */
      if (ct == MATERIAL)
      {
        int mat_id;
        ret = twiddleRead(&mat_id, sizeof(int), 1, in);
        all_cells[key]->setFill(all_materials[mat_id]);
      }
      else if (ct == FILL)
      {
        int univ_id;
        ret = twiddleRead(&univ_id, sizeof(int), 1, in);
        fill_cell_universes[key] = univ_id;
      }

      /* Read cell rotations */
      bool rot;
      ret = twiddleRead(&rot, sizeof(bool), 1, in);
      if (rot)
      {
        double rotation[3];
        ret = twiddleRead(rotation, sizeof(double), 3, in);
        all_cells[key]->setRotation(rotation, 3, "radians");
      }

      /* Read cell translations */
      bool trans;
      ret = twiddleRead(&trans, sizeof(bool), 1, in);
      if (trans)
      {
        double translation[3];
        ret = twiddleRead(translation, sizeof(double), 3, in);
        all_cells[key]->setTranslation(translation, 3);
      }

      /* Read ring / sector information */
      int num_rings, num_sectors;
      ret = twiddleRead(&num_rings, sizeof(int), 1, in);
      ret = twiddleRead(&num_sectors, sizeof(int), 1, in);
      all_cells[key]->setNumRings(num_rings);
      all_cells[key]->setNumSectors(num_sectors);

      /* Read parent cell */
      bool has_parent;
      ret = twiddleRead(&has_parent, sizeof(bool), 1, in);
      if (has_parent)
      {
        int parent_id;
        ret = twiddleRead(&parent_id, sizeof(int), 1, in);
        cell_parent[key] = parent_id;
      }

      /* Read bounding surfaces */
      // FIXME, read regions
      int num_cell_surfaces;
      ret = twiddleRead(&num_cell_surfaces, sizeof(int), 1, in);
      for (int s = 0; s < num_cell_surfaces; s++)
      {
        int surface_id;
        int halfspace;
        ret = twiddleRead(&surface_id, sizeof(int), 1, in);
        ret = twiddleRead(&halfspace, sizeof(int), 1, in);
        all_cells[key]->addSurface(halfspace, all_surfaces[surface_id]);
      }

      /* Check that the key and ID match */
      if (key != id)
      {
        std::string str = all_cells[key]->toString();
        log::ferror("Cell key %d does not match its corresponding ID "
                    "%d for cell:\n%s",
                    key, id, str.c_str());
      }
    }

    /* Read all universe information */
    int num_universes;
    ret = twiddleRead(&num_universes, sizeof(int), 1, in);
    for (int i = 0; i < num_universes; i++)
    {

      /* Get key, value pair and general universe information */
      int key, id;
      int length;
      const char *name;
      ret = twiddleRead(&key, sizeof(int), 1, in);
      ret = twiddleRead(&id, sizeof(int), 1, in);
      ret = twiddleRead(&length, sizeof(int), 1, in);
      if (length > 0)
      {
        char *str = new char[length + 1];
        ret = twiddleRead(str, sizeof(char), length, in);
        str[length] = '\0';
        name = str;
      }
      else
      {
        name = "";
      }
      universeType ut;
      ret = twiddleRead(&ut, sizeof(universeType), 1, in);

      if (ut == SIMPLE)
      {

        /* Read all cells in the universe */
        all_universes[key] = new Universe(id, name);
        int num_universe_cells;
        ret = twiddleRead(&num_universe_cells, sizeof(int), 1, in);
        for (int c = 0; c < num_universe_cells; c++)
        {
          int cell_id;
          ret = twiddleRead(&cell_id, sizeof(int), 1, in);
          all_universes[key]->addCell(all_cells[cell_id]);
        }
      }
      else if (ut == LATTICE)
      {

        /* Read lattice information */
        int num_x, num_y, num_z;
        double width_x, width_y, width_z;
        double offset[3];
        ret = twiddleRead(&num_x, sizeof(int), 1, in);
        ret = twiddleRead(&num_y, sizeof(int), 1, in);
        ret = twiddleRead(&num_z, sizeof(int), 1, in);
        ret = twiddleRead(&width_x, sizeof(double), 1, in);
        ret = twiddleRead(&width_y, sizeof(double), 1, in);
        ret = twiddleRead(&width_z, sizeof(double), 1, in);
        ret = twiddleRead(offset, sizeof(double), 3, in);

        std::vector<double> widths_x(num_x), widths_y(num_y), widths_z(num_z);
        bool non_uniform = false;
        ret = twiddleRead(&non_uniform, sizeof(bool), 1, in);

        /* Read widths vectors if the lattice is non-uniform */
        if (non_uniform)
        {
          ret = twiddleRead(&widths_x[0], sizeof(double), num_x, in);
          ret = twiddleRead(&widths_y[0], sizeof(double), num_y, in);
          ret = twiddleRead(&widths_z[0], sizeof(double), num_z, in);
        }

        /* Create lattice */
        RecLattice *new_lattice = new RecLattice(id, name);
        all_universes[key] = new_lattice;
        new_lattice->setNumX(num_x);
        new_lattice->setNumY(num_y);
        new_lattice->setNumZ(num_z);
        if (non_uniform)
        {
          new_lattice->setWidths(widths_x, widths_y, widths_z);
          new_lattice->setWidth(1, 1, 1);
        }
        else
          new_lattice->setWidth(width_x, width_y, width_z);

        new_lattice->setNonUniform(non_uniform);
        new_lattice->setOffset(offset[0], offset[1], offset[2]);

        /* Get universes */
        lattice_universes[key] = new int[num_x * num_y * num_z];
        for (int j = 0; j < num_x * num_y * num_z; j++)
        {
          int universe_id;
          ret = twiddleRead(&universe_id, sizeof(int), 1, in);
          lattice_universes[key][j] = universe_id;
        }
      }
      if (strcmp(name, "") != 0)
        delete[] name;

      /* Check that the key and ID match */
      if (key != id)
      {
        std::string str = all_universes[key]->toString();
        log::ferror("Universe key %d does not match it's corresponding ID "
                    "%d for surface:\n%s",
                    key, id, str.c_str());
      }
    }

    /* Set universe fills in cells */
    std::map<int, int>::iterator id_iter;
    for (id_iter = fill_cell_universes.begin();
         id_iter != fill_cell_universes.end(); ++id_iter)
      all_cells[id_iter->first]->setFill(all_universes[id_iter->second]);

    /* Set parent cells */
    for (id_iter = cell_parent.begin(); id_iter != cell_parent.end(); ++id_iter)
      all_cells[id_iter->first]->setParent(all_cells[id_iter->second]);

    /* Set lattice universes */
    for (auto &l : lattice_universes)
    {
      int id = l.first;
      int *array = l.second;
      RecLattice *lattice = static_cast<RecLattice *>(all_universes[id]);
      int size = lattice->getNumLatticeCells();
      auto universes = new Universe *[size];
      for (int i = 0; i < size; i++)
      {
        universes[i] = all_universes[array[i]];
      }
      int nz = lattice->getNumZ();
      if (lattice->getLatticeType() == latticeType::Rectangle)
      {
        int nx = lattice->getNumX();
        int ny = lattice->getNumY();
        lattice->setUniverses(nz, ny, nx, universes);
        delete[] array;
      }
      else
        log::ferror("Fix me: read Hexagon lattices");
    }

    /* Set root universe */
    _root_universe = all_universes[root_id];

    /* Close the input file */
    fclose(in);

    log::fverbose("Status of reading Geometry and Materials from file: %d", ret);
    log::finfo("Read complete");
  }

  /**
   * @brief Read an integer array from file.
   * @param ptr the integer array to fill with the data read
   * @param size the size of each element to read (here size(int))
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(int *ptr, size_t size, size_t nmemb,
                               FILE *stream)
  {
    size_t ret = fread(ptr, size, nmemb, stream);
    if (_twiddle)
      for (size_t i = 0; i < nmemb; i++)
        ptr[i] = __builtin_bswap32(ptr[i]);
    return ret;
  }

  /**
   * @brief Read a boolean array from file.
   * @param ptr the boolean array to fill with the data read
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(bool *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    size_t ret = fread(ptr, size, nmemb, stream);
    return ret;
  }

  /**
   * @brief Read an array of universeType from file.
   * @param ptr the array to fill with the data read
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(universeType *ptr, size_t size, size_t nmemb,
                               FILE *stream)
  {
    size_t ret = fread(ptr, size, nmemb, stream);
    int *arr = reinterpret_cast<int *>(ptr);
    if (_twiddle)
      for (size_t i = 0; i < nmemb; i++)
        arr[i] = __builtin_bswap32(arr[i]);
    return ret;
  }

  /**
   * @brief Read an array of cellType from file.
   * @param ptr the array to fill with the read data
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(cellType *ptr, size_t size, size_t nmemb,
                               FILE *stream)
  {
    size_t ret = fread(ptr, size, nmemb, stream);
    int *arr = reinterpret_cast<int *>(ptr);
    if (_twiddle)
      for (size_t i = 0; i < nmemb; i++)
        arr[i] = __builtin_bswap32(arr[i]);
    return ret;
  }

  /**
   * @brief Read an array of surfaceType from file.
   * @param ptr the array to fill with the data read
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(surfaceType *ptr, size_t size, size_t nmemb,
                               FILE *stream)
  {
    size_t ret = fread(ptr, size, nmemb, stream);
    int *arr = reinterpret_cast<int *>(ptr);
    if (_twiddle)
      for (size_t i = 0; i < nmemb; i++)
        arr[i] = __builtin_bswap32(arr[i]);
    return ret;
  }

  /**
   * @brief Read an array of boundaryType from file.
   * @param ptr the array to fill with the data read
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(boundaryType *ptr, size_t size, size_t nmemb,
                               FILE *stream)
  {
    size_t ret = fread(ptr, size, nmemb, stream);
    int *arr = reinterpret_cast<int *>(ptr);
    if (_twiddle)
      for (size_t i = 0; i < nmemb; i++)
        arr[i] = __builtin_bswap32(arr[i]);
    return ret;
  }

  /**
   * @brief Read an array of char from file.
   * @param ptr the array to fill with the data read
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(char *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    size_t ret = fread(ptr, size, nmemb, stream);
    return ret;
  }

  /**
   * @brief Read an array of double from file.
   * @param ptr the array to fill with the data read
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(double *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    long *arr = reinterpret_cast<long *>(ptr);
    size_t ret = fread(arr, size, nmemb, stream);
    if (_twiddle)
      for (size_t i = 0; i < nmemb; i++)
        arr[i] = __builtin_bswap64(arr[i]);
    return ret;
  }

  /**
   * @brief Read an array of long int from file.
   * @param ptr the array to fill with the data read
   * @param size the size of each element to read
   * @param nmemb the number of elements to read
   * @param stream the file to read from
   * @return the return status of the read operation
   */
  size_t Geometry::twiddleRead(long *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    long *arr = ptr;
    size_t ret = fread(arr, size, nmemb, stream);
    if (_twiddle)
      for (size_t i = 0; i < nmemb; i++)
        arr[i] = __builtin_bswap64(arr[i]);
    return ret;
  }

  //----------------------------------------------------------------------
  // Spatial domain decomposition
  //----------------------------------------------------------------------

#ifdef ENABLE_MPI_
  /// \brief Domain decomposes the Geometry with MPI and modular ray tracing
  void Geometry::setDomainDecomposition()
  {

    /* Check that the root universe has been set */
    if (_root_universe == nullptr)
      log::ferror("The root universe must be set before domain "
                  "decomposition.");

    if (!mpi::isMPIInitialized())
    {
      log::ferror("MPI communicators must be initialized before setting "
                  "spatial domains for geometry");
    }

    const auto &nx = mpi::getNumDomainsX();
    const auto &ny = mpi::getNumDomainsY();
    const auto &nz = mpi::getNumDomainsZ();

    /* Check that the Geometry needs to be decomposed */
    if (mpi::getNumUniqueDomains() > 1)
    {

      /* Set the bounds in a 1x1x1 Lattice object */
      _domain_bounds = new RecLattice();
      _domain_bounds->setNumX(1);
      _domain_bounds->setNumY(1);
      _domain_bounds->setNumZ(1);
      double width_x = getWidthX();
      double width_y = getWidthY();
      double width_z = getWidthZ();
      _domain_bounds->setWidth(width_x, width_y, width_z);
      double offset_x = width_x / 2.0 + getMinX();
      double offset_y = width_y / 2.0 + getMinY();
      double offset_z = width_z / 2.0 + getMinZ();
      _domain_bounds->setOffset(offset_x, offset_y, offset_z);
      _domain_bounds->computeSizes();
      log::finfo("Successfully set %d x %d x %d spatial domain decomposition",
                 nx, ny, nz);
    }
  }

  /// \brief Returns the MPI communicator to communicate between MPI domains
  /// \return The MPI communicator
  MPI_Comm Geometry::getMPICart()
  {
    if (!isDomainDecomposed())
      log::ferror("Tried to get MPI Cart but domain decomposition has not "
                  "yet been set");

    // This communicator contains processes which have unique spatial domains
    return mpi::getCommUniqueDomains();
  }

  /// \brief See mpi::getNeighborDomain
  int Geometry::getNeighborDomain(int offset_x, int offset_y, int offset_z)
  {
    return mpi::getNeighborDomain(offset_x, offset_y, offset_z);
  }

  /// \brief Counts the number of FSRs in each MPI domain
  void Geometry::countDomainFSRs()
  {

    /* Gather the number of FSRs into an array */
    int num_domains = mpi::getNumUniqueDomains();
    long num_domains_array[num_domains];
    long my_fsrs = getNumFSRs();
    MPI_Allgather(&my_fsrs, 1, mpi::getDatatype<long>(),
                  num_domains_array, 1, mpi::getDatatype<long>(),
                  getMPICart());

    /* Convert to a vector */
    _num_domain_FSRs.resize(num_domains);
    for (int i = 0; i < num_domains; i++)
      _num_domain_FSRs.at(i) = num_domains_array[i];
    _domain_FSRs_counted = true;
  }

  /**
   * @brief Finds the local FSR ID and domain rank from a global FSR ID
   * @param global_fsr_id The global unique identifier of an FSR
   * @param local_fsr_id The unique identifier of an FSR on its current domain
   * @param domain The rank of the domain containing the FSR
   */
  void Geometry::getLocalFSRId(long global_fsr_id, long &local_fsr_id,
                               int &domain)
  {

    /* Count FSRs on each domain if not already counted */
    if (!_domain_FSRs_counted)
      countDomainFSRs();

    /* Determine the local domain where the global FSR resides */
    long cum_fsrs = 0;
    domain = -1;
    for (int i = 0; i < mpi::getNumUniqueDomains(); i++)
    {
      if (cum_fsrs + _num_domain_FSRs.at(i) > global_fsr_id)
      {
        domain = i;
        break;
      }
      else
      {
        cum_fsrs += _num_domain_FSRs.at(i);
      }
    }

    /* Ensure a domain was found with the FSR ID */
    if (domain == -1)
      log::ferror("No domain was found with the global FSR ID %ld. The "
                  "total number of FSRs in the Geometry is %ld.",
                  global_fsr_id,
                  getNumTotalFSRs());

    local_fsr_id = global_fsr_id - cum_fsrs;
  }

#endif

  /**
   * @brief Returns the FSR centroid of a global FSR
   * @param global_fsr_id The global unique identifier of the FSR
   * @return A vector contianing the coordinates of the FSR centroid
   */
  std::vector<double> Geometry::getGlobalFSRCentroidData(long global_fsr_id)
  {
    double xyz[3];
    Point *centroid = getFSRCentroid(global_fsr_id);
    xyz[0] = centroid->getX();
    xyz[1] = centroid->getY();
    xyz[2] = centroid->getZ();

#ifdef ENABLE_MPI_
    if (isDomainDecomposed())
    {

      /* Determine the domain and local FSR ID */
      long fsr_id;
      int domain;
      getLocalFSRId(global_fsr_id, fsr_id, domain);

      /* Get the FSR centroid in the correct domain */
      int rank = mpi::getRankUniqueDomains();
      MPI_Comm_rank(getMPICart(), &rank);
      double temp_xyz[3];
      if (rank == domain)
      {
        Point *centroid = getFSRCentroid(fsr_id);
        temp_xyz[0] = centroid->getX();
        temp_xyz[1] = centroid->getY();
        temp_xyz[2] = centroid->getZ();
      }
      else
      {
        temp_xyz[0] = 0;
        temp_xyz[1] = 0;
        temp_xyz[2] = 0;
      }

      /* Broadcast the centroid */
      MPI_Allreduce(temp_xyz, xyz, 3, MPI_DOUBLE, MPI_SUM, getMPICart());
    }
#endif

    /* Convert centroid data into a vector */
    std::vector<double> data(3);
    for (int i = 0; i < 3; i++)
      data.at(i) = xyz[i];
    return data;
  }

  //----------------------------------------------------------------------
  // Periodic track decomposition
  //----------------------------------------------------------------------

#ifdef ENABLE_MPI_
  /// \brief Reduce all ExtrudedFSRs/FSRs by merging all of the LocalCoords.
  /// \details This method defines a custom AllReduce procedure for FSRs.
  ///          In fact, FSRs are too heavy to synchronize between processes
  ///          so that we have to reduce all of the effective coordinates of
  ///          LocalCoords and then re-construct each of the FSRs.
  /// \param extruded Reduce extruded FSRs or 3D FSRs.
  void Geometry::allReduceFSRs(const std::vector<double> &z_coords, bool extruded)
  {

    if (!mpi::isPrdTrackDecomposed())
      return;

    std::string fsr_name = extruded ? "ExtrudedFSR" : "FSR";

    const auto comm = mpi::getCommSharedDomain();
    const int rank = mpi::getRankSharedDomain();
    const int np = mpi::getNumProcsSharedDomain();

    // Barrier processes
    mpi::mpiBarrierSharedDomain();

    log::info("Reducing {}s over {} processes...", fsr_name, np);

    // Create custom types for communication
    struct SimplePoint
    {
      double x, y, z;
    };

    MPI_Datatype PointMPI;
    MPI_Type_contiguous(3, MPI_DOUBLE, &PointMPI);
    MPI_Type_commit(&PointMPI);

    MPI_Aint mpi_point_lb;
    MPI_Aint mpi_point_size;
    MPI_Type_get_extent(PointMPI, &mpi_point_lb, &mpi_point_size);

    if (sizeof(SimplePoint) != mpi_point_size)
      log::error("Size of SimplePoint must equal to size of mpi datatype");

    // Max depth starting from 0. This level is computed such that
    // 2^(L-1) - 1 < np - 1 <= 2^L - 1
    // That is, the minimum required level which contains at least np nodes.
    int max_tree_level = (np == 1 ? 0 : std::ceil(std::log2(np)));
    log::verbose_once("Reducing {}s, max tree level = {}", fsr_name, max_tree_level);

    // Compute the merge level for myself. My FSRs will be merged into
    // another process at this level. If I'm root, the merge level will be
    // max_tree_level, which means there is nothing to merge all the time.
    // Otherwise, the merging happens at a level less than max_tree_level.
    int merge_level = 0;
    int leaf = rank;
    while (merge_level < max_tree_level && leaf % 2 == 0)
    {
      leaf /= 2;
      merge_level++;
    }

    log::debug("Reducing {}s, merge level = {}", fsr_name, merge_level);

    // A function to get FSR map size
    auto fsr_map_size = [&]()
    {
      if (extruded)
        return _extruded_FSR_keys_map.size();
      else
        return _FSR_keys_map.size();
    };

    // A function to fill a buffer with points
    auto fill_points = [&](SimplePoint *buffer)
    {
      ExtrudedFSR **extruded_fsrs = nullptr;
      FSRData **fsrs = nullptr;
      long num_points;

      if (extruded)
      {
        extruded_fsrs = _extruded_FSR_keys_map.values();
        num_points = _extruded_FSR_keys_map.size();
      }
      else
      {
        fsrs = _FSR_keys_map.values();
        num_points = _FSR_keys_map.size();
      }

      for (long i = 0; i < num_points; ++i)
      {
        Point *p;
        if (extruded)
          p = extruded_fsrs[i]->_coords->getPoint();
        else
          p = fsrs[i]->_point;

        buffer[i].x = p->getX();
        buffer[i].y = p->getY();
        buffer[i].z = p->getZ();
      }

      delete[] extruded_fsrs;
      delete[] fsrs;
    };

    // A function for merging fsrs
    auto merge_extruded_fsr = [&](const SimplePoint *buffer, const long i)
    {
      const auto &p = buffer[i];
      LocalCoords point(p.x, p.y, p.z, true);
      point.setUniverse(_root_universe);

      // Make a linked list
      findCellContainingCoords(&point);

      long region_id;
      if (extruded)
        region_id = createUniqueExtrudedFSR(point, z_coords);
      else
        region_id = findFSRId(&point);

      // Check that the region ID is valid
      if (region_id == -1)
        log::error("Failed to find a valid FSR when re-constructing {} "
                   "from point {}",
                   fsr_name, point.toString());
    };

    // Receive FSRs from adjacent processes.
    // The size of ExtrudedFSR map remains unknown until a process merges received
    // points into its ExtrudedFSR map.
    // This relies on the fact that the difference between the numbers of
    // two adjacent sub-trees is 2^L, where L is the number of nodes in a
    // sub-tree. That is,
    //     The ordinal n of process p1 at level L is
    //     p1 = 2^L * (2*n)
    //     The ordinal n of process p2 adjacent to p1 at level L is
    //     p2 = 2^L * (2*n + 1)
    for (int level = 0; level < merge_level; ++level)
    {
      int sender = rank + (1 << level);
      if (sender < np)
      {

        // Receive the data size
        int tag = rank;
        int recv_size;
        MPI_Status status;
        MPI_Probe(sender, tag, comm, &status);
        MPI_Get_count(&status, PointMPI, &recv_size);

        log::debug("Reducing {}s, receiving {} points from {} at level {}",
                   fsr_name, recv_size, sender, level);

        // Receive data
        auto recv_buf = new SimplePoint[recv_size];
        MPI_Recv(recv_buf, recv_size, PointMPI, sender, tag, comm, &status);

        // Merge data. An ExtrudedFSR is reconstructed everytime the associated
        // LocalCoords is merged.
        size_t num_before_merge = fsr_map_size();
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < recv_size; ++i)
          merge_extruded_fsr(recv_buf, i);

        delete[] recv_buf;

        log::debug("Reducing {}s, {} points received, {} merged at level {}",
                   fsr_name, recv_size, fsr_map_size() - num_before_merge, level);
      }
    }

    // Send data to adjacent processes. Process 0 will not participate.
    if (rank)
    {
      int receiver = rank - (1 << merge_level);
      int tag = receiver;
      log::debug("Reducing {}s, send {} points to {} at level {}",
                 fsr_name, fsr_map_size(), receiver, merge_level);

      // Get all of the FSRs
      size_t num_points = fsr_map_size();

      // Fill buffer
      // TODO: eliminate the send buffer by indexed type
      auto send_buf = new SimplePoint[num_points];
      fill_points(send_buf);

      // Send the data size
      MPI_Send(send_buf, num_points, PointMPI, receiver, tag, comm);
      delete[] send_buf;
    }

    // Synchronize ExtrudedFSR map size
    long num_points = fsr_map_size();
    MPI_Bcast(&num_points, 1, mpi::getDatatype<long>(), 0, comm);

    // Broadcast points
    auto sendrecv_buf = new SimplePoint[num_points];
    if (rank == 0)
    {
      fill_points(sendrecv_buf);
    }
    MPI_Status status;
    MPI_Request handle;
    MPI_Ibcast(sendrecv_buf, num_points, PointMPI, 0, comm, &handle);

    // Clear FSRs and re-construct them to make them consistent
    // among all of the processes
    if (extruded)
      _extruded_FSR_keys_map.clear();
    else
      _FSR_keys_map.clear();

    MPI_Wait(&handle, &status);
    for (long i = 0; i < num_points; ++i)
      merge_extruded_fsr(sendrecv_buf, i);

    delete[] sendrecv_buf;

    log::fverbose_once("Reducing ExtrudedFSRs, %d points have been broadcasted "
                       "and reconstructed",
                       num_points, np - 1);

    // Cleanup
    MPI_Type_free(&PointMPI);

    _extruded_FSRs_reduced = true;
  }

#endif

  //----------------------------------------------------------------------
  // Profiling helpers
  //----------------------------------------------------------------------

  /// \brief Print the number of extruded FSR
  /// \details If spatial domain decomposition takes effect, processes will
  ///          be barriered. This method should be called after ExtrudedFSRs
  ///          are created.
  void Geometry::printNumExtrudedFSRs()
  {
    long num_ext_fsrs = _extruded_FSR_keys_map.size();
    long total_num_ext_fsrs = num_ext_fsrs;
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
    {
      MPI_Allreduce(&num_ext_fsrs, &total_num_ext_fsrs, 1, mpi::getDatatype<long>(),
                    MPI_SUM, getMPICart());
    }
    log::fverbose("Number of extruded FSRs in domain = %ld", num_ext_fsrs);
#endif
    log::finfo("Total number of extruded FSRs = %ld", total_num_ext_fsrs);
  }

  /// \brief Print the number of FSRs
  /// \details If spatial domain decomposition takes effect, processes will
  ///          be barriered.
  void Geometry::printNumFSRs()
  {
    auto num_fsrs = _FSR_keys_map.size();
    auto total = getNumTotalFSRs();

    if (total <= 0)
      log::error("Invalid total number of FSRs: {}", total);

    log::info("Total number of FSRs = {}", total);
    if (mpi::isSpatialDecomposed())
      log::verbose("Number of FSRs in domain = {}", num_fsrs);
  }

  /// \brief Output the extruded FSR storage requirement
  void Geometry::printMemUsageExtrudedFSRs(long num_fsrs_in_stack)
  {
    float size = num_fsrs_in_stack * (sizeof(double) + sizeof(long) + sizeof(Material *)) +
                 _extruded_FSR_keys_map.size() * (sizeof(_extruded_FSR_keys_map.keys()[0]) + sizeof(ExtrudedFSR) + (LOCAL_COORDS_LEN + 1) * sizeof(LocalCoords));

    float max_size = size;
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
      MPI_Allreduce(&size, &max_size, 1, mpi::getDatatype<float>(), MPI_MAX, getMPICart());
#endif

    log::finfo("Max storage for extruded FSRs per domain = %.2f MiB",
               max_size / (1 << 20));
  }

  /// \brief Output the FSR storage requirement (no ExtrudedFSR)
  void Geometry::printMemUsageFSRs()
  {
    float size = _FSR_keys_map.size() // ParallelHashMap
                     * (sizeof(FSRData) + sizeof(FSRData *) + sizeof(std::string) + sizeof(omp_lock_t)) +
                 _FSRs_to_keys.capacity() * sizeof(decltype(_FSRs_to_keys)::value_type) + sizeof(_FSRs_to_keys) + _FSRs_to_centroids.capacity() * (sizeof(Point *) + sizeof(Point)) + sizeof(_FSRs_to_centroids) + _FSRs_to_material_IDs.capacity() * sizeof(decltype(_FSRs_to_material_IDs)::value_type) + sizeof(_FSRs_to_material_IDs) + _FSRs_to_CMFD_cells.capacity() * sizeof(decltype(_FSRs_to_CMFD_cells)::value_type) + sizeof(_FSRs_to_CMFD_cells);

    float max_size = size;
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
      MPI_Allreduce(&size, &max_size, 1, mpi::getDatatype<float>(), MPI_MAX, getMPICart());
#endif

    log::finfo("Max storage for FSRs per domain = %.2f MiB",
               max_size / (1 << 20));
  }

  /// \brief Output the Material storage requirement (no material name)
  void Geometry::printMemUsageMaterials()
  {
    const auto &materials = getAllMaterials();
    const auto &num_groups = getNumEnergyGroups();
    const auto &m = materials.begin();

    unsigned count_sigma = static_cast<unsigned>(m->second->getSigmaT() != nullptr) + static_cast<unsigned>(m->second->getSigmaA() != nullptr) + static_cast<unsigned>(m->second->getSigmaF() != nullptr) + static_cast<unsigned>(m->second->getNuSigmaF() != nullptr) + static_cast<unsigned>(m->second->getChi() != nullptr);

    // Only 0-order scattering
    unsigned count_sigma_s = static_cast<unsigned>(m->second->getSigmaS() != nullptr);

    float size = sizeof(materials) + materials.size() // Material caching map
                                         * (sizeof(*m) + count_sigma * num_groups * sizeof(FP_PRECISION) + count_sigma_s * num_groups * num_groups * sizeof(FP_PRECISION));

    float max_size = size;
#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
      MPI_Allreduce(&size, &max_size, 1, mpi::getDatatype<float>(), MPI_MAX, getMPICart());
#endif

    log::finfo("Max storage for Materials per domain = %.2f MiB",
               max_size / (1 << 20));
  }

} /* namespace antmoc */
