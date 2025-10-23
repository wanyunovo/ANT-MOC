/**
 * @file Geometry.h
 * @brief The Geometry class.
 * @date January 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <map>
#include <string>
#include <vector>

#include "antmoc/enum_types.h"
#include "antmoc/FSRDataTypes.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/ParallelHashMap.h"

namespace antmoc
{

/** Forward declarationss */
class Cell;
class Cmfd;
class RecLattice;
class LocalCoords;
class Material;
class Point;
class Surface;
class Track;
class Track3D;
class Universe;



void reset_auto_ids();


/**
 * @class Geometry Geometry.h "include/Geometry.h"
 * @brief The master class containing references to all geometry-related
 *        objects - Surfaces, Cells, Universes and Lattices - and Materials.
 * @details The primary purpose for the geometry is to serve as a collection
 *          of all geometry-related objects, as well as for ray tracing
 *          of characteristic tracks across the Geometry and computing
 *          FSR-to-cell offset maps.
 */
class Geometry {

private:

  /** An map of FSR key hashes to unique FSRData structs */
  ParallelHashMap<std::string, FSRData*> _FSR_keys_map;
  ParallelHashMap<std::string, ExtrudedFSR*> _extruded_FSR_keys_map;

  /** An vector of FSR key hashes indexed by FSR ID */
  std::vector<std::string> _FSRs_to_keys;

  /** An vector of FSR centroids indexed by FSR ID */
  std::vector<Point*> _FSRs_to_centroids;

  /** A boolean indicating whether any centroids have been set */
  bool _contains_FSR_centroids;

  /** A vector of Material IDs indexed by FSR IDs */
  std::vector<int> _FSRs_to_material_IDs;

  /** A vector of ExtrudedFSR pointers indexed by extruded FSR ID 由轴向挤出的FSR ID索引的轴向挤出FSR指针的Vector数组*/
  std::vector<ExtrudedFSR*> _extruded_FSR_lookup;

  /** An vector of CMFD cell IDs indexed by FSR ID */
  std::vector<int> _FSRs_to_CMFD_cells;

  /* The Universe at the root node in the CSG tree */
  Universe* _root_universe;

  /** A CMFD object pointer */
  Cmfd* _cmfd;

  /** An optional axial mesh overlaid on the Geometry 这个变量不是很理解 文档的解释是“额外定义的一层网，如果这个网轴向有 nz 层，每层都算上”是轴向网吗？*/
  RecLattice* _overlaid_mesh;

  /* A map of all Material in the Geometry for optimization purposes */
  std::map<int, Material*> _all_materials;

  /* A vector containing allocated strings for key generation */
  std::vector<std::string> _fsr_keys;

  /* A boolean to know whether FSRs were counted */
  bool _domain_FSRs_counted;
  
  /* Lattice object of size 1 that contains the local domain */
  RecLattice* _domain_bounds;
  
  /* Number of FSRs in each domain */
  std::vector<long> _num_domain_FSRs;

  /* Number of modular track laydowns within a domain, in the X, Y and 
     Z directions */
  int _num_modules_x;
  int _num_modules_y;
  int _num_modules_z;

  /* Symmetries in X,Y and Z axis used to restrict the computation domain */
  std::vector<bool> _symmetries;

  /* Wether to read bites backwards or forward (for BGQ) */
  bool _twiddle;
  
  /* Whether the geometry was loaded from a .geo file or created in the input 
   * file. This matters for memory de-allocation purposes. */
  bool _loaded_from_file;

  ///< A boolean to know whether energy groups of all materials are consistent
  bool _energy_group_checked;

  ///< A boolean to know whether ExtrudedFSRs were reduced
  bool _extruded_FSRs_reduced;

  /* Function to find the cell containing the coordinates */
  Cell* findFirstCell(LocalCoords* coords, double azim, double polar=M_PI_2);
  
public:

  Geometry();
  virtual ~Geometry();

  /* Handle number of ray-tracing modules in a domain */
  void setNumDomainModules(int num_x, int num_y, int num_z);
  virtual int getNumXModules();
  virtual int getNumYModules();
  virtual int getNumZModules();

  /* Get parameters */
  virtual double getWidthX();
  virtual double getWidthY();
  virtual double getWidthZ();
  virtual double getMinX();
  virtual double getMaxX();
  virtual double getMinY();
  virtual double getMaxY();
  virtual double getMinZ();
  virtual double getMaxZ();
  virtual double getGlobalMinX();
  virtual double getGlobalMaxX();
  virtual double getGlobalMinY();
  virtual double getGlobalMaxY();
  virtual double getGlobalMinZ();
  virtual double getGlobalMaxZ();
  virtual boundaryType getMinXBoundaryType();
  virtual boundaryType getMaxXBoundaryType();
  virtual boundaryType getMinYBoundaryType();
  virtual boundaryType getMaxYBoundaryType();
  virtual boundaryType getMinZBoundaryType();
  virtual boundaryType getMaxZBoundaryType();
  Universe* getRootUniverse();
  long getNumFSRs();
  long getNumExtrudedFSRs();
  long getNumTotalFSRs();
  int getNumEnergyGroups();
  int getNumMaterials();
  int getNumCells();
  std::map<int, Surface*> getAllSurfaces();
  std::map<int, Material*>& getAllMaterials();
  std::map<int, Cell*> getAllCells();
  std::map<int, Cell*> getAllMaterialCells();
  std::map<int, Universe*> getAllUniverses();
  std::vector<double> getUniqueZHeights(bool include_overlaid_mesh = false);
  std::vector<double> getUniqueZPlanes();
  bool isDomainDecomposed();
  bool isRootDomain();
  void getDomainIndexes(int* indexes);
  void getDomainStructure(int* structure);

  /* Assign root universe to geometry */
  void setRootUniverse(Universe* root_universe);

  /* More complex getter methods */
  Cmfd* getCmfd();
  std::vector<std::string>& getFSRsToKeys();
  std::vector<int>& getFSRsToMaterialIDs();
  std::vector<Point*>& getFSRsToCentroids();
  std::vector<int>& getFSRsToCMFDCells();
  std::vector<ExtrudedFSR*>& getExtrudedFSRLookup();
  std::map<int, int> getMaterialIDsToIndices();
  long getFSRId(LocalCoords* coords, bool err_check=true);
  long getGlobalFSRId(LocalCoords* coords, bool err_check=true);
  Point* getFSRPoint(long fsr_id);
  Point* getFSRCentroid(long fsr_id);
  bool containsFSRCentroids();
  int getCmfdCell(long fsr_id);
  ExtrudedFSR* getExtrudedFSR(int extruded_fsr_id);
  std::string getFSRKey(LocalCoords* coords);
  void getFSRKeyFast(LocalCoords* coords, std::string& key);
  void printToString(std::string& str, int& index, int value);
  int getNumDigits(int number);
  ParallelHashMap<std::string, FSRData*>& getFSRKeysMap();
  ParallelHashMap<std::string, ExtrudedFSR*>& getExtrudedFSRKeysMap();

  /* Setter methods */
  void setCmfd(Cmfd* cmfd);
  void setFSRCentroid(long fsr, Point* centroid);
  void setOverlaidMesh(double axial_mesh_height, int num_x=0,
                       int num_y=0, int num_radial_domains=0,
                       int* radial_domains=NULL);
  void clearBoundaries();

  /* Find methods */
  Cell* findCellContainingCoords(LocalCoords* coords);
  Material* findFSRMaterial(long fsr_id);
  long findFSRId(LocalCoords* coords);
  int findExtrudedFSR(LocalCoords* coords);
  Cell* findCellContainingFSR(long fsr_id);
  Cell* findNextCell(LocalCoords* coords, double azim, double polar=M_PI_2);


  /* Other worker methods */
  void retrieveAllMaterials();
  void reserveKeyStrings(int num_threads);
  void subdivideCells();
  void initializeAxialFSRs(std::vector<double> global_z_mesh);
  void reorderFSRIDs();
  void initializeFlatSourceRegions();
  bool getSymmetry(int axis);
  void segmentize2D(Track* track, double z_coord,
                          bool FSRonly = false);
  void segmentize3D(Track3D* track, bool setup=false);
  void segmentizeExtruded(Track* flattened_track,
                          std::vector<double> z_coords,
                          bool FSRonly = false);
  int createUniqueExtrudedFSR(LocalCoords &start,
                              const std::vector<double> &z_coords);
  void fixFSRMaps();
  void initializeFSRVectors();
  void computeFissionability(Universe* univ=NULL);
  void manipulateXS();

  /* Obtain or print information about the geometry */
  std::vector<long> getSpatialDataOnGrid(std::vector<double> dim1,
                                         std::vector<double> dim2,
                                         double offset,
                                         const char* plane,
                                         const char* domain_type);
  std::string toString();
  void printString();
  void printFSRsToFile(const char* plane="xy", int gridsize=1000, 
                       double offset=0.0, double* bounds_x = NULL, 
                       double* bounds_y = NULL, double* bounds_z = NULL);

  void initializeCmfd();
  void initializeSpectrumCalculator(Cmfd* spectrum_calculator);
  bool withinBounds(LocalCoords* coords);
  bool withinGlobalBounds(LocalCoords* coords);

  std::vector<double> getGlobalFSRCentroidData(long global_fsr_id);
  int getDomainByCoords(LocalCoords* coords);
  std::map<Cell*, std::vector<long> > getCellsToFSRs();

  /* Input/output of geometries from/to .geo files */
  void dumpToFile(std::string filename);
  void loadFromFile(std::string filename, bool twiddle=false);
  size_t twiddleRead(int* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(bool* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(char* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(universeType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(cellType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(surfaceType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(boundaryType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(double* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(long* ptr, size_t size, size_t nmemb, FILE* stream);


#ifdef ENABLE_MPI_
  /* Set up domain decomposition */
  void setDomainDecomposition();
  MPI_Comm getMPICart();
  int getNeighborDomain(int offset_x, int offset_y, int offset_z);
  void countDomainFSRs();
  void getLocalFSRId(long global_fsr_id, long &local_fsr_id, int &domain);

  void allReduceFSRs(const std::vector<double> &z_coords, bool extruded = true);
#endif

  /* Profiling helpers */
  void printNumExtrudedFSRs();
  void printNumFSRs();
  void printMemUsageExtrudedFSRs(long num_fsrs_in_stack);
  void printMemUsageFSRs();
  void printMemUsageMaterials();
};


} /* namespace antmoc */

#endif /* GEOMETRY_H_ */
