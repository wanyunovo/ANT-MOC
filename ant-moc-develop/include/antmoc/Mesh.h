/// \file Mesh.h
/// \brief The Mesh class for the alternative C++ build.
/// \date November, 2016
/// \author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
///         Wang An, USTB (wangan@xs.ustb.edu.cn)

#ifndef MESH_H
#define MESH_H

#include <memory>
#include <set>
#include <unordered_map>

#include "antmoc/enum_types.h"
#include "antmoc/FSRDataHandler.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/Point.h"
#include "antmoc/PyVector.h"
#include "antmoc/tally_utils.h"
#include "antmoc/xml_utils.h"

namespace antmoc
{

// Forward declarationss
class Solver;
class Lattice;
class TrackGenerator;

using SolverPtr = std::shared_ptr<Solver>;
using TrackGeneratorPtr = std::shared_ptr<TrackGenerator>;

// A Vector3D is simply a 3-dimensional std::vector of floats
typedef std::vector<std::vector<std::vector<FP_PRECISION> > > Vector3D;


/// \class Mesh
/// \brief A Mesh is a lattice overlaid on the Geometry across which reaction
///        rates can be tallied from converged scalar fluxes in a Solver
class Mesh : public FSRDataHandler {

  /* The lattice defining the zones across which reaction rates are tallied */
  Lattice* _lattice;

  /* A flag indicating whether a lattice has been allocated internally */
  bool _lattice_allocated;

  ///< The number of points outside the mesh
  size_t _num_outsides;


  std::unordered_map<std::string, XMLDocument*> _docs;
  std::unordered_map<std::string, XMLElement*> _cur_elements;

public:

  Mesh(SolverPtr solver, Lattice* lattice=nullptr);
  virtual ~Mesh();

  void createLattice(int num_x, int num_y, int num_z=1);
  void createLattice(const PyVector<double> &widths_x,
                     const PyVector<double> &widths_y,
                     const PyVector<double> &widths_z);
  void createLattice(int num_r, double width_r,
                     const PyVector<double> &widths_z,
                     std::string orientation = "y",
                     const Point &offset = Point{0, 0, 0});

  void setLattice(Lattice* lattice);
  std::vector<FP_PRECISION> getMeshData(TallyType type,
                                        bool volume_average=false,
                                        const int group = -1);
  Vector3D getFormattedMeshData(TallyType type, bool volume_average=false,
                                const int group = -1);
  Vector3D getNonUniformFormattedMeshData(std::vector<std::vector<double> >
                                          widths_offsets, TallyType rx,
                                          bool volume_average=false);


  //--------------------------------------------------------------------
  // Utility
  //--------------------------------------------------------------------
  static std::string getPieceName(const std::string &, int);
  void writeXMLDocToFile(std::string doc_name, std::string file);

  /// \brief Print an existing document to string
  /// FIXME: replace string with enumeration
  std::string printDocToString(std::string doc_name = "mesh");

  /// \brief Return the number of points outside the mesh.
  /// \return The number of points.
  size_t getNumOutsides()         { return _num_outsides; }

  /// \brief Reset the number of points outside the mesh.
  /// \param n The number of outside points.
  void resetNumOutsides(size_t n) { _num_outsides = n; }

  /// \brief Increase the number of points outside the mesh by 1.
  /// \return The updated number of points.
  size_t increaseNumOutsides()    { return _num_outsides++; }


  //--------------------------------------------------------------------
  // Dumping mesh data
  //--------------------------------------------------------------------
  /// \brief Dump reaction rates to a given file (interface).
  /// \param file path of the output file
  /// \param types types of tallied mesh data
  /// \param energy_groups selected energy groups (defaults to all)
  void dumpMeshDataToFile(std::string file,
                          std::set<TallyType> types,
                          std::set<int> energy_groups = std::set<int>());

  /// \brief Print mesh data with given types and groups to an XML document
  ///        in memory.
  /// \param file path of the output file
  /// \param types types of tallied mesh data
  /// \param energy_groups selected energy groups (defaults to all)
  void printMeshDataToXML(std::string doc_name,
                          std::set<TallyType> types,
                          std::set<int> energy_groups = std::set<int>());

  /// \brief Initializes an XML document for the unstructured grid.
  /// \param doc_name the name of the XML document
  void initializeUnstructuredDoc(std::string doc_name);

  /// \brief Print mesh data per type-group to an XML document in memory.
  /// \param doc_name name of the XML document
  /// \param data_name name of the data array
  /// \param type The type of reaction to tally
  /// \param volume_average whether the reaction rates should be volume averaged
  /// \param group energy group id starting from 1
  /// \param precision precision of fixed floating-point numbers
  void printSingleMeshDataToXML(std::string doc_name, std::string data_name,
                                TallyType type,
                                bool volume_average=false,
                                const int group = -1,
                                unsigned precision = 8);

  //--------------------------------------------------------------------
  // Dumping tracks
  //--------------------------------------------------------------------
  void initializePPolyDataDoc(std::string doc_name, std::string prefix);
  void printTracksMeshToXML(TrackGeneratorPtr tg, std::string doc_name,
                            std::string data_name, TallyType t);
  void dumpTrackMeshToDir(TrackGeneratorPtr tg, std::set<TallyType> types,
                          std::string dir);
};


//----------------------------------------------------------------------
// Inlining functions
//----------------------------------------------------------------------

/// \brief Construct a piece file name
inline
std::string Mesh::getPieceName(const std::string &prefix, int rank) {
  return prefix + "." + std::to_string(rank) + ".vtp";
}


} // namespace antmoc

#endif
