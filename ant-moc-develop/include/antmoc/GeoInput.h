/// \file include/GeoInput.h
/// \brief Handling geometry input
/// \date June 13, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef GEOINPUT_H_
#define GEOINPUT_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace antmoc
{

// Forward declarationss
class Geometry;
class MaterialHandler;

using MaterialHandlerPtr = std::shared_ptr<MaterialHandler>;


///---------------------------------------------------------------------
/// \class GeoInput
/// \brief The base class of Geometry parsers
///---------------------------------------------------------------------
class GeoInput
{
protected:
  Geometry *_geometry;
  MaterialHandlerPtr _matinput;

  // The Geometry has been read and built
  bool _contained_geometry;

  /// \brief A parameter used to overwrite node 'refines'
  std::vector<int> _global_refines;

  /// \brief A parameter used to overwrite node 'sectors'
  int _global_sectors;

  /// \brief A parameter used to overwrite node 'rings'
  int _global_rings;

public:
  GeoInput(Geometry *geometry = nullptr, MaterialHandlerPtr matinput = nullptr)
    :_geometry(geometry), _matinput(matinput),
     _contained_geometry(false),
     _global_sectors(-1),
     _global_rings(-1) { }

  virtual ~GeoInput() = default;

  // Setters
  void setGeometry(Geometry *geometry);
  void setMaterialHandler(MaterialHandlerPtr matinput);
  void setGlobalRefines(const std::vector<int> &);
  void setGlobalSectors(const int);
  void setGlobalRings(const int);

  // Getters
  Geometry *getGeometry()   { return _geometry; }
  MaterialHandlerPtr getMaterialHandler() { return _matinput; }
  std::vector<int> getGlobalRefines() { return _global_refines; }
  int getGlobalSectors() const { return _global_sectors; }
  int getGlobalRings() const { return _global_rings; }

  // Workers
  virtual void clear() = 0;
  virtual int eraseUnusedPrimitives() = 0;
  virtual int eraseUnusedMaterials() = 0;

  /// \brief Read geometry from a file.
  /// \brief file Path to the file.
  virtual Geometry *readGeometryFromFile(const std::string &file) = 0;

  /// \brief Read geometry primitives from a file.
  /// \brief file Path to the file.
  virtual Geometry *readGlobalPrimitives(const std::string &file) = 0;

  /// \brief Write geometry settings into files under a directory.
  /// \param directory Directory to store geometry files.
  virtual void dumpSettings(const std::string &directory) = 0;

  bool containedMaterials();
  bool containedGeometry() { return _contained_geometry; }

  /// \brief Print the number of energy groups and materials.
  /// \details This method is supposed to be called after the geometry
  ///          initializes its materials so that the map all_materials
  ///          is not empty.
  void printReport();

};

} // namespace antmoc

#endif  // GEOINPUT_H_
