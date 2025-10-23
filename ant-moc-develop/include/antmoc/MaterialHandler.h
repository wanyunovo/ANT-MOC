/// \file MaterialHandler.h
/// \brief Handling materials input
/// \date June 13, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \author Gen Wang, USTB

#ifndef MATERIAL_HANDLER_H_
#define MATERIAL_HANDLER_H_

#include <string>
#include <unordered_map>

#include "antmoc/enum_types.h"

namespace antmoc
{

// Forward declarationss
class Geometry;
class Material;

// A map for materials
using MaterialsMap = std::unordered_map<std::string, Material *>;


///---------------------------------------------------------------------
/// \class MaterialHandler
/// \brief Methods for reading materials
///---------------------------------------------------------------------
class MaterialHandler {

protected:

  ///< A map storing retrieved materials
  MaterialsMap _materials;

  ///< Number of energy groups to be read in
  int _num_groups;

  ///< Number of energy groups in the file
  int _num_groups_in_file;

public:

  MaterialHandler() = default;
  virtual ~MaterialHandler() = default;

  //------------------------------------------------------------
  // Interfaces to be implemented by derived classes
  //------------------------------------------------------------
  /// \brief It is assumed that it will do nothing if materials
  ///        don't exist anywhere.
  virtual void readAllMaterials();

  /// \brief It is assumed that it will do nothing if a material
  ///        doesn't exist anywhere.
  virtual void readMaterial(std::string name);

  /// \brief Check the material map
  virtual bool hasMaterials();

  /// \brief Check the existence of a specific material
  /// \param name Name of the material.
  virtual bool hasMaterial(std::string name);

  /// \brief Return the number of retrieved materials
  virtual size_t getNumRetrievedMaterials() { return _materials.size(); }

  /// \brief Return the total number of materials
  virtual size_t getTotalNumMaterials()     { return _materials.size(); }

  //------------------------------------------------------------
  // Interfaces
  //------------------------------------------------------------

  /// \brief Return the map of materials.
  MaterialsMap& getAllMaterials();

  /// \brief Return a material by name.
  /// \param name Name of the material.
  /// \return A pointer to the material, or nullptr if not found.
  Material* getMaterial(std::string name);

  /// \brief Save a material to the handler by name.
  /// \param name Name of the material.
  /// \param material A pointer to the material.
  void setMaterial(std::string name, Material *material);

  /// \brief Get the number of energy group.
  /// \details This value could be set by the user or read from a file.
  int getNumEnergyGroups()                { return _num_groups; }

  /// \brief Set the number of energy group.
  void setNumEnergyGroups(int num_groups);

  /// \brief Get the number of energy group found in the HDF5 file.
  int getNumEnergyGroupsInFile()          { return _num_groups_in_file; }

  /// \brief Check whether the number of energy group is valid and reset it.
  /// \details The number of energy groups must be no more than the value
  ///          read from the HDF5 file.
  void validateNumEnergyGroups();

  //------------------------------------------------------------
  // Clean up
  //------------------------------------------------------------

  /// \brief Clean up all of the unused materials in a geometry.
  /// \return The number of erased materials.
  int eraseUnusedMaterials(Geometry *geometry);

};

} // namespace antmoc

#endif  // MATERIAL_HANDLER_H_
