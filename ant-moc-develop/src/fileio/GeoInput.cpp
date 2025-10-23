/// \file src/GeoInput.cpp

#include "antmoc/GeoInput.h"
#include "antmoc/Geometry.h"
#include "antmoc/log.h"
#include "antmoc/MaterialHandler.h"

namespace antmoc
{

/// \brief Set the geometry to be built
void GeoInput::setGeometry(Geometry *geometry) {
  _geometry = geometry;
}

/// \brief Set the material input object
void GeoInput::setMaterialHandler(MaterialHandlerPtr matinput) {
  _matinput = matinput;
}

/// \brief Set the vector of global refines
void GeoInput::setGlobalRefines(const std::vector<int> &refines) {
  _global_refines = refines;
}


void GeoInput::setGlobalSectors(const int sectors) {
  _global_sectors = sectors;
}


void GeoInput::setGlobalRings(const int rings) {
  _global_rings = rings;
}


/// \brief Check if materials exists or not
bool GeoInput::containedMaterials() {
  return (_matinput && _matinput->hasMaterials());
}


//----------------------------------------------------------------------
// Reports
//----------------------------------------------------------------------
void GeoInput::printReport() {
  // Print the number of materials
  log::info("Number of energy groups required      = {}",
             _matinput->getNumEnergyGroups());
  log::info("Number of energy groups in H5 file    = {}",
             _matinput->getNumEnergyGroupsInFile());
  log::info("Total number of materials in H5 file  = {}",
             _matinput->getTotalNumMaterials());
  log::info("Number of retrieved materials         = {}",
             _matinput->getNumRetrievedMaterials());
  log::info("Number of materials in geometry       = {}",
             _geometry->getNumMaterials());
}


} // namespace antmoc
