#include "antmoc/MaterialHandler.h"
#include "antmoc/Geometry.h"
#include "antmoc/log.h"
#include "antmoc/Material.h"

namespace antmoc {


MaterialsMap& MaterialHandler::getAllMaterials() {
  readAllMaterials();
  return _materials;
}


/// \details If the material does not exist, the method will try to
///          read it from a file or somewhere else. If it cannot find
///          the material anywhere, nullptr will be returned.
Material* MaterialHandler::getMaterial(std::string name) {
  readMaterial(name);
  if (hasMaterial(name))
    return _materials.at(name);
  else
    return nullptr;
}


bool MaterialHandler::hasMaterials() {
  return !_materials.empty();
}


/// \details This method checks the underlying file for a material.
bool MaterialHandler::hasMaterial(std::string name) {
  return _materials.count(name) > 0;
}


void MaterialHandler::setMaterial(std::string name, Material *material) {
  _materials[name] = material;
}


void MaterialHandler::readAllMaterials() {
  log::fwarn("Base class MaterialHandler will do nothing during "
                      "readAllMaterials");
  return;
}


void MaterialHandler::readMaterial(std::string name) {
  log::fwarn("Base class MaterialHandler will do nothing during "
                      "readMaterial(%s)", name);
  return;
}


/// \details This method relies on material map of the geometry. If the map
///          has not been generated yet, it does nothing.
int MaterialHandler::eraseUnusedMaterials(Geometry *geometry) {

  int count = 0;

  if (geometry) {
    auto &used_materials = geometry->getAllMaterials();
    auto m_it = _materials.begin();

    while (m_it != _materials.end()) {
      const auto &name     = m_it->first;
      const auto &material = m_it->second;

      // If the material ID could not be found in the geometry
      if (used_materials.count(material->getId()) == 0) {
        log::profile("Cleaning up unused material '{}'", name);

        delete m_it->second;
        m_it = _materials.erase(m_it);

        // Counting erased materials
        ++count;
      }
      else
        ++m_it;
    }
  }
  else
    log::warn("Could not erase materials: geometry is nullptr");

  return count;
}


void MaterialHandler::setNumEnergyGroups(int num_groups) {
  if (num_groups <= 0) {
    log::error("Cannot set the number of energy group to a non-positive "
               "value: {}", num_groups);
  }
  else {
    _num_groups = num_groups;
  }
}


void MaterialHandler::validateNumEnergyGroups() {
  if (getNumEnergyGroups() > getNumEnergyGroupsInFile()) {
    setNumEnergyGroups(getNumEnergyGroupsInFile());
    log::warn("The number of energy groups ({}) is greater than the value read from "
              "the file ({}), so we reset it to {}", getNumEnergyGroups(),
              getNumEnergyGroupsInFile(), getNumEnergyGroupsInFile());
  }
}


} // namespace antmoc
