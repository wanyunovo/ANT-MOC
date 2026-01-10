#include "antmoc/MaterialHandler.h"
#include "antmoc/Geometry.h"
#include "antmoc/log.h"
#include "antmoc/Material.h"

namespace antmoc
{

  MaterialsMap &MaterialHandler::getAllMaterials()
  {
    readAllMaterials();
    return _materials;
  }

  /// \details If the material does not exist, the method will try to
  ///          read it from a file or somewhere else. If it cannot find
  ///          the material anywhere, nullptr will be returned.
  Material *MaterialHandler::getMaterial(std::string name)
  {
    readMaterial(name);
    if (hasMaterial(name))
      return _materials.at(name);
    else
      return nullptr;
  }

  bool MaterialHandler::hasMaterials()
  {
    return !_materials.empty();
  }

  /// \details This method checks the underlying file for a material.
  bool MaterialHandler::hasMaterial(std::string name)
  {
    return _materials.count(name) > 0;
  }

  void MaterialHandler::setMaterial(std::string name, Material *material)
  {
    _materials[name] = material;
  }

  void MaterialHandler::readAllMaterials()
  {
    log::fwarn("Base class MaterialHandler will do nothing during "
               "readAllMaterials");
    return;
  }

  void MaterialHandler::readMaterial(std::string name)
  {
    log::fwarn("Base class MaterialHandler will do nothing during "
               "readMaterial(%s)",
               name);
    return;
  }

  /// \details This method relies on material map of the geometry. If the map
  ///          has not been generated yet, it does nothing.
  /*
  功能：遍历 MaterialHandler 中存储的所有材料，检查它们是否存在于 Geometry 对象的使用列表中。如果某个材料未被几何体使用，则将其从内存中释放并从列表中移除。
  参数：Geometry *geometry —— 指向几何对象的指针，该对象包含了当前模拟中实际用到的所有材料信息
  */
  int MaterialHandler::eraseUnusedMaterials(Geometry *geometry)
  {

    int count = 0;

    if (geometry)
    {
      auto &used_materials = geometry->getAllMaterials();
      // _materials 是一个 std::map<std::string, Material*>，键是材料名，值是材料对象指针
      auto m_it = _materials.begin();

      while (m_it != _materials.end())
      {
        const auto &name = m_it->first;
        const auto &material = m_it->second;

        // used_materials.count(...) 检查这个 ID 是否存在于几何体使用的材料列表中。
        // 如果 count 为 0，说明几何体没用到这个材料。
        if (used_materials.count(material->getId()) == 0)
        {
          log::profile("Cleaning up unused material '{}'", name);

          delete m_it->second;
          m_it = _materials.erase(m_it); // 从 map 中移除该条目，并让迭代器指向下一个元素

          // Counting erased materials
          ++count;
        }
        else
          ++m_it;
      }
    }
    else
      log::warn("Could not erase materials: geometry is nullptr");

    return count; // 返回被删除的材料总数
  }

  void MaterialHandler::setNumEnergyGroups(int num_groups)
  {
    if (num_groups <= 0)
    {
      log::error("Cannot set the number of energy group to a non-positive "
                 "value: {}",
                 num_groups);
    }
    else
    {
      _num_groups = num_groups;
    }
  }

  void MaterialHandler::validateNumEnergyGroups()
  {
    if (getNumEnergyGroups() > getNumEnergyGroupsInFile())
    {
      setNumEnergyGroups(getNumEnergyGroupsInFile());
      log::warn("The number of energy groups ({}) is greater than the value read from "
                "the file ({}), so we reset it to {}",
                getNumEnergyGroups(),
                getNumEnergyGroupsInFile(), getNumEnergyGroupsInFile());
    }
  }

} // namespace antmoc
