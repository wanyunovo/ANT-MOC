/// \file src/GeoInputXml.cpp

#include "antmoc/GeoInputXml.h"
#include "antmoc/BoolParser.h"
#include "antmoc/Cell.h"
#include "antmoc/constants.h"
#include "antmoc/file_utils.h"
#include "antmoc/Geometry.h"
#include "antmoc/log.h"
#include "antmoc/Material.h"
#include "antmoc/MaterialHandlerHDF5.h"
#include "antmoc/math_utils.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/Lattice.h"
#include "antmoc/string_utils.h"
#include "antmoc/Surface.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <stdexcept>

namespace antmoc
{


//----------------------------------------------------------------------
// Read the geometry and primitives
//----------------------------------------------------------------------

/// \brief Read geometry from a file
/// \param file geometry file path
Geometry *
GeoInputXml::readGeometryFromFile(const std::string &file) {

  log::info("Reading geometry from file...");

  // Read geometry file and element
  xmlutils::loadFileXML(_doc_geometry, file, true);
  auto geo_xml = _doc_geometry.FirstChildElement("geometry");

  if (!geo_xml) {
    log::error("Failed to read 'geometry' element from file {}", file);
  }

  // Remove unnecessary elements
  auto e = _doc_geometry.FirstChildElement();
  while (e) {
    if (e->Value() != geo_xml->Value())
      _doc_geometry.DeleteChild(e);
    e = e->NextSiblingElement();
  }

  buildGeometry(geo_xml);
  _contained_geometry = true;

  return _geometry;
}


///  \brief   Read extra global primitives from a file
///  \details This method reads primitives from a file, which can be the
///           same file as geometry input. If there is a 'global' element
///           in the file, it will be used to build global primitives.
///  \param   file file path
Geometry *
GeoInputXml::readGlobalPrimitives(const std::string &file) {

  // Read global primitives file and element
  xmlutils::loadFileXML(_doc_primitives, file, false);

  auto global_xml = _doc_primitives.FirstChildElement("global");

  if (global_xml) {
    log::info("Reading geometry primitives from file...");
    buildGlobalPrimitives(global_xml);

    // Remove unnecessary elements
    auto e = _doc_primitives.FirstChildElement();
    while (e) {
      if (e->Value() != global_xml->Value())
        _doc_primitives.DeleteChild(e);
      e = e->NextSiblingElement();
    }

  }
  else if (!stringutils::isSpaces(file)) {
    // If the file is specified but no 'global' was found in it
    log::error("Failed to read 'global' element from file '{}'", file);
  }

  return _geometry;
}


/// \brief Write geometry back to files under the directory.
void GeoInputXml::dumpSettings(const std::string &directory) {

  log::verbose_once("Printing geometry files...");

  fileutils::createDirectory(directory);
  std::string file_geometry   = directory + "/geometry.xml";
  std::string file_primitives = directory + "/primitives.xml";

  if (mpi::isMPIRoot()) {
    // Print geometry file
    xmlutils::saveFileXML(_doc_geometry, file_geometry);

    // Print geometry primitives file as needed
    if (!_doc_primitives.NoChildren())
      xmlutils::saveFileXML(_doc_primitives, file_primitives);
  }

}


/// \brief Build the whole geometry
/// \details This method parses a geometry element from the top.
///          The geometry is represented as an CSG tree.
///          If there is a global element, it will be parsed before
///          parsing the whole geometry.
/// \param root_xml the geometry element
void GeoInputXml::buildGeometry(XMLElement *root_xml)
{
  auto root_cell_xml = queryFirstChild(root_xml, "cell", false);
  auto root_univ_xml = queryFirstChild(root_xml, "universe", false);

  auto global_xml = queryFirstChild(root_xml, "global", false);
  if (global_xml) {
    // if the geometry element has a child global element
    log::info("Found child global element of the geometry, reading primitives...");
    buildGlobalPrimitives(global_xml);
  }

  Universe *root_universe;
  // There must be a root universe.
  if (root_univ_xml) {
    log::verbose_once("Reading the root universe...");
    root_universe = buildUniverse(root_univ_xml);
  }
  // If there isn't a root universe, we have to create a new one.
  else {
    log::verbose_once("Creating the the root universe...");
    root_universe = new Universe();

    if (root_cell_xml) {
      log::verbose_once("Reading the root cell...");
      Cell *root_cell = buildCell(root_cell_xml);
      root_universe->addCell(root_cell);
    }
  }
  _geometry->setRootUniverse(root_universe);
}


//----------------------------------------------------------------------
// Build the geometry and primitives
//----------------------------------------------------------------------

/// \brief Build global geometry primitives
/// \details This method is used to build primitives prior to building
///          the whole geometry. Global primitves can be accessed at any
///          level of XML nodes.
///          For example,
///
///            <global>
///              <surface id="1" type="z-plane" value="+1.0" direction="-1"/>
///              <surface id="2" type="z-plane" value="-1.0" direction="+1"/>
///
///              <cell id="1" name="cell 1" material="UO2" region="1 2"/>
///            </global>
///
///          There are 2 surfaces and 1 cell, which will be inserted into
///          the global surface map and the global cell map respectively.
///
///          For the moment, primitives have to appear in a certain order
///            <global>
///              ...surfaces...
///              ...cells...
///              ...universes...
///              ...lattices...
///            </global>
void GeoInputXml::buildGlobalPrimitives(XMLElement *global_xml) {
  _building_primitives = true;

  // Read primitives
  log::debug("Reading global surfaces from the XML file...");
  extractSurfaceNodes(global_xml, _global_surfaces);

  log::debug("Reading global cells from the XML file...");
  extractCellNodes(global_xml, _global_cells);

  log::debug("Reading global universes from the XML file...");
  extractUniverseNodes(global_xml, _global_universes);

  log::debug("Reading global lattices from the XML file...");
  extractLatticeNodes(global_xml, _global_universes);

  // Lazy computing
  log::debug("Performing lazy computing for global primitives...");
  processLazyCellsParents();
  processLazyLatticesLayouts();
  processLazyCellsFills();

  _building_primitives = false;
}


//----------------------------------------------------------------------
// Utility for parsing
//----------------------------------------------------------------------

/// \brief Replicate a number many times
/// \param number an integer to be replicated
/// \param times the times the integer appears in the result
/// \return a string of repeated integers
std::string GeoInputXml::replicateANumber(int number, size_t times)
{
  std::ostringstream result;
  for (size_t i = 1; i < times; i++)
    result << number << ' ';
  result << number;
  return result.str();
}


///  \brief Substitute replicated strings for sub-strings in the form of "A [N]X"
///  \details A and [N] are assumed to be integers stored in strings. X is at
///          least one non-digit character, such as "R". The string will be
///          extended to a string with A being replicated N times. If N is not
///          specified, A will be replicated once. If there are adjacent "[N]X"s,
///          they are extended to the same A.
///          For example, '9 2R R' will be extended to '9 9 9 9'.
///  \param line a space-splited string to be extended
///  \return an extended string of repeated numbers
std::string GeoInputXml::extendR(std::string line)
{
  std::string result;
  auto words = stringutils::splitString(line);

  auto pos = std::string::npos;
  std::string e = "0123456789";

  result.clear();
  std::string number; // The number to be replicated
  for (auto &s : words) {
    if ((pos = s.find_first_not_of(e)) == std::string::npos) {
      // No letter in s, thus it is assumed to be a number
      number = s;
      result += " " + number;
    }
    else {
      // A letter is found, which will be extended to digits
      if (number == "") {
        log::warn("There must be a number before the string to be extended!");
        return "";
      }
      // Get digits in front of the matched letter
      std::string digits = s.substr(0, pos);
      // The digit is set to 1 if not specified
      int abs_value = (digits.size() == 0)?
                       1 : std::abs(stoi(digits));
      // Replicate and clear the number
      result += " " + replicateANumber(stoi(number), abs_value);
      // Comment it out to allow adjacent symbols
      // number.clear();
    }
  }
  return result.substr(result.find_first_not_of(" "));
}


//----------------------------------------------------------------------
//
// Instantiate primitives
// See GeoInputXMLUtils.cpp
//
//----------------------------------------------------------------------


//----------------------------------------------------------------------
// Build surfaces, cells, universes and lattices
//----------------------------------------------------------------------

/// \brief Extract cell elements and build Cell with exerted region and surfaces
/// \details Exerted surfaces are inherited from the outer Lattice. A cell element
///          may be filled by materials or a universe. A material cell must be indicated
///          by the attribute 'material'. A filled cell must be indicated by the
///          attribute 'fill'. Cells may have child elements to indicate surfaces and
///          universes it owns. A cell can own only 1 universe or lattice. If the
///          attribute 'fill' is set, The universe whose id equals the value of 'fill'
///          must be found in either child elements or the global universe map.
///
///          For example, a simple material cell looks like
///
///            <cell id="1" name="cell 1" material="UO2" region="-1 +2">
///              <surface id="1" type="z-plane" value="+1.0"/>
///              <surface id="2" type="z-plane" value="-1.0"/>
///            </cell>
///
///          The 'region' attribute contains ids of surfaces to be added to the cell.
///          The ids may come from surfaces in the global surface map or surfaces
///          belonging to the current cell.
///
///          Note that attributes 'id', 'name', 'material', 'region', 'temperature',
///          etc. can also be sub-elements.
///
///            <cell id="1" name="cell 1">
///              <material> UO2 </material>
///              <region>   -1 +2 </region>
///              <surface id="1" type="z-plane" value="+1.0"/>
///              <surface id="2" type="z-plane" value="-1.0"/>
///            </cell>
///
///          A simple filled cell looks like
///
///            <cell id="2" name="cell 2" fill="1"/>
///              <lattice id="1" name="local lattice" type="rectangle">
///                <...>
///              </lattice>
///
///          The cell filled with a lattice will take the local lattice in preference
///          to the global one.
///
/// \param cell_xml cell element of the XML file
/// \param exerted_surfaces surfaces inherited from the outer Lattice
/// \return instantiated Cell object
Cell* GeoInputXml::buildCell(XMLElement *cell_xml,
                             std::string exerted_region,
                             SurfaceMap *exerted_surfaces) {

  // Instantiate a Cell object
  auto cell      = newCell(cell_xml);
  auto input_id  = cell->getInputId();
  auto cell_name = cell->getName();

  // Build the region for the cell. Exerted region will be at the very beginning
  // of the cell region string. Surfaces may be from either local or global
  // maps of surfaces.
  buildCellRegion(cell_xml, cell, exerted_region, exerted_surfaces);


  /*读取变化模块数据，并进行元的变换(旋转和平移)*/
  /*这段还需要补充*/
  //XMLElement* transform = cell_xml->FirstChildElement("transform");

  // Parse 'universes' attribute/sub-element. This attrubutes indicates parent
  // universes of the cell. For global primitives, if the specified universes
  // doesn't exist yet, they are probably to be built soon. So that the lazy
  // computing will be used.
  auto universe_ids = queryNodeString(cell_xml, "universes", false, cell_name);

  // For now, only for global primitives. Normal building procedure requires a
  // root universe or a root cell, that is, there is no need to provide the
  // element 'universes' to specify parent universes.
  if (_building_primitives && universe_ids) {
    std::ostringstream lazy;
    std::istringstream ss(universe_ids);

    int id;
    while (ss >> id) {
      auto it = _global_universes.find(id);

      if (it != _global_universes.end()) {
        // if the universe exist, add this cell to it
        it->second->addCell(cell);
      }
      else {
        // if the universe doesn't exist, add the id to a list
        lazy << id << ' ';
      }
    }

    if (lazy.str() != "") {
      // if there is any undefined universes
      _lazy_cells_parents.insert({cell->getInputId(), lazy.str()});
    }
  }
  else if (!_building_primitives && universe_ids) {
    log::ferror("For now, attribute 'universes' is supposed to "
                "be used only in building global primitives");
  }

  // Parse 'fill' and 'material' attributes. The cell should be filled by
  // either materials or a universe. Empty cell is not allowed.
  auto fill     = queryNodeInt(cell_xml, "fill", false);
  auto mat_name = queryNodeString(cell_xml, "material", false, cell_name);

  bool is_fill = (fill != defaultQueriedInt());
  bool is_material = mat_name;

  if (!is_fill && !is_material) {
    log::error("Failed to build cell '{}'(id={}): not filled by anything",
               cell_name, input_id);
  }
  else if (is_fill && is_material) {
    log::error("Failed to build cell '{}'(id={}): cannot has 'fill' "
               "and 'material' at the same time", cell_name, input_id);
  }
  else if (!is_fill && is_material) { // this is a material cell

    // Set the number of sectors and rings
    int num_sectors = queryNodeInt(cell_xml, "sectors", false, cell_name);
    int num_rings   = queryNodeInt(cell_xml, "rings", false, cell_name);

    if (num_sectors == defaultQueriedInt())
      num_sectors = 0;
    else if (_global_sectors >= 0)
      num_sectors = _global_sectors;

    if (num_rings == defaultQueriedInt())
      num_rings = 0;
    else if (_global_rings >= 0)
      num_rings = _global_rings;

    cell->setNumSectors(num_sectors);
    cell->setNumRings(num_rings);

    // Check whether materials exists
    auto material = getMaterial(mat_name);
    if (!material) {
      log::error("Failed to build cell '{}'(id={}): requires an undefined "
                  "material named '{}'", cell_name, input_id, mat_name);
    }
    cell->setFill(material);

    // Debugging
    log::debug("Filling in cell '{0}'(id={1}) with material '{2}'(id={3}), "
               "this cell will be refined with {4} ring(s) and {5} sector(s)",
               cell_name, input_id, material->getName(), material->getId(),
               num_rings, num_sectors);
  }
  else {
    // The cell is filled with a universe.
    // It is assumed that the universe needed by this cell is either in the
    // global universe map or among sub-nodes of the cell.
    // We have to read all universes of this cell to determine if the XML
    // file is ill-formed.
    UniverseMap universes;
    extractUniverseNodes(cell_xml, universes);
    extractLatticeNodes(cell_xml, universes);

    bool in_local = universes.count(fill);
    bool in_global = _global_universes.count(fill);

    if (!in_local && !in_global) {  // cannot find the universe/lattice
      if (!_building_primitives) {
        log::error("Cell (id={}) requires an undefined universe "
                   "(id={})", cell->getInputId(), fill);
      }
      // lazy computing
      else {
        _lazy_cells_fills.insert({cell->getInputId(), fill});
        // Debugging
        log::debug("Filling in cell '{}'(id={}) with universe {} (lazy)",
                   cell_name, input_id, fill);
      }
    }
    // found the universe/lattice in local or global maps
    else {
      // Prefer the local universe
      Universe *u = (in_local ? universes.at(fill) : _global_universes.at(fill));
      cell->setFill(u);

      // Debugging
      log::debug("Filling in cell '{}'(id={}) with universe '{}'(id={},{})",
                 cell_name, input_id, u->getName(), fill, u->getId());
    }
  }
  return cell;
}


/// \brief Read an element from the XML file and build a Universe
/// \details A universe is simply a container of cells. Therefore, universes
///          can be explicit or implicit. As for now, implicit universes
///          can only be defined in global primitives.
///          For example, a simple definition of a universe looks like
///
///            <universe id='1' name='Universe 1'/>
///
///          To specify its child cells, try
///
///            <universe id='2' name='Universe 2'>
///              <cell id='10' name='a' material='A' region='1 2'/>
///              <cell id='11' name='b' material='B' region='1 2'/>
///            </universe>
///
///          Cells appeared inside a universe will be added to this universe.
///          To take global cells, try
///
///            <universe id='3' name='Universe 3'>
///              <cells> 1 2 3 </cells>
///            </universe>
///          or
///            <universe id='3' name='Universe 3' cells='1 2 3'/>
///
///          This style of definition can be mixed with the above style.
///
/// \param univ_xml the 'universe' element
/// \param exerted_surfaces the surfaces to be exerted
/// \return a dynamic Universe object
Universe*
GeoInputXml::buildUniverse(XMLElement *univ_xml,
                           std::string exerted_region,
                           SurfaceMap *exerted_surfaces) {

  auto univ     = newUniverse(univ_xml);
  auto input_id = univ->getInputId();
  auto name     = univ->getName();

  // Create a buffer to store extracted cells, which are belonging to this
  // universe. Both the XML element and the global map will be searched.
  CellMap cells;
  extractCellNodes(univ_xml, cells, exerted_region, exerted_surfaces);
  readXmlObjCells(univ_xml, cells);

  // If there are no cells, throw an error
  if (cells.size() == 0) {
    log::ferror("Empty universe is not permitted: id " "= %d, name = %s",
                 input_id, name);
  }

  for (auto &c : cells)
    univ->addCell(c.second);

  return univ;
}


/// \brief Read an element from the XML file and build a Lattice
/// \details This method is designed for both rectangle and hexagon lattice.
///          Originally, a lattice may have attributes id, name, type and offset.
///
///            <lattice id="1" name="L1" type="rectangle" offset="1 1 1"/>
///
///          If the offset is omitted, "0 0 0" will be used.
///
///          Type : RecLattice
///
///          A non-trivial lattice must contain a layout, either 2D or 3D.
///          Basically, a layout is built according to a set of elements
///          belonging to the lattice. That is, 'widths', 'refines', and 'row's.
///
///            <lattice id="2" name="L2" type="rectangle">
///              <widths x="[0.1] + 0.1" y="[0.2] * 2" z="[0.3] * 2"/>
///
///              <row>  1  1 </row>
///              <row>  2  2 </row>
///
///              <row>  1  2 </row>
///              <row>  2  1 </row>
///
///            </lattice>
///
///          The example above shows how to write a simple 3-D layout.
///          The layout of 'row's is written to resemble the lattice itself.
///          That is,
///            - the first 2 'row's represent the top view of the top z-section
///            - the last 2 'row's represent the top view of the bottom z-section
///
///          The number of lattice cells is 2*2*2=8. Note that this is a
///          non-uniform lattice, which indicates a bounding box of
///          0.2*0.4*0.6 cm^3.
///          (Note that if the 'widths' element contains expressions, square
///          brackets '[]' are the only way you could use to specify a vector.
///          This feature resembles 'list' of Python.)
///
///          Note that the numbers in a 'row' element are ids of universes, which
///          can be either global or local. If 'row's are omitted, ids should be
///          inside a 'universes' element.
///
///              <universes>
///                1  1
///                2  2
///
///                1  2
///                2  1
///              </universes>
///
///          Another facility is called 'refines', which refines the lattice with
///          the layout unchanged. If a lattice is refined, the number of lattice
///          cells and the widths of the lattice will change. For example,
///
///            <lattice id="3" name="L3" type="rectangular">
///              <widths x="[0.1] * 2" y="[0.2] * 2" z="[0.3] * 2"/>
///              <refines z="2"/>
///
///              <row>  1  1 </row>
///              <row>  2  2 </row>
///
///              <row>  1  2 </row>
///              <row>  2  1 </row>
///
///            </lattice>
///
///          The 'refines' element requires only a z-refine, which will give a
///          lattice with 2*2*(2*2)=16 lattice cells. The width in the z direction
///          becomes 0.3/2=0.15.
///
///          If there is a 'format' element, region and surfaces included in the
///          element will be exerted on each universe of the lattice, that is, on
///          each cell of the child universes of the lattice.
///
///          There are several means to save layout typing.
///
///          1. Repeated rows
///
///            <row repeat="2">  1  1 </row>
///          gives
///            <row>  1  1 </row>
///            <row>  1  1 </row>
///
///          2. Extended rows
///
///            <row>  1  2R  2 </row>
///          gives
///            <row>  1  1  1  2 </row>
///
///          3. Extruded layout
///             This is used to extend a 2-D layout.
///
///            <lattice ... axial="extruded">
///              <widths x="[0.1] * 2" y="[0.2] * 2" z="[0.3] * 2"/>
///              <row>  1  1 </row>
///              <row>  1  1 </row>
///            </lattice>
///          equals
///            <lattice ... >
///              <widths x="[0.1] * 2" y="[0.2] * 2" z="[0.3] * 2"/>
///              <row>  1  1 </row>
///              <row>  1  1 </row>
///
///              <row>  1  1 </row>
///              <row>  1  1 </row>
///            </lattice>
///
///          Type : HexLattice
///
///          A hexagon lattice has different attributes/elements from
///          rectangle lattices. For now, 'n_rings' is used to specify
///          the number of lattice cells along radius, 'r' and 'z'
///          are used to specify the pitches. The non-uniform feature
///          is implemented for z-axis.
///
///            <lattice id='1' type='Hexagon' n_rings='2'>
///              <widths r='0.2' z='[1.4] * 2'/>
///              <universes>
///                    1
///                  2   2
///                    3
///                  2   2
///                    1
///              </universes>
///            </lattice>
///
///          The layout of a hexagon lattice can be only written in
///          a 'universes' element/attribute.
///
/// \param lat_xml the 'lattice' element
/// \return a dynamic Lattice object
Lattice* GeoInputXml::buildLattice(XMLElement *lat_xml) {
  // Immediately return if a null pointer is passed in
  if (lat_xml == nullptr)
    return nullptr;

  // Instantiate the lattice
  auto lattice = newLattice(lat_xml);

  // Get region and surfaces elements to be exerted on child cells.
  // To achieve this, a 'region' attribute/element must defined for
  // the 'format' element. The exerted surfaces can be either global
  // or local to the 'format' element.
  auto format_xml = queryFirstChild(lat_xml, "format", false, lattice->getName());
  std::string exerted_region;
  SurfaceMap surfaces;
  if (format_xml) {
    auto region = queryNodeString(format_xml, "region", false, lattice->getName());
    if (region) {
      exerted_region = std::string(region);
      extractSurfaceNodes(format_xml, surfaces);
      readXmlObjSurfaces(format_xml, surfaces);
    }
  }


  // Read, set and check the layout of the lattice
  auto layout = buildLatticeLayout(lat_xml, lattice);
  buildLatticeRefines(lat_xml, lattice, layout);

  // Extract and instantiate child universes (local)
  UniverseMap unique_universes;
  extractUniverseNodes(lat_xml, unique_universes, exerted_region, &surfaces);

  // Map these universes to lattice cells
  mapUniversesByLayout(lattice, unique_universes, layout);

  return lattice;
}


/// \brief Extract lattice elements and instantiate Lattice
/// \details Lattices and Universes share unique ids and are stored in the
///          same map.
/// \param elem_xml the father element
void
GeoInputXml::extractLatticeNodes(XMLElement *elem_xml,
                                 UniverseMap &universes) {
  auto lattice_xml = elem_xml->FirstChildElement("lattice");

  while (lattice_xml){

    // Each global lattice will be built only once
    Universe *u = nullptr;
    auto id = queryNodeInt(lattice_xml, "id", false);
    if (_building_primitives && _global_universes.count(id)) {
      log::ferror("Detected duplicated global lattice (id=%d)", id);
    }
    else {
      // Instantiate lattices and exert widths on their child universes
      u = buildLattice(lattice_xml);
    }

    if (!u->anonymous()) {
      if (id != u->getInputId()) {
        log::ferror("Lattice id (%d) changed after it is built (%d)",
                     id, u->getInputId());
      }

      // Insert into the global map on demand
      if (_building_primitives) {
        auto result = _global_universes.insert({id, u});
        if (!result.second) {
          log::ferror("Failed to insert lattice %d into the "
                      "global map: duplicated id", id);
        }
      }
    }
    else {
      id = u->getId(); // get the unique id
    }

    // Insert into the result map
    if (&universes != &_global_universes) {
      auto result = universes.insert({id, u});
      if (!result.second)
        log::fwarn_once("Failed to insert lattice %d: id exists", id);
    }

    lattice_xml = lattice_xml->NextSiblingElement("lattice");
  }
}


/// \brief Extract universe elements and instantiate Universe
/// \param elem_xml the father element
void
GeoInputXml::extractUniverseNodes(XMLElement *elem_xml,
                                  UniverseMap &universes,
                                  std::string exerted_region,
                                  SurfaceMap *exerted_surfaces) {
  auto univ_xml = elem_xml->FirstChildElement("universe");

  while (univ_xml){

    // Each global universe will be built only once
    Universe *u = nullptr;
    auto id = queryNodeInt(univ_xml, "id", false);

    if (_building_primitives && _global_universes.count(id)) {
      log::ferror("Detected duplicated global universe (id=%d)", id);
    }
    else {
      // Instantiate the universe and exert a region on its child cells
      u = buildUniverse(univ_xml, exerted_region, exerted_surfaces);
    }

    if (!u->anonymous()) {
      if (id != u->getInputId()) {
        log::ferror("Universe id (%d) changed after it is built (%d)",
                     id, u->getInputId());
      }

      // Insert into the global map on demand, not for anonymous universe
      if (_building_primitives) {
        auto result = _global_universes.insert({id, u});
        if (!result.second) {
          log::ferror("Failed to insert universe %d into the "
                      "global map: duplicated id", id);
        }
      }
    }
    else {
      id = u->getId(); // get the unique id
    }

    // Insert into the result map
    if (&universes != &_global_universes) {
      auto result = universes.insert({id, u});
      if (!result.second) {
        log::fwarn_once("Failed to insert universe %d: id exists", id);
      }
    }

    univ_xml = univ_xml->NextSiblingElement("universe");
  }
}


/// \brief Extract cell nodes and instantiate Cell objects
/// \details XML nodes are extracted by this method. After instantiated,
///          a cell will be inserted into the out SurfaceMap. If we are
///          building global primitives, the cell will be inserted into
///          the global cell map at the same time.
/// \param elem_xml the father element
/// \param cells a map to store return values
void
GeoInputXml::extractCellNodes(XMLElement *elem_xml,
                              CellMap &cells,
                              std::string exerted_region,
                              SurfaceMap *exerted_surfaces) {
  if (!elem_xml) {
    log::ferror("Unable to extract cells from nullptr");
  }

  auto cell_xml = queryFirstChild(elem_xml, "cell", false);

  while (cell_xml){

    // When building global primitives, each cell will be built only once
    Cell *c = nullptr;
    auto id = queryNodeInt(cell_xml, "id", false);

    if (_building_primitives && _global_cells.count(id))
      log::ferror("Detected duplicated global cell (id=%d)", id);
    else // Instantiate the cell and exert a region on its child cells
      c = buildCell(cell_xml, exerted_region, exerted_surfaces);

    if (!c->anonymous()) {
      if (id != c->getInputId())
        log::ferror("Cell id (%d) changed after it is built (%d)",
                     id, c->getInputId());

      // Insert into the global map on demand
      if (_building_primitives) {
        auto result = _global_cells.insert({id, c});
        if (!result.second) {
          log::ferror("Failed to insert cell %d into the "
                      "global map: duplicated id", id);
        }
      }
    }
    else {
      id = c->getId(); // get the unique id
    }

    // Insert into the result map
    if (&cells != &_global_cells) {
      auto result = cells.insert({id, c});
      if (!result.second) {
        log::fwarn_once("Failed to insert cell %d: id exists", id);
      }
    }

    cell_xml = cell_xml->NextSiblingElement("cell");
  }
}


/// \brief Read global cells from the 'cells' attribute/sub-element
/// \details Global cells are written in a 'cells' element
///          as plain text. This method reads and parses the plain
///          text and insert any valid cell into a map.
/// \param elem_xml the father element
void
GeoInputXml::readXmlObjCells(XMLElement *elem_xml, CellMap &cells) {
  if (!elem_xml) {
    log::ferror("Unable to extract cells from nullptr");
  }

  auto cells_ids = queryNodeString(elem_xml, "cells", false);

  if (cells_ids) {
    std::istringstream ss(cells_ids);
    int id;
    while (ss >> id) {
      auto it = _global_cells.find(id);
      if (it != _global_cells.end()) {
        auto result = cells.insert(*it);
        if (!result.second) {
          log::fwarn_once("Failed to insert cell %d: id exists", it->first);
        }
      }
      else
        log::ferror("Unable to find global cell %d", id);
    }
  }
}


/// \brief Extract surface elements
/// \details Surfaces will be instantiated and returned
/// \param elem_xml the father element
/// \return a map contained surfaces
void
GeoInputXml::extractSurfaceNodes(XMLElement *elem_xml,
                                 SurfaceMap &surfaces) {
  if (!elem_xml) {
    log::error("Unable to read 'surface' element from nullptr");
  }

  auto sf_xml = elem_xml->FirstChildElement("surface");

  // Extract surfaces from the XML element
  while (sf_xml)
  {
    auto id = queryNodeInt(sf_xml, "id", false);

    // If we are building primitives, the id should be unique.
    if (_building_primitives && _global_surfaces.count(id))
      log::error("Detected duplicated global surface (id={})", id);

    auto surface = newSurface(sf_xml);

    // If we are building primitives, every surface is made global and unique
    // even if the surface is a sub-element of some other primitive.
    if (_building_primitives && (&surfaces != &_global_surfaces)) {
      auto result = _global_surfaces.insert({surface->getId(), surface});
      if (!result.second) {
        log::error("Failed to insert surface {} into the global map: duplicated id", surface->getId());
      }
    }

    // Insert into the result map
    auto result = surfaces.insert({surface->getId(), surface});
    if (!result.second)
      log::warn_once("Failed to insert surface {}: id exists", surface->getId());

    sf_xml = sf_xml->NextSiblingElement("surface");
  } // end while
}


/// \brief Extract global surfaces from the 'surfaces' attribute/sub-element
/// \details Global surfaces are written in a 'surfaces' element
///          as plain text. This method reads and parses the plain
///          text and insert any valid surface into a map.
/// \param elem_xml the father element
void
GeoInputXml::readXmlObjSurfaces(XMLElement *elem_xml,
                                SurfaceMap &surfaces) {
  if (!elem_xml) {
    log::ferror("Unable to read 'surfaces' from nullptr");
  }

  auto surface_ids = queryNodeString(elem_xml, "surfaces", false);

  if (surface_ids) {
    std::istringstream ss(surface_ids);
    int id;
    while (ss >> id) {
      auto it = _global_surfaces.find(id);
      if (it != _global_surfaces.end()) {
        auto result = surfaces.insert(*it);
        if (!result.second) {
          log::fwarn_once("Failed to insert global surface %d: "
                          "id exists", it->first);
        }
      }
      else
        log::ferror("Unable to find global surface %d", id);
    }
  }
}


/// \brief Build a region and add surfaces to a cell
/// \details The region is represented by a string. A region can be either
///          passed from the outer lattice or defined inside a cell.
///          Exerted regions will be considered first. However, exerted surfaces
///          could be overrided by local ones.
///          <lattice ... >
///            <format region='-1 +2'>
///              <surface id='1' ... >
///              <surface id='2' ... >
///            </format>
///            <universe ... >
///              <cell id='1' region='+1 -3' ...>
///                <surface id='1' ... >
///                <surface id='3' ... >
///              </cell>
///            </universe>
///          </lattice>
///
///          The case above defined a cell with region '-1 +2 +1 -3' and surfaces
///          1 (local one), 2, 3.
void GeoInputXml::buildCellRegion(XMLElement *cell_xml, Cell *cell,
                                  const std::string &exerted_region,
                                  const SurfaceMap *exerted_surfaces) {

  // Initialize the cell with exerted region and surfaces
  std::string region = exerted_region;
  SurfaceMap surfaces;
  if (exerted_surfaces)
    surfaces = *exerted_surfaces;

  // Extract surface elements from cell elements.
  // This will give us a map of local surfaces
  SurfaceMap local_sfmap;
  extractSurfaceNodes(cell_xml, local_sfmap);

  // Read attribute 'region', which contains ids of surfaces
  auto region_lcl = queryNodeString(cell_xml, "region", false, cell->getName());

  if (region_lcl) {
    // Look up surfaces in local and global maps.
    // Local surfaces have higher priority than the global ones. If two
    // surfaces in different maps have the same id, pick the local one.
    region = region + " " + region_lcl;  // append the local region

    auto surface_ids = getSurfaceIdsFromRegion(region);
    findSurfacesByRegion(surface_ids, local_sfmap, surfaces);
    findSurfacesByRegion(surface_ids, _global_surfaces, surfaces);

    if (surface_ids.size() > surfaces.size()) {
      log::ferror("Undefined ids of surfaces used when parsing "
                  "region '%s' of cell %d", region, cell->getInputId());
    }
    else if(surface_ids.size() < surfaces.size()) {
      log::fwarn("Unused surfaces detected when parsing region '%s' of "
                  "cell %d", region, cell->getInputId());
    }

    // Add the region and associated surfaces to the cell.
    // If region is empty, surfaces won't be added to the cell as well.
    cell->addRegionSurfaces(region, surfaces);

  }
  else if (!local_sfmap.empty()) {  // extra local surfaces detected
    log::fwarn_once("Cell (id=%d) has local surfaces but no "
               "'region', maybe it is missing?");
  }

  if (!stringutils::isSpaces(region)) {
    // Debugging
    log::fdebug("Building region for cell '%s'(id=%d): '%s'",
                 cell->getName(), cell->getInputId(), region);
  }
}


/// \brief Parse the region str and store the tokens in a PyVector
PyVector<int>
GeoInputXml::getSurfaceIdsFromRegion(const std::string &region) {

  BoolParser parser;
  auto surface_ids = parser.tokenize(region);
  parser.removeOperators(surface_ids);

  return surface_ids;
}


/// \brief Find all surfaces in maps according to region string
/// \details This method is used to find surfaces among either
///          local or global surfaces. Surfaces will be appended
///          to the out map.
/// \param region a string indicating the region of a cell
/// \param source a map of surfaces from which to read surfaces
/// \param surfaces a map to store surfaces
/// \return the number of found surfaces
int
GeoInputXml::findSurfacesByRegion(const PyVector<int> &surface_ids,
                                  const SurfaceMap &source,
                                  SurfaceMap &surfaces) {
  int count = 0;

  for (auto id : surface_ids) {
    auto it = source.find(std::abs(id));
    if (it != source.end()) {
      // If the id is found, save the surface as return value
      surfaces.insert(*it);
      count++;
      // auto result = surfaces.insert(*it);
      // if (!result.second) {
      //   log::fwarn_once("Failed to insert existing surface (id=%d). "
      //              "The ids of local and global surfaces may conflict. "
      //              "Current size of region = %d", it->first, surfaces.size());
      // }
    }
  }
  return count;
}


/// \brief Process the lazy cells parents
/// \details This method loops over pairs of the lazy cells map and determine
///          if the mentioned universes are defined. If there is an undefined
///          universe, the method will instantiate it and add it to the global
///          universes.
void GeoInputXml::processLazyCellsParents() {
  for (auto &c_u : _lazy_cells_parents) {
    auto cell = _global_cells.at(c_u.first); // the cell must exist
    std::istringstream ss(c_u.second);
    int id;
    while (ss >> id) {
      Universe *u;
      auto it = _global_universes.find(id);
      if (it != _global_universes.end())
        u = it->second; // if the universe exist, add this cell to it
      else {  // instantiate a universe
        u = new Universe(id, cell->getName()); // use the cell's name
        _global_universes.insert({id, u});
      }
      u->addCell(cell);
    }
  }
  _lazy_cells_parents.clear();
}


/// \brief Process the lazy cells fills
/// \details This method must be called after building all universes
///          because it doesn't create any new universe.
void GeoInputXml::processLazyCellsFills() {
  for (auto &c_u : _lazy_cells_fills) {
    auto cell_id = c_u.first;
    auto cell = _global_cells.at(cell_id); // the cell must exist
    auto fill = c_u.second;

    // If the universe exists, fill the cell with it
    auto it = _global_universes.find(fill);
    if (it != _global_universes.end())
      cell->setFill(it->second);
    else {
      log::ferror("Cell (id=%d) requires an undefined universe "
                  "(id=%d)", cell_id, fill);
    }
  }
  _lazy_cells_fills.clear();
}


/// \brief Process the lazy lattice layouts
/// \details A lattice and its layout are processed in this method. which
///          means all universes used by this lattice must exist before
///          this method. These universes will be inserted into the lattice
///          according to the lattice layout.
void GeoInputXml::processLazyLatticesLayouts() {

  for (auto &e : _lazy_lat_layouts) {
    auto lat_id = e.first;
    auto u = _global_universes.at(e.first);
    auto lattice = dynamic_cast<Lattice*>(u);  // it must be a lattice
    if (!lattice) {
      log::ferror("Lattice (id=%d) not found for lazy computing", lat_id);
    }

    auto &layout = e.second;

    // Get all existing universes and map them again
    auto unique_universes = lattice->getUniqueUniverses();
    mapUniversesByLayout(lattice, unique_universes, layout);
  }
  _lazy_lat_layouts.clear();
}


/// \brief Retrieve widths of lattice cells
/// \param width_x retrieved x-direction width
/// \param width_y retrieved y-direction width
/// \param width_z retrieved z-direction width
void
GeoInputXml::readLatticeCellWidths(XMLElement *elem_xml,
                                   WidthVec &widths_x,
                                   WidthVec &widths_y,
                                   WidthVec &widths_z) {
  if (!elem_xml) {
    log::ferror("Unable to extract lattice widths from nullptr");
  }

  auto wx_str = queryNodeString(elem_xml, "x", true);
  auto wy_str = queryNodeString(elem_xml, "y", true);
  auto wz_str = queryNodeString(elem_xml, "z", false);

  // Parse widths along x-axis and y-axis
  widths_x = _parser.parse(wx_str);
  widths_y = _parser.parse(wy_str);

  // Widths_z will be parsed if it exists
  if (wz_str)
    widths_z = _parser.parse(wz_str);
  else
    widths_z = WidthVec();  // An empty vector
}


/// \brief Read the axial pitch and z widths for HexLattice
void
GeoInputXml::readLatticeCellWidths(XMLElement *elem_xml,
                                   double &width_r,
                                   WidthVec &widths_z) {
  width_r = queryNodeDouble(elem_xml, "r", true);

  // Widths_z will be parsed if it exists
  auto wz_str = queryNodeString(elem_xml, "z", false);
  if (wz_str)
    widths_z = _parser.parse(wz_str);
  else
    widths_z = WidthVec();  // An empty vector
}


/// \brief Check exerted, computed and extracted lattice widths
/// \details This method must be called after the lattice widths
///          and lattice cells are set.
void
GeoInputXml::checkLatticeWidths(const Lattice *lattice,
                                std::array<double,3> box) {
  if (box[0] <= 0. || box[1] <= 0. || box[2] <= 0.) {
    log::ferror("Negative widths detected when building lattice "
                "'%s': box = (%f, %f, %f)", lattice->getName(),
                box[0], box[1], box[2]);
  }

  auto wx = lattice->getMaxX() - lattice->getMinX();
  auto wy = lattice->getMaxY() - lattice->getMinY();
  auto wz = lattice->getMaxZ() - lattice->getMinZ();

  auto gt = [](double a, double b)
    { return std::abs(a - b) >= FLT_EPSILON && a > b; };

  if ( gt(wx, box[0]) || gt(wy, box[1]) || gt(wz, box[2]) ) {
    log::ferror("The widths of lattice '%s' are (%f, %f, %f), "
                "but the bounding box from the outer cell is "
                "(%f, %f, %f)",
                lattice->getName(),
                wx, wy, wz,
                box[0], box[1], box[2]);
  }
}


/// \brief Parse the type of boundary conditions
boundaryType GeoInputXml::parseBoundaryType(std::string bc) {

  // Transform the lower character to the upper character
  stringutils::trim(bc);
  stringutils::toUpper(bc);
  boundaryType type = VACUUM;

  if      (bc == "VACUUM")     { type = VACUUM; }
  else if (bc == "REFLECTIVE") { type = REFLECTIVE; }
  else if (bc == "PERIODIC")   { type = PERIODIC; }
  else if (bc == "INTERFACE")  { type = INTERFACE; }
  else {
    log::ferror("Undefined boundary type: %s", bc);
  }

  return type;
}


/// \brief Map unique universes to Lattice layout
/// \details The layout is saved as the same shape as it looks like.
///          That is, the first universe is at the lower left corner.
///                   x
///                 0 1 2
///                ------
///             2 | 7 8 9
///           y 1 | 4 5 6
///             0 | 1 2 3
void
GeoInputXml::mapUniversesByLayout(Lattice *lattice,
                                  UniverseMap &unique_universes,
                                  const LayoutVec &layout) {
  checkLatticeLayout(lattice, layout);

  // Get the number of lattice cells
  int lat_id = lattice->getInputId();
  size_t num_univs = lattice->getNumLatticeCells();
  auto universes = new Universe *[num_univs];

  // Check if this is during lazy computing
  auto &_lazy = _lazy_lat_layouts;
  bool during_lazy = _lazy.find(lat_id) != _lazy.end();

  // For lazy computing
  bool be_lazy = false;

  for (size_t i = 0; i < num_univs; ++i) {
    // Find the universe to be mapped
    auto id = layout[i];

    auto u = getUniverseFromMaps(unique_universes,
                                _global_universes,
                                id);
    if (!u) {
      if (!_building_primitives) {
        log::ferror("Universe (id=%d) cannot be found in local "
                    "and global universe maps", id);
      }
      else if (during_lazy) {
        log::ferror("Lattice (id=%d) requires an undefined universe (id=%d)",
                    lat_id, id);
      }
      else { // enable lazy computing
          be_lazy = true;
      }
    }
    // the layout starts from the top left corner
    universes[i] = u;
  }

  // Determine whether to enable lazy computing
  if (be_lazy && !during_lazy)
    _lazy.insert({lat_id, layout});

  // Fill the lattice with universes
  auto rec_lat = dynamic_cast<RecLattice*>(lattice);
  if (rec_lat) {
    int num_x = rec_lat->getNumX();
    int num_y = rec_lat->getNumY();
    int num_z = rec_lat->getNumZ();
    rec_lat->setUniverses(num_z, num_y, num_x, universes);
  }
  else {
    auto hex_lat = dynamic_cast<HexLattice*>(lattice);
    int num_r = hex_lat->getNumR();
    int num_z = hex_lat->getNumZ();
    hex_lat->setUniverses(num_z, num_r, universes);
  }

  delete [] universes;
}


/// \brief Get a universe from either local or global map
/// \details If the universe is found in localmap, it will be returned.
///         Otherwise, the global map will be searched.
/// \param local local map
/// \param global global map
Universe *
GeoInputXml::getUniverseFromMaps(const UniverseMap &local,
                                 const UniverseMap &global,
                                 const int id) {
  Universe *u = nullptr;

  auto it = local.find(id);
  if (it != local.end()) // among the local universes
    u = it->second;
  else { // among the global universes
    auto it2 = global.find(id);
    if (it2 != global.end())
      u = it2->second;
  }
  return u;
}


/// \details The layout consists of contiguous elements representing the ids of
///          universes. There are two means to input a layout:
///
///          1. From rows, which means each row of the layout should contain
///             the same number of lattice cells.
///
///             <row> 1 2 1 1 1 1 </row>
///             <row> 2 1 1 1 1 1 </row>
///             <row> 2 1 1 1 1 1 </row>
///
///             In this case, rows can be 'compressed' in 2 ways.
///               <row> 1 2 1 3R </row>
///               <row repeat='2'> 2 1 4R </row>
///
///          2. From 'universes' element. In this case, we can only check the
///             total number of lattice cells rather than the number of lattice
///             cells along a specified direction.
///             Only 'R' compression takes effect in this case.
LayoutVec GeoInputXml::buildLatticeLayout(XMLElement *lattice_xml, Lattice *lattice) {

  if (!lattice_xml)
    log::error("Cannot build a layout for an empty element");

  LayoutVec layout;

  if (xmlutils::existNode(lattice_xml, "universes")) {
    // Build from 'universes'
    try {
      layout = newLayoutFromUniverses(lattice_xml, "universes");
    }
    catch (std::exception &e) {
      log::error("An error occurred when reading 'universes' for lattice '{}'(id={}): {}",
                  lattice->getName(), lattice->getId(), e.what());
    }
  }
  else if (dynamic_cast<RecLattice*>(lattice)) {
    // Check the alternative for RecLattice
    try {
      layout = newLayoutFromRows(lattice_xml);
    }
    catch (std::exception &e) {
      log::error("An error occurred when reading 'row's for lattice '{}'(id={}): {}",
                  lattice->getName(), lattice->getId(), e.what());
    }
  }

  // Determine if it is an extruded lattice and then build extruded layout for it.
  buildExtrudedLayout(lattice_xml, lattice, layout);

  return layout;
}


/// \brief Set the offset in global coordinates
void GeoInputXml::setLatticeOffset(Lattice *lattice,
                                   std::string offset) {
  double x, y, z;
  std::istringstream digits(offset);
  digits >> x >> y >> z;

  lattice->setOffset(x, y, z);

  log::debug("Set offset for lattice '{}'(id={}) = ({:.2f}, {:.2f}, {:.2f})",
             lattice->getName(), lattice->getId(), x, y, z);
}


/// \brief Set widths and the number of lattice cells for a RecLattice
void GeoInputXml::setLatticeCells(RecLattice *lattice,
                                  WidthVec &widths_x,
                                  WidthVec &widths_y,
                                  WidthVec &widths_z) {
  lattice->setNumX(widths_x.size());
  lattice->setNumY(widths_y.size());
  lattice->setNumZ(widths_z.size());

  // Contain axial refines
  lattice->setWidths(widths_x, widths_y, widths_z);
  lattice->computeSizes();
}


/// \brief Read and set lattice cells for a HexLattice
void GeoInputXml::setLatticeCells(HexLattice *lattice,
                                  int num_r,
                                  double width_r,
                                  WidthVec &widths_z) {
  lattice->setNumR(num_r);
  lattice->setNumZ(widths_z.size());

  lattice->setWidths(width_r, widths_z);
  lattice->computeSizes();
}


/// \brief Read and set lattice refines
void GeoInputXml::buildLatticeRefines(XMLElement *lat_xml,
                                      Lattice *lattice,
                                      LayoutVec &layout) {

  if (!lattice) {
    log::error("Unable to set refines for a null lattice pointer");
  }

  auto rec_lat = dynamic_cast<RecLattice*>(lattice);
  auto hex_lat = dynamic_cast<HexLattice*>(lattice);

  if (rec_lat) // Refine the lattice on demand
    buildLatticeRefines(lat_xml, rec_lat, layout);
  else
    buildLatticeRefines(lat_xml, hex_lat, layout);
}


/// \brief Read and set lattice refines a RecLattice object.
/// \details This method can only be performed on lattices on
///          the lowest level
void GeoInputXml::buildLatticeRefines(XMLElement *lat_xml,
                                      RecLattice *lattice,
                                      LayoutVec &layout) {
  auto refines_xml = queryFirstChild(lat_xml, "refines",
                                     false,
                                     lattice->getName());
  // Return if no refines to build
  if (!refines_xml) return;

  auto rx = queryNodeInt(refines_xml, "x", false);
  auto ry = queryNodeInt(refines_xml, "y", false);
  auto rz = queryNodeInt(refines_xml, "z", false);

  // Check global refines
  auto &gr = _global_refines;
  if (gr.size() > 0) {
    rx = (gr[0] == 0) ? rx : gr[0];
    ry = (gr[1] == 0) ? ry : gr[1];
    rz = (gr[2] == 0) ? rz : gr[2];
  }

  // If refines doesn't exist, set to 1
  auto inf = xmlutils::defaultQueriedInt();
  if (rx == inf && ry == inf && rz == inf)
    return;
  if ((rx != inf && rx <= 0) ||
      (ry != inf && ry <= 0) ||
      (rz != inf && rz <= 0)) {
    log::ferror("Non-positive number of refines detected "
                "(%d, %d, %d) when parsing Lattice %s",
                rx, ry, rz, lattice->getName());
  }

  WidthVec wx(lattice->getWidthsX());
  WidthVec wy(lattice->getWidthsY());
  WidthVec wz(lattice->getWidthsZ());
  // We can safely extend each row of the layout
  if (wx.size() > 0 && rx > 1) {
    auto &vec = layout;
    auto n = vec.size();
    auto new_n = n * rx;
    vec.resize(new_n);
    for (size_t i = 0; i < n; ++i) {
      auto idx = n - 1 - i;  // Loop over existing elements
      // Fill the vector from the end
      for (auto k = new_n - rx - i*rx; k < new_n - i*rx; ++k)
        vec[k] = vec[idx];
    }
    refineWidthVector(wx, rx); // Refines the width vector
  }

  if (wy.size() > 0 && ry > 1) {
    // Save a copy of the layout
    LayoutVec origin = std::move(layout);
    layout.clear();

    // Replicate each row in the layout
    auto yz = wy.size() * wz.size();
    for (size_t r = 0; r < yz; ++r) {
      auto it = origin.begin() + r*wx.size();
      LayoutVec row(it, it + wx.size());
      layout += row * ry;
    }
    refineWidthVector(wy, ry);
  }

  if (wz.size() > 0 && rz > 1) {
    // Save a copy of the layout
    LayoutVec origin = std::move(layout);
    layout.clear();

    // Replicate each z-section in the layout
    auto xy = wx.size() * wy.size();
    for (size_t z = 0; z < wz.size(); ++z) {
      auto it = origin.begin() + z*xy;
      LayoutVec sec(it, it + xy);
      layout += sec * rz;
    }
    refineWidthVector(wz, rz);
  }

  // Reset lattice cells
  setLatticeCells(lattice, wx, wy, wz);
}


/// \brief Build HexLattice refines along z-axis
void GeoInputXml::buildLatticeRefines(XMLElement *lat_xml,
                                      HexLattice *lattice,
                                      LayoutVec &layout) {

  auto refines_xml = queryFirstChild(lat_xml, "refines",
                                     false,
                                     lattice->getName());

  // Return if no refines to build
  if (!refines_xml) return;

  auto rz = queryNodeInt(refines_xml, "z", false);

  // Check global refines
  auto &gr = _global_refines;
  if (gr.size() > 0) {
    rz = (gr[2] == 0) ? rz : gr[2];
  }

  // If refines doesn't exist, set to 1
  auto inf = xmlutils::defaultQueriedInt();
  if (rz == inf)
    return;
  if (rz != inf && rz <= 0) {
    log::ferror("Non-positive number of refines detected (%d) "
                "when parsing HexLattice %s",
                rz, lattice->getName());
  }

  WidthVec wz(lattice->getWidthsZ());

  if (wz.size() > 0 && rz > 1) {
    // Save a copy of the layout
    LayoutVec origin = std::move(layout);
    layout.clear();

    // Replicate each z-section in the layout
    auto xy = lattice->getNumLatticeCells() / wz.size();
    for (size_t z = 0; z < wz.size(); ++z) {
      auto it = origin.begin() + z*xy;
      LayoutVec sec(it, it + xy);
      layout += sec * rz;
    }
    refineWidthVector(wz, rz);
  }

  // Reset lattice cells
  setLatticeCells(lattice, lattice->getNumR(),
                  lattice->getWidthR(), wz);
}


/// \brief Read and set lattice type
latticeType GeoInputXml::readLatticeType(XMLElement *lat_xml) {

  // Element type must exist
  auto type = queryNodeString(lat_xml, "type", true);

  // Check lattice type
  latticeType lat_type = latticeType::Rectangle;
  std::string s = stringutils::toUpper(stringutils::trim(type));

  if (s == "RECTANGLE")
    ;
  else if (s == "HEXAGON")
    lat_type = latticeType::Hexagon;
  else
    log::ferror("Undefined lattice layout type: %s", type);

  return lat_type;
}


//----------------------------------------------------------------------
// Clean up objects
//----------------------------------------------------------------------
void GeoInputXml::clear() {
  // Erase unused primitives
  eraseUnusedPrimitives();

  _global_surfaces.clear();
  _global_cells.clear();
  _global_universes.clear();
  _lazy_cells_parents.clear();
  _lazy_cells_fills.clear();
  _lazy_lat_layouts.clear();
}


int GeoInputXml::eraseUnusedPrimitives() {

  // Erase the surface if it is not in use
  int s_count = 0;
  const auto used_surfaces = _geometry->getAllSurfaces();
  auto s_it = _global_surfaces.begin();
  while (s_it != _global_surfaces.end()) {
    const auto id = s_it->first;

    if (used_surfaces.count(id) == 0) {

      log::profile("Cleaning up unused surface '{}' (id={})", s_it->second->getName(), id);

      delete s_it->second;
      s_it = _global_surfaces.erase(s_it);

      // Counting erased surfaces
      ++s_count;
    }
    else
      ++s_it;
  }
  log::verbose_once("Number of surfaces cleaned = {}", s_count);

  // Erase the cell if it is not in use
  int c_count = 0;
  const auto used_cells = _geometry->getAllCells();
  auto c_it = _global_cells.begin();
  while (c_it != _global_cells.end()) {
    const auto id = c_it->first;

    if (used_cells.count(id) == 0) {

      log::profile("Cleaning up unused cell '{}' (id={})", c_it->second->getName(), id);

      delete c_it->second;
      c_it = _global_cells.erase(c_it);

      // Counting erased cells
      ++c_count;
    }
    else
      ++c_it;
  }
  log::verbose_once("Number of cells cleaned = {}", c_count);

  // Erase the universe if it is not in use
  int u_count = 0;
  const auto used_universes = _geometry->getAllUniverses();
  auto u_it = _global_universes.begin();
  while (u_it != _global_universes.end()) {
    const auto id = u_it->first;

    if (used_universes.count(id) == 0) {

      log::profile("Cleaning up unused universe '{}' (id={})", u_it->second->getName(), id);

      delete u_it->second;
      u_it = _global_universes.erase(u_it);

      // Counting erased universes
      ++u_count;
    }
    else
      ++u_it;
  }
  log::verbose_once("Number of universes cleaned = {}", u_count);

  // Log the results
  auto count = s_count + c_count + u_count;
  log::info("Number of primitives cleaned = {}", count);

  return count;
}


int GeoInputXml::eraseUnusedMaterials() {
  auto count = _matinput->eraseUnusedMaterials(_geometry);
  log::info("Number of materials cleaned = {}", count);
  return count;
}


//----------------------------------------------------------------------
// MaterialHandler interface
//----------------------------------------------------------------------

/// \brief Get a Material by id
/// \param mat_id the id of the material
/// \return a pointer to the material
Material *GeoInputXml::getMaterial(std::string mat_name) {
  if (_matinput)
    return _matinput->getMaterial(mat_name);
  else
    return nullptr;
}


} // namespace antmoc
