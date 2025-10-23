/// \file src/GeoInputXMLUtils.cpp
/// \brief Simple utilities for GeoInputXML

#include "antmoc/GeoInputXml.h"
#include "antmoc/Cell.h"
#include "antmoc/Lattice.h"
#include "antmoc/log.h"
#include "antmoc/string_utils.h"
#include "antmoc/Surface.h"

namespace antmoc {


//----------------------------------------------------------------------
// Instantiate primitives
//----------------------------------------------------------------------

Surface* GeoInputXml::newSurface(XMLElement *sf_xml) {
  // Read attributes
  auto id     = queryNodeInt(sf_xml, "id", false);
  auto name   = queryNodeString(sf_xml, "name", false);
  auto type   = queryNodeString(sf_xml, "type", true);
  auto bc     = queryNodeString(sf_xml, "boundary", false);

  if (id < 0)
    id = 0;

  if (name == defaultQueriedString())
    name = "unnamed";

  if (type == defaultQueriedString())
    log::error("Failed to build surface '{}'(id={}): no type", name, id);

  std::string type_str(type);
  stringutils::trim(type_str);

  // Debugging
  log::debug("Building surface '{}': user-defined id = {}, type = {}, boundary = {}",
              name, id, type, (bc ? bc : std::string()));

  // Create a surface
  Surface *surface = nullptr;
  if (type_str == "hex-lattice-prism") {
    auto center      = queryNodeString(sf_xml, "center", true);
    auto pitch       = queryNodeDouble(sf_xml, "pitch", true);
    auto n_rings     = queryNodeInt(sf_xml, "n_rings", false);
    auto orientation = queryNodeString(sf_xml, "orientation", false);

    // Defaults
    if (n_rings == defaultQueriedInt())
      n_rings = 1;
    if (orientation == defaultQueriedString())
      orientation = "y";

    // Debugging
    log::fdebug("Surface type is HexLatticePrism: center = %s, pitch = %f, "
                "n_rings = %d, orientation = %s",
                center, pitch, n_rings, orientation);

    surface = newSurfaceHexLatticePrism(center, pitch, n_rings,
                                        orientation, id, name);
  }
  else {
    // Construct a new quadratic surface object
    auto coeffs = queryNodeString(sf_xml, "coeffs", true);

    // Debugging
    log::fdebug("Surface type is quadratic: coeffs = '%s'", coeffs);

    surface = newSurfaceQuadratic(type, coeffs, id, name);
  }

  // Parse boundary condition
  if (bc) {
    auto bc_type = parseBoundaryType(std::string(bc));
    surface->setBoundaryType(bc_type);
  }

  return surface;
}


Surface* GeoInputXml::newSurfaceQuadratic(const char *type,
                                          const char *coeffs,
                                          const int id,
                                          const char *name) {
  Surface *surf = nullptr;
  std::string coeffs_str(coeffs);
  std::istringstream values(coeffs_str);
  std::string type_str(type);
  stringutils::trim(type_str);

  if (type_str == "plane") {
    double A, B, C, D;
    values >> A >> B >> C >> D;
    surf = new Plane(A, B, C, D, id, name);
  }
  else if (type_str == "x-plane") {
    double A;
    values >> A;
    surf = new XPlane(A, id, name);
  }
  else if (type_str == "y-plane") {
    double A;
    values >> A;
    surf = new YPlane(A, id, name);
  }
  else if (type_str == "z-plane") {
    double A;
    values >> A;
    surf = new ZPlane(A, id, name);
  }
  else if (type_str == "z-cylinder") {
    double x, y, r;
    values >> x >> y >> r;
    surf = new ZCylinder(x, y, r, id, name);
  }
  else{
    log::ferror("Undefined type '%s' found when building surface %s", type, name);
  }
  return surf;
}


Surface* GeoInputXml::newSurfaceHexLatticePrism(const char *center,
                                                const double pitch,
                                                const int n_rings,
                                                const char *orientation,
                                                const int id,
                                                const char *name) {

  std::string str(center);
  std::istringstream center_ss(str);

  double x, y;
  center_ss >> x >> y;

  return new HexLatticePrism(x, y, pitch, n_rings,
                             orientation, id, name);
}


Cell* GeoInputXml::newCell(XMLElement *cell_xml) {
  auto input_id  = queryNodeInt(cell_xml, "id", false);
  auto cell_name = queryNodeString(cell_xml, "name", false);
  // Temperature
  auto cell_temp = queryNodeDouble(cell_xml, "temperature", false);
  // Entity (pin, assembly, core, etc.)
  auto cell_phy  = queryNodeString(cell_xml, "type", false);

  if (input_id < 0) input_id = Cell::NOID;
  if (!cell_name) cell_name = "unnamed";
  if (!cell_phy) cell_phy = "ELSE";

  // Instantiate the Cell
  auto cell = new Cell(input_id, cell_name);
  cell->setPhy(cell_phy);
  cell->setTemperature(cell_temp);

  // Debugging
  log::debug("Building cell '{}': user-defined id = {}, unique id = {}",
              cell_name, input_id, cell->getId());

  return cell;
}


Universe* GeoInputXml::newUniverse(XMLElement *univ_xml) {

  auto input_id = queryNodeInt(univ_xml, "id", false);
  auto name     = queryNodeString(univ_xml, "name", false);

  if (input_id < 0) input_id = Universe::NOID;

  // Universe name defaults to unnamed
  if (!name) name = "unnamed";

  auto universe = new Universe(input_id, name);

  // Debugging
  log::debug("Building universe '{}': user-defined id = {}, unique id = {}",
              name, input_id, universe->getId());

  return universe;
}


Lattice* GeoInputXml::newLattice(XMLElement *lat_xml) {
  auto id       = queryNodeInt(lat_xml, "id", false);
  auto lat_name = queryNodeString(lat_xml, "name", false);

  if (id < 0) id = Universe::NOID;
  if (!lat_name) lat_name = "unnamed";

  // Instantiate a lattice of the proper type
  Lattice *lattice;
  auto type = readLatticeType(lat_xml);
  if (type == latticeType::Rectangle) {

    lattice = new RecLattice(id, lat_name);
    log::debug("Building RecLattice '{}': user-defined id = {}, unique id = {}",
                lat_name, id, lattice->getId());
  }
  else {

    lattice = new HexLattice(id, lat_name);
    log::debug("Building HexLattice '{}': user-defined id = {}, unique id = {}",
                lat_name, id, lattice->getId());
  }

  lattice->setLatticeType(type);

  // Retrieve the offset from a string on demand
  if ( auto offset = queryNodeString(lat_xml, "offset", false, lat_name) )
    setLatticeOffset(lattice, offset);

  auto rec_lat = dynamic_cast<RecLattice*>(lattice);
  auto hex_lat = dynamic_cast<HexLattice*>(lattice);

  // Read and set widths, lattice cells
  auto widths_xml = queryFirstChild(lat_xml, "widths", true, lat_name);
  if (rec_lat) {
    WidthVec widths_x, widths_y, widths_z;  ///< widths of nonuniform lattice
    readLatticeCellWidths(widths_xml, widths_x, widths_y, widths_z);
    // Set the widths of each Lattice cell
    setLatticeCells(rec_lat, widths_x, widths_y, widths_z);
  }
  else {
    // Read widths of a HexLattice
    double wr;          ///< pitch
    WidthVec widths_z;  ///< widths of nonuniform lattice
    readLatticeCellWidths(widths_xml, wr, widths_z);

    // Read number of radial tiles of a HexLattice
    int nr = queryNodeInt(lat_xml, "n_rings", true, lat_name);
    setLatticeCells(hex_lat, nr, wr, widths_z);

    // Read orientation of a HexLattice
    auto orientation = queryNodeString(lat_xml, "orientation", false, lat_name);
    if (orientation)
      hex_lat->setOrientation(orientation);
  }

  // FIXME, Check consistency between widths from different sources
  //checkLatticeWidths(lattice, box);

  return lattice;
}


//----------------------------------------------------------------------
// Build lattice components
//----------------------------------------------------------------------

LayoutVec GeoInputXml::newLayoutFromUniverses(XMLElement *parent,
                                              std::string child_name) {

  // Read the XML node
  auto univs = queryNodeString(parent, child_name, false);

  // The layout is stored just as it looks like
  LayoutVec layout;

  if (univs) {
    // Extend non-digits of the string
    //auto u_str = extendR(univs);
    auto words = stringutils::splitString(univs);
    for (auto &w : words) {
      try {
        int id = std::stoi(w);
        layout.push_back(id);
      }
      catch (std::invalid_argument &e) {
        _parser.parse(w);
        auto ids = _parser.getOperand().toIntVector();
        layout += ids;
      }
    }
  }

  return layout;
}


LayoutVec GeoInputXml::newLayoutFromRows(XMLElement *parent) {

  // Read the XML node
  auto row_xml = parent->FirstChildElement("row");

  // The layout is stored just as it looks like
  LayoutVec layout;

  while (row_xml) {
    // Extend non-digits of the string
    auto row_str = extendR(row_xml->GetText());

    // Convert the string to a vector
    std::istringstream row_ss(row_str);
    LayoutVec row;
    int item;

    while (row_ss >> item)
      row.push_back(item);

    // Parse 'repeat' attribute to determine the number of copies
    auto repeat = queryNodeInt(row_xml, "repeat", false);

    // if there is no attribute or failed to query
    if (repeat == std::numeric_limits<int>::min())
      repeat = 1;
    else if (repeat < 0)
      log::error("Negetive value of 'repeat' found when building lattice layout");

    // Insert the row as separated items
    for (int r = 0; r < repeat; ++r)
      layout += row;

    row_xml = row_xml->NextSiblingElement("row");
  }

  return layout;
}


bool GeoInputXml::isExtrudedLatticeNode(XMLElement *node) {

  auto axial = queryNodeString(node, "axial", false);
  if (axial)
    return stringutils::toUpper(stringutils::trim(std::string(axial))) == "EXTRUDED";
  else
    return false;
}


void GeoInputXml::buildExtrudedLayout(XMLElement *node,
                                      const Lattice *lattice,
                                      LayoutVec &layout) {
  if (isExtrudedLatticeNode(node)) {
    size_t nz  = lattice->getNumZ();
    size_t nxy = lattice->getNumLatticeCells() / nz;

    if (nxy != layout.size()) {
      log::error("Lattice {} is defined as extruded, but the number of lattice "
                 "cells({}) doesn't equal the layout size({})",
                  lattice->getName(), nxy, layout.size());
    }

    const LayoutVec layout_copy(layout);
    for (size_t i = 1; i < nz; ++i)
        layout += layout_copy;
  }
}


void GeoInputXml::checkLatticeLayout(const Lattice *lattice,
                                     const LayoutVec &layout) {
  size_t n = lattice->getNumLatticeCells();

  if (layout.size() != n) {
    log::error("When building lattice '{}'(id={}), the number of elements({}) in the "
               "lattice layout does not equal the number of lattice cells({})",
                lattice->getName(), lattice->getId(), layout.size(), n);
  }
}


/// \brief Refine a width vector in place.
/// \details Each element of the vector will be replaced by multiple elements
///          with smaller values. E.g.
///              Vector [1.0, 2.0, 1.0] can be refined with number 2 to be
///              [0.5, 0.5, 1.0, 1.0, 0.5, 0.5]
/// \param vec A width vector to be refined in place.
/// \param r The parameter to refine the vector.
void GeoInputXml::refineWidthVector(WidthVec &vec, int r) {

  // Compute the size of the new vector
  auto n = vec.size();
  auto new_n = n * r;
  vec.resize(new_n);

  // Loop over existing elements
  for (size_t i = 0; i < n; ++i) {
    auto idx = n - 1 - i;       // Get the reverse iterator
    auto value = vec[idx] / r;  // Compute the new value for each element
    // Fill the vector from the end
    for (auto k = new_n - r - i*r; k < new_n - i*r; ++k)
      vec[k] = value;
  }
}

} // namespace antmoc
