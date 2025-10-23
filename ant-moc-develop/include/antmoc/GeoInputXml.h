/// \file include/GeoInputXml.h
/// \brief Handling XML format.
/// \date January 2, 2019
/// \author Ya Fang, USTB (fangya201388@163.com)
///         An Wang, USTB (wangan@xs.ustb.edu.cn)
///         Haodong Shan, CIAE (shanhaodong@126.com)
///         Gen Wang, USTB (s20180686@xs.ustb.edu.cn)

#ifndef GEOINPUTXML_H_
#define GEOINPUTXML_H_

#include <set>
#include <string>
#include <vector>
#include <unordered_map>

#include "antmoc/enum_types.h"
#include "antmoc/GeoInput.h"
#include "antmoc/VecParser.h"
#include "antmoc/lattice_utils.h"
#include "antmoc/xml_utils.h"


namespace antmoc
{

// Forward declarationss
class Cell;
class Geometry;
class Material;
class Lattice;
class RecLattice;
class HexLattice;
class Surface;
class Universe;

// a vector of universe ids
using LayoutVec     = PyVector<int>;
// a map of lattices
using LatticeMap    = std::map<int, Lattice *>;
// a map of rectangle lattices
using RecLatticeMap = std::map<int, RecLattice *>;
// a map of hexagon lattices
using HexLatticeMap = std::map<int, HexLattice *>;
// a map of universes
using UniverseMap   = std::map<int, Universe *>;
// a map of cells
using CellMap       = std::map<int, Cell *>;
// a map of surfaces
using SurfaceMap    = std::unordered_map<int, Surface *>;


///---------------------------------------------------------------------
/// \class GeoInputXml
/// \brief A reader for Geometry information in XML format
///---------------------------------------------------------------------
class GeoInputXml : public GeoInput
{

private:

  // XML reader
  XMLDocument _doc_geometry;
  XMLDocument _doc_primitives;

  // Expression parser
  VecParser _parser;

  // Global primitives
  SurfaceMap _global_surfaces;
  CellMap _global_cells;
  UniverseMap _global_universes;

  // Keep track of used primitives
  std::set<int> _global_surfaces_used;
  std::set<int> _global_cells_used;
  std::set<int> _global_universes_used;

  // Lazy list
  std::map<int, std::string> _lazy_cells_parents;
  std::map<int, int>         _lazy_cells_fills;
  std::map<int, LayoutVec>   _lazy_lat_layouts;

  // Indicates whether we are building global primitives
  bool _building_primitives = false;

public:
  GeoInputXml(Geometry *geometry = nullptr, MaterialHandlerPtr matinput = nullptr)
    : GeoInput(geometry, matinput) { }

  ~GeoInputXml() = default;

  // Getters
  Material* getMaterial(std::string mat_name);
  SurfaceMap& getGlobalSurfaces()   { return _global_surfaces; }
  UniverseMap& getGlobalUniverses() { return _global_universes; }
  CellMap& getGlobalCells()         { return _global_cells; }

  //--------------------------------------------------------------------
  // Public interfaces
  //--------------------------------------------------------------------
  Geometry* readGeometryFromFile(const std::string &) override;
  Geometry* readGlobalPrimitives(const std::string &) override;
  void dumpSettings(const std::string &directory) override;

  //--------------------------------------------------------------------
  // Helpers for parsing XML
  //--------------------------------------------------------------------
  std::string replicateANumber(int number, size_t times);
  std::string extendR(std::string line);

  //-------------------------------------------------------------------
  // Read objects from XML document
  //-------------------------------------------------------------------
  void extractLatticeNodes(XMLElement *, UniverseMap &);
  void extractUniverseNodes(XMLElement *,
                            UniverseMap &,
                            std::string exerted_region = "",
                            SurfaceMap *exerted_surfaces = nullptr);
  Universe* getUniverseFromMaps(const UniverseMap &, const UniverseMap &, const int);

  void extractCellNodes(XMLElement *,
                        CellMap &,
                        std::string exerted_region = "",
                        SurfaceMap *exerted_surfaces = nullptr);
  void readXmlObjCells(XMLElement *, CellMap &);

  void extractSurfaceNodes(XMLElement *, SurfaceMap &);
  void readXmlObjSurfaces(XMLElement *, SurfaceMap &);

  void buildCellRegion(XMLElement *, Cell *, const std::string &, const SurfaceMap *);
  PyVector<int> getSurfaceIdsFromRegion(const std::string &);
  int findSurfacesByRegion(const PyVector<int> &, const SurfaceMap &, SurfaceMap &);

  //-------------------------------------------------------------------
  // Lazy-computing methods for interaction between global and local
  // objects
  //-------------------------------------------------------------------
  void processLazyCellsParents();
  void processLazyCellsFills();
  void processLazyLatticesLayouts();

  //-------------------------------------------------------------------
  // Build geometry objects and their components
  //-------------------------------------------------------------------
  void buildGeometry(XMLElement *geo_xml);
  void buildGlobalPrimitives(XMLElement *);

  /// \brief Allocate and construct a surface object.
  /// \param type The type of the surface.
  /// \param coeffs Coefficients in the surface equation.
  /// \param id The unique id of the surface.
  /// \param name The name of the surface.
  /// \return An object of Surface type.
  Surface* newSurfaceQuadratic(const char *type, const char *coeffs,
                               const int id, const char *name);

  /// \brief Allocate and construct a HexLatticePrism object.
  /// \param center X and y-coordinates of the center.
  /// \param pitch Radial pitch.
  /// \param n_rings The number of radial cells.
  /// \param orientation Lattice orientation.
  /// \param id The unique id of the surface.
  /// \param name The name of the surface.
  /// \return An object of Surface type.
  Surface* newSurfaceHexLatticePrism(const char *center, const double pitch, const int n_rings,
                                     const char *orientation, const int id, const char *name);

  /// \brief Instantiate a Surface object and set properties for it.
  /// \details Each of the properties is necessary.
  /// \param node A surface element in the XML file.
  /// \return A new Surface object.
  Surface* newSurface(XMLElement *node);

  /// \brief Instantiate a Cell object and set basic properties for it.
  /// \details Basic properties include id, name, temperature, and type.
  ///          None of the recursive properties will be set.
  /// \param node A cell element in the XML file.
  /// \return A new Cell object.
  Cell* newCell(XMLElement *node);

  /// \brief Instantiate a Universe object and set basic properties for it.
  /// \details Basic properties include id and name.
  ///          None of the recursive properties will be set.
  /// \param node A universe element in the XML file.
  /// \return A new Universe object.
  Universe* newUniverse(XMLElement *node);

  /// \brief Instantiate a Lattice object and set basic properties for it.
  /// \details Basic properties include id, name, type, offset, orientation,
  //           and widths. The layout and recursive properties are not set.
  /// \param node A lattice element in the XML file.
  /// \return A new Universe object.
  Lattice* newLattice(XMLElement *node);

  Cell* buildCell(XMLElement *,
                  std::string exerted_region = "",
                  SurfaceMap *exerted_surfaces = nullptr);

  Universe* buildUniverse(XMLElement *,
                          std::string exerted_region = "",
                          SurfaceMap *exerted_surfaces = nullptr);

  Lattice* buildLattice(XMLElement *);

  void mapUniversesByLayout(Lattice *lattice,
                            UniverseMap &unique_universes,
                            const LayoutVec &layout);

  /// \brief Build a lattice layout from 'universes' or 'row's.
  /// \param node An element which has 'universes' or 'row's.
  /// \param lattice A lattice object to save the layout.
  /// \return A layout vector.
  LayoutVec buildLatticeLayout(XMLElement *node, Lattice *lattice);

  /// \brief Create a lattice layout from a 'universes' element/attribute.
  /// \param parent An element which has a child element/attribute.
  /// \param child_name The element/attribute from which we read a layout.
  /// \return A layout vector (empty if there is no 'universes').
  LayoutVec newLayoutFromUniverses(XMLElement *parent, std::string child_name);

  /// \brief Create a lattice layout from multiple 'row' elements (RecLattice only).
  /// \param parent An element which has 'row' elements.
  /// \return A layout vector (empty if there is no 'row').
  LayoutVec newLayoutFromRows(XMLElement *parent);

  /// \brief Check whether a lattice element is defined as extruded.
  /// \param node A lattice element.
  bool isExtrudedLatticeNode(XMLElement *node);

  /// \brief Replicate the layout of an lattice if it is defined as exturded.
  /// \details The lattice widths must be set before the layout replication.
  /// \param node A lattice element.
  /// \param lattice A lattice object with widths.
  /// \param layout The lattice layout to be replicated.
  void buildExtrudedLayout(XMLElement *node, const Lattice *lattice, LayoutVec &layout);

  /// \brief Check the consistency between the layout and the number of cells.
  /// \details This method must be called after the Lattice has built its widths.
  /// \param lattice A lattice object with widths.
  /// \param layout The lattice layout to be checked.
  void checkLatticeLayout(const Lattice *lattice, const LayoutVec &layout);

  boundaryType parseBoundaryType(std::string);

  latticeType readLatticeType(XMLElement *);
  void setLatticeOffset(Lattice *, std::string);

  void checkLatticeWidths(const Lattice *, std::array<double,3>);

  void readLatticeCellWidths(XMLElement *, WidthVec &, WidthVec &, WidthVec &);
  void readLatticeCellWidths(XMLElement *, double &, WidthVec &);

  void setLatticeCells(RecLattice *, WidthVec &, WidthVec &, WidthVec &);
  void setLatticeCells(HexLattice *, int, double, WidthVec &);

  /// \brief Refine a width vector in place.
  /// \details Each element of the vector will be replaced by multiple elements
  ///          with smaller values. E.g.
  ///              Vector [1.0, 2.0, 1.0] can be refined with number 2 to be
  ///              [0.5, 0.5, 1.0, 1.0, 0.5, 0.5]
  /// \param vec A width vector to be refined in place.
  /// \param r The parameter to refine the vector.
  void refineWidthVector(WidthVec &vec, int r);

  void buildLatticeRefines(XMLElement *, Lattice *, LayoutVec &);
  void buildLatticeRefines(XMLElement *, RecLattice *, LayoutVec &);
  void buildLatticeRefines(XMLElement *, HexLattice *, LayoutVec &);


  //------------------------------------------------------------
  // Clean up
  //------------------------------------------------------------

  /// \brief Remove all of the unused global primitives.
  /// \details This method should be called ASAP. If the geometry has a complex
  ///          hierarchy, it will be a huge cost to erase unused primitives
  ///          because we have to retrieve all of the CSG objects.
  /// \return The number of erased primitives.
  int eraseUnusedPrimitives() override;

  /// \brief Remove unused materials for the geometry.
  /// \details To erase unused materials, we need a material map from the geometry.
  ///          So the return value of Geometry:getAllMaterials() must be non-empty.
  /// \return The number of erased materials.
  int eraseUnusedMaterials() override;

  /// \brief Remove unused global primitives and clear temporary containers.
  /// \details Be careful with this method because we will lose track of all of
  ///          the global primitives.
  void clear() override;
};

} // namespace antmoc

#endif // GEOINPUTXML_H_
