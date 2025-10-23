/// \file include/Cell.h
/// \brief The Cell class.
/// \date January 18, 2012
/// \author William Boyd, MIT, Course 22 (wboyd@mit.edu)


#ifndef CELL_H_
#define CELL_H_

#ifdef __cplusplus
#include <limits>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <vector>

#include "antmoc/BoolParser.h"
#include "antmoc/enum_types.h"
#include "antmoc/math_utils.h"
#endif

namespace antmoc
{

// Forward declarationss to resolve circular dependencies
class Material;
class RecLattice;
class LocalCoords;
class Point;
class Surface;
class Universe;
class CellIter;

int cell_id();
void reset_cell_id();
void maximize_cell_id(int cell_id);

using SurfaceMap = std::unordered_map<int, Surface*>;

struct Halfspace {
  int halfspace;
  Surface *surface;
};

/// \class Cell Cell.h "src/Cell.h"
/// \brief Represents a Cell inside of a Universe.
class Cell {

private:

  ///< A static counter for the number of Cells
  static int _n;

  ///< A bool parser
  static BoolParser parser;

  ///< A monotonically increasing unique ID for each Cell created
  int _uid;

  ///< A user-defined ID for each Cell created
  int _id;

  ///< A user-defined ID
  int _input_id;

  ///< A user-defined name for the Surface
  char* _name;

  ///< The type of Cell (ie MATERIAL or FILL)
  cellType _cell_type;

  ///< A pointer to the Material or Universe filling this Cell
  void* _fill;

  ///< The volume / area of the Cell computed from overlapping segments
  double _volume;

  ///< The total number of instances of this Cell in the Geometry
  int _num_instances;

  ///< A boolean indicating whether to cell is rotated
  bool _rotated;

  ///< An array with angles in degrees for rotations about x, y, and z
  double _rotation[3];

  ///< A rotation matrix defined in terms of the rotation angles
  double _rotation_matrix[9];

  ///< A boolean indicating whether to cell is translated
  bool _translated;

  ///< An array with translations in x, y and z
  double _translation[3];

  ///< The number of rings sub-dividing this Cell
  int _num_rings;

  ///< The number of sectors sub-dividing this Cell
  int _num_sectors;

  ///< A parent Cell if cloned by another Cell
  Cell* _parent;

  ///< A boolean indicating this is an intersection of a set of halfspaces
  bool _simple;

  ///< An abstract syntax tree of the region
  NodePtr _ast;

  ///< A reverse polish notation of the region which can be generated from _ast
  PyVector<int> _rpn;

  ///< Map of bounding Surface IDs with pointers
  SurfaceMap _surfaces;

  ///< Vector of neighboring Cells
  std::vector<Cell*> _neighbors;

  //材料Cell拥有温度
  float _temperature;
  //Cell表示的物理含义（棒、组件、堆芯等）
  cellPhy _cell_phy;

  void ringify(std::vector<Cell*>& subcells, double max_radius);
  void sectorize(std::vector<Cell*>& subcells);

public:
  static constexpr int NOID = -1;
  bool anonymous() { return _input_id == NOID; }

  //-----------------------------------------------------
  // Iterators
  //-----------------------------------------------------
  friend class CellIter;
  CellIter begin();
  CellIter end();

  Cell(int id=NOID, const char* name="");
  virtual ~Cell();
  int getUid() const;
  int getId() const;
  int getInputId();
  char* getName() const;
  cellType getType() const;
  Material* getFillMaterial();
  Universe* getFillUniverse();
  double getVolume();
  int getNumInstances();
  bool isRotated();
  bool isTranslated();
  double getPhi(std::string units="degrees");
  double getTheta(std::string units="degrees");
  double getPsi(std::string units="degrees");
  double* getRotationMatrix();
  double* getTranslation();
  void retrieveRotation(double* rotations, int num_axes,
			std::string units="degrees");
  void retrieveTranslation(double* translations, int num_axes);
  int getNumRings();
  int getNumSectors();
  double getMinX();
  double getMaxX();
  double getMinY();
  double getMaxY();
  double getMinZ();
  double getMaxZ();
  boundaryType getMinXBoundaryType();
  boundaryType getMaxXBoundaryType();
  boundaryType getMinYBoundaryType();
  boundaryType getMaxYBoundaryType();
  boundaryType getMinZBoundaryType();
  boundaryType getMaxZBoundaryType();
  int getNumSurfaces() const;
  SurfaceMap getSurfaces() const;
  std::vector<Cell*> getNeighbors() const;
  bool hasParent();
  Cell* getParent();
  Cell* getOldestAncestor();

  int getNumZCylinders();

  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setName(const char* name);
  void setFill(Material* fill);
  void setFill(Universe* fill);
  void setVolume(double volume);
  void incrementVolume(double volume);
  void setNumInstances(int num_instances);
  void incrementNumInstances();
  void setRotation(double* rotation, int num_axes, std::string units="degrees");
  void setTranslation(double* translation, int num_axes);
  void setNumRings(int num_rings);
  void setNumSectors(int num_sectors);
  void setParent(Cell* parent);
  void addSurface(int halfspace, Surface* surface);
  void removeSurface(Surface* surface);
  void addNeighborCell(Cell* cell);

  std::string getRegion() const { return parser.toString(_rpn); }
  PyVector<int> getRPN() const  { return _rpn; }
  void setSurfaces(const SurfaceMap surfaces) { _surfaces = surfaces; }
  void addRegionSurfaces(const std::string region, const SurfaceMap &surfaces);
  void generateRPN();

  /// \brief Indicate whether a token represents a halfspace
  bool isValidHalfspace(int token)
    { return !parser.isOperator(token) && _surfaces.count(abs(token)); }

  bool isFissionable();
  bool containsPoint(Point* point);
  bool containsSimple(Point* point);
  bool containsComplex(Point* point);
  bool containsCoords(LocalCoords* coords);
  double minSurfaceDist(Point* point, double azim, double polar);
  double minSurfaceDist(LocalCoords* coords);

  Cell* clone(bool clone_region=true);
  void subdivideCell(double max_radius);
  void buildNeighbors();

  // Extra parameters and complement
  void setTemperature(float);
  float getTemperature();
  void setPhy(const char*);
  cellPhy getPhy();

  std::string toString();
  void printString();
};


/// \class An iterator over cell halfspaces
/// \details This iterator only iterates over valid tokens. Neither operators
///          nor deleted surfaces will be iterated.
///          Note that the iterator will be invalid if the RPN of the cell
///          is changed inside the iteration.
class CellIter {
public:
  CellIter(Cell &cell, size_t pos = 0)
    : _cell(cell), _pos(pos) {
    // Start from the first valid halfspace
    while (_pos < _cell._rpn.size()) {
      int token = _cell._rpn[_pos];
      if (_cell.isValidHalfspace(token))
        break;
      ++_pos;
    }
  }

  /// \brief Two CellIters equal iff they are at the same valid position
  bool operator==(const CellIter &rhs) { return (_pos == rhs._pos); }

  bool operator!=(const CellIter &rhs) { return !(*this == rhs); }

  /// \brief Dereference the iterator
  /// \return a halfspace
  Halfspace operator*() {
      int token = _cell._rpn[_pos];
      auto surface = _cell._surfaces[abs(token)];
      return {copysign(1, token), surface};
  }

  /// \brief Move to the next valid token
  CellIter& operator++()
  {
    // Find the next valid index
    while (_pos < _cell._rpn.size() - 1) {
      ++_pos;
      int token = _cell._rpn[_pos];
      if (_cell.isValidHalfspace(token))
        return *this;
    }
    // The index past the last element
    _pos = _cell._rpn.size();
    return *this;
  }

  /// \brief Return the current position of the iterator
  size_t getPos() const { return _pos; }

protected:
  Cell &_cell;  ///< The cell to be iterated
  size_t _pos;  ///< Index of the current halfspace
};


} // namespace antmoc

#endif  // CELL_H_
