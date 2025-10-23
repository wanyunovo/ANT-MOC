/**
 * @file Universe.h
 * @brief The Universe class.
 * @date January 9, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#ifdef __cplusplus
#include <limits>
#include <map>
#include <set>
#include <vector>

#include "antmoc/constants.h"
#include "antmoc/enum_types.h"
#endif


namespace antmoc
{

/** Forward declarations */
class Cell;
class LocalCoords;
class Material;
class Point;


int universe_id();
void reset_universe_id();
void maximize_universe_id(int universe_id);


/**
 * @class Universe Universe.h "include/Universe.h"
 * @brief A Universe represents an unbounded space in 3D.
 * @details A Universe contains cell which are bounded subspaces in 3D
 *          which together form the Universe. Universes allow
 *          for complex, repeating (i.e. lattices) geometries to be simply
 *          represented with as few data structures as possible.
 */
class Universe {

protected:

  /** A static counter for the number of Universes */
  static int _n;

  /** A monotonically increasing unique ID for each Universe created */
  int _uid;

  /** A user-defined id for each Universe created */
  int _id;

  /** a univese id in the XML file or user-defined */
  int _input_id;

  /** A user-defined name for the Surface */
  char* _name;

  /** The type of Universe (ie, SIMPLE or LATTICE) */
  universeType _type;

  /** A collection of Cell IDs and Cell pointers in this Universe */
  std::map<int, Cell*> _cells;

  /** A boolean representing whether or not this Universe contains a Material
   *  with a non-zero fission cross-section and is fissionable */
  bool _fissionable;

  /** The extrema of the Universe */
  double _min_x;
  double _max_x;
  double _min_y;
  double _max_y;
  double _min_z;
  double _max_z;

  /** A flag for determining if boundaries are up to date */
  bool _boundaries_inspected;

  /** The boundaryTypes of the universe */
  boundaryType _min_x_bound;
  boundaryType _max_x_bound;
  boundaryType _min_y_bound;
  boundaryType _max_y_bound;
  boundaryType _min_z_bound;
  boundaryType _max_z_bound;

public:
  static constexpr int NOID = -1;
  bool anonymous() { return _input_id == NOID; }

  Universe(const int id=NOID, const char* name="");
  virtual ~Universe();
  int getUid() const;
  int getId() const;
  char* getName() const;
  universeType getType();
  int getNumCells() const;
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

  Cell* getCell(int cell_id);
  std::map<int, Cell*> getCells() const;
  std::map<int, Cell*> getAllCells();
  std::map<int, Material*> getAllMaterials();
  std::map<int, Universe*> getAllUniverses();
  bool isFissionable();

  void resetBoundaries();
  void calculateBoundaries();
  void setName(const char* name);
  void setType(universeType type);
  void addCell(Cell* cell);
  void removeCell(Cell* cell);

  bool containsPoint(Point* point);
  Cell* findCell(LocalCoords* coords);
  void setFissionability(bool fissionable);
  void subdivideCells(double max_radius = 
                      std::numeric_limits<double>::infinity());
  void buildNeighbors();
  int getInputId();

  virtual std::string toString();
  void printString();

  Universe* clone();
};

  
/**
 * @brief A helper struct for the Universe::findCell() method.
 * @details This is used to insert a Universe's Cells to the back of a vector
 *          of neighbor Cells in Universe::findCell() routine. This works in
 *          symbiosis with the pair_second method template defined below.
 */
template<typename tPair>
struct second_t {
  typename tPair::second_type operator()(const tPair& p) const {
    return p.second;
  }
};


/**
 * @brief A helper routine for the Universe::findCell() method.
 * @details This is used to insert a Universe's Cells to the back of a vector
 *          of neighbor Cells in Universe::findCell() routine. This works in
 *          symbiosis with the second_t struct template defined above.
 * @param map a std::map iterator
 * @return the second element in the iterator (e.g., map value)
 */
template<typename tMap>
second_t<typename tMap::value_type> pair_second(const tMap& map) {
  return second_t<typename tMap::value_type>();
}

} /* namespace antmoc */

#endif /* UNIVERSE_H_ */
