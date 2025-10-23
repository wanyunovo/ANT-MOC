#include "antmoc/Universe.h"
#include "antmoc/Cell.h"
#include "antmoc/Lattice.h"
#include "antmoc/LocalCoords.h"
#include "antmoc/log.h"
#include "antmoc/Material.h"
#include "antmoc/Surface.h"


namespace antmoc
{

int Universe::_n = 0;

static int auto_id = DEFAULT_INIT_ID;
static std::set<int> used_ids;

/**
 * @brief Returns an auto-generated unique Universe ID.
 * @details This method is intended as a utility method for user's writing
 *          geometry input files. The method makes use of a static Universe
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined Universe IDs greater
 *          than or equal to 10000 is prohibited.
 */
int universe_id() {
  int id = auto_id;
  auto_id++;
  while (used_ids.find(id) != used_ids.end()) {
    id = auto_id;
    auto_id++;
  }
  return id;
}


/**
 * @brief Resets the auto-generated unique Universe ID counter to 10000.
 */
void reset_universe_id() {
  auto_id = DEFAULT_INIT_ID;
  used_ids.clear();   // Clear used ids
}


/**
 * @brief Maximize the auto-generated unique Universe ID counter.
 * @details This method updates the auto-generated unique Universe ID
 *          counter if the input parameter is greater than the present
 *          value. This is useful for the OpenMC compatibility module
 *          to ensure that the auto-generated Universe IDs do not
 *          collide with those created in OpenMC.
 * @param universe_id the id assigned to the auto-generated counter
 */
void maximize_universe_id(int universe_id) {
  if (universe_id > auto_id)
    auto_id = universe_id;
}


/**
 * @brief Constructor assigns a unique and user-specified ID for the Universe.
 * @param id the user-specified optional Universe ID
 * @param name the user-specified optional Universe ID
 */
Universe::Universe(const int id, const char* name) {

  /* If the user did not define an optional ID, create one */
//  if (id == -1)
//    _id = universe_id();

  /* Use the user-defined ID */
//  else
//    _id = id;
 /* create a universe unique id  */

  // If the user-defined ID was used,
  // create a universe unique id
  if (id < 0 || used_ids.count(id) > 0)
    // create a universe unique id
    _id = universe_id();
  else
    _id = id;

  _input_id = id;

  _uid = _n;
  _n++;

  /* Add the ID to the used set */
  used_ids.insert(_id);

  _name = NULL;
  setName(name);

  _type = SIMPLE;

  _boundaries_inspected = false;

  /* By default, the Universe's fissionability is unknown */
  _fissionable = false;
}


/**
 * @brief Destructor clears the Cell pointers container.
 */
Universe::~Universe() {

  delete [] _name;

  /* Clear the map of Cells, cells are deallocated either with the Geometry
     or automatically at the end of your input file */
  _cells.clear();
}


/**
 * @brief Returns the Universe's unique ID.
 * @return the Universe's unique ID.
 */
int Universe::getUid() const {
  return _uid;
}


/**
 * @brief Return the user-specified ID for this Universe.
 * @return the user-specified Universe ID
 */
int Universe::getId() const {
  return _id;
}


/**
 * @brief Return the user-defined name of the Universe.
 * @return the Universe name
 */
char* Universe::getName() const {
  return _name;
}


/**
 * @brief Return the Universe type (SIMPLE or LATTICE).
 * @return the Universe type
 */
universeType Universe::getType() {
  return _type;
}


/**
 * @brief Return the number of Cells in this Universe.
 * @return the number of Cells
 */
int Universe::getNumCells() const {
  return _cells.size();
}


/**
 * @brief Returns the minimum reachable x-coordinate in the Universe.
 * @return the minimum reachable x-coordinate
 */
double Universe::getMinX() {

  if (!_boundaries_inspected)
    calculateBoundaries();

  return _min_x;
}


/**
 * @brief Returns the maximum reachable x-coordinate in the Universe.
 * @return the maximum reachable x-coordinate
 */
double Universe::getMaxX() {

  if (!_boundaries_inspected)
    calculateBoundaries();

  return _max_x;
}


/**
 * @brief Returns the minimum reachable y-coordinate in the Universe.
 * @return the minimum reachable y-coordinate
 */
double Universe::getMinY() {

  if (!_boundaries_inspected)
    calculateBoundaries();

  return _min_y;
}


/**
 * @brief Returns the maximum reachable y-coordinate in the Universe.
 * @return the maximum reachable y-coordinate
 */
double Universe::getMaxY() {

  if (!_boundaries_inspected)
    calculateBoundaries();

  return _max_y;
}


/**
 * @brief Returns the minimum reachable z-coordinate in the Universe.
 * @return the minimum reachable z-coordinate
 */
double Universe::getMinZ() {

  if (!_boundaries_inspected)
    calculateBoundaries();

  return _min_z;
}


/**
 * @brief Returns the maximum reachable z-coordinate in the Universe.
 * @return the maximum reachable z-coordinate
 */
double Universe::getMaxZ() {

  if (!_boundaries_inspected)
    calculateBoundaries();

  return _max_z;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum
 *        reachable x-coordinate in the Universe.
 * @return the boundary conditions at the minimum reachable x-coordinate
 */
boundaryType Universe::getMinXBoundaryType() {

  if (!_boundaries_inspected)  //universe 边界是否是最新的
    calculateBoundaries();

#ifdef ONLYVACUUMBC
  if (_min_x_bound != VACUUM && _min_x_bound != INTERFACE)
    log::error("The code was compiled specially for cases with only "
               "vacuum boundary conditions and a reflective or periodic "
               "boundary condition was found in universe {}.", _id);
#endif

  return _min_x_bound;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum
 *        reachable x-coordinate in the Universe.
 * @return the boundary conditions at the maximum reachable x-coordinate
 */
boundaryType Universe::getMaxXBoundaryType() {

  if (!_boundaries_inspected)
    calculateBoundaries();

#ifdef ONLYVACUUMBC
  if (_max_x_bound != VACUUM && _max_x_bound != INTERFACE)
    log::error("The code was compiled specially for cases with only "
               "vacuum boundary conditions and a reflective or periodic "
               "boundary condition was found in universe {}.", _id);
#endif

  return _max_x_bound;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum
 *        reachable y-coordinate in the Universe.
 * @return the boundary conditions at the minimum reachable y-coordinate
 */
boundaryType Universe::getMinYBoundaryType() {

  if (!_boundaries_inspected)
    calculateBoundaries();

#ifdef ONLYVACUUMBC
  if (_min_y_bound != VACUUM && _min_y_bound != INTERFACE)
    log::error("The code was compiled specially for cases with only "
               "vacuum boundary conditions and a reflective or periodic "
               "boundary condition was found in universe {}.", _id);
#endif

  return _min_y_bound;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum
 *        reachable y-coordinate in the Universe.
 * @return the boundary conditions at the maximum reachable y-coordinate
 */
boundaryType Universe::getMaxYBoundaryType() {

  if (!_boundaries_inspected)
    calculateBoundaries();

#ifdef ONLYVACUUMBC
  if (_max_y_bound != VACUUM && _max_y_bound != INTERFACE)
    log::error("The code was compiled specially for cases with only "
               "vacuum boundary conditions and a reflective or periodic "
               "boundary condition was found in universe {}.", _id);
#endif

  return _max_y_bound;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum
 *        reachable z-coordinate in the Universe.
 * @return the boundary conditions at the minimum reachable z-coordinate
 */
boundaryType Universe::getMinZBoundaryType() {

  if (!_boundaries_inspected)
    calculateBoundaries();

#ifdef ONLYVACUUMBC
  if (_min_z_bound != VACUUM && _min_z_bound != INTERFACE)
    log::warn("The code was compiled specially for cases with only "
              "vacuum boundary conditions and a reflective or periodic "
              "boundary condition was found in universe {}.", _id);
#endif

  return _min_z_bound;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum
 *        reachable z-coordinate in the Universe.
 * @return the boundary conditions at the maximum reachable z-coordinate
 */
boundaryType Universe::getMaxZBoundaryType() {

  if (!_boundaries_inspected)
    calculateBoundaries();

#ifdef ONLYVACUUMBC
  if (_max_z_bound != VACUUM && _max_z_bound != INTERFACE)
    log::warn("The code was compiled specially for cases with only "
              "vacuum boundary conditions and a reflective or periodic "
              "boundary condition was found in universe {}.", _id);
#endif

  return _max_z_bound;
}


/**
 * @brief Returns a Cell in this universe.
 * @param cell_id the integer the cell_id
 * @return Returns the cell pointer.
 */
Cell* Universe::getCell(int cell_id) {

  if (_cells.find(cell_id) == _cells.end()) {
    log::ferror("Unable to return Cell with ID = %d from Universe with "
               "ID = %d since it does not contain this Cell", cell_id, _id);
  }

    return _cells.at(cell_id);
}


/**
 * @brief Return the container of Cell IDs and Cell pointers in this Universe.
 * @return std::map of Cell IDs
 */
std::map<int, Cell*> Universe::getCells() const {
  return _cells;
}


/**
 * @brief Returns the std::map of Cell IDs and Cell pointers in this Universe
 *        at all nested Universe levels.
 * @return std::map of Cell IDs and pointers
 */
std::map<int, Cell*> Universe::getAllCells() {

  std::map<int, Cell*> cells;
  std::map<int, Cell*>::iterator iter;

  /* Add this Universe's Cells to the map */
  cells.insert(_cells.begin(), _cells.end());

  /* Append all Cells in each Cell in the Universe to the map */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    std::map<int, Cell*> nested_cells = iter->second->getAllCells();
    cells.insert(nested_cells.begin(), nested_cells.end());
  }

  return cells;
}


/**
 * @brief Returns the std::map of all IDs and Material pointers filling
 *        this Universe.
 * @return std::map of Material IDs and pointers
 */
std::map<int, Material*> Universe::getAllMaterials() {

  std::map<int, Cell*> cells = getAllCells();
  std::map<int, Material*> materials;

  for (const auto &cell : cells) {
    if (cell.second->getType() == MATERIAL) {
      auto material = cell.second->getFillMaterial();
      materials[material->getId()] = material;
    }
  }

  return materials;
}


/**
 * @brief Returns the std::map of all nested Universe IDs and Universe pointers
 *         filling this Universe.
 * @return std::map of Universe IDs and pointers
 */
std::map<int, Universe*> Universe::getAllUniverses() {

  /* Get all Cells in this Universe */
  std::map<int, Cell*> cells = getAllCells();

  std::map<int, Universe*> universes;
  universes[_id] = this;
  std::map<int, Cell*>::iterator iter;
  Cell* cell;

  /* Append all Universes containing each Cell to the map */
  for (iter = cells.begin(); iter != cells.end(); ++iter) {
    cell = iter->second;
    std::map<int, Universe*> nested_universes = cell->getAllUniverses();
    universes.insert(nested_universes.begin(), nested_universes.end());
  }

  return universes;
}


/**
 * @brief Returns true if the Universe contains a Cell filled by a fissionable
 *        Material and false otherwise.
 * @details This method should not be called prior to the calling of the
 *          Geometry::computeFissionability() method.
 * @return true if contains a fissionable Material
 */
bool Universe::isFissionable() {
  return _fissionable;
}


/**
 * @brief Sets the name of the Universe.
 * @param name the Universe name string
 */
void Universe::setName(const char* name) {
  int length = strlen(name);

  if (_name != NULL)
    delete [] _name;

  /* Initialize a character array for the Universe's name */
  _name = new char[length+1];

  /* Copy the input character array Universe name to the class attribute name */
  for (int i=0; i <= length; i++)
    _name[i] = name[i];
}


/**
 * @brief Sets the Universe type to SIMPLE or LATTICE.
 * @param type the Universe type
 */
void Universe::setType(universeType type) {
  _type = type;
}


/**
 * @brief Sets whether or not this Universe contains a fissionable Material
 *        with a non-zero fission cross-section.
 * @details This method is called by the Geometry::computeFissionability()
 *          class method.
 * @param fissionable true if the Universe contains a fissionable Material;
 *        false otherwise
 */
void Universe::setFissionability(bool fissionable) {
  _fissionable = fissionable;
}


/**
 * @brief Adds a Cell to this Universe.
 * @details Stores the user-specified Cell ID and Cell pointer in a std::map
 *          along with all of other Cells added to this Universe.
 * @param cell the Cell pointer
 */
void Universe::addCell(Cell* cell) {

  try {
    _cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
    log::fdebug("Added Cell with ID = %d to Universe with ID = %d",
               cell->getId(), _id);
  }
  catch (std::exception &e) {
    log::ferror("Unable to add Cell with ID = %d to Universe with"
               " ID = %d. Backtrace:\n%s", _id, e.what());
  }

  _boundaries_inspected = false;
}


/**
 * @brief Removes a Cell from this Universe's container of Cells.
 * @param cell a pointer to the Cell to remove
 */
void Universe::removeCell(Cell* cell) {
  if (_cells.find(cell->getId()) != _cells.end())
    _cells.erase(cell->getId());

  _boundaries_inspected = false;
}


/**
 * @brief Finds the Cell for which a LocalCoords object resides.
 * @details Finds the Cell that a LocalCoords object is located inside by
 *          checking each of this Universe's Cells. Returns NULL if the
 *          LocalCoords is not in any of the Cells.
 * @param coords a pointer to the LocalCoords of interest
 * @return a pointer the Cell where the LocalCoords is located
 */
Cell* Universe::findCell(LocalCoords* coords) {

  /* Sets the LocalCoord type to UNIV at this level */
  coords->setType(UNIV);

  /* Loop over all Cells */
  std::map<int,Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    Cell* cell = iter->second;

    if (cell->containsCoords(coords)) {

      /* Set the Cell on this level */
      coords->setCell(cell);

      /* MATERIAL type Cell - lowest level, terminate search for Cell */
      if (cell->getType() == MATERIAL)
        return cell;

      /* FILL type Cell - Cell contains a Universe at a lower level
       * Update coords to next level and continue search */
      else if (cell->getType() == FILL) {

        LocalCoords* next_coords = coords->getNextCreate(coords->getX(),
                                                         coords->getY(),
                                                         coords->getZ());

        /* Apply translation to position in the next coords */
        if (cell->isTranslated()) {
          double* translation = cell->getTranslation();
          double new_x = coords->getX() - translation[0];
          double new_y = coords->getY() - translation[1];
          double new_z = coords->getZ() - translation[2];
          next_coords->setX(new_x);
          next_coords->setY(new_y);
          next_coords->setZ(new_z);
        }

        /* Apply rotation to position in the next coords */
        if (cell->isRotated()) {
          double x = next_coords->getX();
          double y = next_coords->getY();
          double z = next_coords->getZ();
          double* matrix = cell->getRotationMatrix();
          double new_x = matrix[0] * x + matrix[1] * y + matrix[2] * z;
          double new_y = matrix[3] * x + matrix[4] * y + matrix[5] * z;
          double new_z = matrix[6] * x + matrix[7] * y + matrix[8] * z;
          next_coords->setX(new_x);
          next_coords->setY(new_y);
          next_coords->setZ(new_z);
        }

        Universe* univ = cell->getFillUniverse();
        next_coords->setUniverse(univ);
        coords->setCell(cell);

        if (univ->getType() == SIMPLE)
          return univ->findCell(next_coords);
        else
          return static_cast<Lattice*>(univ)->findCell(next_coords);
      }
    }
  }

  return NULL;
}


/**
 * @brief Subdivides all of the Material-filled Cells within this Universe
 *        into rings and angular sectors aligned with the z-axis.
 * @param max_radius the maximum allowable radius used in the subdivisions
 */
void Universe::subdivideCells(double max_radius) {

  log::fdebug("Subdividing Cells for Universe ID=%d "
             "with max radius %f", _id, max_radius);

  std::map<int, Cell*>::iterator iter;

  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    /* Cells filled with Materials */
    if (iter->second->getType() == MATERIAL) {
      Cell* cell = iter->second;

      if (cell->getNumRings() > 0 || cell->getNumSectors() > 0)
        cell->subdivideCell(max_radius);
    }

    /* Cells filled with Universes */
    else {
      Universe* fill = iter->second->getFillUniverse();
      if (fill->getType() == SIMPLE)
        fill->subdivideCells(max_radius);
      else
        static_cast<Lattice*>(fill)->subdivideCells(max_radius);
    }
  }
}


/**
 * @brief Builds collections of neighboring Cells for all Cells in this
 *        Universe for optimized ray tracing.
 */
void Universe::buildNeighbors() {

  /* Loop over all of the Universe's Cells and make recursive call */
  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter)
    iter->second->buildNeighbors();
}


/**
 * @brief Convert the member attributes of this Universe to a character array.
 * @return a character array representing the Universe's attributes
 */
std::string Universe::toString() {

  std::stringstream string;
  std::map<int, Cell*>::iterator iter;

  string << "Universe ID = " << _id;
  string << ", input Id = " << _input_id;
  string << ", name = " << _name;

  string << ", type = ";
  if (_type == SIMPLE)
    string << "SIMPLE";
  else
    string << "LATTICE";

  string << ", # cells = " << _cells.size() << ", cell IDs = ";

  for (iter = _cells.begin(); iter != _cells.end(); ++iter)
    string << iter->first << ", ";

  return string.str();
}


/**
 * @brief Prints a string representation of the Universe's attributes to
 *        the console.
 */
void Universe::printString() {
  log::fresult(toString().c_str());
}


/**
 * @brief Clones this Universe and copy cells map.
 * @return a pointer to the Universe clone
 */
Universe* Universe::clone() {

  log::fdebug("Cloning Universe %d", _id);

  /* Instantiate new Universe clone */
  Universe* clone = new Universe(universe_id(), _name);

  /* Loop over Cells in Universe and clone each one */
  std::map<int, Cell*>::iterator iter1;
  for (iter1 = _cells.begin(); iter1 != _cells.end(); ++iter1) {

    /* If the Cell is filled with a Material, clone it */
    if (iter1->second->getType() == MATERIAL) {

      /* Clone the Cell */
      Cell* parent = static_cast<Cell*>(iter1->second);
      Cell* cell_clone = parent->clone();

      /* Add Cell clone to the list */
      clone->addCell(cell_clone);
    }
    /* Throw error message if Cell is FILL type */
    else {
      log::ferror("Unable to clone Universe %d since it contains Cell %d"
                 "which is filled with a Universe rather than a Material");
    }
  }

  return clone;
}


/**
  * @brief  Calculates the boundary locations and conditions (VACUUM or
  *         REFLECTIVE) at the maximum and minimum reachable coordinates in the
  *         Universe.
  * 这个函数的作用是计算universe中最大和最小可达坐标处的边界条件以及边界位置，并把变量_boundaries_inspected设为true
  */
void Universe::calculateBoundaries() {

  /* Calculate the minimum reachable x-coordinate in the geometry and store it
   * in _min_x */
  double min_x = std::numeric_limits<double>::infinity();  //设为无限大
  Surface* surf;
  int halfspace;

  /* Calculate the boundary condition at the minimum reachable x-coordinate in
   * the Universe and store it in _min_x_bound */
  _min_x_bound = BOUNDARY_NONE;   //无类型边界

  /* Check if the universe contains a cell with an x-min boundary */
  for (const auto &c : _cells) {   //遍历root-universe下的所有cell

    double cell_min_x = -std::numeric_limits<double>::infinity();  //无限小
    boundaryType cell_min_x_bound = BOUNDARY_NONE;
    for (const auto &hs : *(c.second)) {  //遍历一个cell下的所有表面和半空间信息
      surf = hs.surface;
      halfspace = hs.halfspace;

      if (surf->getSurfaceType() == XPLANE && halfspace == +1 &&
          surf->getBoundaryType() != BOUNDARY_NONE) {  //如果找到x平面且半空间为正方向的表面
        if (surf->getMinX(halfspace) > cell_min_x) { //找到+1方向上围成这个cell中最里面（求交）的xplane
          cell_min_x = surf->getMinX(halfspace);
          cell_min_x_bound = surf->getBoundaryType();
        }
      }
    }
    if (cell_min_x_bound != BOUNDARY_NONE && cell_min_x < min_x) {  //如果找到这个cell中符合条件的surface
      min_x = cell_min_x;
      _min_x_bound = cell_min_x_bound;
    }
  }

  /* If a x-min boundary was not found, get the x-min from the bounding boxes
   * of the cells */
  if (min_x == std::numeric_limits<double>::infinity()) { 
    //如果 min_x 仍然是无穷大，表示未找到有效的 x-min 边界条件，没有找到有效的半平面，后面就直接从cell的边界处寻找
    double cell_min_x = min_x;
    boundaryType cell_min_x_bound = BOUNDARY_NONE;
    for (const auto &c : _cells) {  //接着检查cell的边界值，找到最小 x 坐标 cell_min_x 和相应的边界条件 cell_min_x_bound。
      cell_min_x = c.second->getMinX();
      cell_min_x_bound = c.second->getMinXBoundaryType();
      if (cell_min_x_bound != BOUNDARY_NONE && cell_min_x < min_x) {
        min_x = cell_min_x;
        _min_x_bound = cell_min_x_bound;
      }
    }
  }

  _min_x = min_x;

  /* Calculate the maximum reachable x-coordinate in the geometry and store it
   * in _max_x */
  double max_x = -std::numeric_limits<double>::infinity();

  /* Calculate the boundary condition at the maximum reachable x-coordinate in
   * the Universe and store it in _max_x_bound */
  _max_x_bound = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an x-max boundary */
  for (const auto &c : _cells) {

    double cell_max_x = std::numeric_limits<double>::infinity();
    boundaryType cell_max_x_bound = BOUNDARY_NONE;
    for (const auto &hs : *(c.second)) {
      surf = hs.surface;
      halfspace = hs.halfspace;

      if (surf->getSurfaceType() == XPLANE && halfspace == -1 &&
          surf->getBoundaryType() != BOUNDARY_NONE) {
        if (surf->getMaxX(halfspace) < cell_max_x) {
          cell_max_x = surf->getMaxX(halfspace);
          cell_max_x_bound = surf->getBoundaryType();
        }
      }
    }
    if (cell_max_x_bound != BOUNDARY_NONE && cell_max_x > max_x) {
      max_x = cell_max_x;
      _max_x_bound = cell_max_x_bound;
    }
  }

  /* If a x-max boundary was not found, get the x-max from the bounding boxes
   * of the cells */
  if (max_x == -std::numeric_limits<double>::infinity()) {
    double cell_max_x = max_x;
    boundaryType cell_max_x_bound = BOUNDARY_NONE;
    for (const auto &c : _cells) {
      cell_max_x = c.second->getMaxX();
      cell_max_x_bound = c.second->getMaxXBoundaryType();
      if (cell_max_x_bound != BOUNDARY_NONE && cell_max_x > max_x) {
        max_x = cell_max_x;
        _max_x_bound = cell_max_x_bound;
      }
    }
  }

  _max_x = max_x;

  /* Calculate the minimum reachable y-coordinate in the geometry and store it
   * in _min_y */
  double min_y = std::numeric_limits<double>::infinity();

  /* Calculate the boundary condition at the minimum reachable y-coordinate in
   * the Universe and store it in _min_y_bound */
  _min_y_bound = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an y-min boundary */
  for (const auto &c : _cells) {

    double cell_min_y = -std::numeric_limits<double>::infinity();
    boundaryType cell_min_y_bound = BOUNDARY_NONE;
    for (const auto &hs : *(c.second)) {
      surf = hs.surface;
      halfspace = hs.halfspace;

      if (surf->getSurfaceType() == YPLANE && halfspace == +1 &&
        surf->getBoundaryType() != BOUNDARY_NONE) {
        if (surf->getMinY(halfspace) > cell_min_y) {
          cell_min_y = surf->getMinY(halfspace);
          cell_min_y_bound = surf->getBoundaryType();
        }
      }
    }
    if (cell_min_y_bound != BOUNDARY_NONE && cell_min_y < min_y) {
      min_y = cell_min_y;
      _min_y_bound = cell_min_y_bound;
    }
  }

  /* If a y-min boundary was not found, get the y-min from the bounding boxes
   * of the cells */
  if (min_y == std::numeric_limits<double>::infinity()) {
    double cell_min_y = min_y;
    boundaryType cell_min_y_bound = BOUNDARY_NONE;
    for (const auto &c : _cells) {
      cell_min_y = c.second->getMinY();
      cell_min_y_bound = c.second->getMinYBoundaryType();
      if (cell_min_y_bound != BOUNDARY_NONE && cell_min_y < min_y) {
        min_y = cell_min_y;
        _min_y_bound = cell_min_y_bound;
      }
    }
  }

  _min_y = min_y;

  /* Calculate the maximum reachable y-coordinate in the geometry and store it
   * in _max_y */
  double max_y = -std::numeric_limits<double>::infinity();

  /* Calculate the boundary condition at the maximum reachable y-coordinate in
   * the Universe and store it in _max_y_bound */
  _max_y_bound = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an y-max boundary */
  for (const auto &c : _cells) {

    double cell_max_y = std::numeric_limits<double>::infinity();
    boundaryType cell_max_y_bound = BOUNDARY_NONE;
    for (const auto &hs : *(c.second)) {
      surf = hs.surface;
      halfspace = hs.halfspace;

      if (surf->getSurfaceType() == YPLANE && halfspace == -1 &&
        surf->getBoundaryType() != BOUNDARY_NONE) {
        if (surf->getMaxY(halfspace) < cell_max_y) {
          cell_max_y = surf->getMaxY(halfspace);
          cell_max_y_bound = surf->getBoundaryType();
        }
      }
    }
    if (cell_max_y_bound != BOUNDARY_NONE && cell_max_y > max_y) {
      max_y = cell_max_y;
      _max_y_bound = cell_max_y_bound;
    }
  }

  /* If a y-max boundary was not found, get the y-max from the bounding boxes
   * of the cells */
  if (max_y == -std::numeric_limits<double>::infinity()) {
    double cell_max_y = max_y;
    boundaryType cell_max_y_bound = BOUNDARY_NONE;
    for (const auto &c : _cells) {
      cell_max_y = c.second->getMaxY();
      cell_max_y_bound = c.second->getMaxYBoundaryType();
      if (cell_max_y_bound != BOUNDARY_NONE && cell_max_y > max_y) {
        max_y = cell_max_y;
        _max_y_bound = cell_max_y_bound;
      }
    }
  }

  _max_y = max_y;

  /* Calculate the minimum reachable z-coordinate in the geometry and store it
   * in _min_z */
  double min_z = std::numeric_limits<double>::infinity();

  /* Calculate the boundary condition at the minimum reachable z-coordinate in
   * the Universe and store it in _min_z_bound */
  _min_z_bound = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an z-min boundary */
  for (const auto &c : _cells) {

    double cell_min_z = -std::numeric_limits<double>::infinity();
    boundaryType cell_min_z_bound = BOUNDARY_NONE;
    for (const auto &hs : *(c.second)) {
      surf = hs.surface;
      halfspace = hs.halfspace;

      if (surf->getSurfaceType() == ZPLANE && halfspace == +1 &&
          surf->getBoundaryType() != BOUNDARY_NONE) {
        if (surf->getMinZ(halfspace) > cell_min_z) {
          cell_min_z = surf->getMinZ(halfspace);
          cell_min_z_bound = surf->getBoundaryType();
        }
      }
    }
    if (cell_min_z_bound != BOUNDARY_NONE && cell_min_z < min_z) {
      min_z = cell_min_z;
      _min_z_bound = cell_min_z_bound;
    }
  }

  /* If a z-min boundary was not found, get the z-min from the bounding boxes
   * of the cells */
  if (min_z == std::numeric_limits<double>::infinity()) {
    double cell_min_z = min_z;
    boundaryType cell_min_z_bound = BOUNDARY_NONE;
    for (const auto &c : _cells) {
      cell_min_z = c.second->getMinZ();
      cell_min_z_bound = c.second->getMinZBoundaryType();
      if (cell_min_z_bound != BOUNDARY_NONE && cell_min_z < min_z) {
        min_z = cell_min_z;
        _min_z_bound = cell_min_z_bound;
      }
    }
  }

  _min_z = min_z;

  /* Calculate the maximum reachable z-coordinate in the geometry and store it
   * in _max_z */
  double max_z = -std::numeric_limits<double>::infinity();

  /* Calculate the boundary condition at the maximum
  * reachable z-coordinate in the Universe and store it in _max_z_bound */
  _max_z_bound = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an z-max boundary */
  for (const auto &c : _cells) {

    double cell_max_z = std::numeric_limits<double>::infinity();
    boundaryType cell_max_z_bound = BOUNDARY_NONE;
    for (const auto &hs : *(c.second)) {
      surf = hs.surface;
      halfspace = hs.halfspace;

      if (surf->getSurfaceType() == ZPLANE && halfspace == -1 &&
          surf->getBoundaryType() != BOUNDARY_NONE) {
        if (surf->getMaxZ(halfspace) < cell_max_z) {
          cell_max_z = surf->getMaxZ(halfspace);
          cell_max_z_bound = surf->getBoundaryType();
        }
      }
    }
    if (cell_max_z_bound != BOUNDARY_NONE && cell_max_z > max_z) {
      max_z = cell_max_z;
      _max_z_bound = cell_max_z_bound;
    }
  }

  /* If a z-max boundary was not found, get the z-max from the bounding boxes
   * of the cells */
  if (max_z == -std::numeric_limits<double>::infinity()) {
    double cell_max_z = max_z;
    boundaryType cell_max_z_bound = BOUNDARY_NONE;
    for (const auto &c : _cells) {
      cell_max_z = c.second->getMaxZ();
      cell_max_z_bound = c.second->getMaxZBoundaryType();
      if (cell_max_z_bound != BOUNDARY_NONE && cell_max_z > max_z) {
        max_z = cell_max_z;
        _max_z_bound = cell_max_z_bound;
      }
    }
  }

  _max_z = max_z;
  _boundaries_inspected = true;

  // Fix Z bounds for 2D cases
  if (_min_z > FLT_INFINITY)   //判断 _min_z 是否未被初始化更新为一个有效的值,因为_min_z的一开始赋值就为无穷大
    _min_z = -std::numeric_limits<double>::infinity();
  if (_max_z < -FLT_INFINITY)  //判断 _max_z 是否未被初始化
    _max_z = std::numeric_limits<double>::infinity();
}


/**
  * @brief Sets _boundaries_not_updated to true so boundaries will be
  *        recalculated if needed.
  */
void Universe::resetBoundaries() {
  _boundaries_inspected = false;
}


/**
 * @brief Determines whether a Point is contained inside a Universe.
 * @details Queries each Cell in the Universe to determine if the Point
 *          is within the Universe. This point is only inside the Universe
 *          if it is inside one of the Cells.
 * @param point a pointer to a Point
 * @returns true if the Point is inside the Universe; otherwise false
 */
bool Universe::containsPoint(Point* point) {

  /* Loop over all Cells */
  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (iter->second->containsPoint(point))
      return true;
  }

  return false;
}


/**
 * @description 获取cell在输入卡中的id
 * @return inputId
 *
 */
int Universe::getInputId(){
  return _input_id;
}

} /* namespace antmoc */
