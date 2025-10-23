#include "antmoc/FSRDataTypes.h"
#include "antmoc/LocalCoords.h"
#include "antmoc/Material.h"
#include "antmoc/Point.h"


namespace antmoc {

FSRData::FSRData() :
  _fsr_id(0),
  _cmfd_cell(0),
  _mat_id(0),
  _point(NULL),
  _centroid(NULL)
{}


FSRData::~FSRData() {
  if (_point != NULL)
    delete _point;

  if (_centroid != NULL)
    delete _centroid;
}


ExtrudedFSR::ExtrudedFSR() :
  _mesh(NULL),
  _fsr_id(0),
  _fsr_ids(NULL),
  _materials(NULL),
  _num_fsrs(0),
  _coords(NULL)
{}


ExtrudedFSR::~ExtrudedFSR() {
  delete [] _mesh;
  delete [] _fsr_ids;
  delete [] _materials;
  delete _coords;
}

} // namespace antmoc
