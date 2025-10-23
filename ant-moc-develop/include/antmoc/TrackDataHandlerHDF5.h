/// \file TrackDataHandlerHDF5.h
/// \brief An I/O class for track data
/// \date March 12, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef TRACKDATA_HANDLER_HDF5_H_
#define TRACKDATA_HANDLER_HDF5_H_


#include <memory>
#include <set>
#include <string>
#include <vector>

#include "antmoc/HDF5Handler.h"


namespace antmoc
{

class TrackGenerator;

///---------------------------------------------------------------------
/// \class TrackDataHandlerHDF5
/// \brief Methods for writing HDF5 files
///---------------------------------------------------------------------
class TrackDataHandlerHDF5 : public HDF5Handler {

protected:

  using TrackGeneratorPtr = std::shared_ptr<TrackGenerator>;

  ///< The trackgenerator from which data is extracted
  TrackGeneratorPtr _track_generator;

public:

  TrackDataHandlerHDF5(std::string file, HDF5Mode mode, _track_generator tg);
  ~TrackDataHandlerHDF5() = default;

  // To be implemented

};


} // namespace antmoc

#endif  // TRACKDATA_HANDLER_HDF5_H_
