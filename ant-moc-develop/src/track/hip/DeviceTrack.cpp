#include "antmoc/hip/DeviceTrack.h"
#include "antmoc/Material.h"
#include "antmoc/Track.h"


namespace antmoc {
namespace hip {

__device__ void DeviceSegment::printString() {
  printf("Segment:\n"
         "length = %f, start = (%f, %f, %f)\n"
         "material_index = %d, region_id = %ld, track_idx = %d\n",
         _length, _start[0], _start[1], _start[2],
         _material_index, _region_id, _track_idx);
}


__host__ __device__
DeviceTrack::DeviceTrack():
  _uid(-1),
  _start{0, 0, 0},
  _phi(0),
  _segments(nullptr),
  _num_segments(0),
  _azim_index(-1),
  _xy_index(-1),
  _link_index(-1),
  _chain_index(-1),
  _track_next_fwd(-1),
  _track_next_bwd(-1),
  _track_prdc_fwd(-1),
  _track_prdc_bwd(-1),
  _track_refl_fwd(-1),
  _track_refl_bwd(-1),
  _next_fwd_is_fwd(false),
  _next_bwd_is_fwd(false),
  _bc_fwd(boundaryType::VACUUM),
  _bc_bwd(boundaryType::VACUUM)
  { }


__host__ __device__
DeviceTrack::~DeviceTrack() {
  if (_segments)
    free(_segments);
}


__device__ void DeviceTrack::printString() {
  printf("DeviceTrack:\n"
         "uid = %ld, phi = %f, num_segments = %d\n"
         "azim_index = %d, xy_index = %ld, link_index = %ld, chain_index = %ld\n",
         _uid, _phi, _num_segments,
         _azim_index, _xy_index, _link_index, _chain_index);
  printf("starting point:\n");
  _start.printString();

  // Print all segments. This is only for tests since the number of segments
  // is very large in real cases.
  if (_num_segments > 0) {
    printf("contained segments:\n");
    for (int i = 0; i < _num_segments; ++i)
      _segments[i].printString();
  }
}


//------------------------------------------------------------------------------
// Helper functions
//------------------------------------------------------------------------------
void trackCopyHtoD(DeviceTrack *track_d, Track *track_h,
                   std::map<int, int> &material_IDs_to_indices) {

  DeviceTrack new_track;

  new_track._uid = track_h->getUid();
  new_track._phi = track_h->getPhi();
  new_track._num_segments = track_h->getNumSegments();

  // Indices of the 2D track
  new_track._azim_index = track_h->getAzimIndex();
  new_track._xy_index = track_h->getXYIndex();
  new_track._link_index = track_h->getLinkIndex();
  new_track._chain_index = track_h->getChainIndex();

  // Track connections
  new_track._track_next_fwd = track_h->getTrackNextFwd();
  new_track._track_next_bwd = track_h->getTrackNextBwd();
  new_track._track_prdc_fwd = track_h->getTrackPrdcFwd();
  new_track._track_prdc_bwd = track_h->getTrackPrdcBwd();
  new_track._track_refl_fwd = track_h->getTrackReflFwd();
  new_track._track_refl_bwd = track_h->getTrackReflBwd();

  new_track._next_fwd_is_fwd = track_h->getNextFwdFwd();
  new_track._next_bwd_is_fwd = track_h->getNextBwdFwd();

  // Boundary conditions
  new_track._bc_fwd = track_h->getBCFwd();
  new_track._bc_bwd = track_h->getBCBwd();

  // Starting point
  new_track._start = DevicePoint(track_h->getStart()->getX(),
                                 track_h->getStart()->getY(),
                                 track_h->getStart()->getZ());

  // Copy host DeviceTrack to device.
  HIP_ASSERT( hipMemcpyHtoD(track_d, &new_track, sizeof(DeviceTrack)) );

  // Prepare segments for the track
  const int num_segments = track_h->getNumSegments();
  const size_t size_segments = num_segments * sizeof(DeviceSegment);
  DeviceSegment *host_segments = new DeviceSegment[num_segments];

  for (int s=0; s < num_segments; s++) {
    auto curr = track_h->getSegment(s);
    host_segments[s]._length    = curr->_length;
    host_segments[s]._region_id = curr->_region_id;
    host_segments[s]._track_idx = curr->_track_idx;
    host_segments[s]._start[0] = curr->_starting_position[0];
    host_segments[s]._start[1] = curr->_starting_position[1];
    host_segments[s]._start[2] = curr->_starting_position[2];
    host_segments[s]._material_index =
        material_IDs_to_indices[curr->_material->getId()];
  }

  // Copy segments from host to device.
  // We didn't set the pointer new_track._segments to dev_segments because it
  // will cause new_track, a host object, trying to free a device pointer in
  // its destructor.
  DeviceSegment *dev_segments;
  HIP_ASSERT( hipMalloc(&dev_segments, size_segments) );
  HIP_ASSERT( hipMemcpyHtoD(dev_segments, host_segments, size_segments) );
  HIP_ASSERT( hipMemcpyHtoD(&track_d->_segments, &dev_segments, sizeof(DeviceSegment*)) );

  delete [] host_segments;
}


} // namespace hip
} // namespace antmoc

