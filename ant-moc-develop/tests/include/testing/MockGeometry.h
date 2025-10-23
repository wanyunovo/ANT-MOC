/**
 * @file testing/MockGeometry.h
 * @brief Mock Geometry
 * @date April 17, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#ifndef MOCKGEOMETRY_H_
#define MOCKGEOMETRY_H_

#include "gmock/gmock.h"
#include "antmoc/Geometry.h"

namespace antmoc
{

class MockGeometry : public Geometry {
 public:
  
  MockGeometry() : Geometry() {}

  /* Get parameters */
  MOCK_METHOD0(getNumXModules, int());
  MOCK_METHOD0(getNumYModules, int());
  MOCK_METHOD0(getNumZModules, int());
  MOCK_METHOD0(getWidthX, double());
  MOCK_METHOD0(getWidthY, double());
  MOCK_METHOD0(getWidthZ, double());
  MOCK_METHOD0(getMinX, double());
  MOCK_METHOD0(getMaxX, double());
  MOCK_METHOD0(getMinY, double());
  MOCK_METHOD0(getMaxY, double());
  MOCK_METHOD0(getMinZ, double());
  MOCK_METHOD0(getMaxZ, double());

  MOCK_METHOD0(getGlobalMinX, double());
  MOCK_METHOD0(getGlobalMaxX, double());
  MOCK_METHOD0(getGlobalMinY, double());
  MOCK_METHOD0(getGlobalMaxY, double());
  MOCK_METHOD0(getGlobalMinZ, double());
  MOCK_METHOD0(getGlobalMaxZ, double());

  MOCK_METHOD0(getMinXBoundaryType, boundaryType());
  MOCK_METHOD0(getMaxXBoundaryType, boundaryType());
  MOCK_METHOD0(getMinYBoundaryType, boundaryType());
  MOCK_METHOD0(getMaxYBoundaryType, boundaryType());
  MOCK_METHOD0(getMinZBoundaryType, boundaryType());
  MOCK_METHOD0(getMaxZBoundaryType, boundaryType());

};

} /* namespace antmoc */

#endif
