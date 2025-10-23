/// \brief Test lattice cell center computing
/// \date March 10, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "testing/geometry_testutils.h"
#include "antmoc/Lattice.h"
#include "antmoc/math_utils.h"

using namespace antmoc;

namespace {

// Test fixture
class test_Lattice_cellCenters : public testing::Test { };


/// \brief Compute lattice cell centers for a uniform RecLattice
TEST_F(test_Lattice_cellCenters, uniformRecLattice2x2x2) {
  ROOT_ONLY();

  // Construct a 1.2*1.2*1.2 lattice
  RecLattice lattice;
  lattice.setWidth(0.6, 0.6, 0.6);
  lattice.setNumX(2);
  lattice.setNumY(2);
  lattice.setNumZ(2);
  lattice.computeSizes();

  // Lattice cells are ordered in x-y-z order
  PyVector<Point> oracles = {
    {-0.3, -0.3, -0.3}, // lower left at z=0
    {+0.3, -0.3, -0.3}, // lower right at z=0
    {-0.3, +0.3, -0.3}, // upper left at z=0
    {+0.3, +0.3, -0.3}, // upper right at z=0
    {-0.3, -0.3, +0.3}, // lower left at z=1
    {+0.3, -0.3, +0.3}, // lower right at z=1
    {-0.3, +0.3, +0.3}, // upper left at z=1
    {+0.3, +0.3, +0.3}, // upper right at z=1
  };

  for (size_t index = 0; index < oracles.size(); ++index) {
    Point &oracle = oracles[index];
    Point center = lattice.getLatticeCellCenter(index);

    EXPECT_EQ(center, oracle)
      << "Center: " << center.toString() << '\n'
      << "Oracle: " << oracle.toString()
      << std::endl;
  }
}


TEST_F(test_Lattice_cellCenters, nonUniformRecLattice2x2x2) {
  ROOT_ONLY();
  
  // Construct a 1.2*1.2*1.2 lattice
  RecLattice lattice;
  DoubleVec widths = {0.3, 0.9};
  lattice.setWidths(widths, widths, widths);
  lattice.setNumX(2);
  lattice.setNumY(2);
  lattice.setNumZ(2);
  lattice.computeSizes();

  // Lattice cells are ordered in x-y-z order
  PyVector<Point> oracles = {
    {-0.45, -0.45, -0.45}, // lower left at z=0
    {+0.15, -0.45, -0.45}, // lower right at z=0
    {-0.45, +0.15, -0.45}, // upper left at z=0
    {+0.15, +0.15, -0.45}, // upper right at z=0
    {-0.45, -0.45, +0.15}, // lower left at z=1
    {+0.15, -0.45, +0.15}, // lower right at z=1
    {-0.45, +0.15, +0.15}, // upper left at z=1
    {+0.15, +0.15, +0.15}, // upper right at z=1
  };

  for (size_t index = 0; index < oracles.size(); ++index) {
    Point &oracle = oracles[index];
    Point center = lattice.getLatticeCellCenter(index);

    EXPECT_EQ(center, oracle)
      << "Center: " << center.toString() << '\n'
      << "Oracle: " << oracle.toString()
      << std::endl;
  }
}


/// \brief Test a y-orientated HexLattice with nr = 2 and nz = 2
TEST_F(test_Lattice_cellCenters, yOrientatedHexLattice2x2) {
  ROOT_ONLY();

  // 2x2 hexagon lattice with radial pitch 1
  HexLattice lattice;
  double width_r = 1.;
  DoubleVec widths_z = {0.3, 0.9};

  lattice.setOrientation("y");
  lattice.setWidths(width_r, widths_z);
  lattice.setNumR(2);
  lattice.setNumZ(2);
  lattice.computeSizes();

  // Lattice cells are ordered in x-y-z order
  //      cell
  // cell      cell
  //      cell
  // cell      cell
  //      cell
  double cos30 = std::cos(M_PI/6);
  double sin30 = std::sin(M_PI/6);

  PyVector<Point> oracles = {
    { -cos30, -sin30-1, -0.45}, // lower left cell at (0,0,0) (invalid index)
    {      0,       -1, -0.45}, // bottom cell at (1,0,0)
    {  cos30,  sin30-1, -0.45}, // lower right cell at (2,0,0)
    { -cos30,   -sin30, -0.45}, // middle left cell at (0,1,0)
    {      0,        0, -0.45}, // center cell at (1,1,0)
    {  cos30,    sin30, -0.45}, // middle right cell at (2,1,0)
    { -cos30, -sin30+1, -0.45}, // upper left cell at (0,2,0)
    {      0,        1, -0.45}, // top cell at (1,2,0)
    {  cos30,  sin30+1, -0.45}, // upper right cell at (2,2,0) (invalid index)

    { -cos30, -sin30-1, +0.15}, // lower left cell at (0,0,1) (invalid index)
    {      0,       -1, +0.15}, // bottom cell at (1,0,1)
    {  cos30,  sin30-1, +0.15}, // lower right cell at (2,0,1)
    { -cos30,   -sin30, +0.15}, // middle left cell at (0,1,1)
    {      0,        0, +0.15}, // center cell at (1,1,1)
    {  cos30,    sin30, +0.15}, // middle right cell at (2,1,1)
    { -cos30, -sin30+1, +0.15}, // upper left cell at (0,2,1)
    {      0,        1, +0.15}, // top cell at (1,2,1)
    {  cos30,  sin30+1, +0.15}, // upper right cell at (2,2,1) (invalid index)
  };

  for (size_t index = 0; index < oracles.size(); ++index) {
    Point &oracle = oracles[index];
    Point center = lattice.getLatticeCellCenter(index);

    EXPECT_EQ(center, oracle)
      << "LatticeCell: " << index << '\n'
      << "Center: " << center.toString() << '\n'
      << "Oracle: " << oracle.toString()
      << std::endl;
  }
}


/// \brief Test an x-orientated HexLattice with nr = 2 and nz = 2
TEST_F(test_Lattice_cellCenters, xOrientatedHexLattice2x2) {
  ROOT_ONLY();

  // 2x2 hexagon lattice with radial pitch 1
  HexLattice lattice;
  double width_r = 1.;
  DoubleVec widths_z = {0.3, 0.9};

  lattice.setOrientation("x");
  lattice.setWidths(width_r, widths_z);
  lattice.setNumR(2);
  lattice.setNumZ(2);
  lattice.computeSizes();

  // Lattice cells are ordered in x-y-z order
  //      cell      cell
  // cell      cell      cell
  //      cell      cell
  double cos30 = std::cos(M_PI/6);
  double sin30 = std::sin(M_PI/6);

  PyVector<Point> oracles = {
    { -sin30-1, -cos30, -0.45}, // lower left cell at (0,0,0) (invalid index)
    {   -sin30, -cos30, -0.45}, // bottom cell at (1,0,0)
    {  1-sin30, -cos30, -0.45}, // lower right cell at (2,0,0)
    {       -1,      0, -0.45}, // middle left cell at (0,1,0)
    {        0,      0, -0.45}, // center cell at (1,1,0)
    {        1,      0, -0.45}, // middle right cell at (2,1,0)
    {  sin30-1,  cos30, -0.45}, // upper left cell at (0,2,0)
    {    sin30,  cos30, -0.45}, // top cell at (1,2,0)
    {  sin30+1,  cos30, -0.45}, // upper right cell at (2,2,0) (invalid index)
                
    { -sin30-1, -cos30, +0.15}, // lower left cell at (0,0,1) (invalid index)
    {   -sin30, -cos30, +0.15}, // bottom cell at (1,0,1)
    {  1-sin30, -cos30, +0.15}, // lower right cell at (2,0,1)
    {       -1,      0, +0.15}, // middle left cell at (0,1,1)
    {        0,      0, +0.15}, // center cell at (1,1,1)
    {        1,      0, +0.15}, // middle right cell at (2,1,1)
    {  sin30-1,  cos30, +0.15}, // upper left cell at (0,2,1)
    {    sin30,  cos30, +0.15}, // top cell at (1,2,1)
    {  sin30+1,  cos30, +0.15}, // upper right cell at (2,2,1) (invalid index)
  };

  for (size_t index = 0; index < oracles.size(); ++index) {
    Point &oracle = oracles[index];
    Point center = lattice.getLatticeCellCenter(index);

    EXPECT_EQ(center, oracle)
      << "LatticeCell: " << index << '\n'
      << "Center: " << center.toString() << '\n'
      << "Oracle: " << oracle.toString()
      << std::endl;
  }
}
} // namespace
