#ifndef GEOMETRY_TESTUTILS_H_
#define GEOMETRY_TESTUTILS_H_

#include <algorithm>
#include <array>
#include <sstream>

#include "testing/test_utils.h"
#include "antmoc/Cell.h"
#include "antmoc/Lattice.h"
#include "antmoc/Material.h"
#include "antmoc/MaterialHandlerHDF5.h"
#include "antmoc/Surface.h"
#include "antmoc/Universe.h"
#include "antmoc/lattice_utils.h"
#include "antmoc/PyVector.h"

namespace antmoc {

using Widths3D = std::vector<WidthVec>; ///< A vector of widths in x,y,z directions

using std::string;

/**
 * @brief Test if a universe is built correctly
 */
inline void testCell(Cell           *cell,
                     const int      id,
                     const string   &name,
                     bool           is_material,
                     const string   &fill_name,
                     const string   *region   = nullptr,
                     const IntPyVec *sf_ids   = nullptr,
                     const int      n_sectors = 0,
                     const int      n_rings   = 0) {
  // Test cell id and name
  EXPECT_EQ   (id,           cell->getInputId());
  EXPECT_STREQ(name.c_str(), cell->getName());

  // Test fill
  if (is_material) {
    const auto m = cell->getFillMaterial();
    ASSERT_TRUE(m);
    EXPECT_STREQ(fill_name.c_str(), m->getName());
  }
  else {
    const auto u = cell->getFillUniverse();
    ASSERT_TRUE(u);
    EXPECT_STREQ(fill_name.c_str(), u->getName());
  }

  // Test the number of sectors and rings
  EXPECT_EQ(n_sectors, cell->getNumSectors());
  EXPECT_EQ(n_rings,   cell->getNumRings());

  // Test the region (optional)
  if (region) {
    EXPECT_EQ(*region, cell->getRegion());
  }

  // Test the ids of cells (optional)
  if (sf_ids) {
    const auto &surfs = cell->getSurfaces();
    ASSERT_EQ(sf_ids->size(), surfs.size());

    for (auto i : *sf_ids)
      EXPECT_TRUE(surfs.find(i) != surfs.end())
        << "  Halfspace id " << i << " not found in cell " << name
        << " (id=" << id << ")\n";
  }
}


/**
 * @brief Test if a universe is built correctly
 */
inline void testUniverse(Universe       *universe,
                         const int      id,
                         const string   &name,
                         const IntPyVec *cell_ids = nullptr) {
  // Test universe id and name
  EXPECT_EQ   (id,           universe->getInputId());
  EXPECT_STREQ(name.c_str(), universe->getName());

  // Test the ids of cells (optional)
  if (cell_ids) {
    const auto &c_ids = *cell_ids;

    const auto &cells = universe->getCells();
    ASSERT_EQ(c_ids.size(), cells.size());

    // Read input ids to a vector
    IntPyVec oracles;
    for (auto &c : cells)
      oracles.push_back(c.second->getInputId());

    for (auto i : c_ids)
      EXPECT_TRUE(std::find(oracles.begin(), oracles.end(), i) != oracles.end())
        << "  Cell id " << i << " not found in universe " << name
        << " (id=" << id << ")\n";
  }
}


/**
 * @brief Test if a lattice is built correctly
 */
inline void testRecLattice(Lattice        *reclattice,
                           const int      id,
                           const string   &name,
                           const IntPyVec3D *layout_oracles = nullptr,
                           const Widths3D *widths_oracles = nullptr,
                           const IntPyVec *universe_ids   = nullptr) {
  auto lattice = dynamic_cast<RecLattice*>(reclattice);
  ASSERT_TRUE(lattice);

  // Test lattice layout type
  ASSERT_EQ(latticeType::Rectangle, lattice->getLatticeType());

  // Test lattice id and name
  EXPECT_EQ   (id,           lattice->getInputId());
  EXPECT_STREQ(name.c_str(), lattice->getName());

  // Test the ids of unique universes (optional)
  if (universe_ids) {
    const auto &u_ids = *universe_ids;

    const auto &universes = lattice->getUniqueUniverses();
    ASSERT_EQ(u_ids.size(), universes.size());

    // Read input ids to a vector
    IntPyVec oracles;
    for (auto &u : universes)
      oracles.push_back(u.second->getInputId());

    for (auto i : u_ids)
      EXPECT_TRUE(std::find(oracles.begin(), oracles.end(), i) != oracles.end())
        << "  Universe id " << i << " not found in lattice " << name
        << " (id=" << id << ")\n";
  }

  // Test lattice layout (optional)
  if (layout_oracles) {
    const auto &layout = *layout_oracles;

    // Get the numbers of lattice cells
    int nz = layout.size();
    ASSERT_TRUE(nz > 0);  // min(nz) = 1

    int ny = layout[0].size();
    ASSERT_TRUE(ny > 0);  // min(ny) = 1

    int nx = layout[0][0].size();
    ASSERT_TRUE(nx > 0);  // min(nx) = 1

    // Check the numbers of lattice cells
    ASSERT_EQ(nx, lattice->getNumX());
    ASSERT_EQ(ny, lattice->getNumY());
    ASSERT_EQ(nz, lattice->getNumZ());

    // Generate strings for the layout
    std::stringstream ss_oracle, ss_output;
    for (int z = 0; z < nz; z++) {
      for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
          Universe *u = nullptr;
          ASSERT_NO_THROW(u = lattice->getUniverse(x, ny-1-y, nz-1-z));
          if (u)
            ss_output << u->getInputId() << '\t';
          else
            ss_output << -1 << '\t';  // nullptr
          ss_oracle << layout[z][y][x] << '\t';
        }
        ss_oracle << '\n';
        ss_output << '\n';
      }
      ss_oracle << "\n\n";
      ss_output << "\n\n";
    }

    // Generate hash values
    std::hash<std::string> hash_f;
    ASSERT_EQ(hash_f(ss_oracle.str()), hash_f(ss_output.str()))
      << nx << "x" << ny << "x" << nz << " layout:\n" << ss_oracle.str()
      << "But ouput:\n" << ss_output.str();
  }

  // Test lattice widths (optional)
  if (widths_oracles) {
    const auto &oracles = *widths_oracles;
    std::vector<std::vector<double>> ws = {
      lattice->getWidthsX(),
      lattice->getWidthsY(),
      lattice->getWidthsZ()
    };

    for (size_t i = 0; i < oracles.size(); ++i) {
      ASSERT_EQ(oracles[i].size(), ws[i].size());
      for (size_t e = 0; e < oracles[i].size(); ++e)
        EXPECT_DOUBLE_EQ(oracles[i].at(e), ws[i].at(e));
    }
  }
}


/**
 * @brief Test if a HexLattice is built correctly
 */
inline void testHexLattice(Lattice        *hexlattice,
                           const int      id,
                           const string   &name,
                           const IntPyVec *layout_oracles = nullptr,
                           const double   *width_r_oracle = nullptr,
                           const WidthVec *widths_z_oracles = nullptr,
                           const int      *num_r_oracle = nullptr,
                           const IntPyVec *universe_ids   = nullptr) {
  auto lattice = dynamic_cast<HexLattice*>(hexlattice);
  ASSERT_TRUE(lattice);

  // Test lattice layout type
  ASSERT_EQ(latticeType::Hexagon, lattice->getLatticeType());

  // Test lattice id and name
  EXPECT_EQ   (id,           lattice->getInputId());
  EXPECT_STREQ(name.c_str(), lattice->getName());

  // Test the ids of unique universes (optional)
  if (universe_ids) {
    const auto &u_ids = *universe_ids;

    const auto &universes = lattice->getUniqueUniverses();
    ASSERT_EQ(u_ids.size(), universes.size());

    // Read input ids to a vector
    IntPyVec oracles;
    for (auto &u : universes)
      oracles.push_back(u.second->getInputId());

    for (auto i : u_ids)
      EXPECT_TRUE(std::find(oracles.begin(), oracles.end(), i) != oracles.end())
        << "  Universe id " << i << " not found in lattice " << name
        << " (id=" << id << ")\n";
  }

  // Test lattice cells (optional)
  if (num_r_oracle) {
    EXPECT_EQ(*num_r_oracle, lattice->getNumR());
  }

  // Test lattice widths_r (optional)
  if (width_r_oracle) {
    const auto oracle = *width_r_oracle;
    EXPECT_DOUBLE_EQ(oracle, lattice->getWidthR());
  }

  // Test lattice widths_z (optional)
  if (widths_z_oracles) {
    const auto &oracles = *widths_z_oracles;
    const std::vector<double> &ws = lattice->getWidthsZ();

    for (size_t i = 0; i < oracles.size(); ++i)
      EXPECT_DOUBLE_EQ(oracles[i], ws[i]);
  }

  // Test lattice layout (optional)
  if (layout_oracles) {
    const auto &layout = *layout_oracles;

    // Generate strings for the layout
    std::stringstream ss_oracle, ss_output;
    for (auto i : layout)
      ss_oracle << i << ' ';

    for (auto u : *lattice) {
      if (u)
        ss_output << u->getInputId() << ' ';
      else
        ss_output << -1 << ' ';  // nullptr
    }
    ss_oracle << '\n';
    ss_output << '\n';

    // Generate hash values
    std::hash<std::string> hash_f;
    ASSERT_EQ(hash_f(ss_oracle.str()), hash_f(ss_output.str()))
      << " layout:\n" << ss_oracle.str()
      << "But ouput:\n" << ss_output.str();
  }
}


} // namespace antmoc

#endif
