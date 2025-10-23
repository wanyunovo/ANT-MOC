#include "testing/test_utils.h"
#include "antmoc/log.h"

class ANTMOCEnvironment : public ::testing::Environment {
  public:
  virtual void SetUp() {
#ifdef ENABLE_MPI_
    char **argv;
    int argc = 0;
    int provided;
    int mpiError = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    ASSERT_FALSE(mpiError);

    // FIXME: this routine doesn't actually perform the domain decomposition
    antmoc::mpi::defaultInitialize();
#endif
    antmoc::log::initialize();
    antmoc::log::set_level("test");
  }

  virtual void TearDown() {
#ifdef ENABLE_MPI_
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
#endif
  }

  virtual ~ANTMOCEnvironment() {}
};

GTEST_API_ int main(int argc, char **argv) {
  testing::InitGoogleMock(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new ANTMOCEnvironment);

  // Parses custom options
  auto &test_options = TestOptions::get();
  test_options.parse(argc, argv);

  return RUN_ALL_TESTS();
}
