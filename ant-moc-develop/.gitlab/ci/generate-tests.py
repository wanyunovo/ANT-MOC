#!/usr/bin/env python3
#-------------------------------------------------------------------------------
# Preparing
#
# This script depends on Spack to validate specs. Jobs will be generated only
# for valid specs.
#-------------------------------------------------------------------------------
import sys, os, getopt
import subprocess

# set up default values
spack_exe = "/opt/spack/bin/spack"

# parse command line arguments
opts, args = getopt.getopt(sys.argv[1:], "h", ["spack="])

for name, value in opts:
  if name == "--spack":
    spack_exe = value

if os.path.exists(spack_exe):
  print("# Spack found: {}".format(spack_exe))
  print("# All specs will be checked. Jobs with missing specs won't be created.")
  check_specs = True
else:
  print("# Spack not found, specs won't be checked")
  check_specs = False

# CXX compiler and compiler specs
compilers = {
  "gcc": ("g++", "%gcc"),
  "clang": ("clang++", "%clang"),
  "hipcc": ("hipcc", "%gcc")
}

# mpi specs
mpi_specs = {
  "serial": "~mpi",
  "mpich": "+mpi^mpich",
  "openmpi": "+mpi^openmpi"
}

#-------------------------------------------------------------------------------
# Generating
#-------------------------------------------------------------------------------
def on_off(bool_value):
  return "ON" if bool_value else "OFF"

print( \
  "stages:\n" \
  "  - build\n" \
  "  - test\n" \
  "  - install\n" \
  "include:\n" \
  "  - local: .gitlab/ci/base.yml")

count_builds = 0
count_tests = 0

for cc in compilers.keys():
  for mpi in mpi_specs.keys():
    # use different config if dependencies will be built internally
    use_spec    = "antmoc {} {}".format(compilers[cc][1], mpi_specs[mpi])
    job_build   = "{}:{}:build".format(cc, mpi)
    job_test    = "{}:{}:test".format(cc,mpi)
    job_install = "{}:{}:install".format(cc,mpi)

    # use Spack to check the current spec
    if check_specs:
      result = subprocess.run([spack_exe, "find", use_spec], stdout=subprocess.PIPE).stdout.decode('utf-8')
      # skip jobs if the spec is missing
      if result.lower().find("no package") >= 0:
        continue

    use_mpi = (mpi != "serial")
    use_hip = (cc == "hipcc")

    # print jobs for building antmoc
    print('{}:'.format(job_build))
    print('  extends: .build-base')
    print('  variables:')
    print('    C_COMPILER: "{}"'.format(cc))
    print('    CXX_COMPILER: "{}"'.format(compilers[cc][0]))
    print('    USE_SPECS: "{}"'.format(use_spec))
    print('    ENABLE_MPI: "{}"'.format( on_off(use_mpi) ))
    print('    ENABLE_HIP: "{}"'.format( on_off(use_hip) ))
    count_builds += 1

    if cc != "hipcc" and (cc != "clang" or mpi != "openmpi"): # skip tests with hip or clang+openmpi for now
      # print jobs for testing antmoc
      print('{}:'.format(job_test))
      print('  extends: .test-base')
      print('  needs: ["{}"]'.format(job_build))
      print('  variables:')
      print('    USE_SPECS: "{}"'.format(use_spec))
      print('    CTEST_RANDOM: "{}"'.format( on_off(not use_mpi) ))
      count_tests += 1

      # print jobs for installing antmoc
      print('{}:'.format(job_install))
      print('  extends: .install-base')
      print('  needs: ["{}"]'.format(job_build))
      print('  variables:')
      print('    USE_SPECS: "{}"'.format(use_spec))

print("# Summary: {} builds, {} tests/installs".format(count_builds, count_tests))

