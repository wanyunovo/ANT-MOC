ANT-MOC
=======

ANT-MOC (or antmoc) is a 3-D neutron transport code based on OpenMOC.

The code is maintained by students in *HPC&DE Lab, University of Science and Technology Beijing (USTB)*.

## Table of contents

- [Prerequisites](#prerequisites)
  - [Supported compilers](#supported-compilers)
  - [Useful links](#useful-links)
- [How to build antmoc](#how-to-build-antmoc)
  - [Python packages (experimental)](#python-packages-experimental)
  - [Environment variables](#environment-variables)
  - [CMake options](#cmake-options)
  - [Build dependencies manually](#build-dependencies-manually)
  - [Host configs](#host-configs)
- [How to run antmoc](#how-to-run-antmoc)
  - [Run cases](#run-cases)
  - [Run tests](#run-tests)
- [How to visualize outputs](#how-to-visualize-outputs)
  - [VTK data in XML](#vtk-data-in-xml)
  - [HDF5 data file](#hdf5-data-file)
- [Scheduling systems](#scheduling-systems)
- [Citing](#citing)
- [License](#License)


## Prerequisites

| Name           | Minimum Requirement          | Optional | Subproject | Info                                                         |
| -------------- | ---------------------------- | -------- | ---------- | ------------------------------------------------------------ |
| C/C++ compiler | with c++11 support           |          |            | see [C++ compiler support](https://en.cppreference.com/w/cpp/compiler_support#cpp11) |
| OpenMP         | 3.1                          |          |            | see [OpenMP Compilers & Tools](https://www.openmp.org/resources/openmp-compilers-tools/) |
| ROCm HIP       | 3.9                          | Yes      |            | see [ROCm HIP](https://github.com/ROCm-Developer-Tools/HIP)  |
| rocThrust      | (shipped with ROCm)          | Yes      |            | see [ROCm Thrust](https://github.com/ROCmSoftwarePlatform/rocThrust) |
| MPI            | 3.0                          | Yes      |            | see [MPI 3.0](https://www.mpi-forum.org/mpi-30/)             |
| CMake          | 3.16                         |          |            | see [Modern CMake](https://cliutils.gitlab.io/modern-cmake/) |
| {fmt}          | 6.0.0                        |          | Yes        | see [fmtlib/fmt](https://github.com/fmtlib/fmt)              |
| TinyXML2       | 7.0.0                        |          | Yes        | see [leethomason/tinyxml2](https://github.com/leethomason/tinyxml2) |
| toml11         | 3.5.0                        |          | Yes        | see [ToruNiina/toml11](https://github.com/ToruNiina/toml11)  |
| HDF5           | 1.10.0 (with optional PHDF5) |          | Yes        | see [The HDF Group](https://www.hdfgroup.org/solutions/hdf5/) |
| Better Enums   | 0.11.2                       |          | Internal   | see [aantron/better-enums](https://github.com/aantron/better-enums) |
| cxxopts        | 2.2.0                        |          | Internal   | see [jarro2783/cxxopts](https://github.com/jarro2783/cxxopts) |
| googletest     | 1.10.0                       | Yes      | Yes        | see [google/googletest](https://github.com/google/googletest) |
| lcov           | 1.14                         | Yes      |            | see [linux-test-project/lcov](https://github.com/linux-test-project/lcov) |
| PyEVTK         | 1.1.1                        | Yes      |            | see [pyscience-projects/pyevtk](https://github.com/pyscience-projects/pyevtk) |

> A subproject is a dependency which will be downloaded and built automatically when it is missing. If a subproject is marked as internal, its source is a part of the ANT-MOC repository.
>
> Please make sure the environment variable `CMAKE_PREFIX_PATH` has been set correctly, or the dependencies may not be found by CMake.

### Supported compilers

ANT-MOC has been tested on several x86_64 platforms. The compilers listed below should be able to build ANT-MOC correctly.

| Compilers/Toolchains | Executables    | Tested Versions   |
| -------------------- | -------------- | ----------------- |
| GCC                  | gcc, g++       | 5, 6, 7, 8, 9, 10 |
| LLVM                 | clang, clang++ | 8, 9, 10          |
| Intel C/C++ Compiler | icc, icpc      | 18.0              |
| HIP-Clang            | hipcc          | 3.9, 3.10, 4.2    |

> MPI wrappers `mpicc/mpicxx` could also be used to build ANT-MOC.

### Useful links

- [USTB HPCer](https://hpcde.github.io/cluster-docs/)
- [Mirror for antmoc dependencies](https://gitee.com/antmoc)


## How to build antmoc

It is highly recommended that you clone the code from git.

```bash
$ git clone git@gitlab.com:HPCer/neutronics/ant-moc.git
```

Once the repository is cloned, you can create a new directory and build the program by CMake,
which is called an *out-of-source* build.

Build and install antmoc, explicitly calling make

```bash
$ mkdir build && cd build
$ cmake ..
$ make
$ make install
```

or implicitly calling make 

```bash
$ cmake -S . -B build
$ cmake --build build
$ cmake --install build --prefix /path/to/install
```

Note that an *in-source* build is not allowed, whereas you can take advantage of out-of-source builds to build antmoc with different configurations.

```bash
# build with testing
$ cmake -S . -B build -DENABLE_TESTS=ON
$ cmake --build build

# build with debugging information
$ cmake -S . -B build -DENABLE_DEBUG=ON -DENABLE_ALL_WARNINGS=ON
$ cmake --build build

# build with MPI and HIP
$ cmake -S . -B build -DENABLE_MPI=ON -DENABLE_HIP=ON -DCMAKE_C_COMPILER=hipcc -DCMAKE_CXX_COMPILER=hipcc
$ cmake --build build
```

### Python packages (experimental)

There are packages for antmoc to facilitate data analysis. To take advantage of them,
please append the package directory to your environment variables.

For example, if your antmoc source tree is `$HOME/ant-moc`, then you can do

```bash
$ ANTMOC_DIR=$HOME/ant-moc
$ export PYTHONPATH=$ANTMOC_DIR/tools/pylib:$PYTHONPATH
```

Packages could now be imported in your scripts. See `./scripts/extract_records.py` for an example.

External packages:

- [antmoc-mgxs](https://gitlab.com/hpcer/neutronics/antmoc-mgxs), a package for MGXS manipulation.

### Environment variables

After the installation, you can put its location into your search paths to make it available everywhere.

```bash
$ export PATH=path/to/antmoc/bin:$PATH
$ export LD_LIBRARY_PATH=path/to/antmoc/lib:$LD_LIBRARY_PATH
```

### CMake options

We provide several options for building antmoc.

> See `CMakeLists.txt` for the full list.

For example, to enable MPI.

```bash
$ cmake -S . -B build -DENABLE_MPI=ON
```

There are several available options (the first one is the default value).

```
BUILD_SHARED_LIBS  = {ON, OFF}  # build shared or static libraries
ENABLE_MPI         = {OFF, ON}  # enable MPI support
ENABLE_HIP         = {OFF, ON}  # enable HIP support
ENABLE_CMFD        = {OFF, ON}  # enable CMFD
ENABLE_CYCLIC      = {ON, OFF}  # enable cyclic track
ENABLE_TESTS       = {OFF, ON}  # build tests
ENABLE_COVERAGE    = {OFF, ON}  # build a target for code coverage
```

### Build dependencies manually

Some third-party dependencies could be built with antmoc. See `CMakeLists.txt` for details.

```bash
cmake --build build -DUSE_INTERNAL_ALL=ON
```

### Host configs

Host configs are used to build antmoc with a specific set of options on specific machines.

For example, to build antmoc with HIP on Xiandao-1,

```bash
cmake --build build -C host-configs/xiandao-1/hip.cmake
```

## How to run antmoc

There are lots of runtime options available. Options may be used as command line arguments or written in a '.toml' file.

```bash
# To check available command line options, run
$ antmoc --help

# To check available options in a TOML file, run
$ antmoc --help-settings
```

You can also run antmoc from the source directory if you have made an out-of-source build.
Assume that you have built antmoc in `./build/`, then you can invoke the executable like this

```bash
$ ./build/antmoc --help
```

### Run cases

Pre-installed cases are located in directory `cases/`. You can run a case by specify the corresponding setting file.

```bash
$ antmoc -c cases/c5g7/simple-lattice/settings.xml
```

Note that you can run cases in your build tree because the option `ANTMOC_BUILD_CASES` defaults to `ON`,
which means the directory `cases/` in your source tree will be copied into your build tree.

### Run tests

To run tests, you have to build antmoc with an extra option enabled. After that, you can call CTest to run all of the tests.

```bash
$ cmake -S . -B build -DENABLE_TESTS=ON
$ cmake --build build
$ cd build
$ ctest
```

Test coverage is provided by an extra target, which could be enabled by an option.

```bash
$ cmake -S . -B build -DENABLE_COVERAGE=ON -DENABLE_TESTS=ON
$ cmake --build build
$ cd build
$ make coverage
```

Extra arguments may be provided to googletest to facilitate testing. This can be done by either providing them to CMake or to individual tests.
Command line options `--update-results` and `--visualize` can be set for regression testing. For example,

```bash
# Update test data
$ ./path/to/test --update-results
# Output visualization data
$ ./path/to/test --visualize
```

which may also be set globally

```bash
$ cmake .. -DENABLE_TESTS=ON -DENABLE_TESTS_UPDATE=ON -DENABLE_TESTS_VIS=ON
```

## How to visualize outputs

There are two types of files generated by ANT-MOC.

### VTK data in XML

If a tally mesh has been specified in the command line input, ANT-MOC will produce an XML file consisting
of VTK data arrays. Each of the reaction rates is stored separately in a XML node. The precision of
floating-point numbers is fixed and cannot be changed by command line input. But it is sufficient for
visualization with VTK.

Two types of tally mesh are supported by ANT-MOC: rectangular tally mesh and hexagonal tally mesh.
A rectangular tally mesh consists of rectangles (2-D) or hexahedrons (3-D), which is defined by three
arguments `widths_x, widths_y, widths_z` indicating its widths along each of the axis. A hexagonal tally
mesh consists of hexagons (2-D) or hexagonal prisms (3-D), which is defined four arguments
`num_r, width_r, widths_z, orientation` indicating the number of cells, the radial pitch, the cell widths
along the z-axis, and the lattice orientation respectively. Please run the program with `-h` or `--help`
for more details.

It is tempting to do some analysis on the VTK file generated by ANT-MOC. A VTK-based visualization toolkit
such as ParaView is highly recommended. In case you want to write your own scripts to manipulate the data,
we provide a Python script  `utils/vtu2xls.py` which shows how you can extract the data you are interested in.

This script extracts data on a plane and then writes it to an XLS file. It expects no more than 2 arguments.
The first argument is the path to a .vtu file, and the second one is the index of a z-section where the data resides.

```bash
$ python tools/vtu2xls.py "reaction_rates.vtu"
$ python tools/vtu2xls.py "reaction_rates.vtu" "0"
$ python tools/vtu2xls.py "reaction_rates.vtu" "0"
```

### HDF5 data file

ANT-MOC has its own format for visualization data. Information of FSRs and tracks are stored in HDF5 files,
and can be converted to VTK files subsequently. For example, if a simulation has been performed and a file
named `fsr_data.h5` has been produced, we can simply run the Python script to complete the conversion.

```bash
$ python tools/data2vtk.py "fsr_data.h5"
```

The format of our data file is as following
```
  group      /
  group      /FSR
  group      /FSR/Points
  dataset    /FSR/Points/X
  dataset    /FSR/Points/Y
  dataset    /FSR/Points/Z
  group      /FSR/Centroids
  dataset    /FSR/Centroids/X
  dataset    /FSR/Centroids/Y
  dataset    /FSR/Centroids/Z
  dataset    /FSR/Volumes
  group      /FSR/Fission RX
  dataset    /FSR/Fission RX/g1
  dataset    /FSR/Fission RX/g2
  dataset    /FSR/Fission RX/sum
  group      /FSR/Scalar Flux
  dataset    /FSR/Scalar Flux/g1
  dataset    /FSR/Scalar Flux/g2
  dataset    /FSR/Scalar Flux/sum
  ...
```

## Scheduling Systems

It's better to put all of your submission scripts in a subdirectory and name them with special suffixes.
In this repository, scripts are kept in `./scripts`. Each SLURM script has the suffix `.slurm`.

```bash
# Build antmoc when you are in the root directory of the source tree.
$ chmod u+x ./scripts/build-antmoc.sh
$ ./scripts/build-antmoc.sh

# Submit the script
$ sbatch ./scripts/submit-antmoc.slurm
```

## Citing

TBD

## License

MIT

