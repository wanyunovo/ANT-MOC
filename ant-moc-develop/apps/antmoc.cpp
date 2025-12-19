#include "antmoc/Cmfd.h"
#include "antmoc/ConfigInputFile.h"
#include "antmoc/debug_utils.h"
#include "antmoc/enum_types.h"
#include "antmoc/Factory.h"
#include "antmoc/file_utils.h"
#include "antmoc/FSRDataHandlerHDF5.h"
#include "antmoc/GeoInputXml.h"
#include "antmoc/MaterialHandlerHDF5.h"
#include "antmoc/Mesh.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/log.h"
#include "antmoc/Timer.h"
#include "antmoc/TrackGenerator3D.h"

#ifdef ENABLE_HIP_
#include "antmoc/hip/hip_info.h"
#endif

using namespace antmoc;

void printTimerReport(Timer &);
std::shared_ptr<Mesh> createSimpleMesh(Factory::ConfInputPtr, Factory::SolverPtr);
void dumpMeshData(Factory::ConfInputPtr, std::shared_ptr<Mesh>);
void dumpDataToH5File(Factory::ConfInputPtr, Factory::SolverPtr);

int main(int argc, char *argv[])
{

  /*
  并行环境初始化 (MPI)，这部分代码是高性能计算（HPC）中非常标准的“固定模式”。
  #ifdef ENABLE_MPI_ ... #endif: 这是一个预处理指令（宏定义）。含义：如果在编译代码时开启了 MPI（多进程并行计算）功能，编译器才会编译中间这几行代码；如果没有开启，这几行代码会被直接忽略，就像不存在一样。这让程序既可以在单机上跑，也可以在超级计算机上跑。
  MPI_Init_thread(...): 这是 MPI 标准库的启动函数。它负责建立进程间的通信环境。
  MPI_THREAD_SERIALIZED：告诉 MPI 系统，虽然我们可能使用多线程，但我们会保证同一时间只有一个线程调用 MPI 函数（这是为了线程安全）。
  mpi::defaultInitialize(): 这是一个项目自定义的辅助函数，用来初始化该软件内部的一些 MPI 相关的全局变量。
  */

#ifdef ENABLE_MPI_
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  // Initialize MPI static variables to default values - 多线程的 MPI 程序初始化 返回 MPI 实际提供的线程支持级别给provided
  mpi::defaultInitialize();
#endif

  //------------------------------------------------------------
  // Initialize the logger - 初始化日志系统
  log::initialize();

  Timer timer;
  timer.startTimer();

  //------------------------------------------------------------
  // Read simulation arguments - 读取仿真参数
  //------------------------------------------------------------
  timer.startTimer();
  // 这段代码中连续调用两次 timer.startTimer() 并不是重复或错误，而是为了实现嵌套计时（Nested Timing）。
  // Timer 类通常使用栈（Stack）结构来管理时间：第一次调用 (第 42 行)，这是为了开始统计整个 "Input" 阶段的总时间。对应的是代码末尾第 194 行的timer.stopTimer("Input");。
  // 第二次调用 (第 49 行) 则是为了开始统计 "Read Settings" 子阶段的时间。它对应的是第 70 行的 timer.stopTimer("Read Settings");。

  /*
  这里使用了工厂模式（Factory Pattern）。含义：你不需要知道对象是如何创建的复杂细节，只需要调用工厂方法，它就会根据输入参数（argc, argv）解析命令行，读取输入文件，并返回一个包含所有配置信息的对象conf
  argc 和 argv 是 C++ main 函数的标准参数，分别代表命令行参数的个数和内容。
  */
  auto conf = Factory::getConfInput<ConfigInputFile>(argc, argv);

  // Print the help message when needed - 在需要时打印帮助信息
  if (argc < 2)
  {
    conf->showHelp();
  }
  else
  {
    conf->showHelpAsNeeded();
  }

  // Set log level - 设置日志级别
  conf->initializeLogger();

  // Print the title - 打印标题
  log::title("ANT-MOC: an advanced neutron transport code from the CVR suite");

  // Report input arguments - 输出输入参数报告
  conf->printArgumentsReport();

  // Validate input parameters - 校验输入参数
  conf->validateArguments();

  timer.stopTimer("Read Settings"); // 意味着“读取配置”这个阶段结束了。

  // Print pid and hostname with log::profile - 使用 log::profile 打印进程号和主机名
#ifdef ENABLE_MPI_
  mpi::mpiBarrier();
  /*
  同步障碍（Barrier）。含义：这是一个非常重要的 MPI 概念。就像旅游团在景点门口集合一样，代码运行到这里时，所有进程必须都到达这一行，才能继续往下执行。作用：防止有的进程跑得太快，有的太慢，导致后面的输出信息乱序或者数据不同步。
  */
#endif
  printHostnames(); // 打印当前进程运行在哪台机器（节点）上。在超级计算机上，一个任务可能跨越几百台机器，这个信息用于确认任务分布情况。

#ifdef ENABLE_HIP_
  // 如果编译时开启了 HIP 支持（通常用于 AMD GPU 加速，类似于 NVIDIA 的 CUDA）。
  hip::printDeviceProperties();
#endif

#ifdef ENABLE_MPI_

  // Report MPI datatypes - 打印程序自定义的 MPI 数据结构信息
  mpi::showMPIDatatypes();

  //------------------------------------------------------------
  // 域分解 (Domain Decomposition) 这部分是并行计算中最关键的步骤之一：如何分工。
  //------------------------------------------------------------
  timer.startTimer();
  auto domains = conf->getNumDomains(); // 读取用户希望如何切分几何空间
  mpi::setDomainDecomposition(MPI_COMM_WORLD,
                              domains[0],
                              domains[1],
                              domains[2],
                              domains[3]);
  /*
  设置空间域分解。
  原理：假设我们要计算一个巨大的反应堆。单台电脑内存不够，算得也慢。
  这个函数负责把整个几何空间（大蛋糕）切成很多小块（小蛋糕），然后分配给不同的 CPU 核心（进程）去计算。
  MPI_COMM_WORLD 是所有进程的集合。
  domains[0]...[3] 指定了在各个方向上切分的数量
  */
  timer.stopTimer("Domain Decomposition");
#endif

  //------------------------------------------------------------
  // Initialize the MaterialHandler for later use - 初始化 MaterialHandler 以供后续使用
  //------------------------------------------------------------
  log::info("Creating geometry and materials...");

  timer.startTimer();

  Factory::MaterialHandlerPtr mat_input; // 定义一个智能指针，指向材料处理器对象。
  mat_input = Factory::getMaterialHandler<MaterialHandlerHDF5>(conf);
  /*
  再次使用了工厂模式。
  <MaterialHandlerHDF5>: 这里明确指定了要创建的是处理 HDF5 格式的处理器。
  HDF5: 是一种专门用于存储海量科学数据的文件格式。在核工程中，材料的截面数据非常庞大，通常都存储在 HDF5 文件中。
  这行代码的作用是：创建一个能读取 HDF5 文件的对象，并把配置信息 conf 传给它（告诉它去读哪个文件）
  */

  // Set the number of energy group for debugging - 设置用于调试的能群数量
  { //: 这是一对大括号构成的一个局部作用域

    auto n = conf->getRequiredEnergyGroups(); // 从配置中获取用户指定的能群数量
    if (n > 0)
      mat_input->setNumEnergyGroups(n);
  }
  // mat_input->readAllMaterials();

  /* Create Non-uniform CMFD mesh - 创建非均匀 CMFD 网格 */
  /*
  log::info("Creating CMFD mesh...");
  int axial_refines = 15;
  Cmfd* cmfd = new Cmfd();
  cmfd->setSORRelaxationFactor(1.5);
  cmfd->setCMFDRelaxationFactor(0.7);
  std::vector<std::vector<double>> widths = {
    {1.26,1.26,1.26},
    {1.26,1.26,1.26},
    {10.71,10.71,10.71,10.71}
  };
  cmfd->setWidths(widths);
  */

  /* Create CMFD mesh 创建CMFD网格 */
  log::fdebug("Creating CMFD mesh...");
  Cmfd *cmfd = new Cmfd();
  cmfd->setSORRelaxationFactor(1.5);   // 设置求解用的SOR迭代松弛因子
  cmfd->setLatticeStructure(3, 3, 30); // 设置CMFD网格的逻辑结构  每一个cmfd网格都是一个lattice

  if (!conf->isHexTallyMesh()) // 判断不是六边形，默认就是四边形
  {                            // 设置MOC能群和CMFD能群的对应关系
    std::vector<std::vector<int>> cmfd_group_structure;
    cmfd_group_structure.resize(2); // cmfd能群压缩，MOC是7个能群，压缩为2个能群，1-3能群压缩为第1个能群，4-7能群压缩第为2个能群
    for (int g = 0; g < 3; g++)
      cmfd_group_structure.at(0).push_back(g + 1); // 这里g+1是为了便于统计总共有多少个MOC能群数，以及便于后续判断CMFD对应的MOC能群数是否是线性递增的
    for (int g = 3; g < 7; g++)
      cmfd_group_structure.at(1).push_back(g + 1);
    cmfd->setGroupStructure(cmfd_group_structure);
  }
  cmfd->setKNearest(3); // 设置更新MOC通量所关联的径向CMFD网格数量，后续会根据这个数量来更新MOC通量

  // set for HexLattice 000000000 - 为六边形设置 000000000
  cmfd->setHexLatticeEnable(conf->isHexTallyMesh());
  // CMFD Lattice aligns with mesh
  if (conf->isHexTallyMesh())
  {
    log::fdebug("Create CMFD Lattice with Hex ...");
    cmfd->setHexGroups(10);
    auto widths = conf->getTallyMesh();
    cmfd->setOrientation(conf->getTallyMeshOrientation());
    cmfd->setNumR(widths[0][0]);
    cmfd->setWidthR(widths[1][0]);
    cmfd->setNumZ(10);
    // cmfd->setWidthsZ(widths[2]);
  }

  //------------------------------------------------------------
  // Create the geometry - 创建几何
  //------------------------------------------------------------
  auto nmods = conf->getNumModules();
  Geometry geometry;
  geometry.setNumDomainModules(nmods[0], nmods[1], nmods[2]);

  auto geo_input =
      Factory::getGeoInput<GeoInputXml>(&geometry, mat_input, conf);

  // Try to read primitives - 尝试读取基本几何原语
  geo_input->readGlobalPrimitives(conf->getGlobalPrimitivesPath());
  // Read the geometry - 读取几何
  geo_input->readGeometryFromFile(conf->getGeoInputPath());

  // Dump geometry settings as needed - 按需导出几何设置
  if (conf->doesDumpSettings())
  {
    geo_input->dumpSettings(conf->getOutputDirectory());
  }

  timer.stopTimer("Read Geometry");

  // Remove unused CSG objects and delete temporary containers - 移除未使用的 CSG 对象并删除临时容器
  timer.startTimer();
  if (conf->doesCleanupPrimitives())
    geo_input->clear();
  timer.stopTimer("Clean Up");

#ifdef ENABLE_MPI_
  timer.startTimer();
  if (mpi::isSpatialDecomposed())
  {
    geometry.setDomainDecomposition();
  }
  timer.stopTimer("Domain Decomposition");
#endif

  // set CMFD for - 为几何设置 CMFD
  geometry.setCmfd(cmfd); // 将cmfd指针传递给geometry  关联几何

  // Initialize FSR - 初始化 FSR
  timer.startTimer();
  geometry.initializeFlatSourceRegions(); // 涉及CMFD的初始化
  timer.stopTimer("Meshing");

  // Remove unused materials - 移除未使用的材料
  timer.startTimer();
  if (conf->doesCleanupMaterials())
    geo_input->eraseUnusedMaterials();
  timer.stopTimer("Clean Up");

  // Print the number of materials and delete the pointer - 打印材料数量并删除指针
  geo_input->printReport();
  geo_input.reset();

  // Print the memory usage of materials - 打印材料的内存占用
  geometry.printMemUsageMaterials();

  // FIXME: HDF5 interfaces must be called before mpi finalizing - FIXME：在 MPI 结束前必须调用 HDF5 接口
  timer.startTimer();
  mat_input.reset();
  timer.stopTimer("Clean Up");

  timer.stopTimer("Input");
  printTimerReport(timer);

  //------------------------------------------------------------
  // Generate tracks - 生成轨迹
  //------------------------------------------------------------
  log::info("Initializing the track generator...");
  auto quad = Factory::getQuadrature(conf);
  auto tg = Factory::getTrackGenerator(&geometry, conf);
  tg->setQuadrature(quad);
  tg->generateTracks(); // 轨迹生成，涉及粗网格与FSR的映射，以及轨迹段与粗网格的映射关系

  //------------------------------------------------------------
  // Run simulation - 运行仿真
  //------------------------------------------------------------
  log::info("Running simulation...");
  SolverPtr solver = Factory::getSolver(tg, conf);

  solver->computeEigenvalue(conf->getMaxIterations()); // 整个计算过程，涉及MOC迭代，CMFD迭代，以及通量更新
  solver->printTimerReport();

  //------------------------------------------------------------
  // Output reaction rates, fluxes and volumes - 输出反应率、通量和体积
  //------------------------------------------------------------
  if (conf->doesDumpVisualizationData())
  {
    auto mesh_simple = createSimpleMesh(conf, solver);
    dumpMeshData(conf, mesh_simple);
  }

  // Dump output data to HDF5 files - 将输出数据写入 HDF5 文件
  if (conf->doesDumpFSRData() && mpi::isInMainCommUniqueDomains())
    dumpDataToH5File(conf, solver);

  log::header("Finished");

#ifdef ENABLE_MPI_
  MPI_Finalize();
#endif
  return 0;
}

//----------------------------------------------------------------------
// Helpers - 辅助函数
//----------------------------------------------------------------------
/// \brief Print a timing report for pre-processing - 预处理阶段计时报告
void printTimerReport(Timer &timer)
{

#ifdef ENABLE_MPI_
  timer.reduceTimer(mpi::getMPIComm());
#endif

  log::header("Timing Report - Input (average)");

  const std::string tot_string = "Input";
  timer.printSplit(tot_string, "Total Input Time");

  timer.printSplit("Read Settings", "Time to read settings",
                   1, tot_string);

  timer.printSplit("Read Geometry", "Time to read geometry and materials",
                   1, tot_string);

#ifdef ENABLE_MPI_
  timer.printSplit("Domain Decomposition", "Time of domain decomposition",
                   1, tot_string);
#endif

  timer.printSplit("Meshing", "Time of meshing and FSRs initialization",
                   1, tot_string);

  timer.printSplit("Clean Up", "Time to remove unused objects",
                   1, tot_string);

  log::separator("-");
}

/// \brief Create a Mesh object with the shape of a simple lattice. - 创建一个简单lattice形状的 Mesh 对象
/// \param conf A pointer to ConfInput object. - 参数 conf：指向 ConfInput 对象的指针
/// \param solver A pointer to Solver object. - 参数 solver：指向 Solver 对象的指针
/// \return A pointer to the new Mesh object. - 返回：指向新建 Mesh 对象的指针
std::shared_ptr<Mesh> createSimpleMesh(Factory::ConfInputPtr conf, Factory::SolverPtr solver)
{

  log::info("Creating simple tally mesh...");

  // Define the mesh object - 定义 mesh 对象
  std::shared_ptr<Mesh> mesh = nullptr;

  // Retrieve the shape of the lattice - 获取lattice形状
  auto widths = conf->getTallyMesh();

  // Initialize the mesh object - 初始化 mesh 对象
  mesh = std::make_shared<Mesh>(solver);

  if (conf->isHexTallyMesh())
  { // Hexagonal 六边形
    int num_r = (int)widths[0][0];
    double width_r = widths[1][0];
    mesh->createLattice(num_r, width_r, widths[2],
                        conf->getTallyMeshOrientation(),
                        conf->getTallyMeshOffset());
  }
  else
  { // Rectangular 矩形
    mesh->createLattice(widths[0], widths[1], widths[2]);
  }

  return mesh;
}

/// \brief Dump visualization data to files. - 将可视化数据导出到文件
/// \details This method dumps rates, cross-sections, volumes and - 该方法按需导出反应率、截面、体积以及
///          tracks data on demand. - 轨迹数据
/// \param conf A pointer to ConfInput object. - 参数 conf：指向 ConfInput 对象的指针
/// \param mesh A mesh object to dump visualization data. - 参数 mesh：用于导出可视化数据的 Mesh 对象
void dumpMeshData(Factory::ConfInputPtr conf, std::shared_ptr<Mesh> mesh)
{

  // Create the directory if it doesn't exist - 如目录不存在则创建
  std::string dump_dir = conf->getOutputDirectory();
  fileutils::createDirectory(dump_dir);

  // Cross-sections are ignored because I don't not how to sum them up across FSRs - 忽略截面，因为目前尚不清楚如何在多个 FSR 上求和
  auto field_types = conf->getDumpRXTypes();
  tallyutils::expandTallyTypeForRX(field_types);

  // Selected energy groups - 选定的能群
  auto energy_groups = conf->getTallyEnergyGroups();
  mesh->dumpMeshDataToFile(dump_dir + "/reaction_rates.vtu", field_types, energy_groups);

  // Dump 2D or 3D tracks - 导出二维或三维轨迹
  // FIXME: need to be refactored to take advantage of PHDF5 - FIXME：需要重构以利用 PHDF5
  // mesh.dumpTrackMeshToDir(tg, conf->getDumpTracksTypes(), dump_dir);
}

/// \brief Dump visualization data to HDF5 files. - 将可视化数据导出到 HDF5 文件
/// \details This method dumps rates, cross-sections, volumes and - 该方法按需导出反应率、截面、体积以及
///          tracks data on demand. - 轨迹数据
/// \param conf a pointer to ConfInput object - 参数 conf：指向 ConfInput 对象的指针
/// \param solver a pointer to Solver object - 参数 solver：指向 Solver 对象的指针
void dumpDataToH5File(Factory::ConfInputPtr conf, Factory::SolverPtr solver)
{

  // Create the directory if it doesn't exist - 如目录不存在则创建
  std::string dump_dir = conf->getOutputDirectory();
  std::string file = dump_dir + "/fsr_data.h5";
  fileutils::createDirectory(dump_dir);

  // Selected field types - 选定的场类型
  auto field_types = tallyutils::combineRXAndXS(conf->getDumpRXTypes(),
                                                conf->getDumpXSTypes());

  // Create a file and write data to it - 创建文件并写入数据
  FSRDataHandlerHDF5 h5_writer(file, HDF5Mode::Truncate, solver);

  // Selected energy groups - 选定的能群
  auto energy_groups = conf->getTallyEnergyGroups();

  h5_writer.writeMetaData();
  h5_writer.dumpFSRData(field_types, energy_groups);
}
