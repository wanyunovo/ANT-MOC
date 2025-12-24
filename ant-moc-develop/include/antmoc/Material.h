/**
 * @file Material.h
 * @brief
 * @date January 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "antmoc/constants.h"

#include <string>

namespace antmoc
{

#ifdef INTEL
/** Aligned memory allocation for Intel's compiler */
#define MM_MALLOC(size, alignment) _mm_malloc(size, alignment)

/** Aligned memory deallocation for Intel's compiler */
#define MM_FREE(array) _mm_free(array)

#else
/** Aligned memory allocation for GNU's compiler */
#define MM_MALLOC(size, alignment) memalign(alignment, size)

/** Aligned memory deallocation for GNU's compiler */
#define MM_FREE(array) free(array)
#endif

  int material_id();
  void reset_material_id();
  void maximize_material_id(int material_id);

  /**
   * @class Material Material.h "src/Material.h"
   * @brief The Material class represents a unique material and its relevant
   *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
   */
  class Material
  {

  private:
    /** A user-defined ID for each Material created */
    int _id;

    /** A user-defined name for the Material */
    char *_name;

    /** The volume / area of the Material computed from overlapping segments */
    double _volume;

    /** The total number of instances of this Material in the Geometry  实例数*/
    int _num_instances;

    /** The number of energy groups */
    int _num_groups;

    /** An array of the total cross-sections for each energy group */
    FP_PRECISION *_sigma_t;

    /** Max total cross section */
    FP_PRECISION _max_sigma_t;

    /** A 2D array of the scattering cross-section matrix from/into each group
     * 散射矩阵（概念图）：矩阵元素 σ(gp→g) 表示：中子从能群 gp 散射到能群 g
     * 标准约定：核物理领域习惯列=源能群，行=目标能群
             源能群 (gp)
          1    2    3    4
目标 1  [σ11  σ21  σ31  σ41]  从能群1/2/3/4散射到能群1 计算σ_a[0]时，减去这一行
能群 2  [σ12  σ22  σ32  σ42]  从能群1/2/3/4散射到能群2 计算σ_a[1]时，减去这一行
(g)  3  [σ13  σ23  σ33  σ43]  从能群1/2/3/4散射到能群3
     4  [σ14  σ24  σ34  σ44]  从能群1/2/3/4散射到能群4

     虽然逻辑上是二维矩阵，但在C++中用一维数组存储，需要索引转换。
     1D数组索引 = 列索引 × 行数 + 行索引
          = gp × _num_groups + g

逻辑上的二维矩阵：             实际的一维数组存储（列主序）：
     gp=0  gp=1  gp=2  gp=3
g=0 [ 0    4     8    12 ]    索引0:  (gp=0, g=0)
g=1 [ 1    5     9    13 ]    索引1:  (gp=0, g=1)
g=2 [ 2    6    10    14 ]    索引2:  (gp=0, g=2)
g=3 [ 3    7    11    15 ]    索引3:  (gp=0, g=3)
                              索引4:  (gp=1, g=0)
_sigma_s[] = {0,1,2,3, 4,5,6,7, 8,9,10,11, 12,13,14,15}
             └──列0──┘└──列1──┘└───列2───┘└────列3───┘
     */
    FP_PRECISION *_sigma_s;

    /** An array of the absorption cross-sections for each energy group */
    FP_PRECISION *_sigma_a;

    /** An array of the fission cross-sections for each energy group */
    FP_PRECISION *_sigma_f;

    /** An array of the fission cross-sections multiplied by 单次裂变发射的平均中子数 nu \f$ \nu \f$
     *  for each energy group */
    FP_PRECISION *_nu_sigma_f;

    /** An array of the 裂变中子能谱chi \f$ \chi \f$ values for each energy group */
    FP_PRECISION *_chi;

    /** A 2D array of the fission matrix from/into each group */
    FP_PRECISION *_fiss_matrix;

    /** A boolean representing whether or not this Material contains a non-zero
     *  fission cross-section and is fissionable */
    bool _fissionable;

    /** A boolean to indicate whether or not the data has been
     * allocated to be vector aligned for SIMD instructions */
    bool _data_aligned;

    /** The number of vector widths needed to fit all energy groups */
    int _num_vector_groups;

  public:
    Material(int id = 0, const char *name = "");
    virtual ~Material();

    int getId() const;
    char *getName() const;
    double getVolume();
    int getNumInstances();
    int getNumEnergyGroups() const;
    FP_PRECISION *getSigmaT();
    FP_PRECISION getMaxSigmaT();
    FP_PRECISION *getSigmaS();
    FP_PRECISION *getSigmaA();
    FP_PRECISION *getSigmaF();
    FP_PRECISION *getNuSigmaF();
    FP_PRECISION *getChi();
    FP_PRECISION *getFissionMatrix();
    FP_PRECISION getSigmaTByGroup(int group);
    FP_PRECISION getSigmaSByGroup(int origin, int destination);
    FP_PRECISION getSigmaAByGroup(int group);
    FP_PRECISION getSigmaFByGroup(int group);
    FP_PRECISION getNuSigmaFByGroup(int group);
    FP_PRECISION getChiByGroup(int group);
    FP_PRECISION getFissionMatrixByGroup(int origin, int destination);
    bool isFissionable();
    bool isDataAligned();
    int getNumVectorGroups();

    void setName(const char *name);
    void setVolume(double volume);
    void incrementVolume(double volume);
    void setNumInstances(int num_instances);
    void incrementNumInstances();
    void setNumEnergyGroups(const int num_groups);

    void setSigmaT(double *xs, int num_groups);
    void setSigmaS(double *xs, int num_groups);
    void setSigmaF(double *xs, int num_groups);
    void setNuSigmaF(double *xs, int num_groups);
    void setSigmaA(double *xs, int num_groups);
    void setChi(double *xs, int num_groups);

    void setSigmaTByGroup(double xs, int group);
    void setSigmaFByGroup(double xs, int group);
    void setNuSigmaFByGroup(double xs, int group);
    void setSigmaSByGroup(double xs, int origin, int destination);
    void setChiByGroup(double xs, int group);
    void setSigmaAByGroup(double xs, int group);

    void buildFissionMatrix();
    void transposeProductionMatrices();
    void alignData(); // 没懂，后面看
    Material *clone();

    static std::int32_t convertStringToId(std::string id_string);
    static std::string convertIdToString(std::int32_t id);

    std::string toString();
    void printString();
  };

} /* namespace antmoc */

#endif /* MATERIAL_H_ */
