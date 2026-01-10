有限差分法（Finite Difference Method，FDM）是一种求解微分方程数值解的近似方法，其主要原理是对微分方程中的微分项进行直接差分近似，从而将微分方程转化为代数方程组求解

差分就是计算函数导数的离散的近似值。${\left( \frac{\partial u}{\partial x}\right) }_{j} = \mathop{\lim }\limits_{{{\Delta x} \rightarrow  0}}\frac{u\left( {x + {\Delta x}}\right)  - u\left( x\right) }{\Delta x} \approx  \frac{{u}_{j + 1} - {u}_{j}}{\Delta x}$

利用导数的定义可以得到最简单的差分格式


# 1. CMFD

在实际的反应堆物理问题中，通常需要借助 CMFD（Coarse Mesh Finite Difference，**粗网有限差分**）加速来实现收敛。在 CMFD 加速过程中，会求解一个与细网格 MOC 问题保持一致的粗网格问题。

由于粗网格问题能够快速求解，因此可以充分收敛并用于更新 MOC 解，从而在更少的迭代次数内实现全局行为的传递。本文档将把 CMFD 作为一种多重网格方法进行探讨，展示 CMFD 方程如何从基本的多群输运方程推导而来，然后讨论应用 CMFD 加速的过程。

> CMFD方法最初由K.S. Smith提出，作为一种非线性迭代加速方法，通过在粗网格上建立近似扩散方程，快速逼近细网格解的全局分布，从而加速特征线法（MOC）迭代过程。

## 1.1 多重网格方法——CMFD 的理论基础

多重网格方法在数值分析中常用于求解微分方程。其基本思想是在粗网格上传递全局信息比在细网格上快得多。基于这一原理，多重网格方法会交替在粗网格和细网格上求解方程组。在粗网格上，问题规模减小，信息传播速度更快；在细网格上，离散化能精确捕捉问题的解。粗网格解与细网格解保持一致性至关重要。在此背景下，一致性指的是在收敛时，粗网格解和细网格解在粗网格上是一致的。

多重网格方法可以有多种不同的结构，但它们通常包含两个重要阶段：限制和延拓。

1. 限制/粗化（Restriction）：将细网格MOC的解压缩为粗网格（CMFD）的一致形式。

2. 延拓（Prolongation）：利用粗网格的解插值修正细网格的解。将收敛后的CMFD粗网格标通量，反推更新MOC细网格通量，为下一轮MOC迭代提供更优的项。

   > 粗化一般指将细网的计算结果传输至粗网中，而延拓则是指通过粗网格求得的解或来更新细网。在每一层中，粗化和延拓将会分别将网格信息转向更粗和更细的网格

一般来说，多重网格方法可能涉及多层网格。在每一层，分别使用限制和延拓来转移到较粗和较细的网格层。OpenMOC 中仅用两层网格——细网格（MOC）和粗网格（CMFD）。下图展示了使用CMFD加速求解 MOC 方程的过程。

<img src="img\image-20251112093425490.png" alt="image-20251112093425490" style="zoom:50%;" />

> - **轴向压缩**：特指Z 轴方向的粗化，即把 MOC 中沿轴向分布的多个连续 “平源区（FSR）” 合并为单个 CMFD 粗单元，本质是空间尺度的 “细→粗” 压缩，确保轴向通量分布的全局特征被保留。
> - **粗化**：广义上包括**空间粗化（含轴向压缩）** 和**能群粗化**两类操作，均通过 “加权平均” 或 “求和归并” 将 MOC 的细尺度物理量（FSR 级、多能群）转化为 CMFD 的粗尺度物理量（单元级、少能群）

<img src="img\image-20251112173813197.png" alt="image-20251112173813197" style="zoom:50%;" />

**CMFD与常规多重网格方法的差异**：常规多重网格在各网格层使用相同形式的方程，而CMFD求解粗网格问题采用一种**类扩散方程**（基于体积平均量，与角方向无关），而非中子输运方程的相同MOC形式（依赖角通量和轨迹离散化）。

总而言之，CMFD在将细网格粗化时，也就是（Restriction）限制过程，其中有个要求就是要保证粗网格求解出来的解，需要和细网格的解保持一致。而CMFD方程是使用解一致但本质上不同的方程来求解的。

| 符号       | 物理含义                                                     | 具体说明                                                     |
| ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| $r$        | 空间位置向量 (Spatial position vector)                       | 三维矢量（如直角坐标系*x*,*y*,*z*），描述中子在反应堆、屏蔽体等介质中的具体位置。 |
| $\Omega$   | 角方向向量 (Angular direction vector)                        | ==描述中子的运动方向的单位矢量，由 “极角*θ*” 和 “方位角*ϕ*” 确定（球面坐标系）== |
| $E$        | 中子能量 (Neutron energy)                                    | 中子的动能（单位通常为 eV 或 MeV），不同能量的中子与原子核(靶核)发生反应（散射、吸收）的概率差异极大（由 “截面” 决定）。 |
| $\Psi$     | 中子角通量 (Angular neutron flux)                            | 引入立体角Ω后，中子通量对应的中子角通量密度（单位时间，单位立体角内通过单位面积的中子数） |
| $k_{eff}$  | 有效中子增殖系数 (Effective neutron multiplication factor)   | 对给定系统，新生一代的中子数和产生它的直属上一代中子数之比。实际上从中子的平衡关系来定义系统的有效增殖因数${\mathrm{k}}_{\mathrm{{eff}}} = \frac{\text{ 系统内中子的产生率 }}{\text{ 系统内中子的总消失(吸收 + 泄漏)率 }}$，当有效增殖因数等于1时，中子的产生率恰好等于消失率，系统处于临界状态，这种系统叫做临界系统。 |
| $\Sigma^T$ | 中子总截面 (Neutron total cross-section, A.1.1 核截面)       | 一个中子在介质中前进单位距离发生核反应的几率（散射截面+吸收截面（包括裂变、俘获）） |
| $\Sigma^S$ | 中子散射截面 (Neutron scattering cross-section)              | 一个中子在介质中前进单位距离发生散射反应的几率（能量与运动方向会发生变化，能量 ${E}^{\prime } \rightarrow  E$，运动方向 ${\mathbf{\Omega }}^{\prime } \rightarrow  \mathbf{\Omega }$） |
| $\Sigma^F$ | 中子裂变截面 (Neutron fission cross-section)                 | 一个中子在介质中前进单位距离发生裂变反应的几率               |
| $\chi$     | 裂变中子能谱 (Energy spectrum for fission neutrons)          | 每次裂变，这些中子在不同能量下产生的概率（可以得到每次裂变产生的中子的能量分布，平均能量是2Mev） |
| $\nu$      | **单次裂变发射的平均中子数** (Average number of neutrons emitted per fission) | 每次裂变平均产生的中子数                                     |

## 1.2 由MOC方程推导CMFD方程：

CMFD 方程可以从多群输运方程中推导出来。**其基本概念是将基于角通量的输运方程转化为基于反应堆内某个体积上的平均标通量的类扩散问题**。在这个过程中，会引入一些近似。然而，所有这些近似在收敛时都不会引入偏差。因此，它们不会影响解的精度。

首先引入 MOC 方程中的基于多群近似的公式B.1（==实际上是从MOC文档中的各项同性假设开始推导==）：

> 假定散射各向同性，也就是从任意方向散射出的中子在角度上均匀分布，这样散射源只与**标通量**$\Phi_g(s)$ 有关，不再显式依赖角度。
>
> 我们额外假定散射是各向同性的, 那么源项(2.23)的计算可以简化为只用标通量来表达总源项：$Q = \frac{q}{4\pi }$ 
> $$
> {Q}_{g}\left( s\right)  = \frac{1}{4\pi }\left( {\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}{\sum }_{{g}^{\prime } \rightarrow  g}^{S}\left( s\right) {\Phi }_{{g}^{\prime }}\left( s\right)  + \frac{{\chi }_{g}\left( s\right) }{{k}_{eff}}\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}\nu {\sum }_{{g}^{\prime }}^{F}\left( s\right) {\Phi }_{{g}^{\prime }}\left( s\right) }\right)
> $$
> 方位角和极角下标 $m$ 的 $p$ 从 ${Q}_{g}\left( s\right)$ 中移除了，因为求标通量 ${\Phi }_{g}\left( s\right)$ 要在角度相空间上积分，而这两个角度已经包含在积分中了。
>
> * $q$ : **标量总源项**，所有方向的源加起来(单位体积、单位时间，总的中子产生率)。
> * ${Q}_{i,g}$：**按立体角分配的各向同性源**（单位是“每单位体积、每单位立体角”）
>
> - ${4\pi }$ : 整个空间的立体角(一个球面总固角)。
> - 因为假设源在角度上各向同性(每个方向一模一样)，所以把总源 q 平均分到 4π 个方向上，每个方向就是 $q/{4\pi }$ 。

$$
\boldsymbol {\Omega} \cdot \nabla \psi_ {g} (\mathbf {r}, \boldsymbol {\Omega}) + \Sigma_ {t} ^ {g} (\mathbf {r}) \psi_ {g} (\mathbf {r}, \boldsymbol {\Omega}) = \frac {1}{4 \pi} \left(\frac {\chi_ {g} (\mathbf {r})}{k} \sum_ {g ^ {\prime} = 1} ^ {G} v _ {g ^ {\prime}} (\mathbf {r}) \Sigma_ {f} ^ {g ^ {\prime}} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r}) + \sum_ {g ^ {\prime} = 1} ^ {G} \Sigma_ {s} ^ {g ^ {\prime} \rightarrow g} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r})\right) \tag {B.1}
$$

> ${\Omega} \cdot \nabla \Psi(\mathbf{r}, \Omega, E)$：**中子泄漏率**，描述中子沿方向${\Omega}$的**空间输运**（因中子运动导致空间分布的变化）而离开 “位置*r*、方向Ω、能量*E*” 的速率（即中子 “跑走” 的损失，离开我们的控制体积）。 
>
> * ${\nabla} \Psi$ 表示中子角通量 $\Psi$的**梯度**（即空间各方向的变化率）；
> * ${\Omega} \cdot ({\nabla} \Psi)$ 是梯度与中子运动方向单位矢量 ${\Omega}$ 的**点积**，对应“中子沿运动方向${\Omega}$的**方向导数**”，描述中子因**定向运动**离开dv体积元对应的表面ds的中子数。
>
> $\Sigma^t(\mathbf{r}, E) \Psi(\mathbf{r}, \Omega, E)$：**中子移出率**，$\Sigma^t$为**总宏观截面**（吸收截面（比如俘获） + 散射截面之和），表示中子因“吸收”或“散射”而消失的速率，因为任何中子间的相互作用都会导致中子改变其能量和运动方向或者消失，所以将其简化为总宏观截面。 
>
> 1. 被吸收截面吸收的中子
> 2. 在散射截面发生散射到其他能群或方向的中子
> 3. $\Sigma^t(\mathbf{r}, E) $与r有关：意味着当穿过反应堆时，可能遇到不同的材料，反应堆不同部分的截面可能不同

对方程在整个  $4\pi$  角空间进行积分（角通量在整个角空间上的积分即为标通量），将上述基于角通量  $\psi_{g}(\mathbf{r},\Omega)$  的方程转换为基于标通量的方程  $\phi_g(\mathbf{r})$  得公式B.2：（去掉 $\Omega$ 和 $4\pi$ ）
$$
\int_ {4 \pi} d \boldsymbol {\Omega} \boldsymbol {\Omega} \cdot \nabla \psi_ {g} (\mathbf {r}, \boldsymbol {\Omega}) + \Sigma_ {t} ^ {g} (\mathbf {r}) \phi_ {g} (\mathbf {r}) = \frac {\chi_ {g} (\mathbf {r})}{k} \sum_ {g ^ {\prime} = 1} ^ {G} v _ {g ^ {\prime}} (\mathbf {r}) \Sigma_ {f} ^ {g ^ {\prime}} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r}) + \sum_ {g ^ {\prime} = 1} ^ {G} \Sigma_ {s} ^ {g ^ {\prime} \rightarrow g} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r}) \tag {B.2}
$$

上述公式中只有==中子流项==（streaming term，第一项）与角通量  $\psi_{g}(\mathbf{r},\Omega)$  有关系，由于角度方向向量  $\Omega$  与空间变量  $\mathbf{r}$  无关，因此可以将梯度带入中子流项积分之外，如式B.3所示。

$$
\nabla \cdot \int_ {4 \pi} d \Omega \boldsymbol {\Omega} \psi_ {g} (\mathbf {r}, \boldsymbol {\Omega}) + \Sigma_ {t} ^ {g} (\mathbf {r}) \phi_ {g} (\mathbf {r}) = \frac {\chi_ {g} (\mathbf {r})}{k} \sum_ {g ^ {\prime} = 1} ^ {G} v _ {g ^ {\prime}} (\mathbf {r}) \Sigma_ {f} ^ {g ^ {\prime}} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r}) + \sum_ {g ^ {\prime} = 1} ^ {G} \Sigma_ {s} ^ {g ^ {\prime} \rightarrow g} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r}) \tag {B.3}
$$

然后，将==净中子流项==（the net current） $J_{g}(\mathbf{r})$  定义为：

> 能群 g （角度积分得到的）通过某个表面的中子流密度

$$
J _ {g} (\mathbf {r}) = \int_ {4 \pi} d \Omega \Omega \psi_ {g} (\mathbf {r}, \Omega). \tag {B.4}
$$

将公式B.4带入B.3，**并对任意体积进行积分**(这个体积当然可以是 $\mathbf{V}_j$ )得到公式B.5：

$$
\begin{array}{l} \int_ {V} d \mathbf {r} \nabla \cdot J _ {g} (\mathbf {r}) + \int_ {V} d \mathbf {r} \Sigma_ {t} ^ {g} (\mathbf {r}) \phi_ {g} (\mathbf {r}) = \\\int_ {V} d \mathbf {r} \left(\frac {\chi_ {g} (\mathbf {r})}{k} \sum_ {g ^ {\prime} = 1} ^ {G} v _ {g ^ {\prime}} (\mathbf {r}) \Sigma_ {f} ^ {g ^ {\prime}} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r}) + \sum_ {g ^ {\prime} = 1} ^ {G} \Sigma_ {s} ^ {g ^ {\prime} \rightarrow g} (\mathbf {r}) \phi_ {g ^ {\prime}} (\mathbf {r})\right)  \tag {B.5} \\ \end{array}
$$

然后将整个几何结构划分为很多 CMFD 单元。为每个 CMFD 单元  $j$  定义一个体积  $\mathbf{V}_j$ ，该单元  $j$  由有限个**不重叠**的 MOC 源区域（FSR）组成，==输运方程可以用体积为  $\mathbf{V}_i$  的 MOC 平源区  $i$  的通量和恒定截面来表示==($\mathop{\sum }\limits_{i}V_i = V_j$ )，如式 B.6 所示：（下面这个方程将整个几何结构的输运方程，**体现到一个 CMFD 单元  $\mathbf{V}_j$  中**）

> 常数截面近似：在每个 FSR 内，不仅源项 Q 常数，材料截面(总截面、裂变截面、散射截面)，裂变中子能谱也取区域平均常数。

$$
\begin{array}{l} \int_ {V _ {j}} d \mathbf {r} \nabla \cdot J _ {g} (\mathbf {r}) + \sum_ {i \in j} \int_ {V _ {i}} d \mathbf {r} \Sigma_ {t} ^ {i, g} \phi_ {g} (\mathbf {r}) =  \\\sum_ {i \in j} \int_ {V _ {i}} d \mathbf {r} \left(\frac {\chi_ {i , g}}{k} \sum_ {g ^ {\prime} = 1} ^ {G} v _ {i, g ^ {\prime}} \Sigma_ {f} ^ {i, g ^ {\prime}} \phi_ {g ^ {\prime}} (\mathbf {r}) + \sum_ {g ^ {\prime} = 1} ^ {G} \Sigma_ {s} ^ {i, g ^ {\prime} \rightarrow g} \phi_ {g ^ {\prime}} (\mathbf {r})\right) \tag {B.6} \\ \end{array}
$$

> 需要注意的是， CMFD 的网格线（粗网格单元的边界面）不能把某个 FSR 从中间切开；理想状态是 每个 FSR 整体完全落在某一个 CMFD 单元里（CMFD 边界要么在 FSR 外面，要么与 FSR 边界重合，但不能穿过它内部）。在实际中，若发现一个 FSR 横跨两个或更多 CMFD 单元，那就把这个 FSR 沿着 CMFD 边界再细分成多个更小的 FSR，使得 分割后的每个小 FSR 都只属于一个 CMFD 单元。

由于 MOC 方程通常是求解平源区内的平均标通量  $\overline{\phi_{i,g}}$ ，因此输运方程可以改写为（近似， $\overline{\phi_{i,g}}V _ {i}$ 代替在FSRi处的体积积分）：
$$
\int_ {V _ {j}} d \mathbf {r} \nabla \cdot J _ {g} (\mathbf {r}) + \sum_ {i \in j} \Sigma_ {t} ^ {i, g} \overline {{\phi_ {i , g}}} V _ {i} = \sum_ {i \in j} \left(\frac {\chi_ {i , g}}{k} \sum_ {g ^ {\prime} = 1} ^ {G} v _ {i, g ^ {\prime}} \Sigma_ {f} ^ {i, g ^ {\prime}} \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i} + \sum_ {g ^ {\prime} = 1} ^ {G} \Sigma_ {s} ^ {i, g ^ {\prime} \rightarrow g} \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i}\right) \tag {B.7}
$$

除中子流项中的净中子流外，此方程中给出的所有变量都存在于 MOC 方程中。CMFD单元  $\mathbf{j}$  的边界表面定义为  $\mathbf{S}_{\mathbf{j}}$  ，其表面法向量为  $\mathbf{n}$  。对中子流项应用**高斯散度定理**，可以将其转换为表面积分，如式B.8所示。

> 高斯散度定理：假设有一个光滑的闭合曲面 S 和一个被该曲面围起来的体积 V，则有以下关系：$\oint \mathbf{F} \cdot  \mathbf{n}  \mathrm{d}\mathbf{s} = \int \nabla  \cdot  \mathbf{F}\left( \mathbf{r}\right) \mathrm{d}V.$
>
> CMFD 粗网格单元通常是三维的长方体，它的边界面 ${S}_{j}$ 由 6 个平面组成，因此边界曲面积分可以拆成各个面的面积分之和：
>
> ${\int }_{{S}_{j}}\mathbf{J} \cdot  \mathbf{n}{dS} = \mathop{\sum }\limits_{{f \in  \{ x\pm ,y\pm ,z \pm  \} }}{\int }_{{S}_{j,f}}\mathbf{J} \cdot  {\mathbf{n}}_{f}{dS}$

$$
\int_ {S \in S _ {j}} d S J _ {g} (\mathbf {r}) \cdot \mathbf {n} + \sum_ {i \in j} \Sigma_ {t} ^ {i, g} \overline {{\phi_ {i , g}}} V _ {i} = \sum_ {i \in j} \left(\frac {\chi_ {i , g}}{k} \sum_ {g ^ {\prime} = 1} ^ {G} v _ {i, g ^ {\prime}} \Sigma_ {f} ^ {i, g ^ {\prime}} \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i} + \sum_ {g ^ {\prime} = 1} ^ {G} \Sigma_ {s} ^ {i, g ^ {\prime} \rightarrow g} \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i}\right) \tag {B.8}
$$

为了使 CMFD 问题的计算量更小，对能群结构进行了粗化。通过**能群压缩**，==将一个或多个 MOC 能群合并形成一个特定的 CMFD 能群 e==。为得到这一关系，需对 CMFD 能群 e 内包含的所有 MOC 能群 g 求和，如式 B.9 所示。

$$
\begin{array}{l} \sum_ {g \in e} \left(\int_ {S \in S _ {j}} d S J _ {g} (\mathbf {r}) \cdot \mathbf {n} + \sum_ {i \in j} \Sigma_ {t} ^ {i, g} \overline {{\phi_ {i , g}}} V _ {i}\right) =\\ \tag {B.9}  \sum_ {i \in j} \left(\frac {\sum_ {g \in e} \chi_ {i , g}}{k} \sum_ {g ^ {\prime} = 1} ^ {G} \nu_ {i, g ^ {\prime}} \Sigma_ {f} ^ {i, g ^ {\prime}} \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i} + \sum_ {g ^ {\prime} = 1} ^ {G} \left(\sum_ {g \in e} \Sigma_ {s} ^ {i, g ^ {\prime} \to g}\right) \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i}\right) \\ \end{array}
$$

> 当你对能群 $g \in  e$ 求和时，括号里那一坨只依赖 ${g}^{\prime }$ 不依赖 $g$ ，因此可以提出去：$\mathop{\sum }\limits_{{g \in  e}}\frac{{\chi }_{i,g}}{k}\left( {\mathop{\sum }\limits_{{g}^{\prime }}\nu {\sum }_{f}^{i,{g}^{\prime }}{\bar{\phi }}_{i,{g}^{\prime }}{V}_{i}}\right)  = \frac{\mathop{\sum }\limits_{{g \in  e}}{\chi }_{i,g}}{k}\left( {\mathop{\sum }\limits_{{g}^{\prime }}\nu {\sum }_{f}^{i,{g}^{\prime }}{\bar{\phi }}_{i,{g}^{\prime }}{V}_{i}}\right)$
>
> 散射源对某个“目标能群” $g$ 是：$\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}{\sum }_{s}^{i,{g}^{\prime } \rightarrow  g}{\bar{\phi }}_{i,{g}^{\prime }}{V}_{i}$
>
> 对 $g \in  e$ 求和就是把“散射到粗能群 $e$ 内所有细群的贡献”都加起来：$\mathop{\sum }\limits_{{g \in  e}}\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}{\sum }_{s}^{i,{g}^{\prime } \rightarrow  g}{\bar{\phi }}_{i,{g}^{\prime }}{V}_{i} = \mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}\left( {\mathop{\sum }\limits_{{g \in  e}}{\sum }_{s}^{i,{g}^{\prime } \rightarrow  g}}\right) {\bar{\phi }}_{i,{g}^{\prime }}{V}_{i}$，只有括号里那个与g有关。
>
> 物理意义：通过对 MOC 能群求和，将多能群的中子平衡方程(中子流、吸收、裂变、散射项)归并为 CMFD少能群的平衡方程，压缩能量维度的计算规模。

能群粗化后，CMFD 单元的截面（裂变、散射、总截面等）需通过 “MOC 截面的通量加权平均” 重新计算（粗网格和粗群结构上的 CMFD 截面是根据细网格 MOC 量来定义的），确保粗化后的截面能反映细网格的核反应特性，**CMFD 截面**用下标 C 表示，其定义见公式 B.10-B.13。
$$
\chi_ {C} ^ {j, e} = \frac {\sum_ {i \in j} \left[ \left(\sum_ {g \in e} \chi_ {i , g}\right) \left(\sum_ {g ^ {\prime} = 1} ^ {G} v _ {i , g ^ {\prime}} \Sigma_ {f} ^ {i , g ^ {\prime}} \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i}\right) \right]}{\sum_ {i \in j} \sum_ {g = 1} ^ {G} v _ {i , g} \Sigma_ {f} ^ {i , g} \overline {{\phi_ {i , g}}} V _ {i}} \tag {B.10}
$$

$$
v _ {C} ^ {j, e} \Sigma_ {C, f} ^ {j, e} = \frac {\sum_ {i \in j} \sum_ {g \in e} v _ {i , g} \Sigma_ {f} ^ {i , g} \overline {{\phi_ {i , g}}} V _ {i}}{\sum_ {i \in j} \sum_ {g \in e} \overline {{\phi_ {i , g}}} V _ {i}} \tag {B.11}
$$

$$
\Sigma_ {C, s} ^ {j, e ^ {\prime} \rightarrow e} = \frac {\sum_ {i \in j} \sum_ {g ^ {\prime} \in e ^ {\prime}} \left(\sum_ {g \in e} \Sigma_ {s} ^ {i , g ^ {\prime} \rightarrow g}\right) \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i}}{\sum_ {i \in j} \sum_ {g \in e} \overline {{\phi_ {i , g}}} V _ {i}} \tag {B.12}
$$

$$
\Sigma_ {C, t} ^ {j, e} = \frac {\sum_ {i \in j} \sum_ {g \in e} \Sigma_ {t} ^ {i , g} \overline {{\phi_ {i , g}}} V _ {i}}{\sum_ {i \in j} \sum_ {g \in e} \overline {{\phi_ {i , g}}} V _ {i}} \tag {B.13}
$$

> 以裂变截面为例：CMFD 裂变截面归并(式 B.11)： ${\nu }_{C}^{j,e}\mathop{\sum }\limits_{{C,f}}^{{j,e}} = \frac{\mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{g \in  e}}{\nu }_{i,g}\mathop{\sum }\limits_{f}^{{i,g}}{\bar{\phi }}_{i,g}{V}_{i}}{\mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{g \in  e}}{\bar{\phi }}_{i,g}{V}_{i}}$ (B.11)
>
> 符号含义:
>
> - ${\nu }_{C}^{j,e}\mathop{\sum }\limits_{{C,f}}^{{j,e}}$ ：CMFD 单元 $j$ 、能群 $e$ 的 "有效裂变截面"(含每次裂变发射中子数 $\nu$ )；
> - ${\nu }_{i,g}\mathop{\sum }\limits_{f}^{{i,g}}$ ：FSR $i$ 、能群 $g$ 的裂变截面与中子发射数乘积；
> - 直观来看，求哪一项，就把分子的其余项除掉，**所以分母也是所谓的权重。**
> - 物理意义：截面归并以 MOC 通量为权重，确保 CMFD 的核反应率(如裂变率、吸收率)与 MOC 细网格的反应率一致，满足 “解一致性” 要求(PDF 1.1 节)。
>   其他截面(如散射截面 $\mathop{\sum }\limits_{{C,s}}^{{j,{e}^{\prime } \rightarrow  e}}$ 、总截面 $\mathop{\sum }\limits_{{C,t}}^{{j,e}}$ )的归并逻辑与式 B.11 一致，分别对应式 B.12、B.13，核心都是“通量加权平均”。

请注意，这些横截面涉及对 MOC 通量的加权。因此，CMFD 解取决于 MOC 标通量。在收敛时，没有近似误差，这使得 CMFD 解与 MOC 解完全一致。CMFD 单元体积  $V_{C}^{j}$  和==单元格平均标通量  $\phi_{C}^{j,e}$== 分别由式B.14和式B.15定义。
$$
V _ {C} ^ {j} = \sum_ {i \in j} V _ {i} \tag {B.14}
$$


$$
\phi_ {C} ^ {j, e} = \frac {\sum_ {i \in j} \sum_ {g \in e} \overline {{\phi_ {i , g}}} V _ {i}}{\sum_ {i \in j} V _ {i}} \tag {B.15}
$$

> CMFD 单元的体积是其包含的所有 FSR 体积之和。
>
> 能群粗化后，==CMFD 单元的平均标通量==通过 “FSR 通量的体积加权平均” 得到，与前面的定义类似。
>
> 符号含义:
>
> - ${\phi }_{C}^{j,e}$ ：第 $j$ 个 CMFD 单元、第 $e$ 个能群的平均标通量；
> - ${\overline{\phi }}_{i,g}$ ：第 $i$ 个FSR、第 $g$ 个 MOC 能群的平均标通量；
> - $g \in  e$ ：第 $g$ 个 MOC 能群属于第 $e$ 个 CMFD 能群(能群粗化)；
> - 分母 ${\sum_ {i \in j} V _ {i}}$ (即式 B.14 第 $j$ 个 CMFD 单元的体积)。
> - 物理意义：某个CMFD单元单位体积的中子通量

==这些定义与 CMFD 截面定义一起，构成了 CMFD 加速中的粗化部分==。需要注意的是，无论是在空间上还是能量上，CMFD 网格都比 MOC 网格粗糙得多。图 B-2 对空间网格进行了比较，左侧为 MOC 网格。所展示的网格是 OpenMOC 最终结果中使用的实际网格尺寸。由于采用了线性源近似，这里的 MOC 网格比典型的 MOC 网格粗糙不少。

<img src="img\image-20251101151953342.png" alt="image-20251101151953342" style="zoom:50%;" />

对于能群压缩，OpenMOC 计算使用了70个能量群，而 CMFD 求解器使用了25个或更少的能群。与 MOC 问题相比，粗网格和粗能量群的组合导致CMFD 问题规模非常小。在细网格上压缩截面后（前面那些定义就是用来粗化截面的），CMFD 输运方程与原始多群输运方程非常相似，如式B.16所示。

> $$
> \begin{array}{l} \sum_ {g \in e} \left(\int_ {S \in S _ {j}} d S J _ {g} (\mathbf {r}) \cdot \mathbf {n} + \sum_ {i \in j} \Sigma_ {t} ^ {i, g} \overline {{\phi_ {i , g}}} V _ {i}\right) =\\ \tag {B.9}  \sum_ {i \in j} \left(\frac {\sum_ {g \in e} \chi_ {i , g}}{k} \sum_ {g ^ {\prime} = 1} ^ {G} \nu_ {i, g ^ {\prime}} \Sigma_ {f} ^ {i, g ^ {\prime}} \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i} + \sum_ {g ^ {\prime} = 1} ^ {G} \left(\sum_ {g \in e} \Sigma_ {s} ^ {i, g ^ {\prime} \to g}\right) \overline {{\phi_ {i , g ^ {\prime}}}} V _ {i}\right) \\ \end{array}
> $$

> 把 B.9 里反复出现的 $\sum {\bar{\phi }}_{i,g}{V}_{i}$ 记成一个符号：$A =  \mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{g \in  e}}{\bar{\phi }}_{i,g}{V}_{i}$，有CMFD 单元体积为${V}_{C}^{j} = \mathop{\sum }\limits_{{i \in  j}}{V}_{i}$
>
> 则 CMFD 单元平均标通量为：${\phi }_{C}^{j,e} = \frac{A}{{V}_{C}^{j}}$
>
> B. 13 的结构本质是：${\sum }_{C,t}^{j,e} = \frac{\mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{g \in  e}}{\sum }_{t}^{i,g}{\bar{\phi }}_{i,g}{V}_{i}}{\mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{g \in  e}}{\bar{\phi }}_{i,g}{V}_{i}} = \frac{\text{ (总反应率) }}{A}$
>
> 两边乘回 $A$ 就得到你要的“代入形式”：$\mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{g \in  e}}{\sum }_{t}^{i,g}{\bar{\phi }}_{i,g}{V}_{i} = {\sum }_{C,t}^{j,e}A = {\sum }_{C,t}^{j,e}{\phi }_{C}^{j,e}{V}_{C}^{j}$
>
> 同理，对B.13是一样的套路，但B.11、B.12略有不同，多了一步变换。
>
> <img src="img\image-20251217214807238.png" alt="image-20251217214807238" style="zoom: 33%;" />
>
> 
>
> 先把“粗群结构”补成一个严格定义：
>
> 能群压缩就是：把一个或多个 MOC 能群合并成一个 CMFD 能群 $e$ 。严格写法是：
>
> - 细群索引集合: $\mathcal{G} = \{ 1,2,\ldots ,G\}$
> - 粗群索引集合: $\mathcal{E} = \{ 1,2,\ldots ,E\}$
> - 每个粗群对应一个细群子集 ${\mathcal{G}}_{e} \subset  \mathcal{G}$
> - 这些子集构成一个划分 (partition)：$\mathcal{G} = \mathop{\bigcup }\limits_{{e = 1}}^{E}{\mathcal{G}}_{e},\;{\mathcal{G}}_{{e}_{1}} \cap  {\mathcal{G}}_{{e}_{2}} = \varnothing \left( {{e}_{1} \neq  {e}_{2}}\right)$
> - 这条“划分”一旦成立，就有一个纯数学恒等式(只是重排求和次序)：$\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}\left( \cdots \right)  = \mathop{\sum }\limits_{{{e}^{\prime } = 1}}^{E}\mathop{\sum }\limits_{{{g}^{\prime } \in  {\mathcal{G}}_{{e}^{\prime }}}}\left( \cdots \right)$
>
> 把 B.9 里的散射项按这个恒等式拆开
>
> B. 9 的散射源(对固定目标粗群 $e$ )结构是：$\mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}\left( {\mathop{\sum }\limits_{{g \in  e}}{\sum }_{s}^{i,{g}^{\prime } \rightarrow  g}}\right) {\bar{\phi }}_{i,{g}^{\prime }}{V}_{i}$
>
> 用上面的恒等式，把 $\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}$ 拆成粗群求和：$\mathop{\sum }\limits_{{{e}^{\prime } = 1}}^{E}\mathop{\sum }\limits_{{i \in  j}}\mathop{\sum }\limits_{{{g}^{\prime } \in  {e}^{\prime }}}\left( {\mathop{\sum }\limits_{{g \in  e}}{\sum }_{s}^{i,{g}^{\prime } \rightarrow  g}}\right) {\bar{\phi }}_{i,{g}^{\prime }}{V}_{i}$
>
> 到这里你就看到两件事：
>
> * $\mathop{\sum }\limits_{{{g}^{\prime } \in  {e}^{\prime }}}$ 是在描述“来自粗群 ${e}^{\prime }$ 的所有细群贡献”
> * $\mathop{\sum }\limits_{{{e}^{\prime } = 1}}^{E}$ 是把所有来源粗群 ${e}^{\prime }$ 全部加起来，等价于原来的 $\mathop{\sum }\limits_{{{g}^{\prime } = 1}}^{G}$

$$
\frac {1}{V _ {C} ^ {j}} \sum_ {g \in e} \left(\int_ {S \in S _ {j}} d S J _ {g} (\mathbf {r}) \cdot \mathbf {n}\right) + \Sigma_ {C, t} ^ {j, e} \phi_ {C} ^ {j, e} = \frac {\chi_ {C} ^ {j , e}}{k} \sum_ {e ^ {\prime} = 1} ^ {E} v _ {C} ^ {j, e ^ {\prime}} \Sigma_ {C, f} ^ {j, e ^ {\prime}} \phi_ {C} ^ {j, e ^ {\prime}} + \sum_ {e ^ {\prime} = 1} ^ {E} \Sigma_ {C, s} ^ {i, e ^ {\prime} \rightarrow e} \phi_ {C} ^ {j, e ^ {\prime}} \quad (B. 1 6)
$$

回到中子流项，CMFD 单元  $\mathbf{j}$  的整个表面  $\mathbf{S}_{\mathbf{j}}$  **被划分为H个表面**，这些表面形成了单元 $\mathbf{j}$  与另一个CMFD单元之间的表面。这使得 CMFD 单元  $\mathbf{j}$  的总净中子流可以根据这些界面上的净中子流之和来定义。如式B.17所示
$$
\frac {1}{V _ {C} ^ {j}} \sum_ {g \in e} \sum_ {h = 1} ^ {H} \left(\int_ {S \in S _ {j, h}} d S J _ {g} (\mathbf {r}) \cdot \mathbf {n}\right) + \Sigma_ {C, t} ^ {j, e} \phi_ {C} ^ {j, e} = \frac {\chi_ {C} ^ {j , e}}{k} \sum_ {e ^ {\prime} = 1} ^ {E} v _ {C} ^ {j, e ^ {\prime}} \Sigma_ {C, f} ^ {j, e ^ {\prime}} \phi_ {C} ^ {j, e ^ {\prime}} + \sum_ {e ^ {\prime} = 1} ^ {E} \Sigma_ {C, s} ^ {i, e ^ {\prime} \rightarrow e} \phi_ {C} ^ {j, e ^ {\prime}} \tag {B.17}
$$

利用式B.4中给出的定义，可以得到式B.18中的关系。

> $$
> J _ {g} (\mathbf {r}) = \int_ {4 \pi} d \Omega \Omega \psi_ {g} (\mathbf {r}, \Omega). \tag {B.4}
> $$
>
> 能群 g （角度积分得到的）通过某个表面的净中子流

$$
\int_ {S _ {j, h}} d S J _ {g} (\mathbf {r}) \cdot \mathbf {n} = \int_ {S \in S _ {j, h}} d S \int_ {4 \pi} d \boldsymbol {\Omega} \psi_ {g} (\mathbf {r}, \boldsymbol {\Omega}) (\boldsymbol {\Omega} \cdot \mathbf {n}) \tag {B.18}
$$

<img src="img\image-20251101152117597.png" alt="image-20251101152117597" style="zoom:50%;" />

> 图B-3：CMFD 表面的示意图，其法向量为  $\mathbf{n}$ ，正被一个横截面积为  $\delta \mathbf{A}t$ 、沿  $\Omega$  方向运动的轨迹穿透。被穿透的表面积为  $\delta \mathbf{A}t / (\Omega \cdot \mathbf{n})$ 。若 $\Omega$  和  $\mathbf{n}$  都是单位矢量， $\Omega \cdot \mathbf{n} = \cos \theta$ 。上图显示了给定 MOC 轨迹和 CMFD 表面之间的几何关系。这表明，轨迹  $t$  在表面  $S_{j,h}$  上穿透的表面积可计算为  $\delta \mathbf{A}t / (\Omega \cdot \mathbf{n})$  。
>
> 而且每条轨迹所代表的空间体积，可表示为“轨迹长度 × 轨迹横截面积 δAt”的乘积。

考虑到这种几何关系，并考虑公式2.16中的角权重，通过界面表面的总净中子流可以使用公式B.19计算。

$$
\int_ {S \in S _ {j, h}} d S J _ {g} (\mathbf {r}) \cdot \mathbf {n} = \sum_ {(t, s) \in S _ {j, h}} w _ {i, t} \psi_ {g} ^ {t, s} \left(s _ {j, h}\right) \tag {B.19}
$$



> 轨迹权重  $\mathbf{w}_{\mathrm{t}}$  计算为轨道截面积  $\delta A_{t}$  与角权重  $\alpha_{t}$  的乘积，如式2.16所示
>$$
> w _ {t} = \delta A _ {t} \alpha_ {t} \tag {2.16}
> $$
> 
> 这里的(t,s)是离散索引对，用来表示“哪些离散轨迹(及其段)穿过界面 ${S}_{j,h}$ ”：
>- $t$ ：离散方向编号(对应 ${\Omega }_{t}$ 和角权重 ${\alpha }_{t}$ )。
> - $s$ : 不是路径长度，而是轨迹编号 (第几条 track) 或轨迹段编号 (某条 track 在某个 cell 里的那一段)。
> - 所以$\mathop{\sum }\limits_{{\left( {t,s}\right)  \in  {S}_{j,h}}}\left( \cdots \right)$意思是：对所有“方向为 $t$ 且该方向下第 $s$ 条(或第 $s$ 段)轨迹”与界面 ${S}_{j,h}$ 有交的那些对象做求和。

与 CMFD 截面类似，这些中子流是由 MOC 计算得出，因此在收敛前只是近似值。将 CMFD 组  $e$  内的所有 MOC 组  $g$  相加，可以得到公式 B.20 中 CMFD 组  $e$  穿过单元  $j$  表面  $h$  的净中子流  $\tilde{J}_{j,h,e}$  。
$$
\tilde {J} _ {j, h, e} = \sum_ {g \in e} \sum_ {(t, s) \in S _ {j, h}} w _ {t} \psi_ {g} ^ {t, s} \left(s _ {j, h}\right) \tag {B.20}
$$

> <img src="img\image-20251112194804666.png" alt="image-20251112194804666" style="zoom:50%;" />
>
>  $\tilde{J}_{j,h,e}$  表示CMFD粗网格在*x*方向上与相邻粗网格之间表面h的净中子流，即该网格面的净流出项（流出为正，流入为负，流出-流入=净流）。
>
> MOC方法通过统计轨迹穿过FSR的角通量，计算角通量的变化量来更新FSR的标通量。

需要注意的是，这些估算依赖于角通量。如附录A所述，由于整个角通量矢量没有显式存储，当在 MOC 输运扫描期间遇到 CMFD 表面时，必须计算沿轨迹的角通量对 CMFD 表面净中子流的贡献。这通常是一个相对简单的操作，不会给传输扫描增加太多工作量。根据计算得到的中子流，==新的输运方程==如式B.21所示。（由B.17-B.20得）
$$
\frac {1}{V _ {C} ^ {j}} \sum_ {h = 1} ^ {H} \tilde {J} _ {j, h, e} + \Sigma_ {C, t} ^ {j, e} \phi_ {C} ^ {j, e} = \frac {\chi_ {C} ^ {j , e}}{k} \sum_ {e ^ {\prime} = 1} ^ {E} v _ {C} ^ {j, e ^ {\prime}} \Sigma_ {C, f} ^ {j, e ^ {\prime}} \phi_ {C} ^ {j, e ^ {\prime}} + \sum_ {e ^ {\prime} = 1} ^ {E} \Sigma_ {C, s} ^ {i, e ^ {\prime} \rightarrow e} \phi_ {C} ^ {j, e ^ {\prime}} \tag {B.21}
$$

通过这种新的中子平衡表示法，可以求解新的标通量。然而，这不是一个特征值问题的形式，由于中子流项不依赖于标通量。因此，==引入了新的项扩散系数，将中子流项与标通量联系起来。==如式B.22所示。

$$
\frac {\tilde {J} _ {j , h , e}}{A _ {j , h}} = - u (j, h) \hat {D} _ {j, e} \left(\phi_ {C} ^ {I (j, h), e} - \phi_ {C} ^ {j, e}\right) - \tilde {D} _ {j, h, e} \left(\phi_ {C} ^ {I (j, h), e} + \phi_ {C} ^ {j, e}\right) \tag {B.22}
$$

函数  $u(j,h)$  为单元  $\mathbf{j}$  上的表面  $\mathbf{h}$  的方向（若表面为正  $\mathbf{x}$ 、 $\mathbf{y}$  或  $\mathbf{z}$  方向表面，则方向值为 +1；若为负  $\mathbf{x}$ 、 $\mathbf{y}$  或  $\mathbf{z}$  方向表面，则方向值为-1）， $A_{j,h}$  为表面  $S_{j,h}$  的面积， $\hat{D}_{j,e}$  为CMFD单元  $\mathbf{j}$，粗群 $\mathbf{e}$ 的扩散系数， $\tilde{D}_{j,h,e}$  为在表面h处的非线性修正扩散系数。

> 涉及  $\hat{D}_{j,e}$  的第一项的灵感来自于扩散理论。扩散理论(Fick 定律)将某个粗能群 $e$ 在沿通量梯度方向上的中子流密度近似为$J = - D\left( \frac{\partial \phi}{\partial X}\right)$（若是一维，就看成沿x方向）。在 CMFD 粗网格里，我们只有单元 j 的平均通量 ${\phi }_{C}^{j,e}$ 与相邻单元的平均通量${\phi }_{C}^{I\left( {j,h}\right) ,e}$ 。为了得到面 $h$ 上的中子流，标准做法是:
>
> 1. 引入面中心通量 ${\phi }_{h}^{e}$ (未知)。
>     假设在每个单元里，从单元中心到界面是线性变化，则单元 $j$ 一侧有$\frac{{J}_{j,h,e}}{{A}_{j,h}} =  - {D}_{j,e}\frac{{\phi }_{h}^{e} - {\phi }_{C}^{j,e}}{\Delta {r}_{h}}.$
> 2. 同一个界面，邻居单元 $I = I\left( {j,h}\right)$ 一侧也能写出(连续性)：$\frac{{J}_{j,h,e}}{{A}_{j,h}} =  - {D}_{I,e}\frac{{\phi }_{C}^{I,e} - {\phi }_{h}^{e}}{\Delta {r}_{h}}.$
> 3. 把 ${\phi }_{h}^{e}$ 消元。两式相等解出 ${\phi }_{h}^{e}$ 并代回，得到$\frac{{J}_{j,h,e}}{{A}_{j,h}} =  - \frac{{D}_{j,e}{D}_{I,e}}{\Delta {r}_{h}\left( {{D}_{j,e} + {D}_{I,e}}\right) }\left( {{\phi }_{C}^{I,e} - {\phi }_{C}^{j,e}}\right) .$
>4. 于是自然得出${\widehat{D}}_{j,h,e} = \frac{{D}_{j,e}{D}_{I\left( {j,h}\right) ,e}}{\Delta {r}_{h}\left( {{D}_{j,e} + {D}_{I\left( {j,h}\right) ,e}}\right) }$
> 
> $\hat{D}_{j,h,e}$ 在CMFD均匀网格假设下使用公式B.23计算
> $$
>  \hat {D} _ {j, h, e} = \frac {D _ {j , e} D _ {I (j , h) , e}}{\Delta \mathbf {r} _ {h} \left(D _ {j , e} + D _ {I (j , h) , e}\right)} \tag {B.23}
> $$
> 
> 其中  $\Delta \mathbf{r}_h$  为CMFD单元质心与界面h之间的距离，==函数  $I(j,h)$  计算单元 j 表面  $\mathbf{h}$   相邻CMFD单元的索引==。由于采用均匀网格假设，单元j的质心与表面h之间的距离与相邻单元  $I(j,h)$  的质心与表面h之间的距离相同。
>
> $D _ {j, e}$在公式B.24中定义。
> 
> $$
> D _ {j, e} = \frac {\sum_ {i \in j} \sum_ {g \in e} \frac {1}{3 \Sigma_ {t} ^ {i , g}} \overline {{\phi_ {i , g}}} V _ {i}}{\sum_ {i \in j} \sum_ {g \in e} \overline {{\phi_ {i , g}}} V _ {i}} \tag {B.24}
> $$
> 
> 

> $u(j,h)$  由式B.25计算，其中1是三维空间中全为1的向量，  $\mathbf{n}_{j,h}$  是表面  $S_{j,h}$  的法向量。对于笛卡尔均匀网格，如果曲面是正的x、y或z曲面，则意义为+1，如果是负的x、y或z曲面，则意义为-1。
>
> $$
> u (j, h) = \frac {\mathbb {1} \cdot \mathbf {n} _ {j , h}}{| \mathbb {1} \cdot \mathbf {n} _ {j , h} |} \tag {B.25}
> $$
> 
>仅靠第一项(纯扩散差分)一般不能重现 MOC 统计到的表面净中子流，所以 B.22 再加一项。 ${\widetilde{D}}_{j,h,e}$ 是用 **"一致性条件"**确定——把 MOC 给出的表面净中子流 ${\widetilde{J}}_{j,h,e}$ 和MOC 给出的粗网格单元平均通量 ${\widetilde{\phi }}_{C}$ 代入 B.22（这两项都是当前MOC迭代时计算出来的），使等式成立，然后对 ${\widetilde{D}}_{j,h,e}$ 做代数求解，就得到 B.26。从 B.22 直接移项就能看出 B.26 的结构。
> 
>在表面h处的非线性修正系数的计算如下：
> 
>
> $$
> \tilde {D} _ {j, h, e} = \frac {- u (j , h) \hat {D} _ {j , h , e} \left(\tilde {\phi} _ {C} ^ {I (j , h) , e} - \tilde {\phi} _ {C} ^ {j , e}\right) - \frac {\tilde {J} _ {j , h , e}}{A _ {j , h}}}{\tilde {\phi} _ {C} ^ {I (j , h) , e} + \tilde {\phi} _ {C} ^ {j , e}} \tag {B.26}
> $$
>
> 其中  $\tilde{\phi}_c^{j,e}$  是由MOC迭代计算得到的CMFD单元格平均标通量，波浪表示该量来自 MOC 计算，且在整个 CMFD 迭代过程中保持不变。

回到平衡方程，并代入净中子流的关系，得到==新平衡方程==，如公式B.27所示。

$$
\begin{array}{l} \frac {1}{V _ {C} ^ {j}} \sum_ {h = 1} ^ {H} A _ {j, h} \left(- u (j, h) \hat {D} _ {j, e} \left(\phi_ {C} ^ {I (j, h), e} - \phi_ {C} ^ {j, e}\right) - \tilde {D} _ {j, h, e} \left(\phi_ {C} ^ {I (j, h), e} + \phi_ {C} ^ {j, e}\right)\right) + \Sigma_ {C, t} ^ {j, e} \phi_ {C} ^ {j, e} = \tag {B.27} \\ \frac {\chi_ {C} ^ {j , e}}{k} \sum_ {e ^ {\prime} = 1} ^ {E} \nu_ {C} ^ {j, e ^ {\prime}} \Sigma_ {C, f} ^ {j, e ^ {\prime}} \phi_ {C} ^ {j, e ^ {\prime}} + \sum_ {e ^ {\prime} = 1} ^ {E} \Sigma_ {C, s} ^ {i, e ^ {\prime} \rightarrow e} \phi_ {C} ^ {j, e ^ {\prime}} \\ \end{array}
$$

## 1.3 求解用于MOC加速的CMFD方程

式B.27中的CMFD平衡方程与 MOC 平衡方程在物理上是等价的系统，二者均可用于求解输运方程。对平衡方程中的各项进行重新排列，并按标通量分组，可得到式B.28所示的关系。

> 从 B.27(对固定的 cell $j$ 、粗群 $e$ )出发：
>
> 第一步：两边同乘 ${V}_{C}^{j}$
> $$
> \mathop{\sum }\limits_{{h = 1}}^{H}{A}_{j,h}\left( \cdots \right)  + {\sum }_{C,t}^{j,e}{V}_{C}^{j}{\phi }_{C}^{j,e} = \frac{{\chi }_{C}^{j,e}}{k}\mathop{\sum }\limits_{{e}^{\prime }}{\nu }_{C}^{j,{e}^{\prime }}{\sum }_{C,f}^{j,{e}^{\prime }}{V}_{C}^{j}{\phi }_{C}^{j,{e}^{\prime }} + \mathop{\sum }\limits_{{e}^{\prime }}{\sum }_{C,s}^{j,{e}^{\prime } \rightarrow  e}{V}_{C}^{j}{\phi }_{C}^{j,{e}^{\prime }}.
> $$
> 
>
> 第二步：把中子流项那一坨展开，并把 ${\phi }_{C}^{j,e}$ 与 ${\phi }_{C}^{I\left( {j,h}\right) ,e}$ 的系数拆开
>
> 对每个面 $h$ :
> $$
> - u\widehat{D}\left( {{\phi }_{I} - {\phi }_{j}}\right)  - \widetilde{D}\left( {{\phi }_{I} + {\phi }_{j}}\right)\\= \underset{\text{属于本单元 }j}{\underbrace{\left( {u\widehat{D} - \widetilde{D}}\right) {\phi }_{j}}} + \underset{\text{属于邻单元 }I\left( {j,h}\right) }{\underbrace{\left( {-u\widehat{D} - \widetilde{D}}\right) {\phi }_{I}}}\\= \left( {u\widehat{D} - \widetilde{D}}\right) {\phi }_{j} - \left( {u\widehat{D} + \widetilde{D}}\right) {\phi }_{I}
> $$
> 
>
> 于是左边的项变成：$\mathop{\sum }\limits_{{h = 1}}^{H}{A}_{j,h}\left\lbrack  {\left( {u{\widehat{D}}_{j,e} - {\widetilde{D}}_{j,h,e}}\right) {\phi }_{C}^{j,e} - \left( {u{\widehat{D}}_{j,e} + {\widetilde{D}}_{j,h,e}}\right) {\phi }_{C}^{I\left( {j,h}\right) ,e}}\right\rbrack $
>
> 第三步：把 $\mathop{\sum }\limits_{{e}^{\prime }}{\sum }_{C,s}^{j,{e}^{\prime } \rightarrow  e}{V}_{C}^{j}{\phi }_{C}^{j,{e}^{\prime }}$ 从右边搬到左边，整理得到B.28：

$$
\begin{array}{l} \left[ \Sigma_ {C, t} ^ {j, e} V _ {C} ^ {j} + \sum_ {h = 1} ^ {H} A _ {j, h} (u (j, h) \hat {D} _ {j, e} - \tilde {D} _ {j, h, e}) \right] \phi_ {C} ^ {j, e} - \sum_ {h = 1} ^ {H} A _ {j, h} (\tilde {D} _ {j, h, e} + u (j, h) \hat {D} _ {j, e}) \phi_ {C} ^ {I (j, h), e} \\ - \sum_ {e ^ {\prime} = 1} ^ {E} \Sigma_ {C, s} ^ {i, e ^ {\prime} \rightarrow e} V _ {C} ^ {j} \phi_ {C} ^ {j, e ^ {\prime}} = \frac {\chi_ {C} ^ {j , e}}{k} \sum_ {e ^ {\prime} = 1} ^ {E} v _ {C} ^ {j, e ^ {\prime}} \Sigma_ {C, f} ^ {j, e ^ {\prime}} V _ {C} ^ {j} \phi_ {C} ^ {j, e ^ {\prime}} \tag {B.28} \\ \end{array}
$$

B.28通过定义==标通量向量==  $\Phi_{C}$  写成矩阵特征值问题的形式，==该向量包含所有CMFD单元格平均标通量  $\phi_{c}^{j,e}$==，损失矩阵 A 和公式B.29中的乘法矩阵 M 表示成B.29。

> B.28 本质上是：对每个 CMFD 单元 $j$ 、每个粗能群 $e$ ，都有一条线性方程，未知量包含：
> - 本单元本群: ${\phi }_{C}^{j,e}$
> - 邻单元同群: ${\phi }_{C}^{I\left( {j,h}\right) ,e}$
> - 本单元其它群 (散射耦合): ${\phi }_{C}^{j,{e}^{\prime }}$
>
> 把所有 cell 与能群的未知量按某个固定顺序排成一个长向量：${\Phi }_{C} = {\left\lbrack  \begin{array}{llllll} {\phi }_{C}^{1,1} & \cdots & {\phi }_{C}^{1,E} & {\phi }_{C}^{2,1} & \cdots & {\phi }_{C}^{J,E} \end{array}\right\rbrack  }^{T}$
>
> 接着把 B.28 的左边所有项(泄漏/去除/散射耦合等)统一记成一个损失矩阵 $A$ ，右边的裂变产生项统一记成乘法矩阵 $M$ ，就得到广义特征值问题B.29。

$$
A \Phi_ {C} = \frac {1}{k} M \Phi_ {C} \tag {B.29}
$$

> <img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218213816988.png" alt="image-20251218213816988" style="zoom:33%;" />
>
> A矩阵副对角线上的元素一定为0
>
> <img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218214201669.png" alt="image-20251218214201669" style="zoom:33%;" />



两边左乘 A 逆，可以和一个==常规特征值问题==联系起来

> 特征值问题： 对于 $n$ 阶方阵 $\mathbf{A}$ ，若有数 $\lambda \alpha$ 和 $n$ 维非零列向量$\alpha$，满足${A\alpha} = {\lambda \alpha}$，则称 $\lambda$ 为 $A$ 的特征值，称非零向量 $\alpha$ 为 $A$ 的对应于特征值 $\lambda$ 的特征向量。注意到 ${A\alpha} = {\lambda \alpha}$ ，即 $\left( {A - {\lambda E}}\right)  \alpha   = 0$ ，而 $\left( {A - {\lambda E}}\right) \alpha = 0$ 有非零解的充要条件是$\det \left( {A - {\lambda E}}\right)  = 0$ ，即：
> $$
> f(\lambda)=\det \left( {A - {\lambda E}}\right)=\left| \begin{matrix} {a}_{11} - \lambda & {a}_{12} & \cdots & {a}_{1n} \\  {a}_{21} & {a}_{22} - \lambda & \cdots & {a}_{2n} \\  \vdots & \vdots & & \vdots \\  {a}_{n1} & {a}_{n2} & \cdots & {a}_{nn} - \lambda  \end{matrix}\right|  = 0
> $$
> 上式是关于 $\lambda$ 的一元 $n$ 次方程，称为方阵 $\mathbf{A}$ 的特征方程。其左端 $\left| {\mathbf{A} - \lambda \mathbf{E}}\right|$ 是关于 $\lambda$ 的 $n$ 次多项式，称为方阵 $\mathbf{A}$ 的特征多项式，记为 $f\left( \lambda \right)$ 。在线性代数中，方阵 $\mathbf{A}$ 的特征值求解可通过方阵 $\mathbf{A}$ 的特征多项式的根来得到，但5次以上的代数多项式是没有求根公式的，它们的求解大部分是相当困难的。所以通过特征多项式来求矩阵特征值不太好。在计算机中，通常使用迭代法。

$$
A ^ {- 1} M \Phi_ {C} = k \Phi_ {C}. \tag {B.30}
$$

> 注意：实际计算里不会显式求 ${A}^{-1}$ ，而是每次需要 ${A}^{-1}\left( \cdot \right)$ 时，通过解线性方程组 ${Ax} = b$ 来实现。

任何常用的特征值求解器都可用于求解该方程组。本文采用了简单的**幂迭代法(幂迭代怎么用的？步骤是?幂迭代中为什么会用到红黑SOR迭代)**。这要求每次迭代都求解一个线性方程组。本文使用**红黑SOR迭代法**来求解该线性方程组。

> SOR迭代法是对Gauss-Seidel迭代法的改进，其基本思想是从 ${\mathbf{x}}^{\left( k\right) }$ 出发，先用Gauss-Seidel 方法计算下一步迭代值 ${x}_{i}^{\left( k + 1\right) }$ 。
> $$
> {\tilde{x}_{i}^{\left( k + 1\right) }} = \frac{1}{{a}_{ii}}\left( {{b}_{i} - \mathop{\sum }\limits_{{j = 1}}^{{i - 1}}{a}_{ij}{x}_{j}^{\left( k + 1\right) } - \mathop{\sum }\limits_{{j = i + 1}}^{n}{a}_{ij}{x}_{j}^{\left( k\right) }}\right) \tag {A}
> $$
> 然后取某个参数 $\omega$ ,将当前迭代值 ${\mathbf{x}}^{\left( k\right) }$ 与 Gauss-Seidel 迭代值 ${\widetilde{\mathbf{x}}}_{i}^{\left( k + 1\right) }$ 做加权平均
> $$
> {x}_{i}^{\left( k + 1\right) } = \omega {\tilde{x}_{i}^{\left( k + 1\right) }} + \left( {1 - \omega }\right) {x}_{i}^{\left( k\right) },i = 1,2,\cdots ,n\tag {B}
> $$
> 式(B)即为 SOR 迭代法。由式(A)和式(B)，SOR 迭代法可变形为
> $$
> {x}_{i}^{\left( k + 1\right) } = \left( {1 - \omega }\right) {x}_{i}^{\left( k\right) } + \frac{\omega }{{a}_{ii}}\left( {{b}_{i} - \mathop{\sum }\limits_{{j = 1}}^{{i - 1}}{a}_{ij}{x}_{j}^{\left( k + 1\right) } - \mathop{\sum }\limits_{{j = i + 1}}^{n}{a}_{ij}{x}_{j}^{\left( k\right) }}\right),i = 1,2,\cdots ,n\tag {C}
> $$
> 其中的实参数 $\omega$ 称为松弛因子， $\omega  > 1$ 称为超松弛， $\omega  < 1$ 为低松弛，当 $\omega  = 1$ 时就是 Gauss-Seidel 迭代法。将式(C) 两边同乘 ${a}_{ii}$ ,得
> $$
> {a}_{ii}{x}_{i}^{\left( k + 1\right) } + \omega \mathop{\sum }\limits_{{j = 1}}^{{i - 1}}{a}_{ij}{x}_{j}^{\left( k + 1\right) } = {a}_{ii}\left( {1 - \omega }\right) {x}_{i}^{\left( k\right) } - \omega \mathop{\sum }\limits_{{j = i + 1}}^{n}{a}_{ij}{x}_{j}^{\left( k\right) } + \omega {b}_{i}
> $$
> 上式的向量形式为
> $$
> \left( {\mathbf{D} + \omega \mathbf{L}}\right) {\mathbf{x}}^{\left( k + 1\right) } = \left\lbrack  {\left( {1 - \omega }\right) \mathbf{D} - \omega \mathbf{U}}\right\rbrack  {\mathbf{x}}^{\left( k\right) } + \omega \mathbf{b}\\{\mathbf{x}}^{\left( k + 1\right) } = {\left( \mathbf{D} + \omega \mathbf{L}\right) }^{-1}\left\lbrack  {\left( {1 - \omega }\right) \mathbf{D} - \omega \mathbf{U}}\right\rbrack  {\mathbf{x}}^{\left( k\right) } + \omega {\left( \mathbf{D} + \omega \mathbf{L}\right) }^{-1}\mathbf{b}
> $$
> 简记为
> $$
> {\mathbf{x}}^{\left( k + 1\right) } = {\mathbf{B}}_{\omega }{\mathbf{x}}^{\left( k\right) } + {\mathbf{d}}_{\omega }\tag {D}
> $$
> 式(D)即为 SOR 迭代法的向量形式，其中 ${\mathbf{B}}_{\omega } = {\left( \mathbf{D} + \omega \mathbf{L}\right) }^{-1}\lbrack \left( {1 - \omega }\right) \mathbf{D} -\left. {\omega U}\right\rbrack$ 称为 SOR 迭代矩阵。
>

通常，为确保稳定性，会对修正扩散系数施加一个松弛因子（==这个松弛因子和SOR的松弛因子一样吗？==）。

引入松弛因子  $\omega$  后，在第  $n + 1$  次迭代中，计算得到的修正扩散系数会通过以下公式进行衰减：

$$
\tilde {D} _ {j, h, e} ^ {n + 1} = \omega \tilde {D} _ {j, h, e} ^ {n + 1 / 2} + (1 - \omega) \tilde {D} _ {j, h, e} ^ {n} \tag {B.31}
$$

其中，“半迭代”（half-iterations）指的是如式B.26所示的、未施加衰减时计算得到的扩散系数。松弛因子  $\omega$  可为用户选定的区间[0,1]内的任意实数。松弛因子越小，稳定性越高，但同时也收敛速度也会越慢。

## 1.4 收敛性判别准则

对于当前的标量通量估计 ${\bar{\phi }}_{i,g}$ ，裂变率的相对变化量(记为 $\mathbf{{eps}}$ -MOC)在第 $n$ 次迭代中可按下式计算：

$$
eps-MOC = \frac{1}{{N}_{f}}\sqrt{\mathop{\sum }\limits_{{i \in  {N}_{f}}}{\left( \frac{{f}_{i}^{n} - {f}_{i}^{n - 1}}{{f}_{i}^{n - 1}}\right) }^{2}}
$$
其中， ${N}_{f}$ 是中子产生率非零的区域数量，${f}_{i}^{n}$ 表示第 $n$ 次迭代区域 $i$ 的中子产生率。当式 (7.17) 中定义的 **eps-MOC** 降至  ${10}^{-4}$，并且相对于上一迭代的特征值估计变化小于 **1 pcm** 时，即认为源迭代收敛。采用 **CMFD 加速**时也使用相同的收敛判据。

当应用 CMFD 加速时，还需要为 CMFD 求解器设置容差（数值计算中判断迭代是否收敛的“误差阈值”，当计算误差小于容差时，认为当前步骤收敛）。由于 CMFD 方程采用幂迭代法求解，**且在每次迭代过程中需通过红黑 SOR 算法求解线性方程组**，因此需要两类 CMFD 容差：**幂迭代的容差**以及**线性求解器的容差**。

如果把容差设为常数，那么当系统离收敛还很远时，可能会浪费大量计算时间。在 MOC 的 CMFD 加速背景下，早期的 MOC 源迭代会在中子源上有较大误差，因此即使采用严格的收敛准则，得到的 CMFD 解也会有较大误差。基于这个原因，收敛判据总是相对于**当前 MOC 解的残差**，或相对于**从初始猜测开始的误差下降幅度**来定义的。在每次迭代中，初始猜测取为先前已收敛的解，因此也就隐含地与 **MOC 残差**相关联。

> 残差：当前迭代的解与前一次迭代解之间的偏差，两次迭代之间的相对变化量，若残差过大，说明解仍在显著变化，需继续迭代；若残差小于阈值，则可停止迭代，避免不必要的计算开销。

在本文中，CMFD 幂迭代的容差选为：相对于第一次迭代残差，误差需降低到其 **0.03 倍**（即误差减少  $97\%$ ）。对于线性求解，要求误差降低到 **0.1 倍**（即误差减少  $90\%$ ），或中子产生残差小于 **0.01 eps-MOC**。无论是幂迭代还是线性求解，都强制至少进行 **25 次迭代**。

## 1.5 延拓（Prolongation）阶段

一旦解出CMFD方程，就用该解去**更新 MOC 通量**，从而为下一次迭代产生新的源项。**对MOC通量的更新**就是 CMFD 过程中的“延拓步骤”（prolongation step）。

通量的更新方式有多种。

一种简单的做法是：对所有落在 CMFD 单元 $j$ 、CMFD 能群 $e$ 所覆盖范围内的 MOC 平源区 $i$ 、MOC 能群 $g$的通量，按式 (B.32) 直接缩放：
$$
\phi_ {i, g} ^ {\text {new}} = \phi_ {i, g} ^ {\text {old}} \frac {\phi_ {C , \text {new}} ^ {j , e}}{\phi_ {C , \text {old}} ^ {j , e}} \tag {B.32}
$$

$\phi_{i,g}^{\mathrm{new}}$  为延拓后更新得到的MOC通量，  $\phi_{i,g}^{\mathrm{old}}$  为延拓前的MOC通量，  $\phi_{C,\mathrm{old}}^{j,e'}$  为CMFD求解开始时的**CMFD单元平均通量**(直接由MOC通量计算，这里是某一个单元，不是前面定义的向量)，  $\phi_{C,\mathrm{new}}^{j,e'}$  为CMFD求解收敛后CMFD单元平均通量。

> 需重点注意的是，新、旧标通量必须以相同方式进行归一化处理。

当 CMFD 网格足够细时，上述延拓方式效果很好；但如果已知通量分布形态，还可以通过**对 CMFD 解做插值处理**来进一步增强 CMFD 加速效果。对轴向方向而言，由于几何轴向结构较简单，预期通量分布会平滑变化，因此将每个单元内的通量分布近似为**二次函数**。CMFD 更新前后的通量分布都用相邻域做**二次插值拟合**；所选插值函数需满足保持三个节点（**当前单元及其两个轴向相邻单元**）中的**平均通量**不变。具体而言，能群 e 下 CMFD 单元 j 的轴向通量分布表示为：
$$
\begin{array}{l} \phi_ {C} ^ {j, e} (z) \approx \phi_ {C} ^ {j - 1, e} \left(\frac {\left(z - z _ {j} ^ {C}\right) ^ {2}}{2} - \left(z - z _ {j} ^ {C}\right) - \frac {1}{2 4}\right) + \phi_ {C} ^ {j, e} \left(- \left(z - z _ {j} ^ {C}\right) ^ {2} + \frac {2 6}{2 4}\right) + \tag {B.33} \\ \phi_ {C} ^ {j + 1, e} \left(\frac {\left(z - z _ {j} ^ {C}\right) ^ {2}}{2} - \left(z - z _ {j} ^ {C}\right) - \frac {1}{2 4}\right) \\ \end{array}
$$

式中，z为轴向高度，  $\mathbf{Z}_{\mathrm{j}}^{\mathrm{C}}$  为CMFD单元格j的质心高度，CMFD单元格  j+1  和 j-1 分别代表上下相邻的轴向单元。对于区域边界处的CMFD单元，采用相邻CMFD单元的展开式进行计算。随后，利用上述二次展开式来更新 MOC 通量：

$$
\phi_ {i, g} ^ {\text {new}} = \phi_ {i, g} ^ {\text {old}} \frac {\phi_ {C , \text {new}} ^ {j , e} \left(z _ {i} ^ {C}\right)}{\phi_ {C , \text {old}} ^ {j , e} \left(z _ {i} ^ {C}\right)} \tag {B.34}
$$

$z_{i}^{C}$  是MOC平源区 i 的质心。如果在该区域上只有两个轴向CMFD单元，则用**线性拟合**代替二次拟合。如果区域上只有一个轴向CMFD单元，则不做拟合，直接采用式 (B.32) 的简单关系更新 MOC 通量。CMFD通量轴向拟合示意图见图 B-4。

<img src="IMG\image-20251101154429594.png" alt="image-20251101154429594" style="zoom:50%;" />

> 图 B-4：用于在 CMFD 加速下更新 MOC 通量的轴向延拓（prolongation）示意图。绿色虚线表示在上面两个轴向 CMFD 单元中使用的展开式，橙色虚线表示在下面两个 CMFD 单元中使用的展开式。黑色虚线表示这些展开式的组合结果。

用 CMFD 解来更新 MOC 通量可以显著加快收敛速度，同时快速捕获在粗网格上的通量分布形态。此外，与仅进行一次 MOC输运扫描的计算成本相比，使 CMFD 解完全收敛的计算代价很小。因此，按照本章开头图 B-1 所给出的流程，每次MOC输运扫描后，都会构建并求解CMFD方程。采用 CMFD 加速来求解 MOC 中子输运特征值问题的具体过程在算法 B-1 中给出。

## 1.6 cmfd算法流程

<img src="img\image-20251112174357004.png" alt="image-20251112174357004" style="zoom:50%;" />



<img src="IMG\image-20251101141635108.png" alt="image-20251101141635108" style="zoom:50%;" />

变量说明：

- $\phi$ : MOC 细网格 (FSR×能群) 的标通量向量。
- $\psi$ ：MOC 的角通量向量。
- ${\Phi }_{C}$ ：CMFD 粗网格 (CMFD cell×CMFD 能群) 的标通量向量。它的一个分量才是你看到的 ${\phi }_{C}^{j,e}$ : 表示“第 $j$ 个 CMFD 单元、第 $e$ 个 CMFD 能群”的通量(一个标量)。 
- $A,M$ ：CMFD 的两个稀疏矩阵

### 步骤详解

#### 1. 初始化（Algorithm B-1：第 1–6 行）

1. `procedure COMPUTEEIGENVALUE(geometry, tracks, N)`定义算法主函数，输入参数包括：
   - `geometry`：反应堆几何模型（用于确定网格、材料截面等基础参数）；
   - `tracks`：MOC 计算所需的 “特征线轨迹”（角通量传输的路径，是 MOC 离散角空间的核心）；
2. 从几何模型中隐式定义矩阵 / 算子 $  F, S, H, D, T  $
  * $  F  $：裂变矩阵（描述中子裂变产生新中子的过程）；$  S  $：散射矩阵（描述中子从一个能群 / 方向散射到另一个的过程）；$  H, D, T  $：MOC 输运相关算子（分别对应角度权重、源项分解、轨迹传输，用于角通量*ψ*的计算）。
3. $\phi  \leftarrow  1$ : 标通量初值设为全 1 (常见初值)。
4. $\psi  \leftarrow  0$ ：角通量初值设为 0。（全 0 向量，角通量描述不同方向的中子通量，初始无中子时设为 0）
5. $k \leftarrow  1$ ：特征值 $\left( {k}_{\mathrm{{eff}}}\right)$ 初值设为 1。（表示中子增殖与损失平衡，是反应堆临界状态的初始假设）
6. 归一化: $\phi  \leftarrow  \phi /\left( {{\mathbf{1}}^{T}{F\phi }}\right)$ 。这是用总裂变源把通量规模固定住，避免数值迭代中通量幅值无意义增长。

#### 2.  外层循环：MOC 源迭代（第 7–22 行）

外层` while not converged` 就是在做”用当前 $\phi ,k$ 生成源项 $\rightarrow$ 输运扫描 $\rightarrow$ 得到新 $\phi$ ”的迭代；CMFD 加速被插在每次输运扫描之后。

1. 计算中子源项（第 8 行）：$\mathbf{q} \leftarrow  \frac{1}{4\pi }\left( {{S\phi } + \frac{1}{k}{F\phi }}\right)$

   * $  S{\phi}  $：散射源；
     * $  \frac{1}{k}F{\phi}  $：裂变源。

2. MOC 输运扫描 (第 9-10 行)

   * 角通量 $  {\psi} \leftarrow T^{-1}HD^{-1}{q}  $：本质是“沿特征线解输运方程”得到角通量 $\psi$；
     * 加权量 $  {p} \leftarrow W{\psi}  $：将角通量按角度权重求和，转化为与标量通量相关的量。

   - 第 10 行计算标通量 ${\phi} \leftarrow D^{-1}{q} + D^{-1}{p} $ ：
   - $  D^{-1}{q}  $：未散射中子贡献的通量；$  D^{-1}{p}  $：散射中子贡献的通量，两者之和为总标通量。
     - 此过程中会**统计 CMFD 表面的净中子流**，为后续 CMFD 矩阵构建提供数据。

#### 3. CMFD “粗化” 步骤（第 11–14 行）

这一步把 MOC 的细网格信息压缩成 CMFD 粗网格一致的量。

* 11：构建 CMFD 矩阵 $  A  $和 $  M  $，定义见式 B.28、B.29。

* 12：计算 CMFD 粗网格标通量向量 $  {\Phi}_C  $，定义见式 B.15（本质是 MOC 通量在 CMFD 单元内的体积加权平均）。

* 13：对 $  {\Phi}_C  $ 按**CMFD 总裂变率归一化**（$  {\Phi}_C \leftarrow {\Phi}_C / (\mathbb{1}^T M {\Phi}_C)  $）。

* 14：保存 “CMFD 求解前的粗网格通量” $  \tilde{{\Phi}}_C \leftarrow {\Phi}_C  $，用于后续 “新旧通量比” 的计算。

#### 4. 内层循环：解 CMFD 特征值问题（第 15–19 行）

CMFD 要解的是粗网格广义特征值问题$A \Phi_ {C} = \frac {1}{k} M \Phi_ {C} \tag {B.29}$


  * 16：用线性求解器（高斯 - 塞德尔）计算 $  {\Phi}_C \leftarrow A^{-1}M{\Phi}_C  $；

> 这一步在实现上不是显式求 ${A}^{-1}$ ，而是每次都解线性方程组：${Ax} = M{\Phi }_{C},$ 令 $x = {A}^{-1}M{\Phi }_{C}$


  * 17：更新特征值 $  k \leftarrow \mathbb{1}^T M {\Phi}_C  $；

  * 18：再次归一化 $ {\Phi}_C  $，直到 CMFD 收敛（残差满足文档设定的容差）。

> **CMFD 的两层容差**：因为“外面是幂迭代、里面每步要解线性方程组”，所以需要 **幂迭代容差** + **线性求解容差**。
>
> 文档还给了经验设置：幂迭代误差降到首次残差的 0.03 倍；线性求解误差降到 0.1 倍或产生残差 < 0.01 eps-MOC；并且两者都强制至少 25 次迭代。

#### 5. CMFD “延拓” 步骤（第 20 行）

这一步把 CMFD  粗网格的收敛信息反馈到 MOC 细网格通量上，达到加速 MOC 收敛的目的：

- 20：对所有落在 CMFD 单元 $j$ 、CMFD 能群 $e$ 所覆盖范围内的 MOC 平源区 $i$ 、MOC 能群 $g$的通量，按新旧 CMFD 通量比缩放：${\phi }_{i,g}^{new} = {\phi }_{i,g}^{old}\frac{{\phi }_{C,{new}}^{j,e}}{{\phi }_{C,{old}}^{j,e}}$
- 还可以做轴向插值 ；若只有 1 个轴向 CMFD 单元则不拟合，就按上式缩放。

#### 6. 收敛判断与输出（第 21–23 行）

* 21：再次对 $  {\phi}  $ 进行**总裂变源归一化**，保证物理意义。

* 22：判断外层是否收敛：
  * 裂变率的相对变化量（$  \text{eps-MOC} \leq 10^{-4}  $）

  * 且特征值 $  k  $ 的变化（<1 pcm）。

*  23：输出收敛的 MOC 标通量 $  {\phi}  $ 和特征值 $  k  $。

## 1.7 CMFD面临的挑战

为了最大限度地减少传输解的迭代，CMFD通常基于燃料棒级计算CMFD，在求解方程时，燃料棒单元作为粗网格单元导致CMFD矩阵过大，从而影响CMFD加速效果。

<img src="img\image-20251112204232121.png" alt="image-20251112204232121" style="zoom:33%;" />

# 2. 核心代码

## antmoc.cpp

### void Cmfd::setLatticeStructure

该函数用于设置 CMFD（粗网有限差分）网格的三维结构，即指定 X/Y/Z 三个方向的CMFD网格数，通过`_width_x`/`_num_x`等函数计算单个 CMFD 单元的几何宽度。**整个CMFD网格系统用一个Lattice对象来表示，该Lattice包含3×3×30个单元。**

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251220091108715.png" alt="image-20251220091108715" style="zoom:50%;" />

```c++
void Cmfd::setLatticeStructure(int num_x, int num_y, int num_z) {
    setNumX(num_x);
    setNumY(num_y);
    setNumZ(num_z);
}
void Cmfd::setNumX(int num_x)
{

    if (num_x < 1)
        log::ferror("The number of lattice cells in the x direction "
                    "must be > 0. Input value: %i",
                    num_x);

    _num_x = num_x;        // 整个几何上x方向上CMFD网格的数量
    _local_num_x = _num_x; // 本进程（MPI并行中的单个进程）负责的X方向网格数（默认等于全局）

    /*
    _domain_communicator：MPI通信器对象（用于多进程协调）
    如果启用了MPI并行：
    _num_domains_x：X方向被分成几个子域（进程数）
    将全局网格数均分给各个进程
    例如：100个网格分给4个进程 → 每个进程25个
    */
    if (_domain_communicator != NULL)
        _local_num_x = _num_x / _domain_communicator->_num_domains_x;

    // 计算单个网格宽度
    if (_width_x != 0.)
        _cell_width_x = _width_x / _num_x; // _width_x整个几何上X方向的宽度/整个几何上x方向上CMFD网格的数量，就是一个cmfd网格在x方向上的长度
}
```

###  void Geometry::initializeFlatSourceRegions

这部分先划分 FSR，再设置 CMFD 网格的边界条件、几何尺寸（各个方向的单元数和每个单元的宽度）

<img src="img\image-20251101194715670.png" alt="image-20251101194715670" style="zoom:50%;" />

```c++
  /*
  该方法应在开始源迭代前调用：先通过 Geometry::subdivideCells() 细分所有单元，再初始化 CMFD 对象
   */
  void Geometry::initializeFlatSourceRegions()
  {

    log::finfo("Initializing flat source regions...");

    /* Subdivide Cells into sectors and rings - 将单元细分为扇区与环形 */
    subdivideCells();

    /*  拿到所有cell里用到的材料，存放在_all_materials里，键是材料id，值是材料指针。该方法应在任何需要调用getAllMaterials()之前执行 */
    retrieveAllMaterials();
      
    auto &all_materials = getAllMaterials();

    /* Initialize absorption XS if CMFD absent - 若CMFD模块为空，则初始化材料的吸收截面（Σₐ） */
    if (_cmfd == nullptr)
    {
      /*
      const auto &pair：
      const：不修改元素
      auto：自动推导类型（这里 pair 是键值对）
      &：引用，避免复制
      */
      for (const auto &pair : all_materials) // 遍历：遍历材料映射表中的每个键值对
        pair.second->getSigmaA();            // 初始化时预先计算每个材料在所有能群下的吸收截面
    }

    if (_cmfd != nullptr) /* 若CMFD模块已存在，则初始化CMFD */
      initializeCmfd();
  }
```

####  void Geometry::initializeCmfd()

前面这部分是六边形的初始化，后续补充。

```c++
  void Geometry::initializeCmfd()
  {

	if (_cmfd->GetHexLatticeEnable()){....}// 这部分是HexLattice的CMFD初始化代码
    else
    { // RecLattice

      /* Get the global Geometry boundary conditions
       获取几何模型六个面的边界条件
      */
      boundaryType min_x_bound = _root_universe->getMinXBoundaryType();
      boundaryType max_x_bound = _root_universe->getMaxXBoundaryType();
      boundaryType min_y_bound = _root_universe->getMinYBoundaryType();
      boundaryType max_y_bound = _root_universe->getMaxYBoundaryType();
      boundaryType min_z_bound = _root_universe->getMinZBoundaryType();
      boundaryType max_z_bound = _root_universe->getMaxZBoundaryType();

      /* Get the global Geometry boundaries
       获取几何模型六个边界的坐标值
      */
      double min_x = _root_universe->getMinX();
      double max_x = _root_universe->getMaxX();
      double min_y = _root_universe->getMinY();
      double max_y = _root_universe->getMaxY();
      double min_z = _root_universe->getMinZ();
      double max_z = _root_universe->getMaxZ();
      
           // Set CMFD mesh boundary conditions
      _cmfd->setBoundary(SURFACE_X_MIN, min_x_bound);
      _cmfd->setBoundary(SURFACE_Y_MIN, min_y_bound);
      _cmfd->setBoundary(SURFACE_Z_MIN, min_z_bound);
      _cmfd->setBoundary(SURFACE_X_MAX, max_x_bound);
      _cmfd->setBoundary(SURFACE_Y_MAX, max_y_bound);
      _cmfd->setBoundary(SURFACE_Z_MAX, max_z_bound);

      /* Set CMFD mesh dimensions 设置CMFD网格在三个方向的总宽度，把universe几何体的物理尺寸传递给CMFD */
      _cmfd->setWidthX(max_x - min_x);
      _cmfd->setWidthY(max_y - min_y);
      _cmfd->setWidthZ(max_z - min_z);

      /* Initialize the CMFD lattice 初始化 CMFD lattice */
      Point offset; // 选择几何中心作为参考点
      offset.setX(min_x + (max_x - min_x) / 2.0);
      offset.setY(min_y + (max_y - min_y) / 2.0);
      if (std::abs(min_z + (max_z - min_z) / 2.0) < FLT_INFINITY) // 如果Z中心坐标有效（不是无穷大），使用计算值；否则设为0
        offset.setZ(min_z + (max_z - min_z) / 2.0);
      else
        offset.setZ(0.);

      _cmfd->setGeometry(this);

      _cmfd->initializeLattice(&offset); //  设置表示整个CMFD的lattice的沿各方向的单元宽度以及累积宽度，并设置此Lattice的几何中心在父Universe坐标系中的位置

#ifdef ENABLE_MPI_
      if (isDomainDecomposed()) // 检查是否进行了域分解，域分解: 将计算区域分割给多个进程并行计算
      {
        /*
         mpi::getNumDomainsX(): 获取X方向的域数量（MPI进程数）
         mpi::getDomainIndexX():  获取当前进程负责的子区域在X方向的索引
         */
        _cmfd->setNumDomains(mpi::getNumDomainsX(),
                             mpi::getNumDomainsY(),
                             mpi::getNumDomainsZ());

        _cmfd->setDomainIndexes(mpi::getDomainIndexX(),
                                mpi::getDomainIndexY(),
                                mpi::getDomainIndexZ());
      }
#endif
      /* Initialize CMFD Maps  只是给一个域中的所有CMFD单元对应的FSR数组分配内存*/
      _cmfd->initializeCellMap();
    }
  }
```

####  void Cmfd::initializeLattice

<img src="img\image-20251220224301091.png" alt="image-20251220224301091" style="zoom:50%;" />

同样前面这部分是六边形的初始化。

```c++
/**
   * @brief Initialize the CMFD lattice.
   */
  void Cmfd::initializeLattice(Point *offset)
  {

    if (_hexlattice_enable){...}
    else
    { // Create RecLattice - 创建四边形lattice
      if (_non_uniform)
      { // _non_uniform：布尔标志，表示网格是否非均匀；非均匀：每个网格宽度不同，例如 [2, 1, 3]
        // 从已设置的宽度数组获取网格数量
        setNumX(_cell_widths_x.size());
        setNumY(_cell_widths_y.size());
        setNumZ(_cell_widths_z.size());
      }
      else
      {
        // 均匀：每个方向的单元宽度 = 总宽度 / 单元数
        _cell_width_x = _width_x / _num_x;
        _cell_width_y = _width_y / _num_y;
        _cell_width_z = _width_z / _num_z;

        _cell_widths_x.resize(_num_x, _cell_width_x); // 如果是均匀的CMFD，则将_cell_widths_x数组的所有宽度设为一样的大小
        _cell_widths_y.resize(_num_y, _cell_width_y);
        _cell_widths_z.resize(_num_z, _cell_width_z);
      }

      // 计算三个方向的累计宽度（前缀和），用于定位单元边界与做宽度一致性校验
      _accumulate_x.resize(_num_x + 1, 0.0);
      // 将 _accumulate_x 调整为长度 _num_x + 1 的向量，并把所有元素初始化为 0.0 。num_x 个单元对应 num_x + 1 个边界（含起点和终点）。例如3个格子需要4个边界位置。
      _accumulate_y.resize(_num_y + 1, 0.0);
      _accumulate_z.resize(_num_z + 1, 0.0);

      //_accumulate_x = [0.0, 2.0, 4.0, 6.0...]
      for (int i = 0; i < _num_x; i++)
        _accumulate_x[i + 1] = _accumulate_x[i] + _cell_widths_x[i]; // X 方向累计宽度

      for (int i = 0; i < _num_y; i++)
        _accumulate_y[i + 1] = _accumulate_y[i] + _cell_widths_y[i]; // Y 方向累计宽度

      for (int i = 0; i < _num_z; i++)
        _accumulate_z[i + 1] = _accumulate_z[i] + _cell_widths_z[i]; // Z 方向累计宽度

      // 一致性校验：检查累积宽度是否等于总宽度（允许微小误差）
      if (fabs(_width_x - _accumulate_x[_num_x]) > FLT_EPSILON ||
          fabs(_width_y - _accumulate_y[_num_y]) > FLT_EPSILON ||
          fabs(_width_z - _accumulate_z[_num_z]) > FLT_EPSILON)
      {
        log::ferror("The sum of non-uniform mesh widths are not consistent "
                    "with geometry dimensions. width_x = %20.17E, width_y = %20.17E"
                    ", width_z = %20.17E, sum_x = %20.17E, sum_y = %20.17E, sum_z ="
                    " %20.17E, diff_x = %20.17E, diff_y = %20.17E, diff_z = %20.17E"
                    ", FLT_EPSILON = %20.17E",
                    _width_x, _width_y, _width_z,
                    _accumulate_x[_num_x], _accumulate_y[_num_y],
                    _accumulate_z[_num_z], fabs(_width_x - _accumulate_x[_num_x]),
                    fabs(_width_y - _accumulate_y[_num_y]),
                    fabs(_width_z - _accumulate_z[_num_z]), FLT_EPSILON);
      }

      /* Delete old lattice if it exists - 若存在旧 lattice 则删除 */
      if (_lattice != NULL)
        delete _lattice;

      /* Initialize the lattice */
      _lattice = new RecLattice();
      _lattice->setNumX(_num_x); // 设置 x 方向单元数
      _lattice->setNumY(_num_y); // 设置 y 方向单元数
      _lattice->setNumZ(_num_z); // 设置 z 方向单元数

      if (_non_uniform)
        _lattice->setWidths(_cell_widths_x, _cell_widths_y, _cell_widths_z); // 非均匀：传入每个方向的单元宽度数组
      else
        _lattice->setWidth(_cell_width_x, _cell_width_y, _cell_width_z); // 均匀：传入每个方向的单元宽度即可

      _lattice->setOffset(offset->getX(), offset->getY(), offset->getZ()); // 设置lattice偏移量，offset是整个几何的中心点
      _lattice->computeSizes();                                            // 计算lattice三个方向的累计宽度（前缀和），用于定位单元边界，与前面的_accumulate_x的逻辑一样
    }
  }
```

### void TrackGenerator::generateTracks

这部分生成 MOC 特征线轨迹并进行射线追踪（Ray Tracing，其实算分段），这是一个非常耗时的步骤。

确定每条轨迹分段所在的 CMFD 单元和表面，建立 FSR 与 CMFD 单元的双向映射关系。

- 内部流程：
  1. **initializeQuadrature()**: 确保求积组已正确初始化。
  2. **checkBoundaryConditions()**: 检查几何边界条件（反射、真空、周期性）。
  3. initializeTracks():
     - 直观上来讲是在几何体上按指定的方位角和间距“铺设”特征线（==不确定==）。
     - 初始化轨迹(2D和3D)，包括 ：
       - 分配内存空间。
       - 修正轨迹的角度和间距。
       - 计算轨迹信息,如起点、终点、link 号等。
       - 设置轨迹的边界条件。  
       - 统计轨迹数量，如每角度、每轨迹堆中的轨迹数量。
  4. segmentize() (射线追踪):
     - 让每一条轨迹穿过几何体，计算它穿过了哪些平源区（FSR, Flat Source Regions）。
     - 生成 **轨迹段（Segments）**：记录每一段的长度、所属的材料、所属的 FSR ID。
     - **CMFD 关联**：如果启用了 CMFD，这里还会计算轨迹段属于哪个 CMFD 粗网格，建立 MOC 细网格与 CMFD 粗网格的映射关系（用于后续的加速计算）。CMFD只与轨迹分段相关，与轨迹生成无关，CMFD主要在3D轨迹分段中使用。

<img src="img\image-20251102141250289.png" alt="image-20251102141250289" style="zoom:50%;" />

进入到generateTracks()，在该函数里调用TrackGenerator3D中的segmentize()来进行3D轨迹分段，该函数的逻辑如下:  

* 若不采用 Explicit 3D 策略：使用轴向挤出几何，在重叠平面上分段(segmentizeExtruded)，分段结束。

<img src="img\image-20251210165831629.png" alt="image-20251210165831629" style="zoom:50%;" />

<img src="img\image-20251210170008434.png" alt="image-20251210170008434" style="zoom:50%;" />

#### 重叠平面上的轨迹分段

##### TrackGenerator3D::segmentizeExtruded()

该函数的逻辑为:

* 确定不同轴向层的高度(即 z 坐标,见getUniqueZPlanes)  

* 遍历每条轨迹，执行重叠平面上的轨迹分段(Geometry::segmentizeExtruded(Track *flattened_track,

  ​                  std::vector<double> z_coords, bool FSRonly)) 

<img src="img\image-20251210165430983.png" alt="image-20251210165430983" style="zoom:50%;" />

* 初始化 FSR 的相关数据结构，包括：
  * Geometry::initializeAxialFSRs(std::vector<double> global_z_mesh) ：在已有“2D FSR 和全局/局部轴向网”的基础上，把 3D FSR 的材料、id、网格信息初始化出来。
  * Geometry::initializeFSRVectors()

<img src="img\image-20251102153455189.png" alt="image-20251102153455189" style="zoom:50%;" />

* 统计轨迹上的线段数量

##### Geometry::segmentizeExtruded()

这个函数接受一条 2D 轨迹对象，返回分段后的结果，每一段都存储在该轨迹中。该函数的主要逻辑如下:

* 初始化起点和终点，找到起点所在的Cell。

* 进入轨迹分段的循环，从轨迹起点所在Cell 开始，只要终点 end 仍然落在几何里的某个 Cell 中，就继续分下一段

* 在重叠平面上求最短线段

  * 先初始化线段最小长度(重叠平面上的)，用变量记录局部点所在的 FSR 的 id 和 z_coords索引，这个索引记录最短线段所在的轴向层高度。
  * 每次把起点设置为上一次的终点，然后遍历所有轴向层高度，即对 **每一个可能的 z 高度 z_coords[i]** ，调用findNextCell 尝试找下一次与几何边界的交点，选出**所有 z 里最近的那个交点**，把这一小段当成一个 segment。
  * 到这里为止： 我们已经知道了，从当前起点出发，在所有轴向层里，哪个 z 层的面最先被撞到，以及对应的最短距离 min_length 和层索引 min_z_ind。

* 创建或查找一个 **轴向挤出 FSR（ExtrudedFSR）**，得到 `region_id`。

  ![image-20260102171958579](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20260102171958579.png)

  *  start 是当前轨迹段中、落在这个挤出区里的一个局部坐标点；

  *  先把局部点设定在最短线段所在的轴向层。

  * 遍历局部点的version_num，最大值是个宏MAX_VERSION_NUM。

    *  给局部坐标点（start）设定version_num，尝试创建局部点对应的 FSR 标识串(多线程情况下，可能已经被其他线程创建了)。  
    *  用局部点对应的 FSR 标识串去 map 里找对应 FSR 的特征点(可能就是当前局部点，也可能是其他线程的局部点)。 
    *  遍历所有轴向层，比较这两个点的 CSG 层次信息。如果在某一轴向层上，两者信息不相同，增加version_num，继续循环；如果在每一轴向层上，两者信息都相同，跳出循环。 

  *  如果version_num 已经遍历完，但每一对局部点和特征点的 CSG 层次信息都不完全相同，打印一条错误信息。  

    > 平面上一个 FSR 在重叠平面中可能对应于多个 FSR，因为所有轴向层中的面都会叠加起来。两个点是否位于不同 FSR 中，无法单独在某一个轴向层中判断，只能遍历所有轴向层，在每一层中都比较两个点生成的 FSR 标识串。如果两点确实位于重叠平面上的同一 FSR 中，意味着在任一轴向层生成的 FSR 标识串都相等，也就是说在任一轴向层它们都位于该层的同一 FSR 中；反之，只要有一层不同，就能说明两个点在重叠平面中是分开的。
    > 在version_num 的循环中，每次都会创建新的 FSR 标识串和 id。也就是说，第一次迭代时，如果局部点和 FSR 特征点的 CSG 层次信息完全一致，则跳出循环，此时只创建了一个 FSR 标识串和 id，实例化一个ExtrudedFSR 结构。否则会有多个，且这多个 FSR 的标识串只有version_num  的值不同。

  ![image-20260102172050755](img\image-20260102172050755.png)

*  创建线段，把 `segment._length = 最短长度`、`segment._region_id = 这个 ExtrudedFSR 的 id` 填进去，然后添加线段到2D轨迹中。找下一条线段的终点，继续循环。

* 统计穿过 CMFD 粗网格边界的中子流，以便进行 CMFD 加速。确定当前线段分别穿过了哪个 CMFD 网格的“进”面和“出”面。

> 所以：**一条 2D 轨迹 → 多个 segment，每个 segment 对应一个 ExtrudedFSR。**
>
> 这一步只在 x-y 里分段；后面用 `initializeAxialFSRs` 再沿 z 把每个 ExtrudedFSR 划成一条“3D 柱子”。

```c++
/*
   flattened_track：要分段的这条 2D 轨迹
   z_coords：所有轴向层的 z 高度。
   FSRonly：如果为 true，只建 FSR 不真正往轨迹上加 segment（一般为 false）。
*/ 
void Geometry::segmentizeExtruded(Track *flattened_track,
                                    std::vector<double> z_coords, bool FSRonly)
  {

    /* Track starting Point coordinates and azimuthal angle */
    double x0 = flattened_track->getStart()->getX();
    double y0 = flattened_track->getStart()->getY();
    double z0 = z_coords[0];
    double phi = flattened_track->getPhi();
    double delta_x, delta_y; // 后面处理 CMFD 时需要做一个微小平移（nudging）
    // unused
    // double delta_z;

    log::fdebug("segmentizeExtruded starts at (%f, %f, %f) along "
                "angle %f",
                x0, y0, z0, phi);

    /* Length of each segment */
    double length;
    int min_z_ind; // 对于当前段来说，哪个 z 层给出的交点最近；最短距离对应的z_coords索引
    int region_id; // 找到的 ExtrudedFSR 的全局 id（也就是 FSR id）

    /* Use a LocalCoords for the start and end of each segment */
    LocalCoords start(x0, y0, z0, true);
    LocalCoords end(x0, y0, z0, true);
    start.setUniverse(_root_universe);
    end.setUniverse(_root_universe);

    /* Find the Cell containing the Track starting Point*/
    Cell *curr = findFirstCell(&end, phi); // 起点所在的最底层 Cell，这里的 &end 不是“按引用传参”，而是取对象 end 的地址

    /* If starting Point was outside the bounds of the Geometry */
    if (curr == nullptr)
    {
      int dom = mpi::getDomainUid();
      log::ferror("Could not find a Cell containing the start Point "
                  "of this Track on domain %d with bounds [%f, %f] x [%f, %f]"
                  " x [%f, %f]. Track: %s. Current coords: %s",
                  dom, getMinX(),
                  getMaxX(), getMinY(), getMaxY(), getMinZ(), getMaxZ(),
                  flattened_track->toString().c_str(),
                  end.getPoint()->toString().c_str());
    }

    /* While the end of the segment's LocalCoords is still within the Geometry,
     * move it to the next Cell, create a new segment, and add it to the
     * Geometry */
    int find_cell_count = 0;
    while (curr != nullptr) // 大循环：沿着轨迹一步步往前走,只要终点 end 仍然落在几何里的某个 Cell 中，就继续分下一段
    {

      /* Check if stuck in loop */
      find_cell_count++; // find_cell_count 防止无限循环
      if (find_cell_count > 1e6)
        log::ferror("Caught in infinite loop finding next cell");

      /* Records the minimum length to a 2D intersection */
      double min_length = std::numeric_limits<double>::infinity();
      region_id = -1;
      min_z_ind = -1; 

      end.copyCoords(&start); // 会把 end 的坐标和整条 LocalCoords 链复制给 start（把终点作为起点），从而开启另一段最小线段的寻找

      /* Loop over all z-heights to find shortest 2D intersection */
      for (size_t i = 0; i < z_coords.size(); i++) // 内层循环：对每一个 z 层试一遍交点，选最近的那个
      {

        /* Change z-height and copy starting coordinates to end */
        // 不同z对应的点所在层不同，得到的Cell也可能不同
        start.setZ(z_coords[i]);
        start.prune(); // 只保留这一级 LocalCoords，把下面的链都删掉，避免上一次循环的 CSG 信息干扰现在这一次；

        findCellContainingCoords(&start);
        start.copyCoords(&end); // 把start的值赋给end，重新寻找最小的线段长度

        /* Find the next Cell along the Track's trajectory */
        curr = findNextCell(&end, phi); // 找从当前点出发，沿方向 phi，最近会碰到哪一个几何边界，并把 end 移动到那个交点；同时返回交点所在的 Cell。

        /* Checks that segment does not have the same start and end Points 起点和终点不能相同*/
        if (fabs(start.getX() - end.getX()) < FLT_EPSILON &&
            fabs(start.getY() - end.getY()) < FLT_EPSILON)
          log::ferror("Created segment with same start and end "
                      "point: x = %f, y = %f",
                      start.getX(), start.getY());
        if (end.getX() != end.getX() ||
            end.getY() != end.getY())
        {
          log::ferror("Nan is found in points when loop over z-heights "
                      "during segmenting: x = %f, y = %f, z_index = %d",
                      end.getX(), end.getY(), i);
        }

        /* Find the segment length and extruded FSR */
        length = double(end.getPoint()->distanceToPoint(start.getPoint())); // 计算特定Z轴高度下的最小的长度

        /* Check if the segment length is the smallest found */
        if (length < min_length)
        { // 把最短线段的长度和所在轴向层高度的索引记录下来
          min_length = length;
          min_z_ind = i;
        }
      }
      // 到这里为止： 我们已经知道了，从当前起点出发，在所有轴向层里，哪个 z 层的面最先被撞到，以及对应的最短距离 min_length 和层索引 min_z_ind。

      // 把 start 的 CSG 链清掉；把 z 设置成刚才找到的那一层；再根据这个坐标重新找 Cell，建好 LocalCoords 链。
      start.prune();
      start.setZ(z_coords[min_z_ind]);
      findCellContainingCoords(&start);

      // Create a unique ExtrudedFSR by comparing ExtrudedFSRs with the same CSG
      region_id = createUniqueExtrudedFSR(start, z_coords); // 根据最小线段的Z高度来创建轴向挤出FSR

      /* Move the coordinates to the next intersection */
      start.copyCoords(&end);
      curr = findNextCell(&end, phi, M_PI_2); // 在选定的层上重新走一遍 findNextCell，移动end到终点，并返回其所在的Cell

#ifdef ENABLE_DEBUG_
      log::fdebug("segment start x = %f, y = %f; end x = %f, y = %f",
                  start.getX(), start.getY(), end.getX(), end.getY());
#endif

      /* Check that the region ID is valid */
      if (region_id == -1)
        log::ferror("Failed to find a valid FSR during axial extruded "
                    "segmentation");

      /* Create a new 2D Track segment with extruded region ID  给 segment 填数据*/
      segment new_segment;
      if (!FSRonly)
      { // FSRonly 默认为false
        new_segment._length = min_length;
        new_segment._region_id = region_id; // 这一段所在的 ExtrudedFSR 编号
      }

           /* Save indices of CMFD Mesh surfaces that the Track segment crosses 保存轨迹段所穿过的CMFD网格surface面的索引*/
      if (!FSRonly && _cmfd != nullptr)
      {

        /* Find cmfd cell that segment lies in */
        int cmfd_cell = _cmfd->findCmfdCell(&start); // 计算线段起点所在的CMFD cell，如果起点不在任何 CMFD 网格内（返回 -1），说明该线段不需要贡献给 CMFD 统计，直接使用 goto 跳转到 add_segment 标签，跳过后续 CMFD 处理。

        if (cmfd_cell == -1)
        {
          log::fdebug("segment not in CMFD Cell start(%lf,%lf), end(%lf,%lf)",
                      start.getX(), start.getY(), end.getX(), end.getY());
          goto add_segment;
        }

        if (_cmfd->GetHexLatticeEnable())
        {                                          // 检查是否启用了六角形网格模式。
          if (!_cmfd->CellinHexLattice(cmfd_cell)) // 检查找到的 cmfd_cell 是否属于有效的六角形晶格结构。如果不在有效晶格内，同样跳过。
          {
            log::fdebug("segment not in Hex CMFD Cell start(%lf,%lf), end(%lf,%lf)",
                        start.getX(), start.getY(), end.getX(), end.getY());
            goto add_segment;
          }
        }

        /* Reverse nudge from surface to determine whether segment start or end
         * points lie on a CMFD surface. */
        delta_x = cos(phi) * TINY_MOVE;
        delta_y = sin(phi) * TINY_MOVE;
        start.adjustCoords(-delta_x, -delta_y, 0); // 为了避免起点/终点正好落在 CMFD 网格的顶点或面上，会先往反方向退一点点（TINY_MOVE）； 这有助于判断“是穿过哪一个面”
        end.adjustCoords(-delta_x, -delta_y, 0);

        /* Calculate CMFD surfaces */
        int cmfd_surfaces[2];

        // cmfd_surfaces[0] = _cmfd->findCmfdSurface(cmfd_cell, &end, phi, M_PI_2);
        // cmfd_surfaces[1] = _cmfd->findCmfdSurface(cmfd_cell, &start, phi, M_PI_2);

        cmfd_surfaces[0] = _cmfd->findCmfdSurface(cmfd_cell, &end, phi, M_PI_2);          //  从终点沿正方向 phi 延伸，看它是从哪个 CMFD 面出界；
        cmfd_surfaces[1] = _cmfd->findCmfdSurface(cmfd_cell, &start, phi + M_PI, M_PI_2); // 从起点沿反方向 phi+π 看，是从哪个 CMFD 面进来的。

        /* Ensure surfaces are x-y surfaces (no z-crossings) */
        /* Note: this code takes advantage of the numeric representation of
           surfaces to find a mapping that removes z-surfaces   删除z表面的映射，转为XY平面映射*/
        _cmfd->GetXYSurfaces(cmfd_surfaces);

        /* Save CMFD surfaces */
        // 将找到的“前向面”（出射面）和“后向面”（入射面）索引保存到 new_segment 结构体中。后续求解器在扫描这条线段时，会利用这两个索引将中子流贡献累加到对应的 CMFD 面上。
        new_segment._cmfd_surface_fwd = cmfd_surfaces[0];
        new_segment._cmfd_surface_bwd = cmfd_surfaces[1];

        /* Re-nudge segments from surface. */
        start.adjustCoords(delta_x, delta_y, 0);
        end.adjustCoords(delta_x, delta_y, 0);
      }

    add_segment: // goto标签
      if (!FSRonly)
      {
        /* Add the segment to the 2D track */
        flattened_track->addSegment(new_segment); // 把这段加入到 2D 轨迹的 segments list
        // log::fdebug("segment start(%lf,%lf), end(%lf,%lf), length:%lf, fwd:(%d,%d), bwd:(%d,%d)",
        //             start.getX(), start.getY(), end.getX(), end.getY(),
        //             new_segment._length,new_segment._cmfd_surface_fwd,new_segment._cmfd_surface_fwd % HEX_NUM_SURFACES,
        //             new_segment._cmfd_surface_bwd,new_segment._cmfd_surface_bwd % HEX_NUM_SURFACES);
      }
    }

    /* Truncate the linked list for the LocalCoords  把 LocalCoords 链表后面的节点都删掉，释放内存。*/
    start.prune();
    end.prune();
  }
```

主要负责 **将 MOC（特征线法）的轨迹线段（Track Segment）与 CMFD（粗网有限差分）网格的表面进行关联**。

在 MOC 计算中，我们需要统计穿过 CMFD 粗网格边界的中子流，以便进行 CMFD 加速。这段代码的作用就是给每个线段（Segment） 标注：在 CMFD 网格里，它从哪个面进、哪个面出。

<img src="img\image-20251102152727573.png" alt="image-20251102152727573" style="zoom:50%;" />

```c++
/*
 * 根据给定的空间坐标，找到它所属的 CMFD（粗网有限差分）网格单元的索引。
 */
int Cmfd::findCmfdCell(LocalCoords *coords)
{
    Point *point = coords->getHighestLevel()->getPoint();
    int global_cmfd_cell = _lattice->getLatticeCell(point); // 计算这个 point 落在 Lattice 的哪一个格子（Cell）里，并返回该格子的全局一维索引。

    if (_hexlattice_enable) // 如果是六角形网格模式（_hexlattice_enable 为真），直接返回刚才算出的全局索引。
        return global_cmfd_cell;
    else // 处理笛卡尔网格 (Cartesian/Rectangular)：普通的矩形网格，将全局索引转换为当前进程视角的本地索引。
    {
        /**
       	* 背景：在大规模并行计算（MPI）中，整个CMFD被切分成多个“域”（Domain），每个 CPU 核心只负责计算其中一块区域。
       	* getLocalCMFDCell将全局 CMFD 单元索引转换为当前 MPI 进程的本地索引。如果归我管，就要把它映射到本地数组的索引。本地数组的大小是 _local_num_x * _local_num_y * _local_num_z。
       	*/
        int local_cmfd_cell = getLocalCMFDCell(global_cmfd_cell);
        return local_cmfd_cell;
    }
```

##### _geometry->initializeFSRVectors();

```c++
/**
   * 在前面的分段过程中，FSR 是这样被记录的：
   * _FSR_keys_map：key → FSRData*
   * key：一个字符串，唯一描述一个 FSR 在 CSG 结构中的位置；
   * FSRData 结构体里有 _fsr_id、_mat_id、_point（质心）等。
   *
   * 但是后续求解时，我们通常是通过“fsr_id”去查：
   * 这个编号对应的 key 是啥？（有时要 debug 用）
   * 他的特征点centroids是哪个？
   * 它的材料 id 是多少？
   * 它在 CMFD 的哪个粗网格 cell 里？
   * initializeFSRVectors 就是把FSRData 结构体的数据，整理成几个“数组/向量”，让你可以通过 fsr_id O(1) 访问到对应信息。
   */
void Geometry::initializeFSRVectors()
{

    /* Get keys and values from map */
    log::info("Initializing FSR lookup vectors");
    auto key_list = _FSR_keys_map.keys();
    auto value_list = _FSR_keys_map.values();

    /* Allocate vectors */
    /* 分配向量 */
    size_t num_FSRs = _FSR_keys_map.size();
    _FSRs_to_keys = std::vector<std::string>(num_FSRs);           // FSR_ID与key的对应关系
    _FSRs_to_centroids = std::vector<Point *>(num_FSRs, nullptr); // FSR_ID与centroids的对应关系,不在这里赋值
    _FSRs_to_material_IDs = std::vector<int>(num_FSRs);           // FSR_ID与material_IDs的对应关系
    _FSRs_to_CMFD_cells = std::vector<int>(num_FSRs);             // FSR_ID与CMFD_cells的对应关系 它在 CMFD 的哪个粗网格 cell 里？
    // 创建四个向量，分别存储FSR ID到key、质心点、材料ID和CMFD单元的映射关系

    #pragma omp parallel for
    // 是一个 OpenMP 指令，告诉编译器“把接下来的 for 循环拆成多个线程同时执行”。编译时如果启用了 OpenMP，这条指令会让循环体（通常是紧跟着的 for (...) { ... }）自动分配给多个线程，提高并行性能；如果没启用 OpenMP，它会被当作普通注释忽略，不影响串行执行。
    for (size_t i = 0; i < num_FSRs; i++) // 填充_FSRs_to_keys, _FSRs_to_material_IDs
    {
        auto key = key_list[i];
        auto fsr = value_list[i];
        long fsr_id = fsr->_fsr_id;
        _FSRs_to_keys.at(fsr_id) = key;
        _FSRs_to_material_IDs.at(fsr_id) = fsr->_mat_id;
    }

    /* Add cmfd information serially */
    /* 串行添加CMFD信息 */
    if (_cmfd != nullptr)
    {
        for (size_t i = 0; i < num_FSRs; i++)
        {
            auto fsr = value_list[i];
            auto fsr_id = fsr->_fsr_id;
            // 每个 FSRData 里有 _cmfd_cell（所在的粗网格单元）和 _point（特征点）；
            auto point = fsr->_point;
            /*
        下面两行对CMFD方法非常重要：建立FSR与CMFD单元之间的双向映射关系
        1. 将FSR添加到对应的CMFD单元中
        2. 记录每个FSR所属的CMFD单元
        std::vector<std::vector<long>> _cell_fsrs;
        这个是二维vector数组，第一维确定是具体的哪个CMFD单元，第二个是单个CMFD中所关联的所有FSR
        */
            _cmfd->addFSRToCell(fsr->_cmfd_cell, fsr_id);     // 让 CMFD 知道这个 cell 里有哪些 FSR
            _FSRs_to_CMFD_cells.at(fsr_id) = fsr->_cmfd_cell; // fsr_id与CMFD cell的对应关系
            log::fdebug("FSR %ld Point x = %.2f, y = %.2f, z = %.2f", fsr_id, point->getX(), point->getY(), point->getZ());
            log::fdebug("cmfd cell is %d, fsr id is %d", fsr->_cmfd_cell, fsr_id);
        }
        if (_cmfd->GetHexLatticeEnable()) // 00000000
            _cmfd->findEmptyCmfdCells();
    }

    //_cmfd->printCellFSRs();

    /* Check if extruded FSRs are present 检查是否存在轴向挤出FSR*/
    size_t num_extruded_FSRs = _extruded_FSR_keys_map.size();
    if (num_extruded_FSRs > 0)
    {

        /* _extruded_FSR_keys_map 是 key → ExtrudedFSR* 的 map；
      同样，我们希望通过 extruded_fsr_id 直接拿到对应的 ExtrudedFSR*/
        _extruded_FSR_lookup = std::vector<ExtrudedFSR *>(num_extruded_FSRs);
        auto extruded_value_list = _extruded_FSR_keys_map.values();

        #pragma omp parallel for
        for (size_t i = 0; i < num_extruded_FSRs; i++)
        {
            long fsr_id = extruded_value_list[i]->_fsr_id;
            _extruded_FSR_lookup[fsr_id] = extruded_value_list[i];
        }

        delete[] extruded_value_list;
    }

    /* Delete key and value lists */
    delete[] key_list;
    delete[] value_list;

    // Output storage requirement of FSRs
    // 输出FSR的存储需求
    printNumFSRs();
    printMemUsageFSRs();

    log::verbose_once("Finished initialization of FSR lookup arrays");
}
```

<img src="img\image-20260102181454633.png" alt="image-20260102181454633" style="zoom:50%;" />

```c++

  /**
   * @brief Add an FSR ID to a vector that contains all the FSR IDs
   *        contained within a CMFD mesh cell.
   * @param The CMFD cell ID.
   * @param The FSR ID.
   */
  void Cmfd::addFSRToCell(int cmfd_cell, long fsr_id)
  {
    if (_hexlattice_enable)
    { // 0000000000000
      if (CellinHexLattice(cmfd_cell))
      {
        _cell_fsrs.at(cmfd_cell).push_back(fsr_id);
      }
    }
    else
    {
      /** _cell_fsrs是二维vector数组，第一维确定是具体的哪个CMFD单元，第二个是单个CMFD中所关联的所有FSR
       * std::vector<std::vector<long>> _cell_fsrs;
       * */
      _cell_fsrs.at(cmfd_cell).push_back(fsr_id); // 把fsr_id数据添加到具体的cmfd_cell单元相关联
    }
  }
```

### solver->computeEigenvalue

整个程序求解，其中涉及CMFD的求解

<img src="img\image-20251102160108685.png" alt="image-20251102160108685" style="zoom:50%;" />

<img src="img\image-20251102160507749.png" alt="image-20251102160507749" style="zoom:50%;" />

进入到**Solver::initializeCmfd()**，这里是==初始化CMFD相关的矩阵和矢量对象(这部分是重点)==

<img src="img\image-20251102161027891.png" alt="image-20251102161027891" style="zoom:50%;" />

<img src="img\image-20251102161339580.png" alt="image-20251102161339580" style="zoom:50%;" />

回到最初的solver->computeEigenvalue

<img src="img\image-20251102161734278.png" alt="image-20251102161734278" style="zoom:50%;" />

选择基于CPU的输运扫描算法。

<img src="img\image-20251102162013493.png" alt="image-20251102162013493" style="zoom:50%;" />

<img src="img\image-20251102162342292.png" alt="image-20251102162342292" style="zoom:50%;" />

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251102163001603.png" alt="image-20251102163001603" style="zoom: 50%;" />

<img src="img\image-20251102163819493.png" alt="image-20251102163819493" style="zoom:50%;" />

后续结合MOC代码和文档来看

<img src="img\image-20251102164814797.png" alt="image-20251102164814797" style="zoom:50%;" />

<img src="img\image-20251102165115531.png" alt="image-20251102165115531" style="zoom:50%;" />

下面这个函数通过遍历轨迹分段求解 MOC 方程，累加各分段对 CMFD 表面的中子流贡献。

<img src="img\image-20251102165843470.png" alt="image-20251102165843470" style="zoom:50%;" />

```c++
  /**
   * @brief Tallies the current contribution from this segment across the
   *        the appropriate CMFD mesh cell surface.
   *用于记录当前段（segment）在适当的CMFD网格单元表面上的中子流贡献
   * @param curr_segment a pointer to the Track segment of interest 指向轨迹分段的指针
   * @param azim_index the azimuthal index for this segment 此分段的方位角索引
   * @param polar_index the polar index for this segmenbt 此分段的极角索引
   * @param track_flux a pointer to the Track's angular flux 指向轨迹角通量的指针
   * @param fwd boolean indicating direction of integration along segment 沿段方向（前进或后退）
   */
  void CPUSolver::tallyCurrent(segment *curr_segment, int azim_index,
                               int polar_index, float *track_flux,
                               bool fwd)
  {

    /* Tally surface currents if CMFD is in use */
    /* 如果使用CMFD，则统计表面中子流 */
    if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    {
      if (_cmfd->GetHexLatticeEnable()) // 00000000000
        _cmfd->tallyHexCurrent(curr_segment, track_flux, azim_index, polar_index, fwd);
      else
        _cmfd->tallyCurrent(curr_segment, track_flux, azim_index, polar_index, fwd); // 四边形
    }
  }
```

<img src="img\image-20251102171245371.png" alt="image-20251102171245371" style="zoom:50%;" />

回到solver->computeEigenvalue ，判断是否开启CMFD。

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251102171614835.png" alt="image-20251102171614835" style="zoom:50%;" />

进入到 ==_cmfd->computeKeff(重点)==，这里先通过`collapseXSO`将 MOC 能群截面归并为 CMFD 能群截面（压缩 MOC 细网格数据到 CMFD 粗网格），然后构建矩阵`A`和`M`，采用幂迭代法求解特征值方程，得到有效增殖系数`k_eff`。

<img src="img\image-20251102172605836.png" alt="image-20251102172605836" style="zoom:50%;" />

这部分使用幂迭代法求解特征值问题，与前面1.3节原理相对应。<img src="img\image-20251102173142075.png" alt="image-20251102173142075" style="zoom:50%;" />

# 问题

1. void Geometry::segmentizeExtruded(Track *flattened_track,

​                  std::vector<double> z_coords, bool FSRonly)



1. ![image-20251210212550953](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251210212550953.png)

按道理从终点应该也是向反方向移动才对，应该终点已经在下一个cell中了。

// 为了避免起点/终点正好落在 CMFD 网格的顶点或面上，会先往反方向退一点点（TINY_MOVE）； 这有助于判断“是穿过哪一个面”。**这一步应该没有必要，写了是不是起了反效果？**，因为在findnextcell方法中会把终点向给定方向移动最小距离，肯定跑到下一个cell中去了。

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218135233695.png" alt="image-20251218135233695" style="zoom:33%;" />

除非位置关系是这样的。

![image-20260102122938080](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20260102122938080.png)

2.

void Geometry::segmentizeExtruded(Track *flattened_track,

​                  std::vector<double> z_coords, bool FSRonly)

这个version_num的用处不理解



<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218140124550.png" alt="image-20251218140124550" style="zoom:33%;" />

![image-20251218140417930](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218140417930.png)

按道理来说retrieved_coords = _extruded_FSR_keys_map.at(fsr_key)->_coords;它获得的局部点应该是start点才对。文档中说其他线程的局部点？只有传递给 Geometry::createUniqueExtrudedFSR的函数的start参数不同才可能，但是又是通过fsr_key在extruded_FSR_keys_map中取，应该不可能取到其他点吧？

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218141534603.png" alt="image-20251218141534603" style="zoom:33%;" />

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251212141008897.png" alt="image-20251212141008897" style="zoom:50%;" />

图5.45的问题

1. 最后想要的分段效果，应该是重叠平面的分段效果？
2. 但thread1为什么少了一部分点？没叠加下来
3. 还是前面那个点，结合这个图该怎么理解？thread 0 和 thread 1 在不同时间、从不同的 2D 轨迹到达了同一块投影区域，它们各自在那一刻调用 createUniqueExtrudedFSR，所以 _coords 可能是 thread 0 的点，而当前 start 是 thread 1 的点。但是fsr_key相同又属于同一个区域，他们是csg层次肯定是一样的吧





3.

over_mesh与轴向网的区别？

![image-20251212145736516](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251212145736516.png)

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218133610318.png" alt="image-20251218133610318" style="zoom:50%;" />





4.

这里是六个面的和吗？

> 从符号上讲, ${S}_{j}$ 是 CMFD 单元 $j$ 的“边界闭合曲面”，对应高斯散度定理里的那个“把体积${V}_{j}$ 围起来的曲面”。

![image-20251217195008971](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251217195008971.png)

若是按前面的理论，这里的 **H=6**，这些表面形成单元 j 与另一个 CMFD 单元之间的表面就说的通了。

![image-20251218142405822](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218142405822.png)





5.

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218145007603.png" alt="image-20251218145007603" style="zoom:33%;" />

这个面是cmfd单元的一个面吧？这里这个图不能直观的体现他们的关系，按道理该是三维的。这个面把它看成cmfd单元的右侧面好理解一点





6.

这一步的推导

![image-20251218152938578](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218152938578.png)

B. 19 本质上就是把 B. 18 里的**"面元积分 + 角度积分"，用 MOC 的离散角度和轨迹 **来做数值近似，最后把几何因子 $\left( {\Omega  \cdot  n}\right)$ 消掉。

从 B.18 出发(对某个界面 ${S}_{j,h}$ )：${\int }_{{S}_{j,h}}{dS}{\mathbf{J}}_{g}\left( \mathbf{r}\right)  \cdot  \mathbf{n} = {\int }_{{S}_{j,h}}{dS}{\int }_{4\pi }{d\Omega }{\psi }_{g}\left( {\mathbf{r},\Omega }\right) \left( {\Omega  \cdot  \mathbf{n}}\right) .$

1. 角度积分离散 (求积公式)
    MOC 用离散方向 ${\Omega }_{t}$ 和角权重 ${\alpha }_{t}$ 近似：${\int }_{4\pi }{d\Omega }{\psi }_{g}\left( {\mathbf{r},\Omega }\right) \left( {\Omega  \cdot  \mathbf{n}}\right)  \approx  \mathop{\sum }\limits_{t}{\alpha }_{t}{\psi }_{g}\left( {\mathbf{r},{\Omega }_{t}}\right) \left( {{\Omega }_{t} \cdot  \mathbf{n}}\right) .$
    代回去：${\int }_{{S}_{j,h}}{dS}{\mathbf{J}}_{g} \cdot  \mathbf{n} \approx  \mathop{\sum }\limits_{t}{\alpha }_{t}{\int }_{{S}_{j,h}}{dS}{\psi }_{g}\left( {\mathbf{r},{\Omega }_{t}}\right) \left( {{\Omega }_{t} \cdot  \mathbf{n}}\right) .$

2. 面积分用“轨迹穿透的面积”离散

  对固定方向 ${\Omega }_{t}$ ，穿过界面 ${S}_{j,h}$ 的轨迹段记为==$\left( {t,s}\right)  \in  {S}_{j,h}$==。图 B-3 给出几何关系：轨迹的横截面积为$\delta {A}_{t}$ ，它在该界面上“覆盖/穿透”的面积元可取
  $\delta {S}_{t,s} = \frac{\delta {A}_{t}}{\left| {\Omega }_{t} \cdot  \mathbf{n}\right| }.$

  因此${\int }_{{S}_{j,h}}{dS}{\psi }_{g}\left( {\mathbf{r},{\Omega }_{t}}\right) \left( {{\Omega }_{t} \cdot  \mathbf{n}}\right)  \approx  \mathop{\sum }\limits_{{\left( {t,s}\right)  \in  {S}_{j,h}}}{\psi }_{g}^{t,s}\left( {s}_{j,h}\right) \left( {{\Omega }_{t} \cdot  \mathbf{n}}\right) \delta {S}_{t,s}.$

  > 这里的(t,s)是离散索引对，用来表示“哪些离散轨迹(及其段)穿过界面 ${S}_{j,h}$ ”：
  >
  > - $t$ ：离散方向编号(对应 ${\Omega }_{t}$ 和角权重 ${\alpha }_{t}$ )。
  > - $s$ : 不是路径长度，而是轨迹编号 (第几条 track) 或轨迹段编号 (某条 track 在某个 cell 里的那一段)。
  > - 所以$\mathop{\sum }\limits_{{\left( {t,s}\right)  \in  {S}_{j,h}}}\left( \cdots \right)$意思是：对所有“方向为 $t$ 且该方向下第 $s$ 条(或第 $s$ 段)轨迹”与界面 ${S}_{j,h}$ 有交的那些对象做求和。
  >
  > 沿方向 ${\Omega }_{t}$ 的特征线可写成$\mathbf{r}\left( s\right)  = {\mathbf{r}}_{0} + {\Omega }_{t}s$。于是角通量在这条线上就是 ${\psi }_{g}\left( {\mathbf{r}\left( s\right) ,{\Omega }_{t}}\right)$ ，常简写成 ${\psi }_{g}^{t}\left( s\right)$ 。
  >
  > **界面处的通量就是“在交点处取值”，当这条轨迹穿过界面 ${S}_{j,h}$ 时，会有一个交点 ${\mathbf{r}}_{j,h}^{t,s}$ 。这个交点对应的沿特征线的坐标就是 ${s}_{j,h}$** 。所以${\psi }_{g}^{t,s}\left( {s}_{j,h}\right)  \equiv  {\psi }_{g}\left( {{\mathbf{r}}_{j,h}^{t,s},{\Omega }_{t}}\right)$

  把 $\delta {S}_{t,s}$ 代入就得到关键的“消去”：$\left( {{\Omega }_{t} \cdot  \mathbf{n}}\right) \delta {S}_{t,s} = \left( {{\Omega }_{t} \cdot  \mathbf{n}}\right) \frac{\delta {A}_{t}}{\left| {\Omega }_{t} \cdot  \mathbf{n}\right| } = \delta {A}_{t}\operatorname{sign}\left( {{\Omega }_{t} \cdot  \mathbf{n}}\right) .$

  如果你的集合 $\left( {t,s}\right)  \in  {S}_{j,h}$ 按外法向只统计流出方向，那 ${\Omega }_{t} \cdot  n > 0$ ，上式就直接变成 $({\Omega }_{t}n)\delta {S}_{t,s} = \delta {A}_{t}$ ；若同时统计进/出，则符号由 $\operatorname{sign}\left( {{\Omega }_{t} \cdot  n}\right)$ 体现“净”流。

  于是整体变为${\int }_{{S}_{j,h}}{dS}{\mathbf{J}}_{g} \cdot  \mathbf{n} \approx  \mathop{\sum }\limits_{t}{\alpha }_{t}\mathop{\sum }\limits_{{\left( {t,s}\right)  \in  {S}_{j,h}}}\delta {A}_{t}{\psi }_{g}^{t,s}\left( {s}_{j,h}\right) .$

3. 把权重合并成 ${w}_{t,s}$

  文中用式 (2.16) 定义轨迹权重${w}_{t} = \delta {A}_{t}{\alpha }_{t}$，于是就得到 B.19：${\int }_{{S}_{j,h}}{dS}{\mathbf{J}}_{g}\left( \mathbf{r}\right)  \cdot  \mathbf{n} \approx  \mathop{\sum }\limits_{{\left( {t,s}\right)  \in  {S}_{j,h}}}{w}_{t,s}{\psi }_{g}^{t,s}\left( {s}_{j,h}\right) $

  > 注意这里内层的集合写成 $\left( {t,s}\right)  \in  {S}_{j,h}$ ，其实就已经在暗示：对每个方向 $t$ 只把该方向下穿过界面 ${S}_{j,h}$ 的那些轨迹段 $s$ 拿来求和。
  >
  > 把权重合并 ${w}_{t} = {\alpha }_{t}\delta {A}_{t}$ 后：$\mathop{\sum }\limits_{t}\mathop{\sum }\limits_{{s \in  {\mathcal{S}}_{j,h}\left( t\right) }}{w}_{t}{\psi }_{g}^{t,s}\left( {s}_{j,h}\right)$
  >
  > 现在如果我们把“所有穿过该界面的轨迹段”统一收集成一个大集合${\mathcal{T}}_{j,h} = \left\{  {\left( {t,s}\right)  \mid  \text{ 轨迹段 }\left( {t,s}\right) \text{ 穿过 }{S}_{j,h}}\right\}  ,$
  >
  > 那么上面的双重求和就等价于一个“单重求和”：$\mathop{\sum }\limits_{{\left( {t,s}\right)  \in  {\mathcal{T}}_{j,h}}}{w}_{t}{\psi }_{g}^{t,s}\left( {s}_{j,h}\right)$

**为什么这里要变成界面处的通量(就是“在交点处取值”)**

因为你现在算的是穿过界面 ${S}_{j,h}$ 的面电流：${\int }_{{S}_{j,h}}{dS}{\psi }_{g}\left( {\mathbf{r},{\Omega }_{t}}\right) \left( {{\Omega }_{t} \cdot  \mathbf{n}}\right)$

这里的被积函数 ${\psi }_{g}\left( {\mathbf{r},{\Omega }_{t}}\right) \left( {{\Omega }_{t} \cdot  n}\right)$ 是定义在“面上”的量：你是在对所有位于该界面上的点 $\mathbf{r} \in  {S}_{j,h}$ 积分，所以离散化时自然要用界面上的通量值。更具体地说，MOC 做 B.19 这一步时隐含了一个标准的数值近似:

1. 用轨迹把界面分成很多小“面片”
    每条方向为 ${\Omega }_{t}$ 的轨迹束(横截面积 $\delta {A}_{t}$ )穿过界面时，对应界面上的一小块面积$\delta {S}_{t,s} \approx  \frac{\delta {A}_{t}}{\left| {\Omega }_{t} \cdot  n\right| }$，这等价于把整个界面 ${S}_{j,h}$ 划分/覆盖成很多个小面积片 $\delta {S}_{t,s}$ 。

2. 在每个小面积片上把通量当作“近似常数”，对某个小面片 $\delta {S}_{t,s}$ ，${\int }_{\delta {S}_{t,s}}{dS}{\psi }_{g}\left( {\mathbf{r},{\Omega }_{t}}\right) \left( {{\Omega }_{t} \cdot  n}\right)  \approx  {\psi }_{g}\left( {{\mathbf{r}}_{t,s}^{\text{hit }},{\Omega }_{t}}\right) \left( {{\Omega }_{t} \cdot  n}\right) \delta {S}_{t,s}$。==也就是用一个代表点(通常就选轨迹与界面的交点 ${\mathbf{r}}_{t,s}^{\mathrm{{hit}}}$ )==来做“面片常值近似”(collocation/中点法思想)。于是就出现你说的“在交点处取值”：${\psi }_{g}^{t,s}\left( {s}_{j,h}\right)  \equiv  {\psi }_{g}\left( {{\mathbf{r}}_{t,s}^{\text{hit }},{\Omega }_{t}}\right) .$

3. 物理上也必须是界面值
    CMFD 要的就是“这个面上穿过去多少粒子”(净流/漏泄)，那决定因素就是粒子到达界面时的角通量(以及方向与法向的夹角)。用体内某点的通量替代界面通量，会把“还没走到界面时的状态”拿来算穿面电流，几何和物理意义都不对。



7.

![image-20251218194852341](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251218194852341.png)

1. 对函数  $u(j,h)$  的理解
2. 补充的扩散理论的正确性





8.



<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251223204415817.png" alt="image-20251223204415817" style="zoom:50%;" />





9.

组件内栅元级粗网是什么意思？

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20251219132408546.png" alt="image-20251219132408546" style="zoom: 67%;" />







