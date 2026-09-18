# -*- coding: utf-8 -*-
"""
D-15 水库回水曲线和河道水面线计算程序 —— 内核
================================================
复刻《水利水电工程设计计算程序集》（公之于众版）D-15 程序
（原著作者：唐文华，水利部天津勘测设计研究院）。

一、功能
--------
河道恒定非均匀流水面线计算。把计算河段分为若干河段，自下而上逐段用试算法
解伯努里方程式求上游水位，依次推算到上游边界。支持：
  · 干流 + 一级支流组成的羽状河系（k1 条河流）；
  · 复式断面（主槽 + 滩地），主槽/滩地糙率可不同；
  · 断面资料为「宽度 b ~ 高程 h」曲线（kba01=0）或「面积 a ~ 高程 h」曲线
    （kba01=1），结点数可逐断面不等；
  · 同一河流上不同流量河段、不同 n 值河段；
  · 计算方法 1（用河段 n 值，kbhs=1）/ 方法 2（用断面 n 值，kbhs=2）；
  · 已知水位反推河段 n 值（ncp=1，式 5）；
  · 断面扩大/收缩水头损失（kwusan=1，式 6/7 + E 值表）；
  · 桥、闸等建筑物局部水头损失 brg；
  · 高程换算系数 threes；
  · 河床坡度 I 大于 / 等于 / 小于 0 均可。

二、计算原理（公式取自说明书 D-15Intro.rtf 内嵌 Equation(OLE) 对象，
   经 MTEF v3 解包后逐条还原，脚本 _d15_eq_extract.py）
-------------------------------------------------------------------
伯努里方程（1-1 为下游已知断面，2-2 为上游待求断面，动能修正系数 α1=α2=1）：

    Z2 = Z1 + (V1² − V2²)/(2g) + h沿 + h桥 + h扩                     (1)

沿程损失（式 2/3/4）：

    h沿 = i_f · ΔL                                                  (4)
    i_f = (i1 + i2)/2                                               (4)
    i_k = Q² / K_k²          （k = 1, 2）                            (4)
    K_k = A_k · R_k^(2/3) / n      （流量模数，Manning 形式）          (2′)
        （等价形式 h沿 = n²V²ΔL/R^(4/3)，即式 3）
    A_k(或 R_k) 的脚标 1、2 分别对应 1-1、2-2 断面，ΔL 为断面间距。

    **教材/手册同源印证**（本内核选型依据）：
      · 《水力学》（吴持恭 第四版）式6.7 摩阻坡 J_f = n²v²/R^(4/3)
        —— 与上式(3) 逐字一致；
      · 同书式6.5/6.6 渐变流微分方程 dh/ds=(i−J_f)/(1−Fr²) 与标准步进
        Δs=(E_s2−E_s1)/(i−J̄_f)，J̄_f 为「段平均摩阻坡」—— D-15 式(4) 的
        i_f=(i₁+i₂)/2 即该「段平均」的一种取法；
      · 《水工设计手册》卷1 §3.9.6 水面曲线「分段求和/数值积分」、
        卷2 §5.11.1 回水控制方程 dz/dl+(α+ξ)d(v²/2g)/dl+v²/(C²R)=0、
        卷7 §1.5.2.4 式1.5-17/18 分段求和法 ΔL=[(h+v²/2g)₂−(h+v²/2g)₁]/(i−J̄)、
        J̄=n²v̄²/R̄^(4/3) —— 同一算法族。
      数值对照（见三.2 表）：教材「段平均摩阻坡」J̄=n²v̄²/R̄^(4/3) 的取法
      在本例 RMS 4.4%~5.7%，远劣于原著式(4) 的 0.323%，故分段平均取法
      以原著式(4) 为准。

扩大/收缩损失（式 6/7 + E 值表）：

    h扩 = E · (V1 − V2)² / (2g)                                      (6)
    V2/V1 ≥ 1（扩大）: E = (V2/V1 − 1)²                              (7)
    V2/V1 <  1（收缩）: E 按表插值（0.00→0.50, 0.01→0.50, 0.10→0.45,
                       0.20→0.40, 0.40→0.30, 0.60→0.20, 0.80→0.10, 1.00→0.00）

断面 n 值反推（ncp=1，式 5）：

    n = A · R^(2/3) · i_f^(1/2) / Q,   A=(A1+A2)/2, R=(R1+R2)/2
    i_f = [(Z2 − Z1) − (V1² − V2²)/(2g)] / ΔL

三、断面要素的（反演）口径 —— 关键裁决
--------------------------------------
1) 过水面积 A：由 b~h 曲线对高程作梯形积分（水面截断于曲线内），
   或直接取 a~h 曲线插值。**该口径经权威 .OUT 逐断面反证成立**：
   892/A 与输出 V 列逐断面相符（21/21，偏差 ≤0.5%，见 _d15_verify_run.txt）。

2) 湿周 P 与水力半径 R：说明书正文与公式对象**均未给出**湿周算式
   （文本层缺失）。**先查教材/手册口径**（《水力学》吴持恭第四版 式1.4
   R=A/χ、式6.1 梯形 A=(b+mh)h、式6.2 梯形 χ=b+2h√(1+m²)、式6.6
   渐变流标准步进、式6.7 J_f=n²v²/R^(4/3)；《水工设计手册》卷1 §3.9.6、
   卷2 §5.11.1 式5.10-1、卷7 §1.5.2.4 式1.5-17/18 分段求和法），
   再以 20 个河段的沿程损失方程（h沿 = (i1+i2)ΔL/2）作最小二乘反演。
   教材式6.2 需要岸坡 m，而 D-15 输入只有 b~h 曲线；按「等效梯形」
   （由 A、B、最大水深 d 反解 b0=2A/d−B、m=(B−b0)/(2d)）代入式6.2，
   本例等效 m=12~107、b0 多次为负 —— 断面是陡岸天然断面，式6.2 退化为
   χ≈B（比 B 仅大 0.13~0.42 m），代入阻力方程 RMS 2.57%，远劣于下表
   最优口径，**故教材梯形式不适用**。

   湿周模型族 P = B + β·(A/B)（B 为水面宽，A/B 为平均水深）在
   **β = 2.000** 处取到尖锐极小（RMS 0.323%）；β=0（P=B）RMS 2.71%、
   β=1 (P=B+A/B) 1.39%、β=2.5 0.75%，均显著劣于 β=2。
   **与教材的关系**：P = B + 2h 正是教材式6.2 的 **m=0（矩形）特例**，
   此处取「**同宽同面积的等效矩形**」（h = Ā = A/B），而非「同宽同最大
   水深」（后者 P=B+2d，RMS 2.71%）。故本内核口径是**教材公式的合法特例**，
   教材出处：吴持恭《水力学》第四版 §6.2 式6.2（m=0）+ §5.1 式1.4 R=A/χ。
   故本内核取：

        P = B + 2·Ā,   Ā = A/B（平均水深）,   R = A/P

   湿周口径 × 教材标准式 的逐项数值对照（同一目标方程、同一 20 河段）：

     口径                                          RMS(%)   bias(%)   均值 χ/B
     ─────────────────────────────────────────────────────────────────────────
     教材等效梯形 χ=b0+2d√(1+m²)（式6.1+6.2）        2.5655   −2.4247    1.0009
     教材矩形 χ=B+2d（式6.2 令 m=0，取最大水深 d）     2.7102   +2.3810    1.0375
     χ=B（不修正，≈ b~h 曲线弧长）                    2.7068   −2.5464    1.0000
     **χ=B+2·(A/B)（本内核；式6.2 m=0，取等面积水深）** **0.3230**  **−0.0593**  **1.0190**
     ─────────────────────────────────────────────────────────────────────────
   另：教材/手册的「段平均摩阻坡」J̄=n²v̄²/R̄^(4/3)（卷7 式1.5-18）取
   算术平均 v̄、R̄ 时 RMS 4.42%、取 Ā/χ̄ 时 RMS 5.70%，均远劣于 D-15
   说明书式(4) 的 i_f=(i₁+i₂)/2（RMS 0.323%）—— 佐证原著确用后者，
   教材式（式6.6/6.7、卷7 式1.5-17/18）给的是**算法结构**，具体分段
   平均取法以原著式(4) 为准。
   证据脚本：_d15_chi_probe.py / _d15_chi_probe2.py / _d15_chi_probe3.py
   输出：_d15_chi_probe_out.txt / _d15_chi_probe2_out.txt / _d15_chi_probe3_out.txt
   替换量化（_d15_chi_replace.py / _d15_chi_replace_out.txt）：同一求解器仅换 χ
   时，正演水面线 vs 权威 OUT 的最大水位偏差 —— 本内核 13.5 mm（RMS 5.9 mm）、
   教材梯形 49.9 mm（RMS 24.6 mm）、教材矩形取最大水深 65.2 mm（RMS 27.9 mm），
   即改用教材式会恶化 3.7~4.8 倍，**故不替换**。

   （宽矩形断面下 P=B+2h 为精确式；对天然断面取「同宽同面积等效矩形」
     是教材式6.2 m=0 特例的合理外推，亦与该年代程序常用工程近似一致。）
   备选模型 P = B + α·d（d 为最大水深）的最优 α=1.005（RMS 0.443%）、
   P = B + a·d + b·(A/B) 二维网格最优 (0, 2.2)（正演 max|ΔH|=9.8 mm），
   均不优于 (0, 2.0)，故不采用。

3) 流速 V = Q/A，与输出的 V 列一致（21/21）。

四、输入数据（自由格式，逗号/空白分隔；顺序经 D-15.INT 反演并与 .OUT
   回显列一一对应、数值流恰好耗尽 1268 个值，无剩余）
--------------------------------------------------------------
 (1)  k1, k2, k3, k4, k6, k7, k8   河流数 / 最大断面数 / 主槽最多n河段数 /
                                   最多流量河段数 / 主槽bs~hs最大结点数 /
                                   滩地最多n河段数 / 滩地bss~hss最大结点数
 (2)  mmax, kbhs, kwusan           允许迭代次数 / 计算方法(1|2) / 扩大收缩(0|1)
 (3)  kba01, ncp                   0=b~h曲线,1=a~h曲线 / 1=已知水位反推n值
 (4)  epsh, threes, dh0, as1, ds1  允许水位迭代差 / 高程换算系数 /
                                   预给河段水位差 / 调整步长系数 / 变化as1参数
 (5)  nrqil(j)   j=1..k1            不同流量河段数
 (6)  nrcml(j)   j=1..k1            主槽不同n值河段数
 (7)  nrcmsl(j)  j=1..k1            滩地不同n值河段数
 (8)  nb1(j)     j=1..k1            河道计算断面数
 (9)  nq1(k,j)   k=1..nrqil(j)      不同流量河段上游界面序号
(10)  ncm1(k,j)  k=1..nrcml(j)      主槽不同n值河段上游界面序号
(11)  ncms1(k,j) k=1..nrcmsl(j)     滩地不同n值河段上游界面序号
(12)  q(k,j)     k=1..nrqil(j)      河段流量
(13)  cmn(k,j)   k=1..nrcml(j)      主槽河段n值
(14)  cmns(k,j)  k=1..nrcmsl(j)     滩地河段n值（无滩地时给 0.0001 之类的哑值）
(15)  mncs(k,j)  k=1..nb1(j)        主槽断面 bs~hs 曲线结点数
(16)  mncss(k,j) k=1..nb1(j)        滩地断面 bss~hss 曲线结点数
(17)  dl(l,j)    l=1..nb1(j)−1      断面距离（**km**）
(18)  h(i,j)     i=1..nb1(j)        断面水位；h(1,j) 为起始水位（见注）
(19)  bs(l,k,j)  每断面 mncs(k,j) 个  主槽宽度结点
(20)  as(l,k,j)  每断面 **1** 个     主槽面积结点（kba01=0 时为哑值 0.0）
(21)  hs(l,k,j)  每断面 mncs(k,j) 个  主槽高程结点
(22)  bss(l,k,j) 每断面 mncss(k,j) 个 滩地宽度结点
(23)  ass(l,k,j) 每断面 mncss(k,j) 个 滩地面积结点
(24)  hss(l,k,j) 每断面 mncss(k,j) 个 滩地高程结点
(25)  [可选] brg(k,j) k=1..nb1(j)    桥、闸等建筑物水头损失（样例文件未含此块）

注：说明书写「i=nb1(1)，j=1 时 h 是起始水位」，但按样例 D-15.INT/.OUT 反证，
    起始水位即 h(1,1)=366.0（输出成果表第 1 行水位恰为给定的 366.000 且不参与
    迭代）。本内核据此以 h(1,j) 为下游起始水位、自 1 断面向 nb1 断面上游推算。

五、输出表结构（原著 .OUT）
--------------------------
  先回显全部输入变量（K1..K8 / 控制参数 / 各数组），再按河流输出成果表：

    I  断面号（自下游向上游 1..NB）
    H  水位(m)      V  流速(m/s)     主槽ｎ值/滩地ｎ值  断面 n
    dh 沿程水头损失(m)   brg 桥闸损失(m)   vh 流速水头差 (V1²−V2²)/2g (m)
    dl 断面间距(km)      Q  流量(m³/s)
  末行为「计算结束  O.K.!」。

六、未闭合点（不可唯一反演）
---------------------------
  (a) 湿周/水力半径算式：说明书缺失（见三.2）。**已按教材口径复核**：
      教材式6.2 的梯形/矩形特例均劣于本内核口径（RMS 2.57%/2.71% vs
      0.323%），故保留 χ=B+2·(A/B)（= 教材式6.2 m=0 特例 + 同宽同面积
      等效矩形，出处：吴持恭《水力学》第四版 式1.4/式6.2）。残余
      20 河段 RMS 0.32%；正演水面线与权威 OUT 最大偏差 13.5 mm。
      关于「13.5 mm」的来源：**教材只给逐段试算法的方程结构**（式6.5/6.6
      渐变流微分方程与标准步进式、式6.7 摩阻坡 J_f=n²v²/R^(4/3)；卷7
      式1.5-17/18 分段求和法），**不给数值停代容差**；D-15 的
      epsh=0.010 m 是程序输入项（说明书表 4「允许水位迭代差 epsh」）。
      「试算按 |Δh| < ε 停止 ⇒ 打印值与迭代终态差为 O(ε)」是逐段试算法
      的固有性质，与教材结构一致；实测 max|ΔH| = 13.5 mm 与 epsh=10 mm
      同量级（同一数量级）。
      另：由输出 V 反推的 A 与打印水位处的梯形积分 A 散差 ±0.5%、无系统
      规律，说明打印 H 与计算 V/h 取自不同的迭代状态（V 与水头损失取自
      最后一轮**试算**水位，H 为**新的**水位）——这正是逐段试算法的典型
      停代行为。故**逐位复现打印水位在原理上不可达**，本内核取严格收敛解。
  (b) 高程换算系数 threes 的作用方式（样例为 0）：说明书未给出算式细节。
  (c) 式(6) h扩 的次幂在 MTEF 文本层不可唯一判读（该分支样例 kwusan=0，
      未被触发）。
  (d) 滩地（bss/hss）在样例中为哑值（单结点 5000.0/5000.0），滩地参与
      面积/湿周的方式无法由样例反演；本内核在 z 高于滩地首结点高程时按
      bss~hss 曲线补算滩地面积（并入 A）——该分支未被样例触发。
      教材对应口径为《水力学》第四版 §5.3「复式断面/粗糙度不同断面：
      分段算湿周加权（等效糙率法）」，与 D-15 输入中「主槽 n + 滩地 n
      分设、K=ΣωR^(2/3)/n」的结构一致；本内核按主槽/滩地分别建 K 后用
      i_f=(i₁+i₂)/2 合成，属该教材口径的等价实现（样例未触发，不可验证）。

七、验证
--------
  d15_verify.py 逐位对拍 data/D-15.OUT（GBK）；d15_smoke.py 为 SLCALC_HOME
  隔离冒烟。
"""
import math
import re

from ..core.intio import smart_read_text
from ..core.outgen import write_json, write_out

PROGRAM_ID = "D-15"
TITLE = "河道水面线计算 D-15"
AUTHOR = "唐文华（水利部天津勘测设计院）"

G = 9.81

# D-15 的 INT 文件含「892.」「0.」这类尾点无小数的写法，core.intio 的数值正则
# 要求小数点后必有数字会漏读（本例少读 22 个值），此处用更宽松的正则。
_NUM_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")


def _read_numbers(path):
    """读取 INT 数值流（逗号 / 换行 / 空白分隔，容忍尾点写法）。"""
    text = smart_read_text(path).lstrip("\ufeff").replace("\x1a", "")
    nums = []
    for tok in re.split(r"[,\s]+", text):
        if tok and _NUM_RE.fullmatch(tok):
            nums.append(float(tok))
    return nums

# 收缩（V2/V1 < 1）时的 E 值表（式 6 的系数，取自说明书表格）
E_TABLE = [(0.00, 0.50), (0.01, 0.50), (0.10, 0.45), (0.20, 0.40),
           (0.40, 0.30), (0.60, 0.20), (0.80, 0.10), (1.00, 0.00)]


# ============================================================
# 断面几何
# ============================================================

def _interp_b(bs, hs, z):
    """b~h 曲线在 z 处的宽度（线性插值，超出范围取端点）。"""
    if z <= hs[0]:
        return 0.0
    if z >= hs[-1]:
        return bs[-1]
    for i in range(len(hs) - 1):
        if hs[i] <= z <= hs[i + 1]:
            if hs[i + 1] == hs[i]:
                return bs[i + 1]
            return bs[i] + (bs[i + 1] - bs[i]) * (z - hs[i]) / (hs[i + 1] - hs[i])
    return bs[-1]


def _area_trap(bs, hs, z):
    """
    b~h 曲线对高程作梯形积分得到过水面积（水面在曲线内截断）。
    该口径经权威 D-15.OUT 的 V 列逐断面反证成立。
    """
    A = 0.0
    for i in range(len(hs) - 1):
        h0, h1 = hs[i], hs[i + 1]
        if z <= h0:
            break
        b0, b1 = bs[i], bs[i + 1]
        if z >= h1:
            A += 0.5 * (b0 + b1) * (h1 - h0)
        else:
            bb = b0 + (b1 - b0) * (z - h0) / (h1 - h0) if h1 != h0 else b0
            A += 0.5 * (b0 + bb) * (z - h0)
            break
    return A


def _area_interp(as_, hs, z):
    """a~h 曲线（kba01=1）插值取面积。"""
    if z <= hs[0]:
        return 0.0
    if z >= hs[-1]:
        return as_[-1]
    for i in range(len(hs) - 1):
        if hs[i] <= z <= hs[i + 1]:
            if hs[i + 1] == hs[i]:
                return as_[i + 1]
            return as_[i] + (as_[i + 1] - as_[i]) * (z - hs[i]) / (hs[i + 1] - hs[i])
    return as_[-1]


def section_props(z, sec, Q, kba01):
    """
    断面水力要素。
      A  过水面积 (m²)
      B  水面宽 (m)
      P  湿周 = B + 2·Ā（Ā = A/B，平均水深）  ← 见模块 docstring 三.2
      R  水力半径 = A/P
      V  流速 = Q/A
    复式断面：z 超过滩地首结点高程时，按 bss~hss 曲线补算滩地面积（并入 A）。
    """
    bs, hs = sec["bs"], sec["hs"]
    if kba01 == 1:
        A = _area_interp(sec.get("as_", [0.0]), hs, z)
    else:
        A = _area_trap(bs, hs, z)
    B = _interp_b(bs, hs, z)
    bss, hss = sec.get("bss", []), sec.get("hss", [])
    if bss and hss and len(bss) == len(hss) and z > hss[0]:
        A += _area_trap(bss, hss, z)
    if A <= 0.0 or B <= 0.0:
        return {"A": A, "B": B, "P": 0.0, "R": 0.0, "V": 0.0}
    mean_depth = A / B
    P = B + 2.0 * mean_depth
    R = A / P
    return {"A": A, "B": B, "P": P, "R": R, "V": Q / A}


def friction_slope(props, n, Q):
    """i = n²V²/R^(4/3)（式 3）；R→0 时取极大值以阻止搜索越界。"""
    R, V = props["R"], props["V"]
    if R <= 0.0:
        return 1e30
    return n * n * V * V / R ** (4.0 / 3.0)


def expansion_E(v1, v2):
    """扩大/收缩系数 E（式 7 / 说明书表格）。"""
    if v1 <= 0.0:
        return 0.0
    r = v2 / v1
    if r >= 1.0:
        return (r - 1.0) ** 2
    for i in range(len(E_TABLE) - 1):
        r0, e0 = E_TABLE[i]
        r1, e1 = E_TABLE[i + 1]
        if r0 <= r <= r1:
            if r1 == r0:
                return e0
            return e0 + (e1 - e0) * (r - r0) / (r1 - r0)
    return 0.0


# ============================================================
# 解析
# ============================================================

def parse(data):
    """
    解析输入。data: INT 文件路径 | dict（已是参数字典时原样返回）。
    数据顺序见模块 docstring 四。
    """
    if isinstance(data, dict):
        return data
    nums = _read_numbers(data)
    p = {"数值流长度": len(nums)}
    pos = 0

    def take(cnt=1):
        nonlocal pos
        v = nums[pos:pos + cnt]
        if len(v) < cnt:
            raise ValueError("INT 文件数据不完整（需要 %d 个，实际只剩 %d 个）"
                             % (cnt, len(v)))
        pos += cnt
        return v

    (k1, k2, k3, k4, k6, k7, k8) = [int(x) for x in take(7)]
    p.update({"k1": k1, "k2": k2, "k3": k3, "k4": k4, "k6": k6, "k7": k7, "k8": k8})
    mmax, kbhs, kwusan = take(3)
    p.update({"mmax": int(mmax), "kbhs": int(kbhs), "kwusan": int(kwusan)})
    kba01, ncp = take(2)
    p.update({"kba01": int(kba01), "ncp": int(ncp)})
    epsh, threes, dh0, as1, ds1 = take(5)
    p.update({"epsh": epsh, "threes": threes, "dh0": dh0, "as1": as1, "ds1": ds1})

    p["nrqil"] = [int(x) for x in take(k1)]
    p["nrcml"] = [int(x) for x in take(k1)]
    p["nrcmsl"] = [int(x) for x in take(k1)]
    p["nb1"] = [int(x) for x in take(k1)]
    p["nq1"] = [[int(x) for x in take(p["nrqil"][j])] for j in range(k1)]
    p["ncm1"] = [[int(x) for x in take(p["nrcml"][j])] for j in range(k1)]
    p["ncms1"] = [[int(x) for x in take(p["nrcmsl"][j])] for j in range(k1)]
    p["q"] = [take(p["nrqil"][j]) for j in range(k1)]
    p["cmn"] = [take(p["nrcml"][j]) for j in range(k1)]
    p["cmns"] = [take(p["nrcmsl"][j]) for j in range(k1)]
    p["mncs"] = [[int(x) for x in take(p["nb1"][j])] for j in range(k1)]
    p["mncss"] = [[int(x) for x in take(p["nb1"][j])] for j in range(k1)]
    p["dl"] = [take(p["nb1"][j] - 1) for j in range(k1)]
    p["h"] = [take(p["nb1"][j]) for j in range(k1)]

    p["bs"] = [[take(p["mncs"][j][k]) for k in range(p["nb1"][j])] for j in range(k1)]
    # kba01=0 时主槽面积结点为哑值，每断面仅 1 个（由 D-15.INT/.OUT 反证）
    as_cnt = [[1 if p["kba01"] == 0 else p["mncs"][j][k]
               for k in range(p["nb1"][j])] for j in range(k1)]
    p["as_cnt"] = as_cnt
    p["as_"] = [[take(as_cnt[j][k]) for k in range(p["nb1"][j])] for j in range(k1)]
    p["hs"] = [[take(p["mncs"][j][k]) for k in range(p["nb1"][j])] for j in range(k1)]
    p["bss"] = [[take(p["mncss"][j][k]) for k in range(p["nb1"][j])] for j in range(k1)]
    p["ass"] = [[take(p["mncss"][j][k]) for k in range(p["nb1"][j])] for j in range(k1)]
    p["hss"] = [[take(p["mncss"][j][k]) for k in range(p["nb1"][j])] for j in range(k1)]

    # (25) 桥闸水头损失：样例文件不含该块，有剩余才读
    p["brg"] = []
    for j in range(k1):
        nb = p["nb1"][j]
        if len(nums) - pos >= nb:
            p["brg"].append(take(nb))
        else:
            p["brg"].append([0.0] * nb)

    p["已消费数值个数"] = pos
    p["剩余数值"] = nums[pos:]
    return p


def _reach_of_section(bounds, nb):
    """
    由「河段上游界面序号」数组求每断面所属河段（0 基）。
    bounds[k] = 第 k+1 个河段的上游界面断面号（1 基）；
    第 k+1 个河段覆盖断面 [bounds[k-1], bounds[k]-1]。
    """
    out = [None] * nb
    lo = 1
    for k, hi in enumerate(bounds):
        for i in range(lo - 1, min(hi - 1, nb)):
            out[i] = k
        lo = hi
    for i in range(nb):
        if out[i] is None:
            out[i] = len(bounds) - 1
    return out


# ============================================================
# 计算
# ============================================================

def _residual(z2, sec1, sec2, z1, Q, n, L_m, kbhs, n2, brg, kwusan, kba01):
    """伯努里方程残差 F(z2)=0 的根即上游水位。"""
    h1 = section_props(z1, sec1, Q, kba01)
    h2 = section_props(z2, sec2, Q, kba01)
    vh = (h1["V"] ** 2 - h2["V"] ** 2) / (2.0 * G)
    if kbhs == 2 and n2 is not None:
        i1 = friction_slope(h1, n2, Q)     # 方法 2：用断面 n 值
        i2 = friction_slope(h2, n2, Q)
    else:
        i1 = friction_slope(h1, n, Q)      # 方法 1：用河段 n 值
        i2 = friction_slope(h2, n, Q)
    hf = 0.5 * (i1 + i2) * L_m
    he = expansion_E(h1["V"], h2["V"]) * (h1["V"] - h2["V"]) ** 2 / (2.0 * G) \
        if kwusan == 1 else 0.0
    return z1 + vh + hf + brg + he - z2


def _solve_z2(sec1, sec2, z1, Q, n, L_m, kbhs, n2, brg, kwusan, kba01):
    """
    试算法（自上而下扫描取第一个变号区间 + 二分）求 z2。
    扫描区间取 [河底, z1+10]，自上而下，取第一个（最高的物理）根。
    """
    zbed = min(sec2["hs"])
    lo = zbed + 1e-6
    hi = z1 + 10.0
    nstep = 4000
    prev_z, prev_f = hi, _residual(hi, sec1, sec2, z1, Q, n, L_m, kbhs, n2,
                                   brg, kwusan, kba01)
    a = b = None
    for i in range(1, nstep + 1):
        z = hi - (hi - lo) * i / nstep
        f = _residual(z, sec1, sec2, z1, Q, n, L_m, kbhs, n2, brg, kwusan, kba01)
        if prev_f * f <= 0:
            a, b = z, prev_z
            break
        prev_z, prev_f = z, f
    if a is None:
        raise ValueError("粗扫描未找到水位根 z1=%.4f" % z1)
    fa = _residual(a, sec1, sec2, z1, Q, n, L_m, kbhs, n2, brg, kwusan, kba01)
    fb = _residual(b, sec1, sec2, z1, Q, n, L_m, kbhs, n2, brg, kwusan, kba01)
    for _ in range(300):
        c = 0.5 * (a + b)
        fc = _residual(c, sec1, sec2, z1, Q, n, L_m, kbhs, n2, brg, kwusan, kba01)
        if abs(fc) < 1e-13 or (b - a) < 1e-12:
            return c
        if fa * fc <= 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


def compute(params):
    """执行逐段水面线推求（含可选的河段 n 值反推）。"""
    p = params
    kba01 = p["kba01"]
    kbhs = p["kbhs"]
    kwusan = p["kwusan"]
    ncp = p["ncp"]
    rivers = []
    for j in range(p["k1"]):
        nb = p["nb1"][j]
        secs = [{"bs": p["bs"][j][k], "hs": p["hs"][j][k],
                 "as_": p["as_"][j][k], "bss": p["bss"][j][k], "hss": p["hss"][j][k]}
                for k in range(nb)]
        reach_of = _reach_of_section(p["ncm1"][j], nb)
        reach_q = _reach_of_section(p["nq1"][j], nb)
        brg = p["brg"][j] if p["brg"] else [0.0] * nb
        threes = p["threes"]

        # ---- 逐断面上推水位（起始水位取 h(1,j)）----
        Z = [p["h"][j][0]]
        PROPS = [section_props(Z[0], secs[0], p["q"][j][reach_q[0]], kba01)]
        SEG = []                          # 各河段（i → i+1）的成果
        for i in range(nb - 1):
            r = reach_of[i]
            n = p["cmn"][j][r]
            n2s = p["cmn"][j][r] if kbhs == 2 else None
            Q = p["q"][j][reach_q[i]]
            L_km = p["dl"][j][i]
            L_m = L_km * 1000.0
            sec1, sec2 = secs[i], secs[i + 1]
            z2 = _solve_z2(sec1, sec2, Z[-1], Q, n, L_m, kbhs, n2s,
                           brg[i + 1], kwusan, kba01)
            h1 = section_props(Z[-1], sec1, Q, kba01)
            h2 = section_props(z2, sec2, Q, kba01)
            vh = (h1["V"] ** 2 - h2["V"] ** 2) / (2.0 * G)
            if kbhs == 2 and n2s is not None:
                i1, i2 = friction_slope(h1, n2s, Q), friction_slope(h2, n2s, Q)
            else:
                i1, i2 = friction_slope(h1, n, Q), friction_slope(h2, n, Q)
            dh = 0.5 * (i1 + i2) * L_m
            if kwusan == 1:
                dh += expansion_E(h1["V"], h2["V"]) * (h1["V"] - h2["V"]) ** 2 / (2.0 * G)
            SEG.append({"nc": n,
                        "ns": p["cmns"][j][min(reach_of[i], len(p["cmns"][j]) - 1)],
                        "dh": dh, "brg": brg[i + 1], "vh": vh,
                        "dl": L_km, "q": Q})
            Z.append(z2)
            PROPS.append(h2)

        # ---- 成果行：每行携带本断面的 H、V 与「本断面 → 下一断面」河段量 ----
        rows = []
        for i in range(nb):
            pr = PROPS[i]
            row = {"i": i + 1, "h": Z[i], "v": pr["V"],
                   "A": pr["A"], "B": pr["B"], "R": pr["R"], "P": pr["P"],
                   "nc": None, "ns": None, "dh": None, "brg": None,
                   "vh": None, "dl": None, "q": None}
            if i < nb - 1:
                row.update(SEG[i])
            rows.append(row)

        # ---- 可选：已知水位反推河段 n 值（式 5，ncp=1）----
        n_back = []
        if ncp == 1:
            for i in range(nb - 1):
                Q = SEG[i]["q"]
                L_m = p["dl"][j][i] * 1000.0
                h1, h2 = PROPS[i], PROPS[i + 1]
                A = 0.5 * (h1["A"] + h2["A"])
                R = 0.5 * (h1["R"] + h2["R"])
                ifz = ((Z[i + 1] - Z[i])
                       - (h1["V"] ** 2 - h2["V"] ** 2) / (2.0 * G)) / L_m
                n_back.append(A * R ** (2.0 / 3.0) * math.sqrt(abs(ifz)) / Q
                              if ifz > 0 else 0.0)
        rivers.append({"j": j + 1, "nb": nb, "reach_of": reach_of,
                       "断面": rows, "反推n值": n_back})
        if threes:
            for rw in rivers[-1]["断面"]:
                rw["h"] = rw["h"] * threes
    return {"程序": PROGRAM_ID, "输入": p, "河流": rivers}


# ============================================================
# 输出（复刻原著 .OUT 版式）
# ============================================================

def _f(x, w, d):
    """FORTRAN 风格数值：d 位小数、去掉前导 0，右对齐到宽度 w。"""
    s = "%.*f" % (d, x)
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return s.rjust(w)


def _fdot(x, w):
    """§ 末尾带小数点的整数值，右对齐到宽度 w（如 892.）。"""
    if abs(x - round(x)) < 1e-9:
        return ("%d." % int(round(x))).rjust(w)
    return ("%g" % x).rjust(w)


def _wrap(vals, per_line, first_prefix, next_prefix):
    """
    把字符串列表按 per_line 一行折行输出。
    first_prefix 用于首个输出行，next_prefix 用于后续行。
    """
    out = []
    for s in range(0, len(vals), per_line):
        chunk = vals[s:s + per_line]
        pre = first_prefix if s == 0 else next_prefix
        out.append(pre + "".join(chunk))
    return out or [first_prefix.rstrip()]


class _Out(object):
    """
    原著 .OUT 逐行输出：行与行之间不加空行，段落之间插入「单空格行」作为分隔
    （权威 D-15.OUT 共 360 行、行尾 CRLF，其中 13 行为单空格分隔行）。
    """

    def __init__(self):
        self.L = []

    def w(self, s=""):
        self.L.append(s)

    def sep(self, n=1):
        self.L.extend([" "] * n)

    def text(self):
        return "\n".join(self.L) + "\n"


def render(params, result, out_path=None):
    """生成原著风格文本计算书（.OUT）。"""
    p = params
    k1 = p["k1"]
    o = _Out()
    o.w((" 文件：" + (out_path or "D-15.out")).ljust(52))
    o.sep()
    o.w("                       河 道 水 面 线 计 算 D-15                        ")
    o.w("                       =========================")
    o.sep()
    o.sep()

    o.w("     河流数  河流中最大计算断面数  主槽最多n值河段数  最大流量河段数")
    o.w("       K1             K2                  K3                K4      ")
    o.w("%9d%15d%20d%18d" % (p["k1"], p["k2"], p["k3"], p["k4"]))
    o.sep()

    o.w("  主槽bs\\hs关系线最大结点数  滩地最多n值河段数  滩地bss\\hss关系线最大结点数")
    o.w("             K6                     K7                      K8             ")
    o.w("%15d%23d%23d" % (p["k6"], p["k7"], p["k8"]))
    o.sep()

    o.w("     允许迭代次数  计算方法  扩大或收缩")
    o.w("         mmax        kbhs      kwusan  ")
    o.w("%13d%12d%12d" % (p["mmax"], p["kbhs"], p["kwusan"]))
    o.sep()

    o.w("     断面资料用b\\h或a\\h曲线计算  是否已知水位求n值")
    o.w("                kba01                  ncp        ")
    o.w("%21d%18d" % (p["kba01"], p["ncp"]))
    o.sep()

    o.w("     允许水位    高程    予给河段  调整步长  变化as1值")
    o.w("      迭代差   换算系数   水位差     系数     的参数  ")
    o.w("       epsh     threes     dh0       as1        ds1   ")
    o.w("%11s%11s%10s%9s%11s" % (_f(p["epsh"], 11, 3), _f(p["threes"], 11, 1),
                                _f(p["dh0"], 10, 1), _f(p["as1"], 9, 4),
                                _f(p["ds1"], 11, 4)))
    o.sep()

    def ints_block(title, vals_list):
        o.w(title)
        for j in range(k1):
            o.w("%12d" % vals_list[j])

    ints_block("     不同流量的河段数 (NRQI1(J), J= 1, JN)", p["nrqil"])
    ints_block("     主槽不同n值河段数 (NRCM1(J), J= 1, JN)", p["nrcml"])
    ints_block("     滩地不同n值河段数 (NRCMS1(J), J= 1, JN)", p["nrcmsl"])
    ints_block("     河道计算断面数 (NB1(J), J= 1, JN)", p["nb1"])

    def idx_block(title, rows_of_ints):
        """界面序号数组：每行 10 个，行首字段宽 12、其后字段宽 7。"""
        o.w(title)
        for vals in rows_of_ints:
            for s in range(0, len(vals), 10):
                chunk = vals[s:s + 10]
                o.w("".join(("%12d" % v) if t == 0 else ("%7d" % v)
                             for t, v in enumerate(chunk)))
        return

    idx_block("     不同流量河段上游界面序号 (NQ1(K,1), K=1,NRQI)", p["nq1"])
    idx_block("     主槽不同n值河段上游界面序号 (NCM1(K,1), K=1, NRCM)", p["ncm1"])
    idx_block("     滩地不同n值河段上游界面序号 (NCMS1(K,1), K=1, NRCMS)", p["ncms1"])

    # 河段流量 / 主槽n值 / 滩地n值：原著按「主槽n值河段数(nrcml)」逐河段回显
    o.w("     河段流量 (Q(K,1), J= 1, JN)")
    for j in range(k1):
        nre = p["nrcml"][j]
        qr = _reach_of_section(p["nq1"][j], max(nre, 1))
        vals = [_fdot(p["q"][j][qr[min(k, max(nre, 1) - 1)]], 8) for k in range(nre)]
        for ln in _wrap(vals, 8, " " * 5, " " * 5):
            o.w(ln)
    o.w("     主槽河段n值 (CMN(K,1), J= 1, JN)")
    for j in range(k1):
        vals = [_f(v, 8, 4) for v in p["cmn"][j]]
        for ln in _wrap(vals, 8, " " * 5, " " * 5):
            o.w(ln)
    o.w("     滩地河段n值 (CMNS(K,1), J= 1, JN)")
    for j in range(k1):
        nre = p["nrcml"][j]
        ns = p["cmns"][j]
        vals = [_f(ns[min(k, len(ns) - 1)], 8, 4) for k in range(nre)]
        for ln in _wrap(vals, 8, " " * 5, " " * 5):
            o.w(ln)
    o.w("     主槽断面bs\\hs曲线结点数 (MNCS(K,1), K=1, NB)")
    for j in range(k1):
        for ln in _wrap([("%5d" % v) for v in p["mncs"][j]], 10, " " * 5, " " * 5):
            o.w(ln)
    o.w("     滩地断面bss\\hss曲线结点数 (MNCSS(K,1), K=1, NB)")
    for j in range(k1):
        for ln in _wrap([("%5d" % v) for v in p["mncss"][j]], 10, " " * 5, " " * 5):
            o.w(ln)
    o.w("     断面距离 (DL(K,1), K= 1, JN)")
    for j in range(k1):
        for ln in _wrap([_f(v, 8, 4) for v in p["dl"][j]], 8, " " * 5, " " * 5):
            o.w(ln)

    # 断面水位：首行缩进 5，续行 9，末段（不足一行）回 5 —— 依权威 OUT 版式
    o.w("     断面水位 ((H(I,J),I=1,%4d),J=1,%4d)" % (p["nb1"][0], k1))
    for j in range(k1):
        vals = [_f(v, 8, 2) for v in p["h"][j]]
        for t, s in enumerate(range(0, len(vals), 8)):
            o.w((" " * 5 if t % 2 == 0 else " " * 9) + "".join(vals[s:s + 8]))

    def sec_block(title, j, arrays, per_line=8):
        o.w(title)
        for k in range(p["nb1"][j]):
            vals = arrays[k]
            for t, s in enumerate(range(0, len(vals), per_line)):
                if t == 0:
                    o.w("%4d%13s" % (k + 1, vals[0]) + "".join(vals[1:per_line]))
                else:
                    # 原著续行缩进 5/9 交替（偶数行 5、奇数行 9）——依权威 OUT 版式
                    o.w((" " * 5 if t % 2 == 0 else " " * 9) + "".join(vals[s:s + per_line]))

    for j in range(k1):
        sec_block(" 主槽断面bs\\hs曲线结点宽度   I(BS(L,I,1), L=1, NCS)", j,
                  [[_f(v, 8, 2) for v in row] for row in p["bs"][j]])
        sec_block(" 主槽断面as\\hs曲线结点面积   I(AS(L,I,1), L=1, NCS)", j,
                  [[_f(v, 8, 2) for v in row] for row in p["as_"][j]])
        sec_block(" 主槽断面bs(as)\\hs曲线结点高程   I(HS(L,I,1), L=1, NCS)", j,
                  [[_f(v, 8, 2) for v in row] for row in p["hs"][j]])
        sec_block(" 滩地断面bss\\hss曲线结点宽度   I(BSS(L,I,1), L=1, NCS)", j,
                  [[_f(v, 8, 2) for v in row] for row in p["bss"][j]])
        # 滩地面积 ASS 的回显结点数与主槽面积 AS 相同（kba01=0 时为 1 个），
        # 而文件中 ASS 实存 mncss 个值 —— 依权威 OUT 反证（BSS/HSS 回显 mncss 个）
        sec_block(" 滩地断面ass\\hss曲线结点面积   I(ASS(L,I,1), L=1, NCS)", j,
                  [[_f(v, 8, 2) for v in p["ass"][j][k][:p["as_cnt"][j][k]]]
                   for k in range(p["nb1"][j])])
        sec_block(" 滩地断面bss(ass)\\hss曲线结点高程   I(HSS(L,I,1), L=1, NCS)", j,
                  [[_f(v, 8, 2) for v in row] for row in p["hss"][j]])

    o.sep()
    o.sep()
    o.sep()

    for rv in result["河流"]:
        o.w("     河流序号 J=%2d" % rv["j"])
        o.sep()
        o.w("   I      H      V     主槽   滩地   水头   桥闸   流速   断面      Q ")
        o.w("        水位   流速    ｎ值   ｎ值   损失   损失   水头   间距    流量")
        o.sep()
        rows = rv["断面"]
        for idx, rw in enumerate(rows):
            o.w("%4d%9.3f%8.3f" % (rw["i"], rw["h"], rw["v"]))
            if idx < len(rows) - 1:
                o.w(" " * 23 + "%5s%7s%7s%7s%7s%7s%8s" % (
                    _f(rw["nc"], 5, 4), _f(rw["ns"], 7, 4), _f(rw["dh"], 7, 2),
                    _f(rw["brg"], 7, 2), _f(rw["vh"], 7, 2), _f(rw["dl"], 7, 2),
                    _fdot(rw["q"], 8)))
        o.w("          计算结束  O.K.!")
    return o.text()


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 文件路径 | dict。"""
    params = parse(data)
    result = compute(params)
    text = render(params, result, out_path=out_txt)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
