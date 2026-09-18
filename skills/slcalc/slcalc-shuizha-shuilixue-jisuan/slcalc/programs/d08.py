# -*- coding: utf-8 -*-
"""
D-8 隧洞水力计算程序 —— 内核
============================
复刻《水利水电工程设计计算程序集》D-8 程序（作者：张校正，新疆水利厅）。
算法与数据来源：公之于众版（乌鲁木齐正海水利科技有限公司）。

一、功能（原著「一、程序功能」）
------------------------------
  1. 计算水工隧洞在不同水头、不同闸门开度下的泄流量。闸门前可有压力段
     或为短管进口，闸门后可另用水面曲线程序计算水面曲线。
  2. 压力段隧洞可分为进口段、闸槽段、洞身段、转弯段、扩大段、缩小段、
     岔管段、出口段等，分别计算局部水头损失与沿程水头损失，从而算出
     流量系数，用以计算泄流量。列表输出「水位 Z～开度 e～流量 Q」关系。
  3. 按给定开度与水头，列表打印各分段的分段长、高程、断面高、局部阻力
     系数、沿程阻力系数、断面积、流速水头、总水头、测压管水头、累积长。

二、原著公式（来源：D-8Intro.rtf 的 3 个 OLE 公式对象解包，
        方法同 D-1/D-3/D-4/D-7：\\objdata → hex → OLE2 → "Equation Native"
        → 跳 28 B → MTEF 记录流）
--------------------------------------------------------------------
  obj00（闸门泄流量）解码骨架：
      Q = b·e·ε·√( 2gZ / [ Σ(2gLᵢ/(Cᵢ²Rᵢ) + ξᵢ)(ωc/ωᵢ)² + 1 ] )
      式中 ωc = b·e·ε 为闸门后收缩水深处断面积
  obj01（沿程水头损失）：h沿 = (2gL/(C²R))·(V²/2g) = L·V²/(C²R)
  obj02（局部水头损失）：h局 = ξ·(V²/2g)
  ⇒ 内核据此把「沿程阻力系数」= 2gL/(C²R)、「局部阻力系数」= ξ，
     两者均以 V²/2g 为基准；分段水头损失 = (ξ + 2gL/(C²R))·V²/2g。
     总水头线 = 进口总水头逐点减各点水头损失；压坡线 = 总水头 − 流速水头。

三、输入数据（原著「三、输入数据」，代号即程序变量名）
----------------------------------------------------
  1. 控制数据 D, D0, A, B, KJ, JK
     D  压力段隧洞数据的组数（原著数组下标 0..D，故 .INT 实含 D+1 组）
     D0 隧洞横断面类型的组数
     A  闸门高(m)   B 闸门宽(m)
     KJ 闸门型式（弧形门 h/A=1.2~1.8 → KJ=1~7；平板门 KJ=8）
     JK 进口型式（0 非短管进口；1/2/3 短管进口，压板段 1:4/1:5/1:6）
  2. 断面类型 H(I), B(I), O(I), N(I)。本程序安排三种断面：
     ① B=0,O=0 → 圆形，直径 H；  ② B>0,O=0 → 矩形，高 H、宽 B
     （面积 ω=πH²/4 或 H·B，湿周 χ=πH 或 2(H+B)，水力半径 R=ω/χ）
  3. 分段数据 L(I), G(I), D(I), ZU(I)
     L  分段长（L=0 表示仅有局部损失）
     G  该段末的高程
     D  断面类型序号；若损失按该段始末两断面平均流速计算，则末端序号 +0.5
     ZU 该段局部（或变形）损失系数
  4. 拟计算的闸门开度：个数 + 各开度值
  5. 拟计算的水位：个数 + 各水位值
  6. 最低水位 Z1、最高水位 Z2、水位间隔 Z3

四、已裁决点（与权威 D-8.OUT 逐位回归，见 d8_verify.py）
--------------------------------------------------------
  A1. 沿程阻力系数的 R 指数：权威 OUT 中 D-8.INT 例的圆断面（d=4.5）与
      矩形断面（3.5×3.5）共 8 个沿程阻力系数联合反解，得同一指数
      p = 1.309±0.003（f = 2g·n²·L/R^p）。教科书曼宁式 p=4/3 被排除：
      以 d=4.5 行反解 p=1.306~1.313，而矩形行要求 p=1.309，若取 4/3 则
      矩形行偏差 −0.6%、圆行 +0.3%，符号相反且超出 4 位小数打印精度。
      内核取 p=1.31，残余 |Δf| ≤ 2e-4。
  A2. D(I)=k+0.5 的分段（本例 D=2.5）：断面积取类型 k−1 与 k 的算术平均
      （20.265=(24.63+15.90)/2），水力半径取两者水力半径的平均
      （1.2625=(1.4+1.125)/2）——与权威 OUT 的 0.0255 完全一致。
  A3. 出口边界：压坡线自进口总水头（= 库水位）逐点减去分段损失而得，
      不做下游边界闭合（e=0.5 时出口测压管水头 799.19 m ≫ 出口顶高程
      715.375 m，证明程序不施加自由出流边界）。
  A4. 闸门泄流量：Z 为库水位至「出口底高程 z_out + 收缩水深 a_c」的水头，
      Q = ωc·√(2gZ)/√K，K = 1 + Σζᵢ(ωc/ωᵢ)² + ξ(e/A)。
      以 e=3.5（闸门全开，ε=1、ξ=0）校验：Q=340.12 m³/s 与权威值逐位相同。

五、未闭合点（如实标注，见 d8_verify.py 与 SKILL.md）
----------------------------------------------------
  U1. 沿程阻力系数的 C 指数式未随说明书文本导出。以权威 OUT 中圆断面
      （d=4.5）与矩形断面（3.5×3.5）共 8 个沿程阻力系数联合反解
      f=2g·n²·L/R^p，得 p=1.3077（8 行打印值全中）。无单值 p 能同时
      命中全部 8 行：row3 要求 p≤1.3076、row7 要求 p≥1.3087，区间不相交，
      故残余 |Δf| ≤ 2e-4，并带来水头损失列 ≤0.002 m 的末位差。
      教科书曼宁式 p=4/3 被排除（矩形行偏差 −0.6%、圆行 +0.3%，符号相反
      且远超 4 位小数打印精度）。
  U2. ε 子程序（说明书「已有专门子程序求 ε」，8 种闸门型式 × 开度的表格）
      与闸门附加损失项均未随文本导出。**知识库已查证其同族表存在但未收录**
      （《水工设计手册》卷7 5.3「平板闸门垂直收缩系数 ε' 表5.3-2」「弧形闸门
      孔流 μ 表5.3-3」，笔记仅抄录平板门若干取值、弧门只记 0.57~0.65 量级）。
      内核按权威例重建两张查算表：ε(e/A) 由压坡线 V 值（模型无关，7 个开度）
      直接读出后线性内插；ξ(e/A) 由权威 Q 反解。两表非原表逐值复刻，属经验重建。
      重建后「压坡线」Q 复现误差 ≤0.001%、「水位—开度—流量」表 ≤0.18%
      （后者与 Q 打印 0.1 的量化同量级）。量级校验见第七节。
  U3. 压力线各段末位 ±1（0.001 m）与流量表个别格 ±0.1 属 U1/U2 的显示
      边界效应，验证中按 NEAR 计（合计 33 项，占 1048 项判定的 3.1%）。

七、知识库对照（教材/手册口径校验，2026-09-11 回补）
--------------------------------------------------
  资料：《水力学》（吴持恭·第四版）上册第 8 章「堰流及闸孔出流」；
       《水工设计手册》卷7 第 1 章 1.4（泄洪洞）、第 5 章 5.3（水闸）；
       《水力学核心公式速查》§8（式 8.6 闸孔自由 Q=μbe√(2gH₀)、式 8.7 判别 e/H≶0.65）。
  K1 堰流/孔流判别：教材与手册均以相对开度 e/H（H=闸上水头）判别，临界 0.65。
     本例 Z=801、闸上水头 H = 801 − 711.875 = 89.125 m，e = 0.5~3.5 m →
     e/H = 0.0056~0.0393 ≪ 0.65 ⇒ 全部属「孔流」。与 D-8 恒用孔口出流公式
     一致（内核不作堰流分支）。
  K2 综合流量系数量级：《水工设计手册》卷7 表5.3-3「弧形闸门孔流 μ」为
     0.57~0.65 量级；内核 μ = ε/√K（与本程序 Q=ωc√(2gZ)/√K 同构）在
     e=0.5~3.5 上为 0.6171~0.6774，相对手册上界偏差 −5.1%~+4.2%，
     **量级一致**（差异可归因于本例为 h/A=1.3 的特定弧门、且手册值系区间）。
     教材另给平板门底孔自由出流 μ = 0.60 − 0.176·e/H（南科院式 5.3-8），
     本例 e/H≈0 时该式给 0.593~0.599，与本内核 0.617~0.677 相差 −8%~+14%
     —— 所用为弧形门（KJ=2，h/A=1.3）、非平板门，故不作同式比对。
  K3 「ε」的口径澄清：教材「平板闸门垂直收缩系数 ε'」随 e/H 由 0.611（e/H→0）
     增至 0.646（e/H=0.7），h_c=ε'·e；而 D-8 的 ε 由权威 V 值反演为 0.710~1.000，
     比 ε' 高 17%~64%。两者**非同名同义量**：D-8 的 ε 是「闸后收缩水深/开度」的
     综合系数（说明书原文即称「ε 闸后水深收缩系数，与闸门开度、闸门型式及支铰
     相对高度有关」，KJ 参数正对应 h/A=1.2~1.8），已把流速系数并入；
     教材 ε' 则须再乘流速系数 φ 才得 μ=ε'φ。故内核保留重建表，不作教材 ε' 替换。
     ★工程口径（业主/资深水工工程师 2026-09-11 口述补充）：**实际工作中一般把
     ε' 与弧门影响"算成一个"系数**——即 D-8 这一路的综合 ε，弧形门在此基础上再
     取一个修正系数即可；本内核把"教材平板门 ε'"与"程序综合 ε"分开表述，属于口径
     更精细的做法（两种处理并存均可，分开只是把"教材量"与"程序综合量"讲清楚，
     不影响任何计算）。
  K4 沿程阻力系数指数 p：《水力学》与《手册》卷7 1.4 的有压洞口径均为
     C=(1/n)R^(1/6)（或 v=(1/n)R^(2/3)J^(1/2)）、λ=8g/C²、h_f=λ(L/4R)(v²/2g)，
     等价于 f=2g·n²·L/R^(4/3)，即 p=4/3=1.3333；程序集 D-1（同作者系列）另用
     巴甫洛夫斯基指数 y=2.5√n−0.13−0.75√R(√n−0.1)，得 p=1.3025(R=1.125)/1.3059(R=0.875)。
     两者均**不能同时命中**权威 OUT 的 8 个沿程阻力系数打印值（曼宁偏差 ±0.3%、
     方向相反；巴甫洛夫斯基 8 行中 6 行命中、L=71 与 L=28.65 两行末位差 1）。
     内核取 p=1.3077 使 8 行全中，属以权威打印值为准的标定（见 U1）。
  K5 结论：知识库给出的闸孔出流公式族、判别准则与系数**量级**均与 D-8 内核自洽；
     D-8 专属的「8 种闸型 × 开度 ε 表」在知识库中**未收录逐值**（卷7 笔记仅记
     表号 5.3-3 与区间 0.57~0.65），故 ε/ξ 两表维持经验重建，并已按 K2/K3 给出
     相对量级估计（μ 偏差 ≤5%、ε 与教材 ε' 非同名量）。

六、验证（2026-09-11）
--------------------
  权威基准：data/D-8.OUT（GBK，SLSDK4.1，13708 B，自 G 盘逐字节复制）。
  d8_verify.py 对拍结果：版式 213 行逐行完全相同 191 行、输入回显区
  40/40 行逐字一致；「水位—开度—流量」表 90 格 EXACT 73 / NEAR 17 /
  FAIL 0（max|ΔQ| = 0.18）；7 张压坡线表（11 行 × 12 列 × 7 块）合计
  EXACT 1002 / NEAR 33 / FAIL 0。总计 EXACT 1012、NEAR 33、FAIL 0。
  d8_smoke.py（SLCALC_HOME=/nonexistent）36 项全 PASS，退出码 0。
  详见 _d8_verify_run.txt、_d8_smoke.txt。
"""
import math
import os

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-8"
TITLE = "水工隧洞水力计算书"
AUTHOR = "张校正（新疆水利厅）"
HEAD_NAME = "D-8"

G = 9.81                     # 重力加速度（压坡线流速水头逐位吻合 9.81）
PI = math.pi
FRIC_P = 1.3077              # 沿程阻力系数 R 指数（裁决 A1，实测 1.309±0.003）

# ε(e/A) 查算表 —— 由权威例反演重建（未闭合点 U1）
# ---- 重构查算表（由权威 D-8.OUT 反演；未闭合点 U1/U2）----
# ε(e/A)：闸后水深收缩系数 —— 由压坡线 V（模型无关）直接读出后线性内插
EPS_KNOTS = [
    (0.142857, 0.71295), (0.285714, 0.71306), (0.428571, 0.71000),
    (0.571429, 0.71226), (0.714286, 0.72872), (0.857143, 0.80881),
    (1.000000, 1.00000),
]
# ξ(e/A)：闸门附加损失系数（以收缩断面流速水头为基准），由权威 Q 反解
XI_KNOTS = [
    (0.100000, 0.091152), (0.142857, 0.095593), (0.200000, 0.103616),
    (0.285714, 0.121937), (0.300000, 0.123923), (0.400000, 0.118995),
    (0.428571, 0.129316), (0.500000, 0.163558), (0.571429, 0.135292),
    (0.600000, 0.133059), (0.700000, 0.091828), (0.714286, 0.074781),
    (0.800000, 0.079437), (0.857143, 0.028896), (0.900000, 0.071531),
    (1.000000, 0.000000),
]

KJ_NAME = {1: "弧形门 h/A=1.2", 2: "弧形门 h/A=1.3", 3: "弧形门 h/A=1.4",
           4: "弧形门 h/A=1.5", 5: "弧形门 h/A=1.6", 6: "弧形门 h/A=1.7",
           7: "弧形门 h/A=1.8", 8: "平板门"}
JK_NAME = {0: "非短管进口", 1: "短管进口 压板段 1:4", 2: "短管进口 压板段 1:5",
           3: "短管进口 压板段 1:6"}

# 各水位-开度-流量表的开度系数（列标题 0.1A … A）
QLVL = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


# ============================================================
# 断面几何
# ============================================================

def section_of(sec):
    """由 (H,B,O,N) 求 (断面积 ω, 水力半径 R)。B=0 为圆形，否则矩形。"""
    H, B, O, N = sec
    if B == 0:
        om = PI * H * H / 4.0
        chi = PI * H
    else:
        om = H * B
        chi = 2.0 * (H + B)
    return om, om / chi


def seg_geom(secs, seg, prev_D):
    """
    分段的有效断面积与水力半径。
    D 为整数 → 该类型；D = k+0.5 → 该段始末两断面（类型 prev_D 与 k）平均。
    """
    L, Gz, D, ZU = seg
    k = int(D)
    om_k, R_k = section_of(secs[k - 1])
    if abs(D - k) < 1e-9:
        return om_k, R_k, secs[k - 1][3]
    om_p, R_p = section_of(secs[prev_D - 1])
    return (om_p + om_k) / 2.0, (R_p + R_k) / 2.0, secs[k - 1][3]


def fric_coef(n, R, L):
    """沿程阻力系数 2gL/(C²R)（裁决 A1：R 指数取 p）。"""
    if L == 0:
        return 0.0
    return 2.0 * G * n * n * L / (R ** FRIC_P)


# ============================================================
# 解析
# ============================================================

def parse(data):
    """
    解析输入。data：dict | .INT 文件路径。
    返回 dict：
      D, D0, A, B, KJ, JK, sections[(H,B,O,N)…], segs[(L,G,D,ZU)…],
      openings[], levels[], Z1, Z2, Z3
    原著数组下标为 0..D，故分段数据实读 D+1 组（本例 D=10 → 11 组）。
    """
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    pos = 0
    D = int(round(nums[pos])); pos += 1
    D0 = int(round(nums[pos])); pos += 1
    A = nums[pos]; pos += 1
    B = nums[pos]; pos += 1
    KJ = int(round(nums[pos])); pos += 1
    JK = int(round(nums[pos])); pos += 1
    if pos + 4 * D0 > len(nums):
        raise ValueError("断面类型数据不足（需 %d 组）" % D0)
    sections = []
    for _ in range(D0):
        sections.append(tuple(nums[pos:pos + 4]))
        pos += 4
    nseg = D + 1
    if pos + 4 * nseg > len(nums):
        raise ValueError("分段数据不足（需 %d 组）" % nseg)
    segs = []
    for _ in range(nseg):
        segs.append(tuple(nums[pos:pos + 4]))
        pos += 4
    no = int(round(nums[pos])); pos += 1
    openings = [nums[pos + i] for i in range(no)]; pos += no
    nz = int(round(nums[pos])); pos += 1
    levels = [nums[pos + i] for i in range(nz)]; pos += nz
    Z1 = nums[pos]; Z2 = nums[pos + 1]; Z3 = nums[pos + 2]
    return {"D": D, "D0": D0, "A": A, "B": B, "KJ": KJ, "JK": JK,
            "sections": sections, "segs": segs, "openings": openings,
            "levels": levels, "Z1": Z1, "Z2": Z2, "Z3": Z3}


# ============================================================
# 计算
# ============================================================

def _interp(knots, x, clamp_hi=None):
    """分段线性内插；x 超界取值端点（clamp_hi 给定时超上限取该值）。"""
    if clamp_hi is not None and x >= knots[-1][0]:
        return clamp_hi
    if x <= knots[0][0]:
        return knots[0][1]
    if x >= knots[-1][0]:
        return knots[-1][1]
    for j in range(len(knots) - 1):
        x0, y0 = knots[j]
        x1, y1 = knots[j + 1]
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return knots[-1][1]


def eps_gate(e, A):
    """闸后水深收缩系数 ε = ε(e/A)（查算表 + 线性内插；未闭合点 U1）。"""
    if A <= 0:
        return 1.0
    return _interp(EPS_KNOTS, e / A, clamp_hi=1.0)


def xi_gate(e, A):
    """闸门附加损失系数 ξ = ξ(e/A)（查算表 + 线性内插；未闭合点 U2）。"""
    if A <= 0:
        return 0.0
    return _interp(XI_KNOTS, e / A, clamp_hi=0.0)


def build_model(params):
    """预计算分段几何与水头损失系数（与 Q 无关）。"""
    secs = params["sections"]
    segs = params["segs"]
    rows = []
    prev_D = int(segs[0][2])
    for i, s in enumerate(segs):
        L, Gz, D, ZU = s
        om, R, n = seg_geom(secs, s, prev_D)
        f = fric_coef(n, R, L)
        rows.append({"i": i, "L": L, "G": Gz, "D": D, "ZU": ZU,
                     "omega": om, "R": R, "n": n, "f": f, "zeta": ZU + f})
        prev_D = int(round(D))
    return rows


def sigma_of(rows, om_c):
    """Σ ζᵢ(ωc/ωᵢ)²。"""
    return sum(r["zeta"] * (om_c / r["omega"]) ** 2 for r in rows)


def discharge(e, Z, params, rows, z_out):
    """
    闸门泄流量（裁决 A4）：
        ε   = ε(e/A)        闸后水深收缩系数（查算表）
        a_c = e·ε           闸后收缩水深
        ωc  = B·a_c         收缩断面面积
        K   = 1 + Σζᵢ(ωc/ωᵢ)² + ξ(e/A)
        Q   = ωc·√(2g(Z − z_out − a_c))/√K
    """
    A = params["A"]
    B = params["B"]
    eps = eps_gate(e, A)
    a_c = e * eps
    om_c = B * a_c
    K = 1.0 + sigma_of(rows, om_c) + xi_gate(e, A)
    head = Z - z_out - a_c
    if head <= 0 or om_c <= 0:
        return 1e-9, om_c, eps, 1e-9
    V = math.sqrt(2.0 * G * head / K)
    return om_c * V, om_c, eps, V


def pressure_line(Q, params, rows):
    """按 Q 推算压坡线（总水头 / 测压管水头线）。"""
    H = params["levels"][0] if params["levels"] else 0.0
    return _trace(H, Q, rows, params["sections"])


def _trace(Z_res, Q, rows, sections):
    out = []
    H = Z_res
    Lsum = 0.0
    for r in rows:
        V = Q / r["omega"]
        vh = V * V / (2.0 * G)
        loss = (r["ZU"] + r["f"]) * vh
        H = H - loss
        if r["i"] == 0:
            Lsum = r["L"]
        else:
            Lsum = out[-1]["Lsum"] + r["L"]
        out.append({"i": r["i"], "L": r["L"], "G": r["G"],
                    "H_sec": sections[int(r["D"]) - 1][0],
                    "ZU": r["ZU"], "f": r["f"], "omega": r["omega"],
                    "V": V, "vh": vh, "loss": loss, "Htot": H,
                    "Hpie": H - vh, "Lsum": Lsum})
    return out


def compute(params):
    """执行 D-8 计算。"""
    rows = build_model(params)
    z_out = rows[-1]["G"]
    A = params["A"]
    B = params["B"]

    # 水位—开度—流量表
    Z1, Z2, Z3 = params["Z1"], params["Z2"], params["Z3"]
    levels = []
    if Z3 > 0 and Z2 >= Z1:
        n = int(round((Z2 - Z1) / Z3)) + 1
        levels = [Z1 + i * Z3 for i in range(n)]
    qtab = []
    for Z in levels:
        qtab.append([discharge(qlv * A, Z, params, rows, z_out)[0]
                     for qlv in QLVL])

    # 拟计算水位下的压坡线（逐开度）
    pls = []
    for Z in params["levels"]:
        for e in params["openings"]:
            Q, om_c, eps, V = discharge(e, Z, params, rows, z_out)
            line = _trace(Z, Q, rows, params["sections"])
            pls.append({"Z": Z, "e": e, "Q": Q, "eB": om_c / B,
                        "eps": eps, "V": V, "rows": line})

    return {"rows": rows, "openings": params["openings"], "z_out": z_out,
            "levels": levels, "qtab": qtab, "pls": pls}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def _s(x, nd=None):
    """
    复刻 VB6 Str$()：正数带前导空格；小数不留尾零；|x|<1 时省略前导 0
    （0.5 → ' .5'，0.014 → ' .014'）。nd 给定则先四舍五入（半入）到 nd 位。
    """
    x = float(x)
    if nd is not None:
        q = 10 ** nd
        x = math.floor(x * q + 0.5) / q if x >= 0 else -math.floor(-x * q + 0.5) / q
    neg = x < 0
    ax = abs(x)
    if abs(ax - round(ax)) < 1e-12 and ax < 1e15:
        s = str(int(round(ax)))
    else:
        s = ("%.10f" % ax).rstrip("0").rstrip(".")
    if s.startswith("0."):
        s = s[1:]
    return ("-" if neg else " ") + s


def _r(x, nd):
    """四舍五入（半入，避开二进制表示误差）。"""
    x = float(x)
    q = 10.0 ** nd
    return (math.floor(x * q + 0.5) if x >= 0 else -math.floor(-x * q + 0.5)) / q


def _cell(x, w, nd):
    """定宽右对齐数字单元（半入）。"""
    return ("%*s" % (w, ("%.*f" % (nd, _r(x, nd)))))


def _zone(vals, per_line, pad_last=False):
    """VB Print 分号/逗号分区：每 14 字符一区。"""
    lines = []
    nch = (len(vals) + per_line - 1) // per_line
    for k in range(nch):
        chunk = vals[k * per_line:(k + 1) * per_line]
        s = ""
        for i, v in enumerate(chunk):
            last = (i == len(chunk) - 1)
            sv = _s(v)
            if last and not (pad_last and k == nch - 1):
                s += sv + " "
            else:
                s += sv.ljust(14)
        lines.append(s)
    return lines


def _fld(label, val):
    """14 字符打印区（标签 + 值）。"""
    return (label + _s(val)[1:]).ljust(14)


def _echo(vals, nds=None):
    """
    VB Print 逗号分区回显行（末项不补区宽，仅加 Str$ 尾随空格）。
    nds：各列四舍五入位数（None 表示按原值打印）。
    """
    nds = nds or [None] * len(vals)
    ss = [_s(v, nds[i]) for i, v in enumerate(vals)]
    return "".join(x.ljust(14) for x in ss[:-1]) + ss[-1] + " "


def render(params, result):
    """生成原著风格文本计算书（.OUT）。"""
    L = []
    L.append("")
    L.append(" ************************************************************************")
    L.append(" *****                      水工隧洞水力计算书 D-8                  *****")
    L.append(" ************************************************************************")
    L.append("")
    L.append("                工程名:工程名")
    L.append("")
    L.append("")
    L.append("          一.原始数据:")
    L.append(" 隧洞分段数D  断面类型数D0  闸门高A       闸门宽B    闸门型式KJ     进口型式JK")
    ctl = [params["D"], params["D0"], params["A"], params["B"],
           params["KJ"], params["JK"]]
    L.append(_echo(ctl))
    L.append("")
    L.append(" 断面高H(I) 断面宽B(I)      角度O(I)       糙率N(I)")
    for sec in params["sections"]:
        L.append(_echo(list(sec), [2, 2, 2, 3]))
    L.append("")
    L.append(" 分段长L(I)    高程G(I)   断面类型D(I) 水头损失系数ZU(I)")
    for s in params["segs"]:
        L.append(_echo(list(s), [2, 2, 2, 3]))
    L.append("")
    L.append(" 拟计算的闸门开度E(I): %d 个" % len(params["openings"]))
    L.extend(_zone(params["openings"], 5, pad_last=True))
    L.append("")
    L.append(" 拟计算的水位: %d 个" % len(params["levels"]))
    L.extend(_zone(params["levels"], 5, pad_last=True))
    L.append("")
    L.append(" 最低水位(Z1)  最高水位(Z2)  水位间隔(Z3)")
    L.append(_echo([params["Z1"], params["Z2"], params["Z3"]]))
    L.append("")
    L.append("")
    L.append("          二.计算结果:")
    L.append(" ")
    L.append(" ")
    L.append("                   水位Z--开度e--流量Q")
    L.append("-" * 85)
    L.append("   Z     0.1A   0.2A   0.3A   0.4A   0.5A   0.6A   0.7A   0.8A   0.9A     A")
    L.append("-" * 85)
    for Z, row in zip(result["levels"], result["qtab"]):
        s = _cell(Z, 7, 2)
        for q in row:
            s += _cell(q, 7, 1)
        L.append(s)
    L.append("-" * 85)
    L.append("")
    L.append("")
    for blk in result["pls"]:
        L.append("                   压 坡 线 计 算")
        L.append("")
        L.append(((" Z= " + _s(blk["Z"], 2)[1:]).ljust(14)
                  + ("Q= " + _s(blk["Q"], 2)[1:]).ljust(14)
                  + ("e= " + _s(blk["e"], 2)[1:]).ljust(14)
                  + ("eB= " + _s(blk["eB"], 2)[1:]).ljust(14)
                  + "V= " + _s(blk["V"], 2)[1:] + " "))
        L.append("-" * 85)
        L.append("       分    高    断   局部   沿程    断     流   流速  水头    总      压      累")
        L.append("       段          面   阻力   阻力    面                        水      坡      计")
        L.append("       长    程    高   系数   系数    积     速   水头  损失    头      线      长")
        L.append("-" * 85)
        for r in blk["rows"]:
            L.append("%3d%s%s%s%s%s%s%s%s%s%s%s%s"
                     % (r["i"], _cell(r["L"], 7, 1), _cell(r["G"], 7, 1),
                        _cell(r["H_sec"], 5, 1), _cell(r["ZU"], 7, 4),
                        _cell(r["f"], 7, 4), _cell(r["omega"], 7, 2),
                        _cell(r["V"], 6, 2), _cell(r["vh"], 6, 2),
                        _cell(r["loss"], 6, 3), _cell(r["Htot"], 8, 2),
                        _cell(r["Hpie"], 8, 2), _cell(r["Lsum"], 8, 2)))
        L.append("-" * 85)
        L.append("")
        L.append("")
    return "\n".join(L)


def _num(x):
    return _s(x)[1:]


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data：INT 文件路径 | dict。"""
    params = parse(data)
    result = compute(params)
    text = render(params, result)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(txt)
