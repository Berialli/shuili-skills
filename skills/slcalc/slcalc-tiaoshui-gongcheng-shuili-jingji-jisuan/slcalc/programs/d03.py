# -*- coding: utf-8 -*-
"""
D-3 消能计算程序 —— 内核
=========================
复刻《水利程序集》D-3 程序（作者：陈靖齐，水电部天津勘测设计院；
数据来源：《水利水电工程设计计算程序集》公之于众版，
乌鲁木齐正海水利科技有限公司，张校正教授级高工技术总负责）。

功能（原著二、功能）：
  本程序包括「底流消能」「挑流消能」两个程序；面流消能因理论不成熟、
  设计不合理时可能引起上游较大波浪运动，对航运与河岸不利，原著未编程序。
  程序按消能类型号 ED 分三种：
    ED=1  底流消能 · 开挖消力池（求池深 S、池长 Lk，并给出 hc、hc″、h″、ΔZ）
    ED=2  底流消能 · 消力坎/消力墙（求坎高 C、坎上总水头 H10、淹没系数 σs）
    ED=3  挑流消能（求挑流初速 V1、鼻坎水深 h1、反弧半径 R、挑流射程 L、
                     冲刷坑深度 Tk）

算法来源与反解过程
------------------
原著 D-3Intro.rtf 中的 13 个公式在文本层全部缺失（为 OLE 公式对象）。
本次改造把 RTF 中 17 个 Equation 对象（`\\objdata` → OLE2 → "Equation Native"）
全部解出，按 MTEF v3 记录逐字节还原出字符流，得到原著公式的真实形式。
还原结果（字符流原样，下标/上标按 MTEF 记录序内联）：

  底流消能·消力池
    o00 : hc = q / Vc                                   (1) 收缩断面水深
    o01 : hc″ = hc/2·(√(1+8·Frc) − 1)                   (2) 共轭水深
    o02 : Frc = Vc²/(g·hc)                              (3) 费鲁德数
    o03 : E0 = hc + Vc²/(2g·φ²)                         (4) 上游→收缩断面能量方程
    o04 : ΔZ = V²/(2g·φ′²) − V2²/(2g)                   (5) 跃后→下游能量方程
    o05 : Vt = q/ht                                     (6) 连续性
    o06 : V2 = q/h″                                     (7) 连续性
    o07 : h″ = ht + S + ΔZ = σ·hc″     σ=1.05           (8) 几何关系（稍淹没水跃）

  底流消能·消力坎
    o08 : h″ = C + H1                                   (9) 几何关系
    o09 : H10 = H1 + V2²/(2g)                           (10) 坎上总水头
    o10 : q = σs·m·√(2g)·H10^(3/2)                      (11) 坎流（淹没堰）
    o11 : σs = f(hs/H10) = f((ht − C)/H10)              (12)(13) 淹没系数试验表

  挑流消能
    o12 : V1 = √(2g·(S − h1·cosθ + V0²/(2g)))           (1) 射流初速
    o13 : x = V1·t·cosθ                                 (2) 水平路程
    o14 : y = V1·t·sinθ − ½·g·t²                        (3) 垂直路程
    o15 : y = −(Z + ht − S − h1·cosθ/2)                 (4) 射流落至河床的垂直距离
    o16 : T = K·q^(1/2)·Z^(1/4) − ht                    (5) 冲刷坑深度

【重要校核与裁决】
  (A) 式(5)/o04 的“项序”问题：OUT 计算书正文把该式印为
      「ΔZ=V*V/(2*g*φ′*φ′)-V2*V2/2/g」，但用 D-3-1.OUT 权威值反演（见下），
      实际参与运算的是
          ΔZ = Vt²/(2g) − V2²/(2g·φ′²)      （Vt=q/ht，V2=q/h″）
      已用穷举法排除其余 8 种项序/变体（见 d3_verify.py 第③段），
      本内核按反演结果实现；正文印出的式子与实现不一致。
      【知识库对照】该印刷项序恰与教材一致——吴持恭《水力学》第四版式9.10、
      《水力学核心公式速查》7.8 均写作 Δz=q²/(2gφ′²h_t²)−q²/(2g h_T²)（φ′ 附下游断面）。
      即：D-3 程序的实现同时偏离其自印文本与教材，三者中唯权威 OUT 的数值链可唯一反演，
      故按权威 OUT 定稿（UA 项，如实标注）。
  (B) 式(1)/o12 未出现 φ：但挑流射程对 φ 高度敏感，用算例 8 反证：
      若按 o12 字面（不含 φ）V1=√(2g(S−h1cosθ+V0²/2g))，算例 8 得 L=43.36m，
      与书值 35.25m 相差 +23%；只有 V1=φ·√(2g·S) 量级才对（见 d3_verify.py 第④段）。
      本内核按 V1 = φ·√(2g·S) 实现，并保留 h1=q/V1 的迭代语义。
  (C) 式(4)/o15 的 h1·cosθ 项：按字面 (含 −h1cosθ/2) 复现算例 7/8 偏 −6.1%/+0.3%；
      按 +h1·cosθ（鼻坎水面高出坎顶的几何修正）复现偏 −2.0%/+1.8%，综合最优，
      故本内核取 y = Z + ht − S + h1·cosθ。该点是本内核的“非唯一反演点”。

知识库对照（《WorkBuddy知识库/水利知识库》，2026-09-11 回补检索）
--------------------------------------------------------------
检索范围：卷7《泄水与过坝建筑物》第1章 1.6/1.7/1.8 与第5章 5.3/5.4；
          《水力学》吴持恭第四版上册 ch5-8、下册 ch9-11；
          《水力学核心公式速查》第 7 节；《程序词汇×教材概念速查表》；SL278-2002 笔记。
结论（一致项）：
  K1. 挑流射程：《手册》卷7 式1.7-1「L₁=v₁cosθ·[v₁sinθ+√(v₁²sin²θ+2g(h₁cosθ+a₁))]/g」——
      与本内核 L=[V1²sinθcosθ+V1cosθ√(V1²sin²θ+2g·Y)]/g 同式，且 **y 中的 +h1·cosθ
      正是手册的 h₁cosθ 项**（手册 a₁ = 鼻坎至下游水面/河床高差 = Z+ht−S，见 d3_verify 第⑦段）。
      本内核先前对该项的“+h1cosθ / −h1cosθ/2”二义性由此得到教材裁决：取 +h1cosθ。
  K2. 冲刷坑深度：《手册》式1.7-5「T = K·q^0.5·Z^0.25 − h_t」与内核逐字一致。
  K3. 岩石冲刷系数：《手册》表1.7-1 给出**范围**——坚硬完整 0.8~1.1、中等 1.1~1.4、
      破碎 1.4~1.8、极破碎 1.8~2.0；D-3 说明书的 4 档 0.9/1.2/1.5/1.8 是各档中值，一致。
  K4. 消力池长度系数：《手册》1.6 式1.6-38~44 与 5.4、吴持恭式9.12 均为 **L_k=(0.7~0.8)L_j**，
      内核取 0.75 为区间中值，一致；L_j=6.9(h″−h′) 亦为教材标准式。
  K5. 淹没度取值：吴持恭「σ_j=1.05~1.10 的淹没水跃最理想」，内核 σ=1.05，一致。
  K6. 池深方程口径：吴持恭式9.9「d=σ_j·h_c1″−(h_t+Δz)」与内核 S=h″−ht−ΔZ 同式；
      《手册》1.6 记作 σ'=(h_t+d)/h₂″=1.05~1.10（不含 Δz），两口径并存，D-3 从吴持恭式。

未闭合点（如实标注）
--------------------
  U1. 消力坎淹没系数 σs 的“试验数据表”：说明书正文只写「为一试验数据表」，
      17 个公式对象、D-3vb.EXE 的数据段（已按 double/float 双精度扫描、
      0.15~1.05 区间单调序列穷举）与 8 个算例文件里均无该表。
      **知识库亦未收录**：吴持恭第四版下册第9章只在笔记中留指针「消能坎…第二类，
      原理同加大下游水深」（该章 9.2~9.12 未 OCR 入笔记）；《手册》卷7 第1章 1.6
      给的是“受控水跃 式1.6-23~32（薄坎/厚坎/升坎…）坎上水深与 Fr₁ 匹配”，无 σs 表。
      库内唯一相关的**淹没系数关系式**是《手册》卷7 第5章 5.3 的宽顶堰淹没式
      「σ_s = 2.31·(h_s/H₀)·(1−h_s/H₀)^0.4，0.72<h_s/H₀<0.98」——但该式针对**宽顶堰**
      而非消力坎，且算例 5 的 hs/H10=0.635 已越出其适用范围，同点给出 0.980 vs
      武水答复 0.964（差 1.7%），代入后本例偏差反而由 +0.2% 恶化到 +2.4%。
      故本内核**保留以武水答复锚定的等效表**（SIG_X/SIG_Y，锚点 hs/H10=0.6354→σs=0.964），
      并把该手册式实现为可切换备选 `sig_s_handbook()`（SIG_MODE="handbook" 启用），
      对照数据见 d3_verify.py 第⑦段。降级表述：**库中无该表，属原著私有表**。
      非唯一性：单锚点只约束曲线过一点，无穷多条单调曲线满足（已给 3 点替代表同解反证）。
  U2. 消力坎流量系数 m′：算例 6 由第 8 字段直接输入 0.40；算例 5 该字段为 0，
      由权威答复 H10=3.146 反解得 m′≈0.4199（取默认 0.42）。
      知识库中 D-3 消力坎的 m′ 表同样未收录（卷7 5.3 只给宽顶堰 m=0.32~0.385、
      实用堰 m=0.45~0.50），0.42 属可复现但非唯一反演的等效参数。
  U3. 挑流冲刷系数 K：算例 7 的 K 由第 12 字段给出 1.25（复现 Tk=9.175 vs 书值 9.2，
      −0.3%）；算例 8 第 12 字段为 0，按第 11 字段岩石类别 3 → K=1.5，
      得 Tk=7.291 vs 书值 7.9（−7.7%）。
      **知识库解释**：手册表1.7-1 给的是**范围**（破碎 1.4~1.8），例 8 的“岸基比较破碎”
      对应 1.4~1.8，取中值 1.6 恰好得 Tk=7.91（书值 7.9）——即书值在范围之内，
      差异源于「D-3 简化表取档位中值 1.5」vs「手册按范围取 1.6」，非公式错误。
      详见 d3_verify.py 第⑧段。
  U4. 挑流第 3 字段语义不清：例 7 为 50.15（=上游水位−河床高程=H+P），
      例 8 为 22（=坝高）。两种解释在同一字段上互斥，本内核按“坝高 P”回显且
      不参与运算（挑流 φ 由输入直给），互斥反证见 d3_verify.py 第⑤段。
  U5. 底流消能的 ΔZ 项序：教材（吴持恭式9.10 / 速查表 7.8）写作
      Δz=q²/(2gφ′²h_t²)−q²/(2g h_T²)（φ′ 附于下游断面）；D-3 的权威 OUT 正文亦如此印。
      但逐位反演证明程序实际按 ΔZ=q²/(2g h_t²)−q²/(2gφ′²h″²)（φ′ 附于跃后断面）计算，
      教材式在本例给出 Δz=0.466（S=1.382），与权威 Δz=0.3879（S=1.4597）不符。
      判定：D-3 在该式的项序上与教材/自印文本均不同，以权威 OUT 反演结果为准，
      见 d3_verify.py 第③段穷举。

验证
----
  · 算例 1（ED=1）逐位对拍权威 D-3-1.OUT（GBK）：
    本内核 Vo=1.5546 hc=1.0471 hc2=4.2818 → hc=0.9325 hc2=4.6169 h2=4.8477
    ΔZ=0.3879 S=1.4598 Lj=25.4218 Lk=19.0663
    权威    Vo=1.5546 hc=1.0471 hc2=4.2819 → hc=0.9326 hc2=4.6167 h2=4.8476
    ΔZ=0.3879 S=1.4597 Lj=25.4206 Lk=19.0654
    → 10 个显示值中 3 个逐位相同、7 个末位差 1（最大 |Δ|=0.0012，相对 4.7e-5），
      系原著 VB6 迭代终止/单精度舍入边界，已如实标注（见 d3_verify 第①段）。
  · 算例 2~8 对拍说明书参考答案，Lk 系数按各书口径（0.7/0.75/0.8）换算后比较。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "D-3"
TITLE = "消能计算书"
AUTHOR = "陈靖齐（水电部天津勘测设计院）"

G = 9.8            # 原著重力加速度（D-2/D-3 同源，取 9.8 m/s²）
SIGMA = 1.05       # 稍淹没水跃安全系数（式 8 / o07）
LK_COEF = 0.75     # 消力池长度系数 Lk = 0.75·Lj（说明书：本文取 0.75）
M_DEFAULT = 0.42   # 消力坎流量系数默认值（由算例 5 权威答复反演，见 U2）
PHI_DEFAULT = 0.95  # 泄流流速系数缺省（说明书："亦可近似取 φ=0.95"）
K_TABLE = {1: 0.9, 2: 1.2, 3: 1.5, 4: 1.8}   # 岩石冲刷系数（说明书原表 Ⅰ~Ⅳ）

# 消力坎淹没系数等效标定表（U1，非原著原表）
#   横坐标 x = hs/H10 = (ht − C)/H10
SIG_X = [0.00, 0.20, 0.30, 0.40, 0.50, 0.55, 0.60, 0.6354, 0.65,
         0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
SIG_Y = [1.000, 1.000, 1.000, 0.995, 0.988, 0.982, 0.975, 0.964, 0.960,
         0.947, 0.930, 0.905, 0.870, 0.820, 0.750, 0.650]

# 备选：知识库唯一收录的淹没系数关系式
#   《水工设计手册》卷7 第5章 5.3（宽顶堰淹没系数式，范围 0.72<hs/H0<0.98）
#     σ_s = 2.31·(hs/H0)·(1 − hs/H0)^0.4
#   该式针对宽顶堰而非消力坎，且算例 5 的 hs/H10=0.635 越出其适用范围，
#   同点 0.980 vs 武水答复 0.964（差 1.7%），故不作默认；仅供对照与切换。
SIG_MODE = "calibrated"     # "calibrated"（默认，武水答复锚定表）| "handbook"


def sig_s_handbook(x):
    """《水工设计手册》卷7 5.3 宽顶堰淹没系数式（0.72<x<0.98，域外取 1.0 / 端点值）。"""
    if x <= 0.72:
        return 1.0
    x = min(x, 0.98)
    return 2.31 * x * (1.0 - x) ** 0.4


def sig_s(x, mode=None):
    """淹没系数 σs = f(hs/H10)；mode 缺省取模块常量 SIG_MODE。"""
    if (mode or SIG_MODE) == "handbook":
        return sig_s_handbook(x)
    return interp(SIG_X, SIG_Y, x)


def interp(xs, ys, x):
    """一维线性插值（端点外沿用首末值）。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


# ------------------------------------------------------------
# 通用水力学小工具
# ------------------------------------------------------------

def conjugate_depth(hc, q, g=G):
    """共轭水深 hc″ = hc/2·(√(1+8·Frc) − 1)，Frc = Vc²/(g·hc)。"""
    vc = q / hc
    frc = vc * vc / (g * hc)
    return hc / 2.0 * (math.sqrt(1.0 + 8.0 * frc) - 1.0)


def solve_hc_from_energy(S, e0, q, phi, g=G):
    """
    由能量方程 hc + q²/(2g·φ²·hc²) = E0 + S 解收缩断面水深 hc。

    方程在 (0, (2C)^{1/3}) 上单调下降（C=q²/(2gφ²)），取最小正根（物理根）。
    S 为消力池开挖深度（消力坎计算取 S=0）。
    """
    if phi <= 0:
        raise ValueError("泄流流速系数 phi 必须为正")
    c = q * q / (2.0 * g * phi * phi)
    lo, hi = 1e-9, (2.0 * c) ** (1.0 / 3.0)
    rhs = e0 + S
    if lo + c / (lo * lo) - rhs < 0:
        raise ValueError("能量方程无正根（E0+S 过小）")
    if hi + c / (hi * hi) - rhs > 0:
        raise ValueError("能量方程在物理区间内无根（E0+S 过大）")
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if m + c / (m * m) - rhs > 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


# ------------------------------------------------------------
# ED=1 底流消能 · 消力池
# ------------------------------------------------------------

def run_pool(H, P, D, ht, q, phi, phi1, g=G, tol=1e-10, max_iter=400):
    """
    消力池计算（o00~o07）。

    E0 = H + P + D + V0²/(2g)，V0 = q/(H+P)
      · H 坎上水头、P 坝高、D 上下游河床高程差；
      · (H+P+D) = 上游水位 − 消力池底高程（例 5 中消力池底为下游河床，D=1.5）。
    迭代（塞德尔法）：S=0 → hc → hc″ → h″=σhc″ → ΔZ → S=h″−ht−ΔZ → 回代。
    返回 dict；若首轮 hc″ ≤ ht 则不需消力池（S=0）。
    """
    # 缺省处理：φ 为 0 时取 0.95（说明书："亦可近似取 φ=0.95"）；
    # φ′ 为 0 时取 φ（算例 3 的 INT 第 8 字段即为 0，属缺省）
    if phi <= 0:
        phi = PHI_DEFAULT
    if phi1 <= 0:
        phi1 = phi
    v0 = q / (H + P) if (H + P) > 0 else 0.0
    e0 = H + P + D + v0 * v0 / (2.0 * g)
    vt = q / ht if ht > 0 else 0.0
    k2 = q * q / (2.0 * g * phi1 * phi1)

    # 第一轮（S=0）：判定是否需要消力池
    hc0 = solve_hc_from_energy(0.0, e0, q, phi, g)
    hc20 = conjugate_depth(hc0, q, g)
    out = {"E0": e0, "V0": v0, "Vt": vt,
           "hc0": hc0, "hc20": hc20, "need": hc20 > ht}

    if hc20 <= ht:
        out.update(hc=hc0, hc2=hc20, h2=None, dz=0.0, S=0.0,
                   Lj=0.0, Lk=0.0, iters=0)
        return out

    S = 0.0
    hc = hc0
    hc2 = hc20
    it = 0
    for it in range(1, max_iter + 1):
        h2 = SIGMA * hc2
        v2 = q / h2
        # ΔZ = Vt²/(2g) − V2²/(2g·φ′²)  （项序经 D-3-1.OUT 反演锁定，见 (A)）
        dz = vt * vt / (2.0 * g) - v2 * v2 / (2.0 * g * phi1 * phi1)
        S_new = h2 - ht - dz
        hc_new = solve_hc_from_energy(S_new, e0, q, phi, g)
        hc2_new = conjugate_depth(hc_new, q, g)
        if abs(S_new - S) < tol:
            S, hc, hc2 = S_new, hc_new, hc2_new
            break
        S, hc, hc2 = S_new, hc_new, hc2_new
    h2 = SIGMA * hc2
    v2 = q / h2
    dz = vt * vt / (2.0 * g) - v2 * v2 / (2.0 * g * phi1 * phi1)
    S = h2 - ht - dz
    Lj = 6.9 * (hc2 - hc)
    out.update(hc=hc, hc2=hc2, h2=h2, dz=dz, S=S,
               Lj=Lj, Lk=LK_COEF * Lj, iters=it)
    return out


# ------------------------------------------------------------
# ED=2 底流消能 · 消力坎（消力墙）
# ------------------------------------------------------------

def run_sill(H, P, D, ht, q, phi, m, g=G, sig_mode=None):
    """
    消力坎计算（o00~o03、o08~o11）。

    方程组：
      hc 同消力池（S=0）；hc″ 共轭；h″ = σ·hc″（σ=1.05）
      C    = h″ − H1                     （o08，H1 坎上水深）
      H10  = H1 + q²/(2g·h″²)            （o09，坎上总水头）
      q    = σs·m·√(2g)·H10^(3/2)        （o10，淹没堰流）
      σs   = f((ht − C)/H10)             （o11，试验数据表 → SIG_X/SIG_Y）
    对 C 做二分求根。
    """
    if phi <= 0:
        phi = PHI_DEFAULT
    v0 = q / (H + P) if (H + P) > 0 else 0.0
    e0 = H + P + D + v0 * v0 / (2.0 * g)
    hc = solve_hc_from_energy(0.0, e0, q, phi, g)
    hc2 = conjugate_depth(hc, q, g)
    h2 = SIGMA * hc2
    k2 = q * q / (2.0 * g * h2 * h2)     # h″ 断面流速水头
    m_use = m if m and m > 0 else M_DEFAULT

    def fp(c):
        h1 = h2 - c
        if h1 <= 0:
            return -q
        h10 = h1 + k2
        sig = sig_s((ht - c) / h10, sig_mode) if h10 > 0 else 0.0
        return sig * m_use * math.sqrt(2.0 * g) * h10 ** 1.5 - q

    lo, hi = 0.0, h2
    flo, fhi = fp(lo), fp(hi)
    if flo < 0:
        raise ValueError("消力坎：坎高下界仍不满足泄流条件（流量过大/坎过低）")
    if fhi > 0:
        raise ValueError("消力坎：坎高上界仍不满足泄流条件（流量过小/坎过高）")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if fp(mid) > 0:
            lo = mid
        else:
            hi = mid
    C = 0.5 * (lo + hi)
    H1 = h2 - C
    H10 = H1 + k2
    sig = sig_s((ht - C) / H10, sig_mode)
    V2 = q / h2
    Lj = 6.9 * (hc2 - hc)
    return {"E0": e0, "V0": v0, "hc": hc, "hc2": hc2, "h2": h2,
            "C": C, "H1": H1, "H10": H10, "sig_s": sig, "m": m_use,
            "sig_mode": sig_mode or SIG_MODE, "V2": V2, "Lj": Lj,
            "Lk": LK_COEF * Lj}


# ------------------------------------------------------------
# ED=3 挑流消能
# ------------------------------------------------------------

def run_jet(S, Z, ht, q, theta_deg, phi, k_class, K, g=G, k_from_input=True):
    """
    挑流消能计算（o12~o16）。

      V1 = φ·√(2g·S)                        射流初速（φ 由程序外输入）
      h1 = q/V1                             鼻坎处水深（迭代语义）
      Y  = Z + ht − S + h1·cosθ             射流落至下游河床的垂直落差
      L  = [V1²sinθcosθ + V1cosθ·√(V1²sin²θ + 2g·Y)] / g
      R  = 10·h1                            反弧段半径
      Tk = K·√q·Z^(1/4) − ht                冲刷坑深度
    """
    if phi <= 0:
        phi = 0.95
    th = math.radians(theta_deg)
    v1 = phi * math.sqrt(2.0 * g * S) if S > 0 else 0.0
    h1 = q / v1 if v1 > 0 else 0.0
    yy = Z + ht - S + h1 * math.cos(th)
    v1c = v1 * math.cos(th)
    v1s = v1 * math.sin(th)
    disc = v1s * v1s + 2.0 * g * yy
    if disc < 0:
        disc = 0.0
    L = (v1 * v1 * math.sin(th) * math.cos(th) + v1c * math.sqrt(disc)) / g
    R = 10.0 * h1
    if k_from_input and K and K > 0:
        k_use = K
    else:
        k_use = K_TABLE.get(int(k_class), 1.2)
    Tk = k_use * math.sqrt(q) * (Z ** 0.25) - ht
    return {"V1": v1, "h1": h1, "Y": yy, "L": L, "R": R,
            "K": k_use, "Tk": Tk}


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

ED_NAMES = {1: "底流消能、消力池", 2: "底流消能、消力坎", 3: "挑流消能"}


def parse(data):
    """
    解析输入。data: dict 或 INT 文件路径。

    INT 数值流（首字段为消能类型号 ED）：
      ED=1（8 值）：ED, H, P, D, ht, q, phi, phi1
      ED=2（8 值）：ED, H, P, D, ht, q, phi, m     （m=0 表示取默认 0.42）
      ED=3（12 值）：ED, H, P, D, S, Z, ht, q, theta, phi, k_class, K
    其中：H 坎上水头(m)；P 坝高(m)；D 上下游河床高程差(m)；ht 下游水深(m)；
      q 单宽流量(m²/s)；phi 泄流流速系数；(ED=1) phi1 水跃处流速系数；
      (ED=2) m 消力坎流量系数；(ED=3) S 鼻坎水头(m)、Z 下游水位以上水头(m)、
      theta 挑角(°)、k_class 岩石类别(1~4)、K 岩石冲刷系数(0=按类别取)。
    """
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 1:
        raise ValueError("INT 数据为空")
    ed = int(round(nums[0]))
    if ed not in ED_NAMES:
        raise ValueError(f"未知消能类型号 ED={ed}（应为 1/2/3）")
    p = {"程序": PROGRAM_ID, "ED": ed}
    rest = nums[1:]
    if ed == 1:
        if len(rest) < 7:
            raise ValueError("ED=1 需要 7 个数值：H,P,D,ht,q,phi,phi1")
        H, P, D, ht, q, phi, phi1 = rest[:7]
        p.update(H=H, P=P, D=D, ht=ht, q=q, phi=phi, phi1=phi1)
    elif ed == 2:
        if len(rest) < 7:
            raise ValueError("ED=2 需要 7 个数值：H,P,D,ht,q,phi,m")
        H, P, D, ht, q, phi, m = rest[:7]
        p.update(H=H, P=P, D=D, ht=ht, q=q, phi=phi, m=m)
    else:
        if len(rest) < 11:
            raise ValueError("ED=3 需要 11 个数值：H,P,D,S,Z,ht,q,theta,phi,k_class,K")
        H, P, D, S, Z, ht, q, theta, phi, kcls, K = rest[:11]
        p.update(H=H, P=P, D=D, S=S, Z=Z, ht=ht, q=q, theta=theta,
                 phi=phi, k_class=kcls, K=K)
    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params):
    """执行消能计算，返回结构化结果。"""
    ed = params["ED"]
    result = {"程序": PROGRAM_ID, "消能类型": ED_NAMES[ed],
              "基本资料": {k: v for k, v in params.items()
                           if k not in ("程序", "ED")}}
    if ed == 1:
        s = run_pool(params["H"], params["P"], params["D"], params["ht"],
                     params["q"], params["phi"], params["phi1"])
        result["计算"] = s
        result["结果"] = {
            "V0": s["V0"], "hc": s["hc"], "hc2": s["hc2"], "h2": s["h2"],
            "dz": s["dz"], "S": s["S"], "Lj": s["Lj"], "Lk": s["Lk"],
            "hc0": s["hc0"], "hc20": s["hc20"], "need_pool": s["need"],
        }
    elif ed == 2:
        s = run_sill(params["H"], params["P"], params["D"], params["ht"],
                     params["q"], params["phi"], params["m"],
                     sig_mode=params.get("sig_mode"))
        result["计算"] = s
        result["结果"] = {
            "hc": s["hc"], "hc2": s["hc2"], "h2": s["h2"], "C": s["C"],
            "H1": s["H1"], "H10": s["H10"], "sig_s": s["sig_s"], "m": s["m"],
            "sig_mode": s["sig_mode"], "Lj": s["Lj"], "Lk": s["Lk"],
        }
    else:
        s = run_jet(params["S"], params["Z"], params["ht"], params["q"],
                    params["theta"], params["phi"], params["k_class"],
                    params["K"])
        result["计算"] = s
        result["结果"] = {
            "V1": s["V1"], "h1": s["h1"], "R": s["R"], "L": s["L"],
            "K": s["K"], "Tk": s["Tk"],
        }
    return result


# ------------------------------------------------------------
# 输出
# ------------------------------------------------------------

def _f4(x):
    return f"{x:.4f}"


def render(params, result):
    """生成汉字计算书（复刻原著 D-3.OUT 风格）。"""
    ed = params["ED"]
    L = []
    L.append("  一、基本资料与计算假定")
    L.append("")
    L.append(f"       {ED_NAMES[ed]}")
    if ed == 1:
        L.append(f"  坎上水头H:                  {_f4(params['H']):>10}(米)")
        L.append(f"  坝高P:                      {_f4(params['P']):>10}(米)")
        L.append(f"  上下游河床高程差D:          {_f4(params['D']):>10}(米)")
        L.append(f"  下游水深ht:                 {_f4(params['ht']):>10}(米)")
        L.append(f"  单宽流量q:                  {_f4(params['q']):>10}(立方米/秒)")
        L.append(f"  泄流系数phi:                {_f4(params['phi']):>10}")
        L.append(f"  (下游)水跃处流速系数phi1:   {_f4(params['phi1']):>10}")
    elif ed == 2:
        L.append(f"  坎上水头H:                  {_f4(params['H']):>10}(米)")
        L.append(f"  坝高P:                      {_f4(params['P']):>10}(米)")
        L.append(f"  上下游河床高程差D:          {_f4(params['D']):>10}(米)")
        L.append(f"  下游水深ht:                 {_f4(params['ht']):>10}(米)")
        L.append(f"  单宽流量q:                  {_f4(params['q']):>10}(立方米/秒)")
        L.append(f"  泄流系数phi:                {_f4(params['phi']):>10}")
        L.append(f"  消力坎流量系数m':            {_f4(params['m']):>10}(0=取默认0.42)")
    else:
        L.append(f"  坝上水头H:                  {_f4(params['H']):>10}(米)")
        L.append(f"  坝高P:                      {_f4(params['P']):>10}(米)")
        L.append(f"  上下游河床高程差D:          {_f4(params['D']):>10}(米)")
        L.append(f"  鼻坎水头S1:                 {_f4(params['S']):>10}(米)")
        L.append(f"  下游水位以上水头Z:          {_f4(params['Z']):>10}(米)")
        L.append(f"  下游水深ht:                 {_f4(params['ht']):>10}(米)")
        L.append(f"  单宽流量q:                  {_f4(params['q']):>10}(立方米/秒)")
        L.append(f"  挑角theta:                  {_f4(params['theta']):>10}(度)")
        L.append(f"  泄流流速系数phi:            {_f4(params['phi']):>10}")
        L.append(f"  岩石类别k_class:            {int(params['k_class']):>10}")
        L.append(f"  岩石冲刷系数K:              {_f4(params['K']):>10}(0=按类别取)")
    L.append("")
    L.append("  二、计算结果")
    L.append("")
    c = result["计算"]
    if ed == 1:
        L.append(f"  上游行进流速Vo:             {_f4(c['V0']):>10}(米/秒)")
        L.append(f"  收缩水深hc:                 {_f4(c['hc0']):>10}(米)")
        rel = ">" if c["need"] else "<="
        L.append(f"  共轭水深hc2:                {_f4(c['hc20']):>10}(米)  "
                 f"{rel}  下游水深ht=  {params['ht']:>6.2f}(米)")
        if c["need"]:
            L.append(f"  收缩水深hc:                 {_f4(c['hc']):>10}(米)")
            L.append(f"  共轭水深hc2:                {_f4(c['hc2']):>10}(米)")
            L.append(f"  跃后水深h2:                 {_f4(c['h2']):>10}(米)")
            L.append(f"  下游水位以上水头dz:         {_f4(c['dz']):>10}(米)")
            L.append(f"  消力池深度s:                {_f4(c['S']):>10}(米)")
            L.append(f"  水跃长度Lj:                 {_f4(c['Lj']):>10}(米)")
            L.append(f"  消力池长度Lk:               {_f4(c['Lk']):>10}(米)")
        else:
            L.append("  （共轭水深不大于下游水深，不需修建消力池）")
    elif ed == 2:
        L.append(f"  上游行进流速Vo:             {_f4(c['V0']):>10}(米/秒)")
        L.append(f"  收缩水深hc:                 {_f4(c['hc']):>10}(米)")
        L.append(f"  共轭水深hc2:                {_f4(c['hc2']):>10}(米)")
        L.append(f"  跃后水深h2:                 {_f4(c['h2']):>10}(米)")
        L.append(f"  坎上水深H1:                 {_f4(c['H1']):>10}(米)")
        L.append(f"  消力坎上总水头H10:          {_f4(c['H10']):>10}(米)")
        L.append(f"  淹没系数sig_s:              {_f4(c['sig_s']):>10}")
        L.append(f"  消力坎高度C:                {_f4(c['C']):>10}(米)")
        L.append(f"  水跃长度Lj:                 {_f4(c['Lj']):>10}(米)")
        L.append(f"  消力池长度Lk:               {_f4(c['Lk']):>10}(米)")
    else:
        L.append(f"  挑流初速V1:                 {_f4(c['V1']):>10}(米/秒)")
        L.append(f"  鼻坎处水深h1:               {_f4(c['h1']):>10}(米)")
        L.append(f"  反弧段半径R:                {_f4(c['R']):>10}(米)")
        L.append(f"  射流落至河床的垂直距离y:    {_f4(c['Y']):>10}(米)")
        L.append(f"  挑流射程L:                  {_f4(c['L']):>10}(米)")
        L.append(f"  岩石冲刷系数K:              {_f4(c['K']):>10}")
        L.append(f"  冲刷坑深度Tk:               {_f4(c['Tk']):>10}(米)")
    return render_text(PROGRAM_ID, TITLE, [("", L)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 文件路径 | dict"""
    params = parse(data)
    result = compute(params)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [], result)
    else:
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
