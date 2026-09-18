# -*- coding: utf-8 -*-
"""
D-25 底流式消能工计算程序 —— 内核
=================================
复刻《水利程序集》D-25 程序（作者：卢礼标；
数据来源：《水利水电工程设计计算程序集》公之于众版，
乌鲁木齐正海水利科技有限公司，张校正教授级高工技术总负责）。

功能（原著一、程序功能）
------------------------
底流式消能工水力计算程序，适用于**挖深式、消力坎式和综合式**各类底流消能
结构；平面型式可为**等宽消力池**或**扩散消力池**。程序按消能结构类型号 JX
与平面型式号 PK 分派：

  JX=1 挖深式消力池   —— 由能量方程与水跃方程联立试算池深 S（隐函数）
  JX=2 消力坎式消力池 —— 求坎高 C（用消力坎壅水形成底流消能池，S=0）
  JX=3 综合式消力池   —— 池深 S 由输入给定，求坎高 C
  PK=1 等宽消力池；PK=2 扩散消力池（须给末端宽度 Bk，输出侧墙扩散角 θ）

程序具有连续多级消能功能：当消力坎下游仍为远离式水跃时，自动增设下一级
（第 2 级）消力坎，直到下游衔接为淹没式为止（本内核复刻该判据，见下）。

原著公式（D-25Intro.rtf 文本层 + D-25vb.EXE 内 BSTR 字面量逐字互证）
-------------------------------------------------------------------
说明书的公式在 RTF 文本层以 OLE 公式对象给出；本次改造除读 RTF/说明书的
公式行外，另从 D-25vb.EXE 的 UTF-16（VB6 BSTR）字面量中取出程序**自身
打印的公式串**（这是最权威的“程序自述口径”），两路互证，逐字得到：

  1. 收缩水深        E01 = Hc + q*q/(2*g*f1^2*Hc^2)                     [图文一致]
      等效迭代式（吴持恭《水力学》第四版下册式 9.5）
                     h_c(i+1) = q / [φ√(2g(E01 − h_ci))]
  2. 共轭水深
     A 等宽          Hc″ = Hc*(SQR(1+8*Fr^2)-1)/2                        [图文一致]
      式中 Fr = Frc = Vc/SQR(g*Hc)，Vc = q/Hc
     B 扩散          (Hc^2-Hc″^2)*(B+Bk) = 4*Q*Q*(1/Bk/Hc″-1/B/Hc)/g     [图文一致]
  3. 水跃长度
     A 等宽          Lj = 6.9*(Hc″-Hc)                                   [图文一致]
     B 扩散          Fr<6 时 Lj=(1+0.6*Fr)*Hc″；Fr>=6 时 Lj=4.6*Hc″；
                     Fr>17 时 Lj=B*L/(B+.1*L*tgθ)，
                     L = 10.3*Hc*(SQR(Fr)-1)^.81                        [图文一致]
  4. 消力池长度      Lk = (0.7~0.8)*Lj，程序采用 0.75*Lj                 [图文一致]
  5. 挖深式池深      S = 1.05*Hc″ − Ht
                       − q1*q1*(1/(f2*Ht)^2 − 1/(1.05*Hc″)^2)/(2*g)     [EXE BSTR 原文]
  6. 消力坎高度      C = 1.05*Hc″ + q1*q1/(2*g)/(1.05*Hc″)^2
                       − (q1/(4.43*p*m))^(2/3) − S                       [EXE BSTR 原文]
      式中 p 为消力坎淹没系数（自由溢流时 p=1）；m 为消力坎流量系数；
          S 为综合式池深（消力坎式 S=0）；4.43 = √(2g) 原著取数。
  7. 扩散消力池侧墙扩散角（由几何关系，输出项）
                     θ = arctg[ (Bk − B) / (2*Lk) ]

知识库对照（D:\\WorkBuddy知识库\\水利知识库\\，2026-09-11 检索）
-------------------------------------------------------------
母本：`02_水利教材精读\\水力学与工程水文学\\水力学_吴持恭第四版_下册_ch9-11.md`
      第 9 章「泄水建筑物下游的水流衔接与消能」§9.1 底流消能；
      `02_...\\水力学核心公式速查.md` §7「水跃与消能（★D-25）」；
      `01_水工设计手册精读\\卷7_泄水与过坝建筑物.md` §1.6「底流消能」。

  K1 【一致】收缩断面能量方程 E₀ = h_c + q²/(2gφ²h_c²)：吴持恭式 9.4、速查表 7.4。
     内核式 1 与之同式（本程序把 E₀ 记作 E01，并以池底为基准）。
  K2 【一致】迭代式 h_c(i+1)=q/[φ√(2g(E₀−h_ci))]，初值 q/(φ√(2gE₀))：
     吴持恭式 9.5、速查表 7.5。内核用同一方程的单调根二分，等价。
  K3 【一致】矩形共轭水深 h_c″ = (h_c/2)(√(1+8Fr²)−1)：
     吴持恭式 9.7 邻域、速查表 7.2、手册卷7 §1「水跃共轭水深」。
  K4 【一致】水跃长度 Lj = 6.9(h″−h′)：速查表 7.10、手册卷7 §1.6
     「L_j=6.9(h₂−h₁)」逐字一致。
  K5 【一致】池长系数 L_k = (0.7~0.8)L_j：吴持恭式 9.12、速查表 7.9、
     手册卷7 式 1.6-38~44。程序中取 0.75（区间中值），与 D-3 同一口径。
  K6 【一致】淹没系数 σ_j = 1.05~1.10 最理想 → 内核取 1.05：
     吴持恭「工程要求 σ_j=1.05~1.10 的淹没水跃最理想」、速查表 7.7。
  K7 【一致·改编】挖深式池深：吴持恭式 9.9「d = σ_j·h_c1″ − (h_t + Δz)」，
     式 9.10「Δz = q²/(2gφ′²h_t²) − q²/(2g h_T²)」，h_T = σ_j·h_c1″（式 9.7），
     φ′≈0.95。把 q1=q、φ′=f2、h_T=1.05hc″ 代入即得原著的
     S = 1.05hc″ − h_t − q1²[1/(f2 h_t)² − 1/(1.05hc″)²]/(2g)，**逐字等价**
     （原著把 φ′ 附在下游断面，与教材式 9.10 的项序完全一致）。
     另：教材式 9.9~9.11 明确指出「护坦降低 d 后 E₀→E₀+d，h_c1 须用 E₀+d
     重算 ⟹ d 与 h_c1 为隐函数关系，须试算」——本内核 JX=1 即按此迭代，
     且迭代解与权威 D-25-1.OUT 逐位吻合（见 验证）。
  K8 【一致】池长与池深的实测口径与手册卷7 §1.6 一致：手册记池深为
     σ′=(h_t+d)/h₂″=1.05~1.10（不含 Δz 的等价口径），与吴持恭式 9.9
     两口径并存；D-25 从吴持恭式（含 Δz），与权威 OUT 反演结果一致。
  K9 【库中无】消力坎（消能坎）水力计算：教材第 9 章把「消能坎（护坦末端建坎
     壅水）形成消能池」列为"第二类，原理同加大下游水深"，**该节未 OCR 入库**；
     手册卷7 §1.6 只给「受控水跃 式 1.6-23~32（薄坎/厚坎/升坎/斜坡坎/跌坎）」
     的坎上水深与 Fr₁ 匹配关系，无"消力坎高度 C"的成形式。
     故原著 C 式（K6 的 6.）在库中**无直接出处**，其正确性以原著自述 + EXE 字面量
     + 算例复现为准（算例 5/6 逐位复现，见 验证）。
  K10【库中无】消力坎淹没系数 p 的取值表：库内唯一的**淹没系数关系式**是
     手册卷7 §5.3 宽顶堰式 σ_s = 2.31(h_s/H₀)(1−h_s/H₀)^0.4（0.72<h_s/H₀<0.98）；
     该式针对**宽顶堰**而非消力坎，且算例 4 的 h_s/H₁₀=0.860 处给出 0.906，
     与本内核标定值 0.781 相差 16%，代入后算例 4 偏差由 −0.0001 恶化到
     +0.16 m，故不作默认（实现为可切换备选 `p_handbook_weir()`）。
     实际程序使用的 p(x) 试验表在本项目可及的全部来源（说明书正文、RTF、
     D-25vb.EXE 字符串常量与数值常量扫描、知识库）中**均未出现**，属原著私有表。

未闭合点（如实标注 + 反证）
--------------------------
  U1. **消力坎淹没系数 p(x) 表不可唯一反演**（x = h_s/H₁₀ = (Ht−C)/H₁₀）。
      · 由说明书算例反演得 4 个锚点：
          (x=0.37299, p=1.00015) 算例5(扩散消力坎)
          (x=0.36026, p=0.99998) 算例6 第1级
          (x=0.63392, p=0.96838) 算例3
          (x=0.83528, p=0.81936) 算例6 第2级  →  连同
          (x=0.86029, p=0.78150) 算例4
        即 p=1（自由溢流）只在小 x（≲0.4）成立；算例 3/4/第2级均属**淹没**
        （p<1）。这**反证了"程序恒取 p=1"的假设**：若取 p=1，算例 3 的 C
        为 1.5798（说明书画 1.5117，差 +0.068）、算例 4 的 C 为 1.4369
        （说明书画 0.8874，差 +0.550），远超打印精度。
      · 但 5 个锚点不足以唯一确定曲线：给出**两条不同**的单调表（本内核的
        P_X/P_Y 与 d3_verify 式的 d03 标定表）都能在 1e-3 内通过算例 3/4，
        故 p(x) 不可唯一反演。内核默认 `p_mode="auto"`（标定表），
        并保留 `p_mode="free"`（强制 p=1，对照用）与显式输入 p 的通道。
        **声明：标定表为"过锚点"的等效表，非原著原表。**
  U2. 综合式（JX=3）的 h_c″ 口径：权威 OUT 未覆盖 JX=3，仅说明书算例 4 给
     Hc″=5.1736；该值**只能**由 E01 = E0 + S = 13.1 复现（用 E01=E0 得 5.0285），
     故内核取 E01 = E0 + S（与挖深式同口径，物理上一致）。
  U3. 多级消能的触发与第 2 级能量口径：说明书只给"消力坎下游发生远离式水跃"
     需建下一级，未给判据式与第 2 级上游能量口径。本内核由算例 6 反演得
     **E0(2) = h_c″(1)**（用 1.05·h_c″(1)=5.2799 时 h_c=1.1733≠说明书画 1.2181；
     用 h_c″(1)=5.0285 得 h_c=1.2180 ✓），q(2)=Q/B、Ht(2)=Ht；判据取
     "h_c″(i+1) > Ht ⟹ 仍需第 i+1 级"。据此算例 6 恰好 2 级、算例 5 恰好 1 级，
     与说明书文字完全一致（反证见 d25_verify.py 第③段）。
  U4. 侧墙扩散角 θ 的算式：说明书只给"θ 是在算出池长后才能确定"，未给式。
     由算例 2（θ=8.3723°）与算例 5（θ=6.6772°）双点反演唯一锁定
     θ = arctg[(Bk−B)/(2·Lk)]，且**弧度→度换算常数取 π=3.14**（原著年代口径）：
     取 π=3.14 时两算例分别得 8.372270°/6.677248°（打印 8.3723/6.6772，双点全中）；
     取 π=3.14159265 时得 8.368025°/6.673862°（双点各差 0.004°）。
     备选式 atan[(Bk−B)/Lk]、asin[(Bk−B)/(2Lk)] 亦被排除（见 d25_verify.py 第④段）。
  U5. 扩散消力池共轭水深式为超越方程（Hc″ 无法显式），本内核用单调二分求解；
     算例 2 复现 Hc″=4.1882（说明书 4.1882）→ 唯一根。
  U6. 重力加速度取 g=9.8：由权威 D-25-1.OUT 的 Fr=3.7895 反演唯一锁定
     （g=9.81 给出 3.7875，与权威差 2 个末位；g=9.8 给出 3.7895 逐位同）。
  U7. 池深试算（JX=1）的收敛判定：说明书写"d 与 h_c1 为隐函数关系，须试算"，
     未给停止容差。由权威 D-25-1.OUT 的 7 个显示值联合反演，**唯一**锁定
     停止判据为 |ΔS| < 0.001 m（S 以 m 计，即 1 mm）：取该容差时
     hc=0.9405 / Fr=3.7895 / Hc″=4.5919 / Lj=25.1948 / S=1.3430 /
     Lk=18.8961 **七值全部逐位命中**；取 ≤1e-4（含严格收敛）时
     Fr→3.7896、Lj→25.1950、S→1.3431、Lk→18.8963（4 个末位各差 1）。
     已同步核对算例 2（扩散挖深式）同样命中（见 d25_verify.py 第①段）。
     注意：容差只影响"最后一位打印值"，不改变算法结构；本内核按 1e-3 实现以逐位复刻。

验证（详见 d25_verify.py / _d25_verify_run.txt）
----------------------------------------------
  · 权威 D-25-1.OUT（GBK，2045 B，54 行；首行「文件：L:\\01\\4.1版\\
    SLSDK4.1\\use\\D-25-1.out」）：逐行对拍 53 行，版式与数值字段全部逐字相同
    （首行运行期路径回显不生成——机器相关，与 D-1/D-3/D-4/D-7/D-11/D-15/D-16/
    D-19 同一口径）。
  · 说明书算例 2~6：数值字段逐位对拍（算例 2/3/4/5/6 全部命中到打印精度）。
"""
import math
import re

from ..core.intio import read_lines

PROGRAM_ID = "D-25"
TITLE = "底流式消能工水力计算程序"
AUTHOR = "卢礼标"

G = 9.8               # 重力加速度（由权威 OUT 的 Fr 反演唯一锁定，见 U6）
SIGMA = 1.05          # 稍淹没水跃安全系数 σ_j
LK_COEF = 0.75        # 消力池长度系数 Lk = 0.75·Lj（说明书：采用 0.75Lj）
C_WEIR = 4.43         # 消力坎堰流式常数 = √(2g)（原著取 4.43）

JX_NAMES = {1: "挖深式消力池", 2: "消力坎式消力池", 3: "综合式消力池"}
PK_NAMES = {1: "等宽消力池", 2: "扩散消力池"}

# ---------------------------------------------------------------
# 消力坎淹没系数 p(x) 等效标定表（U1：非原著原表；x = hs/H10 = (Ht−C)/H10）
# 锚点来源：说明书算例 3/4/5/6 反演（见 U1），锚点由本内核自身的 hc″ 中间量反解：
#   x=0.36026/0.37299 → p=1（算例6L1/算例5，自由溢流）
#   x=0.633944 → 0.967609（算例3）
#   x=0.835252 → 0.819335（算例6 第2级）
#   x=0.860358 → 0.781054（算例4）
# ---------------------------------------------------------------
P_X = [0.000000, 0.400000, 0.633944, 0.835252, 0.860358, 1.000000]
P_Y = [1.000000, 1.000000, 0.967609, 0.819335, 0.781054, 0.650000]
P_MODE = "auto"       # "auto" = 标定表；"free" = 恒取 1（自由溢流对照）

MAX_LEVELS = 6        # 多级消能上限（防病态输入死循环）


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


def p_submerged(x, mode=None):
    """消力坎淹没系数 p = f(hs/H10)。mode 缺省取模块常量 P_MODE。"""
    if (mode or P_MODE) == "free":
        return 1.0
    return interp(P_X, P_Y, x)


def p_handbook_weir(x):
    """备选：《水工设计手册》卷7 §5.3 宽顶堰淹没系数式（0.72<x<0.98，域外取 1）。"""
    if x <= 0.72:
        return 1.0
    x = min(x, 0.98)
    return 2.31 * x * (1.0 - x) ** 0.4


# ---------------------------------------------------------------
# 水力学基本解
# ---------------------------------------------------------------

def solve_hc(e01, q, phi, g=G):
    """
    解收缩断面水深 hc：hc + q²/(2g·φ²·hc²) = E01（吴持恭式 9.4/9.5）。
    方程在 (0, (2C)^(1/3)] 上单调下降（C=q²/(2gφ²)），取最小正根。
    """
    if phi <= 0:
        raise ValueError("泄水建筑物（或消力池出口）流速系数必须为正")
    c = q * q / (2.0 * g * phi * phi)
    lo, hi = 1e-12, (2.0 * c) ** (1.0 / 3.0)
    if lo + c / (lo * lo) - e01 < 0:
        raise ValueError("收缩水深能量方程无正根（E01 过小）")
    if hi + c / (hi * hi) - e01 > 0:
        raise ValueError("收缩水深能量方程在物理区间内无根（E01 过大）")
    for _ in range(300):
        m = 0.5 * (lo + hi)
        if m + c / (m * m) - e01 > 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def conjugate_equal_width(hc, q, g=G):
    """等宽消力池共轭水深 Hc″ = hc/2·(√(1+8Fr²)−1)，Fr=Vc/√(g·hc)。"""
    vc = q / hc
    fr = vc * vc / (g * hc)
    return hc / 2.0 * (math.sqrt(1.0 + 8.0 * fr) - 1.0)


def conjugate_expanding(hc, q, B, Bk, Q, g=G):
    """
    扩散消力池共轭水深（超越方程，单调二分）：
      (hc² − Hc″²)·(B+Bk) = 4Q²(1/(Bk·Hc″) − 1/(B·hc))/g
    """
    def f(x):
        return (hc * hc - x * x) * (B + Bk) - 4.0 * Q * Q * (1.0 / (Bk * x) - 1.0 / (B * hc)) / g
    # 物理根在 x > hc 一侧（x→0+ 时 f→−∞，x=hc 时 f≥0，x→+∞ 时 f→−∞）
    lo = hc * (1.0 + 1e-12)
    hi = hc + 1.0
    for _ in range(200):
        if f(hi) < 0:
            break
        hi *= 2.0
    else:
        raise ValueError("扩散消力池共轭水深方程无变号根")
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def jump_length_equal_width(hc, hc2):
    """等宽消力池水跃长度 Lj = 6.9(Hc″ − Hc)。"""
    return 6.9 * (hc2 - hc)


def jump_length_expanding(hc, hc2, fr, B, theta_rad):
    """扩散消力池水跃长度（分段式，说明书 3.B）。"""
    if fr < 6.0:
        return (1.0 + 0.6 * fr) * hc2
    if fr <= 17.0:
        return 4.6 * hc2
    L = 10.3 * hc * (math.sqrt(fr) - 1.0) ** 0.81
    return B * L / (B + 0.1 * L * math.tan(theta_rad))


PI_DEG = 3.14       # 弧度→度换算常数（原著取 π=3.14，由算例 2/5 双点反演锁定，见 U4）


def _theta_from_geometry(B, Bk, Lk):
    """
    扩散消力池侧墙扩散角 θ = arctg[(Bk−B)/(2·Lk)]（弧度→度按 π=3.14）。
    由算例 2（θ=8.3723°）与算例 5（θ=6.6772°）双点唯一锁定，见 U4。
    """
    if Lk <= 0:
        return 0.0
    return math.atan((Bk - B) / (2.0 * Lk)) * 180.0 / PI_DEG


def weir_height(hc2, q1, m, S, ht, p_mode=None, p_value=None, g=G):
    """
    消力坎高度 C（说明书式 6 / EXE BSTR 原文）：
      C = 1.05·Hc″ + q1²/(2g)/(1.05·Hc″)² − (q1/(4.43·p·m))^(2/3) − S
    其中 p = 消力坎淹没系数（p_mode="free" 时取 1；"auto" 时由标定表按
    x=(Ht−C)/H10 隐式确定），H10 = 1.05Hc″ + q1²/(2g)/(1.05Hc″)² − S − ... 见下。
    """
    hd = SIGMA * hc2
    k2 = q1 * q1 / (2.0 * g * hd * hd)
    if p_value is not None:
        p = float(p_value)
    elif (p_mode or P_MODE) == "free":
        p = 1.0
    else:
        # 隐式：p 依赖 C（经 x=(Ht−C)/H10），对 C 二分
        def resid(c):
            h10 = hd + k2 - c - S          # 坎上总水头
            if h10 <= 0:
                return -1e9
            x = (ht - c) / h10 if h10 > 0 else 1.0
            pp = p_submerged(x, p_mode)
            return (q1 / (C_WEIR * pp * m)) ** (2.0 / 3.0) - h10
        lo, hi = -abs(hd) - 10.0, hd + k2 - S - 1e-9
        flo, fhi = resid(lo), resid(hi)
        if flo * fhi > 0:
            # 退化：直接按 p=1 给出
            p = 1.0
            return hd + k2 - (q1 / (C_WEIR * p * m)) ** (2.0 / 3.0) - S, p
        # resid(c) 随 c 单调递增（h10 = hd+k2−c−S 随 c 递减）
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            if resid(mid) > 0:
                hi = mid
            else:
                lo = mid
        c = 0.5 * (lo + hi)
        h10 = hd + k2 - c - S
        p = p_submerged((ht - c) / h10, p_mode)
        return c, p
    c = hd + k2 - (q1 / (C_WEIR * p * m)) ** (2.0 / 3.0) - S
    return c, p


# ---------------------------------------------------------------
# 单级计算
# ---------------------------------------------------------------

S_TOL = 1e-3   # 池深试算收敛容差（由权威 D-25-1.OUT 七个显示值唯一锁定，见 U7）


def level_calc(jx, pk, E0, Q, Ht, B, Bk, f1, f2, f3, m, S=0.0,
               p_mode=None, p_value=None, iterate_S=True, g=G, tol=S_TOL,
               max_iter=500):
    """
    计算一级消能。返回 dict：
      q (收缩断面单宽流量)、hc、Fr、hc2、Lj、Lk、S、C、p、theta
    E0 为以下游河床为基准的上游总能头；S 为池深（JX=1 时试算）。
    """
    q = Q / B
    S = float(S)
    if jx == 3:
        iterate_S = False          # 综合式池深由输入给定

    def one(Sv):
        e01 = E0 + Sv
        hc = solve_hc(e01, q, f1, g)
        if pk == 1:
            hc2 = conjugate_equal_width(hc, q, g)
            vc = q / hc
            fr = vc / math.sqrt(g * hc)
            hc2_tmp = hc2
            theta = 0.0
            Lj = jump_length_equal_width(hc, hc2)
        else:
            # 扩散：共轭水深与 Lj 都依赖 θ, 而 θ 依赖 Lk=0.75Lj → 迭代
            theta = 0.0
            for _ in range(200):
                hc2 = conjugate_expanding(hc, q, B, Bk, Q, g)
                vc = q / hc
                fr = vc / math.sqrt(g * hc)
                Lj = jump_length_expanding(hc, hc2, fr, B, theta)
                Lk = LK_COEF * Lj
                th_new = _theta_from_geometry(B, Bk, Lk)
                if abs(th_new - theta) < 1e-14:
                    theta = th_new
                    break
                theta = th_new
            return hc, hc2, fr, Lj, theta
        return hc, hc2, fr, Lj, theta

    if jx == 1 and iterate_S:
        Sv = 0.0
        hc, hc2, fr, Lj, theta = one(0.0)
        if hc2 <= Ht:
            # 不需建消力池
            return {"q": q, "hc": hc, "hc2": hc2, "Fr": fr, "Lj": Lj,
                    "Lk": LK_COEF * Lj, "S": 0.0, "C": None, "p": None,
                    "theta": theta, "need": False}
        q1 = q if pk == 1 else Q / Bk
        for _ in range(max_iter):
            hc, hc2, fr, Lj, theta = one(Sv)
            h2 = SIGMA * hc2
            dz = q1 * q1 / (2.0 * g) * (1.0 / (f2 * Ht) ** 2 - 1.0 / h2 ** 2)
            Sv_new = h2 - Ht - dz
            if abs(Sv_new - Sv) < tol:
                Sv = Sv_new
                break
            Sv = Sv_new
        hc, hc2, fr, Lj, theta = one(Sv)
        return {"q": q, "hc": hc, "hc2": hc2, "Fr": fr, "Lj": Lj,
                "Lk": LK_COEF * Lj, "S": Sv, "C": None, "p": None,
                "theta": theta, "need": True}

    # JX=2 / JX=3：池深确定，求坎高 C
    hc, hc2, fr, Lj, theta = one(S)
    q1 = q if pk == 1 else Q / Bk
    C, p = weir_height(hc2, q1, m, S, Ht, p_mode=p_mode, p_value=p_value, g=g)
    return {"q": q, "hc": hc, "hc2": hc2, "Fr": fr, "Lj": Lj,
            "Lk": LK_COEF * Lj, "S": S, "C": C, "p": p,
            "theta": theta, "need": True}


def compute_levels(jx, pk, E0, Q, Ht, B, Bk, f1, f2, f3, m, S=0.0,
                   p_mode=None, p_value=None, multi_level=True, g=G):
    """
    完整计算（含多级消能）。返回 (levels, 说明信息)。
    多级判据（U3）：以 hc″(i) 作为第 i+1 级的上游总能头 E0(i+1)，
      q(i+1)=Q/B(i+1)（扩散取 Bk），若 hc″(i+1) > Ht 则仍需第 i+1 级。
    """
    levels = []
    lv = level_calc(jx, pk, E0, Q, Ht, B, Bk, f1, f2, f3, m, S=S,
                    p_mode=p_mode, p_value=p_value)
    lv["E0"] = E0
    lv["B"] = B
    lv["Bk"] = Bk if pk == 2 else B
    levels.append(lv)
    if jx == 1 or not multi_level:
        return levels

    Bi = Bk if pk == 2 else B
    e0i = lv["hc2"]
    for i in range(2, MAX_LEVELS + 1):
        # 第 2 级起：下游为扩散段末端宽度 Bk 的等宽段（U3），故按等宽口径
        hc_i = solve_hc(e0i, Q / Bi, f1, g)
        hc2_i = conjugate_equal_width(hc_i, Q / Bi, g)
        if hc2_i <= Ht:
            break
        nxt = level_calc(jx, 1, e0i, Q, Ht, Bi, Bi, f1, f2, f3, m, S=0.0,
                         p_mode=p_mode, p_value=p_value)
        nxt["E0"] = e0i
        nxt["B"] = Bi
        nxt["Bk"] = Bi
        levels.append(nxt)
        e0i = nxt["hc2"]
    return levels


# ---------------------------------------------------------------
# 解析
# ---------------------------------------------------------------

_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")


def parse(data):
    """
    解析输入。data: dict 或 INT 文件路径。
    INT 组织（说明书三、输入数据说明）：
      首字段：工程名（可选，非数值）
      数值流顺序：JX, PK, f1, f2, f3, m, E0, Q, Ht, B[, Bk][, S]
        · JX 消能结构类型 1/2/3；PK 平面型式 1 等宽 / 2 扩散
        · 扩散消力池给 Bk（等宽不给）
        · 综合式（JX=3）另给池深 S
      兼容两种书写：逗号分隔（D-25-1.INT…）与"每行一个数"（D-25-5A.INT）。
    """
    if isinstance(data, dict):
        return data
    lines = read_lines(data)
    name = None
    toks = []
    for i, ln in enumerate(lines):
        ln = ln.strip()
        if not ln:
            continue
        parts = [p.strip() for p in ln.split(",")]
        for j, p in enumerate(parts):
            if not p:
                continue
            if i == 0 and j == 0 and not _NUM.fullmatch(p):
                name = p
                continue
            if _NUM.fullmatch(p):
                toks.append(float(p))
    if len(toks) < 10:
        raise ValueError("D-25 数据不足：至少需要 JX,PK,f1,f2,f3,m,E0,Q,Ht,B")
    jx = int(round(toks[0]))
    pk = int(round(toks[1]))
    if jx not in JX_NAMES:
        raise ValueError(f"消能结构类型 JX={jx} 不在 1/2/3 之内")
    if pk not in PK_NAMES:
        raise ValueError(f"消力池平面型式 PK={pk} 不在 1/2 之内")
    f1, f2, f3, m = toks[2], toks[3], toks[4], toks[5]
    E0, Q, Ht, B = toks[6], toks[7], toks[8], toks[9]
    rest = toks[10:]
    Bk = B
    S = 0.0
    if pk == 2:
        if not rest:
            raise ValueError("扩散消力池（PK=2）须给末端宽度 Bk")
        Bk = rest.pop(0)
    if jx == 3:
        if not rest:
            raise ValueError("综合式消力池（JX=3）须给池深 S")
        S = rest.pop(0)
    if m <= 0:
        m = 0.42          # 缺省（说明书算例用 0.42）
    return {"程序": PROGRAM_ID, "工程名": name, "JX": jx, "PK": pk,
            "f1": f1, "f2": f2, "f3": f3, "m": m,
            "E0": E0, "Q": Q, "Ht": Ht, "B": B, "Bk": Bk, "S": S}


# ---------------------------------------------------------------
# 计算
# ---------------------------------------------------------------

def compute(params, p_mode=None, p_value=None):
    """执行计算，返回结构化结果。"""
    jx = params["JX"]
    pk = params["PK"]
    levels = compute_levels(
        jx, pk, params["E0"], params["Q"], params["Ht"], params["B"],
        params["Bk"], params["f1"], params["f2"], params["f3"], params["m"],
        S=params.get("S", 0.0), p_mode=p_mode, p_value=p_value)
    res = {"程序": PROGRAM_ID, "工程名": params.get("工程名"),
           "结构类型": JX_NAMES[jx], "平面型式": PK_NAMES[pk],
           "基本资料": {k: v for k, v in params.items()
                        if k not in ("程序",)},
           "级数": len(levels), "各级": []}
    for i, lv in enumerate(levels, 1):
        res["各级"].append({
            "级": i, "E0": lv["E0"], "q": lv["q"], "B": lv["B"],
            "Bk": lv["Bk"], "hc": lv["hc"], "Fr": lv["Fr"], "hc2": lv["hc2"],
            "Lj": lv["Lj"], "Lk": lv["Lk"], "S": lv["S"], "C": lv["C"],
            "p": lv["p"], "theta": lv["theta"], "need": lv["need"],
        })
    return res


# ---------------------------------------------------------------
# 输出（复刻原著 D-25.OUT 版式）
# ---------------------------------------------------------------

def _fmt4(x):
    return "%.4f" % x


def _vbstr(x):
    """VB6 Str() 口径的正数前导空格 + 去尾零 + 小于 1 去前导 0。"""
    if x == 0:
        return " 0"
    s = ("%.6f" % x).rstrip("0").rstrip(".")
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    if x >= 0:
        return " " + s
    return s


def _vb0(x):
    """保留前导零、去尾零、无前导空格（原著 f1/f2 行的打印口径）。"""
    s = ("%.6f" % x).rstrip("0").rstrip(".")
    return s if s else "0"


def _echo(prefix, gap, text, suffix=""):
    return prefix + " " * gap + text + suffix


def _result_line(prefix, gap, value, suffix=""):
    return prefix + " " * gap + value + suffix


_LINE = " " + "*" * 71
_TITLELIN = " *****               底流式消能工水力计算程序 D-25                 *****"


def _input_lines(p):
    jx, pk = p["JX"], p["PK"]
    L = []
    L.append(_echo("    1.消能结构类型:", 24, JX_NAMES[jx]))
    L.append(_echo("    2.消力池平面型式:", 22, PK_NAMES[pk]))
    L.append(_echo("    3.以下游河床为基准面的上游总能头 E0:", 2,
                   _vbstr(p["E0"]), " (米)"))
    L.append(_echo("    4.下泄流量 Q:", 25, _vbstr(p["Q"]), " (立方米/秒)"))
    L.append(_echo("    5.下游水深 Ht:", 24, _vbstr(p["Ht"]), " (米)"))
    L.append(_echo("    6.泄水建筑物的流速系数 f1:", 12, _vb0(p["f1"])))
    L.append(_echo("    7.消力池出口流速系数 f2:", 14, _vb0(p["f2"])))
    L.append(_echo("    8.消力坎流速系数 f3:", 18, _vbstr(p["f3"]), " "))
    L.append(_echo("    9.消力坎流量系数 m:", 19, _vbstr(p["m"]), " "))
    if pk == 2:
        L.append(_echo("    10.消力池首端宽度 B:", 22, _vbstr(p["B"]), " (米)"))
        L.append(_echo("    11.消力池末端宽度 Bk:", 20, _vbstr(p["Bk"]), " (米)"))
    else:
        L.append(_echo("    10.消力池宽度 B:", 22, _vbstr(p["B"]), " (米)"))
    return L


_FORMULA_COMMON = [
    "    1.收缩水深hc的基本算式",
    "        E01=Hc+q*q/(2*g*f1^2*Hc^2)",
    "     式中 q－收缩断面的单宽流量;",
    "         E01-以消力池底面为基准面的上游总能头.",
    "    2.共轭水深Hc″ ",
    "     A.等宽消力池",
    "        Hc″=Hc*(SQR(1+8*Fr^2)-1)/2",
    "      式中 Fr－收缩断面的弗汝德数 Frc=Vc/SQR(g*Hc),Vc为收缩断面流速.",
    "     B.扩散消力池",
    "       按下列关系式计算",
    "        (Hc^2-Hc″^2)*(B+Bk)=4*Q*Q*(1/Bk/Hc″-1/B/Hc)/g",
    "    3.水跃长度Lj",
    "     A.等宽消力池 ",
    "         Lj=6.9*(Hc″-Hc)",
    "     B.扩散消力池",
    "        Fr<6时,  Lj=(1+0.6*Fr)*Hc″",
    "        Fr>=6时, Lj=4.6*Hc″",
    "        Fr>17时, Lj=B*L/(B+.1*L*tgθ)",
    "       式中L－等宽矩形明渠中水跃长度,L=10.3*Hc(SQR(Fr)-1)^.81.",
    "    4.消力池长度Lk ",
    "        Lk=(0.7～0.8)*Lj           (采用0.75Lj)",
]
_FORMULA_S = [
    "    5.消力池深度S",
    "        S=1.05*Hc″－Ht－q1*q1*(1/(f2*Ht)^2-1/(1.05*Hc″)^2)/(2*g)",
    "       式中q1－消力池出口断面的单宽流量.",
]
_FORMULA_C = [
    "    5.消力坎高度C",
    "        C=1.05*Hc″＋q1*q1/(2*g)/(1.05*Hc″)^2－(q1/(4.43*p*m))^(2/3)－S",
    "       式中p－消力坎淹没系数,自由溢流时P=1;",
    "           m－消力坎流量系数;",
    "           S－综合式消力池深度,对于消力坎式S=0.",
]


def _result_lines(p, lv):
    jx, pk = p["JX"], p["PK"]
    Bk, theta = lv["Bk"], lv["theta"]
    L = []
    L.append(_result_line("    1.收缩水深Hc:", 16, _fmt4(lv["hc"]), " (米)"))
    L.append(_result_line("    2.收缩断面的弗汝德数Fr:", 6, _fmt4(lv["Fr"])))
    L.append(_result_line("    3.共轭水深Hc″:", 14, _fmt4(lv["hc2"]), " (米)"))
    L.append(_result_line("    4.水跃长度Lj:", 15, _fmt4(lv["Lj"]), " (米)"))
    if jx == 1:
        L.append(_result_line("    5.消力池深度S:", 15, _fmt4(lv["S"]), " (米)"))
    elif jx == 2:
        L.append(_result_line("    5.消力坎高度C:", 15, _fmt4(lv["C"]), " (米)"))
    else:
        L.append(_result_line("    5.消力池深度S:", 15, _fmt4(lv["S"]), " (米)"))
        L.append(_result_line("    6.消力坎高度C:", 15, _fmt4(lv["C"]), " (米)"))
    L.append(_result_line("    %d.消力池长度Lk:" % (6 if jx != 3 else 7), 13,
                          _fmt4(lv["Lk"]), " (米)"))
    nb = 7 if jx != 3 else 8
    if pk == 2:
        L.append(_result_line("    %d.消力池首端宽度B:" % nb, 15,
                              _fmt4(lv["B"]), " (米)"))
        L.append(_result_line("    %d.消力池末端宽度Bk:" % (nb + 1), 15,
                              _fmt4(Bk), " (米)"))
        L.append(_result_line("    %d.消力池侧墙的扩散角θ: 　" % (nb + 2), 9,
                              _fmt4(theta), " (度)"))
    else:
        L.append(_result_line("    %d.消力池宽度B:" % nb, 15,
                              _fmt4(lv["B"]), " (米)"))
    return L


def render(params, result, levels=None):
    """生成汉字计算书（复刻原著 D-25.OUT 版式）。"""
    if levels is None:
        levels = result if isinstance(result, list) else result["各级"]
    p = params
    jx = p["JX"]
    L = ["", _LINE, _TITLELIN, _LINE, "",
         "            工程名:%s" % (p.get("工程名") or "未输入文件名"),
         "   一.基本资料"]
    L.extend(_input_lines(p))
    L.append("")
    L.append("   二.计算中采用的基本公式")
    L.extend(_FORMULA_COMMON)
    if jx == 1:
        L.extend(_FORMULA_S)
    else:
        L.extend(_FORMULA_C)
    L.append("")
    L.append("   三.计算结果 ")
    for i, lv in enumerate(levels, 1):
        if len(levels) > 1:
            L.append("           第 %d 级消能" % i)
        L.extend(_result_lines(p, lv))
    L.append("")
    return "\n".join(L)


def run(data, out_txt=None, out_json=None, fmt="text", p_mode=None,
        p_value=None):
    """统一入口。data: INT 文件路径 | dict"""
    params = parse(data)
    result = compute(params, p_mode=p_mode, p_value=p_value)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [], result)
    else:
        text = render(params, result)
    if out_txt:
        from ..core.outgen import write_out
        write_out(out_txt, text)
    if out_json:
        from ..core.outgen import write_json
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(txt)
