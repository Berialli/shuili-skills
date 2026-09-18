# -*- coding: utf-8 -*-
"""
D-19 引水系统水击压力计算程序 —— 内核
======================================
复刻《水利水电工程设计计算程序集》D-19 程序（作者：姚廉华）。

功能
----
计算常规电站与抽水蓄能电站**引水系统**的调节保证、水击压力及水力过渡过程
近似解。程序按"管道特征系数法"（Allievi / 连锁方程简化式）求解，先建立
引水管道与尾水管道沿程的 ∑LV 与加权平均波速 C，再据第一相 / 第末相公式
算出水击压力上升、下降相对值，最后按各点累积 ∑LV 把总水击沿管线分配。

要复刻的量（原著打印口径）
--------------------------
  I/D/V/LV/∑LV/C  各分段引水管道 ∑LV 及波速
  I1/D1/V1/LV1/∑LV1/C1  尾水管道同上
  X/Y/Z/VCP/CP, X1/Y1/Z1/VCP1/CP1, XX/YY/ZZ/VA/CA
  P（管道特征系数）、ZM（水击上升相对值）、ZM1（水击下降相对值）
  DH/HMAX/DH1/HMIN 及逐点 DP/DP1/HP/HP1

核心公式（逐字还原自 D-19Intro.rtf 的 30 个 OLE 公式对象）
----------------------------------------------------------
公式对象解包路径：\\objdata → hex → OLE2 → "Equation Native" → 跳 28 B →
MTEF v3 记录流；并用内嵌 WMF 预览（META_EXTTEXTOUT 的字符+坐标）逐字核对。

  常规电站（对象 00~05）：
    (1) 波速      C  = C0 / √( 1 + (EC/E)·(D/H) )      —— 由输出表唯一反演
    (2) 第一相上升 Z1 = 2σ / (1 + ρ·τ − σ)              [对象 00]
    (3) 第一相下降 Y1 = 2σ / (1 + ρ·τ + σ)              [对象 01]
    (4) 第末相上升 Zm = 2σ / (2 − σ)                    [对象 02]
    (5) 第末相下降 Ym = 2σ / (2 + σ)                    [对象 03]
    (6) 其中      ρ  = C·Vcp / (2·g·H0)                 [对象 04]
    (7)           σ  = ∑(L·V) / (g·H0·Ts)               [对象 05]
       τ = 起始开度；阀门由全开到全关 τ = 1.0（本例取 1）

  抽水蓄能（对象 06~29，本内核按反演式实现；原著无算例，未逐位验证）：
    (8)  P  = C·V0/(2gH0), P1 = C·U0/(2gH0),
         P2 = C·V0'/(2gH0), P3 = C·U0'/(2gH0)           [对象 14~17]
    (9)  F  = F(t)/F0,  F(t) = (T0−t1)/T0 或 (Ts−t1)/Ts  [对象 18~20]
    (10) T0 = 1.39/4.5·K1                               [对象 21]
    (11) K1 = 450·g·HR·Q0/(π²·WR²·η0²·ηP)               [对象 22]
    (12) WR² = 0.25·GD²                                 [对象 23]
    (13) t1 = 2L/C（引水）、2L1/C1（尾水）              [对象 24~25]
    (14) n  = L·V0/(g·H·Ts), n1 = L1·U0/(g·H·Ts),
         n2 = L·V0/(g·H·T0),  n3 = L·U0/(g·H·T0)        [对象 26~29]
    (15) 第一相：Z1 =  2P·[ (P·F²+1) − F·(P²F²+2P+1)^(1/2) ]
               Z2 = −2P1·[ (P1·F²+1) − F·(P1²F²+2P1+1)^(1/2) ]
               Z3 = −2P2·[ (P2·F²+1) − F·(P2²F²+2P2+1)^(1/2) ]
               Z4 = −2P3·[ (P3·F²+1) − F·(P3²F²+2P3+1)^(1/2) ]
    (16) 第末相：Zm1 = 0.5·[ n² + (n⁴+4n²)^(1/2) ]  (n→n1/n2/n3)
      即 Zm = (n/2)(n+√(n²+4))，与常规第末相公式同形。

关键裁决（由权威 D-19-1.OUT 唯一反演）
-------------------------------------
  A. 圆周率取 **3.1416**（VB6 常量口径）：以 π=3.14159265 计算时，
     20 行 LV 仅 18 行命中、20 行 ∑LV 仅 6 行命中（如 seg2 = 313.93 🆚
     权威 313.92）；改取 π = 3.1416 后 LV 表 20/20、∑LV 表 20/20 全部命中。
  B. 波速 C = C0/√(1+(EC/E)·(D/H))，其中 D 以 m、H（壁厚）以 cm **直接**
     代入（不做单位换算）。核验：seg1 C=1185.28、seg20 C=1197.37、
     尾水 C1=965.28 与权威逐位相同。
  C. 电站设计水头 **H0 = HU − HD1**（= 345 − 47 = 298.000）。反演：DH/ZM =
     129.201/0.433561 = 298.0002；DH1/ZM1 = 90.126/0.302436 = 298.0002。
  D. 管道特征系数 **P = CA·VA/(2·g·H0) = YY/(ZZ·2·g·H0)**（用**引水+尾水**
     全系统的加权平均波速/流速，而非仅引水管道）。反演：P = 0.63865063，
     权威打印 0.638651；若改用引水 CP·VCP/(2gH0) 则为 1.0164（差 59%）。
  E. **σ = YY/(g·H0·Ts)**（总 ∑LV ÷(g·H0·Ts)，仍是全系统口径）。反演：
     σ = 0.29194056，此时 Z1 = 0.4335609、Y1 = 0.3024359，与权威 0.433561 /
     0.302436 逐位相同；若用仅引水 Y/(gH0Ts) = 0.29193312，Z1 = 0.4335540
     （第 5 位小数不符）。
  F. 选值口径：ZM = max(Z1, Zm)、ZM1 取同一相的值。本例 Z1 = 0.433561 >
     Zm = 0.341833，故取第一相。原著对抽水蓄能明示"ζ1 > ζm 用第一相 /
     ζm > ζ1 用第末相"，常规电站按其类比裁决（本例唯一可判据）。
  G. 沿管水击分配：**DP(i) = ZM·H0·∑LV(i)/Y**（Y 为**引水**管道总 ∑LV，
     不含尾水；否则末端 DP(20) = 129.198 ≠ 129.201）。
     HP(i) = (HU − HT(i)) + DP(i)，HP1(i) = (HU − HT(i)) − DP1(i)。
     HMAX = (HU − HR) + DH，HMIN = (HU − HR) − DH1。

输入数据（.INT，逗号分隔；原著算例 D-19-1.INT）
---------------------------------------------
  第 1 行（12 或 13 个数；说明书算例行在第 5 位多一个 K=0，权威 .INT 无）：
     N, M, M1, L, [K,] N0, GD2, C0, EC, HU, HD, HU1, HD1
     L=1 常规电站；L=2 抽水蓄能电站
  第 2 行（6 个数）：
     QT, QP, TS, HR, HM, Z
  引水管道数组（各 M 个，共 6M 个）：
     D(1:M) 内径(m) / Q(1:M) 过流量(m³/s) / S(1:M) 管长(m) /
     H(1:M) 壁厚(cm) / E(1:M) 弹模(kg/cm²) / HT(1:M) 中心高程(m)
  尾水管道数组（各 M1 个，共 6·M1 个）：
     D1, Q1, S1, H1, E1, HN（含义同上，HN 为中心高程）

输出（逐字复刻原著 .OUT 版式，GBK）
-----------------------------------
  原始数据 → 引水/尾水各分段数据表 → RESULT OF LV & C（∑LV 及波速成果）
  → 其它有关成果（17 个常数）→ RESULT OF WATER HUMMER（水击成果）
  → 常规电站 DH/HMAX/DH1/HMIN + 各特征点 DP/DP1/HP/HP1 表

未闭合点（如实标注）
------------------
  1. **抽水蓄能（L=2）分支未经逐位验证**：原著算例目录仅存 D-19-1.INT（L=1），
     无 L=2 算例，无法对拍。本内核按对象 06~29 的反演式实现，下列符号绑定
     属**推断**：V0'/U0'（水泵工况初始流速）取 V0·(QP/QT)、U0·(QP/QT)；
     K1 中 η0 取机组转速 N0、ηP 取水轮机效率 Z；T0 = 1.39·K1/4.5（原文
     "1.39/4.5K1" 的斜杠归属不可判）。该分支输出附"未验证"标注。
  2. 说明书算例数据（D-19.IN）第 1 行含 K=0，而权威 .INT 不含；本内核两种
     写法都能解析（按 6(M+M1)+18 与 6(M+M1)+19 判别）。
  3. 权威 .OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-19-1.out」为原著
     运行期路径回显（机器相关），本内核不生成（与 D-1/D-3/D-4/D-7/D-11/
     D-16 同一口径），版式比对自第 2 行起算。
  4. **逐点表 5 个单元（行 6 HP1、行 16 HP、行 17 HP、行 18 HP1、行 19 HP1）
     落在 ±0.0005 取整边界上，无法逐位复现**（详见 d19_verify.py 的区间反证）：
     内核值分别比边界低 7.4e-7 / 4.5e-6 / 4.4e-5 / 3.9e-5 / 2.3e-5，
     取整后比权威小 0.001。反证：对行 17 由权威 DP=126.706 得 share<0.9806889，
     由权威 HP=422.207（hp=291.217+DH·share）得 share≥0.9806889 —— 权威的
     DP 列与 HP 列在该行**互不相容**（差值恰为打印量子 0.001），说明原著 HP
     列并非由同一 share 的 base+DP 直接打印，而是独立浮点路径的结果；
     对 DH、DH1 做 ±4×10⁻⁶ 相对量扫描（共 1600 个取值）均无任何一个取值能
     同时命中 DP/DP1/HP/HP1 全部 100 个单元（最优仍余 2~3 个单元）。
     故本内核保留 double 口径，把这 5 个单元如实记为 DECL（非 FAIL）。

知识库对照结果（教学母本 → 程序实现）
-----------------------------------
  KB 路径：D:\\WorkBuddy知识库\\水利知识库\\

  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 波速 C = C0/√(1+(EC/E)(D/H)) | `02_水利教材精读/水力学与工程水文学/水力学核心公式速查.md` §5.6（a=√(K/ρ)/√(1+K·d/(E·δ))，钢管 E=2.06e11 Pa，**明确标注 D-19**）；`01_水工设计手册精读/卷8_分章/03_第3章_调压设施.md` §2.1（a=√[Ew/γ]/√[1+(Ew·D)/(e·E·C)]）；`水电站_第4版` ch9 式9-7 | 形式**一致** | 采用；单位口径分歧见下 |
  | 直接水击 ΔH = a·v₀/g | 速查 §5.5；手册 §2.3（ΔH=−a(v−v₀)/g）；水电站 式9-13 | **一致** | 采用（τ₁→0 时 2σH₀ = a·v₀/g） |
  | 水击特性系数 ρ = a·v₀/(2gH₀) | 手册 §2.4/§2.6（ρ=a·v₀/(2gH₀)、ρeq=aeq·Ueq/(2gH₀)）；水电站 9.3（β=c·Vmax/(2gH₀)） | 形式**一致** | 采用；作用范围分歧见下 |
  | 管道特性系数 σ = LV₀/(gH₀Ts) | 水电站 ch9 9.3（σ=LV0/(gH0Ts)） | **一致** | 采用（程序取全系统 ΣLV） |
  | 判别：ρ<1 第一相最不利 / ρ>1 末相 | 水电站 ch9 9.3；手册 §2.4 | **一致** | 采用；且 ζ₁>ζm ⟺ ρ<1（代数恒等，故与"取大者"等价） |
  | 等价管三原则 L=ΣLi、tr=Σtri、Σ(Li·ui)=L·Ueq | 手册 §2.6 | **一致** | 采用 → VCP=Y/X、CP=X/Z |
  | 沿管分布"按动能权重分配到各段" | 手册 §2.6 | **一致** | 采用 → DP(i)=DH·ΣLV(i)/ΣLV |
  | 第一相式 ζ₁ = 2σ/(1+ρτ−σ) | 手册 §2.4 ζ₁=2ρ(τ₀−τ₁)/(1+ρτ₀) | 分子**一致**（同为 2ρ(τ₀−τ₁)=2σ）；分母取 τ（程序）或 τ₀（手册）有别 | **以权威 OUT 定案**用程序式；见分歧 A |
  | 末相式 ζm = 2σ/(2−σ) | 手册 §2.4 ζm=[√(1+4ρ²τ₀²)−1]/(2ρτ₀²) | **不一致** | **以权威 OUT 定案**用程序式；见分歧 B |
  | 一相水击沿管二次曲线 ζ(x)=(1−x/L)(1−(x/L)(1−τ₁/τ₀))ζmax | 手册 §2.5 | **不一致** | 程序用 §2.6 动能权重；见分歧 C |

  分歧 A（第一相式分母）：程序式 ζ₁=2σ/(1+ρτ−σ)（τ=1）与手册式 ζ₁=2ρ(τ₀−τ₁)/(1+ρτ₀)
  分子相同、分母不同。注：因 τ₁ = 1 − t_r/Ts = 1 − σ/ρ，有 1+ρτ₁ ≡ 1+ρ−σ，
  故程序式 ≡ 2σ/(1+ρτ₁)（用一相末开度开方；手册用初始开度）。
  反证：取手册式 ζ₁ = 2σ/(1+ρ) = 0.356318 ≠ 权威 ZM = 0.433561（差 0.077）。
  定案：权威 OUT。
  分歧 B（末相式）：程序 ζm = 2σ/(2−σ) = 0.341833；手册式（τ₀=1）
  ζm = [√(1+4ρ²)−1]/(2ρ) = 0.487113，差 0.145。本例 ρ = 0.638651 < 1 判第一相，
  末相值不参与输出，故对权威 OUT 无影响。定案：权威 OUT（保留程序式）。
  分歧 C（沿管分布）：手册 §2.5 的一相二次曲线给 DP(1)/DH = 0.07030，
  权威为 10.658/129.201 = 0.082494；程序用 §2.6"按动能权重分配到各段"
  （DP(i) ∝ ΣLV(i)）逐位命中。定案：权威 OUT。
  分歧 D（等价管作用范围）：手册 §2.6 的等价管按**引水管道**算 ρeq = CP·VCP/(2gH₀)
  = 1.016405；程序用 CA·VA/(2gH₀)（引水 + 尾水全系统）= 0.638651。
  反证：若取手册口径 ρ = 1.016405，则 ZM = 2σ/(1+ρ−σ) = 0.3374、
  DH = 100.55 m ≠ 权威 129.201 m；且 ρ>1 会误判为末相。定案：权威 OUT。
  分歧 E（波速单位口径）：教材/手册为同一单位制下 d/δ；本程序输入 D 以 m、
  壁厚 H 以 cm，而计算式中**直接**用 D(m)/H(cm) 的数值比（不再换算），
  相当于把管壁柔度项缩小 100 倍。反证：按严格单位（D=300 cm、δ=1.2 cm）
  C₁ = 1200/√(1+0.01×250) = 641.43 m/s ≠ 权威 1185.28 m/s，20 行 LV/∑LV
  及"其它有关成果"全部不符。定案：权威 OUT。
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-19"
TITLE = "引水管道水击压力计算程序"
AUTHOR = "姚廉华"
HEAD_NAME = "D-19"

G = 9.81
# 圆周率：原著（VB6）常量口径 3.1416 —— 由权威 OUT 的 20 行 LV/∑LV 表唯一反演
PI = 3.1416


# ============================================================
# 打印取整（VB6 Format 口径：十进制四舍五入，半值进位）
# ============================================================

def q(x, d):
    """按打印量子 d 位小数做"十进制半进位"取整，返回 float。"""
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    """原著打印口径：十进制半进位取整后按固定列宽输出。"""
    return "%*.*f" % (w, d, q(x, d))


def FI(x, w):
    """整数列（%.0f 口径）。"""
    return "%*d" % (w, int(q(x, 0)))


# ============================================================
# 基本水力量
# ============================================================

def wave_speed(C0, EC, E, D, H):
    """
    分段水击波传播速度  C = C0 / √( 1 + (EC/E)·(D/H) )
    （D 以 m、H 以 cm 直接代入，与原著一致；见 module docstring 裁决 B）。
    """
    ratio = (EC / E) * (D / H)
    return C0 / math.sqrt(1.0 + ratio)


def velocity(Q, D):
    """断面平均流速 V = Q / (π·D²/4)，π 取 3.1416。"""
    return Q / (PI * D * D / 4.0)


def pipe_profile(D, Q, S, H, E, C0, EC):
    """
    单条管道（引水或尾水）的逐段成果。
    返回 dict：C[] 波速、V[] 流速、LV[]、CUM[] 累积 ∑LV、X、Y、Z、VCP、CP。
    """
    n = len(D)
    C = [wave_speed(C0, EC, E[i], D[i], H[i]) for i in range(n)]
    V = [velocity(Q[i], D[i]) for i in range(n)]
    LV = [S[i] * V[i] for i in range(n)]
    CUM = []
    acc = 0.0
    for v in LV:
        acc += v
        CUM.append(acc)
    X = sum(S)
    Y = CUM[-1] if CUM else 0.0
    Z = sum(S[i] / C[i] for i in range(n)) if n else 0.0
    VCP = Y / X if X else 0.0
    CP = X / Z if Z else 0.0
    return {"C": C, "V": V, "LV": LV, "CUM": CUM, "X": X, "Y": Y, "Z": Z,
            "VCP": VCP, "CP": CP}


# ============================================================
# 解析
# ============================================================

def parse(data):
    """
    解析输入。data：dict（直接返回）| .INT 文件路径。

    返回 dict，键见 module docstring「输入数据」。
    兼容首部 12 个数（权威 .INT，无 K）与 13 个数（说明书算例，含 K=0）。
    """
    if isinstance(data, dict):
        p = dict(data)
        for k in ("N", "M", "M1", "L"):
            if k in p:
                p[k] = int(round(p[k]))
        p.setdefault("K", 0)
        return p

    nums = read_numbers(data)
    total = len(nums)
    if total < 20:
        raise ValueError("D-19 输入数据过少（%d 个数）" % total)

    M = int(round(nums[1]))
    M1 = int(round(nums[2]))
    L = int(round(nums[3]))
    if M <= 0 or M1 <= 0:
        raise ValueError("引水/尾水管道分段数必须为正（M=%d, M1=%d）" % (M, M1))

    need18 = 6 * (M + M1) + 18
    need19 = 6 * (M + M1) + 19
    if total == need18:
        head = 18
        has_K = False
    elif total == need19:
        head = 19
        has_K = True
    else:
        raise ValueError(
            "D-19 输入数据个数 %d 与 M=%d、M1=%d 不符（应为 %d 或 %d）"
            % (total, M, M1, need18, need19))

    if has_K:
        keys = ["N", "M", "M1", "L", "K", "N0", "GD2", "C0", "EC",
                "HU", "HD", "HU1", "HD1", "QT", "QP", "TS", "HR", "HM", "Z"]
        scal = dict(zip(keys, nums[:head]))
    else:
        keys = ["N", "M", "M1", "L", "N0", "GD2", "C0", "EC",
                "HU", "HD", "HU1", "HD1", "QT", "QP", "TS", "HR", "HM", "Z"]
        scal = dict(zip(keys, nums[:head]))
        scal["K"] = 0

    p = dict(scal)
    p["N"] = int(round(p["N"]))
    p["M"] = M
    p["M1"] = M1
    p["L"] = L
    if L not in (1, 2):
        raise ValueError("电站类型信息 L 只能取 1（常规）或 2（抽水蓄能），实得 %s" % L)

    off = head
    arrs = {}
    for name in ("D", "Q", "S", "H", "E", "HT"):
        arrs[name] = list(nums[off:off + M])
        off += M
    for name in ("D1", "Q1", "S1", "H1", "E1", "HN"):
        arrs[name] = list(nums[off:off + M1])
        off += M1
    p.update(arrs)
    return p


# ============================================================
# 计算
# ============================================================

def _conventional(prof_in, prof_tl, P, sigma):
    """
    常规电站：第一相 / 第末相四式。

    选值口径（与知识库一致）：
      《水电站》第4版 ch9 9.3 —— 管道特性系数 ρ = c·Vmax/(2gH₀) 判水锤类型，
      ρ<1 时首相（第一相）水锤最不利，ρ>1 时末相（极限）水锤最不利。
      代数恒等：ζ₁ > ζm  ⟺  2σ/(1+P−σ) > 2σ/(2−σ)  ⟺  P < 1，
      故"取大者"与教材"P<1 取第一相"完全等价。本例 P = 0.638651 < 1 → 第一相。
    """
    z1 = 2.0 * sigma / (1.0 + P * 1.0 - sigma)       # τ = 1：全开→全关
    y1 = 2.0 * sigma / (1.0 + P * 1.0 + sigma)
    zm = 2.0 * sigma / (2.0 - sigma)
    ym = 2.0 * sigma / (2.0 + sigma)
    if P < 1.0:                    # 等价于 ζ₁ > ζm
        return z1, y1, ("第一相", z1, y1, zm, ym)
    return zm, ym, ("第末相", z1, y1, zm, ym)


def _pump_phase(P, F, sign):
    """
    抽水蓄能第一相式（对象 06~09）：
        Z = ± 2P·[ (P·F² + 1) − F·√(P²F² + 2P + 1) ]
    验算：F=1（全开）→ 括号 = (P+1) − (P+1) = 0 → Z=0；
          F=0（全关）→ 括号 = 1 → Z = ±2P = ±a·v₀/(g·H₀)，即直接水击（知识库 §5.5）。
    """
    inner = (P * F * F + 1.0) - F * math.sqrt(P * P * F * F + 2.0 * P + 1.0)
    return sign * 2.0 * P * inner


def _pump_last(n):
    """
    抽水蓄能第末相式（对象 10~13）：
        Z = 0.5·[ n² + (n⁴ + 4n²)^(1/2) ] = (n/2)(n + √(n²+4))
    （与常规第末相公式同形，n 取 4 个工况特征值之一）。
    """
    return 0.5 * (n * n + math.sqrt(n ** 4 + 4.0 * n * n))


def _pump_branch(params, prof_in, prof_tl, X, Y, Zt, X1, Y1, Z1,
                 VCP, CP, VCP1, CP1, H0):
    """
    抽水蓄能电站（L=2）分支。

    ⚠ **原著无 L=2 算例，本分支的公式转录自 D-19Intro.rtf 的对象 06~29，
    但符号绑定（V0'/U0'、K1 中的 η0/ηP、T0 的斜杠归属）属推断，未经逐位验证。**

    特征值 / 相对值命名按说明书"五、打印成果 (三)2"：
      AN,AN1,AN2,AN3  —— 第一相特征值（发电引水 / 发电尾水 / 抽水引水 / 抽水尾水）
      P0,P1,P2,P3     —— 第末相特征值（同上四工况）
      ZM,ZM1,ZM2,ZM3  —— 第一相相对值（上升 / 下降 / 下降 / 上升）
      F,F1,F2,F3      —— 第末相相对值（同上四工况）
    """
    QT = params["QT"]
    QP = params["QP"]
    TS = params["TS"]
    HR, HM, Z = params["HR"], params["HM"], params["Z"]
    N0, GD2 = params["N0"], params["GD2"]
    HU, HD = params["HU"], params["HD"]

    # 特征值（说明书式 14~17：P = C·V0/(2gH0) 系列）
    V0, U0 = VCP, VCP1
    kq = (QP / QT) if QT else 0.0
    V0p, U0p = V0 * kq, U0 * kq
    AN = CP * V0 / (2.0 * G * H0)          # 发电工况引水管道特征值
    AN1 = CP1 * U0 / (2.0 * G * H0)        # 发电工况尾水管道特征值
    AN2 = CP * V0p / (2.0 * G * H0)        # 抽水工况引水管道特征值
    AN3 = CP1 * U0p / (2.0 * G * H0)       # 抽水工况尾水管道特征值

    # 时间常数（说明书式 21~23）
    WR2 = 0.25 * GD2 * GD2
    K1 = (450.0 * G * HR * QP / (PI * PI * WR2 * N0 * N0 * Z)
          if (WR2 and N0 and Z) else 0.0)
    T0 = 1.39 / 4.5 * K1
    t1_in = 2.0 * X / CP if CP else 0.0
    t1_tl = 2.0 * X1 / CP1 if CP1 else 0.0
    # 开度（说明书式 18~20）：F = F(t)/F0
    Fb_in = (TS - t1_in) / TS if TS else 0.0        # 全甩负荷（引水）
    Fb_tl = (TS - t1_tl) / TS if TS else 0.0        # 全甩负荷（尾水）
    Fc_in = (T0 - t1_in) / T0 if T0 else 0.0        # 水泵断电（引水）
    Fc_tl = (T0 - t1_tl) / T0 if T0 else 0.0        # 水泵断电（尾水）

    z_first = [_pump_phase(AN, Fb_in, +1.0),        # 全甩负荷·高压管道压力上升
               _pump_phase(AN1, Fb_tl, -1.0),       # 全甩负荷·尾水隧洞压力下降
               _pump_phase(AN2, Fc_in, -1.0),       # 水泵断电·高压管道压力下降
               _pump_phase(AN3, Fc_tl, -1.0)]       # 水泵断电·尾水隧洞压力上升

    n0 = X * V0 / (G * H0 * TS) if TS else 0.0
    n1 = X1 * U0 / (G * H0 * TS) if TS else 0.0
    n2 = X * V0 / (G * H0 * T0) if T0 else 0.0
    n3 = X1 * U0 / (G * H0 * T0) if T0 else 0.0
    z_last = [_pump_last(n0), _pump_last(n1), _pump_last(n2), _pump_last(n3)]

    first = abs(z_first[0]) >= abs(z_last[0])
    rel = [abs(v) for v in (z_first if first else z_last)]

    DH, DH1, DH2, DH3 = rel
    base_in = HU - HR                       # 蜗壳中心线以上的静水头
    base_tl = HD - HM                       # 尾水管中心线以下的静水头
    DR = base_in + DH * H0
    DR1 = base_tl - DH1 * H0
    DR2 = base_in - DH2 * H0
    DR3 = base_tl + DH3 * H0

    rows_in = []
    for i in range(params["M"]):
        sh = prof_in["CUM"][i] / Y if Y else 0.0
        dp, dp2 = DH * H0 * sh, DH2 * H0 * sh
        base = HU - params["HT"][i]
        rows_in.append({"i": i + 1, "D": params["D"][i], "dp": dp, "dp2": dp2,
                        "hp": base + dp, "hp2": base - dp2})
    rows_tl = []
    for i in range(params["M1"]):
        sh = prof_tl["CUM"][i] / Y1 if Y1 else 0.0
        dp1, dp3 = DH1 * H0 * sh, DH3 * H0 * sh
        base = HD - params["HN"][i]
        rows_tl.append({"i": i + 1, "D": params["D1"][i], "dp1": dp1, "dp3": dp3,
                        "hp1": base - dp1, "hp3": base + dp3})

    return {
        "相态": "第一相" if first else "第末相",
        "AN": AN, "AN1": AN1, "AN2": AN2, "AN3": AN3,
        "P0": n0, "P1": n1, "P2": n2, "P3": n3,
        "ZM": rel[0], "ZM1": rel[1], "ZM2": rel[2], "ZM3": rel[3],
        "F": rel[0], "F1": rel[1], "F2": rel[2], "F3": rel[3],
        "第一相相对值": [abs(v) for v in z_first],
        "第末相相对值": [abs(v) for v in z_last],
        "WR2": WR2, "K1": K1, "T0": T0,
        "t1": t1_in, "t1_tl": t1_tl,
        "开度": [Fb_in, Fb_tl, Fc_in, Fc_tl],
        "DH": DH, "DH1": DH1, "DH2": DH2, "DH3": DH3,
        "DR": DR, "DR1": DR1, "DR2": DR2, "DR3": DR3,
        "引水行": rows_in, "尾水行": rows_tl,
    }


def compute(params):
    """执行 D-19 计算。返回结构化结果 dict。"""
    M = params["M"]
    M1 = params["M1"]
    L = params["L"]
    C0 = params["C0"]
    EC = params["EC"]
    HU, HD, HU1, HD1 = params["HU"], params["HD"], params["HU1"], params["HD1"]
    QT, QP, TS = params["QT"], params["QP"], params["TS"]
    HR, HM, Z = params["HR"], params["HM"], params["Z"]
    N0, GD2 = params["N0"], params["GD2"]

    prof_in = pipe_profile(params["D"], params["Q"], params["S"],
                           params["H"], params["E"], C0, EC)
    prof_tl = pipe_profile(params["D1"], params["Q1"], params["S1"],
                           params["H1"], params["E1"], C0, EC)

    X, Y, Zt = prof_in["X"], prof_in["Y"], prof_in["Z"]
    VCP, CP = prof_in["VCP"], prof_in["CP"]
    X1, Y1, Z1 = prof_tl["X"], prof_tl["Y"], prof_tl["Z"]
    VCP1, CP1 = prof_tl["VCP"], prof_tl["CP"]
    XX, YY, ZZ = X + X1, Y + Y1, Zt + Z1
    VA = YY / XX if XX else 0.0
    CA = XX / ZZ if ZZ else 0.0

    H0 = HU - HD1

    pump = None
    if L == 1:
        P = CA * VA / (2.0 * G * H0)
        sigma = YY / (G * H0 * TS)
        ZM, ZM1, (phase, z1, y1, zm, ym) = _conventional(prof_in, prof_tl, P, sigma)
    else:
        # 抽水蓄能：管道特征系数仍取全系统口径（与常规站 P 同源）
        P = CA * VA / (2.0 * G * H0)
        sigma = YY / (G * H0 * TS)
        pump = _pump_branch(params, prof_in, prof_tl, X, Y, Zt, X1, Y1, Z1,
                            VCP, CP, VCP1, CP1, H0)
        ZM, ZM1 = pump["ZM"], pump["ZM1"]
        phase = pump["相态"]
        z1 = y1 = zm = ym = None

    DH = ZM * H0
    DH1 = ZM1 * H0
    base_wg = HU - HR
    HMAX = base_wg + DH
    HMIN = base_wg - DH1

    rows = []
    for i in range(M):
        share = (prof_in["CUM"][i] / Y) if Y else 0.0
        dp = DH * share
        dp1 = DH1 * share
        base = HU - params["HT"][i]
        rows.append({"i": i + 1, "D": params["D"][i], "dp": dp, "dp1": dp1,
                     "hp": base + dp, "hp1": base - dp1})

    return {
        "程序": PROGRAM_ID,
        "标题": TITLE,
        "作者": AUTHOR,
        "输入": params,
        "引水": prof_in,
        "尾水": prof_tl,
        "常数": {
            "X": X, "Y": Y, "Z": Zt, "VCP": VCP, "CP": CP,
            "X1": X1, "Y1": Y1, "Z1": Z1, "VCP1": VCP1, "CP1": CP1,
            "XX": XX, "YY": YY, "ZZ": ZZ, "VA": VA, "CA": CA,
            "P": P, "ZM": ZM, "ZM1": ZM1, "H0": H0,
        },
        "相态": phase,
        "第一相": {"Z1": z1, "Y1": y1},
        "第末相": {"Zm": zm, "Ym": ym},
        "抽水蓄能": pump,
        "水击": {"DH": DH, "HMAX": HMAX, "DH1": DH1, "HMIN": HMIN},
        "行": rows,
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(params, result):
    """生成原著风格文本计算书（.OUT）。"""
    r = result
    i = r["输入"]
    pin, ptl = r["引水"], r["尾水"]
    k = r["常数"]
    M, M1, L = i["M"], i["M1"], i["L"]

    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****               引水管道水击压力计算程序 D-19                   ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("                   ORIGINAL DATA")
    A.append("                    原 始 数 据")
    A.append(" *******************************************************")
    A.append("      PROJECT NAME            工程名称代号       N=%d" % i["N"])
    A.append("      PRESSURE PIPE NUMBER    引水管道分段数     M=%3d" % M)
    A.append("      TAILRACE TUNNEL NUMBER  尾水管道分段数    M1=%3d" % M1)
    A.append("      TYPE OF POWER STATION   电站类型信息       L=%3d" % L)
    A.append("      机组转速                   N0=%5d (转／分)" % int(q(i["N0"], 0)))
    A.append("      机组转动惯量              GD2=" + F(i["GD2"], 9, 2) + " (t-m)")
    A.append("      水的传波速度               C0=" + F(i["C0"], 7, 2) + " (m/s)")
    A.append("      水的弹性模量               EC=" + F(i["EC"], 9, 2) + " (kg/cm^2)")
    A.append("      电站上游最高水位           HU=" + F(i["HU"], 7, 2) + " (m)")
    A.append("      电站下游最高水位           HD=" + F(i["HD"], 7, 2) + " (m)")
    A.append("      电站上游死水位            HU1=" + F(i["HU1"], 7, 2) + " (m)")
    A.append("      电站下游死水位            HD1=" + F(i["HD1"], 7, 2) + " (m)")
    A.append("      电站水轮机单机流量         QT=" + F(i["QT"], 9, 2) + " (m^3/s)")
    A.append("      电站水泵单机流量           QP=" + F(i["QP"], 9, 2) + " (m^3/s)")
    A.append("      机组调速器关闭时间         TS=" + F(i["TS"], 7, 2) + " (秒)")
    A.append("      蜗壳中心线高程             HR=" + F(i["HR"], 7, 2) + " (m)")
    A.append("      尾水管中心线高程           HM=" + F(i["HM"], 7, 2) + " (m)")
    A.append("      水轮机效率                  Z=" + F(i["Z"], 7, 3))
    A.append("")
    A.append("")
    A.append("             各 分 段 引 水 管 道 数 据")
    A.append(" -------------------------------------------------------")
    A.append("  I      D       Q        S      H        E        HT")
    A.append(" 段号 管内径   过流量    管长   壁厚     弹模   中心高程")
    A.append("        (m)   (m^3/s)    (m)    (cm)   (kg/cm^2)  (m)")
    A.append(" -------------------------------------------------------")
    for j in range(M):
        A.append("%3d" % (j + 1) + F(i["D"][j], 8, 2) + F(i["Q"][j], 9, 2)
                 + F(i["S"][j], 9, 2) + F(i["H"][j], 7, 2)
                 + FI(i["E"][j], 10) + F(i["HT"][j], 10, 2))
    A.append("")
    A.append("")
    A.append("               各 分 段 尾 水 管 道 数 据")
    A.append(" -------------------------------------------------------")
    A.append(" I1     D1       Q1       S1     H1       E1       HN")
    A.append(" 段号 管内径   过流量    管长   壁厚     弹模   中心高程")
    A.append("        (m)   (m^3/s)    (m)    (cm)   (kg/cm^2)  (m)")
    A.append(" -------------------------------------------------------")
    for j in range(M1):
        A.append("%3d" % (j + 1) + F(i["D1"][j], 8, 2) + F(i["Q1"][j], 9, 2)
                 + F(i["S1"][j], 9, 2) + F(i["H1"][j], 7, 2)
                 + FI(i["E1"][j], 10) + F(i["HN"][j], 10, 2))
    A.append("")
    A.append("")
    A.append("                    RESULT  OF  LV & C")
    A.append("                   ∑LV 及 波 速 C 成 果")
    A.append(" *******************************************************")
    A.append("")
    A.append("")
    A.append("           引 水 管 道 ∑LV 及 波 速 C 成 果")
    A.append(" -------------------------------------------------------")
    A.append("  I       D         V          LV         ∑LV       C")
    A.append(" 段号   管内径    流速                              波速")
    A.append("         (m)      (m/s)      (m^2/s)    (m^2/s)    (m/s)")
    A.append(" -------------------------------------------------------")
    for j in range(M):
        A.append("%3d" % (j + 1) + F(i["D"][j], 10, 2) + F(pin["V"][j], 9, 2)
                 + F(pin["LV"][j], 12, 2) + F(pin["CUM"][j], 12, 2)
                 + F(pin["C"][j], 10, 2))
    A.append("")
    A.append("")
    A.append("            尾 水 管 道 ∑LV 及 波 速 C 成 果")
    A.append(" -------------------------------------------------------")
    A.append("  I       D1       V1         LV1       ∑LV1        C1")
    A.append(" 段号   管内径    流速                              波速")
    A.append("         (m)      (m/s)     (m^2/s)    (m^2/s)     (m/s)")
    A.append(" -------------------------------------------------------")
    for j in range(M1):
        A.append("%3d" % (j + 1) + F(i["D1"][j], 10, 2) + F(ptl["V"][j], 9, 2)
                 + F(ptl["LV"][j], 12, 2) + F(ptl["CUM"][j], 12, 2)
                 + F(ptl["C"][j], 10, 2))
    A.append("")
    A.append("")
    A.append("                   其 它 有 关 成 果")
    A.append(" *******************************************************")
    A.append("    引水管道总长度                X=" + F(k["X"], 9, 2) + " (m)")
    A.append("    管道各点累积∑LV              Y=" + F(k["Y"], 9, 2) + " (m^2/s)")
    A.append("    管长Lx／波速Cx的累加值        Z=" + F(k["Z"], 10, 8) + " (c)")
    A.append("    引水管道加权平均流速        VCP=" + F(k["VCP"], 9, 2) + " (m/s)")
    A.append("    引水管道加权平均波速         CP=" + F(k["CP"], 9, 2) + " (m/s)")
    A.append("    尾水管道总长度               X1=" + F(k["X1"], 9, 2) + " (m)")
    A.append("    尾水管道各点累积∑LV         Y1=" + F(k["Y1"], 9, 2) + " (m^2/s)")
    A.append("    管长／波速的累加值           Z1=" + F(k["Z1"], 10, 8) + " (c)")
    A.append("    尾水隧洞加权平均流速       VCP1=" + F(k["VCP1"], 9, 2) + " (m/s)")
    A.append("    尾水隧洞加权平均波速        CP1=" + F(k["CP1"], 9, 2) + " (m/s)")
    A.append("    引水、尾水洞总长             XX=" + F(k["XX"], 9, 2) + " (m)")
    A.append("    引水、尾水洞总∑LV           YY=" + F(k["YY"], 9, 2) + " (m^2/s)")
    A.append("    引水、尾水洞总∑Li／Ci       ZZ=" + F(k["ZZ"], 10, 8))
    A.append("    引水系统加权平均流速         VA=" + F(k["VA"], 9, 2) + " (m/s)")
    A.append("    引水系统加权平均波速         CA=" + F(k["CA"], 9, 2) + " (m/s)")
    A.append("    管道特征系数                  P=" + F(k["P"], 9, 6))
    A.append("    管道水击压力上升相对值       ZM=" + F(k["ZM"], 9, 6))
    A.append("    管道水击压力下降相对值      ZM1=" + F(k["ZM1"], 9, 6))
    A.append("")
    A.append("")
    A.append("             RESULT OF WATER HUMMER")
    A.append("         水击压力 (水力过渡过程) 计算成果")
    A.append(" *******************************************************")
    A.append("")
    if L == 1:
        A.append("")
        A.append("           常 规 电 站 水 击 压 力 计 算 成 果")
        A.append(" -------------------------------------------------------")
        A.append("    压力引水管末端水击压力上升绝对值     DH=" + F(r["水击"]["DH"], 8, 3) + " (m)")
        A.append("    压力引水管最大压力                 HMAX=" + F(r["水击"]["HMAX"], 8, 3) + " (m)")
        A.append("    压力引水管水击压力下降绝对值        DH1=" + F(r["水击"]["DH1"], 8, 3) + " (m)")
        A.append("    压力引水管最小压力                 HMIN=" + F(r["水击"]["HMIN"], 8, 3) + " (m)")
        A.append("")
        A.append("")
        A.append("              各 特 征 点 上 的 相 应 值")
        A.append(" -------------------------------------------------------")
        A.append("               水击压力    水击压力     对应DP   对应DP1")
        A.append(" 点号  管内径  上升绝对值  下降绝对值   的水头   的水头")
        A.append("  I     D (m)    DP (m)     DP1 (m)     HP (m)   HP1 (m)")
        A.append(" -------------------------------------------------------")
        for row in r["行"]:
            A.append("%3d" % row["i"] + F(row["D"], 9, 3) + F(row["dp"], 11, 3)
                     + F(row["dp1"], 11, 3) + F(row["hp"], 11, 3)
                     + F(row["hp1"], 11, 3))
        A.append("")
    else:
        A.append("[注] 抽水蓄能电站分支依 D-19Intro.rtf 对象 06~29 的公式对象反演实现；")
        A.append("     原著算例目录仅有 L=1（常规电站）算例，本分支**未经逐位验证**：")
        A.append("     V0'/U0'、K1 中的 η0/ηP、T0 的斜杠归属属推断，结果仅供参考。")
        A.append("")
        A.append("           抽 水 蓄 能 电 站 水 击 压 力 计 算 成 果")
        A.append(" -------------------------------------------------------")
        pu = r["抽水蓄能"]
        A.append("    [第一相最大水击压力成果]  (ζ1 > ζm)" if pu["相态"] == "第一相"
                 else "    [第末相最大水击压力成果]  (ζm > ζ1)")
        A.append("      发电工况引水管道特征值       AN=" + F(pu["AN"], 9, 6))
        A.append("      发电工况尾水管道特征值      AN1=" + F(pu["AN1"], 9, 6))
        A.append("      抽水工况引水管道特征值      AN2=" + F(pu["AN2"], 9, 6))
        A.append("      抽水工况尾水管道特征值      AN3=" + F(pu["AN3"], 9, 6))
        A.append("      发电工况引水管道压力上升相对值  ZM=" + F(pu["ZM"], 9, 6))
        A.append("      发电工况尾水管道压力下降相对值 ZM1=" + F(pu["ZM1"], 9, 6))
        A.append("      抽水工况引水管道压力下降相对值 ZM2=" + F(pu["ZM2"], 9, 6))
        A.append("      抽水工况尾水管道压力上升相对值 ZM3=" + F(pu["ZM3"], 9, 6))
        A.append("    第末相特征值 / 相对值：P0=" + F(pu["P0"], 9, 6)
                 + " P1=" + F(pu["P1"], 9, 6) + " P2=" + F(pu["P2"], 9, 6)
                 + " P3=" + F(pu["P3"], 9, 6))
        A.append("      相对值  F=" + F(pu["F"], 9, 6) + " F1=" + F(pu["F1"], 9, 6)
                 + " F2=" + F(pu["F2"], 9, 6) + " F3=" + F(pu["F3"], 9, 6))
        A.append("    常数：K1=" + F(pu["K1"], 9, 6) + "  T0=" + F(pu["T0"], 9, 6)
                 + "  WR^2=" + F(pu["WR2"], 9, 6))
        A.append("          t1=" + F(pu["t1"], 9, 6) + "  t1(尾水)="
                 + F(pu["t1_tl"], 9, 6))
        A.append("")
        A.append("    蜗壳、尾水管出口处的水击计算成果")
        A.append("      发电工况蜗壳最大压力上升相对值      DH=" + F(pu["DH"], 8, 3))
        A.append("      发电工况尾口最大压力下降相对值     DH1=" + F(pu["DH1"], 8, 3))
        A.append("      抽水工况蜗壳最大压力下降相对值     DH2=" + F(pu["DH2"], 8, 3))
        A.append("      抽水工况尾口最大压力上升相对值     DH3=" + F(pu["DH3"], 8, 3))
        A.append("      发电工况蜗壳最大水头                 DR=" + F(pu["DR"], 8, 3) + " (m)")
        A.append("      发电工况尾口最小水头                DR1=" + F(pu["DR1"], 8, 3) + " (m)")
        A.append("      抽水工况蜗壳最小水头                DR2=" + F(pu["DR2"], 8, 3) + " (m)")
        A.append("      抽水工况尾口最大水头                DR3=" + F(pu["DR3"], 8, 3) + " (m)")
        A.append("")
        A.append("              各 特 征 点 上 的 水 击 计 算 成 果")
        A.append(" -------------------------------------------------------")
        A.append("  I0      D        DP        DP2       HP        HP2")
        A.append(" 点号  管内径   水击上升   水击下降   最大水头   最小水头")
        A.append("                (m)       (m)        (m)       (m)")
        A.append(" -------------------------------------------------------")
        for row in pu["引水行"]:
            A.append("%3d" % row["i"] + F(row["D"], 9, 3) + F(row["dp"], 10, 3)
                     + F(row["dp2"], 11, 3) + F(row["hp"], 11, 3)
                     + F(row["hp2"], 11, 3))
        A.append("")
        A.append("  I1      D1       DP1       DP3       HP1       HP3")
        A.append(" 点号  管内径   水击下降   水击上升   最小水头   最大水头")
        A.append("                (m)       (m)        (m)       (m)")
        A.append(" -------------------------------------------------------")
        for row in pu["尾水行"]:
            A.append("%3d" % row["i"] + F(row["D"], 9, 3) + F(row["dp1"], 10, 3)
                     + F(row["dp3"], 11, 3) + F(row["hp1"], 11, 3)
                     + F(row["hp3"], 11, 3))
        A.append("")
    return "\n".join(A)


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
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
