# -*- coding: utf-8 -*-
"""
D-9 调压室水力学计算程序 —— 内核
=================================
复刻《水利水电工程设计计算程序集》D-9 程序（作者：张校正，新疆水利厅）。

原著说明书 `D-9Intro.rtf`（说明书正文 `rtf_text/D-9Intro.txt`）给出五种调压室型式，
其基本微分方程以 **7 个 OLE 公式对象** 内嵌于 RTF，由 `_d9_eq_extract.py`（MTEF v3 令牌流）
与 `_d9_wmf.py`（内嵌 WMF 预览：字符 + 字体 + 坐标）双路解包还原（见 `_d9_eq_formulas.txt` /
`_d9_wmf_formulas.txt`）：

  [00] 差动式（总方程）
       dQ_e/dt = −[ Z_c·g·f/L + φ_c·Q_e²/(2·L·f) ]
       dZ_c/dt = (Q_e − Q − Q_c − Q_B)/F_c
       dZ_p/dt = (Q_c + Q_B)/F_p
  [01] 阻抗孔流量（方向切换流量系数）
       Q_c = M1·F_o·√(2g·(Z_c − Z_p))      或  Q_c = −M2·F_o·√(2g·(Z_p − Z_c))
  [02] 升管顶溢流（堰流，方向切换流量系数）
       Q_B = M3·π·D_c·√(2g)·(Z_c − ZZ)^1.5   或  Q_B = −M4·π·D_c·√(2g)·(ZZ − Z_c)^1.5
  [03] 升管顶以上（Z_c、Z_p 均高于升管顶时）有效面积式： dZ_c/dt = (Q_e − Q)/F_cp
  [04] 阻抗式（圆筒式即 φ_o=0）
       dQ_e/dt = −[ Z·g·f/L + φ_c·Q_e²/(2·L·f) + φ_o·f·(Q_e−Q)²/(2·L·f_o²) ]
       dZ/dt   = (Q_e − Q)/F
  [05] 双室式
       dQ_e/dt = −[ Z·g·f/L + φ_c·Q_e²/(2·L·f) + Σφ·f·(Q_e−Q)²/(2·L·f_0²) ]
       dZ/dt   = (Q_e − Q)/F(z)
  [06] 溢流式（差动式 + 双室式变断面合并）
       dQ_e/dt = −[ Z_c·g·f/L + φ_c·Q_e²/(2·L·f) + Σφ·f·(Q_e−Q)²/(2·L·f_0²) ]
       dZ_c/dt = (Q_e − Q − Q_c − Q_B)/F_c      dZ_p/dt = (Q_c + Q_B)/F_p

其中 f = 引水隧洞断面积 F(输入)、F_c = 升管断面积、F_p = 升管顶以下外室面积、
D_p = 升管顶以上外室面积、F_cp = 升管顶以上调压室有效面积（本内核取 D_p）、
ZZ = 升管溢流口高程 ZS 折算到计算基准面的相对高程。

求解：四阶龙格—库塔（说明书「各种调压室的微分方程均可用四阶龙格－库塔公式求解」）。
步长：常规用大时段 DD；当「升管水位与外室水位恰有一个越过升管顶」时用小时段 TD
（说明书：TD「用于调压井断面突变处，例如上下室与竖井交界处、差动式升管顶部等」）。
注意：小步替换会改变**时间网格**，故输出表的「时段」列是**步号**而非等时距时间。

★ 由权威 OUT 反演并逐位验证的口径（见 d09_verify.py）：
  · 波动周期 T = 2π·√(L·A/(g·f))：g=9.81、π=π（P=1 反演得 136.719286，权威
    136.719279305478，相对差 5e-8；P=3 得 395.3981、P=5 得 177.3184，与权威逐位一致）。
    A 取「最大有效断面」：P=1/P=5 用 D_p（380.134）、P=2/P=3 用 F_c、P=4 用 F_2（竖井）。
  · 初始条件：给定工况 Z_c(0)=Z_p(0)=−K7·Q1²/(2g·F²)、Q_e(0)=Q1；
    丢弃满负荷 −K8·Qo²/(2g·F²)、Q_e(0)=Qo；增加满负荷 Z=0、Q_e=0。
  · 阻抗孔流量 Q_c 与溢流量 Q_B 的**方向切换流量系数**（M1/M2、M3/M4）已被权威 OUT
    的 Q_c/Q_B 列逐位验证（给定工况 Q_c=−93.26 对应 −M2·F_o·√(2gΔZ)，M2=0.8；
    丢弃工况 Q_c=+123.64、Q_B=60.99 对应 M1=0.5、M3=0.33）。
  · 给定工况的 Q(N) 为**节点线性插值**（Q 在第 k 步之值 = 第 k 个数据；两值之间线性过渡，
    稳定流量保持到末时段），由权威 OUT 反演确定（若按每步恒定则残差 1.51，按首值延迟则 3.88）。
  · 输出末行的 Q_c、Q_B 恒为 0.00（原著打印口径：末步之后不再计算孔口/溢流流量）。

未闭合点（如实标注 + 量化）
--------------------------
  1. 给定工况（C=0）与权威 OUT **逐位一致**（31 行 × 3 列，最大偏差 0.01 m³/s = 1 个打印量子；
     水位 0.000 m）；三工况 279 个主列单元 EXACT 132 / NEAR 35 / DECL 112 / **FAIL 0**。
  2. 「突然丢弃满负荷」（C=1）与权威 OUT 存在**最大 0.38 m（水位）/ 6.35 m³/s（流量）**
     的偏差；「突然增加满负荷」（C=2）最大 **0.56 m³/s**。两者共用同一斜坡律
     Q=Q_o·s(To)（说明书「水轮机导叶按直线关闭或直线开启」）与同一套方程，
     即给定工况已把方程、孔口/溢流系数、离散格式全部锁定，剩余差异只能来自
     ①斜坡律的细分口径或 ②升管顶「合并」处理。已做的**数值反证**（见 _d9_probe.py）：
       · 离散格式：RK4 最优（给定工况残差 0.005；Heun 2.29、Euler 3.10）→ 排除低阶格式；
       · 子步加密（DD/2、DD/4、DD/8）：残差反增（0.56→1.05），说明权威解即 DD 步长解，
         而非收敛解 → 排除「程序内部细分」；
       · 步首冻结 Q（每步恒定）：给定工况残差恶化到 2.16 → 排除；
       · 阀特性 √H 律 Q=Q_o·s·√(H/H₀)（H₀ 取 Zo−ZD=220 / Z−ZD=170 / Zo−Z=50）：
         残差恶化到 4.7~20.3 → 排除；
       · 斜坡启动延迟 τ（0.25~1.5 s）：全部恶化（0.56→1.5~12.0）→ 排除；
       · 有效开闭时间 T_e=3.1 s：工况3 残差降至 0.115，但偏离输入 To=3 达 3%，
         与说明书「直线关闭/开启」矛盾 → 不采用；
       · 斜坡指数 s^p：p=1.08 使工况3 残差降至 0.060，同样偏离「直线」表述 → 不采用；
       · 合并（Z_c、Z_p 均高于升管顶）处理：由权威 OUT 的 rows 30→45（水位 6.21→5.37）
         反算有效面积 ≈372 m²，**支持 Fcp = D_p=380.13**（面积口径正确）；再对
         「5 种面积取值 × 3 种切换口径 × 12 组」扫描，取「不切换面积 + Fcp + 体积加权
         强制单水位 + 退出裕度 0.1 m（反向反演常数）」，残差由 8.07 降至 6.35 m³/s。
    结论：残差为**光滑、有限、非发散**（≈3%），且给定工况已逐位命中，属未被反演出的
    某一细分口径，**已量化并排除 8 类假设**，非「不可反演」。
  3. P=2/3/4/5 型式的输出按说明书的微分方程与数据顺序实现；P=3 已用说明书 D-93.OUT
     极值复核（丢弃工况水位 +4.39 / −3.18 m 一致），P=2/P=4/P=5 无权威 OUT，未逐位对拍。

知识库对照结果（教学母本 → 程序实现）
------------------------------------
  KB 路径：D:\\WorkBuddy知识库\\水利知识库\\

  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 连续方程 dZ_c/dt=(Q_e−Q−Q_c−Q_B)/F_c、dZ_p/dt=(Q_c+Q_B)/F_p | `01_水工设计手册精读/卷8_分章/03_第3章_调压设施.md` 式(3.3-8) F·dZ/dt=Q_t−f·u（3.3.4.1） | 形式一致（差一符号约定：手册 Z 向下为正，程序 Z 向上为正） | 采用程序约定 |
  | 动力方程 dQ_e/dt=−(Z_c·g·f/L+φ_c·Q_e²/(2Lf)) | 手册 式(3.3-9) (L/g)·du/dt=Z−h_w−h_c | 一致（同乘 g·f/L 即得） | 采用 |
  | 阻抗孔损失 h_c=φ_o·Q_c²/(2g·f_o²)（方向切换 φ_o） | 手册 式(3.3-11) h_c=Q_c²/(2g·μ²·ω²)，μ 取 0.6~0.8 | 一致（φ_o=1/μ²，f_o=ω）；手册另给方向性：φ_大室→升管 ≥ φ_升管→大室（3.3.6.4） | 采用，方向切换系数由权威 OUT 反演 |
  | 升管顶溢流堰 Q_B=M·π·D_c·√(2g)·H^1.5 | 手册 3.3.6.4 式(3.3-38/39) Q_o=φ_2·ω·√(2g(h_w0+|Z_max|))、溢流层厚 Δh_y=(Q_y/(M·B))^{2/3}（堰流） | 一致（B=π·D_c 为圆形溢流口周长，H=Z_c−ZZ） | 采用（M=M9=.33 / M0=.3） |
  | 升管顶以上有效面积 F_cp | 手册 3.3.6.4「升管顶以上调压井有效面积」；差动式理想设计准则（大室与升管最终同高） | 一致（取 D_p，合并段强制单水位） | 采用 |
  | 波动周期 T=2π√(L·F/(g·f)) | 手册 3.3.5；`02_水利教材精读/水电站输水与厂房/水电站_第4版_精读笔记.md` ch10.6 式10-46 | 一致 | 采用（A 取最大有效断面） |
  | 数值积分（四阶龙格—库塔） | 手册 3.3.6.1 数值法；水电站 ch10.9（式10-55/56 RK4） | 一致 | 采用 |
  | 双室式变断面 F(Z) 分段 | 手册 3.3.6.2 式(3.3-28)（上室 F_u、竖井 F_s）、Jaeger 双室式(3.3-16) | 一致 | 采用（P=4/P=5） |
  | 托马稳定断面 F_T=L·f/(2g·α·H₀) | 手册 式(3.3-13) | 一致 | 本程序不输出，仅记录 |
  | 最高/最低涌波解析式（式3.3-20~27、3.3-37） | 手册 3.3.6.2 | 程序**未用解析式**（全过程数值解） | 记录，不实现 |

知识产权声明：本改造工作基于公之于众版（乌鲁木齐正海水利科技有限公司，张校正教授级高工
技术总负责）的公开算法。改造实现（Python 代码、架构设计、验证数据、自动化流程）为哈胜的原创成果。
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-9"
TITLE = "调压室水力学计算程序"
AUTHOR = "张校正"
HEAD_NAME = "D-9"

G = 9.81
PI = math.pi          # 独立反演：见 docstring（P=1/P=3/P=5 周期逐位一致）

MODES = {1: "差动式", 2: "阻抗式", 3: "园筒式", 4: "双室式", 5: "溢流式"}


# ============================================================
# 打印取整（VB6 Format 口径：十进制四舍五入，半值进位）
# ============================================================

def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    return "%*.*f" % (w, d, q(x, d))


def FI(x, w):
    return "%*d" % (w, int(q(x, 0)))


def vbfmt(x):
    """VB6 Format(x, "0.##") 风格：最多两位小数、去尾零、无前导 0（如 .5 / 38.47 / 1425.5）。"""
    s = "%.2f" % q(x, 2)
    if s.endswith("0"):
        s = s[:-1]
    if s.endswith("."):
        s = s[:-1]
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return s


# ============================================================
# 解析（按 P 分派，数据顺序见说明书「三、输入数据准备」）
# ============================================================

def parse(data):
    if isinstance(data, dict):
        p = dict(data)
        p["P"] = int(round(p.get("P", 1)))
        p["M"] = int(round(p["M"]))
        return p
    n = read_numbers(data)
    if len(n) < 6:
        raise ValueError("D-9 输入数据过少（%d 个数）" % len(n))
    P = int(round(n[0]))
    p = {"P": P, "M": int(round(n[1])), "DD": n[2], "TD": n[3]}
    i = 4
    p["F"], p["K7"], p["K8"], p["K9"], p["L"] = n[i:i + 5]
    i += 5
    if P == 1:
        p["Fc"], p["D"], p["Zs"] = n[i:i + 3]
        i += 3
        p["Fp"], p["Dp"], p["Fo"] = n[i:i + 3]
        i += 3
        p["M1"], p["M2"], p["M3"], p["M4"] = n[i:i + 4]
        i += 4
    elif P == 2:
        p["Fc"], p["Fo"], p["M1"], p["M2"] = n[i:i + 4]
        i += 4
    elif P == 3:
        p["Fc"] = n[i]
        i += 1
    elif P == 4:
        p["F1"], p["F2"], p["F3"] = n[i:i + 3]
        i += 3
        p["Z1"], p["Z2"], p["Z3"] = n[i:i + 3]
        i += 3
        (p["M1"], p["M2"], p["M3"], p["M4"], p["M5"], p["M6"]) = n[i:i + 6]
        i += 6
    elif P == 5:
        p["Fc"], p["D"], p["Zs"] = n[i:i + 3]
        i += 3
        p["Fp"], p["Dp"], p["Fo"] = n[i:i + 3]
        i += 3
        (p["M7"], p["M8"], p["M9"], p["M0"]) = n[i:i + 4]
        i += 4
        p["F1"], p["F2"], p["F3"] = n[i:i + 3]
        i += 3
        p["Z1"], p["Z2"], p["Z3"] = n[i:i + 3]
        i += 3
        (p["M1"], p["M2"], p["M3"], p["M4"], p["M5"], p["M6"]) = n[i:i + 6]
        i += 6
    else:
        raise ValueError("D-9 未知调压室型式 P=%d" % P)
    p["Qo"], p["To"] = n[i:i + 2]
    i += 2
    p["Zo"], p["Z"], p["ZD"] = n[i:i + 3]
    i += 3
    p["A"] = int(round(n[i]))
    i += 1
    p["Qlist"] = n[i:i + p["A"]]
    if len(p["Qlist"]) < p["A"]:
        raise ValueError("D-9 引水流量数据不足（需 %d 个）" % p["A"])
    return p


# ============================================================
# 中间成果常数
# ============================================================

def tank_area(p):
    """波动周期用的「最大有效调压室断面」。"""
    P = p["P"]
    if P in (1, 5):
        return p["Dp"]
    if P in (2, 3):
        return p["Fc"]
    return p["F2"]          # P=4 竖井


def constants(p):
    L, f = p["L"], p["F"]
    A = tank_area(p)
    TT = 2.0 * PI * math.sqrt(L * A / (G * f))
    return {"T": TT, "A_tank": A, "F": f, "L": L}


# ============================================================
# 计算核心
# ============================================================

def _riser_eq(p, Ks, datums, areaswitch, snap, amerg=None, relrm=0.0,
              mnames=("M1", "M2", "M3", "M4"), only=None):
    """
    差动式/溢流式（含升管 + 外室两水位）公共动力学。
    Ks: {case: 沿程+局部损失系数}；datums: {case: 计算基准面高程}。
    返回 dict：Zc/Zp/Qe 序列、小步清单。
    """
    F, L, Fc, Fp, Dp, Fo, D = p["F"], p["L"], p["Fc"], p["Fp"], p["Dp"], p["Fo"], p["D"]
    M, DD, TD = p["M"], p["DD"], p["TD"]
    sp2g = math.sqrt(2.0 * G)

    def ramp_s(t, case):
        To = p["To"]
        if case == "drop":
            return max(0.0, 1.0 - t / To)
        if case == "add":
            return max(0.0, t / To)
        return None

    def Qt(t, case):
        if case == "given":
            Qlist = p["Qlist"]
            k = t / DD
            if k >= len(Qlist) - 1:
                return Qlist[-1]
            i = int(k)
            return Qlist[i] + (Qlist[i + 1] - Qlist[i]) * (k - i)
        return p["Qo"] * min(1.0, ramp_s(t, case))

    def flow(Zc, Zp, Zs):
        dz = Zc - Zp
        if dz >= 0:
            Qc = p[mnames[0]] * Fo * math.sqrt(2.0 * G * dz)
        else:
            Qc = -p[mnames[1]] * Fo * math.sqrt(2.0 * G * (-dz))
        if Zc > Zs and Zp > Zs:            # 两水位均高于升管顶：合并，无溢流
            QB = 0.0
        elif Zc > Zs:
            QB = p[mnames[2]] * PI * D * sp2g * (Zc - Zs) ** 1.5
        elif Zp > Zs:
            QB = -p[mnames[3]] * PI * D * sp2g * (Zp - Zs) ** 1.5
        else:
            QB = 0.0
        return Qc, QB

    def deriv(Zc, Zp, Qe, t, case, Zs, K):
        Q = Qt(t, case)
        Qc, QB = flow(Zc, Zp, Zs)
        mg = (Zc > Zs and Zp > Zs)
        if areaswitch == "dp":
            Ac, Ap = (Dp, Dp) if mg else (Fc, Fp)
        elif areaswitch == "riser":
            Ac, Ap = (Dp, Fp) if mg else (Fc, Fp)
        else:                               # none
            Ac, Ap = Fc, Fp
        dQe = -(G * F / L) * Zc - K * Qe * abs(Qe) / (2.0 * L * F)
        return dQe, (Qe - Q - Qc - QB) / Ac, (Qc + QB) / Ap

    def small_step(Zc, Zp, Zs):
        return (Zc > Zs) != (Zp > Zs)

    res = {}
    for case in ("given", "drop", "add"):
        if case == "given" and not p["Qlist"]:
            continue
        if only is not None and case not in only:
            continue
        K = Ks[case]
        datum = datums[case]
        Zs = p["Zs"] - datum
        if case == "given":
            Q0 = p["Qlist"][0] if p["Qlist"] else 0.0
            Z0 = -K * Q0 * Q0 / (2.0 * G * F * F)
            Qe = Q0
        elif case == "drop":
            Z0 = -K * p["Qo"] * p["Qo"] / (2.0 * G * F * F)
            Qe = p["Qo"]
        else:
            Z0 = 0.0
            Qe = 0.0
        Zc = Zp = Z0
        t = 0.0
        rows = [(0, Zc, Zp, Qe)]
        small = []
        merged = False
        for k in range(1, M + 1):
            if merged:
                # 合并段（两水位同高、无孔口/溢流交换）：单水位，有效面积 Am
                Am = (Dp if amerg is None else amerg)
                h = DD
                def dm(zz, qe, tt):
                    q = Qt(tt, case)
                    return (-(G * F / L) * zz - K * qe * abs(qe) / (2.0 * L * F),
                            (qe - q) / Am)
                k1 = dm(Zc, Qe, t)
                k2 = dm(Zc + h / 2 * k1[1], Qe + h / 2 * k1[0], t + h / 2)
                k3 = dm(Zc + h / 2 * k2[1], Qe + h / 2 * k2[0], t + h / 2)
                k4 = dm(Zc + h * k3[1], Qe + h * k3[0], t + h)
                Zc += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
                Qe += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
                Zp = Zc
                t += h
                if Zc <= Zs - relrm:         # 退出合并：保持同高，恢复常规方程
                    merged = False
                if k % 5 == 0 or k == M:
                    rows.append((k, Zc, Zp, Qe))
                continue
            h = TD if small_step(Zc, Zp, Zs) else DD
            if h == TD:
                small.append(k)
            k1 = deriv(Zc, Zp, Qe, t, case, Zs, K)
            k2 = deriv(Zc + h / 2 * k1[1], Zp + h / 2 * k1[2], Qe + h / 2 * k1[0], t + h / 2, case, Zs, K)
            k3 = deriv(Zc + h / 2 * k2[1], Zp + h / 2 * k2[2], Qe + h / 2 * k2[0], t + h / 2, case, Zs, K)
            k4 = deriv(Zc + h * k3[1], Zp + h * k3[2], Qe + h * k3[0], t + h, case, Zs, K)
            Zc += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            Zp += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
            Qe += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            t += h
            if Zc > Zs and Zp > Zs:          # 进入合并：体积加权强制单水位
                v = (Fc * Zc + Dp * Zp) / (Fc + Dp)
                Zc = Zp = v
                merged = True
            if k % 5 == 0 or k == M:
                rows.append((k, Zc, Zp, Qe))
        res[case] = {"rows": rows, "small": small, "datum": datum, "K": K}
    return res


def _single_eq(p, Ks, datums, phi_laws, only=None):
    """
    阻抗式/圆筒式/双室式（单水位）公共动力学。
    phi_laws: callable(Qe, Q) -> φ（阻抗系数，按流向切换）
    """
    F, L, M, DD, TD = p["F"], p["L"], p["M"], p["DD"], p["TD"]
    Fc = p.get("Fc", p.get("F2", 1.0))
    disp = p.get("_disp")                     # 双室式断面函数
    Qo, To = p["Qo"], p["To"]
    Qlist = p["Qlist"]

    def area_fun(Z):
        return disp(Z) if disp else Fc

    def Qt(t, case):
        if case == "given":
            k = t / DD
            if k >= len(Qlist) - 1:
                return Qlist[-1]
            i = int(k)
            return Qlist[i] + (Qlist[i + 1] - Qlist[i]) * (k - i)
        if case == "drop":
            return Qo * min(1.0, max(0.0, 1.0 - t / To))
        return Qo * min(1.0, max(0.0, t / To))

    def deriv(Z, Qe, t, case, K):
        Q = Qt(t, case)
        Qc = Qe - Q
        phi = phi_laws(Qc)
        f0 = p.get("f0", 0.0)
        ori = (phi * F * Qc * Qc / (2.0 * L * f0 * f0)) if f0 else 0.0
        dQe = -(G * F / L) * Z - K * Qe * abs(Qe) / (2.0 * L * F) - ori
        return dQe, Qc / area_fun(Z)

    def cross(Z, Zprev):
        if disp is None:
            return False
        for e in p.get("_edges", []):
            if (Zprev - e) * (Z - e) < 0:
                return True
        return False

    res = {}
    for case in ("given", "drop", "add"):
        if case == "given" and not Qlist:
            continue
        if only is not None and case not in only:
            continue
        K = Ks[case]
        if case == "given":
            Q0 = Qlist[0] if Qlist else 0.0
            Z0 = -K * Q0 * Q0 / (2.0 * G * F * F)
            Qe = Q0
        elif case == "drop":
            Z0 = -K * Qo * Qo / (2.0 * G * F * F)
            Qe = Qo
        else:
            Z0 = 0.0
            Qe = 0.0
        Z = Z0
        t = 0.0
        rows = [(0, Z, Z, Qe)]
        small = []
        for k in range(1, M + 1):
            h = DD
            k1 = deriv(Z, Qe, t, case, K)
            k2 = deriv(Z + h / 2 * k1[1], Qe + h / 2 * k1[0], t + h / 2, case, K)
            k3 = deriv(Z + h / 2 * k2[1], Qe + h / 2 * k2[0], t + h / 2, case, K)
            k4 = deriv(Z + h * k3[1], Qe + h * k3[0], t + h, case, K)
            Zp0 = Z
            Z += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            Qe += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            t += h
            if cross(Z, Zp0):
                small.append(k)
            if k % 5 == 0 or k == M:
                rows.append((k, Z, Z, Qe))
        res[case] = {"rows": rows, "small": small, "datum": datums[case], "K": K}
    return res


def compute(p):
    P = p["P"]
    k = constants(p)
    # 误差系数：给定工况 K7、突弃 K8、突增 K9（说明书「给定工况/突然丢弃/突然增加」三档）
    Ks = {"given": p["K7"], "drop": p["K8"], "add": p["K9"]}
    # 基准面：丢弃满负荷 → 水库最高水位 Zo；增加满负荷 → 水库最低水位 Z；
    # 给定工况 → 权威 OUT 反演为水库最低水位 Z（增荷型工况）
    datums = {"given": p["Z"], "drop": p["Zo"], "add": p["Z"]}
    if P in (1, 5):
        if P == 5:
            # 溢流式：突然丢弃满负荷用「升管 + 外室」模型（权威版式为 7 列）；
            # 给定/增加工况用「单水位 + 变断面」模型（权威版式为 4 列两栏）
            p["f0"] = p["F2"]
            p["_edges"] = [p["Z1"], p["Z2"], p["Z3"]]

            def disp5(Z):
                if Z > p["Z1"]:
                    return p["F1"]
                if Z > p["Z2"]:
                    return p["F2"]
                return p["F3"]

            def phiw(qc):
                if qc >= 0:
                    return p["M1"] + p["M3"] + p["M5"]
                return p["M2"] + p["M4"] + p["M6"]

            p["_disp"] = disp5
            sim = _single_eq(p, Ks, datums, phiw, only=("given", "add"))
            rd = _riser_eq(p, Ks, datums, "none", True, relrm=0.1,
                           mnames=("M7", "M8", "M9", "M0"), only=("drop",))
            sim.update(rd)
        else:
            sim = _riser_eq(p, Ks, datums, "none", True, relrm=0.1)
    elif P == 3:
        p["f0"] = 0.0
        sim = _single_eq(p, Ks, datums, lambda qc: 0.0)
    elif P == 2:
        p["f0"] = p["Fo"]
        sim = _single_eq(p, Ks, datums, lambda qc: p["M1"] if qc >= 0 else p["M2"])
    else:                                   # P=4 双室式
        p["f0"] = p["F2"]
        p["_edges"] = [p["Z1"], p["Z2"], p["Z3"]]

        def disp(Z):
            if Z > p["Z1"]:
                return p["F1"]
            if Z > p["Z2"]:
                return p["F2"]
            return p["F3"]

        def philaw(qc):
            if qc >= 0:
                return p["M1"] + p["M3"] + p["M5"]
            return p["M2"] + p["M4"] + p["M6"]

        p["_disp"] = disp
        sim = _single_eq(p, Ks, datums, philaw)
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR, "输入": p,
            "常数": k, "模拟": sim, "基准": datums}


# ============================================================
# 输出（复刻原著 .OUT 版式，GBK）
# ============================================================

def _small_block(small):
    A = [" 下列时段按小时段间隔计算:"]
    if small:
        for i in range(0, len(small), 20):
            A.append("".join("%4d" % x for x in small[i:i + 20]))
    A.append(" 共计 %d  个小时段" % len(small))
    A.append("")
    return A


def render(p, r):
    P = p["P"]
    sim = r["模拟"]
    M = p["M"]
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" *****                  调压室水力计算书 D-9X                      *****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("                        工程名:工程名")
    A.append("")
    A.append("")
    name = MODES[P]
    A.append("                     %s调压室水力计算" % name)
    A.append("                     ====================")
    A.append("")
    A.append("                      二. 基  本  数  据")
    A.append("")
    A.append("                           P= %d " % P)
    A.append(" 计算时段数:                              M= %d " % M)
    A.append(" 大的时段长(秒):                         DD= %s " % vbfmt(p["DD"]))
    A.append(" 小的时段长(秒):                         TD= %s " % vbfmt(p["TD"]))
    A.append(" 引水隧洞断面面积:                        F= %s " % vbfmt(p["F"]))
    A.append(" 引水隧洞沿程损失与局部损失系数:")
    A.append("     给定工况的沿程局部损失系数之和:     K7= %s " % vbfmt(p["K7"]))
    A.append("       突然丢弃满负荷的损失系数之和:     K8= %s " % vbfmt(p["K8"]))
    A.append("       突然增加满负荷的损失系数之和:     K9= %s " % vbfmt(p["K9"]))
    A.append(" 引水隧洞长度(米):                        L= %s " % vbfmt(p["L"]))
    if P == 1:
        A.append(" 升管断面面积(平方米):                   Ff= %s " % vbfmt(p["Fc"]))
        A.append(" 升管溢流口直径(米):                      D= %s " % vbfmt(p["D"]))
        A.append(" 升管溢流口高程(米):                     ZS= %s " % vbfmt(p["Zs"]))
        A.append(" 升管顶以下之外室断面面积(平方米):       Fp= %s " % vbfmt(p["Fp"]))
        A.append(" 升管顶以上之外室断面面积(平方米):       Dp= %s " % vbfmt(p["Dp"]))
        A.append(" 阻抗孔口断面面积(平方米):               Fo= %s " % vbfmt(p["Fo"]))
        A.append(" 水流通过阻抗孔口流入外室时的流量系数:   M7= %s " % vbfmt(p["M1"]))
        A.append(" 水流通过阻抗孔口流入升管时的流量系数:   M8= %s " % vbfmt(p["M2"]))
        A.append(" 水流溢入外室时的流量系数:               M9= %s " % vbfmt(p["M3"]))
        A.append(" 水流溢入升管时的流量系数:               Mo= %s " % vbfmt(p["M4"]))
    elif P == 2:
        A.append(" 调压室断面面积(平方米):                  FC= %s " % vbfmt(p["Fc"]))
        A.append(" 阻抗孔口断面面积(平方米):               Fo= %s " % vbfmt(p["Fo"]))
        A.append(" 水流通过阻抗孔口流入调压井的阻抗系数:   M1= %s " % vbfmt(p["M1"]))
        A.append(" 水流通过阻抗孔口流出调压井的阻抗系数:   M2= %s " % vbfmt(p["M2"]))
    elif P == 3:
        A.append(" 调压室断面面积(平方米):                  FC= %s " % vbfmt(p["Fc"]))
    elif P == 4:
        A.append(" 调压室上室断面面积(平方米):             F1= %s " % vbfmt(p["F1"]))
        A.append(" 调压室竖井断面面积(平方米):             F2= %s " % vbfmt(p["F2"]))
        A.append(" 调压室下室断面面积(平方米):             F3= %s " % vbfmt(p["F3"]))
        A.append(" 调压室上室底高程(米):                   Z1= %s " % vbfmt(p["Z1"]))
        A.append(" 调压室下室顶高程(米):                   Z2= %s " % vbfmt(p["Z2"]))
        A.append(" 调压室下室底高程(米):                   Z3= %s " % vbfmt(p["Z3"]))
    else:
        A.append(" 升管断面面积(平方米):                   Ff= %s " % vbfmt(p["Fc"]))
        A.append(" 升管溢流口直径(米):                      D= %s " % vbfmt(p["D"]))
        A.append(" 升管溢流口高程(米):                     ZS= %s " % vbfmt(p["Zs"]))
        A.append(" 升管顶以下之外室断面面积(平方米):       Fp= %s " % vbfmt(p["Fp"]))
        A.append(" 升管顶以上之外室断面面积(平方米):       Dp= %s " % vbfmt(p["Dp"]))
        A.append(" 阻抗孔口断面面积(平方米):               Fo= %s " % vbfmt(p["Fo"]))
        A.append(" 水流通过阻抗孔口流入外室时的流量系数:   M7= %s " % vbfmt(p["M7"]))
        A.append(" 水流通过阻抗孔口流入升管时的流量系数:   M8= %s " % vbfmt(p["M8"]))
        A.append(" 水流溢入外室时的流量系数:               M9= %s " % vbfmt(p["M9"]))
        A.append(" 水流溢入升管时的流量系数:               Mo= %s " % vbfmt(p["M0"]))
        A.append(" 调压室上室断面面积(平方米):             F1= %s " % vbfmt(p["F1"]))
        A.append(" 调压室竖井断面面积(平方米):             F2= %s " % vbfmt(p["F2"]))
        A.append(" 调压室下室断面面积(平方米):             F3= %s " % vbfmt(p["F3"]))
        A.append(" 调压室上室底高程(米):                   Z1= %s " % vbfmt(p["Z1"]))
        A.append(" 调压室下室顶高程(米):                   Z2= %s " % vbfmt(p["Z2"]))
        A.append(" 调压室下室底高程(米):                   Z3= %s " % vbfmt(p["Z3"]))
    A.append(" 电站满负荷时引水隧洞流量(立方米/秒):     Qo= %s " % vbfmt(p["Qo"]))
    A.append(" 水轮机导水叶关闭时间(秒):               To= %s " % vbfmt(p["To"]))
    A.append(" 水库最高水位(米):                        Zo= %s " % vbfmt(p["Zo"]))
    A.append(" 水库最低水位(米):                        Z= %s " % vbfmt(p["Z"]))
    A.append(" 水轮机导水叶安装高程(米):               ZD= %s " % vbfmt(p["ZD"]))
    A.append(" 电站引水流量数据个数:                    A= %d " % p["A"])
    if p["A"] > 0:
        A.append(" 电站引水流量(初始流量到稳定流量 立方米/秒):")
        qs = p["Qlist"]
        seg = vbfmt(qs[0]) + "  " + "  ".join(vbfmt(x) for x in qs[1:2])
        for x in qs[2:]:
            seg += " , " + vbfmt(x)
        seg += " , " + vbfmt(qs[-1]) + " ,......"
        A.append("                                       Q(N)= " + seg)
    A.append("")
    A.append("                      三. 计  算  结  果:")
    A.append("")

    # ---- 各工况结果块 ----
    def qc_qb(p, Zc, Zp, datum):
        Fc, Fp, Dp, Fo, D = p["Fc"], p["Fp"], p["Dp"], p["Fo"], p["D"]
        Zs = p["Zs"] - datum
        dz = Zc - Zp
        if dz >= 0:
            Qc = p["M1"] * Fo * math.sqrt(2.0 * G * dz)
        else:
            Qc = -p["M2"] * Fo * math.sqrt(2.0 * G * (-dz))
        if Zc > Zs and Zp > Zs:
            QB = 0.0
        elif Zc > Zs:
            QB = p["M3"] * PI * D * math.sqrt(2.0 * G) * (Zc - Zs) ** 1.5
        elif Zp > Zs:
            QB = -p["M4"] * PI * D * math.sqrt(2.0 * G) * (Zp - Zs) ** 1.5
        else:
            QB = 0.0
        return Qc, QB

    twocol = P in (2, 3, 4)

    def block(tag, title, sub, rows, datum, two_col):
        B = []
        B.append("                    " + title)
        B.append("                   " + sub)
        rws = rows
        # 极值统计
        iqmax = max(range(len(rws)), key=lambda i: rws[i][3])
        iqmin = min(range(len(rws)), key=lambda i: rws[i][3])
        B.append(" 引水隧洞最大流量:          Qe(%2d)=%9.2f" % (rws[iqmax][0], rws[iqmax][3]))
        B.append(" 引水隧洞最小流量:          Qe(%2d)=%9.2f" % (rws[iqmin][0], rws[iqmin][3]))
        if P in (1, 5) and not two_col:
            izc = max(range(len(rws)), key=lambda i: rws[i][1])
            izcl = min(range(len(rws)), key=lambda i: rws[i][1])
            izp = max(range(len(rws)), key=lambda i: rws[i][2])
            izpl = min(range(len(rws)), key=lambda i: rws[i][2])
            B.append(" 升管最高水位:              Zc(%2d)=%9.2f" % (rws[izc][0], rws[izc][1]))
            B.append(" 升管最低水位:              Zc(%2d)=%9.2f" % (rws[izcl][0], rws[izcl][1]))
            B.append(" 外室最高水位:              Zp(%2d)=%9.2f" % (rws[izp][0], rws[izp][2]))
            B.append(" 外室最低水位:              Zp(%2d)=%9.2f" % (rws[izpl][0], rws[izpl][2]))
            hi = max(rws[izc][1], rws[izp][2])
            lo = min(rws[izcl][1], rws[izpl][2])
            if tag == "drop":
                B.append(" 调压室最高水位高程:%28.2f" % (datum + hi))
            else:
                B.append(" 调压室最低水位高程:%28.2f" % (datum + lo))
        else:
            iz = max(range(len(rws)), key=lambda i: rws[i][1])
            izl = min(range(len(rws)), key=lambda i: rws[i][1])
            B.append(" 调压室最高水位:                   Z(%2d)=%9.2f" % (rws[iz][0], rws[iz][1]))
            B.append(" 调压室最低水位:                   Z(%2d)=%9.2f" % (rws[izl][0], rws[izl][1]))
            if tag == "drop":
                B.append(" 调压室最高水位高程:%28.2f" % (datum + rws[iz][1]))
            else:
                B.append(" 调压室最低水位高程:%28.2f" % (datum + rws[izl][1]))
        B.append("")
        B.append("")
        B.append("                         调压室水位波动过程")
        B.append("                         ==================")
        B.append("                              周期 T= %.12f " % r["常数"]["T"])
        B.append("                         大时段长 DD= %s " % vbfmt(p["DD"]))
        B.append("                         小时段长 TD= %s " % vbfmt(p["TD"]))
        B.append("")
        if P in (1, 5) and not two_col:
            B.append(" 时段     升管水位    外室水位   隧洞流量    隧洞流速   阻抗孔流量     溢流量")
            for i, (idx, Zc, Zp, Qe) in enumerate(rws):
                if i == len(rws) - 1:
                    Qc = QB = 0.0
                else:
                    Qc, QB = qc_qb(p, Zc, Zp, datum)
                B.append(FI(idx, 5) + F(Zc, 11, 2) + F(Zp, 11, 2) + F(Qe, 10, 2)
                         + F(Qe / p["F"], 11, 2) + F(Qc, 11, 2) + F(QB, 11, 2))
        else:
            head = " 时段    水位   隧洞流量  隧洞流速             时段    水位   隧洞流量  隧洞流速"
            B.append(head)
            half = len(rws) // 2
            left, right = rws[:half], rws[half:]
            for i in range(half):
                li = left[i]
                sL = FI(li[0], 5) + F(li[1], 9, 2) + F(li[3], 10, 2) + F(li[3] / p["F"], 10, 2)
                if i < len(right):
                    ri = right[i]
                    sR = "".join(" " * 8) + FI(ri[0], 5) + F(ri[1], 9, 2) + F(ri[3], 10, 2) + F(ri[3] / p["F"], 10, 2)
                else:
                    sR = ""
                B.append(sL + sR)
        B.append("")
        B += _small_block(sim[tag]["small"])
        return B

    orders = []
    if p["A"] > 0:
        orders.append(("given", "给定工况之计算结果", "=================="))
    orders.append(("drop", "电站突然丢弃满负荷计算结果", "=========================="))
    orders.append(("add", "电站突然增加满负荷计算结果", "=========================="))
    for tag, title, sub in orders:
        two_col = twocol or (P == 5 and tag == "given")
        A += block(tag, title, sub, sim[tag]["rows"], r["基准"][tag], two_col)
    return "\n".join(A)


def run(data, out_txt=None, out_json=None):
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
