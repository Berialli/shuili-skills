# -*- coding: utf-8 -*-
"""
D-23 小波动过渡过程计算程序 —— 内核
=====================================
复刻《水利水电工程设计计算程序集》D-23（作者：姚廉华）。
原著说明书明示"计算方法参照王树人《调压井水力计算理论与方法》…"，
其「二、计算方法及公式」仅以 4 个 OLE 公式对象给出 **T1/T2/FTH1/FTH2** 四式
（由 `_d23_mtef.py`（MTEF v3 令牌流）+ `_d23_wmf.py`（内嵌 WMF 预览字符+坐标）
双路还原，见 `_d23_mtef_tokens.txt` / `_d23_wmf_formulas.txt`）：

  [00] T1 = 2π·√( (L1/g)·(F1/f1) )        [MTEF 令牌流可辨 FRAC/ROOT 结构]
  [01] T2 = 2π·√( (L2/g)·(F2/f2) )
  [02] FTH1 = …（含 q0(θ0)、L1、f1、2g、H、ξ1、V10²、F1 等符号）
  [03] FTH2 = …（同构，下标换成 2）

其余打印量（ZS/H/P/A/B/AA/BB/AL/AN/R）说明书未给式，本内核**由权威
`D-23-1.OUT` 唯一算例逐位反演**得到（见下「★已精确复现」与「★未闭合」）。

★已精确复现（权威算例 1111：D1=D2=9, D01=D02=14, S1=948.108, S2=2020.257,
Q0=273, H0=541.8, C1=1.8577, C2=3.1137）
--------------------------------------------------------------------
  f_i = 0.7854·D_i²（隧洞面积）；F_i = 0.7854·D0_i²（调压井面积）
  V_i0 = Q0/f_i；HN = H0 − C1 − C2（净水头 536.8286）
  ω_i = √(g·f_i/(L_i·F_i))；T_i = 2π/ω_i        → 96.086 / 140.260 ✓
  TT = T1/T2                                    → 0.685 ✓
  ZS_i = Q0/(F_i·ω_i)                           → 27.120 / 39.589 ✓
  H_i  = HN/ZS_i                                → 19.794 / 13.560 ✓
  P_i  = C_i/ZS_i                               → 0.068 / 0.079 ✓
  n_i  = 2·g·C_i·HN·F_i²/(L_i·f_i²·V_i0²)       → 6.562 / 5.161 ✓
  FTH_i = F_i/n_i                               → 23.460 / 29.825 ✓
  k_i = 2·g·C_i/(L_i·V_i0)；cq = Q0/HN；c_i = cq/F_i
  A_i = k_i − c_i                               → 5.65484805E-03 / 3.74310188E-03 ✓(9位)
  B_i = ω_i² − k_i·c_i                          → 4.24642966E-03 / 1.98346215E-03 ✓(9位)
  AA_i = −k_i/C_i                               → −4.82230237E-03 / −2.26310982E-03 ✓
  BB_i = −k_i²/C_i                              → −4.32000729E-05 / −1.59473324E-05 ✓
  特征四次方程 x⁴+AL1x³+AL2x²+AL3x+AL4=0：
    AL1 = A1 + A2                               → 9.39794993E-03 ✓(精确)
    AL4 = B1·B2 − BB1·BB2                        → 8.42194357E-06 ✓(精确)
    R   = −(AL1·AL2·AL3 − AL3² − AL1²·AL4)      （Routh–Hurwitz 行列式，符号相反）
  （A_i/B_i 恰为 4×4 状态矩阵两个 2×2 对角块的特征多项式系数：
    x²+A1x+B1 = det(xI−M11)，M11=[[−k1,−g/L1],[f1/F1,c1]]…由此 AL1=tr 关系成立。）

★未闭合（AL2/AL3 的修正项，已量化）
----------------------------------
  由 4×4 分块矩阵 [[M11,M12],[M21,M22]] 展开，特征多项式 = p(x)q(x) − R(x)：
    α1 = A1+A2（与耦合无关，恒成立）
    α2 = B1+B2+A1A2 − X,   α3 = A1B2+A2B1 − Y,   α4 = B1B2 − Z
  权威值 ⇒ X = +4.1952387576E-06, Y = +9.46708768E-06, Z = +6.8893090915E-10。
  本内核取分块模型值 X0 = c1·c2, Y0 = c1·c2·(k1+k2), Z0 = c1·c2·k1·k2：
    X0 = 1.0913358E-05（偏大 2.60 倍 → α2 偏低 0.24%）
    Y0 = 1.7466E-07（偏小 54 倍 → α3 偏低 26.4%）
    Z0 = 6.8896E-10（与 Z 相对差 4.4E-05 → α4 精确）
  已系统排除（数值反证）：① 耦合项符号全组合（8 组）；② 水头含损失变化项
  (2C_i/V_i0)·v_i（4 组符号）；③ 恒定开度 ΔQ=Q0ΔH/(2HN)（α1 差 35%）；
  ④ 十字块 4 参数反解（α2/α3/α4 仅由 a·c、b·c 两个量决定，三方程超定无解）；
  ⑤ A/B/AA/BB 间的 2 项乘积组合全搜索（351² 组，X/Y 均无 1E-8 级命中）。
  根因：说明书未给 AL2/AL3 显式式，缺王树人原书 α1~α4 表达式；本内核按
  分块模型实现，**α2/α3 记为 DECL**（偏差 0.24% / 26.4%，已在 verify 中量化）。

常量口径（本程序独立反演）
--------------------------
  g = 9.81（由 AL1 = A1+A2 的 k_i/ω_i 逐位命中确认；g=9.8 时 AL1 差 1e-3 以上）
  面积系数 0.7854（f_i、F_i 同）；π = math.pi（T_i 打印 3 位小数无法更细分辨，
  但 ω_i、ZS_i、FTH_i 的 6~9 位命中共同确认 π 为双精度）
输出：逐字复刻原著 .OUT 版式（GBK 原件，本内核输出 UTF-8 文本）。
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-23"
TITLE = "小波动过渡过程计算程序"
AUTHOR = "姚廉华"
HEAD_NAME = "D-23"

G = 9.81
K_AREA = 0.7854

_KEYS12 = ["M", "N", "D1", "D2", "D01", "D02", "S1", "S2", "Q0", "H0", "C1", "C2"]
_KEYS11 = ["M", "D1", "D2", "D01", "D02", "S1", "S2", "Q0", "H0", "C1", "C2"]


# ============================================================
# 打印取整（VB6 Format 口径）
# ============================================================

def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    return "%*.*f" % (w, d, q(x, d))


def FE(x, d=8):
    """VB6 科学计数法（8 位小数、两位指数）：5.65484805E-03"""
    return "%0.*E" % (d, float(x))


# ============================================================
# 解析
# ============================================================

def parse(data):
    """解析输入。data：dict | .INT 文件路径。
    G 盘权威 .INT 为 11 个数（无 N）；说明书算例 D-23.IN 为 12 个数（M,N,…）。"""
    if isinstance(data, dict):
        p = dict(data)
        for k in ("M", "N"):
            if k in p:
                p[k] = int(round(p[k]))
        return p
    nums = read_numbers(data)
    n = len(nums)
    if n == 12:
        p = dict(zip(_KEYS12, nums))
    elif n == 11:
        p = dict(zip(_KEYS11, nums))
        p["N"] = 0
    else:
        raise ValueError("D-23 输入数据个数 %d 应为 12（M,N,…）或 11（无 N）" % n)
    p["M"] = int(round(p["M"]))
    p["N"] = int(round(p["N"]))
    return p


# ============================================================
# 计算
# ============================================================

def constants(p):
    D1, D2 = p["D1"], p["D2"]
    D01, D02 = p["D01"], p["D02"]
    L1, L2 = p["S1"], p["S2"]
    Q0, H0 = p["Q0"], p["H0"]
    C1, C2 = p["C1"], p["C2"]
    f1 = K_AREA * D1 * D1
    f2 = K_AREA * D2 * D2
    F1 = K_AREA * D01 * D01
    F2 = K_AREA * D02 * D02
    V10 = Q0 / f1
    V20 = Q0 / f2
    HN = H0 - C1 - C2
    w1 = math.sqrt(G * f1 / (L1 * F1))
    w2 = math.sqrt(G * f2 / (L2 * F2))
    T1 = 2.0 * math.pi / w1
    T2 = 2.0 * math.pi / w2
    TT = T1 / T2
    ZS1 = Q0 / (F1 * w1)
    ZS2 = Q0 / (F2 * w2)
    H10 = HN / ZS1
    H20 = HN / ZS2
    P10 = C1 / ZS1
    P20 = C2 / ZS2
    k1 = 2.0 * G * C1 / (L1 * V10)
    k2 = 2.0 * G * C2 / (L2 * V20)
    n1 = 2.0 * G * C1 * HN * F1 * F1 / (L1 * f1 * f1 * V10 * V10)
    n2 = 2.0 * G * C2 * HN * F2 * F2 / (L2 * f2 * f2 * V20 * V20)
    FTH1 = F1 / n1
    FTH2 = F2 / n2
    cq = Q0 / HN
    c1 = cq / F1
    c2 = cq / F2
    A1 = k1 - c1
    A2 = k2 - c2
    B1 = w1 * w1 - k1 * c1
    B2 = w2 * w2 - k2 * c2
    AA1 = -k1 / C1
    AA2 = -k2 / C2
    BB1 = -k1 * k1 / C1
    BB2 = -k2 * k2 / C2
    # 特征四次方程系数（分块模型；AL1/AL4 精确，AL2/AL3 为 DECL）
    X = c1 * c2
    Y = c1 * c2 * (k1 + k2)
    Z = c1 * c2 * k1 * k2
    AL1 = A1 + A2
    AL2 = B1 + B2 + A1 * A2 - X
    AL3 = A1 * B2 + A2 * B1 - Y
    AL4 = B1 * B2 - Z
    R = -(AL1 * AL2 * AL3 - AL3 * AL3 - AL1 * AL1 * AL4)
    return {"f1": f1, "f2": f2, "F1": F1, "F2": F2, "V10": V10, "V20": V20,
            "HN": HN, "w1": w1, "w2": w2, "T1": T1, "T2": T2, "TT": TT,
            "ZS1": ZS1, "ZS2": ZS2, "H10": H10, "H20": H20, "P10": P10, "P20": P20,
            "k1": k1, "k2": k2, "AN1": n1, "AN2": n2, "FTH1": FTH1, "FTH2": FTH2,
            "A1": A1, "A2": A2, "B1": B1, "B2": B2,
            "AA1": AA1, "AA2": AA2, "BB1": BB1, "BB2": BB2,
            "AL1": AL1, "AL2": AL2, "AL3": AL3, "AL4": AL4, "R": R}


def compute(p):
    k = constants(p)
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR, "输入": p, "常数": k}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, r):
    k = r["常数"]
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****                小波动过渡过程计算程序  D-23                   ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("                   ORIGINAL DATA")
    A.append("                    原 始 数 据")
    A.append(" *******************************************************")
    A.append("      PROJECT NAME  工程名称代号     M=" + FI(p["M"], 5))
    A.append("      引水管道直径                  D1=" + F(p["D1"], 9, 2) + " (m)")
    A.append("      尾水管道直径                  D2=" + F(p["D2"], 9, 2) + " (m)")
    A.append("      上游调压井直径               D01=" + F(p["D01"], 8, 2) + " (m)")
    A.append("      尾水调压井直径               D02=" + F(p["D02"], 8, 2) + " (m)")
    A.append("      上游管道长度                  S1=" + F(p["S1"], 9, 3) + " (m)")
    A.append("      尾水管道长度                  S2=" + F(p["S2"], 9, 3) + " (m)")
    A.append("      管道过流量                    Q0=" + F(p["Q0"], 9, 3) + " (m^3/s)")
    A.append("      上、下游水位差                H0=" + F(p["H0"], 9, 3) + " (m)")
    A.append("      上游管道水头损失值            C1=" + F(p["C1"], 10, 5) + " (m)")
    A.append("      下游管道水头损失值            C2=" + F(p["C2"], 10, 5) + " (m)")
    A.append("")
    A.append("")
    A.append("              RESULT OF CONSTAND DATA")
    A.append("                   中 间 成 果")
    A.append(" =======================================================")
    A.append("    ZS1=" + F(k["ZS1"], 8, 3) + "     ZS2=" + F(k["ZS2"], 8, 3)
             + "     H10=" + F(k["H10"], 8, 3))
    A.append("    H20=" + F(k["H20"], 8, 3) + "     P10=" + F(k["P10"], 8, 3)
             + "     P20=" + F(k["P20"], 8, 3))
    A.append("     A1=" + FE(k["A1"]) + "          A2=" + FE(k["A2"]))
    A.append("     B1=" + FE(k["B1"]) + "          B2=" + FE(k["B2"]))
    A.append("    AA1=" + FE(k["AA1"]) + "         AA2=" + FE(k["AA2"]))
    A.append("    BB1=" + FE(k["BB1"]) + "         BB2=" + FE(k["BB2"]))
    A.append("")
    A.append("")
    A.append("               RESULT OF OSCILATION")
    A.append("                   计 算 成 果")
    A.append(" *******************************************************")
    A.append("      上游调压井水位振荡周期         T1=" + F(k["T1"], 8, 3) + " (秒)")
    A.append("      尾水调压井水位振荡周期         T2=" + F(k["T2"], 8, 3) + " (秒)")
    A.append("      TT  ─ T1/T2                   TT=" + F(k["TT"], 8, 3))
    A.append("      AL1 ─ 相当於公式中的 α1　　　AL1=" + FE(k["AL1"]))
    A.append("      AL2 ─ 相当於公式中的 α2　　　AL2=" + FE(k["AL2"]))
    A.append("      AL3 ─ 相当於公式中的 α3　　　AL3=" + FE(k["AL3"]))
    A.append("      AL4 ─ 相当於公式中的 α4　　　AL4=" + FE(k["AL4"]))
    A.append("      上游调压井的最小稳定断面     FTH1=" + F(k["FTH1"], 8, 3) + " (m^2)")
    A.append("      尾水调压井的最小稳定断面     FTH2=" + F(k["FTH2"], 8, 3) + " (m^2)")
    A.append("      AN1 ─ 相当於公式中的 n1      AN1=" + F(k["AN1"], 8, 3))
    A.append("      AN2 ─ 相当於公式中的 n2      AN2=" + F(k["AN2"], 8, 3))
    A.append("      R ─ 相当於公式中的 Γ　　　　　R=" + FE(k["R"], 9))
    A.append("")
    return "\n".join(A)


def FI(x, w):
    return "%*d" % (w, int(q(x, 0)))


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
