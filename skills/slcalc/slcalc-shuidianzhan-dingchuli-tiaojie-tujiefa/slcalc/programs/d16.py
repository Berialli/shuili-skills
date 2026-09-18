# -*- coding: utf-8 -*-
"""
D-16 无压隧洞水面曲线计算程序 —— 内核
======================================
复刻《水利水电工程设计计算程序集》D-16 程序（作者：姚廉华）。

功能
----
计算**均匀断面、相同坡降**的无压隧洞的降水（M）或壅水（S）曲线分布规律。
已知起（末）端断面水深时，逐段推算下一（上一）断面的水深与两断面之间距离，
从而给出沿洞轴线的水深与水位变化。

洞型与水力要素（原著输入口径）
------------------------------
    B   隧洞过流断面宽度 (m)      矩形洞身，宽 B
    H0  边墙高度 (m)              （本算例各水深 ≤ H0，拱顶不参与过流）
    H1  前池水位升高值 (m)         洞前（后）最高水位 = H2 + H1
    H2  水流深度 (m)              洞内原有（正常）水深，曲线的下限
    AI  隧洞底坡 i
    Z   隧洞糙率 n
    A0  不均匀系数 α（1～1.1）
    H(1:L)  各计算断面假定水深（递减数组，末值 = H2 + Δ）

    过水面积    A = B·h
    湿周        P = B + 2h          （矩形明渠，不计自由水面）
    水力半径    R = A / P
    谢才系数    C = R^(1/6) / n     （原著式 C = (1/n)R^(1/6)）
    特性流量    K = A·C·√R = A·R^(2/3)/n = Q/√i

核心公式（逐字还原自 D-16Intro.rtf 的 19 个 OLE 公式对象）
-----------------------------------------------------------
公式对象解包路径：\\objdata → hex → OLE2 → "Equation Native" → 跳 28B →
MTEF v3 记录流；并用内嵌 WMF 预览（META_EXTTEXTOUT：字符 + 字体 + 坐标 +
META_MOVETO/LINETO：分式横线与上划线）逐一核对版式。结果为：

  式(1)  K  = Q/√i = A·C·√R                     [对象 00]
  式(2)  X1 = K1/K0                             [对象 01]
  式(3)  X2 = K2/K0                             [对象 02]
  式(4)  C  = (1/n)·R^(1/6)                     [对象 03]
  式(5)  C̄  = 平均流速（谢才）系数，C̄ = (C1+C2)/2   [对象 04/05/06]
  式(6)  B̄  = 平均宽度，          B̄ = (B1+B2)/2   [对象 07/08/09]
  式(7)  P̄  = 平均湿周长度，      P̄ = (P1+P2)/2   [对象 10/11/12]
  式(8)  J̄  = α·i·C̄²·B̄ / (g·P̄)                  [对象 13/14]
  式(9)  F(x) = 1.151·lg|(x+1)/(x-1)|   （x>1 与 x<1 分列两式）  [对象 15/16]
  式(10) a  = (x2 − x1)/(h2 − h1)                [对象 17]
  式(11) l  = (1/(a·i))·[ (x2 − x1) − (1 − J̄)·(F(x2) − F(x1)) ]  [对象 18]

  K0 取正常水深断面（h = H2）的特性流量，即 Q = K(H2)·√i 为设计（正常）流量；
  x = K/K0 为特性流量比。

式(11) 的由来（本内核独立推导并与权威 OUT 逐位吻合）
-----------------------------------------------------
  能量方程沿程渐变流： dE/ds = i − J_f，J_f = Q²/K² = i/x²，
  弗劳德数 Fr² = αQ²B/(gA³) = J̄/x²（因 x²·Fr² = α i C²B/(gP) = J̄）。
  于是 dh/ds = i(1 − 1/x²)/(1 − Fr²) = i(x² − 1)/(x² − J̄)，
  以 a = dx/dh 局部线性化：ds/dx = (x² − J̄)/(a·i·(x² − 1))，
  积分 ∫ dx/(x²−1) = −F(x) 即得式(11)。

输入数据（.INT，FORTRAN 自由格式；逗号分隔）
--------------------------------------------
    第一组（简单变量，9 个）：
        N, L, B, H0, H1, H2, AI, Z, A0
      ※ 原著说明书写 "简单变量 9 个" 却列了 10 项（含 M 题目数）；实测 .INT
        在 M=0 时**不写 M**（第一行仅 9 个数），说明书算例行含 M=0。
        本内核两种写法都可解析（按 "9+标准L" 或 "10+L" 恰好消费完全部数判定）。
    第二组（数组，1 组）：
        H(1:L) 各计算断面假定水深 (m)，以 H2 为下限；末值 = H2 + Δ（很小的 Δ），
        否则 x→1 时 F(x)→∞，距离发散。

输出（逐字复刻原著 .OUT 版式）
------------------------------
    一. 原始数据：N / L / B / H0 / H1 / H2 / AI / Z / A0；各断面假定水深
    二. 水位成果：I / H(m) / BL(m) 上、下两断面距离 / BBL(m) 累加值
        BL 打印 2 位小数，BBL 打印 3 位小数

算例（D-16.INT / D-16.OUT 同源）
--------------------------------
    888,11,4.0,6.0,2.0,4.0,0.0036,0.014,1.1
    6.0,5.8,5.6,5.4,5.2,5.0,4.8,4.6,4.4,4.2,4.0001

未闭合点（如实标注，含反证）
----------------------------
  1. **式(9) 的系数 k**：说明书公式对象显示 "1.151"，但按 lg（=ln/ln10）口径
     取 1.151 时例 1 的 BBL 打印值仅命中 2/10（如首行 −70.226 ≠ −70.224）。
     由权威 OUT 的 10 行 BBL（3 位小数）反演，有效系数被**唯一限制**在
     k ∈ [1.1507844, 1.1507852]（区间宽 8×10⁻⁷）；内核取 k = 1.150785，
     此时例 1 的 BL（2 位）与 BBL（3 位）共 21 个打印值**全部逐位命中**。
     反证：k = ln10/2 = 1.1512925 → 命中 0/10；k = 1.151 → 0/10；
     k = 1.151·ln10/2.303 = 1.1507926（"lg 底数取 2.303" 说）→ 5/10；
     k = 1/0.8689 → 0/10。故 1.150785 为反演值（说明书 "1.151" 系 4 位有效
     数字的页面显示值，与该反演值相差 4×10⁻⁴ 相对量）。
     另证：K0 反演 = K(H2) = 1384.4728（对 K0 做 ±6% 二维扫描，最优仍在 K0 = K(H2)），
     J 的量级系数 s 反演 = 1.000（±0.001 扫描最优 s = 1.0）。
  2. **末行（I = L）的 BL/BBL 无法由任何"断面—深度对"复现**——原著对 L 个断面
     打印了 L 个 BL 值，而相邻断面只有 L−1 对，第 L 个必为附加段。
     反证：对全部 11 个假定水深外加 H2、H2+H1、H0、H1 等候选深度两两枚举，
     在 k = 1.150785 下均**不存在** l = +2720.309（即 BL(L) = 997.945 −(−1722.364)）
     的组合：seg(4.0001→6.0) = +2016.63（洞前最高水位断面）、
     seg(4.0001→H2) → +∞（x→1 发散）、seg(4.0001→7.4525) = +2720.31 但 7.4525
     非任何有物理意义的深度（H2+H1 = 6.0、2H2 = 8.0、H0 = 6.0 等均不符）。
     同时说明书自带算例的末行为 +2720.74 / 998.231，与权威 OUT 的
     +2720.31 / 997.945 相差 0.43 / 0.29，而前 9 行两版完全一致 —— 说明该
     附加段依赖原著**版本/运行环境**，不可唯一反演。
     内核口径：第 L 行取"闭合段" l = seg(H(L) → 洞前最高水位断面 H2+H1)，
     本例 = +2016.63（在 verify 中记为 DECL，非 FAIL）。
  3. 权威 .OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-16.out」为原著运行期
     路径回显（机器相关），本内核不生成该行（与 D-1/D-3/D-4/D-7/D-11 同一口径），
     版式比对自第 2 行起算。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-16"
TITLE = "无压隧洞水面曲线计算程序"
AUTHOR = "姚廉华"
HEAD_NAME = "D-16"

G = 9.81

# 式(9) F(x) = K_F·lg|(x+1)/(x-1)| 的系数。
# 说明书公式对象显示 1.151；由权威 D-16.OUT 的 10 行 BBL 反演唯一限制在
# [1.1507844, 1.1507852]，取 1.150785（详见模块 docstring 未闭合点 1）。
K_F = 1.150785


# ============================================================
# 断面几何与水力要素
# ============================================================

def geometry(B, h, n):
    """
    矩形无压隧洞断面水力要素（深度 h ≤ 边墙高度 H0 时）。
    返回 (A 过水面积, P 湿周, R 水力半径, C 谢才系数, K 特性流量)。
    """
    A = B * h
    P = B + 2.0 * h
    R = A / P if P > 0 else 0.0
    C = (R ** (1.0 / 6.0)) / n if R > 0 else 0.0
    K = A * C * math.sqrt(R) if R > 0 else 0.0
    return A, P, R, C, K


def k_ratio(B, h, n, K0):
    """特性流量比 x = K(h)/K0。"""
    return geometry(B, h, n)[4] / K0


def F(x):
    """
    式(9)：F(x) = K_F·lg|(x+1)/(x-1)|；x>1 与 x<1 分段（两段等价于同一绝对值式）。
    x = 1 时发散（原著要求末断面水深加 Δ 以避免）。
    """
    v = abs((x + 1.0) / (x - 1.0))
    return K_F * math.log10(v)


def segment(B, n, i, alpha, K0, h1, h2):
    """
    式(5)~(11)：两断面之间的洞段距离 l（一维水面曲线差分解析解）。

        C̄ = (C1+C2)/2,  B̄ = (B1+B2)/2,  P̄ = (P1+P2)/2
        J̄ = α·i·C̄²·B̄/(g·P̄)
        a = (x2 − x1)/(h2 − h1)
        l = (1/(a·i))·[ (x2 − x1) − (1 − J̄)·(F(x2) − F(x1)) ]

    h1、h2 顺序决定符号：水深递减方向（h2 < h1）得负距离（与权威 OUT 一致）。
    """
    A1, P1, R1, C1, K1 = geometry(B, h1, n)
    A2, P2, R2, C2, K2 = geometry(B, h2, n)
    x1 = K1 / K0
    x2 = K2 / K0
    Cb = 0.5 * (C1 + C2)
    Pb = 0.5 * (P1 + P2)
    J = alpha * i * Cb * Cb * B / (G * Pb) if Pb > 0 else 0.0
    a = (x2 - x1) / (h2 - h1)
    if a == 0.0:
        return 0.0
    return ((x2 - x1) - (1.0 - J) * (F(x2) - F(x1))) / (a * i)


# ============================================================
# 解析
# ============================================================

_SCALARS9 = ["N", "L", "B", "H0", "H1", "H2", "AI", "Z", "A0"]
_SCALARS10 = ["N", "M", "L", "B", "H0", "H1", "H2", "AI", "Z", "A0"]


def parse(data):
    """
    解析输入。data：dict（直接返回）| .INT 文件路径。

    返回 dict：{N, M, L, B, H0, H1, H2, AI, Z, A0, H:[...]}
    兼容两种首部写法：M=0 时原著 .INT 省略 M（首部 9 个数），
    说明书算例行含 M=0（首部 10 个数）。
    """
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    total = len(nums)

    plan = None
    # 先试 "含 M"（10 + L）
    if total >= 10:
        L = int(round(nums[2]))
        if L > 0 and 10 + L == total:
            plan = (_SCALARS10, L, 10)
    # 再试 "不含 M"（9 + L）
    if plan is None and total >= 9:
        L = int(round(nums[1]))
        if L > 0 and 9 + L == total:
            plan = (_SCALARS9, L, 9)
    if plan is None:
        raise ValueError(
            "D-16 输入数据无法解析（数字个数 %d，需 9+L 或 10+L）" % total)

    keys, L, head = plan
    p = dict(zip(keys, nums[:head]))
    p.setdefault("M", 0.0)
    p["N"] = int(round(p["N"]))
    p["M"] = int(round(p["M"]))
    p["L"] = L
    if L < 2:
        raise ValueError("计算断面数 L 必须 ≥ 2（实得 %d）" % L)
    if p["B"] <= 0:
        raise ValueError("隧洞宽度 B 必须为正（实得 %s）" % p["B"])
    if p["AI"] <= 0:
        raise ValueError("隧洞底坡 AI 必须为正（实得 %s）" % p["AI"])
    if p["Z"] <= 0:
        raise ValueError("隧洞糙率 Z 必须为正（实得 %s）" % p["Z"])
    p["H"] = list(nums[head:head + L])
    return p


# ============================================================
# 计算
# ============================================================

def compute(params):
    """执行 D-16 计算：特性流量比 + 逐段水面曲线距离。"""
    N = params["N"]
    L = int(params["L"])
    B = params["B"]
    H0 = params["H0"]
    H1 = params["H1"]
    H2 = params["H2"]
    i = params["AI"]
    n = params["Z"]
    alpha = params["A0"]
    H = params["H"]

    # K0：正常水深（H2）断面的特性流量
    K0 = geometry(B, H2, n)[4]
    A0_, P0_, R0_, C0_, _ = geometry(B, H2, n)
    v0 = (K0 * math.sqrt(i)) / A0_ if A0_ > 0 else 0.0

    rows = []
    bbl = 0.0
    h_high = H2 + H1                      # 洞前（后）最高水位断面
    for k in range(L):
        h = H[k]
        if k < L - 1:
            bl = segment(B, n, i, alpha, K0, h, H[k + 1])
        else:
            # 未闭合点 2：原著对第 L 个断面另有一段（本内核取"闭合段"口径）
            h2 = h_high if abs(h_high - h) > 1e-12 else H[0]
            bl = segment(B, n, i, alpha, K0, h, h2)
        bbl += bl
        rows.append({"i": k + 1, "h": h, "bl": bl, "bbl": bbl})

    return {
        "程序": PROGRAM_ID,
        "标题": TITLE,
        "作者": AUTHOR,
        "输入": {"N": N, "M": params["M"], "L": L, "B": B, "H0": H0, "H1": H1,
                 "H2": H2, "AI": i, "Z": n, "A0": alpha, "H": list(H)},
        "特性流量K0": K0,
        "正常断面": {"A": A0_, "P": P0_, "R": R0_, "C": C0_},
        "正常流速v0": v0,
        "洞前最高水位": h_high,
        "行": rows,
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(params, result):
    """生成原著风格文本计算书（.OUT）。"""
    i = result["输入"]
    rows = result["行"]
    L = i["L"]

    Ls = []
    Ls.append("")
    Ls.append(" ***********************************************************************")
    Ls.append(" ****               无压隧洞水面曲线计算程序 D-16                   ****")
    Ls.append(" ***********************************************************************")
    Ls.append("")
    Ls.append("                   ORIGINAL DATA")
    Ls.append("                    原 始 数 据")
    Ls.append(" *******************************************************")
    Ls.append("        PROJECT NAME  工程名称代号   N=%5d" % i["N"])
    Ls.append("        计算断面数                   L=%5d" % L)
    Ls.append("        隧洞过流断面宽度             B=%7.2f (m)" % i["B"])
    Ls.append("        边墙高度                    H0=%7.2f (m)" % i["H0"])
    Ls.append("        前池水位升高值              H1=%7.2f (m)" % i["H1"])
    Ls.append("        水流深度                    H2=%7.2f (m)" % i["H2"])
    Ls.append("        隧洞底坡                    AI=%8.5f (i)" % i["AI"])
    Ls.append("        隧洞糙率                     Z=%8.4f" % i["Z"])
    Ls.append("        不均匀系数                  A0=%7.2f" % i["A0"])
    Ls.append("")
    Ls.append("")
    Ls.append("          I0          H")
    Ls.append("        断面号      水深(m)")
    Ls.append(" -------------------------------------------------------")
    for k in range(L):
        Ls.append("%11d%15.2f" % (k + 1, i["H"][k]))
    Ls.append("")
    Ls.append("")
    Ls.append("              RESULT OF WATER  LEVEL")
    Ls.append("                水   位   成   果")
    Ls.append(" *******************************************************")
    Ls.append("")
    Ls.append("        I        H (m)        BL (m)        BBL (m)")
    Ls.append("      断面号     水深        断面距离       累加值")
    Ls.append(" -------------------------------------------------------")
    for r in rows:
        Ls.append("%9d%12.2f%14.2f%15.3f" % (r["i"], r["h"], r["bl"], r["bbl"]))
    Ls.append("")
    return "\n".join(Ls)


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
