# -*- coding: utf-8 -*-
"""
D-11 Ｄ－１１ 侧槽式变量流水面曲线计算程序 —— 内核
==================================================
复刻《水利水电工程设计计算程序集》D-11 程序（作者：谭冬初，新疆农业大学水利系）。

功能
----
侧槽式变量水流广泛用于**侧槽式溢洪道**与**底栏栅式渠首引水廊道**。
本程序自侧槽末端（控制断面）向上游逐段推算水面曲线，输出各断面
距首端的距离 X、流量 Q、水深 H、流速 V，供设计方案比较与校核验算。

输入数据（.INT，逗号分隔，固定顺序 9 个数）
------------------------------------------
    Z, M, I0, N, LM, BU, B0, Q, H0
    Z   侧槽计算分段数（起点在侧槽首端，Z 越大精度越高）
    M   侧槽内边坡系数（梯形断面，m）
    I0  侧槽内底坡
    N   侧槽内衬砌糙率
    LM  侧槽长度 (m)
    BU  侧槽首端槽底宽度 (m)
    B0  侧槽末端槽底宽度 (m)
    Q   侧槽末端处流量（= 侧槽总进流量）(m³/s)
    H0  控制断面起端水深（= 调整段起端水深），H0 = (1.05~1.50)·HK
    例：`6,.5,.024,.018,60,5,15,240,3.8`

断面与水力要素
--------------
    底宽沿程线性变化：      B(x) = (B0-BU)·x/LM + BU
    梯形过水面积：          A(x,h) = B(x)·h + M·h²
    水面宽：                Bw(x,h) = B(x) + 2·M·h
    流速：                  V = Q / A

核心公式（逐字还原自 D-11Intro.rtf 的 OLE 公式对象，见下）
----------------------------------------------------------
式(1) 一般动力微分方程（侧槽变量流，由动量方程推得）：

    dy/dx = [ i₀ − i_f − 2α·Q·(dQ/dx)/(g·A²) ] / [ 1 − α·Q²·B_w/(g·A³) ]

式(2) 取 α = 1、i_f = 0 后的差分方程（程序实际使用；y 为水面高程）：

    Δy = Q₁(V₁+V₂)/(g(Q₁+Q₂)) · [ (V₂−V₁) + V₂·(Q₂−Q₁)/Q₁ ]

    ① 为上游断面（较小流量），② 为下游断面（较已知流量）；
    Q₁ → 0（侧槽首端）时该式退化为极限 Δy = V₂²/g。

式(3) 几何关系（① 上游、② 下游）：

    Δy = h₁ − h₂ + i₀·Δx     ⟺     h₁ = h₂ + Δy − i₀·Δx

式(4) 侧槽底坡限值（调整段起端水深取 (1.05~1.50)h_K，末端为临界水深）：

    i_K = (2/L)·(A_K/B_K)
    A_K、B_K 为与侧槽**末端**临界水深 h_K 对应的过水面积与水面宽；
    梯形断面 A_K = b·h_K + m·h_K²，B_K = b + 2m·h_K（b = B0）。

    程序要求 I0 ≤ I_K，否则报「槽内底坡I0>临界底坡Ik,改变数据重算!」。

临界水深（程序实测口径）
------------------------
h_K 解 α_K·Q²·B_w/(g·A³) = 1（断面取侧槽末端，即 B = B0），
其中 α_K = 1.1（与说明书「取 α = 1」的式(2) 所用 α 不同；见"未闭合点 2"）。
程序以等步长试算（step = 1e-3）取第一个满足该式 ≤ 1 的水深。

计算流程（推算方向：末端 → 首端）
----------------------------------
  1. 由 H0、Q 算出 h_K、I_K，校验 I0 ≤ I_K；
  2. NO = Z（x = LM）处 h = H0；
  3. 对 k = Z, Z−1, …, 1：已知 x₂ = k·Δx 处 h₂、Q₂ = Q·k/Z，
     推求 x₁ = (k−1)·Δx 处 h₁、Q₁ = Q·(k−1)/Z（试算）：

         h₁ := h₂                       （初值：取下游断面水深）
         Δy := 式(2)(h₁ = h₂)
         重复：h₁ ← h₂ + Δy − i₀·Δx ;  Δy' ← 式(2)(h₁)
               若 |Δy′ − Δy| ≤ ε 则结束，否则 Δy ← Δy′
         返回该 h₁

     判据 |Δy′ − Δy| 即说明书「由式(2)求得的 Δy 须逼近由式(3)求得的 Δy」
     （由式(3)，Δy = h₁ − h₂ + i₀Δx 恰为上一轮的 Δy）。
     ε = 0.005（说明书未给出数值，由权威 OUT 反演，见"未闭合点 1"）。
  4. 输出 NO = Z…0 的距离、流量、水深、流速。

输出（逐字复刻原著 .OUT 版式）
------------------------------
  一.原始数据：Z / M / I0 / N / LM / BU / B0 / Q / H0
  二.计算结果：临界水深 Hk、临界底坡 Ik；表 NO / 距离 / 流量 / 水深 / 流速
  （Hk 打印 3 位小数、Ik 打印 4 位小数；表内数值打印 2 位小数）

未闭合点（如实标注）
--------------------
  1. **迭代容差 ε 未在说明书中给出数值**（原文仅「直至满足某精度要求为止」）。
     由权威 D-11.OUT 的 7 行 × (h, v) 共 14 个打印值反演：ε ∈ [0.0034, 0.0076]
     时可**逐位复现**全部 14 个值；ε ≤ 0.0033 时 NO=2 行打印 4.20（原 4.21），
     ε ≥ 0.0077 时 NO=5/3 行偏离。内核取 ε = 0.005（区间内的整齐值），
     例 1 全部 14 个打印值逐位复现。
  2. **临界水深系数 α_K = 1.1（而非说明书式(2)的 α = 1）**。若按 α = 1 求梯形
     断面临界水深，例 1 得 h_K = 2.870，与权威 OUT 的 2.960 相差 0.09 m；
     取 α_K = 1.1 时 2.9594。**但 2.9594 打印为 2.959 而权威 OUT 为 2.960**，
     故程序对 h_K 系等步长试算（本内核取步长 1e-3，首个使 αQ²Bw/(gA³) ≤ 1
     的水深 = 2.960）。反证：α_K = 1.1 + 严格二分 → 2.959（✗）；
     α_K = 1.0 → 2.870（✗）；α_K = 1.1 + 步长 ≤ 7e-4 的扫描 → 2.959（✗）。
     由此推定源码以 ~1e-3 网格「列表试算」h_K，本内核照此复刻。
  3. 权威 OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-11.out」为原著运行期
     路径回显（机器相关），本内核不生成该行（与 D-1/D-3/D-4/D-7 内核同一口径）。
  4. 说明书式(1) 与式(2) 并非同一方程的严格离散：把式(1)（α=1, i_f=0）直接
     差分离散后，例 1 首段 |dy/dx| = 0.031~0.040 而实测水面坡降为 0.038，
     且逐段累计会明显偏离；按式(2)（本内核口径）逐段吻合。故程序实际使用的是
     式(2)，式(1) 仅为方法概述。本内核按式(2) 实现。
  5. 说明书式(2) 分母含 Q₁，侧槽首端 Q₁ = 0 时原式为 0/0；本内核取极限
     Δy = V₂²/g（等价于将式(2) 展开为
     Δy = Q₁(V₁+V₂)(V₂−V₁)/(g(Q₁+Q₂)) + V₂(V₁+V₂)(Q₂−Q₁)/(g(Q₁+Q₂)) 后令 Q₁→0）。
     权威 OUT 的 NO=0 行 h = 4.10 与该极限一致（V₂²/g = 0.1190 ⇒ h = 4.098）。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-11"
TITLE = "侧槽式变量流水面曲线计算书"
AUTHOR = "谭冬初（新疆农业大学水利系）"
HEAD_NAME = "D-11"

G = 9.81                # 重力加速度 (m/s²)
ALPHA_CRIT = 1.1        # 临界水深式 α_K·Q²Bw/(gA³)=1 的 α_K（见未闭合点 2）
EPS_DY = 0.005          # 试算容差 |Δy(式2) − Δy(式3)|（见未闭合点 1）
HK_STEP = 0.001         # 临界水深试算步长（见未闭合点 2）
MAX_ITER = 2000         # 试算迭代上限（保护）


# ============================================================
# 断面几何
# ============================================================

def bottom_width(x, LM, BU, B0):
    """侧槽计算断面底宽（单边/两边扩散型线性变化）。"""
    if LM == 0:
        return BU
    return (B0 - BU) * x / LM + BU


def area(x, h, LM, BU, B0, M):
    """梯形断面过水面积 A = B(x)·h + M·h²。"""
    return bottom_width(x, LM, BU, B0) * h + M * h * h


def top_width(x, h, LM, BU, B0, M):
    """水面宽 B_w = B(x) + 2·M·h。"""
    return bottom_width(x, LM, BU, B0) + 2.0 * M * h


def velocity(x, h, q, LM, BU, B0, M):
    """断面平均流速 V = Q/A；Q = 0 时取 0（首端无水流入）。"""
    if q <= 0:
        return 0.0
    a = area(x, h, LM, BU, B0, M)
    return q / a if a > 0 else 0.0


# ============================================================
# 临界水深 / 临界底坡
# ============================================================

def critical_flow_ratio(h, q, LM, BU, B0, M):
    """α_K·Q²·B_w/(g·A³)；= 1 即临界流（断面取侧槽末端）。"""
    a = area(LM, h, LM, BU, B0, M)
    if a <= 0:
        return float("inf")
    bw = top_width(LM, h, LM, BU, B0, M)
    return ALPHA_CRIT * q * q * bw / (G * a ** 3)


def critical_depth(q, LM, BU, B0, M):
    """
    侧槽末端临界水深：1e-3 等步长网格试算，取**首个**使 αQ²Bw/(gA³) ≤ 1 的水深。
    （见"未闭合点 2"：与原著 1e-3 网格列表试算口径一致）

    说明：为免大流量下逐点扫描过慢，先用步长倍增法框定上界 H
    （满足 ratio(H) ≤ 1 且 ratio(H/2) > 1），再从 H/2 起按 1e-3 网格逐点试算；
    由于网格对齐于 0，所得"首个满足点"与全程逐点扫描完全一致。
    """
    if critical_flow_ratio(HK_STEP, q, LM, BU, B0, M) <= 1.0:
        return HK_STEP
    hi = HK_STEP * 2.0
    while critical_flow_ratio(hi, q, LM, BU, B0, M) > 1.0:
        hi *= 2.0
        if hi > 1.0e7:
            raise ValueError("临界水深试算不收敛（Q=%g）" % q)
    k = int(math.ceil((hi / 2.0) / HK_STEP))
    h = k * HK_STEP
    while h <= hi:
        if critical_flow_ratio(h, q, LM, BU, B0, M) <= 1.0:
            return h
        h += HK_STEP
    return hi


def critical_slope(hk, LM, BU, B0, M):
    """式(4)：i_K = (2/L)·(A_K/B_K)。"""
    if LM <= 0:
        return 0.0
    ak = area(LM, hk, LM, BU, B0, M)
    bk = top_width(LM, hk, LM, BU, B0, M)
    return (2.0 / LM) * ak / bk


# ============================================================
# 式(2)：断面间水面降落差
# ============================================================

def delta_y(x1, h1, q1, x2, h2, q2, LM, BU, B0, M):
    """
    式(2)：Δy = Q₁(V₁+V₂)/(g(Q₁+Q₂))·[(V₂−V₁) + V₂(Q₂−Q₁)/Q₁]
    ① 上游（x1 较小流量 q1）、② 下游（x2 较大流量 q2）。
    q1 = 0（侧槽首端）时退化为极限 Δy = V₂²/g。
    """
    v1 = velocity(x1, h1, q1, LM, BU, B0, M)
    v2 = velocity(x2, h2, q2, LM, BU, B0, M)
    if q1 <= 0:
        return v2 * v2 / G
    first = q1 * (v1 + v2) / (G * (q1 + q2))
    return first * ((v2 - v1) + v2 * (q2 - q1) / q1)


def solve_upstream_depth(h2, q1, q2, x1, x2, dx, i0, eps,
                         LM, BU, B0, M):
    """
    由下游断面 (x2, h2, q2) 试算上游断面水深 h1。

    按说明书：「假设 h1 由式(3)求得 Δy，再由式(2)求得 Δy，须逼近至满足精度」。
    实现为不动点迭代（初值 h1 := h2）：

        h₁ ← h₂ + Δy(式2) − i₀·Δx            （此即式(3)）
        直至 |Δy(式2 于新 h₁) − Δy(式3 于新 h₁)| ≤ ε

    返回该 h₁。
    """
    dy = delta_y(x1, h2, q1, x2, h2, q2, LM, BU, B0, M)
    h1 = h2
    for _ in range(MAX_ITER):
        h1 = h2 + dy - i0 * dx
        dy_new = delta_y(x1, h1, q1, x2, h2, q2, LM, BU, B0, M)
        if abs(dy_new - dy) <= eps:
            return h1
        dy = dy_new
    return h1


# ============================================================
# 解析
# ============================================================

_KEYS = ["Z", "M", "I0", "N", "LM", "BU", "B0", "Q", "H0"]


def parse(data):
    """解析输入。data：dict（直接返回）| .INT 文件路径。"""
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 9:
        raise ValueError("D-11 输入数据不足（需 9 个数：Z,M,I0,N,LM,BU,B0,Q,H0，"
                         "实得 %d）" % len(nums))
    p = dict(zip(_KEYS, nums[:9]))
    p["Z"] = int(round(p["Z"]))
    if p["Z"] <= 0:
        raise ValueError("侧槽分段数 Z 必须为正整数（实得 %s）" % p["Z"])
    if p["LM"] <= 0:
        raise ValueError("侧槽长度 LM 必须为正（实得 %s）" % p["LM"])
    return p


# ============================================================
# 计算
# ============================================================

def compute(params):
    """执行 D-11 计算：临界水深/底坡 + 末端→首端逐段推算水面曲线。"""
    Z = int(params["Z"])
    M = params["M"]
    I0 = params["I0"]
    N = params["N"]
    LM = params["LM"]
    BU = params["BU"]
    B0 = params["B0"]
    Q = params["Q"]
    H0 = params["H0"]

    hk = critical_depth(Q, LM, BU, B0, M)
    ik = critical_slope(hk, LM, BU, B0, M)
    warn = None
    if I0 > ik:
        warn = "槽内底坡I0>临界底坡Ik,改变数据重算!"

    dx = LM / Z
    rows = []
    h = H0
    for k in range(Z, -1, -1):
        x = k * dx
        q = Q * k / Z
        v = velocity(x, h, q, LM, BU, B0, M)
        rows.append({"no": k, "x": x, "q": q, "h": h, "v": v})
        if k > 0:
            x1 = (k - 1) * dx
            q1 = Q * (k - 1) / Z
            h = solve_upstream_depth(h, q1, q, x1, x, dx, I0, EPS_DY,
                                     LM, BU, B0, M)

    return {
        "程序": PROGRAM_ID,
        "标题": TITLE,
        "作者": AUTHOR,
        "输入": {"Z": Z, "M": M, "I0": I0, "N": N, "LM": LM,
                 "BU": BU, "B0": B0, "Q": Q, "H0": H0},
        "临界水深Hk": hk,
        "临界底坡Ik": ik,
        "警告": warn,
        "行": rows,
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(params, result):
    """生成原著风格文本计算书（.OUT）。"""
    i = result["输入"]
    Z = i["Z"]
    rows = result["行"]

    L = []
    L.append("")
    L.append(" ***********************************************************************")
    L.append(" ******                侧槽式变量流水面曲线计算书 D-11            ******")
    L.append(" ***********************************************************************")
    L.append("")
    L.append("                    一.原始数据")
    L.append("")
    L.append(" 侧槽分段数 Z=%3d       侧槽内边坡系数 M=%6.3f   侧槽内底坡 I0=%7.4f"
             % (Z, i["M"], i["I0"]))
    L.append(" 侧槽内糙率 N=%7.4f   侧槽长度 LM=%7.2f       侧槽首端底宽 BU=%6.2f"
             % (i["N"], i["LM"], i["BU"]))
    L.append(" 侧槽末端底宽 B0=%6.2f 侧槽末端处流量 Q=%7.2f  起端水深 H0=%6.2f"
             % (i["B0"], i["Q"], i["H0"]))
    L.append("")
    L.append("")
    L.append("                    二.计算结果")
    L.append("")
    L.append(" 临界水深 Hk=%6.3f       临界底坡 Ik=%7.4f"
             % (result["临界水深Hk"], result["临界底坡Ik"]))
    L.append("")
    L.append("  NO" + "".join("           " + s for s in ("距离", "流量", "水深", "流速")))
    for r in rows:
        L.append("%4d%15.2f%15.2f%15.2f%15.2f"
                 % (r["no"], r["x"], r["q"], r["h"], r["v"]))
    if result.get("警告"):
        L.append("")
        L.append(" *** " + result["警告"])
    L.append("")
    return "\n".join(L)


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
