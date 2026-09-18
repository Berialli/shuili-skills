# -*- coding: utf-8 -*-
"""
E-6 用时间平方根法和三点法计算固结系数计算程序 —— 内核
========================================================
复刻《水利水电工程设计计算程序集》E-6（作者：邓铭江，新疆水利厅）。

功能
----
把《规程》求 t90 的**图解过程**改为**解析计算**：以量表读数 d(mm) 为纵坐标、
时间平方根 √t(分) 为横坐标，对开始段直线 ot1 作最小二乘拟合；再过同一截距 ds
作横坐标为原直线 1.15 倍的直线 ot2；ot2 与 d~√t 曲线相邻两点 M1M2 连线的交点
横坐标平方即为固结度 90% 所需时间 t90，进而求固结系数 Cv。并同时给出"三点法"
的 Cv 以便比较。

公式（E-6Intro「二、计算原理」+ 权威 E-6.OUT 逐位复现）
-------------------------------------------------------
    1) 时间平方根法
       d = ds + k·√t            （对序号 N1..N2 的点作最小二乘，本例 1..4）
       ot2: d = ds + (k/1.15)·√t
       取 ot2 自曲线上方穿过曲线的一对相邻点 M1、M2，其连线与 ot2 之交点 x* 满足
           √t90 = x* ,   t90 = x*²  (分)
    2) 固结系数   Cv = 0.848 · h̄² / t90(秒)      (cm²/sec)
       0.848 为 U=90% 对应的时间因数（由 Tv=π/4·U² 或
       1−(8/π²)e^(−π²Tv/4)=0.9 逆解 Tv=0.84802）
    3) 三点法   取 t1、t2=4t1、t3=16t1 三个点 (R1,R2,R3)，
       由三点最小二乘直线 d = d0 + m·√t 得理论零点 d0；按早期线性律
       U ∝ √t ⇒ t ∝ (d−d0)²，故
           t90 = t1 · [(R3−d0)/(R1−d0)]²
       Cv = 0.848 · h̄² / t90(秒)

★ 反演 / 未闭合点（逐条量化；详见 _reports/E6_E10_report.md）
----------------------------------------------------------
  U1. **h̄（试样最大排水距离）口径未闭合**。按《规程》标准口径"某级压力下试样起始
      和终了高度平均值之半"（h̄ = [(H0−d_1)+(H0−d_n)]/4/10）得 h̄=0.917350 cm、
      Cv=0.0023171，与权威 0.0023323 相差 **−0.65%**。由权威 Cv 反解所需
      h̄* = 0.9203440 cm（等价"等效终了读数"1.82824 mm，**不落在任何实测点上**，
      介于 t=42.25min 的 1.820 与 t=60min 的 1.836 之间）→ 口径未定位。
      已排除的候选（量化残差）：
        · 直接取 H0/2（忽略变形）h̄=1.000000 → Cv=0.0027536（+18.07%）
        · 取"起始读数/主固结末读数(60min)" h̄=0.920150 → Cv=0.0023315（−0.034%）
        · 取 3 个三点法读数的平均高度之半 h̄=0.923150 → Cv=0.0023468（+0.62%）
        · t90 误以"分"代"秒" → 0.139938（+5900%）
      常数 0.848 已由 **E-6vb.EXE 常量池 offset 88122 的 double 0.848** 直接确证；
      即未闭合的只是 h̄ 的取值口径，本内核按标准口径实现并如实标注 DECL。
  U2. **三点法 t90 亦按同一 h̄ 口径**：结构式 t90=t1·[(R3−d0)/(R1−d0)]²（d0 取三点
      最小二乘截距 1.384001）给出 t90₃=272.5 s、Cv₃=0.0026187，与权威 0.0026120
      相差 **+0.26%**。权威反解 t90₃=273.2 s（若 h̄ 取标准口径）或 275.0 s
      （若 h̄ 取 h̄*），两种 h̄ 口径不能同时命中 U1、U2，说明 三点法 t90 的结构式
      与原著仍有细微差异（未定位）。
  U3. 权威 .OUT 首行「文件：…」为运行期路径回显，本内核不生成。
  U4. 时间列标题为 t(秒) 但数值实为**分**（三点法把同一时间×60 打印为"秒"；
      例 t=0.25 分 → 15 秒），本内核照排该口径。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | Tv = Cv·t/H²；U = 1 − (8/π²)e^(−π²Tv/4) | `01_水工设计手册精读/卷1_基础理论.md` 4.6「固结（单向固结太沙基）Tv=cvt/H²；U=1−(8/π²)e^(−π²Tv/4)」 | **一致**；U=0.9 逆解 Tv=0.84802 ≈ 0.848 | 采用（0.848 亦由 EXE 常量池确证） |
  | Cv = 0.848·h̄²/t90（h̄ 最大排水距离，t90 以秒计） | 卷1 4.6；`03_水利程序集/程序词汇×教材概念速查表.md` 族谱定位 E-6 | **一致** | 采用（h̄ 口径见 U1） |
  | 时间平方根法：延长开始段直线交纵轴于理论零点 ds，再作横坐标 1.15 倍的直线，交点定 t90 | 卷1 4.6；教材固结试验章 | **一致** | 采用（解析化） |
  | 三点法（理论零点/理论终点） | 卷1 4.6 提及"三点法"但未载公式；教材固结试验章亦未载 | **库中无公式** | 由权威 OUT 反演结构与常数（U2） |

输入数据（E-6Intro「数据语句的顺序」）
---------------------------------------
    第 1 行： N, H0, L2, L1      （读数次数、试样原始高度 mm、该级/前一级仪器变形量 mm）
    其后 N 行： t, d             （经过时间（分）、量表读数 mm）
    另需起始/终止序号 N1,N2（原著交互输入；本例 1,4）。

输出（逐字复刻原著 .OUT 版式）
------------------------------
    试验数据表 → SQR t90、t90、Cv → 三点法计算结果（三点数据 + Cv）
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-6"
TITLE = "用时间平方根法计算固结系数"
AUTHOR = "邓铭江"
HEAD_NAME = "E-6"

# U=90% 的时间因数（Tv）；亦由 E-6vb.EXE 常量池 double 0.848 确证
CV_COEF = 0.848
# time-square-root 法 ot2 直线的横坐标放大倍数（规程 1.15）
OT2_FACTOR = 1.15
# 拟合段默认序号（原著交互输入；权威算例 1,4）
DEFAULT_N1N2 = (1, 4)


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    return "%*.*f" % (w, d, q(x, d))


def qstr(x):
    """VB6 CStr(Double) 口径：最多 15 位有效数字。"""
    v = float(x)
    s = repr(v)
    if "e" in s or "E" in s:
        return s
    neg = s.startswith("-")
    if neg:
        s = s[1:]
    ip, _, fp = s.partition(".")
    digits = (ip + fp).lstrip("0")
    if len(digits) > 15:
        s = "%.15g" % abs(v)
    return ("-" if neg else "") + s


# ============================================================
# 解析
# ============================================================

def parse(data):
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("N1N2", DEFAULT_N1N2)
        return p
    nums = read_numbers(data)
    if len(nums) < 4:
        raise ValueError("E-6 输入数据过少（%d 个数）" % len(nums))
    n = int(round(nums[0]))
    H0, L2, L1 = nums[1], nums[2], nums[3]
    rest = nums[4:]
    rows = []
    for i in range(n):
        rows.append((rest[2 * i], rest[2 * i + 1]))
    return {"N": n, "H0": H0, "L2": L2, "L1": L1, "rows": rows,
            "N1N2": DEFAULT_N1N2}


# ============================================================
# 计算
# ============================================================

def _lin(xs, ys):
    """VB6 口径的原点和式最小二乘（与权威 OUT 的 15 位有效数字逐位一致）。"""
    n = len(xs)
    Sx = sum(xs)
    Sy = sum(ys)
    Sxy = sum(x * y for x, y in zip(xs, ys))
    Sxx = sum(x * x for x in xs)
    Syy = sum(y * y for y in ys)
    k = (n * Sxy - Sx * Sy) / (n * Sxx - Sx * Sx)
    b = (Sy - k * Sx) / n
    den = (n * Sxx - Sx * Sx) * (n * Syy - Sy * Sy)
    r = (n * Sxy - Sx * Sy) / math.sqrt(den) if den > 0 else 0.0
    return k, b, r


def compute(p):
    rows = p["rows"]
    H0 = p["H0"]
    n1, n2 = p["N1N2"]
    xs = [math.sqrt(t) for (t, d) in rows]
    ds_list = [d for (t, d) in rows]

    # ot1 最小二乘（序号 N1..N2，1 基）
    k, ds, _r = _lin(xs[n1 - 1:n2], ds_list[n1 - 1:n2])
    k2 = k / OT2_FACTOR

    def gap(i):
        return ds_list[i] - (ds + k2 * xs[i])

    # ot2 自曲线上方穿过曲线：取 (d − ot2) 由正变负的一对相邻点
    pair = None
    for i in range(len(rows) - 1):
        if gap(i) > 0.0 and gap(i + 1) <= 0.0:
            pair = (i, i + 1)
            break
    if pair is None:                      # 兜底：任意变号处
        for i in range(len(rows) - 1):
            if gap(i) * gap(i + 1) < 0:
                pair = (i, i + 1)
                break
    if pair is None:
        raise ValueError("E-6 ot2 直线未与 d~√t 曲线相交（请检查 N1N2 与数据）")
    i, j = pair
    b_m = (ds_list[j] - ds_list[i]) / (xs[j] - xs[i])
    # 以 M1 为起点沿连线走到与 ot2 的交点（与权威 OUT 逐位一致的形式）
    x90 = xs[i] + (ds + k2 * xs[i] - ds_list[i]) / (b_m - k2)
    t90 = x90 * x90                       # 分

    # 试样最大排水距离 h̄：起始与终了高度平均值之半（标准口径，见 U1）
    hbar = ((H0 - ds_list[0]) + (H0 - ds_list[-1])) / 2.0 / 2.0 / 10.0
    Cv = CV_COEF * hbar * hbar / (t90 * 60.0)

    # ---- 三点法：t1、4t1、16t1 ----
    times = [t for (t, d) in rows]
    trip = None
    for t1 in times:
        if t1 <= 0:
            continue
        cand = []
        for m in (1.0, 4.0, 16.0):
            hit = [idx for idx, tv in enumerate(times) if abs(tv - t1 * m) < 1e-9]
            if not hit:
                break
            cand.append(hit[0])
        if len(cand) == 3:
            trip = cand
            break
    if trip is None:
        trip = [1, 2, 4]
    tt = [times[idx] for idx in trip]
    RR = [ds_list[idx] for idx in trip]
    d0 = _lin([math.sqrt(x) for x in tt], RR)[1]
    t90_3 = tt[0] * ((RR[2] - d0) / (RR[0] - d0)) ** 2
    Cv3 = CV_COEF * hbar * hbar / (t90_3 * 60.0)
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": p, "k": k, "ds": ds, "k2": k2, "x90": x90, "t90": t90,
            "hbar": hbar, "Cv": Cv, "三点序号": trip, "三点t": tt,
            "三点R": RR, "d0": d0, "t90_3": t90_3, "Cv3": Cv3}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, res):
    B = []
    B.append("")
    B.append(" ***********************************************************************")
    B.append(" ****                用时间平方根法计算固结系数 E-6                 ****")
    B.append(" ***********************************************************************")
    B.append("")
    B.append("                  固结试验")
    B.append("")
    B.append("  序     号    各次读数经过的时间      各次量表的读数")
    B.append("      N              t(秒)                  d(毫米)")
    for idx, (t, d) in enumerate(p["rows"]):
        B.append(F(idx + 1, 5, 0) + F(t, 23, 3) + F(d, 22, 3))
    B.append(" SQR t90=" + F(res["x90"], 8, 3))
    B.append(" 固结度达90%%所需时间t90= %s " % qstr(res["t90"]))
    B.append(" 固结系数Cv=" + F(res["Cv"], 11, 7))
    B.append("")
    B.append("")
    B.append("                  三点法计算结果")
    B.append("                  试验数据")
    B.append("         量表读数(毫米)       相应的时间(秒)")
    for tv, RR in zip(res["三点t"], res["三点R"]):
        B.append(F(RR, 17, 4) + "        " + "%d" % int(round(tv * 60.0)) + " ")
    B.append(" 固结系数Cv=" + F(res["Cv3"], 10, 7))
    return "\n".join(B) + "\n"


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
