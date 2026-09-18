# -*- coding: utf-8 -*-
"""
D-14B 比降法计算天然河道水面曲线程序 —— 内核
=============================================
复刻《水利水电工程设计计算程序集》**D-14B**（作者：张校正，新疆水利厅）。
原著为 VB6 程序（运行程序/EXE/D-14Bvb.EXE，2010 版），本模块为 Python 复刻。

一、程序功能
------------
已知天然河道各横断面的地形点资料、断面处各点的不同糙率、以及**各段的纵坡 i**，
用均匀流（正常水深）公式求出**每一个断面各自的正常水深（水位）**，
打印桩号、纵坡、水位、过水断面积、断面平均流速。
各断面**相互独立**求解（作者论文原话：「D-14B程序各个断面独立计算，
没有误差积累的问题」），这正是"比降法"与 D-14A"推求法"的根本差别。

二、计算公式（说明书 D-14BIntro.rtf 公式对象解包所得，与教材逐条对照）
--------------------------------------------------------------------
[00] Q = ω·C·√(R·i)                          …… 谢才均匀流公式
[01] R = ω / χ                               …… 水力半径（教材 §1.4）
[02] C = (1/n)·R^(1/6)                       …… 曼宁谢才系数（教材 §4.5）
[03] Q = √i·ω^(5/3) / (n·χ^(2/3))            …… [00]+[01]+[02] 代入结果
[04] F(h) = 1 − Q·n·χ(h)^(2/3) / (√i·ω(h)^(5/3)) = 0
                                             …… 求 h0 的方程
[05] h = (h1+h2)/2                           …… 二分法
[06] F(h2)=J、[07] F(h)=K                    …… J·K>0 → h2=h；J·K<0 → h1=h2,h2=h
[08] S_i = √((X_{i+1}−X_i)²+(Y_{i+1}−Y_i)²)                       …… 全淹段湿周
[09] W_i = ½·(X_{i+1}−X_i)·(Z−Y_i)²/(Y_{i+1}−Y_i)                 …… 水边相交段面积
[10] S_i = (Z−Y_i)·√(ΔX²+ΔY²)/(Y_{i+1}−Y_i)                      …… 水边相交段湿周

教材/知识库出处：
  · 知识库 `水力学核心公式速查.md` §1.4  R = A/χ
  · 知识库 `水力学核心公式速查.md` §4.5  C=(1/n)R^(1/6)；§4.7 Q=(1/n)A·R^(2/3)·i^½
  · 知识库 `水力学核心公式速查.md` §6.7  J_f = n²v²/R^(4/3)（摩阻坡）
  · 吴持恭《水力学》第四版·下册 明渠非均匀流：均匀流 J_f = i ⇒ 上式即曼宁正常水深
  · 作者论文《推求法计算天然河道水面曲线的局限性和解决办法》：
    比降法 = 「用曼宁公式计算出水位～流量关系曲线」，水 位与流量一一对应。

三、口径裁决：说明书/教材公式 vs 权威 OUT 反演（重要）
------------------------------------------------------
说明书 [02] 与教材均给出**曼宁** C=(1/n)R^(1/6)。但以权威 `D-14B.OUT`
（GBK，SLSDK4.1，2010-09-25）反演，曼宁口径**不能复现**：

  断面   OUT水位Z   OUT面积W   曼宁解z    Δz      巴甫洛夫斯基解z   Δz
   1     103.65     1294.7     103.8698  +0.2244   103.6392        −0.0062
   2     103.53      833.34    103.7206  +0.1941   103.4572        −0.0693
   3     103.37      750.68    103.6487  +0.2793   103.3876        +0.0182
   4     102.94      694.58    103.0661  +0.1212   102.9745        +0.0296
   5     102.87      568.26    103.0323  +0.1637   102.9444        +0.0758
  （"z" 由 OUT 打印的 W 依本节断面几何反解，精度优于 ±0.0002 m；
    曼宁 max|Δz| = 0.279 m，巴甫洛夫斯基 max|Δz| = 0.076 m）

同族姊妹程序 **D-14A** 亦有同样现象，且更明确：
  D-14A.OUT 与 C=R^y/n（巴甫洛夫斯基）吻合 max|Δz|=0.019 m；
  与 C=R^(1/6)/n（曼宁）偏差 max|Δz|=0.198 m。
（复算脚本：`slcalc/_d14a_hf.py`、`slcalc/_d14_yscan.py`）

**裁决**：内核采用**谢才-巴甫洛夫斯基** C = R^y/n（y = 2.5√n − 0.13 − 0.75√R·(√n−0.10)），
即与 D-14A 同族一致的口径；说明书 [02] 的曼宁系数作为**备选口径**保留（`CHEZY_MODE='manning'`）。
裁决依据是权威 OUT（本项目既定惯例：公式有分歧时以 OUT 定案），并附上述反证表。

四、未闭合点（如实标注）
------------------------
1. 巴甫洛夫斯基层面仍有 ≤0.076 m 的逐位残差（断面 2 为 −0.069、断面 5 为 +0.076，
   断面 1/3/4 ≤0.030）。**反证**：对「湿周口径」3 种（水边相交段取部分长 / 取全段长 /
   不计）×「零宽竖直段」2 种（计入 / 剔除）×「谢才系数」2 种（曼宁 / 巴甫洛夫斯
   基），共 12 组组合逐组求解 5 个断面（脚本 `slcalc/_d14ab_joint.py`、`_d14b_grid.py`），
   无任何一组能把 5 个断面同时压到 0.005 m 以内；最优组（巴甫洛夫斯基 + 水边相交段
   取部分长 + 计入竖直段）即本内核口径。
   另对固定迭代次数二分（N=1..25）×3 种初值区间（脚本 `slcalc/_d14b_iter.py`）
   亦不能解释残差 —— 说明残差不是收敛容差造成的，本内核按"收敛到根"处理。
2. 原作者 2010 版 OUT 与 2003 版 INT 数据一致（OUT 回显的 19/18/18/18/17 个地形点
   与 INT 逐点相同），故残差不是数据版本问题。
3. OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-14B.out」为原著运行期路径回显
   （机器相关），内核不生成（与 D-1/D-3/D-4/D-7/D-8/D-11/D-16 同一口径）。
4. 原著在两岸地形点高程不足或纵坡太小时「提示出问题的断面序号并停机」；
   内核抛 `ValueError` 并给出断面号。

五、D-14A / D-14B / D-15 口径对照表（同族程序，勿混用）
--------------------------------------------------------
  项                D-14A（推求法）          D-14B（比降法）           D-15（回水/河道水面线）
  计算对象          逐段推求各断面水位        各断面独立正常水深       水库回水/河道水面线
  起始条件          需已知起始断面水位 Z0     不需起始水位             需坝前/下游控制水位
  基本方程          伯努利 Z1+α1V1²/2g        均匀流 Q=ω·C·√(R·i)     能量方程 + 摩阻比降
                    = Z2+α2V2²/2g+hf+hj      ⇒ Q=(1/n)ωR^(2/3)√i
  沿程损失 hf       Q²·ΔL/(K1·K2)（K 几何平均） 不存在（每断面独立）    ΔL·(i1+i2)/2（摩阻比降算术平均）
  谢才系数 C        R^y/n（巴甫洛夫斯基层面） 同左（由 OUT 反演裁决）   见 d15.py
  局部损失 hj       ξ(V1²−V2²)/2g            不存在                   见 d15.py
  输入              桩号/ARF/GG/M + 点         桩号/**纵坡 i**/M + 点    见 d15.py
  二分目标准则      能量残差 K=0               正常水深残差 F(h)=0      见 d15.py
  断面几何          分块累加（梯形/三角形）    同左（同一几何模块）      同左
  输出              水位/面积/流速/横断面图    水位/面积/流速/横断面图   水位/面积/流速

参考：`slcalc/programs/d14a.py`、`slcalc/programs/d15.py`。
"""
import math

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-14B"
TITLE = "天然河道横断面正常水深计算书"
AUTHOR = "张校正（新疆水利厅）"

G = 9.81

#: 谢才系数口径：'pavlovsky'（默认，按权威 OUT 裁决）| 'manning'（说明书 §[02]）
CHEZY_MODE = "pavlovsky"

#: 断面几何分块口径
PER_RULE = "part"      # 水边相交段湿周：'part' 取水下部分 | 'full' 取全段长
SKIP_VERT = False      # 是否剔除零宽竖直段（ΔX=0）对湿周的贡献


# ============================================================
# 断面几何（说明书 (1)(2)(3) 式 / 教材 §1.4）
# ============================================================

def area_wetted(z, xy, per_rule=PER_RULE, skip_vert=SKIP_VERT):
    """
    分块累加过水面积 ω 与湿周 χ（水面 z 截断，梯形 / 三角形）。
    xy: [(x,y), ...] 地形点序列（相邻点连线段）。

    说明书 (1)(2)(3)：
      (1) Z≥Yi 且 Z≥Yi+1：W_i=(Z−(Y_i+Y_{i+1})/2)·ΔX，S_i=√(ΔX²+ΔY²)
      (2) Y_i<Z<Y_{i+1}：W_i=ΔX·(Z−Y_i)²/(2ΔY)，S_i=(Z−Y_i)·√(ΔX²+ΔY²)/ΔY
      (3) Y_i>Z>Y_{i+1}：把 (2) 中 Y_i、Y_{i+1} 互换即可
    返回 (ω, χ)。
    """
    area = 0.0
    perim = 0.0
    for i in range(len(xy) - 1):
        x1, y1 = xy[i]
        x2, y2 = xy[i + 1]
        dx = x2 - x1
        dy = y2 - y1
        sl = math.hypot(dx, dy)
        vert = abs(dx) < 1e-12
        ylo, yhi = (y1, y2) if y1 < y2 else (y2, y1)
        if z <= ylo:
            continue                                  # 整段在岸上
        if z >= yhi:                                  # 情况 (1) 全淹
            area += (z - 0.5 * (y1 + y2)) * dx
            if not (skip_vert and vert):
                perim += sl
        else:                                         # 情况 (2)(3) 水边相交
            t = (z - y1) / dy if dy != 0 else 0.0
            if y1 <= z:
                area += 0.5 * (t * dx) * (z - y1)
                frac = t
            else:
                area += 0.5 * ((1 - t) * dx) * (z - y2)
                frac = 1 - t
            if skip_vert and vert:
                pass
            elif per_rule == "part":
                perim += frac * sl
            elif per_rule == "full":
                perim += sl
    return area, perim


# ============================================================
# 谢才系数
# ============================================================

def pavlovsky_y(n, R):
    """巴甫洛夫斯基指数 y = 2.5√n − 0.13 − 0.75√R·(√n − 0.10)。"""
    sn = math.sqrt(n)
    return 2.5 * sn - 0.13 - 0.75 * math.sqrt(R) * (sn - 0.10)


def chezy_C(n, R, mode=CHEZY_MODE):
    """谢才系数。mode='pavlovsky'：C=R^y/n；mode='manning'：C=R^(1/6)/n。"""
    if R <= 0:
        return 0.0
    if mode == "manning":
        return R ** (1.0 / 6.0) / n
    return R ** pavlovsky_y(n, R) / n


def section_hydraulics(xy, n, z, Q, mode=CHEZY_MODE, per_rule=PER_RULE,
                       skip_vert=SKIP_VERT):
    """断面水力要素：ω、χ、R、V、C、K=ω·C·√R（流量模数）。"""
    W, P = area_wetted(z, xy, per_rule, skip_vert)
    R = W / P if P > 0 else 0.0
    V = Q / W if W > 0 else 0.0
    C = chezy_C(n, R, mode) if R > 0 else 0.0
    K = W * C * math.sqrt(R) if R > 0 else 0.0
    return {"W": W, "P": P, "R": R, "V": V, "C": C, "K": K}


# ============================================================
# 正常水深（比降法）求解
# ============================================================

def normal_depth(sec, Q, mode=CHEZY_MODE, per_rule=PER_RULE, skip_vert=SKIP_VERT):
    """
    二分法求断面正常水深水位 Z，使流量模数满足 Q = K(Z)·√i（说明书式[00]⇒[04]）。

    K(Z) = ω(Z)·C(n,R)·√(R) 随 Z 单调递增，二分区间 [min(Y)+ε, max(Y)+Δ]。
    若在 max(Y)+Δ 处仍 K√i < Q（两岸地形点高程不够 / 纵坡太小），
    按原著"提示出问题的断面序号并停机"，抛 ValueError。
    """
    xy = sec["xy"]
    i = sec["i"]
    n = sec["n"]
    ymin = min(y for _, y in xy)
    ymax = max(y for _, y in xy)

    def Ksqrt(z):
        return section_hydraulics(xy, n, z, Q, mode, per_rule, skip_vert)["K"] * math.sqrt(i)

    lo = ymin + 1e-6
    hi = ymax + 5.0
    guard = 0
    while Ksqrt(hi) < Q and guard < 60:               # 自适应抬高上界
        hi += max(5.0, (hi - lo) * 0.5)
        guard += 1
    if Ksqrt(hi) < Q:
        raise ValueError("第 %s 个断面（桩号 %s）计算不出正常水深："
                         "两岸地形点高程不够或纵坡太小，请检查数据！"
                         % (sec.get("order", "?"), sec.get("ll", "?")))
    if Ksqrt(lo) >= Q:
        raise ValueError("第 %s 个断面（桩号 %s）计算不出正常水深："
                         "最低地形点以下即已满足流量，数据可能有误！"
                         % (sec.get("order", "?"), sec.get("ll", "?")))

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if Ksqrt(mid) < Q:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10:
            break
    return 0.5 * (lo + hi)


def normal_depth_manning(sec, Q):
    """备选口径：曼宁正常水深（说明书式[02] C=R^(1/6)/n），供对照/回归。"""
    return normal_depth(sec, Q, mode="manning")


# ============================================================
# 解析
# ============================================================

def parse(data):
    """
    解析 D-14B 输入。data: INT 文件路径 | dict。

    INT 固定顺序（说明书"三、输入数据顺序"）：
      HE LIOU, MM, Q                    第一行：河流名称、横断面个数、流量
      以下每断面：
        LL, II, M                       桩号(米)、本段纵坡、地形点个数
        X(k), Y(k), N(k)  ×M            横坐标、绝对高程、该点糙率
    返回 dict：river, n_sections, Q,
              sections=[{ll, i, npts, xy, n_list, n}, ...]
    """
    if isinstance(data, dict):
        return data

    lines = read_lines(data)
    lines = [ln for ln in lines if ln.strip() != ""]
    if not lines:
        raise ValueError("INT 文件为空")

    first = lines[0].split(",")
    river = first[0].strip().strip('"')
    n_sec = int(float(first[1]))
    Q = float(first[2])

    sections = []
    pos = 1
    for s in range(n_sec):
        if pos >= len(lines):
            raise ValueError("INT 文件数据不完整（断面数不足）")
        parts = [p for p in lines[pos].split(",") if p.strip() != ""]
        ll = float(parts[0])
        i = float(parts[1])
        M = int(float(parts[2]))
        pos += 1
        xy, ns = [], []
        for _ in range(M):
            a, b, c = lines[pos].split(",")[:3]
            pos += 1
            xy.append((float(a), float(b)))
            ns.append(float(c))
        sections.append({
            "order": s + 1, "ll": ll, "i": i, "npts": M,
            "xy": xy, "n_list": ns,
            "n": ns[0] if ns else 0.03,       # 主糙率（各点相同时取首点）
        })
    return {"river": river, "n_sections": n_sec, "Q": Q, "sections": sections}


# ============================================================
# 计算
# ============================================================

def compute(params, mode=CHEZY_MODE):
    """逐断面独立求解正常水深（比降法）。"""
    Q = params["Q"]
    sections = params["sections"]

    result = {
        "程序": PROGRAM_ID,
        "标题": TITLE,
        "作者": AUTHOR,
        "河流名称": params.get("river", ""),
        "流量Q(m3/s)": Q,
        "断面数": len(sections),
        "谢才系数口径": mode,
        "断面": [],
    }
    for sec in sections:
        z = normal_depth(sec, Q, mode=mode)
        h = section_hydraulics(sec["xy"], sec["n"], z, Q, mode)
        result["断面"].append({
            "断面号": sec["order"],
            "桩号(m)": sec["ll"],
            "纵坡i": sec["i"],
            "地形点数": sec["npts"],
            "糙率": sec["n"],
            "水位Z(m)": z,
            "过水面积W(m2)": h["W"],
            "湿周P(m)": h["P"],
            "水力半径R(m)": h["R"],
            "平均流速V(m/s)": h["V"],
            "谢才系数C": h["C"],
            "流量模数K": h["K"],
        })
    return result


# ============================================================
# 输出（复刻原著 .OUT 版式）
# ============================================================

def _fmt_trim(x, nd):
    """按 nd 位小数格式化并去掉尾部多余 0（复刻 VB 打印）。"""
    s = ("%%.%df" % nd) % x
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s or "0"


def _fmt_i(x):
    """纵坡：VB Str() 风格（0.00045 → '.00045'）。"""
    s = "%g" % x
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return s


def render(params, result):
    """
    生成原著风格中文计算书——**逐行复刻** D-14B.OUT 版面
    （首行「文件：…」为原著运行期路径回显，本内核不生成）。

    版式要点（由权威 OUT 逐行反演）：
      · 抬头 8 行：空行 + 3 行星号框 + 空行 + 河流名称行 + 2 空行
      · 每断面 12 行 + ceil(M/2) 行地形点（左右双列，列距 12 空格）
        ——M 为**偶数**时地形点行末尾多一个空行（原著 Print 行为）
      · 断面之间 2 个空行；纵坡字段左对齐占宽 17
    """
    secs = params["sections"]
    rows = result["断面"]
    L = []
    L.append("")
    L.append(" **********************************************************************")
    L.append(" *****             天然河道横断面正常水深计算书 D-14B             *****")
    L.append(" **********************************************************************")
    L.append("")
    L.append("              河流名称:%s     流量 Q= %s " % (result["河流名称"],
                                                      _fmt_trim(result["流量Q(m3/s)"], 6)))
    L.append("")
    L.append("")
    for k, r in enumerate(rows):
        sec = secs[k]
        if k > 0:
            L.append("")
            L.append("")
        L.append("                        第 %d 个断面:%s" % (r["断面号"],
                                                          _fmt_trim(r["桩号(m)"], 6)))
        L.append("")
        L.append("              纵坡:  i= %s地形点个数:  M= %d "
                 % (_fmt_i(r["纵坡i"]).ljust(17), r["地形点数"]))
        L.append("")
        L.append("                  %d 个地形点的坐标及横断面处糙率:" % r["地形点数"])
        L.append("                 ==================================")
        L.append(" NO     X       Y      n                NO     X       Y      n")
        xy = sec["xy"]
        ns = sec["n_list"]
        M = len(xy)
        half = (M + 1) // 2
        for a in range(half):
            left = "%3d%8.2f%8.2f%8.4f" % (a + 1, xy[a][0], xy[a][1], ns[a])
            b = a + half
            if b < M:
                right = "%3d%8.2f%8.2f%8.4f" % (b + 1, xy[b][0], xy[b][1], ns[b])
                L.append(left + " " * 12 + right)
            else:
                L.append(left)
        if M % 2 == 0:
            L.append("")
        L.append("                             计算结果")
        L.append("                             ========")
        L.append("")
        L.append(" 水位:  Z= %s     过水断面积:  W= %s     平均流速:  V= %s"
                 % (_fmt_trim(r["水位Z(m)"], 2), _fmt_trim(r["过水面积W(m2)"], 2),
                    _fmt_trim(r["平均流速V(m/s)"], 3)))
        L.append("-" * 69)
    return "\n".join(L)


def run(data, out_txt=None, out_json=None, fmt="text", mode=CHEZY_MODE):
    """统一入口。data：INT 文件路径 | dict。"""
    params = parse(data)
    result = compute(params, mode=mode)
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
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
