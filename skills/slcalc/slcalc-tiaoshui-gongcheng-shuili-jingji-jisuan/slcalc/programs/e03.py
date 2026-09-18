# -*- coding: utf-8 -*-
"""
E-3 用卡式碟式仪测定土体液限计算程序 —— 内核
=============================================
复刻《水利水电工程设计计算程序集》E-3 程序（作者：尤红霞，新疆水利水电学校）。

原著说明书（E-3Intro）：
  「本程序根据卡式碟式仪所测得的试验数据即击次 N 和各击数所对应的含水量 ω。
    计算机以含水量 ω 为纵坐标，击次的对数 logN 为横坐标绘出坐标图…然后计算机
    应用最小二乘法原理拟合各试验点的直线方程。并根据此方程求出击数为 25 次时
    所对应的含水量即液限。这样就将规程《界限含水量试验说明》中的图解过程变成
    了解析法。」

算法（唯一，无分支）
--------------------
  1. 读入 M 组 (N(i), ω(i))，i=1..M（说明书要求 N 由小到大输入）。
  2. 取 x(i) = log10(N(i))（十进制对数），y(i) = ω(i)，作一元线性最小二乘
        k = (n·Σxy − Σx·Σy) / (n·Σx² − (Σx)²)
        c = (Σy − k·Σx) / n
  3. 液限 ω_L 取 N=25 处的拟合值：ω_L = k·log10(25) + c

INPUT（E-3.INT 文本）
---------------------
  M
  N(1),ω(1)
  …
  N(M),ω(M)

OUTPUT（逐字复刻权威 E-3.OUT，GBK）
------------------------------------
  首行「文件：…E-3.out」为原著运行期路径回显（机器相关），本内核不生成。

常量口径（由权威 OUT 反演）
--------------------------
  · 对数底数 = 10（说明书「logN」「LOG(N)」字样 + 权威 OUT 方程系数 −8.889 量级互证）。
  · 无其它常数（纯最小二乘）。

未闭合点（如实标注，量化）
--------------------------
  权威 E-3.OUT 第 16 行方程系数 k=-8.88922、c=58.84124，与本内核按双精度最小二乘
  所得 k=-8.8891215、c=58.8411818 相差 Δk=9.9e-5、Δc=5.8e-5（打印末位单位的
  ~10 倍 / ~6 倍）。反证与排查：
    ① 该二系数**互相不自洽**：若 k=-8.88922，则 c 应为 ȳ−k·x̄ = 58.84133；
       若 c=58.84124，则 k 应为 (ȳ−c)/x̄ = −8.88916。二者差 9e-5，说明原著此处
       并非常规最小二乘的同一组解（或经过单精度/分行刷值）。
    ② 已逐一排除：float32 全程累加（k=-8.889137）、纯偏差式公式、在线递归最小二乘、
       升/降序、x 以 Single 存、对数底数常数扫描 D∈[2.3020,2.3024]（最佳仍差 8.7e-4）、
       x 四舍五入至 3/5/6 位小数、y 整体平移（不影响 k）。均不能同时命中 k 与 c。
    ③ 结论：差异仅出现在**方程系数回显**上；最后的液限 ω_L=46.41 与 3 位有效数字的
       k、c 行（第 13、14 行）在两种口径下**逐字相同**。故 k、c 的 5 位小数回显记 DECL。

知识库对照结果（教学母本 → 程序实现）
------------------------------------
  KB 路径：D:\\WorkBuddy知识库\\水利知识库\\

  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 液限 = ω(N=25) | `02_水利教材精读/土质学与土力学/土质学与土力学_第5版_精读笔记.md` 界限含水量章：卡式碟式仪液限以 25 击为标准（击数-含水量关系图解至 25 击） | **一致** | 采用 |
  | ω 与 logN 呈直线 | 同上「以含水量为纵坐标、击次对数为横坐标绘制关系曲线」 | **一致** | 采用（log10） |
  | 最小二乘拟合直线 | 同上「应用最小二乘法原理拟合各试验点的直线方程」 | **一致** | 采用 |
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-3"
TITLE = "用卡式碟式仪测定土体液限计算"
AUTHOR = "尤红霞（新疆水利水电学校）"
HEAD_NAME = "E-3"

N_STD = 25.0          # 液限标准击次


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    return "%*.*f" % (w, d, q(x, d))


# ============================================================
# 解析
# ============================================================

def parse(data):
    """解析输入。data：dict | .INT 文件路径。"""
    if isinstance(data, dict):
        d = dict(data)
        pts = [(float(a), float(b)) for a, b in d["points"]]
        return {"M": len(pts), "points": pts}
    nums = read_numbers(data)
    if len(nums) < 1:
        raise ValueError("E-3 输入为空")
    m = int(round(nums[0]))
    rest = nums[1:]
    if len(rest) < 2 * m:
        raise ValueError("E-3 输入数据不足：M=%d，数值 %d 个" % (m, len(rest)))
    pts = [(rest[2 * i], rest[2 * i + 1]) for i in range(m)]
    return {"M": m, "points": pts}


# ============================================================
# 计算
# ============================================================

def compute(p):
    pts = p["points"]
    n = len(pts)
    if n < 2:
        raise ValueError("E-3 至少需要 2 组试验点")
    xs = [math.log10(a) for a, _ in pts]
    ys = [b for _, b in pts]
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(a * a for a in xs)
    sxy = sum(a * b for a, b in zip(xs, ys))
    den = n * sxx - sx * sx
    if den == 0.0:
        raise ValueError("E-3 试验点击次全相同，无法拟合直线")
    k = (n * sxy - sx * sy) / den
    c = (sy - k * sx) / n
    n_std = N_STD
    wl = k * math.log10(n_std) + c
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": p, "点数": n, "k": k, "c": c, "n_std": n_std, "液限": wl}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, r):
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****              用卡式碟式仪测定土体液限计算  E-3                ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("         试验记录")
    A.append("  编  号    击次(N)      含水量w(%)   ")
    for i, (nn, ww) in enumerate(p["points"]):
        A.append("%5d%9d%15.3f" % (i + 1, int(round(nn)), q(ww, 3)))
    A.append(" k= %.3f" % q(r["k"], 3))
    A.append(" c= %.3f" % q(r["c"], 3))
    A.append("")
    A.append(" w= %.5f*(LOG(N))+ %.5f" % (q(r["k"], 5), q(r["c"], 5)))
    A.append(" n=%4d" % int(round(r["n_std"])))
    A.append(" w= %.2f(%%)" % q(r["液限"], 2))
    A.append("")
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
