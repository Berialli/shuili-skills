# -*- coding: utf-8 -*-
"""
E-4 光电仪液限、塑限联合试验计算程序 —— 内核
=============================================
复刻《水利水电工程设计计算程序集》E-4 程序（作者：邓铭江，新疆水利厅）。

原著说明书（E-4Intro）：
  「适用光电式液、塑限联合试验仪测定土体液、塑限的计算。计算液、塑限时本程序将
    《规程》的图解过程（图4-1）完全用解析法代替。…程序的打印成果有：试验误差、
    塑限、液限、塑性指标。」
  打印字符：ωb、ωc —— 圆锥入土深度为 2mm 时的两个含水量（图4-1 中 b、c 两点）；
            |ωb−ωc| —— 本试验要求 ≤2%，否则补点试验；
            ωp、ωL —— 塑限、液限；Ip —— 塑性指数。

算法（《土工试验方法标准》GB/T 50123 液塑限联合测定法图解过程的解析化）
----------------------------------------------------------------------
  1. 试验点 (h(i), ω(i))，i=1,2,3（h = 圆锥入土深度 mm，ω = 含水量 %）。
  2. 坐标系：双对数（横轴 log h，纵轴 log ω），三点在双对数坐标上近似呈直线。
  3. 「通过较高含水量的点与其余两点连成两条直线」，在 h=2mm 处查得两个含水量
        ωb = 直线(<最高 ω 点>, <余点1>) 在 h=2 处的 ω
        ωc = 直线(<最高 ω 点>, <余点2>) 在 h=2 处的 ω
  4. 试验误差 |ωb−ωc|（要求 ≤2%）。
  5. 「以该两点含水量的平均值与较高含水量点连成一直线」——
        塑限 ωp = (ωb+ωc)/2 ；该直线在 h=10mm 处查得
        液限 ωL
  6. 塑性指数 Ip = ωL − ωp

INPUT（E-4.INT 文本）：三行，每行 h, ω
OUTPUT（逐字复刻权威 E-4.OUT，GBK）
    首行「文件：…E-4.out」为运行期路径回显，本内核不生成。

常量口径
--------
  · 纵/横轴对数底数 = 10（双对数坐标纸），对直线拟合而言底数可任意（本题与底数无关）。
  · 塑限判据深度 h_p = 2.0 mm；液限判据深度 h_L = 10.0 mm
    （由权威 OUT 逐位反演：h=10 处拟合外推 33.5696→33.57=ωL；h=2 处 ωp=22.4207→22.42）。

未闭合点（如实标注，量化）
--------------------------
  权威 E-4.OUT 中 ωb=22.306（3 位小数）、|ωb−ωc|=0.2294414（7 位有效）与本内核
  ωb=22.306911（→22.307）、|ωb−ωc|=0.2294624 不符，Δ=9.1e-4 与 2.1e-5。
  反证与排查：
    ① 该偏差是 ωb、ωc **共同**的 −4.0e-5 相对下移（二者差 0.22944/0.22946 亦随之偏移），
       与「评估深度 h0 略小于 2mm」等价：由 ωb、ωc 分别反解得 h0=1.99968、1.99967 mm
       （**两条直线一致到 3e-6**，故 h0 口径可分辨，唯一公因子）。
    ② 已逐一排除的替代解释（单参数反解，均两条直线不自治）：
       最高点深度 h1 → 需 9.79966 / 9.79674（不自洽）；
       最高点含水量 ω1 → 超出 [32,35] 搜索域无解；
       评估深度以外各点 h、ω → 或触界或两条线不自洽；
       改用 lin-lin、半对数、ln-ln、x/y 互换、最小二乘全 3 点直线、
       (P1,P2)+(P1,P3)、(P1,P3)+(P2,P3) 等换线组合 → 均无法命中 ωb。
    ③ 结论：ωb（3 位小数）与 |ωb−ωc| 记 DECL；其余 ωc=22.54、ωp=22.42、ωL=33.57、
       Ip=11.15 四项与权威**逐字命中**（ωc 的打印半宽较宽，掩盖了同一相对偏移）。

知识库对照结果（教学母本 → 程序实现）
------------------------------------
  KB 路径：D:\\WorkBuddy知识库\\水利知识库\\

  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 液塑限联合测定：h–ω 关系（双对数）直线 | `02_水利教材精读/土质学与土力学/土质学与土力学_第5版_精读笔记.md` 界限含水量章：光电式液塑限联合测定仪，圆锥入土深度与含水量在双对数坐标上为直线 | **一致** | 采用 |
  | 塑限 = h 2mm 对应含水量 | 同上（76g 锥：塑限对应入土深度 2mm） | **一致** | 采用 h_p=2.0 |
  | 液限 = h 10mm 对应含水量 | 同上（锥式仪液限对应入土深度 10mm） | **一致** | 采用 h_L=10.0 |
  | |ωb−ωc|≤2% 判据 | GB/T 50123 液塑限联合测定法（说明书原文引述《规程》） | **一致** | 采用 |
  | Ip = ωL − ωp | 同上（塑性指数定义） | **一致** | 采用 |
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-4"
TITLE = "光电仪液限.塑限联合试验计算"
AUTHOR = "邓铭江（新疆水利厅）"
HEAD_NAME = "E-4"

H_PLASTIC = 2.0    # 塑限判据深度 mm
H_LIQUID = 10.0    # 液限判据深度 mm


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def vbnum(x, nd=7):
    """VB6 Format(x,"0.#######") 风格：最多 nd 位小数、去尾零。"""
    s = "%.*f" % (nd, q(x, nd))
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s or "0"


# ============================================================
# 解析
# ============================================================

def parse(data):
    if isinstance(data, dict):
        pts = [(float(a), float(b)) for a, b in data["points"]]
        return {"points": pts}
    nums = read_numbers(data)
    if len(nums) < 6:
        raise ValueError("E-4 需要 3 组 (h, ω)，实际数值 %d 个" % len(nums))
    pts = [(nums[2 * i], nums[2 * i + 1]) for i in range(len(nums) // 2)]
    return {"points": pts}


# ============================================================
# 计算
# ============================================================

def w_at(p, q, h0):
    """双对数坐标下过 p、q 两点的直线在入土深度 h0 处的含水量。"""
    (h1, w1), (h2, w2) = p, q
    x1, x2 = math.log10(h1), math.log10(h2)
    y1, y2 = math.log10(w1), math.log10(w2)
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    return 10.0 ** (m * math.log10(h0) + b)


def w_at_xy(x1, y1, x2, y2, h0):
    """任意两点（已取对数）确定的直线在 h0 处的含水量。"""
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    return 10.0 ** (m * math.log10(h0) + b)


def compute(p):
    pts = p["points"]
    if len(pts) != 3:
        raise ValueError("E-4 液塑限联合试验需 3 组 (h, ω)，实际 %d 组" % len(pts))
    # 「较高含水量的点」
    top = max(pts, key=lambda t: t[1])
    others = [t for t in pts if t is not top]
    wb = w_at(top, others[0], H_PLASTIC)     # 直线(top, 余点1) 与 h=2mm 交点
    wc = w_at(top, others[1], H_PLASTIC)     # 直线(top, 余点2) 与 h=2mm 交点
    err = abs(wb - wc)
    wp = (wb + wc) / 2.0                     # 平均点与较高含水量点连线；h=2mm 处即 wp
    wl = w_at_xy(math.log10(H_PLASTIC), math.log10(wp),
                 math.log10(top[0]), math.log10(top[1]), H_LIQUID)
    ip = wl - wp
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": p, "最高点": top,
            "wb": wb, "wc": wc, "误差": err, "塑限": wp,
            "液限": wl, "塑性指数": ip}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, r):
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****             光电仪液限.塑限联合试验计算 E-4                   ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("   园锥入土深度(mm)     对应的含水量(%)   ")
    for h, w in p["points"]:
        A.append("%6.1f%24.2f" % (q(h, 1), q(w, 2)))
    A.append(" wb= %.3f(%%)" % q(r["wb"], 3))
    A.append(" wc= %.2f(%%)" % q(r["wc"], 2))
    A.append(" |wb-wc|=" + vbnum(r["误差"], 7) + "<=2 (%)")
    A.append(" 塑限 wp= %.2f(%%)" % q(r["塑限"], 2))
    A.append(" 液限 wl= %.2f%%" % q(r["液限"], 2))
    A.append(" 塑性指数 Ip= %.2f" % q(r["塑性指数"], 2))
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
