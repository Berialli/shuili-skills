# -*- coding: utf-8 -*-
"""
E-7 压缩试验计算程序 —— 内核
==============================
复刻《水利水电工程设计计算程序集》E-7（作者：邓铭江，新疆水利厅）。

功能
----
土的侧限压缩试验成果整理：由各级荷重 P 下试样总变形 hz 计算压缩后高度 h、
单位沉降量 Si、孔隙比 ei、压缩系数 av、压缩模量 Es、最大排水距离 hp、
压缩指数 Cc，并据 av(1~2) 判定土的压缩性。

公式（E-7Intro「三、显示及打印字符说明」+ 权威 E-7.OUT 逐位反演）
------------------------------------------------------------------
    γd   = r / (1 + ω/100)                  （干容重）
    e0   = Gs·γw / γd − 1  =  Gs(1+ω/100)/r − 1        （初始孔隙比）
    h    = h0 − hz                          (mm)
    Si   = hz/h0 · 1000                     (mm/m)
    ei   = e0 − (1+e0)·hz/h0
    av_i = (e_{i−1} − e_i)/(P_i − P_{i−1})   (cm²/kg；i=1 时 e_0=e0、P_0=0)
    Es_i = (1 + e0)/av_i                    (kg/cm²)
    hp_i = (h_{i−1} + h_i)/4/10             (cm)  ← 试样最大排水距离（=平均高度之半）
    Cc   = −d(e)/d(lg p)（对交互选定的直线段点作最小二乘）

未闭合点 / 口径说明
------------------
  1. **Cv 列为空**：权威 OUT 表头含 Cv(cm²/sec) 列但算例各行皆空 —— 原著仅当
     另给 t90 之外的数据时才输出，本内核照排留空，不做猜测。
  2. **av==0 时 Es 打印 −1.000**：权威 OUT 第 6 行（P: 8→4，卸荷零变形）Es=−1.000，
     即除零保护值；本内核照排。
  3. **压缩性判定阈值**（cm²/kg）：av<0.01 低 / 0.01≤av<0.05 中 / av≥0.05 高。
     权威算例 av(1~2)=0.059 与 E-8 的 0.0647 均判"高压缩性土"，与 ≥0.05 口径一致；
     低/中两档阈值库中无算例标定，按《土工试验规程》旧制登记（**推断，非反演**）。
  4. **Cc 的选点序号**指"加荷支"（P 单调上升段，不含 P=0 初始行）中的先后次序；
     权威算例选 3,4,5 → 对应 P=2,4,8。若按整表行号(3,4,5→P=1,2,4)则 Cc=0.2128，
     与权威 0.2447 不符 → 已排除（相对差 13.0%）。
  5. 权威 .OUT 首行「文件：…」为运行期路径回显，本内核不生成。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 三相换算 e = Gs(1+ω)ρw/ρ − 1 | `01_水工设计手册精读/卷1_基础理论.md` 4.1「基本物理指标 … e=Gs(1+ω)ρw/ρ−1」；`02_水利教材精读/水能规划与岩土结构/土质学与土力学_第5版_精读笔记.md` 表「孔隙比 e=Gsρw(1+w)/ρ−1」 | **一致** | 采用 |
  | 压缩指标 av、Es 与 e-p 曲线 | 卷1 4.6「土的压缩性（压缩试验/压缩指标…）」；教材压缩性章 | **一致** | 采用 |
  | 压缩指数 Cc：e-lgp 曲线直线段斜率 | 卷1 4.6/4.10「e-lgp 曲线」；`卷6_土石坝.md`「压缩系数」段 | **一致** | 采用（最小二乘） |
  | Es=(1+e1)/av 与库中 Es=(1+e0)/av 之别 | 教材 Es=(1+e1)/a | 程序用 **(1+e0)**（起始孔隙比）而非 (1+e1)，差异 | 采用程序口径（反演命中 11 行） |

输入数据（E-7Intro「数据文件为 E-7.INT」）
-------------------------------------------
    第 1 行： ω,Gs,r,H0,M         （ω %、Gs、r kg/cm³、H0 mm、M 加卸荷总次数）
    其后 M+1 行： P,hz,t90        （P kg/cm²、hz mm 总变形、t90 min；未测者填 0）
    另需 Cc 直线段选点序号（原著交互输入；本例 3,4,5）。

输出（逐字复刻原著 .OUT 版式）
------------------------------
    表头 → 参数 → P/hz/h/Si/ei/t90/Es/av/hp/Cv 表 → av(1~2)、压缩性、Cc
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-7"
TITLE = "压缩试验计算程序"
AUTHOR = "邓铭江"
HEAD_NAME = "E-7"

RULE = "_" * 96
DEFAULT_CC_PTS = [3, 4, 5]

# 压缩性判定阈值（cm²/kg，见 docstring 第 3 条）
AV_LOW = 0.01
AV_HIGH = 0.05


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    return "%*.*f" % (w, d, q(x, d))


def G(x):
    """整数则不打小数点（含水量回显口径）。"""
    if abs(x - round(x)) < 1e-9:
        return "%d" % int(round(x))
    return "%g" % x


# ============================================================
# 解析
# ============================================================

def parse(data):
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("Cc选点", DEFAULT_CC_PTS)
        return p
    nums = read_numbers(data)
    if len(nums) < 5:
        raise ValueError("E-7 输入数据过少（%d 个数）" % len(nums))
    w, Gs, r, H0, M = nums[:5]
    rest = nums[5:]
    rows = []
    for i in range(0, len(rest) - 2, 3):
        rows.append((rest[i], rest[i + 1], rest[i + 2]))
    return {"w": w, "Gs": Gs, "r": r, "H0": H0, "M": int(round(M)),
            "rows": rows, "Cc选点": DEFAULT_CC_PTS}


# ============================================================
# 计算
# ============================================================

def compute(p):
    w, Gs, r, H0 = p["w"], p["Gs"], p["r"], p["H0"]
    rows = p["rows"]
    e0 = Gs * (1.0 + w / 100.0) / r - 1.0
    tab = []
    # 第 1 行为初始状态（P=0），其 t90/Es/av/hp 一律留空
    for i, (P, hz, t90) in enumerate(rows):
        h = H0 - hz
        Si = hz / H0 * 1000.0
        ei = e0 - (1.0 + e0) * hz / H0
        tab.append({"P": P, "hz": hz, "h": h, "Si": Si, "ei": ei, "t90": t90,
                    "Es": None, "av": None, "hp": None})
    for i in range(1, len(tab)):
        prev, cur = tab[i - 1], tab[i]
        dP = cur["P"] - prev["P"]
        de = prev["ei"] - cur["ei"]
        av = de / dP if dP != 0 else 0.0
        if av == 0.0:
            av = 0.0                       # 归一化 -0.0 → +0.0（打印口径）
        Es = -1.0 if av == 0.0 else (1.0 + e0) / av
        cur["av"] = av
        cur["Es"] = Es
        cur["hp"] = (prev["h"] + cur["h"]) / 4.0 / 10.0

    # av(1~2)：P=1 → P=2 的压缩系数（即 P=2 那一级的 av）
    av12 = None
    for i, row in enumerate(tab):
        if i >= 1 and abs(row["P"] - 2.0) < 1e-9:
            av12 = row["av"]          # 取首次出现（加荷支）的 P=2 级
            break
    if av12 is None:
        av12 = 0.0
    if av12 < AV_LOW:
        soil = "低压缩性土"
    elif av12 < AV_HIGH:
        soil = "中压缩性土"
    else:
        soil = "高压缩性土"

    # Cc：对加荷支选定序号作 e ~ lg p 最小二乘
    pmax = max(row["P"] for row in tab)
    Cc = None
    r0 = None
    if pmax > 4.0:
        # 加荷支 = 第 2 行到 P 首次取最大值的那一行
        imax = next(i for i, row in enumerate(tab) if abs(row["P"] - pmax) < 1e-12)
        branch = list(range(1, imax + 1))
        sel = [branch[j - 1] for j in p["Cc选点"]
               if 1 <= j <= len(branch)]
        xs = [math.log10(tab[i]["P"]) for i in sel]
        ys = [tab[i]["ei"] for i in sel]
        n = len(xs)
        mx = sum(xs) / n
        my = sum(ys) / n
        sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        sxx = sum((x - mx) ** 2 for x in xs)
        syy = sum((y - my) ** 2 for y in ys)
        slope = sxy / sxx
        Cc = -slope
        r0 = sxy / math.sqrt(sxx * syy)
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": p, "e0": e0, "表": tab, "av12": av12, "压缩性": soil,
            "Cc": Cc, "r0": r0}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, res):
    B = []
    B.append("")
    B.append(" ***********************************************************************")
    B.append(" ****                     压缩试验计算程序 E-7                      ****")
    B.append(" ***********************************************************************")
    B.append("")
    B.append("")
    B.append("    工程名称_____________             试验者___________")
    B.append("")
    B.append("    土样编号_____________             计算者___________")
    B.append("")
    B.append("    试验日期_____________             校核者___________")
    B.append(RULE)
    B.append("")
    B.append("    试样含水量w= %s %%                 试样比重Gs= %s"
             % (G(p["w"]), F(p["Gs"], 0, 3)))
    B.append("    试样容重r= %s(kg/cm^3)         试样原始高度h0= %smm"
             % (F(p["r"], 0, 3), F(p["H0"], 0, 1)))
    B.append(RULE)
    B.append("")
    B.append("    P        hz       h           Si        ei        t90       Es"
             "        av        hp        Cv")
    B.append(" (kg/cm^2)  (mm)     (mm)       (mm/m)     (-)       (min)    (kg/cm^2)"
             " (cm^2/kg)  (cm)    (cm^2/sec)")
    for i, row in enumerate(res["表"]):
        s = F(row["P"], 6, 2) + F(row["hz"], 10, 3) + F(row["h"], 10, 3) \
            + F(row["Si"], 12, 3) + F(row["ei"], 10, 3)
        if i == 0:
            s = s.ljust(87)
        else:
            s += (F(row["t90"], 10, 3) if row["t90"] > 0 else " " * 10)
            s += F(row["Es"], 10, 3) + F(row["av"], 10, 3) + F(row["hp"], 9, 3)
        B.append(s)
    B.append(RULE)
    B.append("")
    B.append(" av(1~2)= %s(cm^2/kg)" % F(res["av12"], 0, 3))
    B.append(" %s" % res["压缩性"])
    if res["Cc"] is not None:
        B.append(" 压缩指数Cc= %s" % F(res["Cc"], 0, 4))
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
