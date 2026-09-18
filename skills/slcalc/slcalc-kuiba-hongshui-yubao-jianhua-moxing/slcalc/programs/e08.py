# -*- coding: utf-8 -*-
"""
E-8 快速压缩试验计算程序 —— 内核
==================================
复刻《水利水电工程设计计算程序集》E-8（作者：孟建平，新疆水利水电学校）。

功能
----
计算试样初始孔隙比 e0、各级荷重下（固结 1h）经"快速法校正"后的总变形量、
单位沉降量 Si、孔隙比 ei、压缩系数 av、压缩模量 Es，并按 av(1−2) 判定土的压缩性。

公式（GB/T 50123「14.3 快速压缩试验」+ E-8Intro + 权威 E-8.OUT 逐位反演）
--------------------------------------------------------------------------
    e0    = Gs(1+ω/100)/r − 1                       （初始孔隙比）
    校正前  ht_i = R_i·0.01 − L_i                   (mm)  （1h 变形减该级仪器变形）
    校正系数 K   = (h_n)_T / (h_n)_t
                  = (Rw·0.01 − L_n) / (R_n·0.01 − L_n)      （最后一级：稳定/1h）
    校正后  hT_i = K · ht_i                         (mm)
    h     = h0 − hT                                 (mm)
    Si    = hT/h0·1000                              (mm/m)
    ei    = e0 − (1+e0)·hT/h0
    av_i  = (e_{i−1} − e_i)/(P_i − P_{i−1})         (cm²/kg)
    Es_i  = (1+e0)·(P_i + P_{i−1})/(e_{i−1} − e_i)  (kg/cm²)
    （i=0 时 e_{−1}=e0、P_{−1}=0）

⚠ **Es 分母口径（由权威 OUT 反演确证）**：原著 Es 的分母用的是 **(P_i + P_{i−1})**
   （"+"而非"−"），与 av 的 (P_i − P_{i−1}) 不同。反证：11 行 Es 仅当取
   (P_i+P_{i−1}) 时逐行命中（例 P:1→2 级 Es=43.54；若取 (P_i−P_{i−1}) 则为 14.51，
   相对差 200%；若取 (1+e_{i−1})/av 则 13.57）。第 1 行 P:0→0.5 时两式等价，故不区分。
   同一族的 E-7 用的是标准 (1+e0)/av（权威 OUT 逐行命中）→ **两程序口径独立，未跨程序复用**。

未闭合点
--------
  1. av==0（卸荷零变形）时 Es 打印 −1.000（除零保护值），权威 OUT 无此情形（E-8 为单纯
     加荷），但 E-7 有；本内核两程序均按 −1.000 处置。
  2. 压缩性判定阈值（cm²/kg）：<0.01 低 / 0.01~0.05 中 / ≥0.05 高；权威算例 av(1−2)=0.0647
     判"高压缩性土"，与 ≥0.05 一致。低/中档阈值库中无算例标定，为**推断**。
  3. 权威 .OUT 首行「文件：…」为运行期路径回显，本内核不生成。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 快速法校正 Δh_i = (h_i)_t·(h_n)_T/(h_n)_t | `01_水工设计手册精读/卷1_基础理论.md` 4.6「土的压缩性（压缩试验/压缩指标）」；`02_水利教材精读/水能规划与岩土结构/土质学与土力学_第5版_精读笔记.md` 压缩性章 | **一致**（与 GB/T 50123 14.3.3-1 同式） | 采用 |
  | 三相换算 e0=Gs(1+ω)ρw/ρ−1 | 卷1 4.1「e=Gs(1+ω)ρw/ρ−1」；教材「孔隙比 e=Gsρw(1+w)/ρ−1」 | **一致** | 采用 |
  | 压缩系数 av=(e1−e2)/(p2−p1)、压缩模量 Es | 教材压缩性章 Es=(1+e1)/a | av **一致**；Es 分母 **不一致**（程序用 P_i+P_{i−1}） | 采用程序口径（已量化反证） |
  | av(1−2) 与压缩性分级 | 卷1 4.6、卷6「压缩系数」段 | 分级阈值库中未载 | 按旧规程登记（推断） |

输入数据（E-8Intro「（三）数据文件顺序」）
-------------------------------------------
    第 1 行： H0,Gs,r,ω,M         （H0 mm、Gs、r kg/cm³、ω %、M 加荷次数）
    其后 M 行： P,R,L             （P kg/cm²、R 量表读数 0.01mm、L 仪器变形 mm）
    末行 1 个数： Rw              （最后一级稳定时的量表读数 0.01mm）

输出（逐字复刻原著 .OUT 版式）
------------------------------
    表头 → 参数 → t/P/R/L/ht/hT/h/Si/ei/av/Es 表 → Pw/Rw → av(1-2)、压缩性
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-8"
TITLE = "快速压缩试验"
AUTHOR = "孟建平"
HEAD_NAME = "E-8"

RULE = "_" * 108
T_HOUR = 60          # 快速法各加压级固结时间 1h = 60 min

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
    if abs(x - round(x)) < 1e-9:
        return "%d" % int(round(x))
    return "%g" % x


# ============================================================
# 解析
# ============================================================

def parse(data):
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 7:
        raise ValueError("E-8 输入数据过少（%d 个数）" % len(nums))
    H0, Gs, r, w, M = nums[:5]
    M = int(round(M))
    rest = nums[5:]
    rows = []
    for i in range(M):
        rows.append((rest[3 * i], rest[3 * i + 1], rest[3 * i + 2]))
    Rw = rest[3 * M]
    return {"H0": H0, "Gs": Gs, "r": r, "w": w, "M": M, "rows": rows, "Rw": Rw}


# ============================================================
# 计算
# ============================================================

def compute(p):
    H0, Gs, r, w = p["H0"], p["Gs"], p["r"], p["w"]
    rows = p["rows"]
    Rw = p["Rw"]
    e0 = Gs * (1.0 + w / 100.0) / r - 1.0
    # 校正前（1h）总变形
    ht = [R * 0.01 - L for (P, R, L) in rows]
    # 校正系数（最后一级稳定/1h）
    hn_t = ht[-1]
    hn_T = Rw * 0.01 - rows[-1][2]
    K = hn_T / hn_t
    tab = []
    for i, (P, R, L) in enumerate(rows):
        hT = K * ht[i]
        tab.append({"t": T_HOUR, "P": P, "R": R, "L": L, "ht": ht[i], "hT": hT,
                    "h": H0 - hT, "Si": hT / H0 * 1000.0,
                    "ei": e0 - (1.0 + e0) * hT / H0, "av": None, "Es": None})
    for i in range(len(tab)):
        prev_e = tab[i - 1]["ei"] if i > 0 else e0
        prev_P = tab[i - 1]["P"] if i > 0 else 0.0
        de = prev_e - tab[i]["ei"]
        dP = tab[i]["P"] - prev_P
        av = de / dP if dP != 0 else 0.0
        if av == 0.0:
            av = 0.0
        denom_p = tab[i]["P"] + prev_P          # ★ 原著口径："+"（见 docstring）
        Es = -1.0 if av == 0.0 else (1.0 + e0) * denom_p / de
        tab[i]["av"] = av
        tab[i]["Es"] = Es
    av12 = None
    for i, row in enumerate(tab):
        if i >= 1 and abs(row["P"] - 2.0) < 1e-9:
            av12 = row["av"]
            break
    if av12 is None:
        av12 = tab[-1]["av"]
    if av12 < AV_LOW:
        soil = "低压缩性土"
    elif av12 < AV_HIGH:
        soil = "中压缩性土"
    else:
        soil = "高压缩性土"
    return {"输入": p, "e0": e0, "K": K, "hn_t": hn_t, "hn_T": hn_T,
            "表": tab, "Pw": rows[-1][0], "Rw": Rw, "av12": av12,
            "压缩性": soil}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, res):
    B = []
    B.append("")
    B.append(" ***********************************************************************")
    B.append(" ****                        快速压缩试验 E-8                       ****")
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
    B.append("    t        P        R          L        ht        hT        h"
             "          Si       ei          av        Es")
    B.append("  (min)   (kg/cm^2)(0.01mm)    (mm)      (mm)      (mm)      (mm)"
             "      (mm/m)    (-)       (cm^2/kg) (kg/cm^2)")
    # 第 1 行（P=0 初始态）：t/P/R/L=0，ht/hT=0，h=h0，Si=0，ei=e0
    z = res["表"][0]
    B.append(F(0, 5, 0) + F(0.0, 10, 1) + F(0.0, 10, 1) + F(0.0, 10, 3)
             + F(0.0, 10, 3) + F(0.0, 10, 3) + F(p["H0"], 11, 3)
             + F(0.0, 12, 3) + F(res["e0"], 6, 3) + " " * 14 + " " * 9)
    for row in res["表"]:
        B.append(F(row["t"], 5, 0) + F(row["P"], 10, 1) + F(row["R"], 10, 1)
                 + F(row["L"], 10, 3) + F(row["ht"], 10, 3) + F(row["hT"], 10, 3)
                 + F(row["h"], 11, 3) + F(row["Si"], 12, 3) + F(row["ei"], 6, 3)
                 + F(row["av"], 14, 4) + F(row["Es"], 9, 2))
    B.append("    Pw=%s  Rw= %s" % (F(res["Pw"], 5, 1), F(res["Rw"], 0, 1)))
    B.append(RULE)
    B.append("")
    B.append(" 压缩系数av(1-2)= %s(cm^2/kg)" % F(res["av12"], 0, 4))
    B.append(" %s" % res["压缩性"])
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
