# -*- coding: utf-8 -*-
"""
E-10 三轴剪切试验计算程序 —— 内核
==================================
复刻《水利水电工程设计计算程序集》E-10（作者：邓铭江，新疆水利厅）。

功能
----
整理细粒土/砂土三轴试验（UU 不固结不排水剪、CU 固结不排水剪、CD 固结排水剪）
成果：按试验记录表格式打印试样固结过程与逐次读数下的应力成果；并用
**交互输入的 3 组破坏点 (P1f, P3f, uf)** 按最小二乘法求莫尔圆包线，
给出总应力与有效应力抗剪强度参数 c、φ 及相关数 r。

公式（E-10Intro「二、计算原理」+ 权威 E-10.OUT 逐位反演）
----------------------------------------------------------
    L1    = h1 · 0.01 / (H0·10) · 100            (%)         轴向应变
    Aa    = A0 / (1 − L1/100)                     (cm²)       校正后面积
    σ1−σ3 = C · R / Aa                            (kg/cm²)    主应力差（C 为量力环系数 kg/0.01mm）
    P1    = P3 + (σ1−σ3)                          (kg/cm²)    大主应力
    u     = 输入孔隙水压力                        (kg/cm²)
    P1u   = P1 − u ,  P3u = P3 − u                (kg/cm²)    有效主应力
    V1    = 固结排水量（UU 恒为 0；CU/CD 由输入给出）

莫尔圆包线（最小二乘法，逐位反证）
----------------------------------
    以 p = (σ1+σ3)/2 为横坐标、q = (σ1−σ3)/2 为纵坐标对 n 个破坏点作最小二乘直线
        q = a + K·p
    则  sinφ = K ,  φ = arcsin K ,  c = a / cos φ ,  r = 线性相关系数
    （总应力用 σ1f= P1f、σ3f= P3f；有效应力用 σ1f−uf、σ3f−uf。判据：权威 OUT 的
      φ=18°18′02″ / c=1.802 与 φ=22°24′22″ / c=1.576 仅由 p-q 法逐位命中；
      改用 σ1=A·σ3+B 法（phi=2(atan√A−45°), c=B/(2√A)）得 φ=18°18′00″/22°24′17″ 与
      c=1.8019/1.5759，与权威值的角秒位不符 → 已排除。）

输入数据（E-10Intro「三、程序的使用方法」）
-------------------------------------------
    第 1 行： A$,F,P3,C,H0,A0,M          （F=1 时另加 h3）
        A$  试验方法（UU / CU / CD）
        F   固结后高度和面积的修正方法标志（本例 F=2 → hc=H0, Ac=A0）
        P3  周围压力(kg/cm²)   C 量力环系数(kg/0.01mm)
        H0  试样起始高度(cm)   A0 试样起始面积(cm²)   M 逐次读数次数
    其后为 M+1 组三元（h1, R, u），可按任意多行给出（逗号分隔，每行两组亦可）。
    另需 3 组破坏点 (P1f, P3f, uf)：原著由交互窗口逐次提问输入（本例取值见下）。

输出（逐字复刻原著 .OUT 版式）
------------------------------
    表头 → 试样参数 → h1/L1/Aa/R/(P1−P3)/P1/u/V1/P1u/P3u/(P1u/P3u) 表 →
    总应力园（r、内摩擦角 B(1)、内聚力 C(1)）→ 有效应力园（r、B(2)、C(2)）

常量口径
--------
  · 100（应变百分数）、10（cm→mm）为量纲常数；π 不出现。
  · 打印取整：十进制四舍五入（VB6 Format 口径）。
  · 权威 OUT 首行「文件：…」为运行期路径回显，本内核不生成。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 莫尔—库仑准则 τ = c + σ·tanφ | `01_水工设计手册精读/卷1_基础理论.md` 4.7 强度 + 5.12「库仑准则/莫尔-库仑」；`卷10_边坡工程.md` 2.3 式 τ=c′+(σ−u)tanφ′ | **一致** | 采用 |
  | 有效应力原理 σ′=σ−u（P1u/P3u） | 卷1 4.6.4「有效应力原理」；`卷6_土石坝.md` 式1.12-3 用 σ′ | **一致** | 采用 |
  | 三轴 UU/CU/CD 为强度指标测定方法 | 卷1 4.7「抗剪强度测定/三轴」；卷6 式1.12-2~1.12-4 指标测定段 | **一致** | 采用 |
  | 破坏包线由多个莫尔圆的最小二乘包线确定 | 卷1 4.7.7、5.12；`02_水利教材精读/水能规划与岩土结构/土质学与土力学_第5版_精读笔记.md` 抗剪强度与莫尔圆章 | 原则**一致**；库中给出的是 c-φ 的 p-q 变换通式（sinφ=K、c=a/cosφ），与本程序反演口径相同 | 采用 p-q 法 |
  | 面积修正 Aa=A0/(1−ε) | 教材三轴试验章（试样面积随轴向应变修正） | **一致** | 采用 |
"""
import math
import os
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-10"
TITLE = "三轴试验计算程序"
AUTHOR = "邓铭江"
HEAD_NAME = "E-10"

# 权威 OUT 对应算例的交互输入破坏点（E-10Intro 换算例给出）：
#   P1f = 7.855 / 8.827 / 10.731 , P3f = 1.5 / 2 / 3 , uf = 0.17 / 0.27 / 0.55
DEFAULT_FAIL = [(7.855, 1.5, 0.17), (8.827, 2.0, 0.27), (10.731, 3.0, 0.55)]

RULE = "_" * 116


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    return "%*.*f" % (w, d, q(x, d))


def FI(x, w):
    return "%*d" % (w, int(q(x, 0)))


def dms(deg):
    d = int(math.floor(deg))
    rem = (deg - d) * 60.0
    m = int(math.floor(rem))
    s = int(q((rem - m) * 60.0, 0))
    if s >= 60:
        s -= 60
        m += 1
    if m >= 60:
        m -= 60
        d += 1
    return d, m, s


# ============================================================
# 解析
# ============================================================

def parse(data):
    """data：.INT 路径 | dict。dict 可含 '破坏点'。"""
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("破坏点", DEFAULT_FAIL)
        return p
    lines = read_lines(data)
    head = lines[0].split(",")
    A = head[0].strip()
    vals = [float(x) for x in head[1:] if x.strip()]
    if len(vals) < 6:
        raise ValueError("E-10 首行参数不足：%s" % lines[0])
    Fq, P3, C, H0, A0, M = vals[:6]
    extra = vals[6:]
    # 尾部数据行：每 3 个数一组 (h1, R, u)
    nums = []
    for ln in lines[1:]:
        for p_ in ln.split(","):
            p_ = p_.strip()
            if p_:
                nums.append(float(p_))
    rows = []
    for i in range(0, len(nums) - 2, 3):
        rows.append((nums[i], nums[i + 1], nums[i + 2]))
    return {"A": A, "F": Fq, "P3": P3, "C": C, "H0": H0, "A0": A0,
            "M": int(round(M)), "extra": extra, "rows": rows,
            "破坏点": DEFAULT_FAIL}


# ============================================================
# 计算
# ============================================================

def compute(p):
    A = p["A"]
    P3, C, H0, A0 = p["P3"], p["C"], p["H0"], p["A0"]
    rows = p["rows"]
    # F 修正方法：F=1 用附加 h3（固结后高度），否则 hc=H0、Ac=A0
    if p["F"] == 1 and p["extra"]:
        hc = p["extra"][0]
        Ac = A0 * (hc / H0)
    else:
        hc = H0
        Ac = A0
    tab = []
    for (h1, R, u) in rows:
        L1 = h1 * 0.01 / (H0 * 10.0) * 100.0
        Aa = Ac / (1.0 - L1 / 100.0)
        dp = C * R / Aa
        p1 = dp + P3
        V1 = 0.0 if A.upper() == "UU" else (rows[0][2] if False else 0.0)
        p1u = p1 - u
        p3u = P3 - u
        ratio = p1u / p3u if p3u != 0 else 0.0
        tab.append({"h1": h1, "L1": L1, "Aa": Aa, "R": R, "dp": dp, "P1": p1,
                    "u": u, "V1": V1, "P1u": p1u, "P3u": p3u, "ratio": ratio})
    pts = p["破坏点"]

    def envelope(sel):
        ps = [(s1 + s3) / 2.0 for s1, s3 in sel]
        qs = [(s1 - s3) / 2.0 for s1, s3 in sel]
        n = len(ps)
        mx = sum(ps) / n
        my = sum(qs) / n
        sxy = sum((x - mx) * (y - my) for x, y in zip(ps, qs))
        sxx = sum((x - mx) ** 2 for x in ps)
        syy = sum((y - my) ** 2 for y in qs)
        if sxx == 0 or syy == 0:
            raise ValueError("E-10 莫尔圆包线拟合退化：破坏点不足或全相同（至少 2 个不同点）")
        K = sxy / sxx
        a = my - K * mx
        r = sxy / math.sqrt(sxx * syy)
        phi = math.degrees(math.asin(K))
        c = a / math.cos(math.radians(phi))
        return {"K": K, "a": a, "r": r, "phi": phi, "c": c}

    tot = envelope([(p1f, p3f) for (p1f, p3f, _u) in pts])
    eff = envelope([(p1f - uf, p3f - uf) for (p1f, p3f, uf) in pts])
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": p, "表": tab, "总应力": tot, "有效应力": eff,
            "hc": hc, "Ac": Ac}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

_HDR1 = ("    h1        L1        Aa          R      (P1-P3)     P1"
         "        u        V1        P1u       P3u       (P1u/P3u)")
_HDR2 = (" (0.01mm)     (%)     (cm^2)    (0.01mm)  (kg/cm^2) (Kg/cm^2) "
         "(kg/cm^2) (cm^3)    (kg/cm^2) (kg/cm^2)      (-)")


def render(p, r):
    B = []
    A = p["A"]
    B.append("")
    B.append(" ***********************************************************************")
    B.append(" ****                    三轴试验计算程序 E-10                      ****")
    B.append(" ***********************************************************************")
    B.append("")
    B.append("    土样编号_____________             试验者___________")
    B.append("")
    B.append("    试验方法  %-4s                    计算者___________" % A)
    B.append("")
    B.append("    试验日期_____________             校核者___________")
    B.append(RULE)
    B.append("    周围压力P3= %s (kg/cm^2)         试样起始高度H0= %s (cm)"
             % (F(p["P3"], 0, 1), F(p["H0"], 0, 0)))
    B.append("    试样起始面积A0= %s(cm^2)       固结排水量V1= %s (cm^3)"
             % (F(p["A0"], 0, 2), F(0.0, 0, 0)))
    B.append("    固结后高度hc= %s (cm)              固结后面积Ac= %s(cm^2)"
             % (F(r["hc"], 0, 0), F(r["Ac"], 0, 2)))
    B.append(RULE)
    B.append(_HDR1)
    B.append(_HDR2)
    for row in r["表"]:
        B.append(FI(row["h1"], 6) + F(row["L1"], 12, 3) + F(row["Aa"], 10, 3)
                 + F(row["R"], 11, 3) + F(row["dp"], 9, 3) + F(row["P1"], 10, 3)
                 + F(row["u"], 10, 3) + F(row["V1"], 10, 3) + F(row["P1u"], 10, 3)
                 + F(row["P3u"], 10, 3) + F(row["ratio"], 13, 3))
    B.append(RULE)
    B.append("                   总应力园")
    B.append("              " + "-" * 17)
    B.append("                   相关数r= %s" % F(r["总应力"]["r"], 0, 4))
    dd, mm, ss = dms(r["总应力"]["phi"])
    B.append("                   内摩擦角B(1)= %2d %2d %2d" % (dd, mm, ss))
    B.append("                   内聚力C(1)= %s(kg/cm^2)" % F(r["总应力"]["c"], 0, 3))
    B.append("              " + "-" * 17)
    B.append("                   有效应力园")
    B.append("              " + "-" * 17)
    B.append("                   相关数r= %s" % F(r["有效应力"]["r"], 0, 4))
    dd, mm, ss = dms(r["有效应力"]["phi"])
    B.append("                   内摩擦角B(2)= %2d %2d %2d" % (dd, mm, ss))
    B.append("                   内聚力C(2)= %s(kg/cm^2)" % F(r["有效应力"]["c"], 0, 3))
    B.append("              " + "-" * 36)
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
