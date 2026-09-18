# -*- coding: utf-8 -*-
"""
D-24 水电站引水渠非恒定流计算程序 —— 内核
============================================
复刻《水利水电工程设计计算程序集》D-24（作者：卢礼标）。
原著说明书明示"主要参考切尔陀乌索夫著《水力学专门教程》提出的明渠非恒定流计算
方法编写"，并给出两套基本假定：
  落波（JQ=1，求最低水位）：① 波峰保持垂直；② 波通过后水面曲线成为直线；
                            ③ 逆落波通过后各渠段水面坡降为常数；
  涌波（JQ=2，求最高水位）：① 波峰保持垂直；② 波通过后水面线成为水平直线；
                            ③ 波进行时无阻力存在。
本程序为梯形/矩形断面引水渠。

★已由权威 D-24-1.OUT 反演验证的部分（算例1：JQ=1,K=3,BD=3,M=1.5,N=.0225,
Q0=10,ΔQ=10，断面 (S,H0,Z0) 见下）
--------------------------------------------------------------------
对每个断面，波前两侧 Q1=Q0、Q2=Q0+ΔQ，断面面积 ω(h)=(BD+M·h)h，
水面宽 B(h)=BD+2M·h，静水压力矩 I(h)=M·h³/3+BD·h²/2。
  波前速度      c = ΔQ/(ω2−ω1)                （连续）
  落波高度      Δh = ΔQ/(B̄·|c|)，B̄=(B1+B2)/2  （=Δω/B̄，与权威 DH0 恒等）
  权威 DH0 = .202812993379937 → 对应 h2 = 3.166187007（Δh=0.202812993）
★波前高度方程（本内核实现口径）
  采用**比力（静矩）式**：f(Q,h)=Q²/(g·ω)+I(h)，
    (f2−f1)·(ω2−ω1) = (ΔQ)²/g
  本内核据此解 h2 → Δh，各断面：0.202217 / 0.220558 / 0.241906 / 0.266138
  （渠末 Δh=0.202217 vs 权威 0.202813，**相对偏差 2.94E-03**）
★未闭合（如实量化）
  在权威 h2 上，各候选波前方程残差：比力式 −2.540E-02、随波坐标动量 +2.540E-02、
  随波坐标能量 +5.531E-04（最接近但非零，相对 1.3E-04）。说明原著波前方程含本
  内核尚未还原的项（疑与初始水面坡降 J=4.375E-05 / 阻力有关）。
  因此 DH0 记为 DECL（2.94E-03 相对偏差）。
  DHT（逆落波至渠首瞬间渠末降）、DHmax（渠末最大降）、TM（总历时）、
  各断面最低水位 ZM —— **未反演出**：已确认 ZM 各断面为**直线**（99.517/99.689/
  99.860/100.032，逐段差 0.172/0.171/0.172），即假定②的直接体现；其直线斜率
  ≈2.1437E-04（由 ZM(3) 反推）与 DHmax 的关系式尚未锁定；体积平衡
  ΔV=∫[ω(h_old)−ω(h_new)]ds 与 ΔQ·TM=11310 m³ 不符（8179 vs 11310，差 28%），
  说明渠首入流亦随时间变化（反射波），需补切尔陀乌索夫原书完整演算流程。
  本内核按已获方程实现 DH0 与波前历时 T，其余量输出——并以 None 占位，
  在 verify 中记为 DECL/未复现。
"""
import math

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-24"
TITLE = "水电站引水渠非恒定流计算"
AUTHOR = "卢礼标"
HEAD_NAME = "D-24"

G = 9.81


def area(h, BD, MST):
    return (BD + MST * h) * h


def width(h, BD, MST):
    return BD + 2.0 * MST * h


def mom(h, BD, MST):
    return MST * h ** 3 / 3.0 + BD * h * h / 2.0


# ============================================================
# 解析
# ============================================================

def parse(data):
    """解析输入。data：dict | .INT 文件路径。
    .INT 结构：首行「工程名,JQ,K,M,BD,N,Q0,DQ」，其后每行「S,H0,Z0」（共 K+1 行，
    允许一行内并列多组）。"""
    if isinstance(data, dict):
        p = dict(data)
        for kk in ("JQ", "K"):
            if kk in p:
                p[kk] = int(round(p[kk]))
        return p
    lines = [ln.strip() for ln in read_lines(data) if ln.strip()]
    head = lines[0].split(",")
    name = head[0].strip()
    nums = [float(x) for x in head[1:] if x.strip()]
    JQ, K, MST, BD, N, Q0, DQ = nums[:7]
    sec = []
    for ln in lines[1:]:
        ps = [float(x) for x in ln.split(",") if x.strip()]
        sec.append(ps)
    flat = []
    for ps in sec:
        flat.extend(ps)
    rows = [flat[i:i + 3] for i in range(0, len(flat), 3)]
    rows = [r for r in rows if len(r) == 3]
    return {"name": name, "JQ": int(round(JQ)), "K": int(round(K)), "M": MST,
            "BD": BD, "N": N, "Q0": Q0, "DQ": DQ, "S": [r[0] for r in rows],
            "H0": [r[1] for r in rows], "Z0": [r[2] for r in rows]}


# ============================================================
# 计算
# ============================================================

def _solve_bore(h1, Q1, Q2, BD, MST):
    """解波前方程 (f2−f1)·(ω2−ω1) = (ΔQ)²/g，返回 (h2, c)。"""
    def eq(h2):
        w1, w2 = area(h1, BD, MST), area(h2, BD, MST)
        f1 = Q1 * Q1 / (G * w1) + mom(h1, BD, MST)
        f2 = Q2 * Q2 / (G * w2) + mom(h2, BD, MST)
        return (f2 - f1) * (w2 - w1) - (Q2 - Q1) ** 2 / G
    lo, hi = h1 * 0.2, h1 * (1 - 1e-12)
    flo, fhi = eq(lo), eq(hi)
    if flo * fhi > 0:
        return None, None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = eq(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    h2 = 0.5 * (lo + hi)
    c = (Q2 - Q1) / (area(h2, BD, MST) - area(h1, BD, MST))
    return h2, c


def compute(p):
    BD, MST = p["BD"], p["M"]
    Q0, DQ = p["Q0"], p["DQ"]
    JQ = p["JQ"]
    K = p["K"]
    S, H0 = p["S"], p["H0"]
    Q1 = Q0
    Q2 = Q0 + DQ if JQ == 1 else Q0 - DQ
    dh = []
    cc = []
    for h1 in H0:
        h2, c = _solve_bore(h1, Q1, Q2, BD, MST)
        dh.append(None if h2 is None else h1 - h2)
        cc.append(c)
    # 波前传播历时（逐段用两端波速的算术平均）
    Tmarch = 0.0
    for i in range(K):
        ca, cb = abs(cc[i]), abs(cc[i + 1])
        if ca > 0 and cb > 0:
            Tmarch += (S[i + 1] - S[i]) / (0.5 * (ca + cb))
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR, "输入": p,
            "波前高度": dh, "波速": cc, "传播历时": Tmarch,
            "DH0": dh[0] if dh else None}


# ============================================================
# 输出（原著版式）
# ============================================================

def render(p, r):
    BD, MST = p["BD"], p["M"]
    K = p["K"]
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" *****               水电站引水渠非恒定流计算  D-24                *****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("            工程名:%s" % p["name"])
    A.append(" 一.原始数据")
    case = ("推求水电站突然增加负荷时引水渠的最低水位" if p["JQ"] == 1
            else "推求水电站突然丢弃负荷时引水渠的最高水位")
    A.append("    计算情况 JQ=%d (%s)" % (p["JQ"], case))
    A.append("    渠道分段数   K= %d " % K)
    A.append("    渠道底宽    BD= %g (m)" % BD)
    A.append("    渠道边坡系数 M= %g " % MST)
    A.append("    渠道糙率     N= %s " % ("%g" % p["N"]).lstrip("0"))
    A.append("    负荷变化前水电站流量 Q0= %g (m3/s)" % p["Q0"])
    A.append("    负荷变化引起水电站流量的改变量(绝对值) ΔQ=　 %g (m3/s)" % p["DQ"])
    A.append("          引水渠各断面的断面水深")
    A.append("          与水面高程(负荷变化前)")
    A.append("")
    A.append("    ----------------------------------------")
    A.append("      断面号  断面距渠末  断面水深  水面高程")
    A.append("        I      距离S(m)      H0(m)     Z0(m)")
    A.append("    ----------------------------------------")
    for i in range(K + 1):
        A.append("%9d%10g%12g%11.3f" % (i, p["S"][i], p["H0"][i], p["Z0"][i]))
    A.append("    ----------------------------------------")
    A.append("")
    A.append(" 二.计算中采用的基本假定")
    if p["JQ"] == 1:
        A.append("    1.在波进行的全过程中,波峰保持垂直;")
        A.append("    2.渠道中的水面曲线,在波通过后成为直线;")
        A.append("    3.在逆落波通过后,各渠段的水面坡降为常数.")
    else:
        A.append("    1.在波进行的全过程中,波峰保持垂直;")
        A.append("    2.渠道中的水面曲线,在波通过后成为水平直线;")
        A.append("    3.在波进行时无阻力存在.")
    A.append("")
    A.append(" 三.计算结果")
    dn = r["DH0"]
    A.append("   渠末第一次落波高度 DH0= %s (m)"
             % ("%s" % (dn if dn is not None else "")).lstrip("0"))
    for tag, key in (("逆落波传至渠首瞬间渠末的水位下降值 DHT", "DHT"),
                     ("渠末最大水位降低值 DHmax", "DHmax")):
        v = r.get(key)
        A.append("   %s= %s (m)" % (tag, "" if v is None else v))
    A.append("   逆落波由渠末传至渠首的历时 T= %.0f (s)" % r["传播历时"])
    v = r.get("TM")
    A.append("   逆落波与反射波在渠道中传播的总历时 TM= %s (s)" % ("" if v is None else v))
    A.append("")
    A.append("           渠道各断面的最低水位")
    A.append("         ----------------------")
    A.append("           断面号I  最低水位ZM(m)")
    A.append("         ----------------------")
    for i in range(K + 1):
        A.append("%15d%14s" % (i, ""))
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
