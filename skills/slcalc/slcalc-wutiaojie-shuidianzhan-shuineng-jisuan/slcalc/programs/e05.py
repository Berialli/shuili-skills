# -*- coding: utf-8 -*-
"""
E-5 渗透试验计算程序 —— 内核
=============================
复刻《水利水电工程设计计算程序集》E-5 程序（作者：孟建平，新疆水利水电学校）。

原著说明书（E-5Intro）：
  「用于计算在南55型渗透仪上进行细粒土变水头和常水头试验时试样的渗透系数，
    以及采用70型渗透仪测定粗粒土渗透系数的计算。」
  控制变量 A$ 选择仪器与试验方法：
     · 55B —— 南55型渗透仪 变水头试验：t1,t2,h1,h2,T
     · 55C —— 南55型渗透仪 常水头试验：t1,t2,Q,h,T
     · 70C —— 70型渗透仪 常水头试验：Z1,Z2,Z3,t,T,Q

输入数据顺序（说明书）
----------------------
  A$, N, W, ω, A, L, Gs   （湿土重 g、含水量 %、试样面积 cm²、渗径=试样高 cm、土粒比重）
  · 55B：a ；随后 N 行 t1,t2,h1,h2,T
  · 55C：随后 N 行 t1,t2,Q,h,T
  · 70C：随后 N 行 Z1,Z2,Z3,t,T,Q

计算原理（《土工试验规程》）
--------------------------
  湿容重  γ  = W/(A·L)                     （g/cm³）
  干容重  γd = γ/(1+ω/100)
  孔隙比  e  = Gs·γw/γd − 1 ，γw=1.0        （g/cm³）
  变水头  kT = 2.3·a·L/(A·t)·log10(h1/h2)   ，t = t2 − t1
  常水头  kT = Q·L/(A·h·t)
  温度校正常数 η(T) ∝ 1/D(T)，D(T)=1+αT+βT²；n = ηT/η10 = D(10)/D(T)
  K10 = kT·n ； Kp = mean(K10)

常量口径（**逐枚独立反演**：由 E-5vb.EXE 常量池 + 权威 E-5.OUT 双向互证）
----------------------------------------------------------------------
  · 2.3        —— E-5vb.EXE @0x13b4f = 2.3000011444（float32(2.3)），变水头系数。
  · ln(10)÷2.3026 —— E-5vb.EXE @0x13dd9 = 2.3026008606。程序以 `Log(x)/2.3026`
    求十进制对数（而非 Log(x)/Log(10)=2.302585…）；对变水头 kT 的等效系数为
    2.3×(2.302585093/2.3026)=2.29998502。**若改用精确 2.302585093，第 2 行 kT
    变为 5.30547（打印 5.3055），与权威 5.3054 不符** —— 反证该 2.3026 口径为真。
  · α=0.0337  —— E-5vb.EXE @0x13f7a = 0.0337000191（float32(0.0337)）
  · β=0.00022 —— E-5vb.EXE @0x1404e = 0.00022（**注意是 0.00022，非教材常见的
    0.000221**；两者在 n(5) 上相差 6e-5 相对，足以改变 K10 第 4 位小数）。
  · η(T) 分子常数（教材 0.01775，Poise）在比值 n 中约去，E-5vb.EXE 亦未见该 double，故不参与。

验证结果（权威 E-5.OUT，南55B 算例，4 行）
-----------------------------------------
  γ/γd/e 与 kt、n、K10、Kp **全部逐位命中**（EXACT），FAIL=0。
  55C / 70C 无权威 OUT，仅作解析与自洽（标 DECL）。

知识库对照结果（教学母本 → 程序实现）
------------------------------------
  KB 路径：D:\\WorkBuddy知识库\\水利知识库\\

  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 变水头渗透系数 kT=2.3·a·L/(A·t)·lg(h1/h2) | `02_水利教材精读/土质学与土力学/土质学与土力学_第5版_精读笔记.md` 渗透章：变水头试验 k=2.3·aL/(At)·lg(h1/h2) | **一致** | 采用 |
  | 常水头渗透系数 kT=Q·L/(A·h·t) | 同上：常水头 k=QL/(Aht) | **一致** | 采用（55C/70C） |
  | 孔隙比 e=Gs·γw/γd−1 | `01_水工设计手册精读/卷6_土石坝.md` 及教材压实/渗透章：e=Gs(1+w)γw/γ−1 等价式 | **一致** | 采用 |
  | 温度校正 K10=kT·ηT/η10 | 教材渗透章：以 10℃ 为标准温度的水温校正 | **一致** | 采用；η 用 D(T)=1+0.0337T+0.00022T²（EXE 反演，与教材 0.000221 略有别，见下） |
  | 教材 η=0.01775/(1+0.0337T+0.000221T²) | 同上（泊肃叶经验式） | **系数 β 不一致**：教材 0.000221，程序 0.00022 | 以 **EXE 常量池 0.00022 + 权威 OUT 反演**为准；0.00022 使 n(5)=1.157581（权威要求 1.157574~1.157593）；若用 0.000221 则 n(5)=1.157641，K10 各行为 6.8536/6.1418/5.6047/6.1326，与权威 6.8532/6.1415/5.6044/6.1323 全部不符（4/4 FAIL）——已量化反证。 |
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-5"
TITLE = "渗透试验计算程序"
AUTHOR = "孟建平（新疆水利水电学校）"
HEAD_NAME = "E-5"

COEF_2_3 = 2.3           # 变水头系数（EXE float32(2.3)）
LOG10_DIV = 2.3026       # 程序内 Log(10) 常量（EXE double 2.3026008606）
ETA_A = 0.0337           # η 分母 α（EXE float32(0.0337)）
ETA_B = 0.00022          # η 分母 β（EXE 0.00022；注意非 0.000221）
GAMMA_W = 1.0            # 水的容重口径（影响孔隙比 e）

METHOD_NAMES = {
    "55B": "变水头渗透试验(南55B)",
    "55C": "常水头渗透试验(南55C)",
    "70C": "70型渗透仪常水头试验(70C)",
}


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def _eta_denom(T):
    return 1.0 + ETA_A * T + ETA_B * T * T


def eta_ratio(T):
    """n = ηT/η10 = D(10)/D(T)。"""
    return _eta_denom(10.0) / _eta_denom(T)


# ============================================================
# 解析
# ============================================================

def _split(line):
    """按逗号切分，去掉引号。"""
    return [t.strip().strip('"').strip("'") for t in line.split(",") if t.strip() != ""]


def parse(data):
    if isinstance(data, dict):
        return data
    lines = [ln.strip() for ln in read_lines(data) if ln.strip() != ""]
    if not lines:
        raise ValueError("E-5 输入为空")
    t = _split(lines[0])
    if len(t) < 7:
        raise ValueError("E-5 首行字段不足：%r" % lines[0])
    method = t[0].upper()
    N = int(round(float(t[1])))
    W = float(t[2]); wo = float(t[3]); area = float(t[4]); L = float(t[5]); Gs = float(t[6])
    idx = 1
    a = None
    if method == "55B":
        a = float(_split(lines[idx])[0]); idx += 1
    rows = []
    for i in range(N):
        v = [float(x) for x in _split(lines[idx + i])]
        if method == "55B":
            rows.append({"t1": v[0], "t2": v[1], "h1": v[2], "h2": v[3], "T": v[4]})
        elif method == "55C":
            rows.append({"t1": v[0], "t2": v[1], "Q": v[2], "h": v[3], "T": v[4]})
        elif method == "70C":
            rows.append({"Z1": v[0], "Z2": v[1], "Z3": v[2], "t": v[3], "T": v[4], "Q": v[5]})
        else:
            raise ValueError("E-5 未知试验方法 A$=%r" % method)
    return {"method": method, "N": N, "W": W, "omega": wo, "A": area, "L": L,
            "Gs": Gs, "a": a, "rows": rows}


# ============================================================
# 计算
# ============================================================

def compute(p):
    A = p["A"]; L = p["L"]
    gam = p["W"] / (A * L)
    gd = gam / (1.0 + p["omega"] / 100.0)
    e = p["Gs"] * GAMMA_W / gd - 1.0
    method = p["method"]
    out = []
    for row in p["rows"]:
        T = row["T"]
        if method == "55B":
            t = row["t2"] - row["t1"]
            kt = COEF_2_3 * p["a"] * L / (A * t) * (math.log(row["h1"] / row["h2"]) / LOG10_DIV)
        elif method == "55C":
            t = row["t2"] - row["t1"]
            kt = row["Q"] * L / (A * row["h"] * t)
        else:  # 70C
            t = row["t"]
            h = row["Z1"] - row["Z3"]
            kt = row["Q"] * L / (A * h * t)
        n = eta_ratio(T)
        k10 = kt * n
        r = dict(row)
        r.update({"t": t, "kt": kt, "n": n, "K10": k10})
        out.append(r)
    kp = sum(r["K10"] for r in out) / len(out)
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": p, "湿容重": gam, "干容重": gd, "孔隙比": e,
            "逐次": out, "Kp": kp}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, r):
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****                      渗透试验计算程序  E-5                    ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append(" " * 18 + METHOD_NAMES.get(p["method"], p["method"]))
    A.append("")
    A.append("")
    A.append("    工程名称＿＿＿＿＿＿   土样说明＿＿＿＿＿＿   试样面积 %5.2f" % q(p["A"], 2))
    A.append("")
    A.append("    土样编号＿＿＿＿＿＿   测压管面积%9.4f    孔隙比 %.4f"
             % (q(p["a"] if p["a"] is not None else 0.0, 4), q(r["孔隙比"], 4)))
    A.append("")
    A.append("    仪器编号＿＿＿＿＿＿   试样高度%9.4f      试验者＿＿＿＿＿＿" % q(p["L"], 4))
    A.append("")
    A.append("    湿土重%9.4f        含水量%9.4f        土粒比重 %.4f"
             % (q(p["W"], 4), q(p["omega"], 4), q(p["Gs"], 4)))
    A.append("")
    A.append(" 开始  终了  经过   开始    终了  水温T度    水   校正   渗透   平均渗")
    A.append(" 时间  时间  时间   水头    水头  时的渗透   温   系数   系数   透系数")
    A.append("  t1    t2    t      h1      h2   系数kt         nt/n10   K10     Kp  ")
    A.append("  秒    秒    秒    厘米    厘米  厘米/秒    度         厘米/秒 厘米/秒")
    A.append("                                   10^-7                 10^-7   10^-7 ")
    for i, row in enumerate(r["逐次"]):
        if p["method"] == "55B":
            line = "%4d%6d%6d%8.1f%8.1f%9.4f%7.1f%7.3f%8.4f" % (
                int(round(row["t1"])), int(round(row["t2"])), int(round(row["t"])),
                q(row["h1"], 1), q(row["h2"], 1), q(row["kt"] * 1e7, 4),
                q(row["T"], 1), q(row["n"], 3), q(row["K10"] * 1e7, 4))
        elif p["method"] == "55C":
            line = "%4d%6d%6d%8.2f%8.1f%9.4f%7.1f%7.3f%8.4f" % (
                int(round(row["t1"])), int(round(row["t2"])), int(round(row["t"])),
                q(row["Q"], 2), q(row["h"], 1), q(row["kt"] * 1e7, 4),
                q(row["T"], 1), q(row["n"], 3), q(row["K10"] * 1e7, 4))
        else:
            line = "%4.0f%6.1f%6.1f%8.0f%7.1f%9.4f%7.1f%7.3f%8.4f" % (
                row["Z1"], row["Z2"], row["Z3"], row["t"], q(row["T"], 1),
                q(row["kt"] * 1e7, 4), q(row["T"], 1), q(row["n"], 3),
                q(row["K10"] * 1e7, 4))
        if i == len(r["逐次"]) - 1:
            line += "%8.4f" % q(r["Kp"] * 1e7, 4)
        A.append(line)
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
