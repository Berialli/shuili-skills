# -*- coding: utf-8 -*-
"""
F-1X 高斯投影计算程序 —— 内核
=================================
复刻《水利水电工程设计计算程序集》F-1X 程序（作者：谢希哲，新疆兵团勘测设计院）。
原著说明书（RTF\\F-1XIntro.RTF）「一. 简介」第 4 条：

    计算方法：
      正算: 先计算相应纬度的子午线弧长,再按展开式计算Ｘ、Ｙ及γ；
      反算: 先按正算中的子午线弧长进行迭代,求得底点纬度 Bf,再按展开式计算B、L及γ；
      换带: 先按 L1 进行大地反算,再按 L2 进行大地正算；

    第 3 条：可适用于 1954 年北京坐标系及 1980 西安坐标系；
    第 5 条：Y 坐标可加亦可不加带号,程序能自动识别；
    「三. 操作注意事项」第 6 条：投影带宽只能选择 3 或 6。

公式（标准高斯—克吕格投影，与《水利工程测量》第5版 §10.2「高斯投影与分带」一致）
-------------------------------------------------------------------------------
  子午线弧长（e 幂级数，取至 e^8）：
      e² = 2f − f²
      X = m₀·[ C₀·B + C₂·sin2B/2 + C₄·sin4B/4 + C₆·sin6B/6 + C₈·sin8B/8 ]
      m₀ = a(1−e²)
      C₀ = 1 + 3/4·e² + 45/64·e⁴ + 175/256·e⁶ + 11025/16384·e⁸
      C₂ = −(3/4·e² + 15/16·e⁴ + 525/512·e⁶ + 2205/2048·e⁸)
      C₄ = 15/64·e⁴ + 105/256·e⁶ + 2205/4096·e⁸
      C₆ = −(35/512·e⁶ + 315/2048·e⁸)
      C₈ = 315/16384·e⁸
  正算（底点纬度展开展开至 l⁶）：
      辅助量 N = a/√(1−e²sin²B)、t = tanB、η² = e′²cos²B、e′² = e²/(1−e²)
      x = X + N/2·sinB·cosB·l²
            + N/24·sinB·cos³B·(5−t²+9η²+4η⁴)·l⁴
            + N/720·sinB·cos⁵B·(61−58t²+t⁴+270η²−330η²t²)·l⁶
      y = N·cosB·l + N/6·cos³B·(1−t²+η²)·l³
            + N/120·cos⁵B·(5−18t²+t⁴+14η²−58η²t²)·l⁵
      γ = l·sinB + l³/3·sinB·cos²B·(1+3η²+2η⁴) + l⁵/15·sinB·cos⁴B·(2−t²)
  反算：以子午线弧长对 x 迭代求底点纬度 Bf，再按
      B = Bf − t/(2N²)(1+η²)y² + t/(24N⁴)(5+3t²+6η²−6η²t²−3η⁴−9η⁴t²)y⁴
            − t/(720N⁶)(61+90t²+45t⁴)y⁶
      l = y/(N·cosB) − (1+2t²+η²)y³/(6N³cosB)
            + (5+28t²+24t⁴+6η²+8η²t²)y⁵/(120N⁵cosB)
      L = L₀ + l，γ 同正算式。

输入数据（.INT，逗号/空白分隔的扁平数值流 + 行结构）
---------------------------------------------------
  第 1 行:  测区名称,计算者,日期,类型,坐标系,K,标志[,带宽,Lo 或 带宽1,L1,带宽2,L2]
      类型 ∈ {+ 大地正算, − 大地反算, T 换带计算}
      坐标系 ∈ {1:1954年北京坐标系（克拉索夫斯基）, 2:1980年西安坐标系（1975国际椭球）}
      K        = Y 坐标加常数（km），一般 500
      标志     = 1（本程序三个算例均为 1；输出 Y 加带号）
      带宽     = 6 或 3（投影带宽；反算时为 0 占位）
      Lo       = 中央子午线（度）；L1/L2 = 换带前/后中央子午线
      ★ 反算文件 F-1X-2.INT 的 Lo 另起一行（`…,1,0\n81\n"WINSU"\n…`），
        本内核按「数值流续读」处理，与 VB6 Input# 的空白/换行同作分隔符一致。
  其后每点两行：点名行、坐标行（正算为 B,L；反算/换带为 X,Y），
      角度坐标以 DD.MMSSss 压缩格式书写（44.351697 = 44°35'16.97″）。

输出（逐字复刻原著 .OUT 版式，GBK；表格用制表符外框）
-----------------------------------------------------
  第 1 行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\F-1X-1.out」为运行期路径回显，不生成。
  自第 2 行起逐行复刻：空行 → 4 行星号题头 → 空行 → 测区行 → 表格（6 行框线）
  → 计算者行。

对拍基准与闭合状态
------------------
  · F-1X-1.OUT（G 盘 RTF\\算例计算结果文件\\F\\，917 B，12 行）为**唯一权威 OUT**。
    本内核 1954 北京坐标系正算结果与权威**逐位相同**：
    X=4939431.543、Y=15451773.622、γ=−0°25′34.8078″ —— EXACT。
  · F-1X-2（大地反算）、F-1X-3（换带计算）**G 盘无权威 OUT，说明书亦未刊印算例结果**
    （F-1XIntro.RTF 全文已逐字提取核对，仅有序言/操作方法/注意事项三节，无结果表）。
    本内核以「正算模型往返闭合」自校：反算 → 正算回代，残差见 verify（≤1e-7 m 量级），
    计 DECL（已量化）。
  · 换带计算的 γ 取换带后中央子午线 L2 处的收敛角（原著未明示，计 DECL）。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  母本：`02_水利教材精读\\水利工程施工与管理\\水利工程测量_第5版_精读笔记.md`
        §2.7 地形图 / §3.3 公式表 10-1、10-2 / §4.5 导线测量内业
  逐条对照见文件末「知识库对照结果」节。
"""
import math
import os

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "F-1X"
TITLE = "高斯投影计算程序"
AUTHOR = "谢希哲"
HEAD_NAME = "F-1X"

# 椭球：(名称, a, 1/f)
ELLIPSOIDS = {
    1: ("1954年北京坐标系", 6378245.0, 1.0 / 298.3),      # 克拉索夫斯基
    2: ("1980年西安坐标系", 6378140.0, 1.0 / 298.257),    # 1975 国际椭球 (IAG75)
}

MODE_CN = {"+": "大地正算", "-": "大地反算", "T": "换带计算"}


# ============================================================
# 角度格式
# ============================================================

def dmss_to_deg(v):
    """DD.MMSSssss 压缩角 → 十进制度（44.351697 → 44.588047222…）。"""
    sgn = -1.0 if v < 0 else 1.0
    s = "%.6f" % abs(float(v))
    ip, fp = s.split(".")
    fp = (fp + "000000")[:6]
    d = int(ip)
    m = int(fp[:2])
    sec = float(fp[2:4] + "." + fp[4:6]) if len(fp) >= 6 else 0.0
    return sgn * (d + m / 60.0 + sec / 3600.0)


def _round4(x):
    """VB6 打印取整口径：十进制四舍五入（半值进位）。"""
    return math.floor(x * 10000.0 + 0.5) / 10000.0


def dms_fields(deg):
    """十进制度 → (度, 分, 秒, 符号)，秒保留 4 位小数并进位。"""
    sgn = "-" if deg < 0 else "+"
    x = abs(deg)
    d = int(x)
    rem = (x - d) * 60.0
    m = int(rem)
    s = _round4((rem - m) * 60.0)
    if s >= 60.0:
        s -= 60.0
        m += 1
    if m >= 60:
        m -= 60
        d += 1
    return d, m, s, sgn


def fmt_dms(deg):
    """'%3d %02d %07.4f' → ' 44 35 16.9700'（列宽 14）。"""
    d, m, s, _ = dms_fields(deg)
    return "%3d %02d %07.4f" % (d, m, s)


def fmt_gamma(deg):
    """'±DD MM SS.ssss' → '- 0 25 34.8078'（列宽 14）。"""
    d, m, s, sgn = dms_fields(deg)
    return "%s%2d %02d %07.4f" % (sgn, d, m, s)


# ============================================================
# 投影正/反算
# ============================================================

def meridian_arc(B, a, f):
    """子午线弧长 X（e 幂级数至 e^8）。"""
    e2 = 2.0 * f - f * f
    C0 = 1.0 + 0.75 * e2 + 45.0 / 64.0 * e2 ** 2 + 175.0 / 256.0 * e2 ** 3 \
        + 11025.0 / 16384.0 * e2 ** 4
    C2 = -(0.75 * e2 + 15.0 / 16.0 * e2 ** 2 + 525.0 / 512.0 * e2 ** 3
           + 2205.0 / 2048.0 * e2 ** 4)
    C4 = 15.0 / 64.0 * e2 ** 2 + 105.0 / 256.0 * e2 ** 3 + 2205.0 / 4096.0 * e2 ** 4
    C6 = -(35.0 / 512.0 * e2 ** 3 + 315.0 / 2048.0 * e2 ** 4)
    C8 = 315.0 / 16384.0 * e2 ** 4
    m0 = a * (1.0 - e2)
    return m0 * (C0 * B + C2 * math.sin(2 * B) / 2.0 + C4 * math.sin(4 * B) / 4.0
                 + C6 * math.sin(6 * B) / 6.0 + C8 * math.sin(8 * B) / 8.0)


def _gamma_eta(l, Brad, eta2):
    """子午线收敛角 γ（弧度）。l、B 均为弧度；用**真实纬度 B**。"""
    sinB, cosB, t = math.sin(Brad), math.cos(Brad), math.tan(Brad)
    return (l * sinB
            + l ** 3 / 3.0 * sinB * cosB ** 2 * (1.0 + 3.0 * eta2 + 2.0 * eta2 ** 2)
            + l ** 5 / 15.0 * sinB * cosB ** 4 * (2.0 - t * t))


def fwd(Bdeg, Ldeg, L0deg, a, f):
    """高斯正算：返回 (x, y, γ/度)。y 为未加常数/带号的原始横坐标。"""
    B = math.radians(Bdeg)
    l = math.radians(Ldeg - L0deg)
    e2 = 2.0 * f - f * f
    ep2 = e2 / (1.0 - e2)
    sinB, cosB, t = math.sin(B), math.cos(B), math.tan(B)
    eta2 = ep2 * cosB * cosB
    N = a / math.sqrt(1.0 - e2 * sinB * sinB)
    X0 = meridian_arc(B, a, f)
    x = (X0 + N / 2.0 * sinB * cosB * l ** 2
         + N / 24.0 * sinB * cosB ** 3 * (5.0 - t * t + 9.0 * eta2 + 4.0 * eta2 ** 2) * l ** 4
         + N / 720.0 * sinB * cosB ** 5
         * (61.0 - 58.0 * t * t + t ** 4 + 270.0 * eta2 - 330.0 * eta2 * t * t) * l ** 6)
    y = (N * cosB * l
         + N / 6.0 * cosB ** 3 * (1.0 - t * t + eta2) * l ** 3
         + N / 120.0 * cosB ** 5
         * (5.0 - 18.0 * t * t + t ** 4 + 14.0 * eta2 - 58.0 * eta2 * t * t) * l ** 5)
    gam = _gamma_eta(l, B, eta2)
    return x, y, math.degrees(gam)


def inv(x, y, L0deg, a, f, max_iter=80):
    """高斯反算：返回 (B/度, L/度, γ/度)。y 为未加常数/带号的原始横坐标。"""
    e2 = 2.0 * f - f * f
    ep2 = e2 / (1.0 - e2)
    # 底点纬度迭代：Bf ← Bf + (x − X(Bf))·(1−e²sin²Bf)^{3/2}/(a(1−e²))
    Bf = x / (a * (1.0 - e2) * (1.0 + 0.75 * e2))
    for _ in range(max_iter):
        d = (x - meridian_arc(Bf, a, f)) \
            * math.sqrt(1.0 - e2 * math.sin(Bf) ** 2) ** 3 / (a * (1.0 - e2))
        Bf += d
        if abs(d) < 1e-15:
            break
    sinB, cosB, t = math.sin(Bf), math.cos(Bf), math.tan(Bf)
    eta2 = ep2 * cosB * cosB
    N = a / math.sqrt(1.0 - e2 * sinB * sinB)
    B = (Bf
         - t / (2.0 * N ** 2) * (1.0 + eta2) * y ** 2
         + t / (24.0 * N ** 4)
         * (5.0 + 3.0 * t * t + 6.0 * eta2 - 6.0 * eta2 * t * t
            - 3.0 * eta2 ** 2 - 9.0 * eta2 ** 2 * t * t) * y ** 4
         - t / (720.0 * N ** 6) * (61.0 + 90.0 * t * t + 45.0 * t ** 4) * y ** 6)
    l = (y / (N * cosB)
         - (1.0 + 2.0 * t * t + eta2) * y ** 3 / (6.0 * N ** 3 * cosB)
         + (5.0 + 28.0 * t * t + 24.0 * t ** 4 + 6.0 * eta2 + 8.0 * eta2 * t * t)
         * y ** 5 / (120.0 * N ** 5 * cosB))
    # γ 用**反算所得的最终纬度 B**（非底点纬度 Bf）代入收敛角级数：
    # 原著「精度可达 0″.0001」，若用 Bf 则 γ 偏差达 0.26″（已反证，见 verify）。
    eta2b = ep2 * math.cos(B) ** 2
    gam = _gamma_eta(l, B, eta2b)
    return math.degrees(B), L0deg + math.degrees(l), math.degrees(gam)


def zone_number(Ldeg, width):
    """投影带号：6°带 N=int(L/6)+1（λ₀=6N−3）；3°带 N'=int((L+1.5)/3)（λ₀=3N'）。"""
    if width == 6:
        return int(Ldeg / 6.0) + 1
    if width == 3:
        return int((Ldeg + 1.5) / 3.0)
    raise ValueError("投影带宽只能为 3 或 6（原著注意事项第 6 条）")


# ============================================================
# 解析
# ============================================================

def _num(tok):
    return float(str(tok).strip())


def parse(data):
    """
    解析 .INT（路径）或 dict。
    dict 形式：{"name","author","date","mode","coord","K","flag",
                "zw","Lo" | "zw1","L1","zw2","L2","points":[(名, 值1, 值2), …]}
    """
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("flag", 1)
        p.setdefault("zw", 0)
        p["points"] = [(str(a), float(b), float(c)) for (a, b, c) in p["points"]]
        p["mode"] = str(p["mode"])
        p["coord"] = int(p["coord"])
        return p

    lines = read_lines(data)
    if not lines:
        raise ValueError("F-1X 输入文件为空")
    head = [t for t in lines[0].split(",")]
    if len(head) < 7:
        raise ValueError("F-1X 首行字段不足 7 个：%r" % lines[0])
    p = {
        "name": head[0].strip().strip('"'),
        "author": head[1].strip().strip('"'),
        "date": head[2].strip().strip('"'),
        "mode": head[3].strip(),
        "coord": int(_num(head[4])),
        "K": _num(head[5]),
        "flag": int(_num(head[6])),
    }
    if p["mode"] not in MODE_CN:
        raise ValueError("F-1X 类型只能为 + / - / T，实为 %r" % p["mode"])
    tail = [_num(t) for t in head[7:] if t.strip()]
    need = 4 if p["mode"] == "T" else 2
    idx = 1
    while len(tail) < need and idx < len(lines):
        tail.extend(_num(t) for t in lines[idx].split(",") if t.strip())
        idx += 1
    if len(tail) < need:
        raise ValueError("F-1X 首部参数不足（%s 模式需 %d 个尾参，实得 %d）"
                         % (p["mode"], need, len(tail)))
    if p["mode"] == "T":
        p["zw1"], p["L1"], p["zw2"], p["L2"] = tail[0], tail[1], tail[2], tail[3]
        p["zw"], p["Lo"] = p["zw1"], p["L2"]
    else:
        p["zw"], p["Lo"] = tail[0], tail[1]
    # 其余行：点名行 + 坐标行
    pts = []
    i = idx
    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            i += 1
            continue
        nm = ln.strip().strip('"')
        i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            raise ValueError("F-1X 点名 %r 之后缺少坐标行" % nm)
        cs = [_num(t) for t in lines[i].split(",") if t.strip()]
        if len(cs) < 2:
            raise ValueError("F-1X 点名 %r 坐标行需 2 个数，实得 %d" % (nm, len(cs)))
        pts.append((nm, cs[0], cs[1]))
        i += 1
    if not pts:
        raise ValueError("F-1X 输入文件无测点")
    p["points"] = pts
    return p


# ============================================================
# 计算
# ============================================================

def compute(p):
    ename, a, f = ELLIPSOIDS[p["coord"]]
    mode = p["mode"]
    K = p["K"]
    out = []
    for nm, u, v in p["points"]:
        rec = {"点名": nm}
        if mode == "+":
            B = dmss_to_deg(u)
            L = dmss_to_deg(v)
            x, y, gam = fwd(B, L, p["Lo"], a, f)
            rec.update({"B": B, "L": L, "x": x, "y": y, "gamma": gam,
                        "B_in": B, "L_in": L})
        elif mode == "-":
            x, y = u, v - K * 1000.0
            B, L, gam = inv(x, y, p["Lo"], a, f)
            # 反算的 Ｘ/Ｙ 两列为输入值回显（原著表格无第二组坐标列）
            rec.update({"B": B, "L": L, "x": u, "y": y, "gamma": gam,
                        "echo_xy": True})
        else:                                   # T 换带
            x1, y1 = u, v - K * 1000.0
            B, L, _g1 = inv(x1, y1, p["L1"], a, f)
            x, y, gam = fwd(B, L, p["L2"], a, f)
            rec.update({"B": B, "L": L, "x": x, "y": y, "gamma": gam,
                        "x_in": u, "y_in": y1})
        # 输出 Y：原始横坐标 + K 加常数（km→m），带号前缀按 flag。
        # 反算模式 Ｘ/Ｙ 为输入值回显（原著表格无第二组坐标列），不加带号与常数。
        if rec.get("echo_xy"):
            y_out = v
            rec["Y_out"] = y_out
            out.append(rec)
            continue
        y_out = y + K * 1000.0
        if p["flag"] and p["zw"] in (3, 6):
            y_out += zone_number(rec["L"], p["zw"]) * 1000000
        rec["Y_out"] = y_out
        out.append(rec)
    return {
        "程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
        "椭球": ename, "椭球a": a, "椭球1/f": 1.0 / f,
        "输入": p, "成果": out,
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

_LINE = " " + "*" * 71
_BANNER = (" ****" + " " * 16 + "高 斯 投 影 计 算" + "  " + "F-1X" + "  "
           + "(99.6版)" + " " * 14 + "****")
_BORDER_TOP = "┌──┬────────┬───────┬───────┬──────┬──────┬───────┐"
_BORDER_MID = "├──┼────────┼───────┼───────┼──────┼──────┼───────┤"
_BORDER_BOT = "└──┴────────┴───────┴───────┴──────┴──────┴───────┘"
_HEADER = ("│ No.│   点      名   │      Ｂ      │      Ｌ      │     Ｘ     │     Ｙ     │"
           "      γ\u3000\u3000  │")


def _head_line(p):
    """测区行：'     测区名称:K' 左对齐 51 列 + 类型名 18 列 + Lo/K 尾串。"""
    left = ("     测区名称:" + str(p["name"])).ljust(51)
    mid = MODE_CN[p["mode"]].ljust(18)
    if p["mode"] == "T":
        tail = ("L1=" + "%3d" % int(round(p["L1"])) + "°   L2="
                + "%3d" % int(round(p["L2"])) + "°   K=" + "%4d" % int(round(p["K"]))
                + " km")
    else:
        tail = ("Lo=" + "%3d" % int(round(p["Lo"])) + "°" + "   "
                + "K=" + "%4d" % int(round(p["K"])) + " km")
    return left + mid + tail


def render(p, r):
    A = []
    A.append("")
    A.append(_LINE)
    A.append(_BANNER)
    A.append(_LINE)
    A.append("")
    A.append(_head_line(p))
    A.append(_BORDER_TOP)
    A.append(_HEADER)
    A.append(_BORDER_MID)
    for i, rec in enumerate(r["成果"]):
        A.append("│" + "%3d " % (i + 1) + "│" + str(rec["点名"]).ljust(16) + "│"
                 + fmt_dms(rec["B"]) + "│" + fmt_dms(rec["L"]) + "│"
                 + "%12.3f" % rec["x"] + "│" + "%12.3f" % rec["Y_out"] + "│"
                 + fmt_gamma(rec["gamma"]) + "│")
    A.append(_BORDER_BOT)
    A.append("    计算者:" + str(p["author"]) + "     日期: " + str(p["date"]) + " ")
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


# ============================================================
# 知识库对照结果（D:\WorkBuddy知识库\水利知识库\）
# ============================================================
#  母本：02_水利教材精读\水利工程施工与管理\水利工程测量_第5版_精读笔记.md
#
#  | 程序公式 | 库中出处 | 是否一致 | 处置 |
#  |---|---|---|---|
#  | 6°带中央子午线 λ₀ = 6N − 3 | 笔记 §2.7 / 公式表 10-1 | **一致** | 采用；
#  |   |  |  | zone_number(L,6)=int(L/6)+1 |
#  | 3°带中央子午线 λ₀ = 3N′ | 笔记 §2.7 / 公式表 10-2 | **一致** | 采用；
#  |   |  |  | zone_number(L,3)=int((L+1.5)/3) |
#  | 横坐标西移 500 km 并加带号 | 笔记 §2.7「横坐标西移 500km 加带号」 | **一致** | 采用；
#  |   |  |  | K 加常数=500（km）+ 带号×10⁶ |
#  | 椭球元素：克拉索夫斯基 a=6378245，1:298.3 | 笔记 §2.1「参考椭球…克拉索夫斯基椭球（1940，a=6378245m，1:298.3）」 | **一致** | 采用（坐标系=1） |
#  | 椭球元素：1975 国际椭球 | 笔记 §2.1「1975 国际椭球」 | **数值库中未给**； | 采用 a=6378140，1/f=298.257；
#  |   |  | 程序说明书称「1980西安坐标系」 | 计 DECL（无权威 OUT 反证） |
#  | 正形投影（角度无变形） | 笔记 §2.7「正形投影（角度无变形）」 | **一致** | 采用标准等角投影级数 |
#  | 高斯投影正算 x,y 级数 | 笔记 §3.3 只给 λ₀ 两式，未给展开式 | 库中**未收** | 取《大地测量学》标准
#  |   |  |  | 展开式（l⁶ / y⁶ 阶），由权威 OUT 逐位验证 |
#  | 子午线弧长级数 | 库中未收 | 库中**未收** | 同上，由权威 OUT 逐位验证 |
#  | 子午线收敛角 γ（两点子午线夹角） | 笔记 §2.4「子午线收敛角（两点子午线夹角）」 | 概念**一致**，展开式库中未收 | 采用 l·sinB 级数，
#  |   |  |  | 由权威 OUT γ=−0°25′34.8078″ 逐位验证 |
#
#  结论：F 类（测量）在知识库中的唯一母本即《水利工程测量》第5版精读笔记。
#  该笔记 §2.7/§3.3 **给出了分带与中央子午线公式（10-1/10-2）与 500km 带号约定**，
#  与本内核完全一致；**未给出**高斯正反算的级数展开式与子午线弧长展开式
#  （教材第10章只讲地形图应用，级数式在《大地测量学》），
#  本内核的级数形式改由权威 F-1X-1.OUT 逐位反演确认（X/Y/γ 三项全中）。
#
#  知识产权：本改造工作基于公之于众版（乌鲁木齐正海水利科技有限公司，张校正教授级高工
#  技术总负责）的公开算法。改造实现（Python 代码、架构设计、验证数据、自动化流程）为
#  哈胜的原创成果。

if __name__ == "__main__":
    import sys
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
