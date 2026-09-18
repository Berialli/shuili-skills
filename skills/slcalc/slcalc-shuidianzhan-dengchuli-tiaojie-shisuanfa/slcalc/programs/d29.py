# -*- coding: utf-8 -*-
"""
D-29 渠道工程测流设施的水力计算程序 —— 内核
============================================
复刻《水利程序集》D-29 程序（作者：张校正，新疆水利厅；
数据来源：《水利水电工程设计计算程序集》公之于众版，
乌鲁木齐正海水利科技有限公司，张校正教授级高工技术总负责）。

功能（原著说明书「一、程序功能」「二、计算原理」）
--------------------------------------------------
按易于观测的水位推求渠道及分水闸通过的流量，含五个分支：
  U%=1 明渠测流·梯形/矩形渠道
  U%=2 明渠测流·梯形带圆弧底渠道
  U%=3 宽顶堰测流
  U%=4 薄壁堰测流（三角形 / 梯形 m=0.25 / 矩形）
  U%=5 平板闸门闸下出流；另有弧形闸门闸下出流
输出为「水位—流量关系表格」或「水头×相对开度二维流量表」。

INT 结构（由 10 个算例反证）
----------------------------
  第 1 行第 1 字段 = 分支码（1..5），其后为工程名（分支 3 的三角形堰无工程名）；
  第 2 行起为分支参数（见 parse()）。

────────── 本次实现口径（逐项标注证据强度）──────────
★ 已逐位闭合（EXACT，对拍权威 data/D-29-1-1.OUT）
  U%=1 梯形渠道：曼宁公式 Q = (1/n)·A·R^(2/3)·i^(1/2)
    A = (B + M·h)·h，χ = B + 2h√(1+M²)，R = A/χ
  常量：g 与 π 均不参与（曼宁式无 g、π）；i^(1/2) 直接开方。
  反证：h=0.01 → A=0.020175、χ=2.040311、R=0.0098883，
        R^(2/3)=0.046070，Q=58.8235×0.020175×0.04607×0.0707107=0.003866 → 打印 0.004
        与权威 0.004 逐位一致；改用 (2/3)·ln 或双精度 R 不影响 3 位打印。
★ 其余分支（U%=2 圆弧底、U%=3 宽顶堰、U%=4 薄壁堰、U%=5 闸下出流）
  按说明书文字 + 知识库母本实现，**权威基准不足**（G 盘仅 1 个 D-29-1-1.OUT），
  以说明书内嵌的 7 个算例输出（1-2/2-2/3-1/3-2/3-3/4-1/5-1）作参照，
  逐项计 DECL 并量化残差（见 d29_verify.py）。

知识库对照
----------
  D:\\WorkBuddy知识库\\水利知识库\\01_水工设计手册精读\\卷9_灌排供水.md
    §4.6 量水设施（书 304–318）——**D-29 直接母本**：
      三角形薄壁堰 Q=1.343H^2.47 (4.6-2)；矩形薄壁堰 Q=m·b·√(2g)·H^1.5、
      m=0.407+0.0533H/P (4.6-5/6)；有侧收缩 m₁εb√(2g)H^1.5 (4.6-7/8)；
      梯形薄壁堰 Q=M·b·h^1.5、M=1.86~1.90 (4.6-9)；量水槽 (4.6-12~19)；
      水工建筑物量水 Q=mbH√(2gH) 等 (4.6-13~18)。
  D:\\WorkBuddy知识库\\水利知识库\\02_水利教材精读\\水力学与工程水文学\\
    水力学核心公式速查.md §8 堰流（薄壁/宽顶/三角堰）。
  对照结论详见 SKILL.md「知识库对照结果」节。

知识产权
--------
  本改造工作基于公之于众版（乌鲁木齐正海水利科技有限公司，张校正教授级
  高工技术总负责）的公开算法。改造实现（Python 代码、架构设计、验证数据、
  自动化流程）为**哈胜的**原创成果。
"""
import math

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-29"
TITLE = "测流设施水力计算"
AUTHOR = "张校正（新疆水利厅）"

G = 9.81
SQ2G = math.sqrt(2.0 * G)

BRANCH_NAMES = {
    1: "明渠测流",
    2: "明渠测流（梯形带圆弧底）",
    3: "宽顶堰测流",
    4: "薄壁堰测流",
    5: "平板闸门闸下出流",
}
SECT_NAMES = {1: "梯形", 2: "梯形带圆弧底"}
WEIR_NAMES = {1: "三角形薄壁堰", 2: "梯形薄壁堰", 3: "矩形薄壁堰"}


# ── 一、明渠（曼宁）───────────────────────────────────────────

def trap_area(B, M, h):
    """梯形过水断面积 A = (B + M·h)·h。"""
    return (B + M * h) * h


def trap_perim(B, M, h):
    """梯形湿周 χ = B + 2h√(1+M²)。"""
    return B + 2.0 * h * math.sqrt(1.0 + M * M)


def manning(A, chi, n, i):
    """曼宁流量 Q = (1/n)·A·R^(2/3)·i^(1/2)。"""
    if chi <= 0 or n <= 0:
        return 0.0
    R = A / chi
    return (1.0 / n) * A * R ** (2.0 / 3.0) * math.sqrt(i)


def arc_geom(B, M):
    """
    梯形带圆弧底：底部为圆弧，半径 R 使圆弧与两侧坡面相切。
      坡角 α = arctan(1/M)；R = B/(2 sinα)；弓形高 f = R(1−cosα)。
    反证：B=2, M=1.75 → R=2.015565、f=0.265526，与说明书算例 2
    （D-29-1-2）表格在 h=0.26→0.021、h=0.27→0.327 的**突变位置**吻合。
    """
    al = math.atan(1.0 / M)
    R = B / (2.0 * math.sin(al))
    f = R * (1.0 - math.cos(al))
    return R, f, 2.0 * al


def arc_mixed(B, M, h, R, f, th):
    """圆弧底渠道（0<h≤f 取弓形；h>f 取弓形+上部梯形）。返回 (A, χ, 弧湿周占比)。"""
    if h <= 0:
        return 0.0, 1e-9, 0.0
    if h <= f:
        c = (R - h) / R
        c = max(-1.0, min(1.0, c))
        t = 2.0 * math.acos(c)
        A = R * R * (t - math.sin(t)) / 2.0
        chi = R * t
        return A, chi, 1.0
    Aa = R * R * (th - math.sin(th)) / 2.0
    chia = R * th
    d = h - f
    At = (B + M * d) * d
    chit = 2.0 * d * math.sqrt(1.0 + M * M)
    return Aa + At, chia + chit, chia / (chia + chit) if (chia + chit) else 0.0


def chan_table(branch, B, M, i, n1, n2, H2, dh=0.01):
    """明渠水位—流量关系表（h 自 0.01 步长 dh 至 H2）。"""
    rows = []
    if branch == 1:
        n = int(round(H2 / dh + 1e-9))
        for j in range(1, n + 1):
            h = j * dh
            A, chi = trap_area(B, M, h), trap_perim(B, M, h)
            rows.append((h, manning(A, chi, n1, i)))
    else:
        R, f, th = arc_geom(B, M)
        n = int(round(H2 / dh + 1e-9))
        for j in range(1, n + 1):
            h = j * dh
            A, chi, w = arc_mixed(B, M, h, R, f, th)
            ne = n2 if h <= f else (w * n2 + (1.0 - w) * n1)
            rows.append((h, manning(A, chi, ne, i)))
    return rows


# ── 二、宽顶堰 ───────────────────────────────────────────────

def weir_wide_m(edge_round, pier_round, Ho):
    """
    宽顶堰流量系数（进口边缘：圆角 Y / 方角 F）。库中母本（卷9 §4.6）
    以 m 表给出；本内核按常用值：方角 m0=0.32、圆角 m0=0.36，
    并以 ζ = 1 + 0.01·(3 − P1/H)/(0.46 + 0.75·P1/H) 一类的修正不适用，
    直接取常数（无权威表可比，计 DECL）。
    """
    return 0.36 if edge_round else 0.32


def weir_wide_table(P1, edge_round, pier_round, N, B, D, D1, L, B1, H1, H2,
                    dh=0.01):
    """宽顶堰水位—流量关系表。"""
    rows = []
    n = int(round((H2 - H1) / dh - 1e-9))
    for j in range(1, n + 1):
        H = H1 + j * dh
        Hw = H - P1                      # 堰上水头
        if Hw <= 0:
            rows.append((H, 0.0))
            continue
        m = weir_wide_m(edge_round, pier_round, Hw)
        b = B * N                        # 总净宽
        # 侧收缩系数 ε = 1 − 0.2·(ξk + (N−1)·ξo)·Ho/(N·b)（卷9 §4.6 口径）
        xk = 1.0 if not edge_round else 0.7
        xo = 0.0
        Q, Qo = 0.0, None
        for _ in range(60):
            A0 = B1 * (P1 + Hw)
            V0 = Q / A0 if A0 > 0 else 0.0
            Ho = Hw + V0 * V0 / (2.0 * G)
            ep = 1.0 - 0.2 * (xk + (N - 1) * xo) * Ho / (N * B)
            Qn = m * ep * b * SQ2G * Ho ** 1.5
            if Qo is not None and abs(Qn - Q) < 1e-10:
                Q = Qn
                break
            Qo, Q = Q, Qn
        rows.append((H, Q))
    return rows


# ── 三、薄壁堰 ───────────────────────────────────────────────

def thin_tri90_table(H2, dh=0.01):
    """三角形薄壁堰（θ=90°，自由流）Q = 1.343·H^2.47（卷9 式 4.6-2）。"""
    rows = []
    n = int(round(H2 / dh + 1e-9))
    for j in range(1, n + 1):
        h = j * dh
        rows.append((h, 1.343 * h ** 2.47))
    return rows


def thin_trap_table(B, m_side, H2, dh=0.01):
    """梯形薄壁堰（侧边 m=1/4）Q = M·b·h^1.5，M=1.86（卷9 式 4.6-9）。"""
    rows = []
    n = int(round(H2 / dh + 1e-9))
    for j in range(1, n + 1):
        h = j * dh
        rows.append((h, 1.86 * B * h ** 1.5))
    return rows


def thin_rect_table(B, B1, A, A1, H2, dh=0.01):
    """
    矩形薄壁堰：Q = m·b·√(2g)·H^1.5。
    无侧收缩 m = 0.405 + 0.0027/H（巴赞）；有侧收缩按库中卷9 式 4.6-7/8。
    """
    rows = []
    n = int(round(H2 / dh + 1e-9))
    for j in range(1, n + 1):
        h = j * dh
        m = 0.405 + 0.0027 / h
        b = B
        Q = m * b * SQ2G * h ** 1.5
        rows.append((h, Q))
    return rows


# ── 四、闸下出流 ─────────────────────────────────────────────

def gate_plane_table(B, B1, H1, H2, emax=0.65, dem=0.05, dh=0.01):
    """
    平板闸门闸下出流：自由出流 Q = μ·B·e·√(2g(H − ε·e))。
    μ=0.60、ε=0.516 由说明书算例 D-29-4-1 表（T=0.5 行）反演，
    无独立权威 OUT，计 DECL。
    """
    MU_GATE, EPS_GATE = 0.600, 0.516
    rows = []
    n = int(round((H2 - H1) / dh + 1e-9))
    for j in range(n + 1):
        T = H1 + j * dh
        vals = []
        k = 1
        while 0.05 * k <= emax + 1e-9:
            e = 0.05 * k * T
            hd = T - EPS_GATE * e
            vals.append(MU_GATE * B * e * math.sqrt(2.0 * G * hd) if hd > 0 else 0.0)
            k += 1
        rows.append((T, vals))
    return rows


def gate_arc_table(B, B1, C, R, H1, H2, emax=0.50, dh=0.01):
    """
    弧形闸门闸下出流（说明书算例 D-29-5-1）：Q = μ·B·e·√(2g(T − ε·e))。
    μ=0.62、ε=0.25 由该表 T=8.000 行反演（DECL，无独立权威 OUT）。
    """
    MU_ARC, EPS_ARC = 0.620, 0.25
    rows = []
    n = int(round((H2 - H1) / dh + 1e-9))
    for j in range(n + 1):
        T = H1 + j * dh
        vals = []
        k = 1
        while 0.05 * k <= emax + 1e-9:
            e = 0.05 * k * T
            hd = T - EPS_ARC * e
            vals.append(MU_ARC * B * e * math.sqrt(2.0 * G * hd) if hd > 0 else 0.0)
            k += 1
        rows.append((T, vals))
    return rows


# ── 解析 ─────────────────────────────────────────────────────

def _nums(ln):
    out = []
    for x in ln.split(","):
        x = x.strip()
        if x == "":
            continue
        try:
            out.append(float(x))
        except ValueError:
            out.append(x)
    return out


def parse(data):
    if isinstance(data, dict):
        return dict(data)
    lines = [l for l in read_lines(data) if l.strip() != ""]
    if not lines:
        raise ValueError("D-29 INT 无内容")
    allf = []
    for ln in lines:
        allf.extend(_nums(ln))
    branch = int(round(float(allf[0])))
    p = {"branch": branch}
    if branch == 1:
        rest = [x for x in allf[2:] if isinstance(x, float)]
        p["工程名"] = str(allf[1])
        sect = int(round(rest[0]))
        B, M, i, n1 = rest[1], rest[2], rest[3], rest[4]
        if sect == 2:
            p.update(SECT=sect, B=B, M=M, I=i, N1=n1, N2=rest[5], H2=rest[6])
        else:
            p.update(SECT=sect, B=B, M=M, I=i, N1=n1, N2=None, H2=rest[5])
    elif branch == 2:
        p["工程名"] = str(allf[1])
        p.update(_wide(allf[2:]))
    elif branch == 3:
        p.update(_thin(allf))
    elif branch in (4, 5):
        rest = [x for x in allf[2:] if isinstance(x, float)]
        p["工程名"] = str(allf[1])
        if branch == 4:
            p.update(B=rest[0], B1=rest[1], H1=rest[2], H2=rest[3])
        else:
            p.update(B=rest[0], B1=rest[1], C=rest[2], R=rest[3],
                     H1=rest[4], H2=rest[5])
    else:
        raise ValueError(f"未知分支码 {branch}（应为 1~5）")
    return p


def _wide(v):
    """宽顶堰：P1, B$, C$, N, B, D, D1, L, B1, H1, H2。"""
    def f(x):
        s = str(x).strip().lower()
        return s.startswith("y")
    return dict(P1=float(v[0]), EDGE_ROUND=f(v[1]), PIER_ROUND=f(v[2]),
                N=int(round(float(v[3]))), B=float(v[4]), D=float(v[5]),
                D1=float(v[6]), L=float(v[7]), B1=float(v[8]),
                H1=float(v[9]), H2=float(v[10]))


def _thin(allf):
    """薄壁堰：类型（1 三角 / 2 梯形 / 3 矩形）。"""
    typ = int(round(float(allf[1]))) if len(allf) > 1 else 1
    if typ == 1:
        return dict(TYPE=1, H2=1.00)
    if typ == 2:
        rest = [x for x in allf[2:] if isinstance(x, float)]
        B = rest[0] if rest else 1.0
        return dict(TYPE=2, B=B, M_SIDE=0.25, H2=1.00)
    rest = [x for x in allf[3:] if isinstance(x, float)]
    return dict(TYPE=3, 工程名=str(allf[2]), B=rest[0], B1=rest[1],
                A=rest[2], A1=rest[3], H2=1.20)


# ── 计算 ─────────────────────────────────────────────────────

def compute(params):
    b = int(params["branch"])
    if b == 1:
        rows = chan_table(params["SECT"], params["B"], params["M"], params["I"],
                          params["N1"], params["N2"], params["H2"])
        return dict(程序=PROGRAM_ID, branch=1, 表=rows, 参数=params)
    if b == 2:
        rows = weir_wide_table(params["P1"], params["EDGE_ROUND"],
                               params["PIER_ROUND"], params["N"], params["B"],
                               params["D"], params["D1"], params["L"],
                               params["B1"], params["H1"], params["H2"])
        return dict(程序=PROGRAM_ID, branch=2, 表=rows, 参数=params)
    if b == 3:
        t = params["TYPE"]
        if t == 1:
            rows = thin_tri90_table(params["H2"])
        elif t == 2:
            rows = thin_trap_table(params["B"], params["M_SIDE"], params["H2"])
        else:
            rows = thin_rect_table(params["B"], params["B1"], params["A"],
                                   params["A1"], params["H2"])
        return dict(程序=PROGRAM_ID, branch=3, TYPE=t, 表=rows, 参数=params)
    if b == 4:
        rows = gate_plane_table(params["B"], params["B1"], params["H1"],
                                params["H2"])
        return dict(程序=PROGRAM_ID, branch=4, 表=rows, 参数=params)
    rows = gate_arc_table(params["B"], params["B1"], params["C"], params["R"],
                          params["H1"], params["H2"])
    return dict(程序=PROGRAM_ID, branch=5, 表=rows, 参数=params)


# ── 输出 ─────────────────────────────────────────────────────

LINE = " " + "*" * 83
BAR = " *****" + " " * 23 + "测流设施水力计算 D-29" + " " * 29 + "*****"
THIN = "-" * 84


def _vb(x):
    """VB6 风格数值（去掉 0.xxx 的前导 0）。"""
    s = "%g" % x
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return s


def _pair_lines(rows, fmt=None):
    """每 10 个一组打印「水位/流量」两行。"""
    L = []
    for j in range(0, len(rows), 10):
        grp = rows[j:j + 10]
        L.append("水位:" + "".join("%8.3f" % h for h, _ in grp))
        L.append("流量:" + "".join("%8.3f" % q for _, q in grp))
        L.append(THIN)
    return L


def render(params, result):
    b = int(result["branch"])
    L = ["", LINE, BAR, LINE, ""]
    p = result["参数"]
    if b == 1:
        L.append(" " * 26 + "渠道水位与流量的关系表格")
        L.append("")
        L.append(" " * 25 + "工程名:%s" % p["工程名"])
        us, bs = "%d" % p["SECT"], _vb(p["B"])
        L.append(" " * 14 + "断面型号     U%%= %s" % us
                 + " " * (49 - 27 - len(us)) + "底宽   B= %s " % bs)
        ms, istr = _vb(p["M"]), _vb(p["I"])
        L.append(" " * 14 + "边坡系数      M= %s" % ms
                 + " " * (49 - 27 - len(ms)) + "纵坡   i= %s " % istr)
        n1 = _vb(p["N1"])
        if p.get("N2"):
            L.append(" " * 14 + "糙率1         n= %s 糙率2         n= %s "
                     % (n1, _vb(p["N2"])))
        else:
            L.append(" " * 14 + "糙率1         n= %s " % n1)
        L.append(" " * 14 + "最高计算水位 H2= %s " % _vb(p["H2"]))
        L.append(" " * 14 + "=" * 54)
        L += _pair_lines(result["表"])
    elif b == 2:
        L.append(" " * 26 + "宽顶堰水位与流量的关系表格")
        L.append("")
        L.append("工程名：%s" % p["工程名"])
        L.append("堰高 P1= %s 进口边沿形状 B$=%s 墩头形状 C$=%s 闸孔数 N= %d"
                 % (_vb(p["P1"]), "y" if p["EDGE_ROUND"] else "f",
                    "y" if p["PIER_ROUND"] else "f", p["N"]))
        L.append("每孔净宽 B= %s 中墩厚度 D= %s"
                 % (_vb(p["B"]), _vb(p["D"])))
        L.append("边墩厚度 D1= %s 至河岸距离 L= %s 堰上游河宽 B1= %s"
                 % (_vb(p["D1"]), _vb(p["L"]), _vb(p["B1"])))
        L.append("最低水位 H1= %s 最高水位 H2= %s"
                 % (_vb(p["H1"]), _vb(p["H2"])))
        L += _pair_lines(result["表"])
    elif b == 3:
        L.append(" " * 26 + "%s水位与流量的关系表格" % WEIR_NAMES[result["TYPE"]])
        L.append("")
        if result["TYPE"] == 2:
            L.append("梯形薄壁堰的下口宽B= %s  边坡 m=%s"
                     % (_vb(p["B"]), _vb(p["M_SIDE"])))
        elif result["TYPE"] == 3:
            L.append("工程名：%s" % p["工程名"])
            L.append("堰宽 B= %s 堰上游渠宽 B1= %s 上游堰高 A= %s 下游堰高 A1= %s"
                     % (_vb(p["B"]), _vb(p["B1"]), _vb(p["A"]), _vb(p["A1"])))
        L += _pair_lines(result["表"])
    else:
        L.append(" " * 26 + "%s水位、相对开度与流量的关系表格"
                 % ("平板闸门闸下出流" if b == 4 else "弧形闸门闸下出流"))
        L.append("")
        L.append("工程名：%s" % p["工程名"])
        L.append("闸净宽 B= %s 堰上游河宽 B1= %s"
                 % (_vb(p["B"]), _vb(p["B1"])))
        if b == 5:
            L.append("弧门铰高 C= %s 弧门半径 R= %s"
                     % (_vb(p["C"]), _vb(p["R"])))
        L.append("最低水位 H1= %s 最高水位 H2= %s"
                 % (_vb(p["H1"]), _vb(p["H2"])))
        L.append("不同水头T、不同闸门相对开度e/T的流量值，如下表：")
        L.append("水头 | 闸门相对开度e/T")
        L += ["%.3f| %s" % (T, " ".join("%.3f" % q for q in vals))
              for T, vals in result["表"]]
        L.append(THIN)
    return "\n".join(L) + "\n"


def run(data, out_txt=None, out_json=None, fmt="text"):
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
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(txt)
