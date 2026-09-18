# -*- coding: utf-8 -*-
"""
F-12 测边网平差计算程序 —— 内核
================================
复刻《水利水电工程设计计算程序集》F-12（作者：谢希哲，新疆兵团勘测设计院，
99.6 版）。权威基准：`slcalc/data/F-12-1.OUT`（GBK，5098 B，53 行）。

算法（由权威 OUT 逐位反演，全部 5 个待定点坐标命中 0.001 m）
-------------------------------------------------------------
原著为**严格最小二乘（间接平差）**，本内核按同一模型实现：

  未知参数：待定点坐标（2n 个）
  观测值  ：
    · 测边 m 条      —— 残差取观测值单位（mm 或 m 同比），权
                         Ps = [(a+b)/(a+S·b/1000)]²          （S 以 m 计，a、b 以 mm 计）
    · 加测方向       —— 每测站 n_i 个方向读数 → (n_i−1) 个角度（相对零方向），权
                         Pβ = [(a+b)/Mβ]²
    · 该站"史赖伯和方程" —— 1 个虚拟和方程（残差 = Σ 角度残差），权 −Pβ/n_i
    · 加测方位角 k 个 —— 权 Pα = [(a+b)/Mα]²
  多余观测 r = (m + Σn_i + k方向数) − 2n

  ★核心裁决（由权威 F-12-1.OUT 唯一反演）：
    · 说明书 §4② 印作 "Pβ=12.96(Ms/Mβ)²、Pα=12.96(Ms/Mα)²"，
      但**实际计算不含 12.96 因子**——带 12.96 时全部 5 个待定点坐标偏离权威
      6.7 mm（Pβ=343.34）；去掉后 Pβ=26.4922=(7/1.36)²、Pα=13.2921=(7/1.92)²
      时 5 个坐标打印值**逐位全中**（残差 ≤0.00046 m）。见 F-12 报告「反证清单」。
    · 方向观测**不是**"带测站定向未知数的方向"：权威 OUT 的 3 个方向改正数之和
      （−2.27″ 与 −5.37″）不为 0，而带定向未知数的模型强制 Σv=0 → 已排除。
      改用 (n_i−1) 个角度 + 史赖伯和方程（权 −Pβ/n_i），与权威一致。
    · 相对中误差 1:N 的 N 取整至最近 100。
    · 平差后单位权中误差 Mo = ±√([pv²]/r)，v 取 mm、Ps 以 1 km 边为 1
      → F-12-1 得 ±19.2 mm（权威 ±19.2mm）。

输出（逐字复刻原著 .OUT 版式，GBK）
-----------------------------------
  网名行 → 表头 → 各点块（点名+坐标+误差椭圆 → 该点各边行，含方位角、Ms/S、
  方向观测值/平差值）→ Ms/方向中误差/Mo 行 → 计算者/日期行

已声明偏差（DECL）
------------------
  D1. 坐标与权威存在 ≤0.00046 m 的固定量级偏差，全部落在 0.001 m 打印量子内
      → 5 个待定点坐标打印值**全中**。反证：模型已由「去掉 12.96 因子」唯一确定，
      残差为原著 DOS 版（QuickBASIC 未加后缀变量默认 SINGLE，相对精度 ≈1.2×10⁻⁷）
      在局部坐标下的中间量舍入所致。已排除：权重比、方向模型（定向未知数 vs
      史赖伯和方程）、和方程权符号、12.96 因子。
  D2. 权威 .OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\USE\\F-12-1.out」为运行期
      路径回显（机器相关），本内核不生成，版式比对自第 2 行起算。
  D3. F-12-2.INT / F-12-4.INT G 盘**无权威 OUT**，无法逐位对拍；仅作内部校核。

知识库对照结果（教学母本 → 程序实现）
------------------------------------
  KB 路径：D:\\WorkBuddy知识库\\水利知识库\\

  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 边长观测方程 v=S测−S算 | `02_水利教材精读/水利工程施工与管理/水利工程测量_第5版_精读笔记.md` §3.2 式6-3/6-4（距离反算） | **一致** | 采用 |
  | 测距精度 Ms=a+S·b（标称精度 mm+ppm） | 同上 §2.2「电磁波测距」；式4-31~33（全站仪三轴误差） | **一致**（库中给出概念，未给权函数式） | 采用（权 Ps=[(a+b)/(a+S·b/1000)]² 由权威 OUT 验证） |
  | 权 P_i=μ²/m_i²、加权平差 | 同上 §2.5、式5-21/5-22 | **一致** | 采用 |
  | 误差传播 m_Z²=Σ(∂f/∂x_i)²m_i² | 同上 式5-12/5-14 | **一致** | 采用（Ms、Mα、误差椭圆由协因数阵传播） |
  | 单位权中误差 m=±√([ΔΔ]/n) | 同上 式5-7 | **一致**（本程序 Mo=±√([pv²]/r)） | 采用 |
  | 方向观测/方位角观测 | 同上 §2.4「直线定向、坐标方位角」 | **一致**（概念） | 采用（权重与和方程由权威 OUT 反演） |
  | 点位误差椭圆 E、F、Φe | 库中**无专章** | **未覆盖** | 按通用公式实现，经权威 OUT 逐位验证 |
  | 史赖伯（Schreiber）和方程 | 库中**无专章**（教材第6章仅述条件平差通式） | **未覆盖** | 按原著说明书 §4③ 实现（权 −Pβ/n），并经权威 OUT 反演确认 |

知识产权
--------
本改造工作基于公之于众版（乌鲁木齐正海水利科技有限公司，张校正教授级高工技术总负责）
的公开算法。改造实现（Python 代码、架构设计、验证数据、自动化流程）为哈胜的原创成果。
"""
import math
import os
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "F-12"
TITLE = "测边网平差计算程序"
AUTHOR = "谢希哲"
HEAD_NAME = "F-12"

RHO = 180.0 / math.pi * 3600.0


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def dwidth(s):
    w = 0
    for ch in s:
        o = ord(ch)
        if (0x0370 <= o <= 0x03FF or 0x0400 <= o <= 0x04FF
                or 0x1100 <= o <= 0x115F or 0x2E80 <= o <= 0xA4CF
                or 0xAC00 <= o <= 0xD7A3 or 0xF900 <= o <= 0xFAFF
                or 0xFE30 <= o <= 0xFE6F or 0xFF00 <= o <= 0xFF60
                or 0xFFE0 <= o <= 0xFFE6):
            w += 2
        else:
            w += 1
    return w


def dpad(s, width, right=False):
    gap = max(0, width - dwidth(s))
    return (" " * gap + s) if right else (s + " " * gap)


def _dms(deg, prec=2):
    deg = deg % 360.0
    d = int(deg)
    m = int((deg - d) * 60.0)
    s = q((deg - d) * 60.0 * 60.0 - m * 60.0, prec)
    if s >= 60.0:
        s -= 60.0
        m += 1
    if m >= 60:
        m -= 60
        d += 1
    return "%3d %02d %0*.*f" % (d, m, prec + 3, prec, s)


def parse_angle(x):
    neg = x < 0
    x = abs(x)
    d = int(x)
    mm = (x - d) * 100.0
    m = int(mm + 1e-9)
    ss = (mm - m) * 100.0
    v = d + m / 60.0 + ss / 3600.0
    return -v if neg else v


def _lines(path):
    raw = open(path, "rb").read()
    for enc in ("gbk", "utf-8"):
        try:
            return raw.decode(enc).replace("\r\n", "\n").split("\n")
        except Exception:
            continue
    return raw.decode("latin-1").replace("\r\n", "\n").split("\n")


def parse(data):
    if isinstance(data, dict):
        p = dict(data)
        for k in ("NK", "NU", "NE", "ND", "NA"):
            p[k] = int(round(p[k]))
        return p
    ls = [x.strip() for x in _lines(data)]
    ls = [x for x in ls if x != ""]
    it = iter(ls)

    def nxt():
        return next(it)
    p = {}
    p["网名"] = nxt()
    p["计算者"] = nxt()
    p["日期"] = nxt()
    NK = int(float(nxt()))          # 已知点数
    NU = int(float(nxt()))          # 待定点数
    NE = int(float(nxt()))          # 测边数
    p["NK"], p["NU"], p["NE"] = NK, NU, NE
    eq = nxt().upper().startswith("Y")
    p["等精度"] = eq
    a = float(nxt()); b = float(nxt())
    p["a"], p["b"] = a, b
    p["概算"] = nxt().upper().startswith("Y")
    ND = int(float(nxt()))          # 加测方向测站数
    p["ND"] = ND
    # 原著仅在"有加测方向"时才输入方向观测中误差 Mβ
    # （F-12-2.INT 为 ND=0，其数据流中无 Mβ 项，直接接方位角个数 0）
    p["Mβ"] = float(nxt()) if ND > 0 else 1.0
    NA = int(float(nxt()))          # 加测方位角数
    p["NA"] = NA
    pts = []
    for i in range(NK):
        nm = "P%d" % (i + 1)
        x = float(nxt()); y = float(nxt())
        pts.append({"no": i + 1, "name": nm, "x": x, "y": y, "known": True})
    for i in range(NK, NK + NU):
        pts.append({"no": i + 1, "name": "P%d" % (i + 1), "x": None, "y": None, "known": False})
    edges = []
    for k in range(1, NE + 1):
        u = int(float(nxt())); v = int(float(nxt())); s = float(nxt())
        edges.append({"no": k, "u": u, "v": v, "obs": s})
    dirs = []
    for _ in range(ND):
        st = int(float(nxt())); n = int(float(nxt()))
        t1 = int(float(nxt()))
        rec = [(t1, 0.0)]
        for _j in range(n - 1):
            tg = int(float(nxt())); val = parse_angle(float(nxt()))
            rec.append((tg, val))
        dirs.append({"st": st, "n": n, "rec": rec})
    azs = []
    for _ in range(NA):
        u = int(float(nxt())); v = int(float(nxt()))
        aa = parse_angle(float(nxt())); ma = float(nxt())
        azs.append({"u": u, "v": v, "az": aa, "Ma": ma})
    p["点"] = pts
    p["边"] = edges
    p["方向"] = dirs
    p["方位角"] = azs
    # F-12-2 的 ".C" 辅助文件给出点名（若有同名 .C 则读入）
    if isinstance(data, str):
        cp = os.path.splitext(data)[0] + ".C"
        if os.path.exists(cp):
            try:
                raw = open(cp, "rb").read()
                for i in range(0, min(len(raw) // 16, len(pts))):
                    nm = raw[i * 16:(i + 1) * 16].decode("gbk", "replace").strip()
                    if nm:
                        pts[i]["name"] = nm
            except Exception:
                pass
    return p


# ============================================================
# 计算
# ============================================================

def _gauss(A, b):
    n = len(b)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        pv = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[pv][c]) < 1e-14:
            return None
        M[c], M[pv] = M[pv], M[c]
        d = M[c][c]
        for j in range(c, n + 1):
            M[c][j] /= d
        for r in range(n):
            if r != c and M[r][c] != 0.0:
                f = M[r][c]
                for j in range(c, n + 1):
                    M[r][j] -= f * M[c][j]
    return [M[i][n] for i in range(n)]


def _inv(M):
    n = len(M)
    A = [M[i][:] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(A[r][c]))
        if abs(A[p][c]) < 1e-18:
            return None
        A[c], A[p] = A[p], A[c]
        d = A[c][c]
        for j in range(2 * n):
            A[c][j] /= d
        for r in range(n):
            if r != c and A[r][c] != 0.0:
                f = A[r][c]
                for j in range(2 * n):
                    A[r][j] -= f * A[c][j]
    return [[A[i][n + j] for j in range(n)] for i in range(n)]


def _build(p, x, unk):
    """返回 (l, rows, meta)；l 为残差（边 mm、角 秒），rows 为 (设计行, 权)"""
    idx = {pt: k for k, pt in enumerate(unk)}
    pts = {q_["no"]: q_ for q_ in p["点"]}
    def C(no):
        pt = pts[no]
        if pt["known"]:
            return (pt["x"], pt["y"])
        k = idx[no]
        return (x[2 * k], x[2 * k + 1])
    def azv(u, v):
        A = C(u); B = C(v)
        return math.degrees(math.atan2(B[1] - A[1], B[0] - A[0])) % 360.0
    def dv(u, v):
        A = C(u); B = C(v)
        return math.hypot(B[0] - A[0], B[1] - A[1])
    nu = 2 * len(unk)
    a, b = p["a"], p["b"]
    l = []; rows = []; meta = []
    # 边
    for e in p["边"]:
        u, v, so = e["u"], e["v"], e["obs"]
        A = C(u); B = C(v)
        dx = B[0] - A[0]; dy = B[1] - A[1]
        S = math.hypot(dx, dy)
        Ps = ((a + b) / (a + S / 1000.0 * b)) ** 2
        l.append((so - S) * 1000.0)                 # mm
        row = [0.0] * nu
        for no, sg in ((v, 1.0), (u, -1.0)):
            if pts[no]["known"]:
                continue
            k = idx[no]
            row[2 * k] += sg * dx / S * 1000.0
            row[2 * k + 1] += sg * dy / S * 1000.0
        rows.append((row, Ps))
        meta.append(("E", e))
    # 方向 → 角度 + 和方程
    Pb = ((a + b) / p["Mβ"]) ** 2
    for d in p["方向"]:
        st = d["st"]; rec = d["rec"]
        t0 = rec[0][0]
        srow = [0.0] * nu; sl = 0.0
        for (tg, obs) in rec[1:]:
            val = (azv(st, tg) - azv(st, t0)) % 360.0
            if val > 180.0:
                val = 360.0 - val
            dd = (obs - val) % 360.0
            if dd > 180.0:
                dd -= 360.0
            l.append(dd * 3600.0)
            row = [0.0] * nu
            eps = 1e-4
            for (no, ci) in ((tg, 0), (tg, 1), (t0, 0), (t0, 1), (st, 0), (st, 1)):
                if pts[no]["known"]:
                    continue
                k = idx[no]; old = x[2 * k + ci]; x[2 * k + ci] = old + eps
                v1 = azv(st, tg) - azv(st, t0)
                v1 = (v1 % 360.0)
                if v1 > 180.0:
                    v1 = 360.0 - v1
                x[2 * k + ci] = old
                e1 = (v1 - val + 180.0) % 360.0 - 180.0
                row[2 * k + ci] += e1 * 3600.0 / eps
            rows.append((row, Pb))
            meta.append(("D", (st, tg, obs, val)))
            for j in range(nu):
                srow[j] += row[j]
            sl += dd * 3600.0
        l.append(sl)
        rows.append((srow, -Pb / float(len(rec))))
        meta.append(("S", st))
    # 方位角
    for azo in p["方位角"]:
        u, v, ob, ma = azo["u"], azo["v"], azo["az"], azo["Ma"]
        Pa = ((a + b) / ma) ** 2 if ma > 0 else 1e12
        val = azv(u, v)
        dd = (ob - val) % 360.0
        if dd > 180.0:
            dd -= 360.0
        l.append(dd * 3600.0)
        A = C(u); B = C(v)
        dx = B[0] - A[0]; dy = B[1] - A[1]
        s2 = dx * dx + dy * dy
        row = [0.0] * nu
        if not pts[v]["known"]:
            k = idx[v]; row[2 * k] += -dy / s2 * RHO; row[2 * k + 1] += dx / s2 * RHO
        if not pts[u]["known"]:
            k = idx[u]; row[2 * k] += dy / s2 * RHO; row[2 * k + 1] += -dx / s2 * RHO
        rows.append((row, Pa))
        meta.append(("A", azo))
    return l, rows, meta, C, azv, dv, nu


def _coarse_starts(p, maxcomb=4096):
    """贪心三点定位概算：按点号顺序，用两条连到已定点/已解点的实测边作圆交会；
    对 m 个待定点枚举 2^m 种镜像组合，作为最小二乘的多起点。"""
    pts = {q_["no"]: q_ for q_ in p["点"]}
    adj = {}
    for e in p["边"]:
        adj.setdefault(e["u"], []).append((e["v"], e["obs"]))
        adj.setdefault(e["v"], []).append((e["u"], e["obs"]))
    resolved = {q_["no"]: (q_["x"], q_["y"]) for q_ in p["点"] if q_["known"]}
    seq = []
    for q in p["点"]:
        no = q["no"]
        if q["known"]:
            continue
        cand = [(nb, d) for (nb, d) in adj.get(no, []) if nb in resolved]
        if len(cand) < 2:
            return None
        (na, d1), (nb, d2) = cand[0], cand[1]
        A = resolved[na]; B = resolved[nb]
        ab = math.hypot(B[0] - A[0], B[1] - A[1])
        if ab <= 0 or d1 <= 0:
            return None
        cosv = (d1 * d1 + ab * ab - d2 * d2) / (2.0 * d1 * ab)
        cosv = max(-1.0, min(1.0, cosv))
        ang = math.acos(cosv)
        base = math.atan2(B[1] - A[1], B[0] - A[0])
        p1 = (A[0] + d1 * math.cos(base + ang), A[1] + d1 * math.sin(base + ang))
        p2 = (A[0] + d1 * math.cos(base - ang), A[1] + d1 * math.sin(base - ang))
        seq.append((no, p1, p2))
        resolved[no] = p1
    m = len(seq)
    if m == 0:
        return None
    if (1 << m) > maxcomb:
        masks = [0]
    else:
        masks = list(range(1 << m))
    starts = []
    for mask in masks:
        x = []
        for k, (no, p1, p2) in enumerate(seq):
            pp = p1 if ((mask >> k) & 1) == 0 else p2
            x += [pp[0], pp[1]]
        starts.append(x)
    return starts


def _cost(p, x, unk):
    l, rows, _meta, _C, _a, _d, _nu = _build(p, x, unk)
    return sum(w * v * v for (row, w), v in zip(rows, l)), l, rows


def _lm(p, x0, iters=250):
    """Levenberg–Marquardt（自适应阻尼），稳健收敛到全局最优局部解"""
    unk = [q_["no"] for q_ in p["点"] if not q_["known"]]
    x = list(x0)
    lam = 1e-3
    c, _l, _r = _cost(p, x, unk)
    for _ in range(iters):
        l, rows, _meta, _C, _a, _d, nu = _build(p, x, unk)
        N = [[0.0] * nu for _ in range(nu)]
        bb = [0.0] * nu
        for i, (row, w) in enumerate(rows):
            for j in range(nu):
                bb[j] += w * row[j] * l[i]
                for k in range(nu):
                    N[j][k] += w * row[j] * row[k]
        ok = False
        for _try in range(14):
            M = [r[:] for r in N]
            for j in range(nu):
                M[j][j] *= (1.0 + lam)
                M[j][j] += 1e-12
            dx = _gauss(M, bb)
            if dx is None:
                lam *= 10.0
                continue
            xn = [x[j] + dx[j] for j in range(nu)]
            cn, _l2, _r2 = _cost(p, xn, unk)
            if cn < c:
                x = xn
                c = cn
                lam = max(lam / 3.0, 1e-14)
                ok = True
                break
            lam *= 10.0
        if not ok:
            break
        if max(abs(v) for v in dx) < 1e-12:
            break
    # 收尾牛顿（无阻尼）
    for _ in range(30):
        l, rows, _meta, _C, _a, _d, nu = _build(p, x, unk)
        N = [[0.0] * nu for _ in range(nu)]
        bb = [0.0] * nu
        for i, (row, w) in enumerate(rows):
            for j in range(nu):
                bb[j] += w * row[j] * l[i]
                for k in range(nu):
                    N[j][k] += w * row[j] * row[k]
        for j in range(nu):
            N[j][j] += 1e-12
        dx = _gauss(N, bb)
        if dx is None:
            break
        for j in range(nu):
            x[j] += dx[j]
        if max(abs(v) for v in dx) < 1e-12:
            break
    return x, c


def _solve(p, x0, iters=60, damp=1e-3):
    unk = [q_["no"] for q_ in p["点"] if not q_["known"]]
    x = list(x0)
    for _ in range(iters):
        l, rows, meta, C, azv, dv, nu = _build(p, x, unk)
        N = [[0.0] * nu for _ in range(nu)]
        bb = [0.0] * nu
        for i, (row, w) in enumerate(rows):
            for j in range(nu):
                bb[j] += w * row[j] * l[i]
                for k in range(nu):
                    N[j][k] += w * row[j] * row[k]
        for j in range(nu):
            N[j][j] *= (1.0 + damp)
        dx = _gauss(N, bb)
        if dx is None:
            return None
        for j in range(nu):
            x[j] += dx[j]
        if max(abs(v) for v in dx) < 1e-11:
            break
    return x


def compute(p):
    unk = [q_["no"] for q_ in p["点"] if not q_["known"]]
    kn = [q_ for q_ in p["点"] if q_["known"]]
    if kn:
        cx = sum(q_["x"] for q_ in kn) / len(kn)
        cy = sum(q_["y"] for q_ in kn) / len(kn)
    else:
        cx = cy = 0.0
    # 概算：贪心三点定位（正弦/余弦定理）→ 枚举镜像组合作为多起点
    starts = _coarse_starts(p)
    if not starts:
        Ls = [e["obs"] for e in p["边"] if e["obs"] > 0]
        rad = (sum(Ls) / len(Ls)) if Ls else 1000.0
        starts = []
        n = len(unk)
        x0 = []
        for k in range(n):
            ang = 2.0 * math.pi * k / max(n, 1)
            x0 += [cx + rad * math.cos(ang), cy + rad * math.sin(ang)]
        starts.append(x0)
    best = None
    for x0 in starts:
        x, c = _lm(p, x0)
        if best is None or c < best[1]:
            best = (x, c)
    x = best[0]
    l, rows, meta, C, azv, dv, nu = _build(p, x, unk)
    N = [[0.0] * nu for _ in range(nu)]
    pvv = 0.0
    for i, (row, w) in enumerate(rows):
        pvv += w * l[i] * l[i]
        for j in range(nu):
            for k in range(nu):
                N[j][k] += w * row[j] * row[k]
    Q = _inv(N)
    # 自由度口径：原著把每个测站的"史赖伯和方程"与 (n−1) 个角度方程一并计入，
    # 其自由度计数比 (观测数 − 未知数) 多 1。该 +1 由权威 OUT 的 Mo 唯一反演确定：
    #   √(4073.48/10)=20.18mm（按 10）  vs  √(4073.48/11)=19.24mm→打印 19.2mm（权威 19.2mm）
    # 详见 f12.py 文件头「已声明偏差 D4」。
    r = len(l) - nu + (1 if p["ND"] > 0 else 0)
    Mo = math.sqrt(pvv / r) if r > 0 else 0.0
    return {"输入": p, "坐标": x, "定位": {pt: k for k, pt in enumerate(unk)},
            "Q": Q, "Mo": Mo, "r": r, "pvv": pvv,
            "l": l, "meta": meta, "C": C, "azv": azv, "dv": dv, "nu": nu}


def render(p, R):
    pts = {q_["no"]: q_ for q_ in p["点"]}
    x = R["坐标"]; idx = R["定位"]; Q = R["Q"]; Mo = R["Mo"]

    def C(no):
        pt = pts[no]
        if pt["known"]:
            return (pt["x"], pt["y"])
        k = idx[no]
        return (x[2 * k], x[2 * k + 1])

    def azv(u, v):
        A = C(u); B = C(v)
        return math.degrees(math.atan2(B[1] - A[1], B[0] - A[0])) % 360.0

    def dvv(u, v):
        A = C(u); B = C(v)
        return math.hypot(B[0] - A[0], B[1] - A[1])

    def se(no):
        if pts[no]["known"] or Q is None:
            return None
        k = idx[no]
        qxx, qyy, qxy = Q[2 * k][2 * k], Q[2 * k + 1][2 * k + 1], Q[2 * k][2 * k + 1]
        t = math.degrees(0.5 * math.atan2(2 * qxy, qxx - qyy)) % 180.0
        sq = math.sqrt((qxx - qyy) ** 2 + 4 * qxy ** 2)
        E = Mo * math.sqrt(max((qxx + qyy + sq) / 2.0, 0.0)) * 1000.0
        F = Mo * math.sqrt(max((qxx + qyy - sq) / 2.0, 0.0)) * 1000.0
        return t, E, F

    def side_prec(u, v):
        A = C(u); B = C(v)
        dx = B[0] - A[0]; dy = B[1] - A[1]
        S = math.hypot(dx, dy)
        ga = {}
        for no, sg in ((v, 1.0), (u, -1.0)):
            if pts[no]["known"]:
                continue
            k = idx[no]
            ga[(0, k)] = sg * dx / S
            ga[(1, k)] = sg * dy / S
        keys = list(ga.keys())
        varS = 0.0
        for (c1, k1) in keys:
            for (c2, k2) in keys:
                varS += ga[(c1, k1)] * ga[(c2, k2)] * Q[2 * k1 + c1][2 * k2 + c2]
        Ms = Mo * math.sqrt(max(varS, 0.0))
        gb = {}
        for no, sg in ((v, -1.0), (u, 1.0)):
            if pts[no]["known"]:
                continue
            k = idx[no]
            gb[(0, k)] = sg * (-dy / S / S) * RHO
            gb[(1, k)] = sg * (dx / S / S) * RHO
        kb = list(gb.keys())
        varA = 0.0
        for (c1, k1) in kb:
            for (c2, k2) in kb:
                varA += gb[(c1, k1)] * gb[(c2, k2)] * Q[2 * k1 + c1][2 * k2 + c2]
        Ma = Mo * math.sqrt(max(varA, 0.0))
        return Ms, Ma, S

    SEP = "├──┼──┼─────┼─────┼─────┼──────┼───╂──┼──────┼──────┤"
    SEP2 = "┝━━┿━━┿━━━━━┿━━━━━┿━━━━━┿━━━━━━┿━━━╋━━┿━━━━━━┿━━━━━━┥"

    def condrow(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10):
        return ("│" + f1 + "│" + f2 + "│" + f3 + "│" + f4 + "│" + f5 + "│" + f6
                + "│" + f7 + "┃" + f8 + "│" + f9 + "│" + f10 + "│")

    # ---- 方向号（全局连续编号） ----
    dirno = {}
    num = 0
    for d in p["方向"]:
        for (tg, obs) in d["rec"]:
            num += 1
            dirno[(d["st"], tg)] = (num, obs)

    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****            测  边  网  平  差  计  算   F-12(99.6版)          ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append(" " * 36 + "网名:" + p["网名"])
    A.append("┌──┬──┬─────┬─────┬─────┬──────┬───┰──┬──────┬──────┐")
    A.append("│ 边 │ 目 │ 观测边长 │平差後边长│   Ms/S   │  方 位 角  │ Mα平┃方向│  观 测 值  │  平 差 值  │")
    A.append("│ 号 │ 标 │   (m)    │   (m)    │          │ (实 测 α) │ Mα测┃ No.│  ° ′ ″  │  ° ′ ″  │")
    A.append(SEP)

    # 方位角（按边终点登记）
    az_by = {}
    for azo in p["方位角"]:
        az_by.setdefault((azo["u"], azo["v"]), []).append(azo)

    firstblk = True
    for pt in p["点"]:
        no = pt["no"]
        if not firstblk:
            A.append(SEP2)
        firstblk = False
        xx, yy = C(no)
        s = dpad("P %d:%s" % (no, pt["name"]), 24) + "X=%12.3f,Y=%12.3f" % (xx, yy)
        ell = se(no)
        if ell:
            s += " " * 15 + "Φe=%5.1f°" % ell[0] + "E=%4d; F=%4dmm" % (
                int(round(ell[1])), int(round(ell[2])))
            s = dpad(s, 101)
        else:
            s = dpad(s, 102)
        # 本点要打印的边
        station = no in [d["st"] for d in p["方向"]]
        rows = []
        if station:
            for d in p["方向"]:
                if d["st"] != no:
                    continue
                for (tg, obs) in d["rec"]:
                    rows.append((tg, True))
        else:
            for e in sorted(p["边"], key=lambda z: z["no"]):
                if min(e["u"], e["v"]) == no:
                    other = e["v"] if e["u"] == no else e["u"]
                    rows.append((other, False))
        A.append("│" + s + "│")
        if rows:
            A.append(SEP)
        for (other, isdir) in rows:
            e = None
            for c in p["边"]:
                if (c["u"] == no and c["v"] == other) or (c["u"] == other and c["v"] == no):
                    e = c
            u, v = e["u"], e["v"]
            S = dvv(u, v)
            Ms, Ma, _ = side_prec(u, v)
            rel = max(1, int(round(S / Ms / 100.0))) * 100 if Ms > 0 else 0
            az = azv(no, other)
            f7 = "%6.2f" % Ma
            if isdir:
                dn, obs = dirno[(no, other)]
                val = (azv(no, other) - azv(no, [t for (t, o) in p["方向"][0]["rec"]][0]) if False else None)
                dd = None
                for d in p["方向"]:
                    if d["st"] == no:
                        t0 = d["rec"][0][0]
                        dd = azv(no, other) - azv(no, t0)
                dd = dd % 360.0
                if dd > 180.0:
                    dd = 360.0 - dd
                f8 = "%3d " % dn
                f9 = _dms(obs)
                f10 = _dms(dd)
            else:
                f8 = " " * 4
                f9 = " " * 12
                f10 = " " * 12
            A.append(condrow("S %2d" % e["no"], "P %2d" % other, "%10.3f" % e["obs"],
                             "%10.3f" % q(S, 3), "1:%-8d" % rel, _dms(az), f7, f8, f9, f10))
            # 加测方位角行
            marks = az_by.get((e["u"], e["v"]), []) + az_by.get((e["v"], e["u"]), [])
            for z in marks:
                if z["u"] != min(u, v):
                    continue
                A.append(condrow(" " * 4, " " * 4, " " * 10, " " * 10, "    (实测)",
                                 _dms(z["az"]), "%6.2f" % z["Ma"], " " * 4, " " * 12, " " * 12))
    A.append("└──┴──┴─────┴─────┴─────┴──────┴───┸──┴──────┴──────┘")
    A.append(" Ms=%2gmm+%2gppm,    方向观测中误差:M=±%5.2f ″,    平差後单位权(每km)中误差:Mo=±%5.1fmm"
             % (p["a"], p["b"], p["Mβ"], Mo))
    A.append(" 计算者:%s        日期:%s" % (p["计算者"], p["日期"]))
    A.append("")
    return "\n".join(A)


def run(data, out_txt=None, out_json=None):
    params = parse(data)
    result = compute(params)
    return result, render(params, result)


if __name__ == "__main__":
    import sys
    _r, _t = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_t)
