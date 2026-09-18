# -*- coding: utf-8 -*-
"""
F-5 任意三角网平差计算程序 —— 内核（含 F-5X 大规模版共用实现）
==============================================================
复刻《水利水电工程设计计算程序集》F-5 / F-5X（作者：谢希哲，新疆兵团勘测设计院）。

功能
----
对以测角为主的任意三角网（导线网）作**方向观测间接平差**：解算各未知点平差坐标、
各方向改正数 V、平差后方向值与方位角、每条边的方位角中误差 Mα 与边长相对中误差、
各未知点的误差椭圆要素，并给出平差后单位权（方向观测）中误差 Mo。

模型（F-5Intro「一. 简介」第 4 条 + 权威 F-5.OUT 逐位反演）
----------------------------------------------------------
  观测方程（每个测站每个方向一条）：
        v_ij = (A_ij − z_i) − L_ij
      A_ij = 由（概略或平差后）坐标反算的 i→j 方位角（″）
      z_i  = 第 i 测站的定向未知数（″）
      L_ij = 归算后方向值（″，以测站零方向为参考，零方向 L=0）
  未知数：各未知点坐标 dx,dy（2/点）＋ 各测站定向 z（1/站）
  定权  ：以方向观测之权为 1（P=1）；加测方位角 Pα=1/Mα²、加测边
          Ps=(M/(a+S·b/1000))²/42545（F-5Intro 第 4 条①）；取固定值时权×1000。
  解算  ：N=AᵗPA → 平方根/高斯法求逆；m₀=±√([pvv]/(n−u))；Q=N⁻¹
  精度  ：m_α = m₀·√(gᵗQg)，g=∂A/∂(坐标)（″/m）；m_S = m₀·√(g_Sᵗ Q g_S)；
          椭圆：q=未知点 2×2 协因数块，E/F=m₀√((q₁+q₂±√((q₁−q₂)²+4q₁₂²))/2)

逐位反演证据（F-5.INT，n=21、u=12、r=9）
-----------------------------------------
  · Mo=√(59.6024/9)=2.5734″ → 打印 2.57 ✓；21 个 V 值**逐字相同** ✓
  · 平差后坐标 P4(4547997.117,14613066.679)、P5(4546901.791,14617148.139)、
    P6(4548248.841,14619238.044) 逐字相同 ✓；11 条边长逐字相同 ✓
  · Mα：2.1973→2.20、1.1958→1.20、1.5215→1.52、1.7787→1.78、1.6912→1.69、
    1.9976→2.00、1.7531→1.75、2.1010→2.10、2.3771→2.38 ✓
  · 边长相对中误差 1:N（N=S/m_S）四舍五入到百位：79300/146300/144000/79500/
    125700/123500/71900/101900/94900 ✓
  · 误差椭圆 P4(166.6°,44,36)、P5(54.2°,23,19)、P6(139.8°,33,18) ✓
  · **方位角秒位为「截尾」而非四舍五入**（关键裁决）：335°07′53.3896″→53.38、
    52°11′19.3680″→19.36、123°17′47.0189″→47.01 —— 逐位命中；
    若按四舍五入则得 53.39/19.37/47.02，与权威不符 → 已排除。

输入数据（.INT，与 F-5 程序一致）
--------------------------------
  第 1 行: "网名","计算者","日期",等级,已知点数,总方向数,M,加测方位角数,加测边数,
           角度取位,坐标取位
  随后为各点三元组： X, Y, 方向数          （点号 1..N，已知点在前）
  再随为各站方向表： 目标点号, 方向值(ddd.mmssss)   （每站 方向数 条，零方向值 0）
  再随为 加测方位角（n_az 组：起,终,值(ddd.mmssss),固定标志）
        与 加测边（n_sd 组：起,a,终,b,S）—— 本算例无（见 DECL）
  点名取自同名的 .C 文件（24 字节定长记录）；缺省时用 P1…Pn。

输出（逐字复刻原著 .OUT 版式，GBK；框线为全角）
--------------------------------------------
  首行「文件：…」为运行期路径回显，本内核不生成。
  表格每行内容按**显示宽度 106 列**补齐（全角字符按 2 列计）。
  F-5X 变体：未知点行「°」后多一个全角空格、末行补 6 个半角空格（F-5Xvb 版式差异）。

基准与闭合状态
--------------
  · 权威 OUT：F-5.OUT（5337 B，52 行）、F-5X.OUT（5350 B，52 行）—— 均自第 2 行起逐行对拍。
  · 加测方位角/边的解析口径（.INT 尾部字段）无权威算例，计 DECL。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  母本：`02_水利教材精读\\水利工程施工与管理\\水利工程测量_第5版_精读笔记.md`
  逐条对照：
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 间接观测（参数）平差 [pvv]=最小 | §2.5「平差：最小二乘原理（[vv]=最小）」；§4.3；§5.1 关联表 | **一致** | 采用 |
  | 方向观测为等权（P=1） | §2.5 + 公式 **5-21**「P_i=μ²/m_i²」（同一测回等精度 → 等权） | **一致** | 采用 |
  | 方位角/边长观测按中误差定权 | 公式 5-21；F-5Intro 第 4 条①给出 Pα=1/Mα²、Ps=(M/(a+S·b/1000))²/42545 | 库中给出定权通式，未给具体系数 | 采用说明书系数，计 DECL（无基准） |
  | 单位权中误差 m₀=±√([pvv]/(n−u)) | §2.5「单位权中误差」；§4.1 | **一致** | 采用 |
  | 精度评定 Mα=m₀√(gᵗQg)、M_S=m₀√(g_SᵗQg_S) | 公式 **5-12**（线性函数误差传播）、**5-14**（一般函数误差传播）；§2.5 | **一致** | 采用 N⁻¹ 协因数阵 |
  | 误差椭圆 E/F/Φe | 库中 §2.5 给出误差传播通式，未列椭圆要素公式 | 库中无对应条款 | 采用标准式 E,F=m₀√((q₁+q₂±√((q₁−q₂)²+4q₁₂²))/2)、Φ=½atan2(2q₁₂,q₁−q₂)，权威逐位命中 |
  | 坐标反算方位角 α=arctan(ΔY/ΔX)、S=√(ΔX²+ΔY²) | 公式 **6-3/6-4**「tanα_AB=Δy/Δx；D=√(Δx²+Δy²)」；§2.4 方位角定义 | **一致** | 采用（本例逐位命中） |
  | 测站定向未知数（零方向归化） | 库中 §4.5「导线测量内业」给出方位角推算与闭合差分配，未展开方向观测间接平差 | 原则**一致**（同源最小二乘） | 采用；属教材深度差异，非分歧 |
  | ddd.mmss 度分秒压缩记号 | §2.4/§3.3 公式 6-5 方位角推算；教材一律用度分秒 | **一致** | 采用 |

  分歧说明：库中教材以「导线（闭合/附合/支）内业 + 角度闭合差分配」为主，
  未展开任意三角网的方向观测间接平差；两者的数学内核（[pvv]=最小、P=μ²/m²、
  误差传播定律）完全同源，差异仅在网形与解算规模 —— 非分歧。
  另：方位角**秒位截尾**（而非四舍五入）为原著打印口径，由权威 OUT 逐位反演确定
  （53.3896→53.38、19.3680→19.36、47.0189→47.01），教材未涉及该打印细节。
"""
import math
import os
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "F-5"
TITLE = "任意三角网平差计算程序"
AUTHOR = "谢希哲(新疆兵团勘测设计院)"
HEAD_NAME = "F-5"

RHO = 206264.80624709636
D2S = 3600.0

# ---------------------------------------------------------------- 版式常量
SEP_TOP = "┌──┬──────┬──────┬───┬──────┬───────┬───┬──────┬──────┐"
HDR1 = ("│ 方 │            │ 归  算  後 │改正数│  平 差 後  │"
        "   平 差 後   │      │  平 差 後  │     Ｍ     │")
HDR2 = ("│ 向 │ 目  标  名 │            │      │            │"
        "  ( 实  测 )  │  Mα │ ( 实  测 ) │    ━━    │")
HDR3 = ("│ 号 │            │ 方  向  值 │  V   │ 方  向  值 │"
        "  方 位 角 α │      │  边    长  │     Ｓ     │")
SEP_HEAVY = "┝━━┷━━━━━━┷━━━━━━┷━━━┷━━━━━━┷━━━━━━━┷━━━┷━━━━━━┷━━━━━━┥"
SEP_LIGHT = "├──┬──────┬──────┬───┬──────┬───────┬───┬──────┬──────┤"
SEP_BOT = "└──┴──────┴──────┴───┴──────┴───────┴───┴──────┴──────┘"
TITLE_LINE = (" ****       任  意  三  角  网  平  差  计  算  F-5  "
              "(99.6版)       ****")
STAR = "*" * 71
COL_W = 106


def disp_width(s):
    return sum(2 if ord(c) > 127 else 1 for c in s)


def pad_disp(s, w):
    n = w - disp_width(s)
    return s + (" " * n if n > 0 else "")


def q(x, d):
    """十进制四舍五入（VB6 Format 口径）。"""
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


# ---------------------------------------------------------------- 角度
def dms_to_sec(v):
    """DDD.MMSSss（度.分秒） → 秒。"""
    neg = v < 0
    v = abs(v)
    d = math.floor(v)
    t = (v - d) * 100.0
    m = math.floor(t)
    s = (t - m) * 100.0
    out = d * 3600.0 + m * 60.0 + s
    return -out if neg else out


def sec_to_dms_trunc(sec):
    """秒 → (d, m, s)，秒位**截尾**至 0.01″（权威口径）。"""
    d = int(math.floor(sec / 3600.0))
    r = sec - d * 3600.0
    m = int(math.floor(r / 60.0))
    s = r - m * 60.0
    s = math.floor(s * 100.0 + 1e-6) / 100.0
    if s >= 60.0:            # 截尾不会上溢，但保留防御
        s -= 60.0
        m += 1
    return d, m, s


def fmt_dir(sec):
    d, m, s = sec_to_dms_trunc(sec)
    return "%3d %02d %05.2f" % (d, m, s)


def az_to_dms(deg):
    return fmt_dir(deg * 3600.0)


# ---------------------------------------------------------------- 输入
def _quote_split(line):
    out, cur, inq = [], [], False
    for c in line:
        if c == '"':
            inq = not inq
            continue
        if c == "," and not inq:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(c)
    out.append("".join(cur).strip())
    return out


def read_point_names(int_path):
    """读取同名 .C 文件（24 字节定长记录）中的点名。"""
    if not int_path:
        return None
    stem = os.path.splitext(str(int_path))[0]
    for ext in (".C", ".c"):
        p = stem + ext
        if os.path.exists(p):
            raw = open(p, "rb").read()
            names = []
            for i in range(0, len(raw), 24):
                chunk = raw[i:i + 24]
                if not chunk.strip():
                    continue
                try:
                    t = chunk.decode("gbk")
                except Exception:
                    t = chunk.decode("latin-1")
                names.append(t.replace("\x00", "").strip())
            return names or None
    return None


def parse(data):
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("网名", "")
        p.setdefault("计算者", "")
        p.setdefault("日期", "")
        p.setdefault("等级", 3)
        p.setdefault("已知点数", 0)
        p.setdefault("M", 0.0)
        p.setdefault("加测方位角", [])
        p.setdefault("加测边", [])
        p.setdefault("点", [])
        p.setdefault("方向", [])
        return p
    lines = read_lines(data)
    head = _quote_split(lines[0])
    name, author, date = head[0], head[1], head[2]
    grade = int(float(head[3]))
    n_known = int(float(head[4]))
    n_dir = int(float(head[5]))
    M = float(head[6])
    n_az = int(float(head[7]))
    n_sd = int(float(head[8]))
    vals = []
    for ln in lines[1:]:
        for t in ln.split(","):
            t = t.strip()
            if not t:
                continue
            try:
                vals.append(float(t))
            except ValueError:
                pass
    # 点数：累计方向数恰为 n_dir 的最短前缀
    pts, k, tot, n_pt = [], 0, 0, None
    while k + 3 <= len(vals):
        x, y, nd = vals[k], vals[k + 1], vals[k + 2]
        pts.append([x, y, int(nd)])
        tot += int(nd)
        k += 3
        if tot == n_dir and 3 * len(pts) + 2 * n_dir <= len(vals):
            n_pt = len(pts)
            break
    if n_pt is None:
        raise ValueError("F-5 无法由 INT 判定点数（累计方向数 %d ≠ 总方向数 %d）"
                         % (tot, n_dir))
    dirs = []
    for i in range(n_pt):
        st = []
        for _ in range(pts[i][2]):
            tgt = int(vals[k])
            v = dms_to_sec(vals[k + 1])
            k += 2
            st.append((tgt, v))
        dirs.append(st)
    azs, sds = [], []
    rest = vals[k:]
    j = 0
    for _ in range(n_az):
        if j + 4 > len(rest):
            break
        azs.append((int(rest[j]), int(rest[j + 1]), dms_to_sec(rest[j + 2]),
                    int(rest[j + 3])))
        j += 4
    for _ in range(n_sd):
        if j + 5 > len(rest):
            break
        sds.append((int(rest[j]), rest[j + 1], int(rest[j + 2]), rest[j + 3],
                    rest[j + 4]))
        j += 5
    names = read_point_names(data) or ["P%d" % (i + 1) for i in range(n_pt)]
    names = (names + ["P%d" % (i + 1) for i in range(n_pt)])[:n_pt]
    return {"网名": name, "计算者": author, "日期": date, "等级": grade,
            "已知点数": n_known, "方向数": n_dir, "M": M,
            "加测方位角": azs, "加测边": sds, "点": pts, "方向": dirs,
            "点号名": names, "源": str(data)}


# ---------------------------------------------------------------- 线性代数
def _inv_solve(N, w):
    """高斯—约当 [N|I|w] → [I|N⁻¹|x]（列归一化以消去量级差）。"""
    m = len(N)
    sc = [math.sqrt(N[a][a]) if N[a][a] > 0 else 1.0 for a in range(m)]
    A = [[N[a][b] / (sc[a] * sc[b]) for b in range(m)] for a in range(m)]
    b = [w[a] / sc[a] for a in range(m)]
    M = [A[i][:] + [1.0 if i == j else 0.0 for j in range(m)] + [b[i]]
         for i in range(m)]
    for c in range(m):
        piv = max(range(c, m), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-13:
            raise ValueError("F-5 法方程奇异（第 %d 列）——网形不足以解算" % (c + 1))
        M[c], M[piv] = M[piv], M[c]
        d = M[c][c]
        M[c] = [v / d for v in M[c]]
        for r in range(m):
            if r != c and M[r][c] != 0.0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    x = [M[i][2 * m] / sc[i] for i in range(m)]
    Q = [[M[i][m + j] / (sc[i] * sc[j]) for j in range(m)] for i in range(m)]
    return x, Q


def _quad(Q, g):
    n = len(g)
    return sum(g[a] * Q[a][b] * g[b] for a in range(n) for b in range(n))


# ---------------------------------------------------------------- 计算
def compute(p):
    pts = p["点"]
    dirs = p["方向"]
    n_pt = len(pts)
    n_known = p["已知点数"]
    if not (0 <= n_known < n_pt):
        raise ValueError("F-5 已知点数 %d 不合理（共 %d 点）" % (n_known, n_pt))
    unk = list(range(n_known, n_pt))
    st_idx = [i for i in range(n_pt) if pts[i][2] > 0]
    n_u, n_s = len(unk), len(st_idx)
    if n_u == 0:
        raise ValueError("F-5 无未知点（全部为已知点）")
    pos = {q_: j for j, q_ in enumerate(unk)}
    spos = {q_: 2 * n_u + j for j, q_ in enumerate(st_idx)}
    m = 2 * n_u + n_s
    X = [q_[0] for q_ in pts]
    Y = [q_[1] for q_ in pts]

    A, L, P, meta = [], [], [], []

    def grad_az(xa, ya, i, j, vec_len):
        """i→j 方位角对未知量的梯度（″/m）。"""
        dx, dy = xa[j] - xa[i], ya[j] - ya[i]
        S2 = dx * dx + dy * dy
        g = [0.0] * vec_len
        if i in pos:
            g[2 * pos[i]] += RHO * dy / S2
            g[2 * pos[i] + 1] += -RHO * dx / S2
        if j in pos:
            g[2 * pos[j]] += -RHO * dy / S2
            g[2 * pos[j] + 1] += RHO * dx / S2
        return g

    def grad_s(xa, ya, i, j, vec_len):
        dx, dy = xa[j] - xa[i], ya[j] - ya[i]
        S = math.hypot(dx, dy)
        g = [0.0] * vec_len
        if i in pos:
            g[2 * pos[i]] += -dx / S
            g[2 * pos[i] + 1] += -dy / S
        if j in pos:
            g[2 * pos[j]] += dx / S
            g[2 * pos[j] + 1] += dy / S
        return g

    # —— 方向观测
    for i in st_idx:
        for (tgt, obs) in dirs[i]:
            j = tgt - 1
            if not (0 <= j < n_pt):
                raise ValueError("F-5 目标点号越界：%d" % tgt)
            dx, dy = X[j] - X[i], Y[j] - Y[i]
            S2 = dx * dx + dy * dy
            az = math.degrees(math.atan2(dy, dx)) * D2S
            g = grad_az(X, Y, i, j, m)
            g[spos[i]] = 1.0
            A.append(g)
            ll = (az - obs + 648000.0) % 1296000.0 - 648000.0
            L.append(ll)
            P.append(1.0)
            meta.append(("dir", i, j, S2 ** 0.5))
    # —— 加测方位角（口径见 DECL）
    for (i, j, val, flag) in p["加测方位角"]:
        i, j = i - 1, j - 1
        dx, dy = X[j] - X[i], Y[j] - Y[i]
        az = math.degrees(math.atan2(dy, dx)) * D2S
        g = grad_az(X, Y, i, j, m)
        A.append(g)
        L.append((az - val + 648000.0) % 1296000.0 - 648000.0)
        P.append(1.0 * (1000.0 if flag else 1.0))
        meta.append(("az", i, j, math.hypot(dx, dy)))
    # —— 加测边（口径见 DECL）
    for (i, a, j, b, S_obs) in p["加测边"]:
        i, j = i - 1, j - 1
        dx, dy = X[j] - X[i], Y[j] - Y[i]
        S = math.hypot(dx, dy)
        g = grad_s(X, Y, i, j, m)
        A.append(g)
        L.append(S_obs - S)
        sd = (a + S_obs * b / 1000.0)
        P.append(((p["M"] / sd) ** 2 / 42545.0) if sd else 1.0)
        meta.append(("S", i, j, S))

    n = len(A)
    if n <= m:
        raise ValueError("F-5 观测数 %d ≤ 未知数 %d，无法平差" % (n, m))
    N = [[0.0] * m for _ in range(m)]
    w = [0.0] * m
    for r in range(n):
        for a_ in range(m):
            if A[r][a_] == 0.0:
                continue
            w[a_] -= P[r] * A[r][a_] * L[r]
            for b_ in range(a_, m):
                if A[r][b_] == 0.0:
                    continue
                N[a_][b_] += P[r] * A[r][a_] * A[r][b_]
    for a_ in range(m):
        for b_ in range(a_):
            N[a_][b_] = N[b_][a_]
    x, Q = _inv_solve(N, w)
    V = [sum(A[r][a_] * x[a_] for a_ in range(m) if A[r][a_] != 0.0) + L[r]
         for r in range(n)]
    pvv = sum(P[r] * V[r] * V[r] for r in range(n))
    r_red = n - m
    m0 = math.sqrt(pvv / r_red) if r_red > 0 else 0.0

    XA = X[:]
    YA = Y[:]
    for j, k in enumerate(unk):
        XA[k] += x[2 * j]
        YA[k] += x[2 * j + 1]

    # —— 逐方向成果
    rows = []
    seen = {}
    for k, (kind, i, j, S0) in enumerate(meta):
        if kind != "dir":
            continue
        dx, dy = XA[j] - XA[i], YA[j] - YA[i]
        S = math.hypot(dx, dy)
        az = math.degrees(math.atan2(dy, dx)) % 360.0
        ekey = (min(i, j), max(i, j))
        first = ekey not in seen
        seen[ekey] = k
        ga = grad_az(XA, YA, i, j, m)
        gs = grad_s(XA, YA, i, j, m)
        qa = _quad(Q, ga) if (i in pos or j in pos) else 0.0
        qs = _quad(Q, gs) if (i in pos or j in pos) else 0.0
        ma = m0 * math.sqrt(max(qa, 0.0))
        ms = m0 * math.sqrt(max(qs, 0.0))
        rel = (S / ms) if ms > 0 else None
        rows.append({"序号": k + 1, "站": i, "目标": j, "S": S, "az": az,
                     "V": V[k], "L": None, "首次": first, "Ma": ma,
                     "mS": ms, "相对N": rel,
                     "有用": (i in pos or j in pos)})
    # 归算後/平差後方向值（以零方向为准）
    idx = 0
    for i in st_idx:
        nd = len(dirs[i])
        v0 = None
        for t_, _o in dirs[i]:
            if _o == 0.0:
                v0 = V[idx]
        if v0 is None:
            v0 = V[idx]
        for t_, _o in dirs[i]:
            rows[idx]["L"] = _o
            rows[idx]["平差方向"] = _o + V[idx] - v0
            idx += 1
    # —— 误差椭圆
    ell = {}
    for j, k in enumerate(unk):
        qxx, qyy = Q[2 * j][2 * j], Q[2 * j + 1][2 * j + 1]
        qxy = Q[2 * j][2 * j + 1]
        t1 = qxx + qyy
        t2 = math.sqrt((qxx - qyy) ** 2 + 4.0 * qxy * qxy)
        E = m0 * math.sqrt(max((t1 + t2) / 2.0, 0.0))
        F = m0 * math.sqrt(max((t1 - t2) / 2.0, 0.0))
        phi = 0.5 * math.degrees(math.atan2(2.0 * qxy, qxx - qyy))
        if phi < 0:
            phi += 180.0
        if phi >= 180.0:
            phi -= 180.0
        ell[k] = {"Phi": phi, "E": E * 1000.0, "F": F * 1000.0,
                  "qxx": qxx, "qyy": qyy, "qxy": qxy}
    # 每点方向数归属：rows 已按站顺序排列
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": {"网名": p["网名"], "计算者": p["计算者"], "日期": p["日期"],
                     "等级": p["等级"], "已知点数": n_known, "点数": n_pt,
                     "方向数": n, "M": p["M"]},
            "点": [{"序号": i + 1, "名": p["点号名"][i], "X": XA[i], "Y": YA[i],
                    "已知": i < n_known, "方向数": pts[i][2]} for i in range(n_pt)],
            "方向表": rows, "椭圆": ell, "m0": m0, "Mo": m0, "pvv": pvv,
            "r": r_red, "未知数": m}


# ---------------------------------------------------------------- 输出
def render(p, res, variant="F-5"):
    fwsp = "\u3000" if variant == "F-5X" else " "
    B = ["", " " + STAR, TITLE_LINE, " " + STAR, "",
         " " * 35 + "网 名:" + res["输入"]["网名"]]
    B.append(SEP_TOP)
    B.append(HDR1)
    B.append(HDR2)
    B.append(HDR3)
    # 逐站
    rows = res["方向表"]
    k = 0
    n_pt = len(res["点"])
    for pi in range(n_pt):
        pt = res["点"][pi]
        B.append(SEP_HEAVY)
        pre = "P%3d:%s" % (pt["序号"], pt["名"])
        inner = pad_disp(pre, 40) + "X=%.3f,Y=%.3f" % (pt["X"], pt["Y"])
        if not pt["已知"]:
            e = res["椭圆"][pi]
            inner += "  " + "Φe=%5.1f°%sE=%4dmm,  F=%4dmm" % (
                q(e["Phi"], 1), fwsp, int(q(e["E"], 0)), int(q(e["F"], 0)))
        B.append("│" + pad_disp(inner, COL_W) + "│")
        if pt["方向数"] <= 0:
            continue
        B.append(SEP_LIGHT)
        for _ in range(pt["方向数"]):
            row = rows[k]
            k += 1
            seg = ["%4d" % row["序号"],
                   pad_disp(res["点"][row["目标"]]["名"], 12),
                   fmt_dir(row["L"]),
                   "%6.2f" % q(row["V"], 2),
                   fmt_dir(row["平差方向"]),
                   " " + az_to_dms(row["az"]) + " "]
            show = row["首次"] and row["有用"]
            seg.append("%6.2f" % q(row["Ma"], 2) if show else " " * 6)
            seg.append("%10.3f  " % q(row["S"], 3))
            if show and row["相对N"]:
                n100 = int(round(row["相对N"] / 100.0)) * 100
                seg.append(("1:%-10s" % n100))
            else:
                seg.append(" " * 12)
            B.append("│" + pad_disp("│".join(seg), COL_W) + "│")
    B.append(SEP_BOT)
    tail = (" M=± %4.2f ″  Mo=±%7.2f″    计算者: %s   日期: %s"
            % (res["输入"]["M"], res["Mo"], res["输入"]["计算者"],
               res["输入"]["日期"]))
    if variant == "F-5X":
        tail += " " * 6
    B.append(tail)
    return "\n".join(B) + "\n"


def run(data, out_txt=None, out_json=None, fmt="text", variant="F-5", **kw):
    params = parse(data)
    result = compute(params)
    text = render(params, result, variant)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    _r, _t = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_t)
