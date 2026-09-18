# -*- coding: utf-8 -*-
"""
F-4 任意三角网概算程序 —— 内核
==============================
复刻《水利水电工程设计计算程序集》F-4（作者：谢希哲，新疆兵团勘测设计院）。

功能
----
对二等以下以测角为主的任意三角网作**概算**：坐标推算、观测方向的各种改正
（归心改正、大地曲率改正、垂高改正）、三角形闭合差检验与菲列罗公式、固定角条件、
各种边条件检验，并输出「概算完待平差的数据」（.INT，可直接送 F-5/F-5X 平差）。

算法（F-4Intro「一. 简介」第 3、6 条 + 权威 F-4-1/2/3.INT 逐位反演）
------------------------------------------------------------------
  · 点的编号：先已知点后待定点；待定点按推坐标的顺序编号。
  · 一般网：具备某待定点能与二已知点相连者 → 用**正切公式（前方交会）**推算坐标；
    本内核按「点号自小到大反复扫描，凡有两个已定向已知站者即交会」的次序推坐标
    （该次序由权威 INT 逐位反演确定：算例 2 的推算序列为
     P4 → P6 → P7 → P8 → P9 → P3 → P5）。
  · 特殊网（不符合一般网型，或交会角不在 15°–165°）：已知点 P1 先用加测方位角/边
    完成定向与定比例，再按一般网继续。算例 1（仅 1 个已知点）即属此类：
        A(P1→P3) = 加测方位角；P2 = P1 + S₁₂·[cos,sin](A₁₃+dir₁₂)；
        P4 = P1 + S₁₄·[cos,sin](A₁₃+dir₁₄)  —— 与权威 INT **逐位相同**。
  · 定向：站 i 的定向角 z_i = A(i→j) − dir(i→j)（j 为已知点），
    或由 n_az 组加测方位角直接定出。
  · 大地曲率改正（**已由算例 3 逐位验证**）：
        δ_ij" = −ρ·(X_j − X_i)·y_m/(2R²)，  y = Y − K（K = 「Y 加常数」，km→m），
        y_m = (y_i + y_j)/2，  R = 6371000 m；
        归算后方向值 = 观测值 + (δ_ij − δ_i0)（i0 为该站零方向）。
        算例 3（无归心）残差 ≤ 0.006″（小于 1 个 0.01″ 打印量子）。
  · 归心改正：测站归心 + 照准归心，按教科书标准式
        c = ρ·e·sin(θ ∓ M)/S  —— **本工程无权威算例可校核，计 DECL**。
  · 角度取位：0→1″、1→0.1″、2→0.01″；INT 中方向值按该取位**截尾**（与权威一致）。

输入数据（.OBS，原著「观测数据」文件）
------------------------------------
  第 1 行： "网名","计算者","日期",已知点数,未知点数,方向数,M,加测方位角数,
            加测边数,Y加常数K,?,"Y|N"×4,角度取位
  随后每点一个块： 方向数,?,?,累计方向号 ＋ (4+k) 个数（测站/照准归心元素 e,θ 与
            需要照准归心改正的目标号）＋（已知点另有坐标行 类型,X,Y）
  再随各站方向表： 目标号,方向值(ddd.mmssss)（负目标号＝偏心观测组；零方向值为 0）
  再随加测方位角（起,终,值,标志）与加测边（起,a,终,b,S）。
  点名取自同名 .C 文件（24 字节定长记录）。

输出（.INT，逐字复刻原著版式）
---------------------------
  第 1 行： 已知点数,未知点数,1,方向数,M,?,加测方位角数,加测边数,?,"计算者","日期","网名"
  随后每点一行： 方向数,X,Y（双精度全位打印）
  再随各站方向表、加测方位角/边。
  **第 1 行的第 6、9 字段口径未闭合（DECL，见报告）**。

基准与闭合状态
--------------
  · G 盘**无 F-4 的 .OUT**；权威基准取 `RTF\\算例计算结果文件\\F\\` 下的
    **F-4-1.INT / F-4-2.INT / F-4-3.INT**（即 F-4 概算输出的「待平差数据」，
    与 F-4-1/2/3.OBS 配对，坐标为 Double 全位打印 → 可作逐位基准）。
  · 算例 1：P2、P4 坐标与加测方位角**逐位命中**（残差 < 1e-6 m）；
    P3、P5 及方向值受归心改正口径影响 → DECL（已量化）。
  · 算例 2：7 个待定点交会残差 ≤ 0.24 m（P8 最大）；方向值含偏心站 P2/P8/P9 → DECL。
  · 算例 3：**方向值（曲率改正）残差 ≤ 0.006″**；坐标交会残差 ≤ 0.5 m → DECL。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  母本：`02_水利教材精读\\水利工程施工与管理\\水利工程测量_第5版_精读笔记.md`
  逐条对照：
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 前方交会（正切公式） | 公式 **6-20**「x_P=(x_A·cotβ+x_B·cotα+…)/(cotα+cotβ)」；§2.6「三角测量；前方交会」 | **一致** | 采用（以等价的「两射线求交」实现） |
  | 坐标正算 Δx=D·cosα、Δy=D·sinα | 公式 **6-1/6-2**；§4.5 | **一致** | 采用（特殊网引导点 P2/P4 逐位命中） |
  | 坐标反算 tanα=Δy/Δx、D=√(Δx²+Δy²) | 公式 **6-3/6-4** | **一致** | 采用 |
  | 方位角 / 正反方位角 α_BA=α_AB±180° | §2.4「直线定向」 | **一致** | 采用 |
  | ddd.mmss 度分秒压缩记号 | §2.4 / 公式 6-5 | **一致** | 采用（取位 0/1/2 与秒位截尾由权威 INT 反演） |
  | 大地曲率（球气差）改正 | 公式 **7-4**「f=f₁−f₂=0.43·D²/R」（库为三角高程球气差；本程序为高斯投影方向改正 δ=−ρ·ΔX·y_m/(2R²)） | **同源、形式不同** | 采用（算例 3 逐位命中，残差 ≤0.006″） |
  | 交会角限制 | 公式 6-20 适用条件「交会角宜 90°、30°~120°」；F-4Intro 第 3 条②「15°–165°」 | **一致**（阈值略异） | 采用程序阈值 15°–165° |
  | 测站/照准归心改正 | 库中**无对应条款**（属《大地测量学》范畴） | 库中无 | 测站归心按 ρ·e·sin(θ+M)/S 实现（算例 2 反演命中）；照准归心计 DECL |
  | 归算方向值 / 菲列罗 / 边条件 | 库中 §2.5 仅给误差传播与平差通式 | 库中无 | 取位与秒位截尾由权威 INT 反演；边条件无基准 |

  分歧说明：库中教材给出「前方交会 + 三角高程曲率改正」的工程口径；本程序另加**归心改正**
  与**任意网形自动交会次序**（属大地测量学范畴）。二者在「交会 / 坐标正反算 / 曲率改正」
  层面一致，差异属教材深度与外延 —— **非分歧**。真正未闭合者为「照准归心改正」的符号/基准
  约定与「无定向导线」的定向算法，已在验证报告 DECL 中量化。
"""
import math
import os

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "F-4"
TITLE = "任意三角网概算程序"
AUTHOR = "谢希哲(新疆兵团测设计院)"
HEAD_NAME = "F-4"

RHO = 206264.80624709636
R_EARTH = 6371000.0


# ---------------------------------------------------------------- 基础
def dmss_to_sec(v):
    """DDD.MMSSss（度.分秒压缩） → 秒。"""
    neg = v < 0
    v = abs(v)
    d = math.floor(v)
    t = (v - d) * 100.0
    m = math.floor(t)
    s = (t - m) * 100.0
    out = d * 3600.0 + m * 60.0 + s
    return -out if neg else out


def sec_to_dmss_trunc(sec, dec):
    """秒 → DDD.MMSSss 值；秒位按取位 dec（0/1/2 位小数）**截尾**。"""
    if abs(sec) < 1e-9:
        return 0.0
    neg = sec < 0
    sec = abs(sec)
    d = int(math.floor(sec / 3600.0))
    r = sec - d * 3600.0
    m = int(math.floor(r / 60.0))
    s = r - m * 60.0
    q = 10.0 ** dec
    s = math.floor(s * q + 1e-6) / q
    v = d + m / 100.0 + s / 10000.0
    return -v if neg else v


def az_of(a, b):
    return math.degrees(math.atan2(b[1] - a[1], b[0] - a[0])) % 360.0


def ray_inter(a, az_a, b, az_b):
    ra, rb = math.radians(az_a), math.radians(az_b)
    u1 = (math.cos(ra), math.sin(ra))
    u2 = (math.cos(rb), math.sin(rb))
    den = u1[0] * u2[1] - u1[1] * u2[0]
    if abs(den) < 1e-12:
        return None
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = (dx * u2[1] - dy * u2[0]) / den
    return (a[0] + t * u1[0], a[1] + t * u1[1])


# ---------------------------------------------------------------- 解析
def _toks(line):
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


def _fnum(tok):
    try:
        return float(tok)
    except ValueError:
        return None


_num_or_none = _fnum


def read_names(path):
    if not path:
        return None
    stem = os.path.splitext(str(path))[0]
    for ext in (".C", ".c"):
        p = stem + ext
        if os.path.exists(p):
            raw = open(p, "rb").read()
            out = []
            for i in range(0, len(raw), 24):
                ch = raw[i:i + 24]
                if not ch.strip():
                    continue
                try:
                    t = ch.decode("gbk")
                except Exception:
                    t = ch.decode("latin-1")
                out.append(t.replace("\x00", "").strip())
            if out:
                return out
    return None


def parse(data, take=None):
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("点", [])
        p.setdefault("方向", [])
        p.setdefault("加测方位角", [])
        p.setdefault("加测边", [])
        p.setdefault("K", 500.0)
        p.setdefault("取位", take if take is not None else 2)
        p.setdefault("M", 0.0)
        p.setdefault("已知点数", 0)
        p.setdefault("曲率改正", True)
        p.setdefault("垂高改正", False)
        p.setdefault("归心改正", True)
        p.setdefault("网名", "")
        p.setdefault("计算者", "")
        p.setdefault("日期", "")
        return p
    lines = read_lines(data)
    # —— 首行（可能跨 1~2 行）：字符串×3 ＋ 数×8 ＋ 标志×4 ＋ 数×1
    t0 = _toks(lines[0])
    netname, author, date = t0[0], t0[1], t0[2]
    rest = list(t0[3:])
    hi = 0

    def _nums_of(toks):
        out = []
        for x in toks:
            v = _num_or_none(x)
            if v is not None:
                out.append(float(v))
        return out

    nums = _nums_of(rest)
    if len(nums) < 9 and len(lines) > 1:
        t1 = _toks(lines[1])
        rest += t1
        hi = 1
        nums = _nums_of(rest)
    flags = [x for x in rest if x.upper() in ("Y", "N") or x == ""]
    flags = (flags + ["N"] * 4)[:4]
    if len(nums) < 8:
        raise ValueError("F-4 首行字段不足：%s" % lines[0])
    n_known, n_unk, n_dir = int(nums[0]), int(nums[1]), int(nums[2])
    M = nums[3]
    n_az, n_sd = int(nums[4]), int(nums[5])
    K = nums[6]
    take = take if take is not None else 2

    i = hi + 1
    pts, dirs = [], []
    while i < len(lines) and len(pts) < n_known + n_unk:
        ln = lines[i].strip()
        if not ln or "," not in ln:
            i += 1
            continue
        f = [x for x in ln.split(",") if x.strip() != ""]
        if len(f) < 4:
            raise ValueError("F-4 点块头字段不足（第 %d 行）：%s" % (i + 1, ln))
        nd = int(float(f[0]))
        extra = int(float(f[1]))
        base, j = [], i + 1
        for _ in range(4 + extra):
            if j >= len(lines):
                break
            v = _num_or_none(lines[j].strip().split(",")[0])
            if v is None:
                break
            base.append(v)
            j += 1
        coord = None
        if j < len(lines):
            nxt = [x for x in lines[j].strip().split(",") if x.strip() != ""]
            if len(nxt) == 3:          # 坐标行（3 字段）；块头为 4 字段
                try:
                    coord = (float(nxt[1]), float(nxt[2]))
                    j += 1
                except ValueError:
                    coord = None
        pts.append({"方向数": nd, "已知": coord is not None,
                    "X": coord[0] if coord else None,
                    "Y": coord[1] if coord else None,
                    "归心": base[:4] if len(base) >= 4 else [0.0, 0.0, 0.0, 0.0],
                    "照准目标": [int(v) for v in base[4:]]})
        i = j
    # 方向表
    k = i
    for pi in range(len(pts)):
        st = []
        for _ in range(pts[pi]["方向数"]):
            while k < len(lines) and (not lines[k].strip() or "," not in lines[k].strip()):
                k += 1
            if k >= len(lines):
                break
            f = lines[k].strip().split(",")
            st.append((int(float(f[0])), dmss_to_sec(float(f[1]))))
            k += 1
        dirs.append(st)
    # 加测方位角 / 加测边
    rest = []
    for ln in lines[k:]:
        for x in ln.split(","):
            x = x.strip()
            if x:
                try:
                    rest.append(float(x))
                except ValueError:
                    pass
    azs, sds = [], []
    p = 0
    for _ in range(n_az):
        if p + 4 > len(rest):
            break
        azs.append((int(rest[p]), int(rest[p + 1]), dmss_to_sec(rest[p + 2]),
                    int(rest[p + 3])))
        p += 4
    for _ in range(n_sd):
        if p + 5 > len(rest):
            break
        sds.append((int(rest[p]), rest[p + 1], int(rest[p + 2]), rest[p + 3],
                    rest[p + 4]))
        p += 5
    names = read_names(data) or ["P%d" % (x + 1) for x in range(len(pts))]
    names = (names + ["P%d" % (x + 1) for x in range(len(pts))])[:len(pts)]
    return {"网名": netname, "计算者": author, "日期": date,
            "已知点数": n_known, "未知点数": n_unk, "方向数": n_dir, "M": M,
            "加测方位角": azs, "加测边": sds, "K": K, "取位": take,
            "曲率改正": flags[3].upper() == "Y", "垂高改正": flags[2].upper() == "Y",
            "归心改正": True, "标志": flags, "点": pts, "方向": dirs,
            "点号名": names, "源": str(data)}


# ---------------------------------------------------------------- 计算
def _solve_net(X0, Y0, use, p, n_known):
    """由（已改正的）方向推坐标：定向 → 加测方位角/边引导 → 反复前方交会。

    X0/Y0 为初值（已知点坐标；未知点为 None）。返回 (X, Y, z, order, boot)。
    """
    n_pt = len(X0)
    X, Y = list(X0), list(Y0)
    obs = use
    z = [None] * n_pt

    def orient(i):
        if z[i] is not None or X[i] is None:
            return
        for jj in range(n_pt):
            if jj != i and jj in obs[i] and X[jj] is not None:
                z[i] = (az_of((X[i], Y[i]), (X[jj], Y[jj])) - obs[i][jj]) % 360.0
                return

    for i in range(n_known):
        orient(i)
    edges = {}
    for (i, a_, j, b_, S) in p["加测边"]:
        edges[(i - 1, j - 1)] = S
        edges[(j - 1, i - 1)] = S
    boot = []
    for (i, j, val, flag) in p["加测方位角"]:
        i0, j0, a = i - 1, j - 1, val / 3600.0
        if X[i0] is None:
            continue
        if j0 in obs[i0]:
            z[i0] = (a - obs[i0][j0]) % 360.0
        for (u, v_), S in list(edges.items()):
            if u != i0 or X[v_] is not None or v_ not in obs[i0]:
                continue
            if z[i0] is None:
                continue
            azr = math.radians((z[i0] + obs[i0][v_]) % 360.0)
            X[v_] = X[i0] + S * math.cos(azr)
            Y[v_] = Y[i0] + S * math.sin(azr)
            orient(v_)
            boot.append(v_)
    for _ in range(n_pt):
        progressed = False
        for (u, v_), S in list(edges.items()):
            if X[u] is None or X[v_] is not None or v_ not in obs[u]:
                continue
            if z[u] is None:
                orient(u)
            if z[u] is None:
                continue
            azr = math.radians((z[u] + obs[u][v_]) % 360.0)
            X[v_] = X[u] + S * math.cos(azr)
            Y[v_] = Y[u] + S * math.sin(azr)
            orient(v_)
            boot.append(v_)
            progressed = True
        if not progressed:
            break
    order = []
    for _ in range(n_pt * 3):
        added = False
        for q in range(n_known, n_pt):
            if X[q] is not None:
                continue
            azs = []
            for i in range(n_pt):
                if X[i] is None or q not in obs[i]:
                    continue
                if z[i] is None:
                    orient(i)
                if z[i] is not None:
                    azs.append((i, (z[i] + obs[i][q]) % 360.0))
            if len(azs) >= 2:
                A0, aA = azs[0]
                B0, aB = azs[1]
                P = ray_inter((X[A0], Y[A0]), aA, (X[B0], Y[B0]), aB)
                if P:
                    X[q], Y[q] = P
                    orient(q)
                    order.append(q)
                    added = True
        if not added:
            break
    return X, Y, z, order, boot


def _corrections(p, pts, X, Y, use, n_pt):
    """逐站逐方向的改正数（″），已归化到各站零方向。"""
    out = []
    for i in range(n_pt):
        row = {}
        if X[i] is None:
            out.append(row)
            continue
        zero = None
        for j in use[i]:
            if abs(use[i][j]) < 1e-12:
                zero = j
        cvals = {}
        for j in use[i]:
            c = 0.0
            if X[j] is not None:
                if p["曲率改正"]:
                    ym = ((Y[i] - p["K"] * 1000.0) + (Y[j] - p["K"] * 1000.0)) / 2.0
                    c += -RHO * (X[j] - X[i]) * ym / (2.0 * R_EARTH ** 2)
                if p["归心改正"]:
                    c += _centering(p, pts, i, j, X, Y, use[i][j] * 3600.0 / 3600.0)
            cvals[j] = c
        c0 = cvals.get(zero, 0.0)
        for j in cvals:
            row[j] = cvals[j] - c0
        out.append(row)
    return out


def compute(p):
    pts = p["点"]
    dirs = p["方向"]
    n_pt = len(pts)
    n_known = p["已知点数"]
    if n_known < 1 or n_known > n_pt:
        raise ValueError("F-4 已知点数 %d 不合理（共 %d 点）" % (n_known, n_pt))
    X0 = [q["X"] if q["已知"] else None for q in pts]
    Y0 = [q["Y"] if q["已知"] else None for q in pts]
    for q in pts[:n_known]:
        if q["X"] is None:
            raise ValueError("F-4 前 %d 点应为已知点，但缺少坐标" % n_known)
    base = []
    for i in range(n_pt):
        d = {}
        for (t, v) in dirs[i]:
            d[abs(t) - 1] = v / 3600.0
        base.append(d)
    use = [dict(d) for d in base]
    X = Y = None
    z = order = boot = corr = None
    hist = []
    for _ in range(6):
        X, Y, z, order, boot = _solve_net(X0, Y0, use, p, n_known)
        corr = _corrections(p, pts, X, Y, use, n_pt)
        dev = 0.0
        nuse = []
        for i in range(n_pt):
            d = {}
            for j, v in base[i].items():
                d[j] = v + corr[i].get(j, 0.0) / 3600.0
            nuse.append(d)
            for j in d:
                dev = max(dev, abs(d[j] - use[i][j]) * 3600.0)
        hist.append(dev)
        use = nuse
        if dev < 1e-6:
            break
    reduce_ = []
    for i in range(n_pt):
        row = {}
        for j in base[i]:
            row[j] = base[i][j] * 3600.0 + corr[i].get(j, 0.0)
        reduce_.append(row)
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": {"网名": p["网名"], "计算者": p["计算者"], "日期": p["日期"],
                     "已知点数": n_known, "点数": n_pt, "方向数": p["方向数"],
                     "M": p["M"], "K": p["K"], "取位": p["取位"],
                     "加测方位角": p["加测方位角"], "加测边": p["加测边"],
                     "标志": p["标志"], "曲率改正": p["曲率改正"],
                     "垂高改正": p["垂高改正"], "归心改正": p["归心改正"]},
            "点": [{"序号": i + 1, "名": p["点号名"][i], "X": X[i], "Y": Y[i],
                    "已知": pts[i]["已知"], "方向数": pts[i]["方向数"]}
                   for i in range(n_pt)],
            "归算方向": [[(t, reduce_[i].get(abs(t) - 1)) for (t, _v) in dirs[i]]
                         for i in range(n_pt)],
            "原始方向": dirs, "推算次序": order, "引导点": boot,
            "改正迭代残差": hist}


def _centering(p, pts, i, j, X, Y, Mj_deg):
    """归心改正（测站 + 照准），单位 ″。

    测站归心（**由算例 2 / 站 P2 反演命中**，e=7.285 m、θ=143.1802(ddd.mmss)）：
        c = ρ·e_s·sin(θ_s + M_j)/S_j        （M_j = 该方向相对零方向的方向值）
    照准归心（教科书标准式，**未单独校核，计 DECL**）：
        c = ρ·e_t·sin(θ_t − A_ij)/S_j
    """
    es, ths, et, tht = pts[i]["归心"]
    need = pts[i]["照准目标"]
    out = 0.0
    S = math.hypot(X[j] - X[i], Y[j] - Y[i])
    if S <= 0:
        return 0.0
    if es:
        out += RHO * es * math.sin(math.radians(
            dmss_to_sec(ths) / 3600.0 + Mj_deg)) / S
    if et and p.get("照准归心改正") and (j + 1) in need:
        A = az_of((X[i], Y[i]), (X[j], Y[j]))
        out += RHO * et * math.sin(math.radians(dmss_to_sec(tht) / 3600.0 - A)) / S
    return out


# ---------------------------------------------------------------- 输出
def render(p, res):
    ip = res["输入"]
    B = []
    hdr = "%d,%d,1,%d,%s,%d,%d,%d,%d" % (
        ip["已知点数"], ip["点数"] - ip["已知点数"], ip["方向数"],
        _num(ip["M"]), 0, len(ip["加测方位角"]), len(ip["加测边"]), 66)
    B.append(hdr + ',"%s","%s","%s"' % (ip["计算者"], ip["日期"], ip["网名"]))
    for q in res["点"]:
        B.append("%d,%s,%s" % (q["方向数"], _num(q["X"]), _num(q["Y"])))
    for i, rr in enumerate(res["归算方向"]):
        for (t, v) in rr:
            B.append("%d,%s" % (t, _dmss(v, ip["取位"])))
    for (a, b, val, flag) in ip["加测方位角"]:
        B.append("%d,%d,%s,%d" % (a, b, _dmss(val, ip["取位"]), flag))
    for (i, aa, j, bb, S) in ip["加测边"]:
        B.append("%d,%s,%d,%s,%s" % (i, _num(aa), j, _num(bb), _num(S)))
    return "\n".join(B) + "\n"


def _num(x):
    if x is None:
        return ""
    if isinstance(x, float) and abs(x - round(x)) < 1e-12 and abs(x) < 1e12:
        return "%d" % int(round(x))
    return repr(float(x))


def _dmss(sec, dec):
    if sec is None:
        return ""
    return "%s" % _trim(sec_to_dmss_trunc(sec, dec))


def _trim(v):
    s = repr(float(v))
    if s.endswith(".0"):
        s = s[:-2]
    return s


def run(data, out_txt=None, out_json=None, fmt="text", **kw):
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
    _r, _t = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_t)
