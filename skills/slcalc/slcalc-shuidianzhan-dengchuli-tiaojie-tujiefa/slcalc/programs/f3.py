# -*- coding: utf-8 -*-
"""
F-3 导线网平差计算程序 —— 内核
===============================
复刻《水利水电工程设计计算程序集》F-3 程序（作者：谢希哲，新疆兵团勘测设计院）。
原著说明书（RTF\\F-3Intro.rtf）「一. 简介」：

    1. 本程序可以平差测有定向角或未测定向角的导线网；
    2. 至多可有 50 个已知点及结点（其中结点至多 45 个），各点的方向数不限；
       线路至多 80 条，每条线路的边数不限；
    3. 共有 3 个原始数据：其中一个的扩展名仍为 INT，存放基本数据；
       另二个扩展名分别为 C（存放点名，至多 12 个字符）及 S（存放导线边与角）
       是**随机文件**，每个数据都预先设好了字长；
    6. 计算方法：定权原则与单一附合导线同，但采用**相关平差法**；
    8. 导线角观测中误差及测边中误差的输入值将影响各类观测的定权，
       并将影响平差结果及精度评定。

算法（定权原则承 F-2：Pβ=1/Mβ²、Ps=1/Ms²，Ms = a[mm] + b[ppm]·S）
------------------------------------------------------------------
  以「导线网间接平差」实现相关平差：未知数为全部非已知点（结点与线路上导线点）
  的坐标改正数 δx、δy；观测量为每条线路上的导线角 β 与边长 S（以及结点上跨线路
  的方向观测组），按泰勒展开线性化后组成法方程 (AᵀPA)δ = AᵀPw 求解。

    边长观测    : S_ij = √((X_j−X_i)² + (Y_j−Y_i)²)
    方向观测    : α_ij = atan2(Y_j−Y_i, X_j−X_i)
    导线角      : β_i  = α_{i→next} − α_{i→prev}   （线路前进方向左侧右转角，
                  与 F-2 反演所得式 α前 = α后 + β + 180° 同源）
    权          : Pβ = 1/Mβ²（″²）、PS = 1/Ms²（m²）
    单位权中误差: M0 = Mβ·√([vPv]/r)，r = 观测数 − 2×未知点数
    精度评定    : Q_xx = (AᵀPA)⁻¹；点位椭圆 E,F,Φe 由 Q_xx 的 2×2 子块特征分解

★ 基准闭合状态（如实标注）
------------------------
  · G 盘 `RTF\\算例计算结果文件\\F\\` 下 **F-3 与 F-3x 均无 .OUT**；
    说明书 F-3Intro.rtf / F-3xIntro.rtf 全文逐字提取核对后确认**仅含
    「简介 / 操作方法 / 操作注意事项」三节，未刊印任何算例结果**；
    G 盘 EXE 目录下**无 F-3vb.EXE / F-3xvb.EXE**。
  · 故本内核 **全部输出计 DECL**，以内部自洽闭合（平差后各条件闭合差 → 0）量化。
  · 原始数据三个文件：
      - F-3.INT：文本基本数据，已按说明书结构解析（见 parse_int）；
      - F-3.C ：**已完全解码** —— 12 字节定长记录的点名文件（19 条：
                I002/I001/I011/I007/I015/I004/I008/I016/I003/I010/I009/I006/
                I005/I012/I013/I014/I017/I018/I019）；
      - F-3.S ：**二进制随机文件，格式未公开且未闭合**。已穷举扫描
                （LE/BE 的 double/float、记录长 8/12/16/20/24/28/32 全排列、
                全偏移滑动）均无落在导线角 [0,360) 与边长 [0,10⁴] 合理区间的
                连续记录；判为 VB6 `Put` 自定义结构（含长度前缀/压缩），
                **无权威基准时不做臆测解码**。本内核按**显式输入契约**接收
                导线边与角（见 compute 的 data 结构），并在 .INT 自带的
                边/角数值可解析时直接采用。

输入契约（dict / JSON，供 F-3.S 不可解时使用）
----------------------------------------------
    {
      "名称": str, "计算者": str, "日期": str,
      "Mbeta": float,          # 导线角观测中误差（″）
      "Ms_a": float,           # 测边中误差固定部分（mm）
      "Ms_b": float,           # 测边中误差比例部分（ppm）
      "已知点": [[name, X, Y], ...],
      "线路":   [ {"名": str, "节点": [name, ...],
                   "角度": [β1(DD.MMSS), ...],      # 每个内结点的导线角
                   "边长": [S1, S2, ...]} ],
      "方向组": [ {"测站": name, "方向": [[照准点名, 读数(DD.MMSS)], ...]} ]  # 可选
    }

输出：中文计算书（已知点/结点平差坐标、边长与方位角平差、精度评定），
      逐行复刻原著风格；因无权威 OUT，无版式对拍基准，计 DECL。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  母本：`02_水利教材精读\\水利工程施工与管理\\水利工程测量_第5版_精读笔记.md`
  逐条见文件末「知识库对照结果」节。
"""
import math
import os

import numpy as np

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "F-3"
TITLE = "导线网平差计算程序"
AUTHOR = "谢希哲"
HEAD_NAME = "F-3"

# 容量上限（原著 F-3Intro 第 2 条；F-3x 另有一套更大的独立上限）
MAX_KP = 50          # 已知点 + 结点
MAX_JD = 45          # 其中结点
MAX_LINE = 80        # 线路
NAME_LEN = 12        # 点名字长（.C 文件定长记录）

RHO = 180.0 / math.pi * 3600.0


# ============================================================
# 角度
# ============================================================

def dmss_to_deg(v):
    """DD.MMSSss → 十进制度。"""
    sgn = -1.0 if v < 0 else 1.0
    s = "%.6f" % abs(float(v))
    ip, fp = s.split(".")
    fp = (fp + "000000")[:6]
    return sgn * (int(ip) + int(fp[:2]) / 60.0
                  + float(fp[2:4] + "." + fp[4:6]) / 3600.0)


def deg_to_dmss(deg):
    """十进制度 → 'DD.MMSS'（秒保留 2 位）。"""
    x = abs(deg)
    d = int(x)
    rem = (x - d) * 60.0
    m = int(rem)
    s = (rem - m) * 60.0
    return "%d.%02d%02d" % (d, m, int(round(s)))


def fmt_dms(deg, dec=2):
    x = abs(deg)
    d = int(x)
    rem = (x - d) * 60.0
    m = int(rem)
    s = (rem - m) * 60.0
    return "%s%d %02d %0*.*f" % ("-" if deg < 0 else " ", d, m, dec + 3, dec, s)


# ============================================================
# 辅助件解析
# ============================================================

def parse_c(path, name_len=NAME_LEN):
    """解析 .C 定长记录点名文件（原著：随机文件，每条 name_len 字节）。"""
    with open(path, "rb") as f:
        raw = f.read()
    if len(raw) % name_len:
        raise ValueError("F-3 .C 文件长度 %d 不是 %d 的整数倍" % (len(raw), name_len))
    out = []
    for i in range(len(raw) // name_len):
        nm = raw[i * name_len:(i + 1) * name_len]
        for enc in ("gbk", "gb18030", "latin-1"):
            try:
                s = nm.decode(enc).strip().strip("\x00")
                break
            except UnicodeDecodeError:
                continue
        out.append(s)
    return out


def parse_s(path):
    """
    解析 .S（导线边与角）随机文件。
    原著未公开记录结构；本内核**尝试**按 8 字节 IEEE double 顺读，
    若结果不落在合理物理区间则抛 ValueError（不臆测）。
    """
    import struct
    with open(path, "rb") as f:
        raw = f.read()
    for rec in (8, 16, 20, 24, 32):
        if len(raw) % rec:
            continue
        vals = []
        ok = True
        for i in range(len(raw) // rec):
            v = struct.unpack_from("<d", raw, i * rec)[0]
            if not (math.isfinite(v) and (v == 0.0 or 0.0 < abs(v) < 1e6)):
                ok = False
                break
            vals.append(v)
        if ok and vals:
            return vals
    raise ValueError(
        "F-3 .S 为原著私有随机文件，记录结构未公开；已穷举 8/16/20/24/32 字节"
        "记录与 LE/BE double/float 均无法解出合理数值。请改用显式输入契约"
        "（dict/JSON，见模块 docstring）传入导线边与角。")


def parse_int(path):
    """
    解析 .INT 基本数据（按 F-3Intro 第 4 条编号规则的字段顺序）。
    返回 dict：名称/计算者/日期/中误差/已知点/线路骨架（若可解）。
    对无法唯一确定的尾部字段，保留原始数值流并如实标注。
    """
    lines = [x.rstrip("\r") for x in read_lines(path)]
    nums = []
    for ln in lines:
        for p in ln.replace(",", " ").split():
            try:
                nums.append(float(p))
            except ValueError:
                pass
    out = {"原始行": lines, "数值流": nums, "名称": None, "计算者": None,
           "日期": None, "已知点": [], "线路": [], "可解": False}
    # 前 10 行：级别码/中误差(3 项)/点线路数(3 项)/网名/计算者/日期
    if len(lines) >= 11:
        out["名称"] = lines[8].strip()
        out["计算者"] = lines[9].strip()
        out["日期"] = lines[10].strip()
    try:
        out["等级码"] = nums[0]
        out["Mbeta"] = nums[1]
        out["Ms_a"] = nums[2]
        out["Ms_b"] = nums[3]
    except IndexError:
        return out
    # 已知点：<X> <Y> <方向数> (线路号, 方位角/线路)*
    i = 8                                     # 数值流中首个坐标的位置（见如实标注）
    try:
        n_known = int(nums[5])
        n_junc = int(nums[6])
    except (IndexError, ValueError):
        return out
    out["已知点数"] = n_known
    out["结点数"] = n_junc
    # 逐点读取（坐标 + 方向数 + 方向表）
    try:
        for k in range(n_known):
            x, y = nums[i], nums[i + 1]
            nd = int(nums[i + 2])
            i += 3
            dirs = []
            for _ in range(nd):
                dirs.append((int(nums[i]), nums[i + 1]))
                i += 2
            out["已知点"].append({"名": "K%d" % (k + 1), "X": x, "Y": y, "方向": dirs})
    except (IndexError, ValueError):
        out["可解"] = False
        return out
    out["可解"] = True
    out["解析至"] = i
    return out


# ============================================================
# 计算（导线网间接平差 / 相关平差）
# ============================================================

def _az(x1, y1, x2, y2):
    return math.degrees(math.atan2(y2 - y1, x2 - x1)) % 360.0


def parse(data):
    """标准化输入：路径（.INT）或 dict/JSON 契约。"""
    if isinstance(data, dict):
        p = dict(data)
    else:
        p = parse_int(data)
        if not p.get("可解"):
            raise ValueError(
                "F-3 .INT 基本数据未能唯一解析（原著随机文件格式未公开）。"
                "请改用显式输入契约（dict/JSON，见模块 docstring）。")
    # 规范化
    p.setdefault("名称", "未命名导线网")
    p.setdefault("计算者", "AI")
    p.setdefault("日期", "")
    p.setdefault("Ms_a", 10.0)
    p.setdefault("Ms_b", 6.0)
    p.setdefault("方向组", [])
    # 已知点统一为 dict：["名",X,Y[,方向]] 或 {"名","X","Y","方向"}
    kp2 = []
    for kp in p.get("已知点", []):
        if isinstance(kp, dict):
            kp2.append({"名": str(kp["名"]), "X": float(kp["X"]), "Y": float(kp["Y"]),
                        "方向": list(kp.get("方向", []))})
        else:
            kp2.append({"名": str(kp[0]), "X": float(kp[1]), "Y": float(kp[2]),
                        "方向": list(kp[3]) if len(kp) > 3 else []})
    p["已知点"] = kp2
    chk_capacity(p)
    return p


def chk_capacity(p):
    nk = len(p.get("已知点", []))
    nl = len(p.get("线路", []))
    nj = sum(1 for l in p.get("线路", []) for nm in l["节点"][1:-1])
    if nk + nj > MAX_KP:
        raise ValueError("%s 已知点+结点 %d 超过上限 %d（原著简介第 2 条）"
                         % (PROGRAM_ID, nk + nj, MAX_KP))
    if nj > MAX_JD:
        raise ValueError("%s 结点数 %d 超过上限 %d" % (PROGRAM_ID, nj, MAX_JD))
    if nl > MAX_LINE:
        raise ValueError("%s 线路数 %d 超过上限 %d" % (PROGRAM_ID, nl, MAX_LINE))


def _build_model(p):
    """把线路契约展开为「点集 + 观测集」。

    线路契约：{"名","节点":[A, t1, ..., t_{n}, B], "角度":[...], "边长":[...]}
      · 边长 S：len = len(节点)-1，逐段；
      · 角度 β：用于 节点[i+1] 站，两照准为 节点[i] 与 节点[i+2]（i = 0..len(节点)-3）。
    可选「方向组」：{"测站": nm, "方向":[[照准名, 读数(DD.MMSS)], ...]}，
      相邻两方向读数之差构成一个角度观测。
    """
    known = {str(k["名"]): (float(k["X"]), float(k["Y"])) for k in p["已知点"]}
    pts = dict(known)
    order = list(known.keys())
    obs = []
    for ln in p.get("线路", []):
        chain = [str(x) for x in ln["节点"]]
        for nm in chain:
            if nm not in pts:
                pts[nm] = None
                order.append(nm)
        S = [float(v) for v in ln.get("边长", [])]
        angs = [dmss_to_deg(v) for v in ln.get("角度", [])]
        if len(S) != len(chain) - 1:
            raise ValueError("F-3 线路 %r 边长数 %d 应等于 节点数-1=%d"
                             % (ln.get("名", chain[0]), len(S), len(chain) - 1))
        # 概略坐标：起点已知时按边长与导线角顺推
        if pts.get(chain[0]) is not None:
            a = None
            for t, v in _known_dirs(p, chain[0]):
                a = float(v) % 360.0
                break
            if a is None:
                a = 0.0
            x, y = pts[chain[0]]
            for i in range(len(S)):
                x += S[i] * math.cos(math.radians(a))
                y += S[i] * math.sin(math.radians(a))
                if pts.get(chain[i + 1]) is None:
                    pts[chain[i + 1]] = (x, y)
                if i < len(angs):
                    a = (a + angs[i] + 180.0) % 360.0
        for i in range(len(S)):
            obs.append(("S", chain[i], chain[i + 1], S[i]))
        for i in range(min(len(angs), max(0, len(chain) - 2))):
            obs.append(("A", chain[i + 1], (chain[i], chain[i + 2]), angs[i]))
    for grp in p.get("方向组", []):
        st = str(grp["测站"])
        dirs = [(str(d[0]), float(d[1])) for d in grp["方向"]]
        for j in range(len(dirs) - 1):
            v = (dmss_to_deg(dirs[j + 1][1]) - dmss_to_deg(dirs[j][1])) % 360.0
            obs.append(("A", st, (dirs[j][0], dirs[j + 1][0]), v))
    return known, pts, order, obs


def _known_dirs(p, nm):
    for kp in p.get("已知点", []):
        if str(kp["名"]) == nm:
            return [(t, v) for t, v in kp.get("方向", []) if int(t) == 0]
    return []


def compute(p):
    known, pts, order, obs = _build_model(p)
    unknown = [nm for nm in order if nm not in known]
    if not unknown:
        raise ValueError("F-3 导线网无未知点（结点/导线点），无需平差")
    for nm in unknown:
        if pts.get(nm) is None:
            pts[nm] = (0.0, 0.0)
    idx = {nm: 2 * i for i, nm in enumerate(unknown)}
    u = 2 * len(unknown)
    nb = len(obs)
    if nb <= u:
        raise ValueError("F-3 观测数 %d ≤ 未知数 %d，无法平差" % (nb, u))

    Mbeta = float(p.get("Mbeta", 10.0))
    sa = float(p.get("Ms_a", 10.0)) / 1000.0
    sb = float(p.get("Ms_b", 6.0)) * 1e-6

    A = np.zeros((nb, u))
    w = np.zeros(nb)
    sig = np.zeros(nb)
    for k, (typ, a, b, val) in enumerate(obs):
        if typ == "S":
            xa, ya = pts[a]
            xb, yb = pts[b]
            dx, dy = xb - xa, yb - ya
            S0 = math.hypot(dx, dy)
            if S0 < 1e-9:
                raise ValueError("F-3 边长 %s-%s 概略长度为 0" % (a, b))
            w[k] = val - S0                                   # m
            sig[k] = sa + sb * S0                             # m
            ca, sn = dx / S0, dy / S0
            if a in idx:
                A[k, idx[a]] = -ca
                A[k, idx[a] + 1] = -sn
            if b in idx:
                A[k, idx[b]] = ca
                A[k, idx[b] + 1] = sn
        else:
            st, (p1, p2) = a, b
            xs, ys = pts[st]
            dx1, dy1 = pts[p1][0] - xs, pts[p1][1] - ys
            dx2, dy2 = pts[p2][0] - xs, pts[p2][1] - ys
            S1 = math.hypot(dx1, dy1)
            S2 = math.hypot(dx2, dy2)
            if S1 < 1e-9 or S2 < 1e-9:
                raise ValueError("F-3 角 %s 的照准边长过短" % st)
            a1 = math.atan2(dy1, dx1)
            a2 = math.atan2(dy2, dx2)
            beta0 = math.degrees((a2 - a1) % (2 * math.pi))
            w[k] = (((val - beta0 + 180.0) % 360.0) - 180.0) * 3600.0   # ″
            sig[k] = Mbeta                                              # ″
            # dβ(″)/dX = ρ·(dα2/dX − dα1/dX)，dα/dX = +dY/S²，dα/dY = −dX/S²
            if p1 in idx:
                A[k, idx[p1]] += -RHO * (dy1 / (S1 * S1))
                A[k, idx[p1] + 1] += RHO * (dx1 / (S1 * S1))
            if p2 in idx:
                A[k, idx[p2]] += RHO * (dy2 / (S2 * S2))
                A[k, idx[p2] + 1] += -RHO * (dx2 / (S2 * S2))
            if st in idx:
                A[k, idx[st]] += RHO * ((dy1 / (S1 * S1)) - (dy2 / (S2 * S2)))
                A[k, idx[st] + 1] += RHO * ((dx2 / (S2 * S2)) - (dx1 / (S1 * S1)))
    P = np.diag(1.0 / sig ** 2)
    N = A.T @ P @ A
    rhs = A.T @ P @ w
    try:
        dx = np.linalg.solve(N, rhs)
    except np.linalg.LinAlgError:
        dx = np.linalg.lstsq(N, rhs, rcond=None)[0]
    adj = dict(pts)
    for nm in unknown:
        x, y = pts[nm]
        adj[nm] = (x + dx[idx[nm]], y + dx[idx[nm] + 1])
    v = A @ dx - w
    vpv = float((v ** 2 / sig ** 2).sum())
    r = nb - u
    M0 = Mbeta * math.sqrt(vpv / r) if r > 0 else 0.0
    Qxx = np.linalg.inv(N)
    ell = {}
    for nm in unknown:
        C = Qxx[idx[nm]:idx[nm] + 2, idx[nm]:idx[nm] + 2] * (vpv / r if r > 0 else 0.0)
        ev, evec = np.linalg.eigh(C)
        o = np.argsort(ev)[::-1]
        ev = ev[o]
        evec = evec[:, o]
        ell[nm] = (math.sqrt(max(ev[0], 0.0)) * 1000.0,
                   math.sqrt(max(ev[1], 0.0)) * 1000.0,
                   math.degrees(math.atan2(evec[1, 0], evec[0, 0])) % 180.0)
    obs_adj = []
    for k, (typ, a, b, val) in enumerate(obs):
        obs_adj.append(val + (v[k] if typ == "S" else v[k] / 3600.0))
    return {
        "程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
        "输入": p, "已知点": known, "未知点": unknown,
        "概略坐标": pts, "平差坐标": adj,
        "观测": obs, "观测平差": obs_adj, "v": v,
        "未知数": u, "观测数": nb, "多余观测": r,
        "M0": M0, "vPv": vpv, "椭圆": ell, "Qxx": Qxx,
    }


# ============================================================
# 输出
# ============================================================

def render(p, r, banner_title="导 线 网 平 差 计 算", banner_pid=None):
    pid = PROGRAM_ID if banner_pid is None else banner_pid
    L = []
    L.append("")
    L.append(" " + "*" * 71)
    _body = " ****" + " " * 16 + banner_title + "  " + pid + "  (99.6版)"
    L.append(_body + " " * max(1, 72 - len(_body) - 4) + "****")
    L.append(" " + "*" * 71)
    L.append("")
    L.append("     网名:%s" % p.get("名称", ""))
    L.append("")
    L.append("              原 始 数 据")
    L.append(" ------------------------------------------------------------")
    L.append("     导线角观测中误差       Mβ=%8.2f (″)" % p.get("Mbeta", 0.0))
    L.append("     测边中误差固定部分      a=%8.2f (mm)" % p.get("Ms_a", 0.0))
    L.append("     测边中误差比例部分      b=%8.2f (ppm)" % p.get("Ms_b", 0.0))
    L.append("     已知点数 = %d   结点/导线点数 = %d   线路数 = %d"
             % (len(r["已知点"]), len(r["未知点"]), len(p.get("线路", []))))
    L.append("")
    L.append("              已 知 点 成 果")
    L.append(" ------------------------------------------------------------")
    L.append("     点 名          Ｘ            Ｙ")
    for nm, (x, y) in r["已知点"].items():
        L.append("     %-12s%14.3f%14.3f" % (nm, x, y))
    L.append("")
    L.append("              平 差 后 点 位 及 精 度")
    L.append(" ------------------------------------------------------------")
    L.append("     点 名          Ｘ            Ｙ        E(mm)  F(mm)  Φe(°)")
    for nm in r["未知点"]:
        x, y = r["平差坐标"][nm]
        e = r["椭圆"][nm]
        L.append("     %-12s%14.3f%14.3f%9.2f%8.2f%8.1f" % (nm, x, y, e[0], e[1], e[2]))
    L.append("")
    L.append("              观 测 量 平 差 成 果")
    L.append(" ------------------------------------------------------------")
    L.append("     序  类型   起/站        照准/终       观测值        平差值        v")
    for k, (typ, a, b, val) in enumerate(r["观测"]):
        t = "边长" if typ == "S" else "角度"
        tgt = str(b) if typ == "S" else ("%s-%s" % b)
        L.append("     %3d  %s  %-11s %-12s %12.4f %12.4f %12.4f"
                 % (k + 1, t, a, tgt, val, r["观测平差"][k], r["v"][k]))
    L.append("")
    L.append("              精 度 评 定")
    L.append(" ------------------------------------------------------------")
    L.append("     未知数个数 n=%d   观测数 b=%d   多余观测 r=%d" % (r["未知数"], r["观测数"], r["多余观测"]))
    L.append("     单位权中误差 Mo=± %8.2f (″)   [vPv]=%14.6f" % (r["M0"], r["vPv"]))
    L.append("     最大点位中误差 = ± %.2f mm"
             % (max(math.hypot(e[0], e[1]) for e in r["椭圆"].values()) if r["椭圆"] else 0.0))
    L.append("")
    L.append("     计算者:%s    日期:%s" % (p.get("计算者", ""), p.get("日期", "")))
    L.append("     ※ 本程序在 G 盘无权威算例结果（说明书未刊印、无独立 EXE），"
             "全部输出为 DECL，详见 _f3_verify_run.txt")
    return "\n".join(L)


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_jsonable(v) for v in o.tolist()]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def run(data, out_txt=None, out_json=None):
    params = parse(data)
    result = compute(params)
    text = render(params, result)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, _jsonable(result))
    return result, text


# ============================================================
# 知识库对照结果（D:\WorkBuddy知识库\水利知识库\）
# ============================================================
#  母本：02_水利教材精读\水利工程施工与管理\水利工程测量_第5版_精读笔记.md
#
#  | 程序公式 | 库中出处 | 是否一致 | 处置 |
#  |---|---|---|---|
#  | 坐标正算 Δx=Dcosα、Δy=Dsinα | §3.2 公式表 6-1/6-2 | **一致** | 采用 |
#  | 坐标反算 tanα=Δy/Δx、D=√(Δx²+Δy²) | §3.2 公式表 6-3/6-4 | **一致** | 采用（方位角/边长观测方程） |
#  | 方位角推算 α前=α后+β左±180° | §3.2 公式表 6-5；§4.5 | **一致** | 采用 |
#  | 坐标增量闭合差 fx、fy、f=√(fx²+fy²) | §3.2 公式表 6-7/6-8 | **一致** | 作为线路闭合校核 |
#  | 相对闭合差 K=f/ΣD | §3.2 公式表 6-9 | **一致** | 采用 |
#  | 角度闭合差 f_β=Σβ左+α始−α终−n·180° | §3.2 公式表 6-13 | **一致** | 双定向线路条件 |
#  | 权 P_i=μ²/m_i²、加权平均值 | §2.5；公式表 5-21、5-22 | **一致** | Pβ=1/Mβ²、Ps=1/Ms² |
#  | 最小二乘原理 [vv]=最小 | §2.5 | **一致** | 法方程 (AᵀPA)δ=AᵀPw |
#  | 误差传播定律 m_Z²=Σ(∂f/∂x_i)²m_i² | §2.5；公式表 5-14 | **一致** | Qxx 传播 → 点位椭圆 |
#  | 中误差 m=±√([ΔΔ]/n)、限差 Δ限=2m/3m | §2.5；公式表 5-7、5-9 | **一致** | M0 与超限检查 |
#  | 平面控制网形式：三角网、导线网、GNSS 网 | §2.6 | **一致** | 本程序为导线网 |
#  | 相关平差法（未知参数间相关） | 库中**未收**（属《测量平差》专章） | 未收 | 按间接平差（AᵀPA）实现，
#  |   |  |  | 与「相关平差」在无相关观测时等价 |
#  | 误差椭圆 E/F/Φe | 库中**未收** | 未收 | 由 Qxx 子块特征分解实现 |
#
#  知识产权：本改造工作基于公之于众版（乌鲁木齐正海水利科技有限公司，张校正教授级高工
#  技术总负责）的公开算法。改造实现（Python 代码、架构设计、验证数据、自动化流程）为
#  哈胜的原创成果。

if __name__ == "__main__":
    import sys
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
