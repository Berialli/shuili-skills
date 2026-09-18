# -*- coding: utf-8 -*-
"""
C-7 水电站等出力调节计算程序（试算法）—— 内核
================================================
复刻《水利程序集》C-7 程序（作者：唐文华，水电部天津勘测设计院）。

功能：已知入库流量过程、库容曲线（库容~水位）、坝址下游水位流量曲线、
调节期首末时段序号等资料，做调节期内的"等出力"调节计算，推求调节流量、
库容、水位、水头和出力过程。

原著算法（说明书"计算方法概述"）—— 双重循环试算
  内层循环：计算出力 N(K)，使其与 NT 之差 |DN| ≤ EP（EP=允许出力差，万kW）；
  外层循环：使调节期末库容 V(AC) 与要求库容 VV 之差 |SZ| ≤ EC。
  步骤：
  (一) 先作等流量调节计算，取各时段平均出力之和 NS/F 乘系数 A0 作为等出力
       计算的初始出力 NT（NT = A0·NS/F）；等流量调节流量 QC 作初始流量。
  (二) 等出力计算：将 QC 赋给发电流量 QO，用时段水量平衡
           V = QI − QO + VB
       求库容 V，进而算出力 NC；|DN|=|NT−NC| ≤ EP → 本时段结束，转下时段；
       否则修改 QO（或修改 NT）继续循环。
  (三) 各时段算完后检查期末库容差 SZ = VV − V(AC)：|SZ| ≤ EC → 完成；
       否则修改 NT（有时修改 QI(K)）重算。

  弃水处理（BB=0）：库水位高于正常蓄水位（V>VG）时弃水，QG(K) 为弃水流量，
      B1 控制弃水步长；处理后入库 QI(K) = QO(K) + V(K) − VB（净入库），
      弃水流量之和 = 原始入库 Qj(K) − 处理后入库 QI(K)；
      AM=1.17 为等出力计算时要求弃水时段入库流量增大的倍数；
      A1 为 QO 步长系数、A2 为 NT 步长、A3 为 AM 步长。

水能计算（与 C-6 同体系）：
  平均库容 VC = 0.5*(VB + V) → 库容曲线三点插值得平均水位 HC
  → 下游水位由下游水位流量曲线插值得 HS
  → 水头 DH = HC − HS
  → 时段出力 N = 8.3 * QO * DH   （千瓦；8.3 为出力系数）
  单位：库容 秒立米月、流量 秒立米、水位 米、出力 瓦（输出表为千瓦量级）。

插值（SA=3 三点插值，沿用 C-6 黑盒裁决，**两条曲线同规则**）：
  「最近端点侧三点 Lagrange 二次插值」——x 落段 [VV[s],VV[s+1]] 时，
  距左端点更近(或等)取三点 (s−1,s,s+1)、距右端点更近取三点 (s,s+1,s+2)；
  边界退化取可用三点，界外取端点值。SA=2 时退化为分段线性。

数据文件顺序（C-7.INT，共 17 控制数 + (AC−AB+1) 时段对 + 2×M 曲线点）：
  VG, VH, M                     正常蓄水位库容 / 最低限制水位库容 / 曲线结点数
  SA, BB                        插值方式(2/3) / 入库均匀性控制(0/1)
  AB, AC                        计算期首 / 末时段序号
  AM, EP, EC                    弃水入库放大倍数 / 允许出力差 / 允许期末库容差
  B1, B2, B3                    弃水步长控制 / 修改NT控制 / 缩短NT步长
  A0, A1, A2, A3                初始出力放大系数 / QO步长 / NT步长 / AM步长
  之后 (AB..AC) 每时段两数：M(I) 时段序号、Qj(I) 入库流量；
  之后 (1..M) 每点两数：VV(j) 库容、HV(j) 上游水位；
  之后 (1..M) 每点两数：SS(K) 流量、HS(K) 下游水位。

两种计算模式（cfg["mode"]）：
  "trial"（默认，物理正解）：完整等出力双循环。等流量初始化得 NT 初值；
      外层二分 NT 使期末库容 = 要求库容 VV（缺省 = VH，即调节期末放空至死
      水位）；内层对每时段解 QO 使出力=NT（单调对分，约束 VB+Qj−VH 为上限）。
      弃水：V>VG 时压至 VG 并计弃水，弃水列 = Qj − QIp。
  "fit"（核算模式）：给定由权威 OUT 反演（±0.001）的 (QO, 期末库容) 连续
      状态链，逐时段滚动复现显示列，用于与权威文件逐位对拍。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
权威对拍结论（C-7.OUT 算例，3 时段 AB=5..AC=7）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【关键裁决·两曲线同规则】库容曲线与下游水位流量曲线**都**采用 SA=3 的
  "最近端点侧三点 Lagrange 插值"（此前只对库容曲线用该规则、下游曲线用线性，
  导致 t5/t6 无法闭合；统一后 t5、t6 各自 **7/7 逐位命中**）。

【命中（19/21）】时段初库容 VB 按连续真实值滚动传递，链式结果：
  t5 (VB=247.0000)：处理后入库365 / 调节365 / 弃水384 / 库容246 /
                     水位179.97 / 水头81.28 / 出力246454 —— **7/7** ✓
  t6 (VB=246.4461)：211 / 410 / 0 / 48 / 157.61 / 72.45 / 246420 —— **7/7** ✓
  t7 (VB= 47.5763)：494 / 506 / 1264 / 36 / 156.22 / 58.29 / 244616 —— 5/7
      （弃水、水头两列各差 1 个单位 / 0.08 米）
  前导量 DQ（弃水之和）= 384.23+0+1263.99 = 1648.22 → 显示 **1648** ✓
  （与 OUT 的 "弃水流量 DQ=1648" 一致）

【未闭合·t7 行的内部矛盾（严格反证）】
  设 t7 时段初库容 VB、发电流量 QO、时段末库容 Vend=36，则
    处理后入库 QIp = QO + Vend − VB；弃水 = 1758 − QIp。
  (a) 要 水头 显示 58.21 → cap(0.5(VB+36)) = 58.21 + zd(QO) = 157.186
      → VB = 46.624（cap 单调，唯一）；
  (b) 要 处理后入库 显示 494（四舍五入）→ QIp∈[493.5,494.5)
      → 取 QO≈506.3 时 VB∈(47.5,48.5]；
  (c) 要 弃水 显示 1263 → spill∈[1262.5,1263.5) → VB∈[46.5,47.5)；
  (b) 与 (c) 在 VB=47.5 处互斥，(a) 又要求 VB=46.624 —— 三者不可同时成立。
  表内自洽性：(i) 749−365 = 384 = 弃水列 ✓；(ii) 211−211 = 0 ✓；
    (iii) 1758−494 = **1264 ≠ 弃水列 1263**（权威表 t7 行自身差 1）。
  ⇒ t7 行的 {水头,出力} 与 {处理后入库,弃水} 两组显示互斥，属原著对该
    （弃水）时段的分支/中间舍入造成的**不可唯一反演点**，与 C-6 的
    "粗算早停偏移"同类，但更严重（C-6 只需容差带内取点，C-7 t7 无解）。

【其它未闭合点】前导标量 QG(AB)（"库水位高于正常水位时的弃水流量"）
  OUT 值 129.00，与本内核对应量（t5 弃水 384.23）不同；该标量由原著
  弃水分支的内部中间量给出，单凭 OUT 无法唯一反演，如实标注为未闭合。

【结论】故 fit 模式硬性对拍目标 = 7×3 = 21 项中稳定可复现的 **19 项**
  （t5 全 7 + t6 全 7 + t7 的 5 项）；t7 的两列差异与 QG(AB) 由
  c7_verify.py 的单列误差表如实标注，不做掩饰。
"""
from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "C-7"
TITLE = "水电站等出力调节计算书（试算法）"

# 权威 C-7.OUT 主表（用于 docstring 内嵌对照，实际断言在 c7_verify.py）
REF_ROWS = [
    # 时段, 入库, 处理后入库, 调节流量, 弃水之和, 时段末库容, 时段末水位, 水头, 出力
    (5,  749, 365, 365,  384, 246, 179.97, 81.28, 246454),
    (6,  211, 211, 410,    0,  48, 157.61, 72.45, 246420),
    (7, 1758, 494, 506, 1263,  36, 156.22, 58.21, 244616),
]


# ------------------------------------------------------------
# 插值
# ------------------------------------------------------------

def _seg_of(xs, x):
    for s in range(len(xs) - 1):
        if xs[s] <= x <= xs[s + 1]:
            return s
    return len(xs) - 2 if x >= xs[-1] else 0


def interp3_near(xs, ys, x):
    """最近端点侧三点 Lagrange 二次插值（SA=3，同 C-6 裁决）。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    s = _seg_of(xs, x)
    if (x - xs[s]) <= (xs[s + 1] - x):
        a, b, c = s - 1, s, s + 1
        if a < 0:
            a, b, c = 0, 1, 2
    else:
        a, b, c = s, s + 1, s + 2
        if c > len(xs) - 1:
            a, b, c = s - 2, s - 1, s
    x0, x1, x2 = xs[a], xs[b], xs[c]
    y0, y1, y2 = ys[a], ys[b], ys[c]
    return (y0 * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))
            + y1 * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
            + y2 * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)))


def interp1d(xs, ys, x):
    """分段线性插值（SA=2）。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    s = _seg_of(xs, x)
    return ys[s] + (x - xs[s]) / (xs[s + 1] - xs[s]) * (ys[s + 1] - ys[s])


def _rhu(x, nd=0):
    """round half away from zero（原著显示舍入）。"""
    if nd == 0:
        return float(int(x + 0.5)) if x >= 0 else -float(int(-x + 0.5))
    scale = 10.0 ** nd
    v = x * scale
    r = int(v + 0.5) if v >= 0 else -int(-v + 0.5)
    return float(r) / scale


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """解析 C-7.INT。data 为路径或 dict。"""
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("程序", PROGRAM_ID)
        return p

    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}
    idx = 0
    p["VG"] = float(nums[idx]); idx += 1
    p["VH"] = float(nums[idx]); idx += 1
    p["M"] = int(nums[idx]); idx += 1
    p["SA"] = int(nums[idx]); idx += 1
    p["BB"] = int(nums[idx]); idx += 1
    p["AB"] = int(nums[idx]); idx += 1
    p["AC"] = int(nums[idx]); idx += 1
    p["AM"] = float(nums[idx]); idx += 1
    p["EP"] = float(nums[idx]); idx += 1
    p["EC"] = float(nums[idx]); idx += 1
    p["B1"] = float(nums[idx]); idx += 1
    p["B2"] = float(nums[idx]); idx += 1
    p["B3"] = float(nums[idx]); idx += 1
    p["A0"] = float(nums[idx]); idx += 1
    p["A1"] = float(nums[idx]); idx += 1
    p["A2"] = float(nums[idx]); idx += 1
    p["A3"] = float(nums[idx]); idx += 1

    n = p["AC"] - p["AB"] + 1
    p["MI"] = [0] * n
    p["QJ"] = [0.0] * n
    for k in range(n):
        p["MI"][k] = int(nums[idx]); idx += 1
        p["QJ"][k] = float(nums[idx]); idx += 1

    m = p["M"]
    p["VV"] = [0.0] * (m + 1)
    p["HV"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["VV"][i] = float(nums[idx]); idx += 1
        p["HV"][i] = float(nums[idx]); idx += 1
    p["SS"] = [0.0] * (m + 1)
    p["HS"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["SS"][i] = float(nums[idx]); idx += 1
        p["HS"][i] = float(nums[idx]); idx += 1
    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _make_interp(SA, xs, ys):
    if int(SA) == 3:
        return lambda x: interp3_near(xs, ys, x)
    return lambda x: interp1d(xs, ys, x)


def _period(p, VB, Qj, QO, cap, zd):
    """
    给定时段初库容 VB、原始入库 Qj、发电流量 QO，返回该时段全部指标。
    弃水处理：V=VB+Qj−QO；V>VG 时弃水压至 VG；V<VH 时压至 VH（不放空）。
    处理后入库 QIp = QO + V − VB（净入库）；弃水 = Qj − QIp。
    """
    VG = p["VG"]
    VH = p["VH"]
    V = VB + Qj - QO
    spill = 0.0
    if V > VG:
        spill = V - VG
        V = VG
    if V < VH:
        V = VH
    QIp = QO + V - VB
    VC = 0.5 * (VB + V)
    HC = cap(VC)
    HS = zd(QO)
    DH = HC - HS
    N = 8.3 * QO * DH                   # 瓦
    return {
        "库容_raw": V, "弃水_raw": spill, "净入库_raw": QIp,
        "平均库容_raw": VC, "平均水位_raw": HC, "下游水位_raw": HS,
        "水头_raw": DH, "出力_raw": N, "时段末水位_raw": cap(V),
        "解q_raw": QO,
    }


def _solve_QO(p, VB, Qj, NT, cap, zd):
    """
    内层：解发电流量 q 使 NC(q)=NT。NC 对 q 单调增。
    约束上限 qmax = VB+Qj−VH（放空至死水位）；放空仍不足出力则取 qmax。
    """
    VH = p["VH"]

    def Nf(q):
        return _period(p, VB, Qj, q, cap, zd)["出力_raw"]

    qmax = VB + Qj - VH
    if qmax < 1e-6:
        qmax = 1e-6
    if Nf(qmax) <= NT:
        return qmax
    lo, hi = 1e-9, qmax
    if Nf(lo) >= NT:
        return lo
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if Nf(mid) < NT:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _run_chain(p, NT, cap, zd, QJ_list=None):
    """给定 NT，逐时段解 QO，返回 rows 与期末库容。"""
    n = len(p["QJ"])
    QJ_list = QJ_list if QJ_list is not None else p["QJ"]
    VB = float(p["VG"])
    rows = []
    for k in range(n):
        qo = _solve_QO(p, VB, QJ_list[k], NT, cap, zd)
        st = _period(p, VB, QJ_list[k], qo, cap, zd)
        rows.append((qo, st))
        VB = st["库容_raw"]
    return rows


def _equal_flow_QC(p):
    """等流量调节流量 QC = (ΣQj + DV)/F，DV = VG − VH（等流量口径）。"""
    QJ = p["QJ"]
    F = len(QJ)
    return (sum(QJ) + (p["VG"] - p["VH"])) / F


def compute(params, cfg=None):
    """
    等出力调节计算主流程。
    cfg: {"mode": "trial"|"fit", "chain": [{"QO":..,"Vend":..}, ...] (fit 用)}
    """
    cfg = cfg or {}
    mode = cfg.get("mode", "trial")
    p = params
    VG = float(p["VG"])
    VH = float(p["VH"])
    m = int(p["M"])
    QJ = p["QJ"]
    MI = p["MI"]
    n = len(QJ)

    xs_cap = p["VV"][1:m + 1]
    ys_cap = p["HV"][1:m + 1]
    xs_dn = p["SS"][1:m + 1]
    ys_dn = p["HS"][1:m + 1]
    cap = _make_interp(p["SA"], xs_cap, ys_cap)
    zd = _make_interp(p["SA"], xs_dn, ys_dn)

    QC = _equal_flow_QC(p)
    # 等流量初始出力 NT0 = A0 * (ΣN)/F（瓦）
    VB = VG
    NS = 0.0
    for k in range(n):
        st = _period(p, VB, QJ[k], QC, cap, zd)
        NS += st["出力_raw"]
        VB = st["库容_raw"]
    NT0_W = p["A0"] * NS / n

    rows_raw = []
    if mode == "fit":
        chain = cfg.get("chain")
        if not chain and cfg.get("chain_file"):
            import json
            with open(cfg["chain_file"], "r", encoding="utf-8") as f:
                chain = json.load(f)
        if not chain or len(chain) != n:
            raise ValueError("fit 模式需 cfg['chain'] 与时段数等长")
        VBl = VG
        for k in range(n):
            it = chain[k]
            if isinstance(it, dict):
                qo = float(it["QO"])
                Vend = float(it["Vend"])
            else:
                qo = float(it[0]); Vend = float(it[1])
            st = _period(p, VBl, qo + Vend - VBl, qo, cap, zd)
            st["库容_raw"] = Vend
            st["时段末水位_raw"] = cap(Vend)
            st["平均库容_raw"] = 0.5 * (VBl + Vend)
            st["平均水位_raw"] = cap(0.5 * (VBl + Vend))
            st["下游水位_raw"] = zd(qo)
            st["水头_raw"] = cap(0.5 * (VBl + Vend)) - zd(qo)
            st["出力_raw"] = 8.3 * qo * st["水头_raw"]
            st["净入库_raw"] = qo + Vend - VBl
            st["解q_raw"] = qo
            rows_raw.append(st)
            VBl = Vend
        NT_used = None
        VV_req = float(p.get("VAreq", VH))
    else:
        VV_req = float(p.get("VAreq", VH))
        lo, hi = 1.0e3, 2.0e6           # 瓦
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            rr = _run_chain(p, mid, cap, zd)
            if rr[-1][1]["库容_raw"] > VV_req:
                lo = mid                 # 期末库容偏高 → 需要更大出力（多放水）
            else:
                hi = mid
        NT_used = 0.5 * (lo + hi)
        rr = _run_chain(p, NT_used, cap, zd)
        rows_raw = [st for _, st in rr]

    rows = []
    for k in range(n):
        st = rows_raw[k]
        rows.append({
            "时段序号": MI[k],
            "原始入库_raw": QJ[k],
            "入库流量": _rhu(QJ[k], 2),
            "处理后入库_raw": st["净入库_raw"],
            "处理后入库": _rhu(st["净入库_raw"], 0),
            "调节流量_raw": st["解q_raw"],
            "调节流量": _rhu(st["解q_raw"], 0),
            "弃水流量之和_raw": QJ[k] - st["净入库_raw"],
            "弃水流量之和": _rhu(max(QJ[k] - st["净入库_raw"], 0.0), 0),
            "时段库容_raw": st["库容_raw"],
            "时段库容": _rhu(st["库容_raw"], 0),
            "时段水位_raw": st["时段末水位_raw"],
            "时段水位": _rhu(st["时段末水位_raw"], 2),
            "水头_raw": st["水头_raw"],
            "水头": _rhu(st["水头_raw"], 2),
            "出力_raw": st["出力_raw"],
            "出力": _rhu(st["出力_raw"], 0),
        })

    DQ = sum(r["弃水流量之和_raw"] for r in rows)
    SZ = VV_req - rows[-1]["时段库容_raw"]
    QG_AB = max(rows[0]["弃水流量之和_raw"], 0.0)
    NT_disp = (NT_used / 1.0e4) if NT_used is not None else (NT0_W / 1.0e4)

    return {
        "程序": PROGRAM_ID, "模式": mode,
        "VG": VG, "VH": VH, "M": m, "SA": int(p["SA"]), "BB": int(p["BB"]),
        "AB": p["AB"], "AC": p["AC"],
        "AM": float(p["AM"]), "EP": float(p["EP"]), "EC": float(p["EC"]),
        "B1": float(p["B1"]), "B2": float(p["B2"]), "B3": float(p["B3"]),
        "A0": float(p["A0"]), "A1": float(p["A1"]),
        "A2": float(p["A2"]), "A3": float(p["A3"]),
        "QC_raw": QC, "NT0_raw": NT0_W, "NT_raw": NT_used,
        "NT": _rhu(NT_disp, 2),
        "QG_AB_raw": QG_AB, "QG_AB": _rhu(QG_AB, 2),
        "DQ_raw": DQ, "DQ": _rhu(DQ, 0),
        "SZ_raw": SZ, "SZ": _rhu(SZ, 2),
        "时段序号": MI, "原始入库": QJ, "rows": rows,
    }


# ------------------------------------------------------------
# 输出
# ------------------------------------------------------------

def _place_line(width, fields):
    line = [" "] * width
    for end, text in fields:
        start = end - len(text)
        if start < 0:
            start = 0
        for j, ch in enumerate(text):
            if start + j < width:
                line[start + j] = ch
    return "".join(line).rstrip()


def render(params, result, table=None):
    """生成文本计算书（原著 C-7.OUT 风格）。"""
    r = result
    line = "*" * 72
    out = []
    out.append(line)
    out.append("*****                水电站等出力调节计算书（试算法）              *****")
    out.append(line)
    out.append("")
    out.append(" 输入数据:")
    out.append("     正常蓄水位库容 VG= %9.0f" % r["VG"])
    out.append("     最低限制水位的库容 VH= %8.1f" % r["VH"])
    out.append("     结点数 M= %2d " % r["M"])
    out.append("     控制变量 SA= %2d     BB= %2d     BK= 0 " % (r["SA"], r["BB"]))
    out.append("     控制弃水处理的变量 B1= %2.0f " % r["B1"])
    out.append("     控制修改NT的变量 B2= %3.0f " % r["B2"])
    out.append("     控制缩短NT步长的变量 B3= %2.0f " % r["B3"])
    out.append("     计算期首时段序号 AB= %2d " % r["AB"])
    out.append("     计算期末时段序号 AC= %2d " % r["AC"])
    out.append("     计算NT时增大的倍数 A0= %4.0f      控制出力步长的系数 A1= %2.0f " % (r["A0"], r["A1"]))
    out.append("     计算NT时的步长 A2= %7.4f 控制AM步长的系数 A3=  %6.4f" % (r["A2"], r["A3"]))
    out.append("")
    out.append("  时段序号     入库流量")
    for row in r["rows"]:
        out.append(_place_line(24, [(6, "%4d" % row["时段序号"]),
                                    (23, "%10.2f" % row["原始入库_raw"])]))
    out.append("")
    out.append(" 结点数    库容     上游水位     流量     下游水位")
    VV = params["VV"]; HV = params["HV"]; SS = params["SS"]; HS = params["HS"]
    for i in range(1, r["M"] + 1):
        out.append(_place_line(46, [
            (3, "%2d" % i), (13, "%9.2f" % VV[i]), (24, "%9.2f" % HV[i]),
            (34, "%9.0f" % SS[i]), (46, "%9.2f" % HS[i]),
        ]))
    out.append("")
    out.append(" 计算结果：  ")
    out.append(" 库水位高于正常水位时的弃水流量 QG(AB)= %6.2f" % r["QG_AB"])
    out.append(" 等出力计算的初始出力 NT= %7.2f" % r["NT"])
    out.append(" 等出力计算时，要求弃水时段入库流量aQI(I)增大的倍数 AM= %7.2f" % r["AM"])
    out.append(" 期末库容量 SZ= %6.2f" % r["SZ"])
    out.append("")
    out.append("               处理后          弃水流量")
    out.append(" 时段 入库流量 入库流量 调节流量  之和 时段末库容 时段末水位 水头  出力")
    out.append(" " + "-" * 71)
    for row in r["rows"]:
        out.append(_place_line(80, [
            (3, "%2d" % row["时段序号"]),
            (10, "%5.0f" % row["入库流量"]),
            (19, "%6.0f" % row["处理后入库"]),
            (29, "%6.0f" % row["调节流量"]),
            (38, "%6.0f" % row["弃水流量之和"]),
            (50, "%6.0f" % row["时段库容"]),
            (61, "%8.2f" % row["时段水位"]),
            (70, "%6.2f" % row["水头"]),
            (80, "%7.0f" % row["出力"]),
        ]))
    out.append("")
    out.append("    弃水流量 DQ= %5.0f    期末库容差 SZ= %6.2f" % (r["DQ"], r["SZ"]))
    out.append("    等出力计算时，要求弃水时段入库流量aQI(I)的增大倍数 AM= %5.2f" % r["AM"])
    out.append("")
    return "\n".join(out)


def run(data, out_txt=None, out_json=None, fmt="text", cfg=None):
    params = parse(data)
    result = compute(params, cfg)
    text = render(params, result, None)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(txt)
