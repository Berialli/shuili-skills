# -*- coding: utf-8 -*-
"""
C-8 水电站定出力调节计算程序（图解法）—— 内核
================================================
复刻《水利程序集》C-8 程序（作者：唐文华，水电部天津勘测设计院）。

功能：已知计算期内各时段的预定出力 NT、入库流量 QI，以及库容曲线
（库容 VV ~ 水位 HV）、下游流量水位曲线（流量 SS ~ 水位 HS）、计算期
消落深度（VG 起、VH 止）与消落范围内的发电水头网格 HH(1..MF)，用
"图解法"（不迭代）逐时段定出力调节，得到各时段发电流量 QO、时段末
库容 V（及水位 H）等水能成果。

────────────────────────────────────────────────────────────
算法（原著"计算方法概述"，公式(1)~(4)）
────────────────────────────────────────────────────────────
把消落深度 VG−VH 划分为 MF 个发电水头 DH（HH 数组，范围略大于
VG−VH），对每个 DH 用 (1) 式算发电流量：

    QO = NT · CC / DH                                    …(1)

  NT 单位万千瓦时，CC = 出力系数倒数的 10000 倍（算例 1204.819277
  = 10000/8.3）。本算例 HS 曲线恒为 18.45 m，故式(1)中的 DQ（下游
  流量，取 QO）对 HS 无影响；一般情形由 QO 在 SS~HS 曲线上查下游
  水位 H 下。

由 DH 与下游水位相加得平均水位，再在库容曲线上反查（水位→库容）
得平均库容 VP，最后按 (2) 式算 FA：

    FA = VP + QO/2                                       …(2)

对 MF 个 DH 算得 MF 组 (FA, QO)，即构成该出力 NT 的"水能计算工作
曲线 FA ~ NT ~ QO"。

逐时段计算：已知时段初库容 VB 与入库流量 QI，由 (3) 式算 FA：

    FA = VB + QI/2                                       …(3)

由该 FA 在工作曲线 FA ~ NT ~ QO 上插值得发电流量 QO，再由 (4) 式算
时段末库容 V：

    V = VB + QI − QO                                     …(4)

同时可由 V 在库容曲线上查得时段末水位 H。若 V ≥ VH 则转入下一时段；
否则 NT 偏大 / 有效库容偏小，须修改 NT 重算（原著处理）。

────────────────────────────────────────────────────────────
黑盒反演裁定的实现口径（对拍权威 C-8.OUT 7 行 × 2 列逐位全中）
────────────────────────────────────────────────────────────
1. 插值（SU=3，三点插值）统一用"最近端点侧三点 Lagrange 二次插值"
    interp3_near：x 落段 [xs[s], xs[s+1]] 时，距左端点更近（或等）取
   三点 (s−1,s,s+1)，距右端点更近取 (s,s+1,s+2)；边界退化、越界取
   端点值。
   · 用途 A：工作曲线构建时"平均水位 → 平均库容"的反查（自变量 HV
     升序、因变量 VV）。
   · 用途 B：逐时段"FA → QO"在工作曲线上的查询。
2. 平均水位（= DH + 下游水位）低于库容曲线最低水位（49 m）时，取
   曲线首两点线性外延（本算例 HH=29 时平均水位 47.45 < 49 触发，
   若简单取端点值会使该结点 FA 显著抬高、曲线在低 FA 段畸形，导致
   末几行 QO 偏差 40 以上）。
3. 工作曲线按 FA 升序排序后再查询（sort=True）：HH 网格低端因顶点
   外推使 FA 不再单调，排序后曲线在查询支上光滑；不排序则在低 FA
   段命中错误的非单调分支。
4. 时段初库容 VB 以"连续真实值"逐时段滚动传递（VB ← V，非显示
   舍入值）；显示列用 half-away-from-zero 舍入（原著 VB 语义）。
   裁决依据：以权威 OUT 的整型 V 作 VB 时 QO 亦 7/7 命中，但以连续
   值滚动时 QO 与 V 两列同时 7/7 命中（oscore=14/14），连续滚动是
   唯一能同时命中两列的口径。
5. 输出 4 个：QO(j) 发电流量、V(j) 时段末库容、N(j) 同 NT(j)、
   DQ 弃水流量之和（本算例全时段无弃水，DQ=0）。

数据文件顺序（C-8.INT，共 10 控制数 + 2×7 时段数 + 10+10 曲线点
+ MF 水头）：
  VG, VH, M, MF, SU, AB, AC, B1, B2, CC
    VG  计算期初库容（=正常蓄水位相应库容，秒立米月）
    VH  计算期末库容（=最低限制水位相应库容）
    M   库容曲线 / 下游流量水位曲线的结点数
    MF  水能计算工作曲线的结点数
    SU  =2 两点（线性）插值；=3 三点插值
    AB  计算期首时段序号；AC 计算期末时段序号
    B1  库容大于 VG 时控制弃水的系数（0.5~5 秒立米月）
    B2  库容小于 VH 时控制溢出的系数（1~10 秒立米月）
    CC  出力系数倒数的 10000 倍（= 10000/8.3 = 1204.819277）
  之后 (AB..AC) 每时段：M(I) 时段序号、QI(I) 入库流量（交错对）
  之后 (AB..AC) 每时段：NT(I) 各时段预定出力（万千瓦）
  之后 (1..M) 每点：VV(j) 库容、HV(j) 水位（交错对）
  之后 (1..M) 每点：SS(K) 流量、HS(K) 下游水位（交错对）
  之后 (1..MF) 每点：HH(j) 按消落深度 VG−VH 计算的发电水头

验证基准：C-8.OUT（算例，7 时段 11,12,1,2,3,4,5 月）
  VG=1975, VH=1215, M=10, MF=10, SU=3, AB=1, AC=7, B1=5, B2=1,
  CC=1204.819277, HH=29..38；
  预定出力 8.43/8.24/7.98/7.63/7.49/7.25/7.16 万千瓦；
  自 VB=VG=1975 起调，全时段无弃水，时段末库容 1869→1725→1568→
  1435→1359→1288→1245 逐时段回落，均高于 VH=1215；
  主表 7 行 × QO/V 两列（14 项）逐位命中（见 c8_verify.py）。
"""
from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "C-8"
TITLE = "水电站定出力调节计算书（图解法）"


# ------------------------------------------------------------
# 插值（原著 SU=3 三点插值：最近端点侧三点 Lagrange）
# ------------------------------------------------------------

def _seg_of(xs, x):
    """返回 x 所在段 s：xs[s] <= x <= xs[s+1]。调用前保证 x 在界内。"""
    for s in range(len(xs) - 1):
        if xs[s] <= x <= xs[s + 1]:
            return s
    if x >= xs[-1]:
        return len(xs) - 2
    return 0


def interp3_near(xs, ys, x):
    """
    最近端点侧三点 Lagrange 二次插值（SU=3 语义，黑盒反演裁定）。
    x ∈ [xs[s], xs[s+1]]：
      距左端点更近(或等) → 三点 (s-1,s,s+1)；距右端点更近 → (s,s+1,s+2)。
    越界取端点 ys；边界处退化为可用三点。
    """
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    s = _seg_of(xs, x)
    if (x - xs[s]) <= (xs[s + 1] - x):      # 距左端点更近(或等)：取左侧三点
        a, b, c = s - 1, s, s + 1
        if a < 0:
            a, b, c = 0, 1, 2
    else:                                   # 距右端点更近：取右侧三点
        a, b, c = s, s + 1, s + 2
        if c > len(xs) - 1:
            a, b, c = len(xs) - 3, len(xs) - 2, len(xs) - 1
    x0, x1, x2 = xs[a], xs[b], xs[c]
    y0, y1, y2 = ys[a], ys[b], ys[c]
    return (y0 * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))
            + y1 * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
            + y2 * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)))


def _lagrange3(x0, x1, x2, y0, y1, y2, x):
    """三点 Lagrange 二次插值（给定三点与 x）。"""
    return (y0 * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))
            + y1 * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
            + y2 * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)))


def interp1d(xs, ys, x):
    """分段线性插值（SU=2 语义）：越界取端点，段内直线内插。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    s = _seg_of(xs, x)
    t = (x - xs[s]) / (xs[s + 1] - xs[s])
    return ys[s] + t * (ys[s + 1] - ys[s])


def _make_interp(SU):
    """按 SU 选插值函数。SU=3 三点（最近端点侧）；SU=2 线性。"""
    if int(SU) == 3:
        return interp3_near
    return interp1d


# ------------------------------------------------------------
# 显示舍入（原著 VB 显示为 half-away-from-zero，非银行家舍入）
# ------------------------------------------------------------

def _rhu(x, nd=0):
    """round half up（远离零），nd 为小数位数。"""
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
    """
    解析输入。
    data: dict 或 INT 文件路径。
    返回标准参数字典。
    """
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("程序", PROGRAM_ID)
        return p

    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}
    idx = 0

    p["VG"] = float(nums[idx]); idx += 1     # 计算期初库容
    p["VH"] = float(nums[idx]); idx += 1     # 计算期末库容
    p["M"] = int(nums[idx]); idx += 1        # 曲线结点数
    p["MF"] = int(nums[idx]); idx += 1       # 工作曲线结点数
    p["SU"] = int(nums[idx]); idx += 1       # 插值控制 2/3
    p["AB"] = int(nums[idx]); idx += 1       # 计算期首时段序号
    p["AC"] = int(nums[idx]); idx += 1       # 计算期末时段序号
    p["B1"] = float(nums[idx]); idx += 1     # 弃水控制系数
    p["B2"] = float(nums[idx]); idx += 1     # 溢出控制系数
    p["CC"] = float(nums[idx]); idx += 1     # 出力系数倒数的 1e4 倍

    n = p["AC"] - p["AB"] + 1
    p["MI"] = [0] * n                        # 时段序号（绝对）
    p["QI"] = [0.0] * n                      # 入库流量
    for k in range(n):
        p["MI"][k] = int(nums[idx]); idx += 1
        p["QI"][k] = nums[idx]; idx += 1

    p["NT"] = [0.0] * n                      # 各时段预定出力（万kW）
    for k in range(n):
        p["NT"][k] = nums[idx]; idx += 1

    m = p["M"]
    p["VV"] = [0.0] * (m + 1)                # 1 基：库容曲线
    p["HV"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["VV"][i] = nums[idx]; idx += 1
        p["HV"][i] = nums[idx]; idx += 1

    p["SS"] = [0.0] * (m + 1)                # 1 基：下游流量水位曲线
    p["HS"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["SS"][i] = nums[idx]; idx += 1
        p["HS"][i] = nums[idx]; idx += 1

    mf = p["MF"]
    p["HH"] = [0.0] * mf                     # 发电水头网格
    for j in range(mf):
        p["HH"][j] = nums[idx]; idx += 1

    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def build_work_curve(p, NT):
    """
    构建出力 NT（万千瓦）的水能计算工作曲线 FA ~ NT ~ QO。
    返回 (FA 升序列表, QO 列表)。

    对每个发电水头 DH = HH[j]：
      QO_j = NT·CC/DH                      …(1)
      HS_j = 下游水位(QO_j)                （本算例恒 18.45）
      HU_j = DH + HS_j                     （平均水位）
      VP_j = 库容反查(HU_j)（水位→库容）    （三点插值；低于曲线下界时线性外延）
      FA_j = VP_j + QO_j/2                 …(2)
    """
    CC = p["CC"]
    HH = p["HH"]
    m = p["M"]
    xs_dn = p["SS"][1:m + 1]
    ys_dn = p["HS"][1:m + 1]
    xs_lv = p["HV"][1:m + 1]         # 水位（升序）
    ys_st = p["VV"][1:m + 1]         # 库容
    f_inv = _make_interp(p["SU"])

    pts = []
    for dh in HH:
        qo = NT * CC / dh                    # (1)
        hu = dh + interp1d(xs_dn, ys_dn, qo)  # 平均水位
        if hu >= xs_lv[0]:
            vp = f_inv(xs_lv, ys_st, hu)     # 水位→库容 反查
        else:
            # 低于库容曲线最低水位：以首两点线性外延（黑盒反演裁定）
            k = (ys_st[1] - ys_st[0]) / (xs_lv[1] - xs_lv[0])
            vp = ys_st[0] + k * (hu - xs_lv[0])
        pts.append((vp + qo / 2.0, qo))      # (2)
    # 按 FA 升序排序（HH 低端外推会使 FA 非单调，排序后查询支光滑）
    pts.sort(key=lambda t: t[0])
    return [t[0] for t in pts], [t[1] for t in pts]


def compute(params, cfg=None):
    """
    图解法定出力调节计算主流程。
    cfg: dict，可选：
      mode: "trial"（默认，正向图解法，逐时段用该时段 NT 建工作曲线）
            | "fit"（核算模式，用给定 QO 序列直接滚动复现权威输出）。
      qo:   fit 模式下逐时段发电流量序列（连续值）。
    返回结果字典（含逐时段明细表 rows 与汇总）。
    """
    cfg = cfg or {}
    mode = cfg.get("mode", "trial")
    p = params
    VG = float(p["VG"])
    VH = float(p["VH"])
    n = len(p["QI"])
    QI = p["QI"]
    MI = p["MI"]
    NT_w = p["NT"]                      # 预定出力，万kW
    f_q = _make_interp(p["SU"])         # 查询插值（FA→QO）

    fit_qo = cfg.get("qo") if mode == "fit" else None
    if mode == "fit":
        if not fit_qo or len(fit_qo) != n:
            raise ValueError("fit 模式需提供与时段数等长的 qo 序列（cfg['qo']）")

    VB = VG                             # 定出力调节自计算期初库容起调
    rows = []
    warnings = []
    for k in range(n):
        NT = NT_w[k]
        FA_c, QO_c = build_work_curve(p, NT)
        fa = VB + QI[k] / 2.0           # (3)
        if mode == "fit":
            qo = float(fit_qo[k])
        else:
            qo = f_q(FA_c, QO_c, fa)     # 工作曲线上查 QO
        v = VB + QI[k] - qo             # (4)
        qg = 0.0
        if v > VG:
            # 库容超过计算期初库容：弃水（原著以 B1 控制弃水步长）
            qg = v - VG
            v = VG
        if v < VH:
            warnings.append(
                f"时段 {MI[k]}: 时段末库容 {v:.2f} 低于计算期末库容 "
                f"VH={VH:.1f}（出力偏大或有效库容偏小，原著须修改 NT 重算；"
                f"本内核照常打印已算成果）")

        h_end = interp1d(p["VV"][1:p["M"] + 1], p["HV"][1:p["M"] + 1], v)

        rows.append({
            "时段序号": MI[k],
            "入库流量_raw": QI[k],
            "预定出力_raw": NT,
            "模式": mode,
            "曲线FA": FA_c,
            "曲线QO": QO_c,
            "查询FA_raw": fa,
            "调节流量_raw": qo,
            "时段末库容_raw": v,
            "时段末水位_raw": h_end,
            "弃水_raw": qg,
            # 显示列（原著 half-away 舍入）
            "入库流量": _rhu(QI[k], 2),
            "调节流量": _rhu(qo, 0),
            "时段库容": _rhu(v, 0),
            "时段水位": _rhu(h_end, 2),
            "预定出力显示": _rhu(NT * 1.0e4, 0),
        })
        VB = v                          # 连续值滚动传递

    dq = sum(r["弃水_raw"] for r in rows)
    result = {
        "程序": PROGRAM_ID,
        "VG": VG, "VH": VH, "M": int(p["M"]), "MF": int(p["MF"]),
        "SU": int(p["SU"]), "AB": int(p["AB"]), "AC": int(p["AC"]),
        "B1": float(p["B1"]), "B2": float(p["B2"]), "CC": float(p["CC"]),
        "HH": list(p["HH"]),
        "时段序号": MI,
        "入库流量": QI,
        "预定出力": NT_w,
        "调节流量": [r["调节流量"] for r in rows],
        "时段末库容": [r["时段库容"] for r in rows],
        "时段末水位": [r["时段水位"] for r in rows],
        "弃水流量之和": dq,
        "rows": rows,
        "warnings": warnings,
    }
    return result


# ------------------------------------------------------------
# 输出
# ------------------------------------------------------------

def _place_line(width, fields):
    """按"字段右端对齐"落字构造一行文本。"""
    line = [" "] * width
    for end, text in fields:
        start = end - len(text)
        if start < 0:
            start = 0
            text = text[-end:] if len(text) > end else text
        for j, ch in enumerate(text):
            if start + j < width:
                line[start + j] = ch
    return "".join(line).rstrip()


def render(params, result, table=None):
    """
    生成文本计算书（原著 C-8.OUT 风格，主表列右端位置逐字对齐）。
    主表列右端位置量取自权威 C-8.OUT：
      时段月序 端5  入库流量 端19  发电流量 端32  时段末库容 端46
      各时段预定出力 端63
    """
    r = result
    line = " " + "*" * 72
    out = []
    out.append(line)
    out.append(" *****               水电站定出力调节计算书（图解法）               *****")
    out.append(line)
    out.append("")
    out.append(" 输入数据:")
    out.append(f"     计算期初库容 VG= {r['VG']:>8.1f}")
    out.append(f"     计算期末库容 VH= {r['VH']:>8.1f}")
    out.append(f"     库容曲线的结点数 M= {r['M']:>1d} ")
    out.append(f"     水能计算工作曲线的结点数 MF= {r['MF']:>1d} ")
    out.append(f"     插值计算控制变量 SU= {r['SU']:>1d} ")
    out.append(f"     计算期首时段序号 AB= {r['AB']:>1d} ")
    out.append(f"     计算期末时段序号 AC= {r['AC']:>1d} ")
    out.append(f"     库容大于VG时控制弃水的系数 B1= {r['B1']:>7.3f}")
    out.append(f"     库容小于VH时控制溢出的系数 B2= {r['B2']:>7.3f}")
    out.append("")
    out.append(" 水能成果表：  ")
    out.append("")
    out.append(" 时段月序  入库流量 QI  发电流量 QO  时段末库容 V   各时段预定出力 N")
    out.append(" " + "-" * 67)
    for row in r["rows"]:
        out.append(_place_line(63, [
            (5, f"{row['时段序号']:>3d}"),
            (19, f"{row['入库流量']:>5.0f}"),
            (32, f"{row['调节流量']:>5.0f}"),
            (46, f"{row['时段库容']:>6.0f}"),
            (63, f"{row['预定出力显示']:>6.0f}"),
        ]))
    out.append("")
    out.append(f" 弃水流量之和 DQ= {r['弃水流量之和']:>7.0f}")
    out.append(f" 计算期初库容 VG= {r['VG']:>7.0f}")
    out.append(f" 计算期末库容 VH= {r['VH']:>7.0f}")
    if r.get("warnings"):
        out.append("")
        for w in r["warnings"]:
            out.append("  警告: " + w)
    out.append("")
    return "\n".join(out)


def run(data, out_txt=None, out_json=None, fmt="text", cfg=None):
    """
    统一入口。
    data: INT 文件路径 | dict
    """
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
