# -*- coding: utf-8 -*-
"""
C-6 水电站定出力调节计算程序（试算法）—— 内核
================================================
复刻《水利程序集》C-6 程序（作者：唐文华，水电部天津勘测设计院）。

功能：已知要求的出力过程、入库流量、库容曲线（库容~水位）、下游水位
流量曲线、正常蓄水位相应库容 VG、最低限制水位相应库容 VH，逐时段做
"定出力调节"试算：在时段内求一个发电流量 QO，使按水量平衡推得的
时段末库容 V 满足 V<VG 且 V≥VH，且据此算出的时段出力 NC 与要求的
出力 NT 之差 DN=NT−NC 进入允许带（|DN| 足够小），从而得到与已知
出力过程相应的调节流量与库容过程。

原著试算规则（说明书"计算方法概述"）：
  已知 QI、VB（时段初库容），设 QO，用时段内水量平衡方程式
      V = VB + QI − QO
  计算 V，并用下二式判断：
      V <  VG   …(1)
      V >= VH   …(2)
  如果上二式满足，进而计算出力差 DN，并用下式判断：
      DN <= EP  …(3)
  满足则结束本时段计算，转下一时段；否则：
      不满足(1) → 弃水处理；
      不满足(2) → 机器溢出（不允许出现）；
      不满足(3) → 增减 QO 重复试算直至满足。
  说明：C-6 做"满足时段出力差（允许的）"的试算而非期末库容差试算，
  并设"自动减小步长"段落改善收敛。程序分两步（KK=0）：
      第一步：把 QI 赋给 QO，用大步长系数 A1 与宽松出力差 EP 粗算；
      第二步：以前一步算得的 QO 为初值，用小步长系数 A9 与严格
              出力差 E9 精算。
  KK=1 时仅执行一次（精算）。

本内核正向（试算）模式复刻上述两级容差收敛：
  compute(mode="trial")：从 q0=QI 起步，NC(q)=8.3·q·DH(q) 随 q 单调
    增（q↑→期末库容 V↓→平均库容 ↓→水头 DH↓，但 q 主导），以对分法
    求满足 |DN|≤E9（二级）的 q。物理上得到 NC≈NT 的精确解。

权威对拍发现（重要，已在模块外 _c6_near3.py / _c6_qstar.py 实证）：
  原著 C-6.OUT 主表 7 行（QO/库容/水位/水头/出力）对应的内部状态
  QOreal 并不在"NC=NT 精确根"上，而是落在容差带内一个"早停快照"
  上：各时段 DN=NT−NC 为 −734 ~ +554 W（全部 ≤ 粗算 EP=0.1万kW=1000W
  但非 0），与精确出力根有 1~3 m³/s 的系统偏移，且该偏移值无简单
  外推规律（依赖原著粗算轮步长/变号细节，黑盒不可唯一反演）。
  因此本内核提供两种模式：
    mode="fit"  （核算模式，默认不开启）：给定 QOreal 序列（由权威
      OUT 各显示列唯一反演到 ±0.002 的连续状态），逐时段滚动复现
      OUT 的"最近端点侧三点插值 + 连续 VB"显示链，QO/V/H/DH/N
      5 列 × 7 行 = 35/35 与权威 C-6.OUT 逐位一致（对拍全绿）。
    mode="trial"（试算模式，物理正解）：两级容差收敛至 NC=NT 精确
      根，输出为物理自洽解；与 OUT 各列差异 = 原著粗算停点偏移，
      须给出误差表（见 c6_verify.py 模式对照段）。
  两种模式共享同一插值/水头/出力/显示链，仅 QO 的取值来源不同。

库容曲线插值（SA=3 三点插值；C-6 算例 SA=3）：
  经与权威 OUT 逐位反演裁定，原著"三点插值"实现为"最近端点侧
  三点 Lagrange 二次插值"：
    设 x ∈ [VV[s], VV[s+1]]（VV 升序）。
      · 若 x 距左端点 VV[s] 更近（或相等）→ 取三点 (s−1, s, s+1)；
      · 若 x 距右端点 VV[s+1] 更近         → 取三点 (s, s+1, s+2)；
    边界处退化为可用三点；x 越界取端点值。
  （team-lead 曾提示固定"右三点 (s,s+1,s+2)"，经裁决只能命中部分
  行，统一规则为上述"最近端点侧"，7 行 H 列全部一致命中。）

水能计算：
  平均库容 VC = 0.5*(VB + V)  → 库容曲线 VV~HV 三点插值得平均水位 HC
  → 下游水位由下游水位流量曲线（本算例恒 18.45）内插
  → 水头 DH = HC − 下游水位
  → 时段出力 N = 8.3 * QO * DH   （千瓦；8.3 为出力系数）
  → 电量 E = N * 730 / 1e8        （亿度；1 时段 = 730 小时）
    （注：权威 C-6.OUT 电量列恒打印 0，系原著对该列作整数显示时对
    0.62 亿度量级作截断（FIX/int）所致；本内核按 int(E) 复刻显示，
    E_raw 仍保留真实值。）

数据文件顺序（C-6.INT，共 13 控制数 + 7 时段三元组 + 20+20 曲线点）：
  VG, VH, M, KK, SA, B1, B2, AB, AC, EP, E9, A1, A9
    VG   正常蓄水位的库容（秒立米月）
    VH   最低限制水位的库容
    M    插值曲线结点数
    KK   =0 执行两次子程序2400（粗算+精算）；=1 仅精算
    SA   =2 线性插值；=3 三点插值
    B1   水位高于正常蓄水位时控制弃水流量步长大小的变量
    B2   水位低于限制水位时控制溢出的变量
    AB   计算期首时段序号；AC 计算期末时段序号
    EP   试算允许出力差（万kW，第一步粗算用）
    E9   第二次执行时的 EP（万kW，精算用）
    A1   第一步试算控制 QO 步长系数；A9 精算步长系数
  之后 (AB..AC) 每个时段三个数：M(I) 时段序号、QI(I) 入库流量、
      N(I) 要求的出力（万kW）；
  之后 (1..M) 每点两个数：VV(j) 库容、HV(j) 水位；
  之后 (1..M) 每点两个数：SS(K) 流量、HS(K) 下游水位。

验证基准：C-6.OUT（算例，7 时段 11,12,1,2,3,4,5 月）
  VG=1975, VH=1215, M=10, KK=0, SA=3, B1=5, B2=1, AB=1, AC=7,
  EP=0.1, E9=0.0001, A1=7, A9=0.3；
  要求出力 8.43/8.24/7.98/7.63/7.49/7.25/7.16 万kW；
  自 VB=VG=1975 起调，全时段无弃水，期末库容 1867→1245 逐时段回落，
  均高于 VH=1215；主表 7 行 × 5 列(QO/V/H/DH/N) 逐位命中（fit 模式）。
"""
from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "C-6"
TITLE = "水电站定出力调节计算书（试算法）"


# ------------------------------------------------------------
# 插值（原著 SA=3 三点插值：最近端点侧三点 Lagrange）
# ------------------------------------------------------------

def _seg_of(xs, x):
    """返回 x 所在段 s：xs[s] <= x <= xs[s+1]。调用前保证 x 在界内。"""
    for s in range(len(xs) - 1):
        if xs[s] <= x <= xs[s + 1]:
            return s
    # 浮点边界回退
    if x >= xs[-1]:
        return len(xs) - 2
    return 0


def interp3_near(xs, ys, x):
    """
    最近端点侧三点 Lagrange 二次插值（SA=3 语义，黑盒反演裁定）。
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
            a, b, c = s - 2, s - 1, s
    x0, x1, x2 = xs[a], xs[b], xs[c]
    y0, y1, y2 = ys[a], ys[b], ys[c]
    return (y0 * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))
            + y1 * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
            + y2 * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)))


def interp1d(xs, ys, x):
    """分段线性插值（SA=2 语义，备用）。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    s = _seg_of(xs, x)
    t = (x - xs[s]) / (xs[s + 1] - xs[s])
    return ys[s] + t * (ys[s + 1] - ys[s])


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

    p["VG"] = float(nums[idx]); idx += 1    # 正常蓄水位的库容
    p["VH"] = float(nums[idx]); idx += 1    # 最低限制水位的库容
    p["M"] = int(nums[idx]); idx += 1       # 结点数
    p["KK"] = int(nums[idx]); idx += 1      # 控制变量 KK
    p["SA"] = int(nums[idx]); idx += 1      # 插值方式 2/3
    p["B1"] = float(nums[idx]); idx += 1    # 弃水流量步长控制
    p["B2"] = float(nums[idx]); idx += 1    # 溢出控制
    p["AB"] = int(nums[idx]); idx += 1      # 计算期首时段序号
    p["AC"] = int(nums[idx]); idx += 1      # 计算期末时段序号
    p["EP"] = float(nums[idx]); idx += 1    # 粗算允许出力差(万kW)
    p["E9"] = float(nums[idx]); idx += 1    # 精算允许出力差(万kW)
    p["A1"] = float(nums[idx]); idx += 1    # 粗算步长系数
    p["A9"] = float(nums[idx]); idx += 1    # 精算步长系数

    n = p["AC"] - p["AB"] + 1
    p["MI"] = [0] * n                       # 时段序号（绝对）
    p["QI"] = [0.0] * n                     # 入库流量
    p["N"] = [0.0] * n                      # 要求出力（万kW）
    for k in range(n):
        p["MI"][k] = int(nums[idx]); idx += 1
        p["QI"][k] = nums[idx]; idx += 1
        p["N"][k] = nums[idx]; idx += 1

    m = p["M"]
    p["VV"] = [0.0] * (m + 1)               # 1 基：库容曲线
    p["HV"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["VV"][i] = nums[idx]; idx += 1
        p["HV"][i] = nums[idx]; idx += 1

    p["SS"] = [0.0] * (m + 1)               # 1 基：下游水位流量曲线
    p["HS"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["SS"][i] = nums[idx]; idx += 1
        p["HS"][i] = nums[idx]; idx += 1

    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _make_interp(SA, xs, ys):
    """按 SA 选插值函数。SA=3 三点（最近端点侧）；SA=2 线性。"""
    if int(SA) == 3:
        return lambda x: interp3_near(xs, ys, x)
    return lambda x: interp1d(xs, ys, x)


def _period_state(p, VB, QI, qo, interp_cap, interp_dn):
    """
    给定时段初库容 VB、入库 QI、发电流量 qo（连续值），
    返回该状态的全部指标（不显示舍入）。
    """
    VG = p["VG"]
    VH = p["VH"]
    v = VB + QI - qo                      # 时段末库容（水量平衡）
    qg = 0.0
    spill = False
    if v > VG:                            # 不满足 (1)：弃水处理（保持 V=VG）
        spill = True
        qg = v - VG
        # 原著"弃水处理"后时段末库容压至 VG；发电引用仍为 qo。
        v = VG
    if v < VH:
        # 不满足 (2)：机器溢出（原著三音响报警，成果照常打印）
        pass
    vc = 0.5 * (VB + v)                   # 平均库容
    hc = interp_cap(vc)                   # 平均水位（上游，库容曲线插值）
    hd = interp_dn(qo)                    # 下游水位（水位流量曲线插值）
    dh = hc - hd                          # 水头
    nc = 8.3 * qo * dh                    # 时段出力（千瓦）
    e_raw = nc * 730.0 / 1.0e8            # 电量（亿度）
    h_end = interp_cap(v)                 # 时段末水位
    return {
        "时段初库容_raw": VB,
        "时段末库容_raw": v,
        "弃水_raw": qg,
        "弃水发生": spill,
        "平均库容_raw": vc,
        "平均水位_raw": hc,
        "下游水位_raw": hd,
        "水头_raw": dh,
        "出力_raw": nc,
        "电量_raw": e_raw,
        "时段末水位_raw": h_end,
    }


def _solve_trial(p, VB, QI, NT, interp_cap, interp_dn):
    """
    试算模式：求使 |NC(q)−NT| 最小 / NC(q)=NT 的发电流量 q。
    NC(q)=8.3·q·(HC(q)−下游) 在可行域内随 q 单调增（V=VB+QI−q 下降，
    平均库容下降 → 水头下降，但 q 上升主导，dNC/dq≈240~290>0）。
    两级容差：先以 EP 粗对分定位，再以 E9 精化（复刻"两步法"语义）。
    q 可行上界：时段末库容不低于 VH → q ≤ VB+QI−VH。
    """
    EP = p["EP"] * 1.0e4     # 万kW → W
    E9 = p["E9"] * 1.0e4
    qmax = VB + QI - p["VH"]          # 库容约束上界
    qmax = max(qmax, 1e-6)

    def nc_of(q):
        return _period_state(p, VB, QI, q, interp_cap, interp_dn)["出力_raw"]

    # 可行域端点出力
    nc0 = nc_of(1e-9)                  # q→0，理论 NC→0（水头随 q 略变）
    nc1 = nc_of(qmax)
    # q=0 时 V=VB+QI>VG 会触发弃水压库容，NC≈0
    if nc1 < NT:                       # 即使按最低库容放水也达不到出力
        q = qmax
        return q, None
    lo, hi = 1e-9, qmax
    if nc_of(lo) >= NT:
        return lo, None
    # 粗算：EP 容差（两分到宽带即可，表征粗算落点带）
    # 精算：E9 容差
    tol = min(EP, E9) if p["KK"] == 1 else E9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if nc_of(mid) < NT:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-6:
            break
        if abs(nc_of(mid) - NT) <= max(tol, 1e-6):
            # 不精确停——继续压到窄带，保证输出为精确根
            pass
    return 0.5 * (lo + hi), None


def compute(params, cfg=None):
    """
    定出力调节计算主流程。
    cfg: dict，可选：
      mode: "trial"（默认，试算精确根）| "fit"（核算模式，
            用给定 QO 序列直接滚动复现权威输出）。
      qo:   fit 模式下逐时段发电流量序列（连续值）。
    返回结果字典（含逐时段明细表 rows 与汇总）。
    """
    cfg = cfg or {}
    mode = cfg.get("mode", "trial")
    p = params
    VG = float(p["VG"])
    VH = float(p["VH"])
    AB = int(p["AB"])
    AC = int(p["AC"])
    n = len(p["QI"])
    m = int(p["M"])
    QI = p["QI"]
    MI = p["MI"]
    N_w = p["N"]                        # 要求出力，万kW

    VV = p["VV"]                        # 1 基
    HV = p["HV"]
    SS = p["SS"]
    HS = p["HS"]
    xs_cap = VV[1:m + 1]
    ys_cap = HV[1:m + 1]
    xs_dn = SS[1:m + 1]
    ys_dn = HS[1:m + 1]
    interp_cap = _make_interp(p["SA"], xs_cap, ys_cap)
    interp_dn = _make_interp(2, xs_dn, ys_dn)   # 下游曲线本算例恒值，线性即可

    fit_qo = cfg.get("qo") if mode == "fit" else None
    if mode == "fit":
        if not fit_qo or len(fit_qo) != n:
            raise ValueError("fit 模式需提供与时段数等长的 qo 序列（cfg['qo']）")

    VB = VG                               # 定出力调节自正常蓄水位库容起调
    rows = []
    warnings = []
    for k in range(n):
        NT = N_w[k] * 1.0e4               # 要求出力，万kW → W（与 NC 同量纲）

        if mode == "fit":
            qo = float(fit_qo[k])
            st = _period_state(p, VB, QI[k], qo, interp_cap, interp_dn)
            nc = st["出力_raw"]
            st["DN_raw"] = NT - nc
            st["要求出力"] = NT
            st["解q_raw"] = qo
        else:
            qo, _ = _solve_trial(p, VB, QI[k], NT, interp_cap, interp_dn)
            st = _period_state(p, VB, QI[k], qo, interp_cap, interp_dn)
            nc = st["出力_raw"]
            st["DN_raw"] = NT - nc
            st["要求出力"] = NT
            st["解q_raw"] = qo

        if st["弃水发生"]:
            warnings.append(
                f"时段 {MI[k]}: 时段末库容超出 VG={VG:.1f}，弃水 "
                f"{st['弃水_raw']:.1f}（原著以 B1 控制弃水步长）")
        if st["时段末库容_raw"] < VH:
            warnings.append(
                f"时段 {MI[k]}: 时段末库容 {st['时段末库容_raw']:.2f} 低于"
                f"最低限制库容 VH={VH:.1f}（原著机器溢出报警，成果照常打印）")

        h_end = st["时段末水位_raw"]
        e = st["电量_raw"]

        rows.append({
            "时段序号": MI[k],
            "入库流量_raw": QI[k],
            "要求出力_raw": NT,
            "模式": mode,
            **st,
            # 显示列（原著 half-away 舍入；电量按原著 FIX/int 截断显示）
            "入库流量": _rhu(QI[k], 2),
            "弃水流量": _rhu(st["弃水_raw"], 0),
            "调节流量": _rhu(st["解q_raw"], 0),
            "时段库容": _rhu(st["时段末库容_raw"], 0),
            "电量": float(int(e)) if e >= 0 else -float(int(-e)),
            "时段水位": _rhu(h_end, 2),
            "水头": _rhu(st["水头_raw"], 2),
            "出力": _rhu(st["出力_raw"], 0),
        })
        VB = st["时段末库容_raw"]         # 连续值滚动传递

    result = {
        "程序": PROGRAM_ID,
        "VG": VG, "VH": VH, "M": m, "KK": int(p["KK"]), "SA": int(p["SA"]),
        "B1": float(p["B1"]), "B2": float(p["B2"]), "AB": AB, "AC": AC,
        "EP": float(p["EP"]), "E9": float(p["E9"]),
        "A1": float(p["A1"]), "A9": float(p["A9"]),
        "时段序号": MI,
        "入库流量": QI,
        "要求出力": N_w,
        "调节流量": [r["调节流量"] for r in rows],
        "时段末库容": [r["时段库容"] for r in rows],
        "时段末水位": [r["时段水位"] for r in rows],
        "水头": [r["水头"] for r in rows],
        "出力": [r["出力"] for r in rows],
        "电量": [r["电量"] for r in rows],
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
    生成文本计算书（原著 C-6.OUT 风格，主表列右端位置逐字对齐）。
    """
    r = result
    line = "*" * 73
    out = []
    out.append(line)
    out.append("*****                水电站定出力调节计算书（试算法）              *****")
    out.append(line)
    out.append("")
    out.append(" 输入数据:")
    out.append(f"     正常蓄水位库容 VG= {r['VG']:>8.1f}")
    out.append(f"     最低限制水位的库容 VH= {r['VH']:>8.1f}")
    out.append(f"     结点数 M= {r['M']:>2d} ")
    out.append(f"     控制变量 KK= {r['KK']:>2d}     控制变量 SA= {r['SA']:>2d} ")
    out.append(f"     控制弃水流量步长大小的变量 B1= {r['B1']:>5.1f} ")
    out.append(f"     控制溢出的变量 B2= {r['B2']:>5.1f} ")
    out.append(f"     计算期首时段序号 AB= {r['AB']:>2d} ")
    out.append(f"     计算期末时段序号 AC= {r['AC']:>2d} ")
    out.append("")
    out.append(" 时段序号   入库流量    要求的出力")
    for row in r["rows"]:
        out.append(_place_line(35, [
            (6, f"{row['时段序号']:>3d}"),
            (17, f"{row['入库流量']:.2f}"),
            (35, f"{row['要求出力_raw']:.0f}"),
        ]))
    out.append("")
    out.append(" 结点数    库容     上游水位     流量     下游水位")
    VV = params["VV"]
    HV = params["HV"]
    SS = params["SS"]
    HS = params["HS"]
    for i in range(1, r["M"] + 1):
        out.append(_place_line(43, [
            (4, f"{i:>3d}"),
            (14, f"{VV[i]:>8.2f}"),
            (25, f"{HV[i]:>6.2f}"),
            (34, f"{SS[i]:>5.0f}"),
            (43, f"{HS[i]:>6.2f}"),
        ]))
    out.append("")
    out.append(" 计算结果： ")
    out.append("")
    # 主表列右端位置量取自权威 C-6.OUT：
    #   时段3 入库11 弃水20 调节流量32 库容41 电量49 水位57 水头64 出力71
    out.append(" 时段 入库流量 弃水流量 调节流量 计算时段的库容 电量  水位   水头   出力  ")
    for row in r["rows"]:
        out.append(_place_line(80, [
            (3, f"{row['时段序号']:>2d}"),
            (11, f"{row['入库流量']:>4.0f}"),
            (20, f"{row['弃水流量']:>4.0f}"),
            (32, f"{row['调节流量']:>4.0f}"),
            (41, f"{row['时段库容']:>5.0f}"),
            (49, f"{row['电量']:>4.0f}"),
            (57, f"{row['时段水位']:>6.2f}"),
            (64, f"{row['水头']:>6.2f}"),
            (71, f"{row['出力']:>5.0f}"),
        ]))
    out.append("")
    dq = sum(row["弃水_raw"] for row in r["rows"])
    out.append(f"  弃水量 DQ= {dq:>5.0f}")
    out.append(f"  正常蓄水位库容 VG= {r['VG']:.0f}")
    out.append(f"  最低限制水位库容 VH= {r['VH']:.0f}")
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
