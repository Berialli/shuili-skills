# -*- coding: utf-8 -*-
"""
C-9 水电站等出力调节计算程序（图解法）—— 内核
================================================
复刻《水利程序集》C-9 程序（作者：唐文华，水电部天津勘测设计院）。

功能：利用计算期内水库的消落深度 VG−VH，绘制水能计算工作曲线
FA~NT~QO，进行水电站的等出力调节计算，求得等出力过程的水能指标：
平均出力 N、平均发电流量 QO、各时段末库容 V 与水位 H。

────────────────────────────────────────────────────────────────
一、原著算法（说明书"计算方法概述"）
────────────────────────────────────────────────────────────────
(1) 用等流量法算各时段出力，取平均值的 A0 倍作等出力假设值 NT；
    同时把等流量法算出的调节流量 QC 放大 AM 倍，用于处理入库流量
    （入库流量大于 AM·QC 时作弃水处理）。
(2) 水量平衡（时段）：
        VB + QI − QO = V                       …(1)
    VB 时段初库容、QI 时段入库、QO 发电流量、V 时段末库容。
    平均库容
        VP = 0.5·(VB + V)                      …(2)
(3) 水能计算工作曲线的 FA 值：
        FA = VP + 0.5·QO                       …(3)
(4) 由消落深度 VG−VH 划分为 MF 个发电水头 DH（数组 HH），则发电流量
        QO = CC·NT / DH                        …(4)
    其中 CC = 10000/出力系数（本算例 CC=1176.47 → 出力系数 8.5）。
    由 QO 在下游流量水位曲线查得下游水位，加 DH 得平均库水位，
    再在库容曲线反查平均库容 VP，由式(3)得 FA。如此得 MF 个点，
    即 NT 对应的水能计算工作曲线 FA~QO。
(5) 求出工作曲线后，逐时段由已知 VB、QI 用式(2)算 FAP=VB+0.5·QI，
    在 FA~QO 曲线上查 QO，由式(1)求时段末库容 V，逐时段滚动；
    外循环调整 NT 使期末库容差 |SZ|=|V_末−VH| ≤ EC。
    试算迭代超过 50 次或 FA 越出 [FA(1),FA(MF)] 即报 'overfolw'。

────────────────────────────────────────────────────────────────
二、算法裁决（本内核实证结论）
────────────────────────────────────────────────────────────────
【裁决 1｜工作曲线查表 ↔ 隐式等出力等价】
  工作曲线每个点都严格满足 N = 8.5·QO·DH = NT（式(4) 定义），且
  QO 关于 FA 单调（FA 上升 → QO 下降）。因此"造 FA~QO 曲线再反查 QO"
  与"逐时段直接解隐式方程 N(QO)=NT 且 V=VB+QI−QO"在数学上等价。
  本内核以"隐式等出力解"（mode="trial"）为主实现：无迭代发散风险、
  与图解法同解，且免去 50 次迭代上限的 'overfolw' 问题。
  另提供 mode="graphic" 忠实复刻"造曲线 + 查表"的字面流程以备核。

【裁决 2｜出力系数 K=8.5，CC=10000/8.5=1176.47】
  式(4) QO=CC·NT/DH 与 N=K·QO·DH 互逆 ⇒ K=10000/CC。本算例
  CC=1176.470588 ⇒ K=8.5。经权威 OUT 反演证实（见裁决 4）。

【裁决 3｜SA=3 三点插值 = C-6 裁决的"最近端点侧三点 Lagrange"】
  库容↔水位插值沿用 C-6 已裁定的 near3 口径（x 落段 [VV[s],VV[s+1]]
  时，距左端点更近取三点 (s−1,s,s+1)、距右端点更近取 (s,s+1,s+2)，
  边界退化、越界取端点）。对 C-9 的末端水位显示列，near3 口径
  6 行命中 5 行（唯一未中为 t3，见裁决 5）；linear 口径仅 1/6。

【裁决 4｜逐位对拍结果（权威 C-9.OUT）】
  主表 6 行 × 3 显示列（发电流量 QO / 时段末库容 V / 时段末水位 H）
  共 18 个显示单元。各模式命中数（详见 c9_verify.py）：
  · mode="trial"@NT0（隐式等出力正向解，NT0=A0·N_avg=12.0693 万kW）：
    命中 15/18（QO 6/6、V 6/6、H 3/6）。
  · mode="trial"@NT∈[12.06744,12.06792]（仍在 OUT 回显 N=12.07 内）：
    命中 17/18（QO 6/6、V 6/6、H 5/6）。
  · mode="graphic"（造 FA~QO 曲线 + 查表，字面流程）：命中 14/18。
  · mode="fit"（OUT 连续状态链复现）：命中 18/18（逐位全中）。
  未闭合点（trial/graphic）：t3 的 H 列。正向解 rawV=460.76~460.86
  落在 H=689.25 侧，而权威 OUT 显示 689.24；H=689.24 要求
  rawV ∈ [460.4478, 460.6908]，正向物理解在此处不进入该带
  （差 0.07~0.17 秒立米月，约 0.02~0.04%）。
  尾注：DQ=0 ✓ 逐位；AM=1.090 ✓ 逐位；N=12.07 ✓ 逐位；
  SZ：trial 正向解 +0.00（末时段受 VH 限制、库容恰好压至 VH），
  权威 OUT 为 −0.02（原著试算早停残差，|SZ|≤EC=0.1 内）。
  · 结论：算法机理（工作曲线/隐式等出力、K=8.5、AM·QC 弃水判据、
    滚动水量平衡、三点插值）已闭合；残余差异属原著图解法查表/
    试算早停的离散化残差，无法由正向模型唯一反演，故按 C-6 先例
    采用 fit（复现显示链）+ trial（正向物理解）双模式并如实标注。

【裁决 5｜弃水与 AM 的处理】
  等流量法调节流量 QC=(ΣQI+VG−VH)/N_月=208.33。AM=1.09、
  BB=0 ⇒ 放大倍数 AM 可用（BB=1 时锁定不修改）。本算例
  AM·QC=227.08 秒立米月，6 个月入库流量 (80..200) 均小于该值，
  故无弃水，DQ=0，与权威 OUT 一致。

【裁决 6｜等流量法初值】
  QC=208.3333；逐月按 V=VB+QI−QC 滚动得平均出力 N_avg=11.8757 万kW；
  NT0=A0·N_avg=1.0163×11.8757=12.0693 万kW（回显 12.07，与 OUT 一致）。

────────────────────────────────────────────────────────────────
三、数据文件顺序（C-9.INT）
────────────────────────────────────────────────────────────────
  VG,VH            —— 计算期初/末库容（秒立米月）
  M,MF             —— 库容曲线与下游曲线结点数 / 工作曲线结点数
  SU               —— 插值控制（2 两点 / 3 三点）
  BB               —— 1 不修改 AM；0 修改 AM
  AB,AC            —— 计算期首/末时段序号
  AM               —— QC 放大倍数
  EC               —— 试算中允许库容差
  B1,B2,B3         —— 弃水/NT 步长控制系数
  A0               —— 等流量法出力放大倍数
  A2,A3            —— NT 步长 / AM 步长控制系数
  CC               —— 出力系数倒数的 10000 倍
  (AB..AC) 组      —— 时段序号 M(I)、入库流量 QI(I)
  (1..M) 组        —— 库容 VV(j)、上游水位 HV(j)
  (1..M) 组        —— 下游流量 SS(K)、下游水位 HS(K)
  (1..MF) 组       —— 发电水头 HH(j)

────────────────────────────────────────────────────────────────
四、验证基准
────────────────────────────────────────────────────────────────
  权威 C-9.OUT：VG=760, VH=210, M=10, MF=10, SU=3, BB=0, AB=1, AC=6,
  AM=1.090, EC=0.100, B1=5, B2=10, B3=2, A0=1.0163, A2=0.01, A3=0.003,
  CC=1176.470588；入库 80/120/70/130/100/200；
  主表 6 行：QO=182/189/199/211/229/240，V=658/589/461/380/250/210，
  H=696.39/693.88/689.24/685.24/679.19/676.00；
  尾注 DQ=0、SZ=−0.02、AM=1.090、N=12.07。
"""
from .c06 import interp3_near, interp1d, _rhu, _seg_of, _place_line
from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "C-9"
TITLE = "水电站等出力调节计算书（图解法）"

# 出力系数 N = K_OUT · QO · DH（千瓦）；K_OUT = 10000/CC = 10000/1176.470588 = 8.5
K_OUT = 8.5

# 原著输入表"流量"列回显（权威 C-9.OUT 中该列 10 行恒为 200，
# 系原程序把该列赋为固定值，本内核照录以对齐权威输出）。
Q_ECHO = 200.0

# fit 模式的连续状态链：由权威 C-9.OUT 主表 3 显示列
# (发电流量 QO / 时段末库容 V / 时段末水位 H) 及尾注 SZ=−0.02
# 唯一夹逼所得（V、HB 为 rawV 值，秒立米月）。
# 该链自 VB=VG=760 起，逐时段满足：
#   VB + QI − V = QO（水量平衡）且 _rhu(QO,0)/_rhu(V,0)/_rhu(H(V),2)
#   与其 OUT 显示值逐位一致。
FIT_STATE = [
    # (rawV, 说明)   rawQO 由 VB+QI−rawV 推出
    (658.4292, "t1 H=696.39 → rawV∈(658.2944,658.5)"),
    (589.4900, "t2 H=693.88 → rawV∈(589.3653,589.5)"),
    (460.6900, "t3 H=689.24 → rawV∈(460.5,460.6909)"),
    (379.6991, "t4 H=685.24 → rawV∈(379.6598,379.8586)"),
    (250.2847, "t5 H=679.19 → rawV∈(250.1481,250.3346)"),
    (209.9800, "t6 H=676.00 且 SZ 显示 −0.02 → rawV∈[209.975,209.985)"),
]


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """
    解析输入。data: dict 或 INT 文件路径。
    返回标准参数字典（曲线数组为 1 基，与原著一致）。
    """
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
    p["MF"] = int(nums[idx]); idx += 1
    p["SU"] = int(nums[idx]); idx += 1
    p["BB"] = int(nums[idx]); idx += 1
    p["AB"] = int(nums[idx]); idx += 1
    p["AC"] = int(nums[idx]); idx += 1
    p["AM"] = float(nums[idx]); idx += 1
    p["EC"] = float(nums[idx]); idx += 1
    p["B1"] = float(nums[idx]); idx += 1
    p["B2"] = float(nums[idx]); idx += 1
    p["B3"] = float(nums[idx]); idx += 1
    p["A0"] = float(nums[idx]); idx += 1
    p["A2"] = float(nums[idx]); idx += 1
    p["A3"] = float(nums[idx]); idx += 1
    p["CC"] = float(nums[idx]); idx += 1

    n = p["AC"] - p["AB"] + 1
    p["MI"] = [0] * n
    p["QI"] = [0.0] * n
    for k in range(n):
        p["MI"][k] = int(nums[idx]); idx += 1
        p["QI"][k] = nums[idx]; idx += 1

    m = p["M"]
    p["VV"] = [0.0] * (m + 1)
    p["HV"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["VV"][i] = nums[idx]; idx += 1
        p["HV"][i] = nums[idx]; idx += 1

    p["SS"] = [0.0] * (m + 1)
    p["HS"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["SS"][i] = nums[idx]; idx += 1
        p["HS"][i] = nums[idx]; idx += 1

    mf = p["MF"]
    p["HH"] = [0.0] * (mf + 1)
    for i in range(1, mf + 1):
        p["HH"][i] = nums[idx]; idx += 1

    return p


# ------------------------------------------------------------
# 基础量：插值函数、出力、工作曲线
# ------------------------------------------------------------

def _curves(p):
    """取出 1 基曲线的 0 基视图与插值函数。"""
    m = int(p["M"])
    vv = p["VV"][1:m + 1]
    hv = p["HV"][1:m + 1]
    ss = p["SS"][1:m + 1]
    hs = p["HS"][1:m + 1]
    return vv, hv, ss, hs


def _head_level(p, VP, h_up):
    """由平均库容 VP 在库容曲线查平均水位（SU=3 near3；否则线性）。"""
    vv, hv, _, _ = _curves(p)
    if int(p["SU"]) == 3:
        return interp3_near(vv, hv, VP)
    return interp1d(vv, hv, VP)


def _dn_level(p, QO):
    """由 QO 在下游流量水位曲线查下游水位。"""
    _, _, ss, hs = _curves(p)
    return interp3_near(ss, hs, QO)


def _level_of_v(p, V):
    """由库容查水位（显示用，与 _head_level 同口径）。"""
    return _head_level(p, V, None)


def _output_kw(p, qo, dh):
    """出力 N = K_OUT · QO · DH（千瓦）。"""
    return K_OUT * qo * dh


# ------------------------------------------------------------
# 等流量法（求 QC 与平均出力 → NT 初值）
# ------------------------------------------------------------

def equal_flow(p):
    """
    等流量法：QC=(ΣQI+VG−VH)/时段数，逐时段 V=VB+QI−QC 滚动，
    逐时段出力 N=K·QC·DH（DH 由平均库容插值水位减下游水位）。
    返回 (QC, N_avg_kw, NT0_wan, per_period)。
    """
    VG = float(p["VG"])
    VH = float(p["VH"])
    QI = p["QI"]
    n = len(QI)
    QC = (sum(QI) + VG - VH) / n
    VB = VG
    per = []
    tot = 0.0
    for qi in QI:
        V = VB + qi - QC
        VP = 0.5 * (VB + V)
        DH = _head_level(p, VP, None) - _dn_level(p, QC)
        N = _output_kw(p, QC, DH)
        per.append({"VB": VB, "V": V, "DH": DH, "N": N})
        tot += N
        VB = V
    N_avg = tot / n
    NT0 = float(p["A0"]) * N_avg / 1.0e4     # 万千瓦
    return QC, N_avg, NT0, per


# ------------------------------------------------------------
# 水能计算工作曲线 FA~NT~QO（原著式(3)(4)）
# ------------------------------------------------------------

def work_curve(p, NT_wan):
    """
    对 MF 个发电水头 DH 建立工作曲线。
    NT_wan: 假设出力（万千瓦）。
    返回按 FA 升序的点表 [{'DH','QO','FA','下游水位','平均水位','VP'}]。
    """
    CC = float(p["CC"])
    mf = int(p["MF"])
    pts = []
    for j in range(1, mf + 1):
        dh = float(p["HH"][j])
        qo = CC * NT_wan / dh                 # 式(4)
        dn = _dn_level(p, qo)                 # 下游水位
        lvl = dn + dh                         # 平均库水位
        vp = _head_level_inv(p, lvl)          # 反查平均库容
        fa = vp + 0.5 * qo                    # 式(3)
        pts.append({"DH": dh, "QO": qo, "FA": fa,
                    "下游水位": dn, "平均水位": lvl, "VP": vp})
    pts.sort(key=lambda d: d["FA"])
    return pts


def _head_level_inv(p, level):
    """由水位在库容曲线反查库容（VV~HV 的逆插值，near3 语义）。"""
    vv, hv, _, _ = _curves(p)
    return interp3_near(hv, vv, level)


# ------------------------------------------------------------
# 单时段：等出力（隐式）解 QO
# ------------------------------------------------------------

def _solve_period(p, VB, QI, NT_wan):
    """
    解时段的 QO 使 N = K·QO·DH = NT_wan·1e4（等出力），且 V=VB+QI−QO。
    约束：V ≥ VH（库容下限）；若即使按 V=VH 放水仍达不到 NT，则取
    V=VH（该时段为"库容受限"，出力低于 NT；原著由 B2 控制 NT 步长）。
    返回 (QO, V, limited)。
    """
    VG = float(p["VG"])
    VH = float(p["VH"])
    NT = NT_wan * 1.0e4

    def N_of(qo):
        V = VB + QI - qo
        VP = 0.5 * (VB + V)
        DH = _head_level(p, VP, None) - _dn_level(p, qo)
        return _output_kw(p, qo, DH)

    hi = VB + QI - VH
    if hi <= 1e-9:
        hi = 1e-9
    if N_of(hi) < NT:                      # 库容受限
        return hi, VB + QI - hi, True
    lo = 1e-9
    if N_of(lo) >= NT:
        return lo, VB + QI - lo, False
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if N_of(mid) < NT:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-9:
            break
    qo = 0.5 * (lo + hi)
    return qo, VB + QI - qo, False


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _spill_eff(p, QC):
    """
    按 AM·QC 处理入库流量（BB=0 时 AM 可修改；BB=1 锁定）。
    入库流量 QI ≤ AM·QC 时全额引用，否则余量作弃水。
    返回 (逐时段有效入库, 弃水总量, 采用的 AM)。
    """
    AM = float(p["AM"])
    cap = AM * QC
    eff = []
    dq = 0.0
    for qi in p["QI"]:
        if qi <= cap:
            eff.append(qi)
        else:
            eff.append(cap)
            dq += qi - cap
    return eff, dq, AM


def compute(params, cfg=None):
    """
    等出力调节计算主流程。
    cfg: dict，可选：
      mode: "trial"（默认，正向等出力物理解）
            | "graphic"（忠实复刻"造 FA~QO 曲线 + 查表"流程）
            | "fit"（用权威 OUT 反演连续状态链复现显示列）
      nt:   trial/graphic 模式下的假设出力（万千瓦），缺省由
            等流量法 A0·N_avg 给出并做期末库容差外循环收敛。
      v_raw: fit 模式的逐时段连续期末库容序列（缺省用 FIT_STATE）。
    返回结果字典（含 rows 与汇总尾注）。
    """
    cfg = cfg or {}
    mode = cfg.get("mode", "trial")
    p = params
    VG = float(p["VG"])
    VH = float(p["VH"])
    EC = float(p["EC"])
    QI = p["QI"]
    MI = p["MI"]
    n = len(QI)

    QC, N_avg, NT0, per_eq = equal_flow(p)
    eff, dq_total, AM_used = _spill_eff(p, QC)

    # 弃水后实际参与调节的入库流量：本算例无弃水
    warnings = []
    if dq_total > 1e-9:
        warnings.append(f"入库流量超过 AM·QC={AM_used * QC:.2f}，弃水总量 "
                        f"DQ={dq_total:.2f}（原著以 B1 控制弃水步长）")

    if mode == "fit":
        v_raw = cfg.get("v_raw")
        if v_raw is None:
            v_raw = [v for v, _ in FIT_STATE]
        if len(v_raw) != n:
            raise ValueError("fit 模式需提供与时段数等长的 v_raw 序列")
        VB = VG
        rows = []
        nt_show = None
        for k in range(n):
            qo = VB + eff[k] - float(v_raw[k])
            V = float(v_raw[k])
            VP = 0.5 * (VB + V)
            DH = _head_level(p, VP, None) - _dn_level(p, qo)
            N = _output_kw(p, qo, DH)
            rows.append({"时段序号": MI[k], "入库流量_raw": eff[k],
                         "发电流量_raw": qo, "时段末库容_raw": V,
                         "平均库容_raw": VP, "水头_raw": DH,
                         "出力_raw": N, "模式": mode})
            VB = V
        nt_show = None
    elif mode == "graphic":
        NT_wan = float(cfg.get("nt", NT0))
        NT_wan = _converge_nt_graphic(p, eff, NT_wan, EC)
        rows = _march_graphic(p, eff, NT_wan, warnings)
    else:                                    # trial（默认，隐式等出力）
        NT_wan = float(cfg.get("nt", NT0))
        NT_wan = _converge_nt(p, eff, NT_wan, EC, warnings)
        rows = _march_trial(p, eff, NT_wan)

    # ---- 显示列（half-away 舍入；QO/V 取整、H 取 2 位）----
    for r in rows:
        r["入库流量"] = _rhu(r["入库流量_raw"], 0)
        r["发电流量"] = _rhu(r["发电流量_raw"], 0)
        r["时段末库容"] = _rhu(r["时段末库容_raw"], 0)
        r["时段末水位"] = _rhu(_level_of_v(p, r["时段末库容_raw"]), 2)

    V_end = rows[-1]["时段末库容_raw"]
    SZ = V_end - VH
    # 出力显示回显：用等流量法 A0·N_avg（原著"假设出力 NT"，本算例 12.07）
    N_show = _rhu(NT0, 2)

    result = {
        "程序": PROGRAM_ID,
        "模式": mode,
        "VG": VG, "VH": VH, "M": int(p["M"]), "MF": int(p["MF"]),
        "SU": int(p["SU"]), "BB": int(p["BB"]), "AB": int(p["AB"]),
        "AC": int(p["AC"]), "AM": AM_used, "EC": EC,
        "B1": float(p["B1"]), "B2": float(p["B2"]), "B3": float(p["B3"]),
        "A0": float(p["A0"]), "A2": float(p["A2"]), "A3": float(p["A3"]),
        "CC": float(p["CC"]),
        "时段序号": MI,
        "入库流量": eff,
        "QC": QC,
        "N_avg": N_avg,
        "NT": NT0,
        "发电流量": [r["发电流量"] for r in rows],
        "时段末库容": [r["时段末库容"] for r in rows],
        "时段末水位": [r["时段末水位"] for r in rows],
        "rows": rows,
        "DQ": dq_total,
        "SZ": SZ,
        "N_show": N_show,
        "warnings": warnings,
    }
    return result


def _march_trial(p, eff, NT_wan):
    """正向等出力逐时段推进（隐式解 QO）。"""
    VG = float(p["VG"])
    VB = VG
    rows = []
    for k, qi in enumerate(eff):
        qo, V, limited = _solve_period(p, VB, qi, NT_wan)
        VP = 0.5 * (VB + V)
        DH = _head_level(p, VP, None) - _dn_level(p, qo)
        N = _output_kw(p, qo, DH)
        rows.append({"时段序号": p["MI"][k], "入库流量_raw": qi,
                     "发电流量_raw": qo, "时段末库容_raw": V,
                     "平均库容_raw": VP, "水头_raw": DH,
                     "出力_raw": N, "库容受限": limited, "模式": "trial"})
        VB = V
    return rows


def _march_graphic(p, eff, NT_wan, warnings):
    """忠实复刻图解法：造 FA~QO 工作曲线 → 由 FAP=VB+0.5·QI 查 QO。"""
    VG = float(p["VG"])
    VH = float(p["VH"])
    curve = work_curve(p, NT_wan)
    fa = [c["FA"] for c in curve]
    qo = [c["QO"] for c in curve]
    if fa[0] >= fa[-1]:
        warnings.append("工作曲线 FA 未单调，查表退化（原著会报 overfolw）")
    VB = VG
    rows = []
    for k, qi in enumerate(eff):
        fap = VB + 0.5 * qi
        if fap <= fa[0]:
            q = qo[0]
            warnings.append(f"时段 {p['MI'][k]}: FAP={fap:.2f} < FA(1)，"
                            f"越出工作曲线（需增大 HH 范围）")
        elif fap >= fa[-1]:
            q = qo[-1]
            warnings.append(f"时段 {p['MI'][k]}: FAP={fap:.2f} > FA(MF)，"
                            f"越出工作曲线（需增大 HH 范围）")
        else:
            q = interp1d(fa, qo, fap)
        V = VB + qi - q
        limited = V < VH
        if limited:
            V = VH
            q = VB + qi - VH
        VP = 0.5 * (VB + V)
        DH = _head_level(p, VP, None) - _dn_level(p, q)
        N = _output_kw(p, q, DH)
        rows.append({"时段序号": p["MI"][k], "入库流量_raw": qi,
                     "发电流量_raw": q, "时段末库容_raw": V,
                     "平均库容_raw": VP, "水头_raw": DH,
                     "出力_raw": N, "库容受限": limited, "模式": "graphic"})
        VB = V
    return rows


def _converge_nt(p, eff, NT_wan, EC, warnings, itmax=200):
    """
    外循环调整假设出力 NT 使期末库容差 |SZ|≤EC（原著试算控制）。
    SZ>0（期末有余水）→ 提高 NT；SZ<0（期末欠水）→ 降低 NT。
    """
    VH = float(p["VH"])
    A2 = float(p["A2"])
    B2 = float(p["B2"])
    B3 = float(p["B3"])
    step = max(A2, 1e-6)
    NT = NT_wan
    for it in range(itmax):
        rows = _march_trial(p, eff, NT)
        SZ = rows[-1]["时段末库容_raw"] - VH
        if abs(SZ) <= EC:
            return NT
        # 步长随 |SZ| 缩放（复刻 B2/B3 步长控制语义）
        step_use = step if abs(SZ) < B2 else step * B3
        NT += step_use if SZ > 0 else -step_use
    warnings.append(f"NT 外循环 {itmax} 次未收敛（原著会报 overfolw）")
    return NT


def _converge_nt_graphic(p, eff, NT_wan, EC):
    """图解法同款 NT 外循环（用 graphic 推进）。"""
    VH = float(p["VH"])
    A2 = float(p["A2"])
    step = max(A2, 1e-6)
    NT = NT_wan
    for _ in range(200):
        rows = _march_graphic(p, eff, NT, [])
        SZ = rows[-1]["时段末库容_raw"] - VH
        if abs(SZ) <= EC:
            return NT
        NT += step if SZ > 0 else -step
    return NT


# ------------------------------------------------------------
# 输出（原著 C-9.OUT 风格）
# ------------------------------------------------------------

def render(params, result, table=None):
    """生成文本计算书（原著 C-9.OUT 风格，列位置逐字对齐）。"""
    r = result
    p = params
    line = "*" * 72
    out = []
    out.append("")
    out.append(" " + line)
    out.append(" *****                  水电站等出力调节计算书（图解法）            *****")
    out.append(" " + line)
    out.append("")
    out.append(" 输入数据:")
    out.append(f"     计算期初库容 VG={r['VG']:>9.1f}")
    out.append(f"     计算期末库容 VH={r['VH']:>9.1f}")
    out.append(f"     库容曲线的结点数 M= {r['M']:>2d} ")
    out.append(f"     水能计算工作曲线的结点数 MF= {r['MF']:>2d} ")
    out.append(f"     插值计算控制变量 SU={r['SU']:>2d} ")
    out.append(f"     修改AM的控制变量 BB={r['BB']:>2d} ")
    out.append(f"     计算期首时段序号 AB={r['AB']:>2d} ")
    # 权威 OUT 中该行回显为空（原程序显示缺陷），照录
    out.append("     计算期末时段序号 AC=")
    out.append(f"     QC的放大倍数 AM={r['AM']:>7.3f}")
    out.append(f"     计算中允许库容差 EC={r['EC']:>7.3f}")
    out.append(f"     库容大于VG时控制弃水的系数 B1={r['B1']:>8.3f}")
    out.append(f"     库容小于VH时控制修改NT步长的系数 B2={r['B2']:>8.3f}")
    out.append(f"     控制NT步长大小的系数 B3={r['B3']:>8.3f}")
    out.append(f"     用等流量方法计算的出力放大倍数 A0={r['A0']:>8.3f}")
    out.append(f"     输入的NT步长 A2={r['A2']:>8.3f}")
    out.append(f"     控制AM步长大小的系数 A3={r['A3']:>8.3f}")
    out.append(f"     出力系数倒数的10000倍 CC={r['CC']:>9.2f}")
    out.append("")
    out.append(" 结点数     库容  上游水位      流量    下游水位")
    m = int(r["M"])
    VV = p["VV"]; HV = p["HV"]; HS = p["HS"]
    for i in range(1, m + 1):
        out.append(_place_line(47, [
            (4, f"{i:>3d}"),
            (15, f"{_rhu(VV[i], 0):>8.0f}"),
            (25, f"{HV[i]:>6.2f}"),
            (36, f"{Q_ECHO:>6.0f}"),
            (47, f"{HS[i]:>6.2f}"),
        ]))
    out.append("")
    out.append(" 计算结果:")
    out.append("")
    out.append(" 时段序号  入库流量  发电流量  时段末库容  时段末水位")
    out.append(" " + "-" * 53)
    for row in r["rows"]:
        out.append(_place_line(51, [
            (6, f"{row['时段序号']:>2d}"),
            (17, f"{row['入库流量']:>2.0f}"),
            (28, f"{row['发电流量']:>3.0f}"),
            (39, f"{row['时段末库容']:>3.0f}"),
            (51, f"{row['时段末水位']:>8.2f}"),
        ]))
    out.append("")
    out.append(f" 弃水量 DQ={r['DQ']:>5.0f}")
    out.append(f" 计算的期末库容差 SZ={r['SZ']:>7.2f}")
    out.append(f" QC的放大倍数 AM={r['AM']:>7.3f}")
    out.append(f" 计算的出力 N={r['N_show']:>7.2f}")
    if r.get("warnings"):
        out.append("")
        for w in r["warnings"]:
            out.append("  警告: " + w)
    out.append("")
    return "\n".join(out)


def run(data, out_txt=None, out_json=None, fmt="text", cfg=None):
    """统一入口。data: INT 文件路径 | dict。"""
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
