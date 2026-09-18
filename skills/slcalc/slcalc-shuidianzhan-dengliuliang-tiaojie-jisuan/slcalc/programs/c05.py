# -*- coding: utf-8 -*-
"""
C-5 水电站蓄水期调节计算程序 —— 内核
=============================================
复刻《水利程序集》C-5 程序（作者：唐文华，水电部天津勘测设计院）。

功能：已知入库流量过程、库容曲线（库容~水位）、下游水位流量曲线、
正常蓄水位相应库容 VG、限制水位相应库容 VT（时段序号≤AF 的最高允许
蓄水位）/VH（最低限制水位）、允许最小调节流量 MQ，进行蓄水期逐时段
水量平衡调节，输出各时段调节流量、时段末库容与水位、出力、电量。

计算方法概要（原著说明书）—— 蓄水期调节规则
  水电站蓄水期的目标：在满足 库容≤V1（V1=VT 当 时段序号≤AF，否则
  =VG）、库容≥VH、调节流量≥MQ 的前提下尽量多蓄水抬高水位。
  第 I 时段水量平衡：QI(I) − QO(I) = V(I) − V(I−1)
  把未知量 QO、V(I) 分别换为 MQ 与 V1 得到判断条件（原著语句 2345）：
       D = QI(I) − MQ − (V1 − V(I−1))
   (i)  D ≤ 0  → 水库能存来水抬高水位：QO = MQ，
                  V(I) = V(I−1) + QI(I) − MQ（≤ V1 自动满足）
   (ii) D > 0  → 库水位已到最高（顶格 V1），来水超过可蓄空间的部分
                  全部用于发电：V(I) = V1，
                  QO(I) = QI(I) − (V1 − V(I−1)) = MQ + D
  蓄水期起始库容取最低限制水位库容 VH（算例闭合验证）。
  若某时段 V(I) < VH，原著以三音响报警并继续打印已算成果（本内核
  以 warnings 列表提示）。

水能计算：
  平均库容 VC = 0.5*(VB + V(I))  → 库容曲线 VV~HV 线性内插得平均水位 HC
  → 下游水位以调节流量 QO 在 SS~HS 下游水位流量曲线上线性内插
  → 水头 DH = HC − 下游水位
  → 时段出力 N = 8.3 * QO * DH   （单位：千瓦，8.3 为出力系数）
  → 电量 E = N * 730 / 1e8        （1 时段=730 小时；单位：亿度）
    （单位说明见原著：库容 秒立米月、流量 秒立米、出力 千瓦、电量 亿度）

数据文件顺序（C-5.INT）：
  VG, VT, VH, M                正常蓄水位/限制/最低限制库容，结点数 M
  MQ, AF                       允许最小调节流量，VT 最高库容的时段最大序号
  AB, AC                       计算期首/末时段序号
  M(I), QI(I)  (AB..AC)        时段序号、入库流量（交错对）
  VV(I), HV(I)  (1..M)         库容曲线结点（库容,水位 交错对）
  SS(I), HS(I)  (1..M)         下游水位流量曲线结点（流量,水位 交错对）

验证基准：C-5.OUT（算例，5 时段）
  VG=247, VT=201.5, VH=99, M=10, MQ=275, AF=4, AB=1..AC=5，
  入库 527/430/550/434/759。时段1 自 VH=99 起蓄：
  D1=252−102.5>0 顶格 VT=201.5 → QO=424.5(显示425)、V=201.5(显示202)，
  时段2~4 库容已满 VT、来水全发 QO=QI，时段5 上限转 VG=247 → QO=713.5
  (显示714)、V=247。主表 5 行 调节流量/库容/水位/出力/电量 逐位命中。
"""
from ..core.intio import read_numbers, smart_read_text
from ..core.outgen import write_out, write_json

PROGRAM_ID = "C-5"
TITLE = "水电站蓄水期调节计算书"


# ------------------------------------------------------------
# 插值
# ------------------------------------------------------------

def interp1d(xs, ys, x):
    """
    分段线性插值，与原著一致（同 C-4）：
      x <= xs(1)   → ys(1)
      x >= xs(M)   → ys(M)
      否则在 [xs(i), xs(i+1)] 段内直线内插。
    xs/ys 为 0 基列表（含全部 M 个结点，升序）。
    """
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    i = 1
    while i < len(xs) - 1 and xs[i] < x:
        i += 1
    t = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
    return ys[i - 1] + t * (ys[i] - ys[i - 1])


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
        p["程序"] = PROGRAM_ID
        return p

    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}
    idx = 0

    p["VG"] = nums[idx]; idx += 1       # 正常蓄水位相应的库容
    p["VT"] = nums[idx]; idx += 1       # 限制时段最高蓄水位相应的库容
    p["VH"] = nums[idx]; idx += 1       # 最低限制水位的库容
    p["M"] = int(nums[idx]); idx += 1   # 曲线结点数 M

    p["MQ"] = nums[idx]; idx += 1       # 蓄水期允许的最小调节流量
    p["AF"] = int(nums[idx]); idx += 1  # 最高库容 VT 的时段最大序号

    p["AB"] = int(nums[idx]); idx += 1  # 计算期首时段序号
    p["AC"] = int(nums[idx]); idx += 1  # 计算期末时段序号

    n = p["AC"] - p["AB"] + 1
    p["MI"] = [0] * n                   # 时段序号（绝对）
    p["QI"] = [0.0] * n                 # 入库流量
    for k in range(n):
        p["MI"][k] = int(nums[idx]); idx += 1
        p["QI"][k] = nums[idx]; idx += 1

    m = p["M"]
    p["VV"] = [0.0] * (m + 1)           # 1 基：库容
    p["HV"] = [0.0] * (m + 1)           # 1 基：水位
    for i in range(1, m + 1):
        p["VV"][i] = nums[idx]; idx += 1
        p["HV"][i] = nums[idx]; idx += 1

    p["SS"] = [0.0] * (m + 1)           # 1 基：下游水位流量曲线
    p["HS"] = [0.0] * (m + 1)
    for i in range(1, m + 1):
        p["SS"][i] = nums[idx]; idx += 1
        p["HS"][i] = nums[idx]; idx += 1

    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params, cfg=None):
    """
    蓄水期调节计算主流程。
    返回结果字典（含逐时段明细表 rows 与汇总）。
    """
    VG = float(params["VG"])
    VT = float(params["VT"])
    VH = float(params["VH"])
    MQ = float(params["MQ"])
    AF = int(params["AF"])
    AB = int(params["AB"])
    AC = int(params["AC"])
    n = len(params["QI"])
    QI = params["QI"]
    MI = params["MI"]
    m = int(params["M"])
    VV = params["VV"]          # 1 基
    HV = params["HV"]
    SS = params["SS"]          # 1 基
    HS = params["HS"]

    A = float(params.get("A", 8.3))      # 出力系数（千瓦/时段流量·米）
    TH = float(params.get("TH", 730.0))  # 1 时段小时数（算例按月 730h）

    xs_cap = VV[1:m + 1]
    ys_cap = HV[1:m + 1]
    xs_dn = SS[1:m + 1]
    ys_dn = HS[1:m + 1]

    VB = VH                              # 蓄水期起始库容 = 最低限制库容
    rows = []
    warnings = []
    for i in range(n):
        # 本时段允许的最高库容：时段序号 ≤ AF 用 VT，其后用 VG
        V1 = VT if MI[i] <= AF else VG

        D = QI[i] - MQ - (V1 - VB)       # 原著语句 2345 判断式
        if D <= 0.0:
            mode = "蓄水"                # (i) 能存来水抬高水位
            qo = MQ
            v = VB + QI[i] - MQ
        else:
            mode = "顶格超发"            # (ii) 库水位已到最高，来水全用于发电
            qo = QI[i] - (V1 - VB)
            v = V1

        if v < VH:
            warnings.append(
                f"时段 {MI[i]}: 时段末库容 {v:.2f} 低于最低限制库容 VH={VH:.2f}"
                f"（原著三音响报警，成果照常打印）")

        VC = 0.5 * (VB + v)              # 平均库容
        HC = interp1d(xs_cap, ys_cap, VC)     # 平均水位
        zdown = interp1d(xs_dn, ys_dn, qo)    # 下游水位（QO 内插）
        DH = HC - zdown                  # 水头
        NC = A * qo * DH                 # 出力（千瓦）
        E = NC * TH / 1.0e8              # 电量（亿度）

        h_end = interp1d(xs_cap, ys_cap, v)   # 时段末水位

        rows.append({
            "时段序号": MI[i],
            "入库流量_raw": QI[i],
            "调节流量_raw": qo,
            "时段末库容_raw": v,
            "时段末水位_raw": h_end,
            "时段初库容_raw": VB,
            "上限库容_raw": V1,
            "判断D_raw": D,
            "模式": mode,
            "平均库容_raw": VC,
            "平均水位_raw": HC,
            "下游水位_raw": zdown,
            "水头_raw": DH,
            "出力_raw": NC,
            "电量_raw": E,
            # 显示列（原著 half-away 舍入）
            "入库流量": round(QI[i], 2),
            "调节流量": _rhu(qo, 0),
            "时段库容": _rhu(v, 0),
            "时段水位": _rhu(h_end, 2),
            "出力": _rhu(NC, 0),
            "电量": _rhu(E, 3),
        })
        VB = v

    result = {
        "程序": PROGRAM_ID,
        "VG": VG, "VT": VT, "VH": VH, "MQ": MQ, "AF": AF,
        "AB": AB, "AC": AC, "M": m,
        "入库流量": QI,
        "调节流量": [r["调节流量"] for r in rows],
        "时段末库容": [r["时段库容"] for r in rows],
        "时段末水位": [r["时段水位"] for r in rows],
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
    """
    按"字段右端对齐"落字构造一行文本。
    fields: [(end_col, text), ...]；text 的最右字符结束于 end_col-1。
    """
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
    生成文本计算书（原著 .OUT 风格）。
    """
    line = "*" * 76
    hdr = [
        line,
        f"{'*':<4}{'水电站蓄水期调节计算书':^64}{'*':>4}",
        line,
        "",
        " 输入数据:",
        f"     正常蓄水位相应的库容 VG= {result['VG']:>8.1f}",
        f"     指定时段以前时段的最高蓄水位相应的库容 VT= {result['VT']:>8.1f}",
        f"     最低限制水位的库容 VH= {result['VH']:>8.1f}",
        f"     结点数 M= {result['M']:>2d} ",
        f"     蓄水期允许的最小调节流量 MQ= {result['MQ']:>8.1f}",
        f"     计算期首时段序号 AB= {result['AB']:>2d} ",
        f"     计算期末时段序号 AC= {result['AC']:>2d} ",
        "",
        "  时段序号     入库流量",
    ]
    for row in result["rows"]:
        hdr.append(_place_line(23, [
            (5, f"{row['时段序号']:>5d}"),
            (22, f"{row['入库流量']:.2f}"),
        ]))
    hdr.append("")
    hdr.append(" 序号I  库容 VV(I)  水位 HV(I)  流量 SS(I)  下游水位 HS(I)")
    for i in range(1, result["M"] + 1):
        hdr.append(_place_line(53, [
            (5, f"{i:>5d}"),
            (17, f"{float(params['VV'][i]):.1f}"),
            (29, f"{params['HV'][i]:.2f}"),
            (41, f"{float(params['SS'][i]):.1f}"),
            (53, f"{params['HS'][i]:.2f}"),
        ]))

    # 主表列右端位置量取自权威 C-5.OUT：
    #   时段 端4  入库流量 端15  调节流量 端25  时段末库容 端38
    #   时段末水位 端51  出力 端64  电量 端75
    body = [" 计算结果: ", ""]
    body.append(" 时段   入库流量  调节流量  时段末库容  时段末水位     出力      电量")
    for row in result["rows"]:
        body.append(_place_line(76, [
            (4, f"{row['时段序号']:>4d}"),
            (15, f"{row['入库流量']:.0f}"),
            (25, f"{row['调节流量']:>4.0f}"),
            (38, f"{row['时段库容']:>8.0f}"),
            (51, f"{row['时段水位']:>8.2f}"),
            (64, f"{row['出力']:>8.0f}"),
            (75, f"{row['电量']:>7.3f}"),
        ]))
    if result.get("warnings"):
        body.append("")
        for w in result["warnings"]:
            body.append("  警告: " + w)
    return "\n".join(hdr + [""] + body + [""])


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
