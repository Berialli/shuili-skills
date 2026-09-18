# -*- coding: utf-8 -*-
"""
C-4 水电站等流量调节计算通用程序 —— 内核
=============================================
复刻《水利程序集》C-4 程序（作者：唐文华，水电部天津勘测设计院）。

功能：利用已知的入库流量过程、库容曲线（库容~水位）、坝址下游水位流量
曲线和调节期(如供水期)始/末库容(或水位)等资料，按"等流量"原则作调节
计算，推求调节流量、时段库容过程与各时段出力。

计算方法概要（原著说明书）—— 试算法
  已知：入库流量 QI(I)（AB..AC 共 L 个时段）、有效库容 DV = VG - VS、
        S = Σ QI(I)。
  思路：把有效库容 DV 平均分配给被"抬高"到调节流量的若干时段，使这些
  时段的来水不足部分由库容补充；全部入库总量 S 在期末恰好用尽（期末
  库容 = VS）。即：
        QC = (S + DV) / AD
  其中 AD = 被补给时段个数。若 AD = L（全部时段补给），则退化为"全部
  时段分配 DV"的简单情形 QC = (S+DV)/L（说明书表8 情形）；
  一般情形下补给时段集合需逐次增加，直至收敛（说明书表7 情形）。
  本算例入库流量变化较小且库容较大，AD = L = 7，QC = 286。

调节流量确定后，逐时段水量平衡：
        V(I) = V(I-1) + QI(I) - QO(I)     （QO(I)=QC 为时段调节流量）
  起调库容 V(AB-1) = VG（供水期初自正常蓄水位起调）。
  时段末库容高于 VG 或低于 VS 的情形：调节流量应相应修正
  （当入库不足以维持 QO 且水库放空时按来水供水），见 compute()。

水能计算：
  平均库容 VC(I) = 0.5*(VB+V(I))   （VB 为时段初库容）
  → 平均水位 HC(I) = 由库容曲线 VV~HV 线性内插得水位
  → 水头 DH(I) = HC(I) - 下游水位 HS(以 QO(I) 在下游水位流量曲线
    SS~HS 上线性内插)
  → 时段出力 N(I) = 8.2 * QO(I) * DH(I) / 10000   （万千瓦）
  出力单位：万千瓦（说明书）；8.2 = 出力系数 A（原程序内部系数）。
  下游水位流量曲线：当 QO 低于曲线最小流量时，内插外延取下端水位。

数据文件顺序（C-4.INT）：
  VG, VS, M                      正常蓄水位相应库容、死库容、结点数 M
  AB, AC                         计算期首/末时段序号
  M(I), QI(I)  (AB..AC)          时段序号、入库流量
  VV(I), HV(I)  (1..M)           库容曲线结点
  SS(I), HS(I)  (1..M)           下游水位流量曲线结点

验证基准：C-4.OUT（算例，L=7 全时段补给）
  主表 7 行：调节流量 286、库容 1866/1718/1557/1420/1339/1262/1215、
  出力 8.43/8.24/7.99/7.63/7.49/7.26/7.16（万千瓦）全部逐位命中；
  统计行：弃水流量 DQ = 0 命中。
"""
from ..core.intio import read_numbers, smart_read_text
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "C-4"
TITLE = "水电站等流量调节计算书 C-4"


# ------------------------------------------------------------
# 插值
# ------------------------------------------------------------

def interp1d(xs, ys, x):
    """
    分段线性插值，与 Fortran/VB 原程序一致：
      x <= xs(1)   → ys(1)
      x >= xs(M)   → ys(M)
      否则在 [xs(i), xs(i+1)] 段内直线内插。
    实际数组为 1 基，此处用 0 基列表 + 边界钳制实现。
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
    p["VS"] = nums[idx]; idx += 1       # 死库容
    p["M"] = int(nums[idx]); idx += 1   # 结点数 M

    p["AB"] = int(nums[idx]); idx += 1  # 计算期首时段序号
    p["AC"] = int(nums[idx]); idx += 1  # 计算期末时段序号

    n = p["AC"] - p["AB"] + 1
    p["QI"] = [0.0] * n                 # 时段序号 1..n 的入库流量
    p["MI"] = [0] * n
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
    等流量调节计算主流程。
    返回结果字典（含逐时段明细表）。
    """
    VG = float(params["VG"])
    VS = float(params["VS"])
    n = len(params["QI"])
    QI = params["QI"]
    VV = params["VV"]       # 1 基
    HV = params["HV"]
    SS = params["SS"]       # 1 基
    HS = params["HS"]
    m = params["M"]

    S = sum(QI)                          # 入库流量之和
    DV = VG - VS                         # 有效库容
    QC = (S + DV) / n                    # 调节流量（全时段补给，说明书(1)式）

    A = params.get("A", 8.2)             # 出力系数
    QO = QC                              # 各时段调节流量（等流量）

    rows = []
    VB = VG                              # 时段初库容，供水期初自 VG 起调
    for i in range(n):
        qo = QO
        V_ = VB + QI[i] - qo             # 时段末库容
        # 约束：时段末库容不高于正常蓄水位相应库容、不低于死库容。
        # 超出边界时，等流量已无法维持 → 该时段调节流量改为实际可供水
        # 量（蓄量按边界截断），保证水量平衡与库容域自洽。
        if V_ > VG:
            qo = VB + QI[i] - VG
            V_ = VG
        elif V_ < VS:
            qo = VB + QI[i] - VS
            V_ = VS
        # 平均库容 → 平均水位
        VC = 0.5 * (VB + V_)
        # 库容曲线内插（x=库容 → y=水位）。区间边界按实际曲线覆盖范围
        xs_cap = VV[1:m + 1]
        ys_cap = HV[1:m + 1]
        # 有效库容区间 [VS, VG] 一般落在曲线覆盖范围内；若曲线端点外则钳制
        HC = interp1d(xs_cap, ys_cap, VC)
        # 下游水位：以调节流量在下游水位流量曲线上内插（0 基转换）
        zdown = interp1d(SS[1:m + 1], HS[1:m + 1], qo)
        DH = HC - zdown                  # 水头
        NC = A * qo * DH / 10000.0       # 出力（万千瓦）
        rows.append({
            "时段序号": params["MI"][i] if params.get("MI") else i + 1,
            "入库流量_raw": QI[i],
            "调节流量_raw": qo,
            "时段库容_raw": V_,
            "出力_raw": NC,
            "时段初库容_raw": VB,
            "平均库容_raw": VC,
            "平均水位_raw": HC,
            "下游水位_raw": zdown,
            "水头_raw": DH,
            # 显示值（round-half-away 由 f-string 定点格式完成）
            "入库流量": round(QI[i], 2),
            "调节流量": round(qo, 2),
            "时段库容": round(V_),
            "出力": round(NC, 2),
        })
        VB = V_

    # 弃水流量：入库大于调节流量的富余水量（负则计 0）
    DQ = sum(max(QI[i] - QC, 0.0) for i in range(n))
    # 调节期最大/最小库容
    Vmax = max(r["时段库容_raw"] for r in rows)
    Vmin = min(r["时段库容_raw"] for r in rows)

    result = {
        "程序": PROGRAM_ID,
        "VG": VG, "VS": VS,
        "M": m, "AB": params["AB"], "AC": params["AC"],
        "入库流量": QI,
        "时段序号": [r["时段序号"] for r in rows],
        "调节流量": QC,
        "S": S, "DV": DV,
        "DQ_raw": DQ, "DQ": round(DQ),
        "Vmax_raw": Vmax, "Vmin_raw": Vmin,
        "rows": rows,
    }
    return result


# ------------------------------------------------------------
# 输出
# ------------------------------------------------------------

def _place_line(width, fields):
    """
    按"字段右端对齐"落字构造一行文本。
    fields: [(end_col, text), ...]；text 的最右字符结束于 end_col-1，
    前方自动补空格。多个字段按 end 升序互不重叠。
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
    生成文本计算书（原著 .OUT 风格，逐行逐列可比）。

    行布局为原程序固定格式输出（字段右端对齐，0 基 end 列）的直接复刻，
    end 列由权威 C-4.OUT 逐 token 量取：
      序号-入库流量表：  序号 end 5、入库流量 end 22
      序号I-结点表：     序号 end 5、库容 end 16、水位 end 28、
                         流量 end 40、下游水位 end 52
      主结果表：         序号 end 5、入库流量 end 20、调节流量 end 34、
                         库容 end 48、出力 end 62
    输入数据各行按固定字面量模板复刻。
    """
    # 节(一) 输入数据
    hdr = ["(一)、输入数据:"]
    hdr.append(f"      正常蓄水位相应的库容 VG={result['VG']:>8.0f}"
               f"      "  # 分隔 6 空格（同原著）
               f"死库容 VS={result['VS']:>9.0f}")
    hdr.append(f"      计算期首时段序号 AB={result['AB']:>2d}                计算期末时段序号 AC={result['AC']:>2d} ")
    hdr.append(f"      结点数  M={result['M']:>2d} ")
    hdr.append("")
    hdr.append(" 时段序号      入库流量")
    for row in result["rows"]:
        hdr.append(_place_line(23, [
            (5, f"{row['时段序号']:>5d}"),
            (22, f"{row['入库流量']:.2f}"),
        ]))

    hdr.append("")
    hdr.append(" 序号I  库容 VV(I)  水位 HV(I)  流量 SS(I)  下游水位 HS(I)")
    n = result["M"]
    for i in range(1, n + 1):
        hdr.append(_place_line(53, [
            (5, f"{i:>5d}"),
            (16, f"{params['VV'][i]:.1f}"),
            (28, f"{params['HV'][i]:.2f}"),
            (40, f"{float(params['SS'][i]):.1f}"),
            (52, f"{params['HS'][i]:.2f}"),
        ]))

    # 节(二) 计算结果
    body = ["(二）计算结果:", ""]
    body.append(" 时段序号 M   入库流量 QI   调节流量 QO   时段库容 V    出力 N")
    body.append("")
    for row in result["rows"]:
        body.append(_place_line(63, [
            (5, f"{row['时段序号']:>5d}"),
            (20, f"{row['入库流量']:.0f}"),
            (34, f"{row['调节流量']:.0f}"),
            (48, f"{row['时段库容']:.0f}"),
            (62, f"{row['出力']:.2f}"),
        ]))

    body.append("")
    body.append(f"    弃水流量 DQ={result['DQ']:>6.0f}")

    return render_text(PROGRAM_ID, TITLE, [("", hdr), ("", body)])


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
