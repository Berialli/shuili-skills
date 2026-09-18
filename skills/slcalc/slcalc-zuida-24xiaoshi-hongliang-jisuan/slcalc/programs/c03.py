# -*- coding: utf-8 -*-
"""
C-3 两种供水保证的水库径流调节计算程序 —— 内核
==============================================
复刻《水利程序集》C-3 程序（作者：王建生，水电部天津勘测设计院）。

原理：时历法径流调节，核心公式：
    V = V0 + R - A - B - E
其中：
    V   ：月末水库蓄量
    V0  ：上月末水库蓄量
    R   ：月径流量
    A   ：高保证率供水月水量（工业用水）
    B   ：低保证率供水月水量（农业用水）
    E   ：水库蒸发渗漏月损失量
若 V < 0，则供水破坏。

验证基准（N=20, VM=7000, DX=50, PA=0.9, PB=0.7, K=1, MO=6）：
    VX=7000, VY=1503: AY=0, AM=0, BY=4, BM=17
    VX=6950, VY=1498: AY=1, AM=1, BY=5, BM=18
    VX=6900, VY=1419: AY=1, AM=2, BY=5, BM=18
    VX=6800, VY=1700: AY=1, AM=1, BY=5, BM=21
    VX=6500, VY=1425: AY=1, AM=2, BY=5, BM=21
    VX=6450, VY=1522: AY=2, AM=6, BY=4, BM=23
"""
from ..core.intio import read_numbers, smart_read_text
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "C-3"
TITLE = "两种供水保证率的水库径流调节计算书 C-3"


def parse(data):
    if isinstance(data, dict):
        return data
    text = smart_read_text(data)
    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}
    idx = 0
    p["N"] = int(nums[idx]); idx += 1
    p["VM"] = nums[idx]; idx += 1
    p["DX"] = nums[idx]; idx += 1
    p["PA"] = nums[idx]; idx += 1
    p["PB"] = nums[idx]; idx += 1
    p["K"] = int(nums[idx]); idx += 1
    p["MO"] = int(nums[idx]); idx += 1
    p["A"] = nums[idx:idx + 12]; idx += 12
    p["E"] = nums[idx:idx + 12]; idx += 12
    years = []
    while idx < len(nums):
        year = int(nums[idx]); idx += 1
        if year == 0 and idx >= len(nums):
            break
        R = []; B = []
        for m in range(12):
            R.append(nums[idx]); idx += 1
            B.append(nums[idx]); idx += 1
        years.append({"year": year, "R": R, "B": B})
    p["years"] = years
    return p


def compute(params, cfg):
    N = params["N"]
    A = params["A"]
    E = params["E"]
    years = params["years"]
    MO = params["MO"]
    VX = cfg["VX"]
    VY = cfg["VY"]

    total_months = N * 12
    R_all = []
    B_all = []
    year_idx = []
    month_idx = []

    for ydata in years:
        year = ydata["year"]
        R = ydata["R"]
        B = ydata["B"]
        for m in range(12):
            adj_m = (m + MO - 1) % 12
            R_all.append(R[m])
            B_all.append(B[m])
            year_idx.append(year)
            month_idx.append(adj_m)

    V = 0.0
    ay_months = []
    am_months = []
    ay_years = set()
    am_years = set()
    storage_points = []

    for i in range(total_months):
        V_new = V + R_all[i] - A[i % 12] - B_all[i] - E[i % 12]
        storage_points.append({
            "年": year_idx[i],
            "月": month_idx[i] + 1,
            "蓄量": round(V_new),
            "径流": round(R_all[i]),
            "高供水": round(A[i % 12]),
            "低供水": round(B_all[i]),
            "损失": round(E[i % 12]),
        })
        if V_new < 0:
            ay_months.append(i)
            ay_years.add(year_idx[i])
            am_months.append(i)
            am_years.add(year_idx[i])
        # 蓄量约束：不能为负，不能超过限制库容
        V_new = max(0.0, min(V_new, float(VY)))
        V = V_new

    AY = len(ay_years)
    AM = len(ay_months)
    BY = len(am_years)
    BM = len(am_months)
    PAY = (N - AY) / N
    PBY = (N - BY) / N
    PAM = (12 * N - AM) / (12 * N)
    PBM = (12 * N - BM) / (12 * N)

    return {
        "AY": AY, "AM": AM, "BY": BY, "BM": BM,
        "PAY": round(PAY, 3), "PBY": round(PBY, 3),
        "PAM": round(PAM, 3), "PBM": round(PBM, 3),
        "storage_points": storage_points,
        "VX": VX, "VY": VY,
    }


def render(params, result, table):
    lines = []
    lines.append(f"径流系列年数 N = {params['N']}")
    lines.append(f"允许最大的兴利库容 VM = {params['VM']:.1f}")
    lines.append(f"兴利库容递减的幅度 DX = {params['DX']}")
    lines.append(f"高保证率(例如工业用水)的供水保证率 PA = {params['PA']:.3f}")
    lines.append(f"低保证率(例如农业用水)的供水保证率 PB = {params['PB']:.3f}")
    k_desc = "K>0 表示保证率以年统计" if params['K'] > 0 else "K<=0 表示保证率以月统计"
    lines.append(f"{k_desc}. K = {params['K']}")
    lines.append(f"每年起始月份 (如日历年 M0=1,水文年一般为6) M0 = {params['MO']}")
    lines.append("")
    lines.append("计算结果:")
    lines.append("  年   月  月径流量  高保证率供水 低保证率供水 蒸发渗漏月损失量")
    for sp in result.get("storage_points", []):
        lines.append(f" {sp['年']:4d} {sp['月']:3d} {sp['径流']:7.0f} {sp['高供水']:9.0f} {sp['低供水']:9.0f} {sp['损失']:9.0f}")
    lines.append("")
    lines.append(f"兴利库容 VX= {result['VX']}")
    lines.append(f"限制库容 VY= {result['VY']}")
    lines.append(f"高保证率供水遭受破坏年数 AY= {result['AY']:3d}   高保证率供水年供水保证率 PAY= {result['PAY']:.3f}")
    lines.append(f"高保证率供水遭受破坏月数 AM= {result['AM']:3d}   低保证率供水年供水保证率 PAM= {result['PAM']:.3f}")
    lines.append(f"低保证率供水遭受破坏年数 BY= {result['BY']:3d}   高保证率供水月供水保证率 PBY= {result['PBY']:.3f}")
    lines.append(f"低保证率供水遭受破坏月数 BM= {result['BM']:3d}   低保证率供水月供水保证率 PBM= {result['PBM']:.3f}")
    lines.append("")
    lines.append("  年   月    水库蓄量")
    for sp in result.get("storage_points", []):
        if sp['蓄量'] < 0 or sp['蓄量'] >= result['VY']:
            lines.append(f" {sp['年']:4d} {sp['月']:3d} {sp['蓄量']:8.0f}")
    lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    params = parse(data)
    result = compute(params, cfg={"VX": params["VM"], "VY": params["VM"]})
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [("一", ["原始数据"]), ("二", ["结果见 JSON"])], result)
    else:
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
