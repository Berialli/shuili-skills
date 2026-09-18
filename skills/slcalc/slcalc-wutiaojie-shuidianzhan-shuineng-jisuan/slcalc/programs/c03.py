# -*- coding: utf-8 -*-
"""
C-3 两种供水保证的水库径流调节计算程序 —— 内核（v2，规则闭合版）
==================================================================
复刻《水利程序集》C-3 程序（作者：王建生，水电部天津勘测设计院）。

【2026-09-05 规则破译】基于权威 OUT 101/101 逐点复现验证（5 组 VX/VY 全部
打印点），蓄量演算真实规则为：

    每月水量平衡（时历法）：
        起调：1958.6 月初 V0 = 0
        若 V0 > VY（限制库容）：低保证供水 B 可供，但以"不使月末蓄量跌破
            VY"为限：  B_eff = max(0, min(B计划, V0+R-A-E-VY))
        若 V0 <= VY          ：B 削减为 0（免供，记低保证供水破坏）
        月末蓄量 V = V0 + R - A - B_eff - E
        若 V > VX（兴利库容）：弃水，V = VX
        若 V < 0            ：高保证供水 A 破坏（负蓄量=破坏深度），
                              下月从 V=0 重新起蓄
    打印规则：仅当 V <= VY（含负值）时打印"年 月 蓄量"

    【破坏统计口径 v7（2026-09-05 收官修正，表2 22/22 + 表3 5/5 全绿）】
        低保证破坏月集 F = {V<VY} ∪ {V==VY 且当月 B 真被压减(0<Be<Bp)}
        ── 月末恰好 V==VY 的"顶格贴线"月本身不算 B 破坏；
           仅当由 B 压减造成(VY 拦截灌溉供水, Be<Bp)才计；若 Bp=0 自然耗至
           VY(如 1961.12)不计。此修正解决旧口径 BM=#{V<=VY} 在 6 组
           (6750/1675 等)多计 1 月的问题。BM=|F|, BY=含F月的年数。

    1961.10 = 1522 的机制（原"2716/129 差值"悬案）：
        1961.9 末 V0 = 1966（> VY=1522），10 月 R=48、A=182、E=57、
        B计划=269。若 B 全额供：V = 1966+48-182-269-57 = 1506 < VY，
        违反"不跌破限制库容"约束，故 B 被削减至
        avail = 1966+48-182-57-1522 = 253 → V = 1522 恰好顶格 VY。
        （此前"2716"系把水文年数据当公历月错位所致，真实 1961.10 R=48 非 41。）

作者：王建生，水电部天津勘测设计院。
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
    p["A"] = nums[idx:idx + 12]; idx += 12   # 从 MO=6 起排
    p["E"] = nums[idx:idx + 12]; idx += 12
    years = []
    while idx < len(nums):
        year = int(nums[idx]); idx += 1
        if year <= 0 and idx >= len(nums):
            break
        R = []; B = []
        for m in range(12):
            R.append(nums[idx]); idx += 1
            B.append(nums[idx]); idx += 1
        years.append({"year": year, "R": R, "B": B})
    p["years"] = years
    return p


def _month_of(mo_idx, year, MO):
    """INT 列序 mo_idx(0..11) -> (公历年, 公历月)，MO 为起调月(1..12)"""
    cm = mo_idx + MO
    while cm > 12:
        cm -= 12
        year += 1
    return year, cm


def compute(params, cfg):
    N = params["N"]
    A = params["A"]
    E = params["E"]
    years = params["years"]
    MO = params["MO"]
    VX = cfg["VX"]
    VY = cfg["VY"]

    V = 0.0  # 起调：1958.6 月初空库
    months = []
    for ydata in years:
        year = ydata["year"]
        for m in range(12):
            cy, cm = _month_of(m, year, MO)
            R = ydata["R"][m]
            B_plan = ydata["B"][m]
            a = A[m]
            e = E[m]
            V0 = V
            # --- 低保证供水 B 削减规则（破译） ---
            if V0 > VY:
                avail = V0 + R - a - e - VY
                B_eff = max(0.0, min(B_plan, avail))
            else:
                B_eff = 0.0
            # --- 月末蓄量 ---
            V_new = V0 + R - a - B_eff - e
            if V_new > VX:
                V_new = float(VX)          # 满兴利库容弃水
            b_cut = B_eff < B_plan
            months.append({
                "年": cy, "月": cm,
                "月初蓄量_raw": V0,
                "蓄量_raw": V_new,
                "月初蓄量": round(V0),
                "蓄量": round(V_new),
                "径流": round(R),
                "高供水": round(a),
                "低供水计划": round(B_plan),
                "低供水实供": round(B_eff),
                "损失": round(e),
                "B削减": b_cut,
                "A破坏": V_new < 0,
            })
            if V_new < 0:
                V = 0.0
            else:
                V = V_new

    # ---- 破坏统计（2026-09-05 v7 口径落地，表2 22/22 + 表3 5/5 权威全绿） ----
    # AY = 年内存在月末蓄量 V<=0 的月 → A 破坏年数
    # AM = #{月末 V<0} ∪ #{月初 V0==0, 首月不计}（V==0 恰好只计 AY 不计 AM）
    # 低保证破坏月集 F = {月末 V<VY} ∪ {月末 V==VY 且当月 B 真被压减(0<Be<Bp)}
    #   BY = F 中年数;  BM = |F|
    #   注：V==VY 顶格贴线中，由 B 压减造成(VY 拦水, Be<Bp)的才算 B 破坏月；
    #       Bp=0 自然耗至 VY(如 1961.12)不计（旧口径 V<=VY 会多计 1 月）。
    # 统计一律用原始浮点值判定（与破译探针一致），round 仅用于显示。
    ay_years, by_years = set(), set()
    am_months = 0
    for i, r in enumerate(months):
        vr = r["蓄量_raw"]
        v0r = r["月初蓄量_raw"]
        # 月末蓄量<=0 → A 破坏年（含恰好为 0）
        if vr <= 0:
            ay_years.add(r["年"])
        # AM：V<0（供水穿底）或 V0==0（空库重启），首月 V0=0 不计
        if vr < 0 or (i >= 1 and v0r == 0):
            am_months += 1
        # B 破坏月判定（v7）：V<VY 一律计；V==VY 仅当 B 实被压减(B_eff<B_plan)
        if vr < VY or (abs(vr - VY) < 1e-9 and r["低供水实供"] < r["低供水计划"]):
            by_years.add(r["年"])
    AY, BY = len(ay_years), len(by_years)
    # BM = B 破坏月数（与 by_years 同一 F 集）
    bm_months = sum(1 for r in months
                    if r["蓄量_raw"] < VY
                    or (abs(r["蓄量_raw"] - VY) < 1e-9
                        and r["低供水实供"] < r["低供水计划"]))

    return {
        "AY": AY, "AM": am_months, "BY": BY, "BM": bm_months,
        "months": months,
        "VX": VX, "VY": VY,
        "N": N,
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
    for r in result.get("months", []):
        lines.append(f" {r['年']:4d} {r['月']:3d} {r['径流']:7.0f} {r['高供水']:9.0f} {r['低供水实供']:9.0f} {r['损失']:9.0f}")
    lines.append("")
    lines.append(f"兴利库容 VX= {result['VX']:.0f}")
    lines.append(f"限制库容 VY= {result['VY']:.0f}")
    lines.append(f"高保证率供水遭受破坏年数 AY= {result['AY']:3d}")
    lines.append(f"高保证率供水遭受破坏月数 AM= {result['AM']:3d}")
    lines.append(f"低保证率供水遭受破坏年数 BY= {result['BY']:3d}")
    lines.append(f"低保证率供水遭受破坏月数 BM= {result['BM']:3d}")
    lines.append("")
    lines.append("  年   月    水库蓄量")
    for r in result.get("months", []):
        if r["蓄量_raw"] <= result["VY"]:
            lines.append(f" {r['年']:4d} {r['月']:3d} {r['蓄量']:8.0f}")
    lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text", cfg=None):
    params = parse(data)
    result = compute(params, cfg or {"VX": params["VM"], "VY": params["VM"]})
    text = render(params, result, None)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    cfg = None
    if len(sys.argv) > 2:
        cfg = {"VX": float(sys.argv[2]), "VY": float(sys.argv[3])}
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None, cfg=cfg)
    print(txt)
