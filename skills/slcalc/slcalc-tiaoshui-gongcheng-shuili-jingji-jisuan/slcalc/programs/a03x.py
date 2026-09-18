# -*- coding: utf-8 -*-
"""
A-3X 水文频率计算程序 —— 内核
===============================
复刻《水利程序集》A-3X 程序（作者：马明，新疆水利水电勘测设计院）。
在 A-3 基础上支持连续/不连续系列，功能：
  1. 系列排队 + 经验频率（连续 P=m/(n+1)；不连续含特大值修正公式）
  2. 矩法统计参数（连续用常规矩法；不连续用修正矩法，分母 N-1）
  3. P-Ⅲ 频率曲线适线优选（固定 Qa，单纯形加速法优选 Cv、Cs）
  4. 理论频率曲线设计值表（21 个设计频率含 P=100%）

数据文件顺序：
  连续系列：  站名,C,系列项数,连续系列第一个年份,系列值...
  不连续系列：站名,D,重现期N,连续系列项数n,特大值项数a,实测中抽出特大值项数L,
              (特大值年份,特大值值)×a 共 a 组..., 连续系列第一个年份,连续系列值...

验证基准：
  A-3X.INT（1096 项连续系列）：矩法 Qa=21.70 Cv=1.198 Cs=2.713；优选 CV=1.226 CS=3.114
  说明书例1（33 项连续）：Qa=107.00 CV=0.143 CS=0.709；优选 CV=0.157 CS=1.023
  说明书例2（不连续 N=90,n=26,a=4,L=2）：Qa=856.00 Cv=1.293；
    优选 CV=1.529 CS=3.142；采用 CV=1.5 CS=3.0 设计值表
"""
import re

from ..core.p3freq import (
    moment_params, moment_params_discontinuous,
    empirical_freqs_continuous, empirical_freqs_discontinuous,
    fitting_objective, p3_quantile,
)
from ..core.numext import nelder_mead
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-3X"
TITLE = "水文频率计算"

# 设计频率表（%），与原著输出一致（含 100%）
_STD_FREQS = [0.01, 0.1, 0.2, 0.5, 1, 2, 3, 5, 10, 20, 30, 40, 50,
              60, 70, 80, 90, 95, 97, 99, 100]


def parse(data):
    """
    解析输入。
    data: INT 文件路径 | 内容字符串 | dict。
    返回 dict：
      连续: station, series_type='C', values, start_year
      不连续: station, series_type='D', N, n, a, l,
              special=[(year,val),...], values(连续系列值含抽出特大值), start_year
    """
    if isinstance(data, dict):
        st = data.get("series_type", "C")
        base = {
            "station": data.get("station", ""),
            "series_type": st,
            "values": [float(v) for v in data.get("values", [])],
            "start_year": data.get("start_year"),
        }
        if st == "D":
            base.update({
                "N": int(data.get("N")),
                "a": int(data.get("a")),
                "l": int(data.get("l", 0)),
                "special": [(int(y), float(v)) for y, v in data.get("special", [])],
            })
        return base

    if isinstance(data, str):
        import os
        if os.path.isfile(data):
            content = open(data, encoding="gbk", errors="replace").read()
        else:
            content = data
    else:
        content = str(data)

    content = content.replace("\x1a", "").strip()

    # 连续系列：站名,C,项数,首年,数值...
    m = re.match(r"^(.*?),([CD]),\s*(\d+),\s*(\d+),(.*)$", content, re.S)
    if not m:
        raise ValueError("无法解析 A-3X 数据格式（需要：站名,C/D,...）")
    station, stype, n_str, year_str, rest = m.groups()
    station = station.strip()

    if stype == "C":
        nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", rest.replace(",", " "))]
        return {
            "station": station,
            "series_type": "C",
            "values": nums,
            "start_year": int(year_str),
        }

    # 不连续系列：站名,D,N,n,a,l, (年份,值)×a, 首年, 值...
    # 头部 5 个数字（N,n,a,l 后接换行，非逗号）
    head_match = re.match(
        r"^(.*?),D,\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)(.*)$", content, re.S)
    if not head_match:
        raise ValueError("无法解析 A-3X 不连续系列格式（需要：站名,D,N,n,a,l,...）")
    station, Ns, ns, as_, ls, rest2 = head_match.groups()
    N, n, a, l = int(Ns), int(ns), int(as_), int(ls)

    # 剩余部分：先 (年份,值) 对 × (a-L)（历史调查/考证特大值），再 首年, 值流
    # 注：数据文件只列出 a-L 个特大值对，从实测抽出的 L 个特大值留在实测系列中
    toks = re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", rest2.replace(",", " "))
    n_hist = a - l  # 历史特大值对数
    if len(toks) < 2 * n_hist + 1:
        raise ValueError(f"不连续系列数据不足：需要 {2*n_hist+1} 个头部数值，实际 {len(toks)}")
    special = []
    for i in range(n_hist):
        y = int(float(toks[2 * i]))
        v = float(toks[2 * i + 1])
        special.append((y, v))
    start_year = int(float(toks[2 * n_hist]))
    values = [float(x) for x in toks[2 * n_hist + 1:]]
    if len(values) != n:
        print(f"⚠ 声明实测项数 {n}，实际数值 {len(values)} 个，以实际为准")

    return {
        "station": station,
        "series_type": "D",
        "N": N, "a": a, "l": l,
        "special": special,
        "values": values,
        "start_year": start_year,
    }


def _discontinuous_split(params):
    """
    不连续系列分组：
      special_vals: 全部 a 个特大值 = (a-l) 个历史/调查特大值 + l 个实测抽出的特大值
      rest_vals:    实测中剩余的 n-l 个普通值（剔除从实测抽出的 l 个特大值）
    实测抽出判定：实测值中最大的 l 个（SL44 惯例：实测中的最大项为特大值）。
    """
    hist_vals = [float(v) for _, v in params["special"]]
    values = params["values"]
    a, l = params["a"], params["l"]

    # 从实测抽出的特大值 = 实测值中最大的 l 个
    extracted = sorted(values, reverse=True)[:l]
    rest = [float(v) for v in values if v not in extracted]
    special_vals = hist_vals + extracted
    return special_vals, rest


def compute(params, fit=True, freqs=None):
    """
    执行水文频率计算。
    返回结构化结果（含矩法参数、优选参数、设计值表）。
    """
    st = params.get("series_type", "C")
    freqs = freqs or _STD_FREQS

    if st == "D":
        special_vals, rest_vals = _discontinuous_split(params)
        if len(special_vals) + len(rest_vals) < 3:
            raise ValueError("数据点过少")
        mom = moment_params_discontinuous(special_vals, rest_vals,
                                          params["N"], params["n"] if "n" in params else len(params["values"]),
                                          params["a"], params["l"])
        emp = [(P, x) for _, x, P in empirical_freqs_discontinuous(
            special_vals, rest_vals, params["N"],
            params["n"] if "n" in params else len(params["values"]),
            params["a"], params["l"])]
    else:
        values = params["values"]
        if len(values) < 3:
            raise ValueError("至少需要 3 个数据")
        mom = moment_params(values)
        emp = [(P, x) for _, x, P in empirical_freqs_continuous(values)]

    Qa, Cv0, Cs0 = mom["Qa"], mom["Cv"], mom["Cs"]

    # 适线优选（固定 Qa，优化 Cv/Cs）
    fit_result = {"Cv": Cv0, "Cs": Cs0, "S": None, "converged": False}
    if fit:
        try:
            obj = lambda cvcs: fitting_objective(cvcs, Qa, emp)
            best, fval = nelder_mead(obj, [Cv0, Cs0], step=max(0.3, Cv0 * 0.5),
                                     tol=1e-10, max_iter=5000)
            if best[0] > 0 and best[1] > 0:
                fit_result = {"Cv": best[0], "Cs": best[1], "S": fval, "converged": True}
        except Exception as e:
            fit_result["error"] = str(e)

    Cv_f, Cs_f = fit_result["Cv"], fit_result["Cs"]

    # 设计值表
    design = []
    for P in freqs:
        try:
            if P == 100:
                # P=100% 时 Qp = 伽马分布下界 δ = Qa*(1-2Cv/Cs)
                delta = Qa * (1.0 - 2.0 * Cv_f / Cs_f)
                Qp, Kp = delta, delta / Qa
            else:
                Qp, Kp = p3_quantile(Qa, Cv_f, Cs_f, P)
            design.append({"P": P, "Kp": round(Kp, 3), "Qp": round(Qp, 3)})
        except Exception as e:
            design.append({"P": P, "error": str(e)})

    return {
        "程序": PROGRAM_ID,
        "站名": params.get("station", ""),
        "系列类型": st,
        "n": len(params["values"]) if st == "D" else len(params["values"]),
        "矩法参数": {"Qa": round(Qa, 3), "Cv": round(Cv0, 3), "Cs": round(Cs0, 3)},
        "优选参数": {"Cv": round(Cv_f, 3), "Cs": round(Cs_f, 3),
                   "R": round(Cs_f / Cv_f, 3) if Cv_f > 0 else None,
                   "S": round(fit_result["S"], 4) if fit_result["S"] else None,
                   "收敛": fit_result["converged"]},
        "设计值表": design,
    }


def render(params, result):
    """生成文本计算书（原著风格）。"""
    st = params.get("series_type", "C")
    lines = []

    lines.append(f"站名: {params.get('station', '')}")
    lines.append("")

    lines.append("一. 基  本  数  据")
    if st == "D":
        lines.append(f"重现期 N= {params['N']:>3d}   实测系列项数 n= {len(params['values']):>3d}")
        lines.append(f"特大值项数 a= {params['a']:>3d}   实测系列中抽出的特大值项数 L= {params['l']:>3d}")
        lines.append("")
        lines.append("  序     原序号        值                序     原序号         值")
        # 特大值区（a-L 个历史特大值，独立编号）
        sp_rows = [(i, y, float(v)) for i, (y, v) in enumerate(params["special"], start=1)]
        for i in range(0, len(sp_rows), 2):
            left = sp_rows[i]
            lstr = f"{left[0]:>3}   {left[1]:>6}   {left[2]:>10.2f}"
            if i + 1 < len(sp_rows):
                right = sp_rows[i + 1]
                rstr = f"{right[0]:>3}   {right[1]:>6}   {right[2]:>10.2f}"
            else:
                rstr = ""
            lines.append(f"{lstr:>30}  {rstr}")
        # 实测区（n 项，独立编号）
        start_year = params.get("start_year")
        for i in range(0, len(params["values"]), 2):
            left = f"{i+1:>3}   {start_year+i:>6}   {params['values'][i]:>10.2f}"
            right = f"{i+2:>3}   {start_year+i+1:>6}   {params['values'][i+1]:>10.2f}" if i + 1 < len(params["values"]) else ""
            lines.append(f"{left:>30}  {right}")
    else:
        lines.append(f"实测系列项数 n= {len(params['values'])}")
        lines.append("")
        lines.append("  序     原序号        值                序     原序号         值")
        start_year = params.get("start_year") or 0
        for i in range(0, len(params["values"]), 2):
            left = f"{i+1:>3}   {start_year+i:>6}   {params['values'][i]:>10.2f}"
            right = f"{i+2:>3}   {start_year+i+1:>6}   {params['values'][i+1]:>10.2f}" if i + 1 < len(params["values"]) else ""
            lines.append(f"{left:>30}  {right}")
    lines.append("")

    lines.append("二. 计  算  结  果")
    lines.append("")
    lines.append("经  验  频  率")
    # 经验频率表
    if st == "D":
        special_vals, rest_vals = _discontinuous_split(params)
        n_all = params["N"]
        emp_list = empirical_freqs_discontinuous(
            special_vals, rest_vals, params["N"],
            params.get("n", len(params["values"])), params["a"], params["l"])
    else:
        emp_list = empirical_freqs_continuous(params["values"])
    lines.append("  序  原序号       值      频率(%)       序  原序号       值      频率(%)")
    pairs = []
    if st == "D":
        special_map = {v: y for y, v in params["special"]}
        for i, (_, x, P) in enumerate(emp_list, start=1):
            if x in special_map:
                pairs.append((i, special_map[x], x, P))
            else:
                # 从 values 找年份（首个匹配）
                for j, v in enumerate(params["values"]):
                    if abs(v - x) < 1e-9:
                        pairs.append((i, params.get("start_year", 0) + j, x, P))
                        break
                else:
                    pairs.append((i, "", x, P))
    else:
        start_year = params.get("start_year") or 0
        sorted_vals = sorted(params["values"], reverse=True)
        for i, (_, x, P) in enumerate(emp_list, start=1):
            idx = params["values"].index(x) if x in params["values"] else i - 1
            pairs.append((i, start_year + idx, x, P))
    for i in range(0, len(pairs), 2):
        lp = pairs[i]
        lstr = f"{lp[0]:>3}  {lp[1]:>6}   {lp[2]:>10.2f}   {lp[3]:>7.2f}"
        if i + 1 < len(pairs):
            rp = pairs[i + 1]
            rstr = f"{rp[0]:>3}  {rp[1]:>6}   {rp[2]:>10.2f}   {rp[3]:>7.2f}"
        else:
            rstr = ""
        lines.append(f"{lstr:>36}  {rstr}")
    lines.append("")

    mom = result["矩法参数"]
    fit = result["优选参数"]
    lines.append("统  计  参  数  值")
    lines.append(f"均值 Qa= {mom['Qa']:>8.3f}   CV= {mom['Cv']:>6.3f}           CS= {mom['Cs']:>6.3f}")
    lines.append("")
    lines.append("参  数  优  选  值")
    R = fit.get("R")
    lines.append(f"                    CV= {fit['Cv']:>6.3f}          CS= {fit['Cs']:>6.3f}          R= {R if R is not None else ''}")
    lines.append("")
    lines.append("理  论  频  率  曲  线  设  计  值  表")
    R_f = fit.get("R")
    lines.append(f"Qa= {mom['Qa']:>8.3f}       CV= {fit['Cv']:>6.3f}          CS= {fit['Cs']:>6.3f}          R= {R_f if R_f is not None else ''}")
    lines.append("频率P(%)      KP        设计值QP       频率P(%)      KP        设计值QP")
    design = result["设计值表"]
    half = (len(design) + 1) // 2
    for i in range(half):
        left = design[i]
        lstr = f"{left['P']:>8.3f}   {left['Kp']:>8.3f}   {left['Qp']:>10.3f}"
        if i + half < len(design):
            right = design[i + half]
            rstr = f"{right['P']:>8.3f}   {right['Kp']:>8.3f}   {right['Qp']:>10.3f}"
        else:
            rstr = ""
        lines.append(f"{lstr:>35}  {rstr}")
    lines.append("")

    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。"""
    params = parse(data)
    result = compute(params)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [("一", [f"站名: {params.get('station','')}"]),
                                                   ("二", ["结果见 JSON"])], result)
    else:
        text = render(params, result)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(txt)
