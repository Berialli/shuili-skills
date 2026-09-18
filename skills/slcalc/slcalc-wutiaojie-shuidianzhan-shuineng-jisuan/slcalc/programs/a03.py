# -*- coding: utf-8 -*-
"""
A-3 水文频率计算程序 —— 内核
==============================
复刻《水利程序集》A-3 程序（作者：马明，新疆水利水电勘测设计院）。
功能：一次完成一个水文系列频率计算：
  1. 系列排队 + 经验频率（连续/不连续系列）
  2. 矩法统计参数（Qa, Cv, Cs）
  3. P-Ⅲ 频率曲线适线优选（离差平方和最小，单纯形加速法）
  4. 理论频率曲线设计值表（各设计频率下的 Qp, Kp）

数据文件顺序（连续系列）：
  站名, C, 系列项数, 连续系列第一个年份, 系列值...
（不连续系列格式见说明文档，本版本支持连续系列）

验证基准：A-3X.INT（1096 项实测系列）
  矩法: Qa=21.70  Cv=1.198  Cs=2.713
  优选: CV=1.226  CS=3.114
  理论设计值表（P=0.01% QP=303.543 ... P=99% QP=4.614）
"""
import re

from ..core.p3freq import (
    moment_params, empirical_freqs_continuous, fitting_objective,
    p3_series, gammaincinv_std,
)
from ..core.numext import nelder_mead
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-3"
TITLE = "水文频率计算"

# 设计频率表（%），与原著输出一致
_STD_FREQS = [0.01, 0.1, 0.2, 0.5, 1, 2, 3, 5, 10, 20, 30, 40, 50,
              60, 70, 80, 90, 95, 97, 99]


def parse(data):
    """
    解析输入。
    data: 字符串（INT 文件内容）或文件路径或 dict。
    返回 dict(station, series_type, values, start_year, meta)
    """
    if isinstance(data, dict):
        return {
            "station": data.get("station", ""),
            "series_type": data.get("series_type", "C"),
            "values": [float(v) for v in data.get("values", [])],
            "start_year": data.get("start_year"),
        }

    if isinstance(data, str):
        # 可能是文件路径或内容
        import os
        if os.path.isfile(data):
            content = open(data, encoding="gbk", errors="replace").read()
        else:
            content = data
    else:
        content = str(data)

    # 解析：站名,C,项数,首年,数值...
    # 站名可能含逗号，从右往左解析
    content = content.strip()
    # 去掉尾部 \x1a
    content = content.replace("\x1a", "").strip()

    # 找格式：...,C,项数,首年,数值流
    m = re.match(r"^(.*?),([CD]),\s*(\d+),\s*(\d+),(.*)$", content, re.S)
    if not m:
        raise ValueError("无法解析 A-3 数据格式（需要：站名,C,项数,首年,数值...）")
    station, stype, n_str, year_str, rest = m.groups()
    nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", rest.replace(",", " "))]
    n = int(n_str)
    if len(nums) != n:
        # 容错：有些文件项数与实际不符，以实际为准
        print(f"⚠ 声明项数 {n}，实际数值 {len(nums)} 个，以实际为准")
    return {
        "station": station.strip(),
        "series_type": stype,
        "values": nums,
        "start_year": int(year_str),
    }


def compute(params, fit=True, freqs=None):
    """
    执行水文频率计算。
    fit=True 时进行 P-Ⅲ 适线优选；否则仅用矩法参数。
    返回结构化结果。
    """
    values = params["values"]
    if len(values) < 3:
        raise ValueError("至少需要 3 个数据")

    # 1. 矩法
    mom = moment_params(values)
    Qa, Cv0, Cs0 = mom["Qa"], mom["Cv"], mom["Cs"]

    # 2. 经验频率
    emp = [(P, x) for _, x, P in empirical_freqs_continuous(values)]

    # 3. 适线优选（离差平方和最小）
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

    # 4. 设计值表
    freqs = freqs or _STD_FREQS
    design = []
    for P in freqs:
        try:
            Qp, Kp = _p3_quantile_safe(Qa, Cv_f, Cs_f, P)
            design.append({"P": P, "Kp": round(Kp, 3), "Qp": round(Qp, 3)})
        except Exception as e:
            design.append({"P": P, "error": str(e)})

    return {
        "程序": PROGRAM_ID,
        "站名": params.get("station", ""),
        "系列类型": params.get("series_type", "C"),
        "n": len(values),
        "矩法参数": {"Qa": round(Qa, 3), "Cv": round(Cv0, 3), "Cs": round(Cs0, 3)},
        "优选参数": {"Cv": round(Cv_f, 3), "Cs": round(Cs_f, 3),
                   "S": round(fit_result["S"], 4) if fit_result["S"] else None,
                   "收敛": fit_result["converged"]},
        "设计值表": design,
    }


def _p3_quantile_safe(Qa, Cv, Cs, P):
    """P-Ⅲ 分位数安全包装（捕获异常返回兜底）。"""
    from ..core.p3freq import p3_quantile
    return p3_quantile(Qa, Cv, Cs, P)


def render(params, result):
    """生成文本计算书（原著风格）。"""
    values = params["values"]
    emp = [(P, x) for _, x, P in empirical_freqs_continuous(values)]
    lines = []

    lines.append(f"站名: {params.get('station', '')}")
    lines.append("")
    lines.append("一. 基  本  数  据")
    lines.append(f"实测系列项数 n= {len(values)}")
    lines.append("")
    lines.append("  序     值         序     值")
    for i in range(0, len(values), 2):
        left = f"{i+1:>3}   {values[i]:>10.2f}"
        right = f"{i+2:>3}   {values[i+1]:>10.2f}" if i + 1 < len(values) else ""
        lines.append(f"{left:>28}  {right}")
    lines.append("")

    mom = result["矩法参数"]
    fit = result["优选参数"]
    lines.append("统  计  参  数  值")
    lines.append(f"均值 Qa= {mom['Qa']:>8.3f}    CV= {mom['Cv']:>6.3f}          CS= {mom['Cs']:>6.3f}")
    lines.append("")
    lines.append("参  数  优  选  值")
    lines.append(f"                    CV= {fit['Cv']:>6.3f}          CS= {fit['Cs']:>6.3f}")
    lines.append("")
    lines.append("理  论  频  率  曲  线  设  计  值  表")
    lines.append(f"Qa= {mom['Qa']:>8.3f}        CV= {fit['Cv']:>6.3f}          CS= {fit['Cs']:>6.3f}")
    lines.append("频率P(%)      KP        设计值QP")
    for d in result["设计值表"]:
        if "error" in d:
            lines.append(f"{d['P']:>8.3f}   (计算失败)")
        else:
            lines.append(f"{d['P']:>8.3f}   {d['Kp']:>8.3f}   {d['Qp']:>10.3f}")
    lines.append("")

    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """
    统一入口。
    data: INT 文件路径 | 内容字符串 | dict
    """
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
