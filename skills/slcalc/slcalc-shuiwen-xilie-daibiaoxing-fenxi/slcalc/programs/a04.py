# -*- coding: utf-8 -*-
"""
A-4 水文系列代表性分析程序 —— 内核
=====================================
复刻《水利程序集》A-4 程序（作者：孙建峰，水电部天津勘测设计院）。

功能（5 项分析，可独立开关）：
  ㈠ 频率计算：按年排队（降序），P = i/(N+1)×100（经验频率）
  ㈡ 统计参数计算：取系列前 m 个（自最近年份起），m=5..N，
     均值 X̄、Cv（无偏 n-1 分母）、Cs（n-3 修正式）、Cs/Cv
  ㈢ 累积平均曲线：k_m = (前 m 个之和/m) / X̄全（模比系数累积均值）
  ㈣ 差积曲线：k_i = Σ(Kj−1)（模比系数减 1 累积）
  ㈤ 滑动平均曲线：k_i = (i 起 HT 年滑动均值) / X̄全

数据文件顺序（原著）：
  系列长,最近年份
  观测值（按时间降序：最近→最早）,...
  做频率计算输入 1 否则 0, 做参数计算输入 1 否则 0,
  做累积计算输入 1 否则 0, 做差积计算输入 1 否则 0,
  做滑动计算输入 1, 滑动年数 HT
  （若控制位缺省，默认全部计算、滑动年数 5——算例 A-4.INT 即此形式）

验证基准（A-4.INT：某站 1951~1980 年 30 年降水）：
  频率表 30 行：排队值 985.40→230.10，P=i/31
  统计参数表 26 行（m=5..30）全部命中
  累积平均 / 差积 / 滑动平均（HT=5）表全部命中
  回归状态：五项输出全基准验证通过 ✓（2026-09-02）
"""
import math
import os
import re

from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-4"
TITLE = "水文系列代表性分析计算书"

# 各分析项的默认开关
_DEFAULTS = {
    "freq": True, "params": True, "cumavg": True,
    "resid": True, "moving": True, "ht": 5,
}


# ------------------------------------------------------------
# 基础计算
# ------------------------------------------------------------

def _fmt(x, nd=2):
    """数字格式化：去掉多余的 0（原著风格 685.60/0.32/1.1）"""
    return f"{x:.{nd}f}"


def freq_table(values, years):
    """
    (一) 频率计算：按年份行输出，排队值=全体降序序列，P=i/(N+1)。
    原著格式：
      年份    值X(该年)   排队值X(%)   频率P(%)
    排队值列 = 降序排列的第 i 个值（与年份无关，i 为行序）。
    返回 [(year, x, rank_i_val, P), ...] 按年份降序。
    """
    n = len(values)
    sorted_desc = sorted(values, reverse=True)
    rows = []
    for i in range(n):
        P = (i + 1) / (n + 1.0) * 100.0
        rows.append((years[i], values[i], sorted_desc[i], P))
    return rows


def stat_params_table(values):
    """
    (二) 统计参数计算：m=5..N。
    均值   = 前 m 个（自最近年份起）的算术平均
    Cv     = sqrt(Σ(x−X̄)²/(m−1)) / X̄
    Cs     = Σ(Ki−1)³ / [(m−3)·Cv³]，Ki = x/X̄
    Cs/Cv  = 比值
    返回 [(m, mean, cv, cs, cs_cv), ...]
    """
    n = len(values)
    rows = []
    for m in range(5, n + 1):
        seg = values[:m]
        mean = sum(seg) / m
        if m > 1:
            var = sum((x - mean) ** 2 for x in seg) / (m - 1.0)
            cv = math.sqrt(var) / mean if mean != 0 else 0.0
        else:
            cv = 0.0
        if m > 3 and cv > 0:
            cs = sum((x / mean - 1.0) ** 3 for x in seg) / ((m - 3.0) * cv ** 3)
        else:
            cs = 0.0
        rows.append((m, mean, cv, cs, cs / cv if cv > 0 else 0.0))
    return rows


def cumavg_series(values):
    """
    (三) 累积平均曲线：k_m = (前 m 个均值)/全体均值。
    原著输出每行含两个年份（左列 1980..1966，右列 1965..1951）。
    返回 [(year, k), ...] 按年份降序。
    """
    n = len(values)
    mean_all = sum(values) / n
    k = []
    s = 0.0
    for i, x in enumerate(values, start=1):
        s += x
        k.append((s / i) / mean_all if mean_all != 0 else 0.0)
    return k


def resid_series(values):
    """
    (四) 差积曲线：k_i = Σ_{j≤i}(Kj − 1)，Kj = x/X̄全。
    返回 [(i, k_i)]（0-based 索引）
    """
    mean_all = sum(values) / len(values)
    out = []
    acc = 0.0
    for x in values:
        acc += (x / mean_all - 1.0) if mean_all != 0 else 0.0
        out.append(acc)
    return out


def moving_avg_series(values, ht):
    """
    (五) 滑动平均曲线（HT 年滑动）：
    k_i = (自 i 起 ht 年滑动均值)/全体均值；不足 ht 年（尾部）用剩余年数均值。
    返回 [(i, k)]。
    """
    n = len(values)
    mean_all = sum(values) / n
    out = []
    for i in range(n):
        seg = values[i:i + ht]
        k = (sum(seg) / len(seg)) / mean_all if mean_all != 0 else 0.0
        out.append(k)
    return out


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """解析输入。data: INT 路径 | 字符串 | dict。"""
    if isinstance(data, dict):
        return data
    if isinstance(data, str) and os.path.isfile(data):
        with open(data, encoding="gbk", errors="replace") as f:
            content = f.read()
    else:
        content = str(data)
    content = content.replace("\x1a", "").strip()

    nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", content)]
    if len(nums) < 4:
        raise ValueError("A-4 数据不足：需要 系列长,最近年份,观测值...,控制位")
    N = int(nums[0])
    last_year = int(nums[1])
    # 观测值：紧接着 N 个（可能是整数年份不参与——按 N 个数值取）
    values = nums[2:2 + N]
    # 剩余为控制位/滑动年数；若剩余 ≥1 则从剩余取
    rest = nums[2 + N:]
    opts = dict(_DEFAULTS)
    if rest:
        # 形式 A: 1,1,1,1,1,5 → 5 个开关 + HT
        if len(rest) >= 5 and all(v in (0.0, 1.0) for v in rest[:5]):
            opts["freq"] = bool(rest[0])
            opts["params"] = bool(rest[1])
            opts["cumavg"] = bool(rest[2])
            opts["resid"] = bool(rest[3])
            opts["moving"] = bool(rest[4])
            if len(rest) >= 6 and rest[5] >= 1:
                opts["ht"] = int(rest[5])
        else:
            # 形式 B: 仅一个滑动年数（算例 A-4.INT: 5）
            if rest[0] >= 1:
                opts["ht"] = int(rest[0])
    years = list(range(last_year, last_year - N, -1))
    if len(values) < N:
        # 容错：数值不足时以实际为准
        print(f"⚠ A-4 声明 {N} 项，实际数值 {len(values)} 个，以实际为准")
        N = len(values)
        years = list(range(last_year, last_year - N, -1))
    return {
        "N": N, "last_year": last_year, "values": values, "years": years,
        "opts": opts,
    }


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params):
    """执行五项水文系列代表性分析。"""
    values = params["values"]
    years = params["years"]
    N = params["N"]
    opts = params["opts"]
    mean_all = sum(values) / N if N else 0.0

    res = {
        "程序": PROGRAM_ID,
        "输入": {
            "系列长N": N, "最近年份": params["last_year"],
            "观测值": [round(v, 4) for v in values],
            "分析项": {k: opts[k] for k in
                      ("freq", "params", "cumavg", "resid", "moving")},
            "滑动年数HT": opts["ht"],
        },
        "全系列均值": round(mean_all, 4),
    }

    if opts["freq"]:
        fr = freq_table(values, years)
        res["频率计算"] = {
            "方法": "P = i/(N+1)",
            "行": [(yr, round(x, 2), round(rv, 2), round(P, 2))
                   for yr, x, rv, P in fr],
        }

    if opts["params"]:
        sp = stat_params_table(values)
        res["统计参数"] = {
            "方法": "均值=前m个; Cv无偏(n-1); Cs=Σ(K−1)³/[(m−3)Cv³]",
            "行": [(m, round(mean, 2), round(cv, 2), round(cs, 2), round(cs / cv, 1) if cv else 0.0)
                   for m, mean, cv, cs, _ in sp],
        }

    if opts["cumavg"]:
        ca = cumavg_series(values)
        res["累积平均曲线"] = {
            "方法": "k_m = (前m个均值)/全系列均值",
            "行": [(yr, round(k, 2)) for yr, k in zip(years, ca)],
        }

    if opts["resid"]:
        rd = resid_series(values)
        res["差积曲线"] = {
            "方法": "k_i = Σ(Kj−1), Kj=x/X̄",
            "行": [(yr, round(k, 2)) for yr, k in zip(years, rd)],
        }

    if opts["moving"]:
        mv = moving_avg_series(values, opts["ht"])
        res["滑动平均曲线"] = {
            "方法": f"k_i = {opts['ht']}年滑动均值/全系列均值",
            "行": [(yr, round(k, 2)) for yr, k in zip(years, mv)],
        }

    return res


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def render(params, result):
    """生成文本计算书（原著 .OUT 风格）。"""
    L = []
    inp = result["输入"]
    L.append(f"系列长 N= {inp['系列长N']}   最近年份= {inp['最近年份']}")
    L.append("")
    vals = inp["观测值"]
    # 回显观测值（10 个一行）
    L.append("观测值：")
    for i in range(0, len(vals), 10):
        L.append("  " + ", ".join(f"{v:>9.2f}" for v in vals[i:i + 10]))
    L.append("")
    L.append(f"滑动年数 HT= {inp['滑动年数HT']}")
    L.append("")

    if "频率计算" in result:
        L.append("(一）频率计算结果：")
        L.append(" 年份       值X     排队值X(%)  频率P(%)")
        for yr, x, rv, P in result["频率计算"]["行"]:
            L.append(f"{yr}    {x:>10.2f}   {rv:>10.2f}   {P:>6.2f}")

    if "统计参数" in result:
        L.append("")
        L.append("(二) 统计参数的计算结果：")
        L.append(" 序号   均值       Cv      Cs    CS/CV")
        for m, mean, cv, cs, r in result["统计参数"]["行"]:
            L.append(f"{m:>4}  {mean:>9.2f}  {cv:>6.2f}  {cs:>6.2f}  {r:>5.1f}")

    # 累积平均 / 差积 / 滑动平均：两列（左 1980..1966，右 1965..1951）
    half = (result["输入"]["系列长N"] + 1) // 2
    for key, title in (("累积平均曲线", "累 积 平 均 曲 线"),
                       ("差积曲线", "差 积 平 均 曲 线"),
                       ("滑动平均曲线", "滑 动 平 均 曲 线")):
        if key in result:
            L.append("")
            L.append(f"        {title}")
            L.append("年份       k        年份      k")
            rows = result[key]["行"]
            for i in range(half):
                y1, k1 = rows[i]
                if i + half < len(rows):
                    y2, k2 = rows[i + half]
                    L.append(f" {y1}    {k1:>7.2f}      {y2}   {k2:>7.2f}")
                else:
                    L.append(f" {y1}    {k1:>7.2f}")

    return render_text(PROGRAM_ID, TITLE, [("", L)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 路径 | dict"""
    params = parse(data)
    result = compute(params)
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
