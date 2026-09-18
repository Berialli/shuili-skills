# -*- coding: utf-8 -*-
"""
A-2 三参数幂函数曲线拟合程序 —— 内核
====================================
复刻《水利程序集》A-2 程序（作者：陈沂，水电部天津勘测设计院）。
以三参数函数 Y = a·X^b + c 为数学模型，采用最小二乘法对一组实验
（或实测）数据进行拟合（如水位~流量、水位~库容曲线）。

方法：
  将 Y = a·X^b + c 改写为 Y = a·x' + c 的直线方程，其中 x' = X^b。
  对 X>0 取 x' = X^b，X=0 时 x' = 0。
  在 [B1, B2] 区间内按步长 E 扫描 b（优选），对每个 b 用最小二乘
  求 a、c，取相关系数平方 r² 最大的 b 为优选值；
  第二段对优选 b 四舍五入到 0.01 步长精度，再求 a、c 为采用成果。

数据文件顺序（INT）：
  数组对数 n,
  成对数组数据 Xi（逗号分隔）,
  成对数组数据 Yi（逗号分隔）,
  b 值优选较小值 B1, b 值优选较大值 B2, b 值优选精度 E

验证基准：A-2.INT（水位~流量 6 点据，B1=0, B2=1, E=0.02）
  说明书输出：
    优选段 b=0.43614, a=0.15520, c=51.99003, r²=0.99988, Dcp=0.01492, Dmax=0.03455
    采用段 b=0.44000, a=0.15041, c=52.00016, r²=0.99990, Dcp=0.01347, Dmax=0.02950
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import render_text, render_markdown, write_out, write_json

PROGRAM_ID = "A-2"
TITLE = "三参数幂函数曲线拟合计算书"


def parse(data):
    """
    解析输入。
    data: 数值流（list[float]）或 dict。
    返回 dict(xs=[...], ys=[...], b1=float, b2=float, e=float)
    """
    if isinstance(data, dict):
        xs = [float(v) for v in data.get("xs", [])]
        ys = [float(v) for v in data.get("ys", [])]
        if not xs or len(xs) != len(ys):
            raise ValueError("xs 与 ys 必须等长且非空")
        return {
            "xs": xs, "ys": ys,
            "b1": float(data.get("b1", 0.0)),
            "b2": float(data.get("b2", 1.0)),
            "e": float(data.get("e", 0.01)),
        }

    # 数值流：n, x1..xn, y1..yn, B1, B2, E
    nums = list(data)
    if len(nums) < 4:
        raise ValueError("数据不足：需要 n, Xi 系列, Yi 系列, B1, B2, E")
    n = int(nums[0])
    need = 1 + 2 * n + 3
    if len(nums) < need:
        raise ValueError(f"数据不足：声明 n={n}，实际仅有 {len(nums)-1} 个数值（应 {need-1}）")
    xs = nums[1:1 + n]
    ys = nums[1 + n:1 + 2 * n]
    b1 = nums[1 + 2 * n]
    b2 = nums[2 + 2 * n]
    e = nums[3 + 2 * n]
    if e <= 0:
        raise ValueError("优选精度 E 必须 > 0")
    return {"xs": xs, "ys": ys, "b1": b1, "b2": b2, "e": e}


def _fit(xs, ys, b):
    """
    对给定 b 做线性回归 y = a·X^b + c。
    返回 (a, c, r2, y_pred, dcp, dmax)
    """
    n = len(xs)
    X = [x ** b if x > 0 else 0.0 for x in xs]
    sx = sum(X)
    sy = sum(ys)
    sxx = sum(x * x for x in X)
    sxy = sum(x * y for x, y in zip(X, ys))
    denom = n * sxx - sx * sx
    if denom == 0:
        raise ValueError("回归方程奇异（数据退化），请检查 X 值")
    a = (n * sxy - sx * sy) / denom
    c = (sy - a * sx) / n
    yp = [a * x + c for x in X]
    ss_res = sum((y - p) ** 2 for y, p in zip(ys, yp))
    mean_y = sy / n
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    dcp = sum(abs(y - p) for y, p in zip(ys, yp)) / n
    dmax = max(abs(y - p) for y, p in zip(ys, yp))
    return a, c, r2, yp, dcp, dmax


def _scan_best(xs, ys, b1, b2, e):
    """
    在 [b1, b2] 区间用黄金分割法求 r² 最大的 b（连续寻优，非格点扫描）。
    e 为收敛精度（说明书中的"优选精度 E"）。
    返回 (b, fit)，fit=(a, c, r2, yp, dcp, dmax)。
    """
    phi = (math.sqrt(5.0) - 1.0) / 2.0  # 黄金分割比 ≈0.618

    def obj(b):
        try:
            return _fit(xs, ys, b)[2]
        except ValueError:
            return -1.0

    lo, hi = float(b1), float(b2)
    # 区间退化保护
    if hi - lo < 1e-12:
        return lo, _fit(xs, ys, lo)

    x1 = hi - phi * (hi - lo)
    x2 = lo + phi * (hi - lo)
    f1 = obj(x1)
    f2 = obj(x2)
    # 保护：E 过小或区间过大时限制迭代次数
    max_iter = 10000
    it = 0
    while abs(hi - lo) > e and it < max_iter:
        if f1 > f2:
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi - phi * (hi - lo)
            f1 = obj(x1)
        else:
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo + phi * (hi - lo)
            f2 = obj(x2)
        it += 1

    b_best = (lo + hi) / 2.0
    return b_best, _fit(xs, ys, b_best)


def compute(params):
    """执行两段式拟合（优选 b + 确定 b 重算）。"""
    xs, ys = params["xs"], params["ys"]
    b1, b2, e = params["b1"], params["b2"], params["e"]

    # 第一段：优选 b
    best_b, (a1, c1, r2_1, yp1, dcp1, dmax1) = _scan_best(xs, ys, b1, b2, e)

    # 第二段：b 四舍五入到 0.01 精度，作为采用成果
    b_use = round(best_b * 100.0) / 100.0
    a2, c2, r2_2, yp2, dcp2, dmax2 = _fit(xs, ys, b_use)

    return {
        "程序": PROGRAM_ID,
        "n": len(xs),
        "b优选段": {
            "b": round(best_b, 5), "a": round(a1, 5), "c": round(c1, 5),
            "r2": round(r2_1, 5), "Dcp": round(dcp1, 5), "Dmax": round(dmax1, 5),
            "Y回代": [round(v, 3) for v in yp1],
        },
        "b确定段": {
            "b": round(b_use, 5), "a": round(a2, 5), "c": round(c2, 5),
            "r2": round(r2_2, 5), "Dcp": round(dcp2, 5), "Dmax": round(dmax2, 5),
            "Y回代": [round(v, 3) for v in yp2],
        },
        "方程": f"Y = {a2:.5f}·X^{b_use:.5f} + {c2:.5f}",
    }


def render(params, result):
    """生成计算书文本（复刻原版三段式）。"""
    xs, ys = params["xs"], params["ys"]
    b1, b2, e = params["b1"], params["b2"], params["e"]
    best_b, (a1, c1, r2_1, yp1, dcp1, dmax1) = _scan_best(xs, ys, b1, b2, e)
    b_use = round(best_b * 100.0) / 100.0
    a2, c2, r2_2, yp2, dcp2, dmax2 = _fit(xs, ys, b_use)

    sections = []
    # (一) 基本数据
    rows = ["   N        X(I)        Y(I)"]
    for i, (x, y) in enumerate(zip(xs, ys)):
        rows.append(f"  {i:>3}   {x:>10.3f}   {y:>10.3f}")
    sections.append(("一", ["基本数据"] + rows))

    # (二) b 优选段
    rows2 = ["*****  输入b值优选范围和精度要求  *****",
             f"幂函数中未知数系数 a = {a1:>10.5f}",
             f"幂函数中未知数指数 b = {best_b:>10.5f}",
             f"幂函数中常数项     c = {c1:>10.5f}",
             f"相关分析中相关系数平方 r*r= {r2_1:>8.5f}"]
    for i, v in enumerate(yp1):
        rows2.append(f"Y({i:>2} )= {v:>8.3f}")
    rows2.append(f"回代计算中平均误差  Dcp= {dcp1:>8.5f}")
    rows2.append(f"回代计算中最大误差 Dmax= {dmax1:>8.5f}")
    rows2.append(f"Y = {a1:>7.5f} X^ {best_b:>7.5f}+ {c1:>8.5f}")
    sections.append(("二", rows2))

    # (三) 确定 b 段
    rows3 = ["****   确定b值重新计算 ，打印采用成果   *****",
             f"幂函数中未知数系数 a = {a2:>10.5f}",
             f"幂函数中未知数指数 b = {b_use:>10.5f}",
             f"幂函数中常数项     c = {c2:>10.5f}",
             f"相关分析中相关系数平方 r*r= {r2_2:>8.5f}"]
    for i, v in enumerate(yp2):
        rows3.append(f"Y({i:>2} )= {v:>8.3f}")
    rows3.append(f"回代计算中平均误差  Dcp= {dcp2:>8.5f}")
    rows3.append(f"回代计算中最大误差 Dmax= {dmax2:>8.5f}")
    rows3.append(f"Y = {a2:>7.5f} X^ {b_use:>7.5f}+ {c2:>8.5f}")
    sections.append(("三", rows3))

    return render_text(PROGRAM_ID, TITLE, sections)


def run(data, out_txt=None, out_json=None, fmt="text"):
    """
    统一入口。
    data: INT 文件路径 | 数值流 | dict
    返回 (result_dict, 计算书文本)
    """
    if isinstance(data, str):
        params = parse(read_numbers(data))
    else:
        params = parse(data)
    result = compute(params)
    if fmt == "markdown":
        text = render_markdown(PROGRAM_ID, TITLE, [("一", ["结果见 JSON"])], result)
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
