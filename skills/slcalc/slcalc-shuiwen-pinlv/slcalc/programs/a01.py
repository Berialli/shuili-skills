# -*- coding: utf-8 -*-
"""
A-1 直线相关计算程序 —— 内核
==============================
复刻《水利程序集》A-1 程序（作者：陈丽棠，水电部天津勘测设计院）。
计算 Xi、Yi 两系列直线相关方程 a、b，相关系数 r、相关系数标准差 Kr、
机误 Er、两系列标准差 Kx、Ky，已知 X 值插补 Y 值，并可自动转换
自变量与倚变量重算。

数据文件顺序（INT）：
  系列项数
  Xi 系列值（逗号分隔）
  Yi 系列值（逗号分隔）

验证基准：A-1.INT（拉萨~唐加 17 年 9 月径流）
  r=0.990  a=77.697  b=0.682  Kr=0.005  Er=0.003
  Kx=205.400  Ky=141.391  meanX=596.706  meanY=484.412
"""
from ..core.numext import linear_regression
from ..core.intio import read_numbers
from ..core.outgen import render_text, render_markdown, write_out, write_json

PROGRAM_ID = "A-1"
TITLE = "直线相关计算书 A-1"

# 说明文档中的插补 X 值（算例固定输出 4 个插补点）
_DEFAULT_INTERP = [502.0, 610.0, 442.0, 526.0]


def parse(data):
    """
    解析输入。
    data: 数值流（list[float]）或 dict。
    返回 dict(xs=[...], ys=[...], x_interp=[...]|None)
    """
    if isinstance(data, dict):
        xs = [float(v) for v in data.get("xs", [])]
        ys = [float(v) for v in data.get("ys", [])]
        interp = data.get("x_interp")
        if interp is not None:
            interp = [float(v) for v in interp]
        if not xs or len(xs) != len(ys):
            raise ValueError("xs 与 ys 必须等长且非空")
        return {"xs": xs, "ys": ys, "x_interp": interp}

    # 数值流：n, x1..xn, y1..yn
    nums = list(data)
    if len(nums) < 3:
        raise ValueError("数据不足：需要 n, Xi 系列, Yi 系列")
    n = int(nums[0])
    if len(nums) < 1 + 2 * n:
        raise ValueError(f"数据不足：声明 n={n}，实际仅有 {len(nums)-1} 个数值")
    xs = nums[1:1 + n]
    ys = nums[1 + n:1 + 2 * n]
    return {"xs": xs, "ys": ys, "x_interp": _DEFAULT_INTERP}


def compute(params):
    """执行直线相关计算（含自变量/倚变量互换两轮）。"""
    xs, ys = params["xs"], params["ys"]
    x_interp = params.get("x_interp")

    r1 = linear_regression(xs, ys, x_interp)
    r2 = linear_regression(ys, xs)  # 互换

    return {
        "程序": PROGRAM_ID,
        "n": r1["n"],
        "第一轮(X→Y)": {
            "r": round(r1["r"], 3), "a": round(r1["a"], 3), "b": round(r1["b"], 3),
            "Kr": round(r1["Kr"], 3), "Er": round(r1["Er"], 3),
            "Kx": round(r1["Kx"], 3), "Ky": round(r1["Ky"], 3),
            "mean_x": round(r1["mean_x"], 3), "mean_y": round(r1["mean_y"], 3),
            "插补": [round(v, 3) for v in r1["y_interp"]] if r1["y_interp"] else None,
        },
        "第二轮(Y→X互换)": {
            "r": round(r2["r"], 3), "a": round(r2["a"], 3), "b": round(r2["b"], 3),
            "Kr": round(r2["Kr"], 3), "Er": round(r2["Er"], 3),
            "Kx": round(r2["Kx"], 3), "Ky": round(r2["Ky"], 3),
            "mean_x": round(r2["mean_x"], 3), "mean_y": round(r2["mean_y"], 3),
        },
    }


def render(params, result):
    """生成计算书文本。"""
    xs, ys = params["xs"], params["ys"]
    x_interp = params.get("x_interp") or []
    r1 = linear_regression(xs, ys, x_interp)
    r2 = linear_regression(ys, xs)

    sections = []
    # (一) 原始资料
    rows = [" I       X(I)        Y(I)       I         X(I)         Y(I)",
            "-" * 62]
    half = (len(xs) + 1) // 2
    for i in range(half):
        j = i + half
        if j < len(xs):
            rows.append(f"{i+1:>3} {xs[i]:>10.3f} {ys[i]:>10.3f}   "
                        f"{j+1:>3} {xs[j]:>10.3f} {ys[j]:>10.3f}")
        else:
            rows.append(f"{i+1:>3} {xs[i]:>10.3f} {ys[i]:>10.3f}")
    sections.append(("一", ["原始资料"] + rows))

    def block(res, label):
        b = [label]
        b.append(f"  相关系数 r= {res['r']:.3f}")
        b.append(f"  回归直线方程 a= {res['a']:>10.3f}          b= {res['b']:.3f}")
        b.append(f"    Y= {res['a']:>10.3f}+ {res['b']:.3f}X")
        b.append(f"  相关系数标准差 Kr= {res['Kr']:.3f}           相关系数的机误 Er= {res['Er']:.3f}")
        b.append(f"  Yi系列标准差   Ky= {res['Ky']:.3f}         Xi系列标准差   Kx= {res['Kx']:.3f}")
        b.append(f"         ** 均  值 **")
        b.append(f"  X= {res['mean_x']:>8.3f}                        Y= {res['mean_y']:>8.3f}")
        if label.startswith("(二)"):
            b.append("         ********     插补 y 值        *******")
            b.append("   X             Y")
            if res["y_interp"]:
                for xv, yv in zip(x_interp, res["y_interp"]):
                    b.append(f"{xv:>10.3f}        {yv:>10.3f}")
        return b

    sections.append(("二", block(r1, "(二) 计算结果:")))
    sections.append(("三", ["********    自变量x与倚变量y互换后，重新计算上述各值。  ********"]))
    sections.append(("四", block(r2, "(二) 计算结果:")))
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
        text = render_markdown(PROGRAM_ID, TITLE, [("一", ["原始资料"]), ("二", ["结果见 JSON"])], result)
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
