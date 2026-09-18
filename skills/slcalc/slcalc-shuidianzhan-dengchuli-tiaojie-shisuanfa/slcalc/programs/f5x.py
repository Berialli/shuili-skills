# -*- coding: utf-8 -*-
"""
F-5X 大规模三角网平差计算程序 —— 内核
====================================
复刻《水利水电工程设计计算程序集》F-5X（作者：谢希哲，新疆兵团勘测设计院）。

F-5XIntro 第 1 条：「与『任意三角网平差计算程序(F-5)』的差别仅在于解算量的扩大」——
至多 450 个控制点、总方向数至多 2500、可加测 5 个方位角及 400 条边（F-5 为
130 点 / 700 方向 / 5 个方位角 / 120 条边）；法方程规模扩大后不整体定义数组，
改用 C 盘临时随机文件，仍按平方根法求解。

因此本内核 = F-5 内核（`slcalc/programs/f5.py`）的同一算法，仅
  · PROGRAM_ID / 注册键名 为 "F-5X"；
  · 渲染走 F-5X 版式变体（未知点行「°」后多一个全角空格、末行补 6 个半角空格）——
    该差异由权威 F-5X.OUT 与 F-5.OUT 逐字节比对确定（其余 51 行完全相同）。

基准：slcalc/data/F-5X.OUT（5350 B，52 行，GB18030）—— 见 f5x_verify.py。
"""
from . import f5

PROGRAM_ID = "F-5X"
TITLE = "F-5X 大规模三角网平差计算"
AUTHOR = "谢希哲(新疆兵团勘测设计院)"
HEAD_NAME = "F-5X"

parse = f5.parse
compute = f5.compute


def render(p, res):
    return f5.render(p, res, variant="F-5X")


def run(data, out_txt=None, out_json=None, fmt="text", **kw):
    params = f5.parse(data)
    result = f5.compute(params)
    text = f5.render(params, result, variant="F-5X")
    if out_txt:
        f5.write_out(out_txt, text)
    if out_json:
        f5.write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    _r, _t = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_t)
