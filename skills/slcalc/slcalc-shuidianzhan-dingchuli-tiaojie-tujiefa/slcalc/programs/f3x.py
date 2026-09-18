# -*- coding: utf-8 -*-
"""
F-3x 大规模导线网平差计算程序 —— 内核
======================================
复刻《水利水电工程设计计算程序集》F-3x 程序（作者：谢希哲，新疆兵团勘测设计院）。
原著说明书（RTF\\F-3xIntro.rtf）「一. 简介」：

    「与『导线网平差计算程序(F-3)』的差别**仅在于解算量的扩大**，除相同之处外，
      还补充本程序扩大解算量所采取的措施。
      2. 至多可有 250 个已知点及结点，各点的方向数不限；线路至多 400 条，
         每条线路的边数不限；
      6. 计算方法：定权原则与单一附合导线同，但采用相关平差法。
         为了对付庞大的数据，F-3 程序即已将最大的法方程阵取其上三角阵并进行了
         变带宽的处理。本程序则更进一步，将法方程阵改在 C 盘建立一个随机文件来
         存放其数据，将原来近 10000 个单元的数组变量所占的内存让出来……
         当结点数在 20 个以内时，仍以采用 F-3 程序为宜。
      三. 操作注意事项：与「导线网平差计算程序（F-3）」**完全一样，数据文件均可兼容**。」

故本内核与 F-3 **同算法、同输入契约、同输出格式**，仅**容量上限独立放大**
（`MAX_KP=250`、`MAX_LINE=400`），并独立声明本枚常量（不跨程序复用）。

★ 基准闭合状态（如实标注）
------------------------
  与 F-3 相同：G 盘 `RTF\\算例计算结果文件\\F\\` 下**无 F-3x .OUT**；
  F-3xIntro.rtf 全文逐字提取核对后确认**未刊印任何算例结果**；
  EXE 目录下**无 F-3xvb.EXE**。故本内核**全部输出计 DECL**，
  以合成网自洽（条件闭合差 → 0）与原始三件（.INT/.C/.S）结构可达性作为量化证据。

知识库对照结果：与 F-3 同源（母本《水利工程测量》第5版精读笔记），
  逐条对照见 `slcalc/programs/f3.py` 文件末「知识库对照结果」节与 `_f3_verify_run.txt`。

知识产权：本改造工作基于公之于众版（乌鲁木齐正海水利科技有限公司，张校正教授级高工
技术总负责）的公开算法。改造实现（Python 代码、架构设计、验证数据、自动化流程）为
哈胜的原创成果。
"""
from . import f3 as _f3

PROGRAM_ID = "F-3x"
TITLE = "大规模导线网平差计算程序"
AUTHOR = "谢希哲"
HEAD_NAME = "F-3x"

# 本枚独立的容量上限（原著 F-3xIntro 第 2 条：250 个已知点及结点、400 条线路）
MAX_KP = 250
MAX_JD = 250
MAX_LINE = 400
NAME_LEN = 12

RHO = _f3.RHO
dmss_to_deg = _f3.dmss_to_deg
deg_to_dmss = _f3.deg_to_dmss
fmt_dms = _f3.fmt_dms
parse_c = _f3.parse_c
parse_s = _f3.parse_s
parse_int = _f3.parse_int


def chk_capacity(p):
    """本枚容量校验（独立阈值）。"""
    nk = len(p.get("已知点", []))
    nl = len(p.get("线路", []))
    nj = sum(max(0, len(l["节点"]) - 2) for l in p.get("线路", []))
    if nk + nj > MAX_KP:
        raise ValueError("%s 已知点+结点 %d 超过上限 %d（原著简介第 2 条）"
                         % (PROGRAM_ID, nk + nj, MAX_KP))
    if nl > MAX_LINE:
        raise ValueError("%s 线路数 %d 超过上限 %d" % (PROGRAM_ID, nl, MAX_LINE))


def parse(data):
    """与 F-3 同一契约；仅容量校验换用本枚阈值。"""
    if isinstance(data, dict):
        p = dict(data)
    else:
        p = parse_int(data)
        if not p.get("可解"):
            raise ValueError(
                "F-3x .INT 基本数据未能唯一解析（原著随机文件格式未公开）。"
                "请改用显式输入契约（dict/JSON，见 f3 模块 docstring）。")
    p.setdefault("名称", "未命名导线网")
    p.setdefault("计算者", "AI")
    p.setdefault("日期", "")
    p.setdefault("Ms_a", 10.0)
    p.setdefault("Ms_b", 6.0)
    p.setdefault("方向组", [])
    kp2 = []
    for kp in p.get("已知点", []):
        if isinstance(kp, dict):
            kp2.append({"名": str(kp["名"]), "X": float(kp["X"]), "Y": float(kp["Y"]),
                        "方向": list(kp.get("方向", []))})
        else:
            kp2.append({"名": str(kp[0]), "X": float(kp[1]), "Y": float(kp[2]),
                        "方向": list(kp[3]) if len(kp) > 3 else []})
    p["已知点"] = kp2
    chk_capacity(p)
    return p


def compute(p):
    """复用 F-3 的间接平差引擎（两程序算法完全一致，原著明示「完全一样」），
    并把结果 dict 中的程序标识改回本枚 F-3x。"""
    r = _f3.compute(p)
    r["程序"] = PROGRAM_ID
    r["标题"] = TITLE
    r["作者"] = AUTHOR
    return r


def render(p, r):
    return _f3.render(p, r, banner_title="大 规 模 导 线 网 平 差", banner_pid=PROGRAM_ID)


def run(data, out_txt=None, out_json=None):
    from ..core.outgen import write_out, write_json
    params = parse(data)
    result = compute(params)
    text = render(params, result)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, _f3._jsonable(result))
    return result, text


if __name__ == "__main__":
    import sys
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
