# -*- coding: utf-8 -*-
"""
D-28 压力供水管道过水能力计算程序 —— 内核
==========================================
复刻《水利水电工程设计计算程序集》公之于众版 D-28（压力供水管道过水能力）。

功能
----
已知管道各段桩号、管内径、管壁类型与**上下游水位差 Zo**，求管道可通过的流量 Q
（即「过水能力」），并逐段打印管段长 L、流量 Q、内径 D、管壁类型 N、流速 V、
流速水头 V²/2g 与水头损失和。

原著算法（由权威 D-28.OUT 反演，含反证）
----------------------------------------
流量 Q 为**解**：使全线沿程水头损失之和恰等于给定水位差 Zo
        Σ_j  i(材料_j, Q, d_j) · L_j = Zo
（D-28.OUT 逐段损失之和 = 341.24 m（按 2 位打印值累加）与 Zo = 341.21 m 一致，
  即**局部水头损失百分比未另行叠加**，或以系数并入 i —— 二者不可分辨，见 docstring 末）。
i 与 D-27 同族：i = C·Q^a·d^(−b)（Q m³/s，d m）。D-28vb.exe 内嵌的指数常量
亦为 double 1.761000633 / 4.761001587（float32(1.761)/float32(4.761)），与 D-27 一致。
C 由权威 D-28.OUT 反演（与 D-27 的 C 略有 0.17% 差异，见「反演与未闭合」）。

版式：权威 D-28.OUT 首行为运行期路径回显（L:\\01\\4.1版\\SLSDK4.1\\use\\D-28.OUT，
机器相关），本内核不生成；比对自第 2 行起算（与 D-1/D-3/D-4/D-7/D-8/D-11/D-14/
D-16/D-17/D-19/D-23~D-27 同一口径）。
"""
import os

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json
from .d27 import material_params

PROGRAM_ID = "D-28"
TITLE = "压力供水管道过水能力"
AUTHOR = "张校正（新疆水利厅）"
HEAD_NAME = "D-28"

#: 圆周率（用于流速 V 与流速水头 V²/2g 显示）
PI = 3.14159
#: 重力加速度 (m/s^2)——V²/2g 打印 2 位小数，9.8/9.81 不可分辨（见 docstring）
G = 9.81
#: 管材 6（UPVC）指数（D-28vb.exe 内嵌 double 常量 1.761000633 / 4.761001587）
A_UPVC = 1.761
B_UPVC = 4.761
#: 管材 6 的系数。与 D-27 取**同一值**（两程序共用同一材料常数）：
#: 由权威 D-28.OUT 反演，命中 41 行的可行区间为 [0.00091842, 0.00091844]（含 D-27 反演值）。
C_UPVC_D28 = 0.000918437


def area(d):
    """圆管断面积 A = π·d²/4 (m²)。"""
    return PI * d * d / 4.0


def velocity(q, d):
    """流速 V = Q/A (m/s)。"""
    return q / area(d)


def fmt_station(s):
    """桩号格式化：0+000.00（km+m，两位小数）。"""
    return "%d+%06.2f" % (int(s // 1000), s - 1000 * int(s // 1000))


# ============================================================
# 解析
# ============================================================

def parse(data):
    """
    解析输入。data：dict（原样返回）| .INT 文件路径。

    .INT 结构（据 D-28E.xls「控制数据 / 节点数据」）：
      工程名,分段数,上下游水位差
      (共「分段数+1」行，按节点)
        桩号（m）,管内径(m),管壁类型,局部水头损失百分比
    """
    if isinstance(data, dict):
        return data
    lines = [ln for ln in read_lines(data) if ln.strip()]
    if not lines:
        raise ValueError("D-28 输入为空")
    head = [x.strip() for x in lines[0].split(",")]
    name = head[0]
    nseg = int(round(float(head[1])))
    zo = float(head[2])
    nodes = []
    for ln in lines[1:2 + nseg]:
        f = [x.strip() for x in ln.split(",")]
        nodes.append({
            "station": float(f[0]),
            "d": float(f[1]),
            "mat": int(round(float(f[2]))),
            "loc": float(f[3]),
        })
    return {"name": name, "nseg": nseg, "zo": zo, "nodes": nodes}


# ============================================================
# 计算
# ============================================================

def _total_loss(p, q, c6):
    """全线沿程水头损失之和（材料 6 用 c6，其余材料按 D-27 名义常数）。"""
    tot = 0.0
    for j in range(1, len(p["nodes"])):
        n0, n1 = p["nodes"][j - 1], p["nodes"][j]
        L = n1["station"] - n0["station"]
        a, b, c = material_params(n1["mat"])
        if int(n1["mat"]) == 6:
            c = c6
        tot += c * (q ** a) * (n1["d"] ** (-b)) * L
    return tot


def compute(p):
    """解 Q 使 Σi·L = Zo（二分法），并给出逐段成果。"""
    lo, hi = 1e-6, 100.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if _total_loss(p, mid, C_UPVC_D28) < p["zo"]:
            lo = mid
        else:
            hi = mid
    q = 0.5 * (lo + hi)

    rows = []
    for j, nd in enumerate(p["nodes"]):
        L = 0.0 if j == 0 else nd["station"] - p["nodes"][j - 1]["station"]
        a, b, c = material_params(nd["mat"])
        if int(nd["mat"]) == 6:
            c = C_UPVC_D28
        loss = c * (q ** a) * (nd["d"] ** (-b)) * L
        rows.append({"I": j, "station": nd["station"], "L": L, "Q": q,
                     "D": nd["d"], "N": nd["mat"], "loss": loss,
                     "V": velocity(q, nd["d"])})
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR, "输入": p,
            "Q": q, "Zo": p["zo"], "分段数": p["nseg"], "成果": rows,
            "总损失": _total_loss(p, q, C_UPVC_D28)}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式，自第 2 行起）
# ============================================================

_LINE = " " + "*" * 71
_TITLE = " **********            压力供水管道过水能力计算书             **********"
_HDR1 = "         终止                                                      水头"
_HDR2 = "节点号            管段长     流量     内径   管壁  流速  流速水头  损失"
_HDR3 = " I       桩号       L          Q        D      N     V    V^2/2g    和 "


def _row_fmt(r):
    """数据行格式（D-28 成果表，71 字符）。"""
    v = r["V"]
    return ("%3d%12s%9.2f%11.6f%8.3f%2d%10.2f%8.2f%8.2f"
            % (r["I"], fmt_station(r["station"]), r["L"], r["Q"], r["D"],
               r["N"], v, v * v / (2.0 * G), r["loss"]))


def render(p, r):
    """生成原著风格文本计算书（.OUT）。返回自第 2 行（空行）起的文本。"""
    L = []
    L.append("")
    L.append(_LINE)
    L.append(_TITLE)
    L.append(_LINE)
    L.append("")
    L.append("              工程名  %s" % p["name"])
    L.append("")
    L.append("              管道分段数 M= %2d              上下游水位差  Zo= %.2f"
             % (r["分段数"], r["Zo"]))
    L.append("")
    L.append("                             分段计算结果")
    L.append("                             ============")
    L.append("")
    L.append(_HDR1)
    L.append(_HDR2)
    L.append(_HDR3)
    L.append("")
    for row in r["成果"]:
        L.append(_row_fmt(row))
    L.append("")
    return "\n".join(L)


def run(data, out_txt=None, out_json=None):
    """统一入口。data：INT 文件路径 | dict。"""
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
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
