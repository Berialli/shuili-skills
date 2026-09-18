# -*- coding: utf-8 -*-
"""
D-26W 管道水力计算程序（W 型：树枝状管网）—— 内核
==================================================
《水利水电工程设计计算程序集》公之于众版 **D-26G / D-26W  管道水力计算程序**
（作者：张校正，新疆水利厅）的 W 型实现。

W 型与 G 型（`programs/d26g.py`）**计算内核完全相同**，区别仅在输入/输出格式：
原著说明书「使用方法（2）」原文：「在树枝状管网的计算中，每个桩号前需要用节点名
加以区别，此时计算可以使用 D-26W 程序，其区别就是输入数据中，每个桩号前，
加一个节点名。」

  · 输入：每行分段数据最前多一个「节点名」（G 型无）；
  · 输出：多一列「距离」（本段长度）与「节点名」，并在**每个以数字开头的
    节点行之后**打印一条 1 空格 + 86 连字符的分隔线；
  · 表格标题栏加宽（88 星号横幅）。

输入数据文件（.INT，逗号分隔，GBK）
-----------------------------------
  第 1 行：工程名 G$, 管道分段数 M, 起始桩号 QZ, 起始水位 Z0, 计算起始信息 K0
  第 2..M+2 行：节点名, 桩号(m), 流量 Q(m^3/s), 管内径 D(m), 管壁材料类型(1~10),
                局部阻力百分比 K, 管底高程 GD(m)

算法、公式、常量口径、裁决、知识库对照、分歧与未闭合点
------------------------------------------------------
与 `programs/d26g.py` 模块 docstring **完全一致**（同一份内核）。本模块只承载
W 型的解析与版式；所有水力要素函数（`unit_loss` / `seg_loss` / `velocity`）、
数值口径（`r2` / `chain`）、常量（`PI` / `G` / `LOCAL_K_COEF` / `PIPE_TABLE`）
与 `compute()` 均直接复用 d26g，保证两枚程序**常量口径完全同源、不重复反演**。

W 型独有的「未闭合点」
----------------------
⑦ **W 型分隔线的判据不可唯一反演**（仅 D-26W.OUT 一个权威算例）：
   权威 11 个数据行共 10 条分隔线，唯一空缺点在第 5 行（节点名 `J4`）之后。
   以「节点名首字符是否为数字」为判据，11 行**全部命中**（0 行差异）；
   「不以 J 开头」判据在本算例上与前者输出**逐字相同**（0 行差异），
   两判据不可区分（另有 6 种判据各 ≥2 行不符，见 d26g.py「未闭合点③」的反证）。
   本模块按「节点名首字符为数字」实现。
⑧ **W 型钳制行（明满流）的括号版式无权威算例可反演**：D-26W.OUT 无钳制行。
   本模块按与 G 型同构的方式处理（`(` 占流速列首分隔位、`)` 占水压高程列首分隔位，
   总宽不变 = 87 字符），标注为**推定版式**；G 型的钳制版式已由 D-26G-1.OUT 第 0 行
   逐位证实（见 d26g.py「裁决 C」）。
⑨ 权威 .OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-26W.OUT」为原著运行期
   路径回显（机器相关），本模块不生成，版式比对自第 2 行起算。
"""
from ..core.intio import read_lines
from ..core.outgen import write_out, write_json
from . import d26g as _g

PROGRAM_ID = "D-26W"
TITLE = "管道水力计算程序"
AUTHOR = "张校正"
HEAD_NAME = "D-26W"

#: 常量口径与 G 型同源（同一份内核，禁止跨程序重复反演）
PI = _g.PI
G = _g.G
LOCAL_K_COEF = _g.LOCAL_K_COEF
PIPE_TABLE = _g.PIPE_TABLE

r2 = _g.r2
chain = _g.chain
unit_loss = _g.unit_loss
seg_loss = _g.seg_loss
velocity = _g.velocity
compute = _g.compute


# ============================================================
# 解析（W 型：每行多一个节点名）
# ============================================================

def parse(data):
    """
    解析输入。data：dict（直接返回）| .INT 文件路径（GBK）。

    返回 dict：{"工程名","M","QZ","Z0","K0",
                "节点":[{"i","name","x","Q","D","N","K","GD"}, ...]}（共 M+1 个）
    """
    if isinstance(data, dict):
        return data

    lines = [ln.strip() for ln in read_lines(data)]
    lines = [ln for ln in lines if ln]
    if len(lines) < 2:
        raise ValueError("D-26W 输入数据不足：至少需 1 行控制数据 + 1 行分段数据")

    head = _g._field(lines[0])
    if len(head) < 5:
        raise ValueError("D-26W 控制数据字段不足（需 工程名,M,QZ,Z0,K0，实得 %d）" % len(head))
    name = head[0]
    try:
        M = int(float(head[1]))
        QZ = float(head[2])
        Z0 = float(head[3])
        K0 = int(float(head[4]))
    except ValueError:
        raise ValueError("D-26W 控制数据数值解析失败：%r" % (lines[0],))
    if K0 not in (1, 2):
        raise ValueError("D-26W 计算起始信息只能为 1（上游→下游）或 2（下游→上游），实得 %s" % K0)
    if M < 1:
        raise ValueError("D-26W 管道分段数 M 必须 ≥ 1，实得 %d" % M)

    nodes = []
    for j, ln in enumerate(lines[1:1 + M + 1]):
        fl = _g._field(ln)
        if len(fl) < 7:
            raise ValueError("D-26W 第 %d 行分段数据字段不足"
                             "（需 节点名,桩号,Q,D,N,K,GD）：%r" % (j + 1, ln))
        nodes.append({
            "i": j,
            "name": fl[0],
            "x": float(fl[1]),
            "Q": float(fl[2]),
            "D": float(fl[3]),
            "N": int(float(fl[4])),
            "K": float(fl[5]),
            "GD": float(fl[6]),
        })
    if len(nodes) != M + 1:
        raise ValueError("D-26W 分段数据行数不符：M=%d 需 %d 行，实得 %d 行" % (M, M + 1, len(nodes)))

    return {"工程名": name, "M": M, "QZ": QZ, "Z0": Z0, "K0": K0, "节点": nodes}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式；横幅 88 星号）
# ============================================================

_HEAD14 = "         终止                            管壁               水头   水压    地面    地面"
_HEAD15 = " 节点名            流量     距离   内径  类型 流速 流速水头 损失   高程    高程    水头"
_HEAD16 = "  I      桩号        Q               D     N    V   V^2/2g   和      Z      GD      GY "
_SEP = " " + "-" * 86


def sep_after(name):
    """是否在本行之后打印分隔线：节点名首字符为数字（见 docstring 未闭合点⑦）。"""
    return bool(name) and name[0].isdigit()


def _row(r):
    s = " %-5s" % r["name"]
    s += chain(r["x"])
    s += "%10.6f" % r["Q"]
    s += "%8.2f" % r2(r["距离"])
    s += "%7.3f" % r["D"]
    s += "%4d" % r["N"]
    if r["钳制"]:
        s += "(%5.2f" % r2(r["V"])
    else:
        s += "%6.2f" % r2(r["V"])
    s += "%7.2f" % r2(r["Vh"])
    s += "%7.2f" % r2(r["损失"])
    if r["钳制"]:
        s += ")%7.2f" % r2(r["Z"])
    else:
        s += "%8.2f" % r2(r["Z"])
    s += "%8.2f" % r2(r["GD"])
    s += "%8.2f" % r2(r["GY"])
    return s


def render(params, result, out_path=None):
    """生成原著风格文本计算书（.OUT）。返回自第 2 行（空行）起的文本。"""
    L = []
    L.append("")
    L.append(" " + "*" * 88)
    L.append(" **********" + " " * 22 + "供水管道水力学计算书" + " " * 23 + "*************")
    L.append(" " + "*" * 88)
    L.append("")
    L.append(" " * 34 + "工程名  " + result["工程名"])
    L.append("")
    L.append(" " * 34 + "管道分段数 M= %d " % result["M"])
    L.append(" " * 14 + "起始桩号:  " + chain(result["QZ"]) + " " * 25
             + "起始水位  Zo=" + ("%.2f" % r2(result["Z0"])))
    L.append("")
    L.append(" " * 34 + "分段计算结果")
    L.append(" " * 34 + "============")
    L.append("")
    L.append(_HEAD14)
    L.append(_HEAD15)
    L.append(_HEAD16)
    L.append("")
    for r in result["节点"]:
        L.append(_row(r))
        if sep_after(r["name"]):
            L.append(_SEP)
    return "\n".join(L)


def run(data, out_txt=None, out_json=None, fmt="text"):
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
