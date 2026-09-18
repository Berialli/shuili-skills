# -*- coding: utf-8 -*-
"""
D-18B Ｙ型分岔管水头损失计算程序 —— 内核
=========================================
复刻《水利水电工程设计计算程序集》公之于众版 D-18B 程序（作者：姚廉华）。

功能
----
计算**Ｙ型分岔管**在发电（分流）与抽水（汇流）两种工况下两条支管的水头
损失。程序参照《水力摩阻》一书编写（说明书「一、程序功能」原文）。

公式（说明书「二、计算公式」内嵌 Equation 对象解包）
------------------------------------------------------
说明书 D-18BIntro.rtf 含 3 个 OLE 公式对象，解包（\\objdata→hex→OLE2→
"Equation Native"→28B 头→MTEF v3 令牌流）得：

  对象[00]  H  = S·V²/(2g)                                   （水头损失）
  对象[01]  S₃₋₁ = 1 − 2·v₁/v₃·cos φ₁ + (v₁/v₃)²             （1# 支管损失系数）
  对象[02]  S₃₋₂ = 1 − 2·v₂/v₃·cos φ₂ + (v₂/v₃)²             （2# 支管损失系数）

MTEF 令牌流（原文逐字）：
  [01] `S 3 - 1 = 1 - 2 v 1 / v 3 cos f 1 + ( v 1 / v 3 ) 2`
  [02] `S 3 - 2 = 1 - 2 v 2 / v 3 cos f 2 + ( v 2 / v 3 ) 2`
其中 v₁、v₂ 为两条支管流速，v₃ 为主管流速，f₁、f₂ 为两条支管的分岔夹角。

程序实现口径（与上式一致）
--------------------------
    断面积        A(D) = π·D²/4                （π = 3.14159）
    发电（分流）  v₁ = AL1·TQ/A₁, v₂ = AL2·TQ/A₂, v₃ = ALM·TQ/A₃
    抽水（汇流）  同上，把 TQ 换成 PQ（速度比值不变，故 X₁₁ ≡ X₁）
    系数          X₁ = 1 − 2(v₁/v₃)cos F2 + (v₁/v₃)²      （1# 支管）
                  X₂ = 1 − 2(v₂/v₃)cos F3 + (v₂/v₃)²      （2# 支管）
    水头损失      H₁ = X₁·v₁²/(2g)                        （m）
                  H₂ = X₂·v₂²/(2g)
                  H₁₁、H₂₁ 为抽水工况改用 PQ 后的对应值
    ※ 结果量在 VB6 中存为 Single（float32），打印用 8 位小数：
      浮点原型经权威 OUT 逐位反演确认（见「裁决 A」）。

输入数据（.INT，逗号分隔，10 个数）
-----------------------------------
    D1  1#支管直径 (m)         D2  2#支管直径 (m)         DM  主管直径 (m)
    AL1 1#支管过流机组数       AL2 2#支管过流机组数       ALM 主管过流机组数
    F2  分岔夹角 (度)          F3  分岔夹角 (度)
    TQ  发电工况(分流)过流量 (m^3/s)   PQ  抽水工况(汇流)过流量 (m^3/s)
原著算例 D-18B-2.INT：
    1.2,1.2,2,1,1,2,45,45,5,5

裁决（由权威 D-18B-2.OUT 逐位反演）
------------------------------------
A. **结果量按 Single(float32) 取值**。以双精度直算 X₁ = 0.964826842383，
   打印 8 位小数为 0.96482684；权威值为 **0.96482682**。取 float32(0.964826842383)
   = 0.964826822281，打印 8 位小数 = 0.96482682，与权威逐位一致。四种等价
   表达式（1−2rc+r²、(1−r)²+2r(1−c)、1−r(2c−r)、(r−c)²+(1−c²)）在 double 下
   同为 0.964826842383、在 float32 下同为 0.964826822281，故表达式形式不可分辨，
   内核取最直白的一种；**唯一性来自 float32 截断**（double 下 4 种全错 1 位）。
   反证：若按 double 打印 8 位 → 0.96482684 ≠ 0.96482682（0/1 命中）。
B. **π = 3.14159、g = 9.81 不能由本算例唯一确定**（记 DECL）：算例两条支管
   直径相同（D1=D2=1.2）、主管 2.0，且 TQ=PQ=5；X₁、X₂ 只依赖面积比
   (DM/D1)² = 25/9（π 完全约去），H₁、H₂ 只打印 2 位小数（0.96），
   π ∈ {3.14, 3.14159, 3.1416, π} × g ∈ {9.8, 9.81, 9.80665} 的全部 12 种组合
   均得 H = 0.96（H 区间 [0.9611, 0.9622]，2 位小数全为 0.96）。
   故内核取程序集 D 族主流口径 π = 3.14159、g = 9.81，并在此显式声明不可反演。
C. 版式：权威 .OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-18B-2.out」为
   原著运行期路径回显（机器相关），本内核不生成（与 D-1/D-3/D-4/D-7/D-11/
   D-15/D-16/D-17/D-19/D-25/D-27/D-28 同一口径），版式比对自第 2 行起算。

未闭合点
--------
① π、g 不可反演（量化见裁决 B）；② 权威 OUT 仅 1 个算例，X₁ 的浮点原型
   仅由该算例锁定（float32 假设与 0.96482682 唯一相容）。
"""
import struct

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-18B"
TITLE = "Ｙ型分岔管水头损失计算程序"
AUTHOR = "姚廉华"
HEAD_NAME = "D-18B"

#: 圆周率（说明书未给；取 D 族主流 5 位口径，见裁决 B）
PI = 3.14159
#: 重力加速度 (m/s^2)（同上）
G = 9.81

_SCALARS = ["D1", "D2", "DM", "AL1", "AL2", "ALM", "F2", "F3", "TQ", "PQ"]


def f32(x):
    """VB6 Single 截断（结果量按 float32 存储，见裁决 A）。"""
    return struct.unpack("<f", struct.pack("<f", float(x)))[0]


def area(d):
    """圆管断面积 A = π·D²/4 (m²)。"""
    return PI * d * d / 4.0


def coef(v_branch, v_main, angle_deg):
    """支管损失系数 S = 1 − 2(v支/v主)cos φ + (v支/v主)²（Single）。"""
    import math
    r = v_branch / v_main
    return f32(1.0 - 2.0 * r * math.cos(math.radians(angle_deg)) + r * r)


# ============================================================
# 解析
# ============================================================

def parse(data):
    """解析输入。data：dict（直接返回）| .INT 文件路径。"""
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 10:
        raise ValueError("D-18B 输入数据不足（需 10 个简单变量，实得 %d）" % len(nums))
    p = dict(zip(_SCALARS, nums[:10]))
    for k in ("D1", "D2", "DM"):
        if p[k] <= 0:
            raise ValueError("%s 必须为正（实得 %g）" % (k, p[k]))
    return p


# ============================================================
# 计算
# ============================================================

def compute(params):
    """执行 D-18B 计算：四组损失系数与水头损失。"""
    D1, D2, DM = params["D1"], params["D2"], params["DM"]
    AL1, AL2, ALM = params["AL1"], params["AL2"], params["ALM"]
    F2, F3 = params["F2"], params["F3"]
    TQ, PQ = params["TQ"], params["PQ"]

    A1, A2, AM = area(D1), area(D2), area(DM)

    out = {}
    for tag, Q in (("发电", TQ), ("抽水", PQ)):
        v1 = AL1 * Q / A1
        v2 = AL2 * Q / A2
        v3 = ALM * Q / AM
        x1 = coef(v1, v3, F2)
        x2 = coef(v2, v3, F3)
        out[tag] = {
            "Q": Q, "v1": v1, "v2": v2, "v3": v3,
            "X1": x1, "X2": x2,
            "H1": f32(x1 * v1 * v1 / (2.0 * G)),
            "H2": f32(x2 * v2 * v2 / (2.0 * G)),
        }
    g, p = out["发电"], out["抽水"]
    return {
        "程序": PROGRAM_ID,
        "标题": TITLE,
        "作者": AUTHOR,
        "输入": params,
        "断面积": {"A1": A1, "A2": A2, "AM": AM},
        "发电": g,
        "抽水": p,
        "结果": {"X1": g["X1"], "X2": g["X2"],
                 "X11": p["X1"], "X21": p["X2"],
                 "H1": g["H1"], "H2": g["H2"],
                 "H11": p["H1"], "H21": p["H2"]},
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

_LINE = " " + "*" * 71
_STAR = " " + "*" * 55
_TITLE = " ****              Ｙ型分岔管水头损失计算程序 D-18B                 ****"


def render(params, result):
    """生成原著风格文本计算书（.OUT）。返回以第 2 行（空行）开始的文本。"""
    i = result["输入"]
    R = result["结果"]
    L = []

    L.append("")
    L.append(_LINE)
    L.append(_TITLE)
    L.append(_LINE)
    L.append("")
    L.append("")
    L.append("                ORIGINAL DATA")
    L.append("                 原 始 数 据")
    L.append(_STAR)
    L.append("")
    L.append("   1#支管直径" + " " * 20 + "D1=%7.2f (m)" % i["D1"])
    L.append("   2#支管直径" + " " * 20 + "D2=%7.2f (m)" % i["D2"])
    L.append("   主管直径" + " " * 22 + "DM=%7.2f (m)" % i["DM"])
    L.append("   1#支管过流机数" + " " * 15 + "AL1=%7.2f" % i["AL1"])
    L.append("   2#支管过流机数" + " " * 15 + "AL2=%7.2f" % i["AL2"])
    L.append("   主管过流机数" + " " * 17 + "ALM=%7.2f" % i["ALM"])
    L.append("   分岔夹角" + " " * 22 + "F2=%7.2f (度)" % i["F2"])
    L.append("   分岔夹角" + " " * 22 + "F3=%7.2f (度)" % i["F3"])
    L.append("   发电工况(分流)过流量" + " " * 10 + "TQ=%12.5f (m^3/s)" % i["TQ"])
    L.append("   抽水工况(汇流)过流量" + " " * 10 + "PQ=%12.5f (m^3/s)" % i["PQ"])
    L.append("")
    L.append("")
    L.append("        RESULT OF Y SHAPE MANIFOLD HEAD LOSSES")
    L.append("         Ｙ 型 岔 管 水 头 损 失 计 算 成 果")
    L.append(_STAR)
    L.append("")
    L.append("       X1 =%12s         X2 =%12s" % ("%.8f" % R["X1"], "%.8f" % R["X2"]))
    L.append("       X11=%12s         X21=%12s" % ("%.8f" % R["X11"], "%.8f" % R["X21"]))
    L.append("")
    L.append("   1#支管发电工况水头损失值" + " " * 10 + "H1=%7.2f (m)" % R["H1"])
    L.append("   2#支管发电工况水头损失值" + " " * 10 + "H2=%7.2f (m)" % R["H2"])
    L.append("   1#支管抽水工况水头损失值" + " " * 9 + "H11=%7.2f (m)" % R["H11"])
    L.append("   2#支管抽水工况水头损失值" + " " * 9 + "H21=%7.2f (m)" % R["H21"])
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
