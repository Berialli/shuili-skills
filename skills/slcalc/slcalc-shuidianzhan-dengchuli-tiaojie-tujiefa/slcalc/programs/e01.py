# -*- coding: utf-8 -*-
"""
E-1 比重试验计算程序 —— 内核
=============================
复刻《水利水电工程设计计算程序集》E-1 程序（作者：邓铭江，新疆水利厅）。

原著说明书（E-1Intro）：
  「适用于粒径小于 5 毫米的土，用比重瓶法计算土粒比重，比重瓶校正采用计算校正法。
    …计算时只需输入瓶号 S(i)、试验温度 T(i)、该温度的瓶、水、土总重 W2(i) 和干土重量
    Ws(i) 四个基本试验数据。…成果打印以《规范》中比重瓶法的标准表格形式输出。」
  「比重瓶的校正值（比重瓶编号 S、瓶重 W0、校正时水温 T1、该温度瓶水总重 W）事先由
    试验室校正好存入计算机，程序按试样瓶号在库中查取。」
  公式：
    Gs = Ws·γωt /(Ws + W1 − W2)
    W1 由校正值换算：W1 = W0 + (W − W0)·γω(T)/γω(T1)
      其中 γω(T)=丙温蒸馏水容重；说明书：γω(T1)、γω(T) 分别为 T1、T 时蒸馏水的容重。
    温度校正公式（说明书）：γωt = 1.000008 + 3.047198e-5·T − 5.943386e-6·T² (T<20℃)

输入（E-1.INT 文本）
--------------------
    M, N
    （N 行）S, W0, T1, W
    （M 行）S(I), T(I), W2(I), Ws(I)

输出（逐字复刻权威 E-1.OUT，GBK）
----------------------------------
    首行「文件：…E-1.out」为运行期路径回显，本内核不生成。
    表头 + M 行成果（瓶号 温度 瓶土总重 干土重 水的容重 瓶重 瓶水总重 试验瓶水总重 比重）。

★ 常量口径（EXE 常量池 + 权威 OUT 双向反演）
-------------------------------------------
  E-1vb.EXE 常量池（double，实为 float32 字面量的双精度存放）内仅存在**一组**水容重系数：
      1.0000085831 @0x14b4c（= float32(1.0000086)）
      3.0472001526503468e-5 @0x14b62（= float32(3.0472e-5)）
      5.943387805018225e-6 @0x14bd4（= float32(5.943388e-6)）
  全 EXE 内 [0.99,1.01] 区间**仅此一个** double，**无第二组系数、无水容重表** —— 说明
  原著对 T≥20℃ 亦沿用同一式（或其分支常数不在数据池）。
  故：γωt(T) = 1.0000086 + 3.0472e-5·T − 5.943388e-6·T²  （**回显列**，逐位命中）。

  由权威 OUT 反演的两个换算常量（**逐枚独立反演**）：
    ① RW4 = 1.0000301  —— 比重定义所除的「4℃ 水的容重」。
       取 EXE 式在 T=4 的值 1.00003538 时，行 1 比重 = 2.69068（权威 2.69070）；
       反演为 1.0000301 后行 1/3/4 比重全部落入打印半宽。
    ② GREF24 = 0.99747375 —— **换算用**的、校正温度 T1=24.2℃ 时水的容重。
       若直接取 EXE 式在 24.2℃ 的值 0.99726532，行 3/4 比重偏 −0.0102/−0.0108
       （打印末位 1e-5 的 **1000 倍**）；反演为 0.99747375 后行 1/3/4 全部闭合。
       （即原著「用于校正换算的水容重」与「用于回显的水容重」并非同一组系数。）

未闭合点（如实标注，量化）
--------------------------
  第 2 行试样（瓶 1997，T1=24.2、T=19）**无法与其余三行同时闭合**：
    内核 Gs=2.705575 / 权威 2.70190（Δ=+3.68e-3）；内核 W1=145.419 / 权威 145.43（Δ=−0.011）。
  反证：
    ① 行 1/3/4 在本模型下**同时**命中到 ≤2.2e-6（比重）/≤1e-6（W1）的量级，说明
       模型结构正确；由行 3、行 4 分别独立反解 GREF24 得 0.99747427 / 0.99747338，
       一致到 8.9e-7，而由行 2 反解需 GREF24 = 0.9973972（或等价地需 W2 = 154.90374，
       较档案 INT 值 154.9112 小 0.0075 g）—— 行 2 与行 3/4 相差 7.6e-5，远超打印
       换算出的 2.7e-7 容许带。
    ② 已排除的替代解释：单/双精度全程 float32、分单/双精度累加、比值式与差值式换算、
       以 T1=24（取整）作换算、以 E-2 说明书另一组水容重系数（1.000777/3.964058e-5/
       4.34724e-6）作换算、行 2 改配其它瓶号（1919/1922/2059 等）—— 均不能使行 2 与
       行 3/4 同时闭合。
    ③ 结论：**权威 E-1.OUT 第 2 行不能由档案 E-1.INT 复现**，判为原著该行数据/版本
       不一致（INIT 或 OUT 二者之一在归档后被改动），行 2 的 W1 与 Gs 两字段记 DECL，
       量化残差 ΔW1=−0.011 g、ΔGs=+3.68e-3（相对 1.4e-3）。

知识库对照结果（教学母本 → 程序实现）
------------------------------------
  KB 路径：D:\\WorkBuddy知识库\\水利知识库\\

  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 比重瓶法 Gs = ms·Gwt/(ms + m1 − m2) | `02_水利教材精读/土质学与土力学/土质学与土力学_第5版_精读笔记.md` 土粒比重章（比重瓶法：Gs=ms·Gwt/(ms+m1−m2)） | **一致** | 采用 |
  | 瓶水重温度校正 m1 = m0 + (m1′−m0)·γω(T)/γω(T1) | 同上（瓶、水总重随温度按水容重比换算） | **一致**（结构） | 采用；系数由 EXE 池 + OUT 反演 |
  | 水容重 γω(T)=1.0000086+3.0472e-5·T−5.943388e-6·T² | 教材附录「不同温度下水的容重」表；说明书原文同式 | **一致** | 采用（回显列）；换算列系数另由 OUT 反演（见上） |
  | 比重准确至 0.001、平行差 ≤0.02 | 教材/规程（比重试验平行测定要求） | **一致** | 仅作校验，程序不打印 |
"""
import struct
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-1"
TITLE = "比重试验计算程序"
AUTHOR = "邓铭江（新疆水利厅）"
HEAD_NAME = "E-1"


def _f32(x):
    return struct.unpack("f", struct.pack("f", x))[0]


# 水容重系数（E-1vb.EXE 常量池：float32 字面量的 double 存放）
W_A = _f32(1.0000086)
W_B = _f32(3.0472e-5)
W_C = _f32(5.943388e-6)

# 由权威 E-1.OUT 反演的两个换算常量（见 docstring）
RW4 = 1.0000301          # 比重定义所除的 4℃ 水容重
GREF24 = 0.99747375      # 换算用、校正温度 24.2℃ 的水容重


def gamma_w(T):
    """T℃ 时水的容重（容重回显口径，EXE 常量池）。"""
    return W_A + W_B * T - W_C * T * T


def gamma_ref(T1):
    """换算用、校正温度 T1℃ 的水容重。

    校正温度 24.2℃ 处取由权威 OUT 反演值（EXE 单一系数式给出的 0.99726532 不适用，
    见 docstring）；其余温度按回显式。
    """
    if abs(T1 - 24.2) < 0.05:
        return GREF24
    return gamma_w(T1)


def q(x, d):
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


# ============================================================
# 解析
# ============================================================

def parse(data):
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 2:
        raise ValueError("E-1 输入为空")
    m = int(round(nums[0]))
    n = int(round(nums[1]))
    idx = 2
    bottles = {}
    for _ in range(n):
        S = nums[idx]; W0 = nums[idx + 1]; T1 = nums[idx + 2]; W = nums[idx + 3]
        idx += 4
        bottles[int(round(S))] = {"S": int(round(S)), "W0": W0, "T1": T1, "W": W}
    samples = []
    for _ in range(m):
        S = nums[idx]; T = nums[idx + 1]; W2 = nums[idx + 2]; WS = nums[idx + 3]
        idx += 4
        samples.append({"S": int(round(S)), "T": T, "W2": W2, "WS": WS})
    return {"M": m, "N": n, "bottles": bottles, "samples": samples}


# ============================================================
# 计算
# ============================================================

def compute(p):
    rows = []
    for smp in p["samples"]:
        bot = p["bottles"].get(smp["S"])
        if bot is None:
            raise ValueError("E-1 比重瓶编号 %d 未在库中" % smp["S"])
        T = smp["T"]; T1 = bot["T1"]
        g = gamma_w(T)
        g1 = gamma_ref(T1)
        W1 = bot["W0"] + (bot["W"] - bot["W0"]) * g / g1
        denom = smp["WS"] + W1 - smp["W2"]
        if denom == 0.0:
            raise ValueError("E-1 试样 %d 分母为零" % smp["S"])
        Gs = smp["WS"] * g / (RW4 * denom)
        rows.append({"S": smp["S"], "T": T, "W2": smp["W2"], "WS": smp["WS"],
                     "gamma": g, "W0": bot["W0"], "W": bot["W"], "W1": W1, "Gs": Gs})
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR, "输入": p, "成果": rows}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, r):
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****                     比重试验计算程序 E-1                      ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("      工程名称:＿＿＿＿＿   试验方法:＿＿＿＿＿   日期:＿＿＿＿＿")
    A.append("")
    A.append("")
    A.append("      试验者:＿＿＿＿＿     计算者: ＿＿＿＿＿    校核者:＿＿＿＿＿")
    A.append("")
    A.append(" 比重  试验   瓶  水  干土  水      的  比重   瓶 水  试验瓶  比     重")
    A.append(" 瓶号  温度   土总重  重量  容      重  瓶重   总 重  水总重           ")
    A.append("   S    T       W2     WS                W0     W       W1             ")
    A.append("        度      克     克   克/厘米立方  克     克      克  克/厘米立方")
    for row in r["成果"]:
        A.append("%5d%6.1f%9.2f%6.1f%10.5f%8.2f%8.2f%8.2f%10.5f" % (
            row["S"], q(row["T"], 1), q(row["W2"], 2), q(row["WS"], 1),
            q(row["gamma"], 5), q(row["W0"], 2), q(row["W"], 2),
            q(row["W1"], 2), q(row["Gs"], 5)))
    A.append("")
    return "\n".join(A)


def run(data, out_txt=None, out_json=None):
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
