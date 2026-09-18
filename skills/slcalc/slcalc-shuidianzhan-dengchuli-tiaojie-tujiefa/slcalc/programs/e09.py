# -*- coding: utf-8 -*-
"""
E-9 直接剪切试验计算程序 —— 内核
=================================
复刻《水利水电工程设计计算程序集》E-9（作者：邓铭江，新疆水利厅）。

功能
----
应变控制式直剪仪上快剪/固结快剪/慢剪试验的成果整理与绘图数据生成：
逐级垂直压力 P 下给出剪应力 f 与剪切位移 dL，取各级最大剪应力作为抗剪强度，
再以库仑强度线（最小二乘）求凝聚力 C 与内摩擦角 B。

公式（E-9Intro「二、计算原理 计算公式」原文）
---------------------------------------------
    剪应力        f  = C · R              (kg/cm²)
    剪切位移      dL = 20·N − R           (0.01 mm)
    库仑强度线    S  = K · P + C0         （最小二乘法）
    内摩擦角      B  = arctg K            （以度-分-秒打印）
    相关数        r  = 线性相关系数

其中 C 为量力环率定系数（kg/cm²/0.01mm）、R 为量力环量表读数（0.01mm）、
N 为手轮转数（N = n·i，n 为「每转 n 周测记一次」，i 为读数序号）。

输入数据（E-9Intro「（三）数据文件顺序」）
------------------------------------------
    y,n,M                      试样个数, 读数间隔转数, 每试样读数次数
    对第 k 个试样（k = 1..y）：
        P(k),C(k)              垂直压力(kg/cm²), 量力环率定系数
        R(1),R(2),...,R(M)     量力环量表读数(0.01mm)（可跨行）

输出（逐字复刻原著 .OUT 版式）
------------------------------
  每个试样一页：表头 → 参数回显 → N/R/L/F 表 → 最大剪应力 Pmax
  末页附：凝聚力 C、内摩擦角 B（度-分-秒）、相关数 r

常量口径
--------
  · 手轮转速固定打印「6转/分」（原著为表头固定文本，非由 n 推算；见
    E-9Intro 记录表样例「手轮转度 6转/分」，与 n=2 无算术关系）
  · 位移系数 20 = 手轮每转对应 0.02mm×1000，即 20×N 单位为 0.01mm（原著公式原文）
  · 打印取整：十进制四舍五入（VB6 Format 口径）

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 库仑强度线 τf = c + σ·tanφ | `01_水工设计手册精读/卷1_基础理论.md` 4.7 强度「库仑定律 τf=c+σtanφ」；`卷6_土石坝.md` 式1.12-2/1.12-3 | **一致** | 采用 |
  | 直剪试验（快剪/固结快剪/慢剪）为 c、φ 的测定方法 | 卷1 4.7.7「不同试验方法对应不同工况，E 类程序输出即此类指标」；`02_水利教材精读/水能规划与岩土结构/土质学与土力学_第5版_精读笔记.md` 抗剪强度章 | **一致** | 采用 |
  | 抗剪强度取峰值 / 稳定值（无峰值时取位移 4mm 处） | 卷1 4.7.7；教材抗剪强度测定节 | **一致**（本程序自动取各级最大剪应力，即峰值口径） | 采用 |
  | 最小二乘拟合直线求 c、φ | 卷10 2.3「数据统计」；教材强度指标整理 | **一致** | 采用（K=Σ(x−x̄)(y−ȳ)/Σ(x−x̄)²，C0=ȳ−K·x̄，r 为标准相关系数） |

未闭合点
--------
  无。权威 E-9.OUT 全部数值（4 试样的 21×4=84 组 N/R/L/F、4 个 Pmax、
  C=0.365、B=26°39′24″、r=0.9976）均可由上述公式逐位复现。
"""
import math
from decimal import Decimal, ROUND_HALF_UP

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "E-9"
TITLE = "直接剪切试验计算程序"
AUTHOR = "邓铭江"
HEAD_NAME = "E-9"

# 手轮位移系数：dL(0.01mm) = 20·N − R（原著式）
WHEEL_FACTOR = 20.0
# 表头固定文本：手轮转速（原记录表固定 6 转/分）
HANDWHEEL_SPEED = "6转/分"

_RULE = "_" * 70


def q(x, d):
    """十进制四舍五入（VB6 Format 口径，半值进位）。"""
    if x is None:
        return 0.0
    quant = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(quant, rounding=ROUND_HALF_UP))


def F(x, w, d):
    return "%*.*f" % (w, d, q(x, d))


def dms(deg):
    """度 → (度, 分, 秒)，秒四舍五入到整数并进位。"""
    d = int(math.floor(deg))
    rem = (deg - d) * 60.0
    m = int(math.floor(rem))
    s = int(q((rem - m) * 60.0, 0))
    if s >= 60:
        s -= 60
        m += 1
    if m >= 60:
        m -= 60
        d += 1
    return d, m, s


# ============================================================
# 解析
# ============================================================

def parse(data):
    """解析 .INT 数据。返回 dict(y, n, M, samples=[{P,C,R:[...]}, ...])。"""
    nums = read_numbers(data) if not isinstance(data, dict) else None
    if isinstance(data, dict):
        return data
    if len(nums) < 3:
        raise ValueError("E-9 输入数据过少（%d 个数）" % len(nums))
    y = int(round(nums[0]))
    n = int(round(nums[1]))
    M = int(round(nums[2]))
    pos = 3
    samples = []
    for _ in range(y):
        if pos + 2 + M > len(nums) + 0 and pos + 2 > len(nums):
            raise ValueError("E-9 数据不足：第 %d 个试样缺少参数" % (len(samples) + 1))
        P = nums[pos]
        C = nums[pos + 1]
        pos += 2
        if pos + M > len(nums):
            raise ValueError("E-9 数据不足：第 %d 个试样缺 R 值" % (len(samples) + 1))
        R = nums[pos:pos + M]
        pos += M
        samples.append({"P": P, "C": C, "R": R})
    return {"y": y, "n": n, "M": M, "samples": samples}


# ============================================================
# 计算
# ============================================================

def _fit(xs, ys):
    """最小二乘直线 y = K·x + C0 及相关系数 r。"""
    m = len(xs)
    mx = sum(xs) / m
    my = sum(ys) / m
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        raise ValueError("E-9 库仑线拟合退化：垂直压力或抗剪强度全相同")
    K = sxy / sxx
    C0 = my - K * mx
    r = sxy / math.sqrt(sxx * syy)
    return K, C0, r


def compute(p):
    n = p["n"]
    out_samples = []
    for smp in p["samples"]:
        C, P, Rs = smp["C"], smp["P"], smp["R"]
        rows = []
        for i, R in enumerate(Rs):
            N = n * (i + 1)
            L = WHEEL_FACTOR * N - R
            f = C * R
            rows.append({"N": N, "R": R, "L": L, "F": f})
        # 原著将各级最大剪应力按打印量子（0.01 kg/cm²）取整后再作库仑线拟合
        pmax = max(q(r["F"], 2) for r in rows)
        out_samples.append({"P": P, "C": C, "rows": rows, "Pmax": pmax,
                            "n_read": len(Rs)})
    xs = [s["P"] for s in out_samples]
    ys = [s["Pmax"] for s in out_samples]
    K, C0, r = _fit(xs, ys)
    B = math.degrees(math.atan(K))
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": p, "试样": out_samples, "K": K, "凝聚C": C0,
            "内摩擦角B": B, "相关数r": r}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def _put(s, token, end):
    """把 token 追加到 s 之后，使其最后一个字符位于 0 基列号 end。"""
    need = end + 1 - len(token)
    if len(s) < need:
        s = s + " " * (need - len(s))
    return s + token


def render(p, r):
    y = p["y"]
    A = []
    A.append("")
    A.append(" ***********************************************************************")
    A.append(" ****                        直剪试验计算书 E-9                     ****")
    A.append(" ***********************************************************************")
    A.append("")
    A.append("")
    for k, smp in enumerate(r["试样"]):
        if k > 0:
            A.append("")
        A.append("    工程名称_____________             试验者___________")
        A.append("")
        A.append("    土样编号_____________             计算者___________")
        A.append("")
        A.append("    试验日期_____________             校核者___________")
        A.append(_RULE)
        A.append("")
        A.append("    试样编号_____________             剪前压缩时间___________")
        A.append("")
        A.append("    仪器编号_____________             剪前压缩沉降量___________")
        A.append("")
        A.append("    量力环系数C(%2d )= %s(kg/cm^2/0.01mm)" % (k + 1, F(smp["C"], 0, 4)))
        A.append("")
        A.append("    手轮转速  %s                  抗剪强度___________" % HANDWHEEL_SPEED)
        A.append("")
        A.append("    垂直压力P(%2d )= %s(kg/cm^2)      剪切历时___________"
                 % (k + 1, F(smp["P"], 0, 1)))
        A.append(_RULE)
        A.append("")
        A.append("              N         R         L         F         V")
        A.append("             (-)      (0.01mm)  (0.01mm) (kg/cm^2)    mm")
        for row in smp["rows"]:
            s = " " * 14 + str(row["N"])
            s = _put(s, F(row["R"], 0, 1), 26)
            s = _put(s, F(row["L"], 0, 1), 37)
            s = _put(s, F(row["F"], 0, 2), 45)
            A.append(s)
        A.append("")
        A.append("  最大剪应力Pmax= %s(kg/cm^2)" % F(smp["Pmax"], 0, 2))
        A.append("")
    A.append("  凝聚力C=%7s(kg/cm^2)" % F(r["凝聚C"], 0, 3))
    dd, mm, ss = dms(r["内摩擦角B"])
    A.append("  内摩擦角B= %2d %2d %2d(度-分-秒)" % (dd, mm, ss))
    A.append("  相关数r= %s" % F(r["相关数r"], 0, 4))
    return "\n".join(A) + "\n"


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
