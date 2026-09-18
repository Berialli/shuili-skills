# -*- coding: utf-8 -*-
"""
D-18A 卜型分岔管水头损失计算程序 —— 内核
=========================================
复刻《水利水电工程设计计算程序集》公之于众版 D-18A 程序（作者：姚廉华）。

功能
----
计算**卜型分岔管**的引水（发电分流 / 抽水汇流）与尾水（发电汇流 / 抽水分流）
两种工况下，各岔管**直通管**与**旁支管**的水头损失表达式、本岔管总表达式，
并给出尾水分岔总系数与全站 4 个岔管的最大总水头损失值。

公式（说明书 D-18AIntro.rtf 内嵌对象解包）
-------------------------------------------
说明书含 9 个 \\objdata 块：3 个 OLE Equation（Y 型）+ 4 个公式图形（卜型）
+ 2 个插图。卜型 4 式经「WMF 预览逐字符（字体+坐标）还原」与「MTEF v3 令牌流
还原」双路解出：

  [3] 分流·旁支管  S_{c-b} = 1 + (q_b·F_c /(q_c·F_b))² − 2(q_b·F_c/(q_c·F_b))·cos φ
  [4] 分流·直通管  S_{c-n} = 1 + (q_n·F_c /(q_c·F_n))² − 2(q_n·F_c/(q_c·F_n))
  [5] 汇流·旁支管  S_{c-b} = 1 + (q_b·F_c /(q_c·F_b))²
                            − 2(q_b/q_c)²(F_c/F_n)·cos φ − 2(1−q_b/q_c)²(F_c/F_n)
  [6] 汇流·直通管  S_{c-n} = (1−q_b/q_c)²(F_c/F_n)² − (q_b·F_c/(q_c·F_b))²

其中 c=主管、b=旁支管、n=直通管；F 为断面积、q 为过流量（机组数），φ 为分岔角。
通用损失表达式（由权威 .OUT 反演）：**Z = S·q_b²/(2g·A²)**，`h = Z × (单机流量)²`。

已**逐位验证**的部分
--------------------
* 尾水·旁支管 发电（汇流式 [5]）、抽水（分流式 [3]）：
  Z41、Z4 三行全部逐位命中（见 verify）。
* 引水·旁支管 发电（分流式 [3]）：Z2 三行全部逐位命中。
* 引水·旁支管 抽水：当 q_b = q_c（末级岔管）时与发电相同，逐位命中。
* π = 3.14159、g = 9.81 由上述逐位命中唯一确定（同一套 A、2gA² 口径）。

未闭合（记 DECL，量化残差 + 反证见 verify 与 _reports/D18_report.md）
-------------------------------------------------------------------
① 引水·**直通管** Z1 / Z11：式 [4]/[6] 在「直通管机数 = 主管机数 − 旁支机数」
   口径下不能复现权威值；受控实验（见下）证明 **Z1 与 F1、DD 无关**，
   **与 DC、DN、机组数有关**，且直通机数为 0 时 Z1 ≠ 0，故 [4] 的
   「(1−X)²·q_n²」结构在直通机数口径上尚有未定环节。
   实测反证：M1=1、J=1、DC=DN=8.0 时 Z1 = 2.0173E-05 = 1/(2g·A(8)²)（S=1）；
   同口径 DC=8.0、DN=5.6 时 Z1 = 3.864E-05（S≠1），说明 S 还含 DN 的其它作用。
② 尾水·**直通管** Z31 / Z3：Z3 恒 ≈ 0 与式 [6] 一致（本算例 F_n = F_b 时 [6] 恒为 0，
   权威值 6E-08 系 Single 浮点残差）；Z31 与 [6]/[4] 均不符。
③ 汇总链 ZA/ZB/ZMAX/ZMAX1/HT1..4/HP1..4：受控实验显示 ZA 在**尾水输入完全相同**
   的两个算例中分别为 0.001 与 0.000，说明该链还耦合引水段中间量，结构未定。

受控实验（仅副本 EXE、隐藏式驱动；见报告）
-------------------------------------------
* H_F90（F1=90）、H_DD5（DD=5.0）：Z1 与基准（F1=60、DD=3.5）**完全相同**
  → Z1 与分岔角、旁支管直径无关，只与 DC/DN/机组数有关（支持式 [4] 无 cos φ）。
* H_M1J1 / P1_K1（M1=1、J=1）：给出单岔管零点对照。

输入数据（.INT，逗号分隔）
-------------------------
  第 1 行  N, M1, M2, J, TQ, PQ, F1, DD
  第 2 行  DC(1:M1)    引水分岔直通管主管直径 (m)
  第 3 行  DN(1:M1)    引水分岔直通管支管直径 (m)
  第 4 行  DC1(1:M2)   尾水分岔直通管主管直径 (m)
  第 5 行  DN1(1:M2)   尾水分岔直通管支管直径 (m)
  第 6 行  CQ1(1:M2)   尾水分岔直通管主管过流量机组数
  第 7 行  ANQ1(1:M2)  尾水分岔直通管支管过流量机组数
  第 8 行  DQ1(1:M2)   尾水分岔旁支管过流量机组数
  第 9 行  DD1(1:M2)   尾水分岔旁支管直径 (m)
  第 10 行 FF(1:M2)    尾水分岔夹角 (度)
原著算例 D-18A-1.INT：
  111,3,3,3.0,68.25,55.56,60.0,3.5 / 8.0,6.852,5.282 / 6.852,5.282,3.5 /
  5.6,5.6,8.0 / 4.0,4.0,5.6 / 2.0,2.0,4.0 / 1.0,1.0,2.0 / 1.0,1.0,2.0 /
  4.0,4.0,5.6 / 65.0,65.0,72.0

版式说明
--------
权威 .OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-18A-1.out」为原著运行期路径
回显（机器相关），本内核不生成（与 D-1/D-3/D-4/D-7/D-11/D-15/D-16/D-17/D-19/
D-25/D-27/D-28 同一口径），版式比对自第 2 行起算。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-18A"
TITLE = "卜型分岔管水头损失计算程序"
AUTHOR = "姚廉华"
HEAD_NAME = "D-18A"

#: 圆周率（由尾水/引水旁支管逐位命中反演确定）
PI = 3.14159
#: 重力加速度 (m/s^2)
G = 9.81


def area(d):
    """圆管断面积 A = π·D²/4 (m²)。"""
    return PI * d * d / 4.0


def _z(s, q, a):
    """通用损失表达式 Z = S·q²/(2g·A²)（h = Z × 单机流量²）。"""
    return s * q * q / (2.0 * G * a * a)


def s_fork_branch(qb, qc, fc, fb, phi_deg):
    """分流·旁支管损失系数（式 [3]，已逐位验证）。"""
    x = (qb * fc) / (qc * fb)
    return 1.0 + x * x - 2.0 * x * math.cos(math.radians(phi_deg))


def s_fork_straight(qn, qc, fc, fn):
    """分流·直通管损失系数（式 [4]）。"""
    x = (qn * fc) / (qc * fn)
    return 1.0 + x * x - 2.0 * x


def s_merge_branch(qb, qc, fc, fb, fn, phi_deg):
    """汇流·旁支管损失系数（式 [5]，已逐位验证）。"""
    x = (qb * fc) / (qc * fb)
    r = qb / qc
    return (1.0 + x * x
            - 2.0 * r * r * (fc / fn) * math.cos(math.radians(phi_deg))
            - 2.0 * (1.0 - r) ** 2 * (fc / fn))


def s_merge_straight(qb, qc, fc, fb, fn):
    """汇流·直通管损失系数（式 [6]）。"""
    x = (qb * fc) / (qc * fb)
    r = qb / qc
    return (1.0 - r) ** 2 * (fc / fn) ** 2 - x * x


# ============================================================
# 解析
# ============================================================

_SCALARS = ["N", "M1", "M2", "J", "TQ", "PQ", "F1", "DD"]
_KEYS = ["DC", "DN", "DC1", "DN1", "CQ1", "ANQ1", "DQ1", "DD1", "FF"]


def parse(data):
    """解析输入。data：dict（直接返回）| .INT 文件路径。"""
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 8:
        raise ValueError("D-18A 输入数据不足（需 8 个简单变量，实得 %d）" % len(nums))
    p = dict(zip(_SCALARS, nums[:8]))
    for k in ("M1", "M2"):
        p[k] = int(round(p[k]))
    M1, M2 = p["M1"], p["M2"]
    if M1 <= 0 or M2 <= 0:
        raise ValueError("M1/M2 必须为正（实得 %d/%d）" % (M1, M2))
    need = 8 + 2 * M1 + 7 * M2
    if len(nums) != need:
        raise ValueError(
            "D-18A 输入数据个数不符：需 8+2M1+7M2 = %d，实得 %d" % (need, len(nums)))
    q = 8
    for key, cnt in (("DC", M1), ("DN", M1), ("DC1", M2), ("DN1", M2), ("CQ1", M2),
                     ("ANQ1", M2), ("DQ1", M2), ("DD1", M2), ("FF", M2)):
        p[key] = nums[q:q + cnt]
        q += cnt
    return p


# ============================================================
# 计算
# ============================================================

def compute(params):
    """执行 D-18A 计算。"""
    M1, M2 = params["M1"], params["M2"]
    J = int(round(params["J"])) if float(params["J"]).is_integer() else params["J"]
    TQ, PQ = params["TQ"], params["PQ"]
    F1, DD = params["F1"], params["DD"]
    DC, DN = params["DC"], params["DN"]
    DC1, DN1, CQ1 = params["DC1"], params["DN1"], params["CQ1"]
    ANQ1, DQ1, DD1, FF = params["ANQ1"], params["DQ1"], params["DD1"], params["FF"]
    Fb = area(DD)

    # ---- 引水卜型 ----
    head_t, head_p = [], []
    for i in range(M1):
        u = M1 - i                 # 该岔主管过流机组数
        qc, qb, qn = u, 1.0, u - 1.0
        Fc, Fn = area(DC[i]), area(DN[i])
        st = s_fork_straight(qn, qc, Fc, Fn)
        sb = s_fork_branch(qb, qc, Fc, Fb, F1)
        mb = s_merge_branch(qb, qc, Fc, Fb, Fn, F1)
        mt = s_merge_straight(qb, qc, Fc, Fb, Fn)
        head_t.append({"i": i + 1, "Z1": _z(st, qn, Fn), "Z2": _z(sb, qb, Fb)})
        head_p.append({"i": i + 1, "Z11": _z(mt, qn, Fn), "Z21": _z(mb, qb, Fb)})
    for r in head_t:
        r["G1"] = r["Z1"] + r["Z2"]
    for r in head_p:
        r["G11"] = r["Z11"] + r["Z21"]

    # ---- 尾水卜型 ----
    tail = []
    for i in range(M2):
        Fc, Fn, Fb2 = area(DC1[i]), area(DN1[i]), area(DD1[i])
        qc, qn, qb = CQ1[i], ANQ1[i], DQ1[i]
        phi = FF[i]
        # 发电 = 汇流；抽水 = 分流
        tail.append({
            "i": i + 1,
            "Z31": _z(s_merge_straight(qb, qc, Fc, Fb2, Fn), qn, Fn),
            "Z41": _z(s_merge_branch(qb, qc, Fc, Fb2, Fn, phi), qb, Fb2),
            "Z3": _z(s_fork_straight(qn, qc, Fc, Fn), qn, Fn),
            "Z4": _z(s_fork_branch(qb, qc, Fc, Fb2, phi), qb, Fb2),
        })

    # ---- 汇总链（式未定，按同一 S·q²/(2gA²) 口径累加；见 DECL ③） ----
    ZA = max([t["Z31"] + t["Z41"] for t in tail] or [0.0])
    ZB = max([t["Z3"] + t["Z4"] for t in tail] or [0.0])
    gcd = [h["Z1"] + h["Z2"] for h in head_t]
    gcp = [h["Z11"] + h["Z21"] for h in head_p]
    ZMAX = max(gcd + [ZA]) if gcd else ZA
    ZMAX1 = max(gcp + [ZB]) if gcp else ZB
    HT = [v * TQ * TQ for v in (gcd + [ZA])][:4]
    HP = [v * PQ * PQ for v in (gcp + [ZB])][:4]

    return {
        "程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
        "输入": params,
        "引水发电": head_t, "引水抽水": head_p, "尾水": tail,
        "汇总": {"ZA": ZA, "ZB": ZB, "ZMAX": ZMAX, "ZMAX1": ZMAX1,
                 "HT": HT, "HP": HP},
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

_LINE = " " + "*" * 71
_THIN = " " + "-" * 55
_STAR = " " + "*" * 55
_TITLE = " ****              卜型分岔管水头损失计算程序 D-18A                 ****"


def render(params, result):
    """生成原著风格文本计算书（.OUT）。返回以第 2 行（空行）开始的文本。"""
    i = result["输入"]
    S = result["汇总"]
    L = []
    L.append("")
    L.append(_LINE)
    L.append(_TITLE)
    L.append(_LINE)
    L.append("")
    L.append("                   ORIGINAL DATA")
    L.append("                    原 始 数 据")
    L.append(_STAR)
    L.append("        PROJECT NAME  工程名称代号   N=%5d" % i["N"])
    L.append("        引水分岔管个数              M1=%3d" % i["M1"])
    L.append("        尾水分岔管个数              M2=%3d" % i["M2"])
    L.append("        机组台数                    AJ=%7.2f" % i["J"])
    L.append("        发电工况单机流量            TQ=%7.2f (m^3/s)" % i["TQ"])
    L.append("        抽水工况单机流量            PQ=%7.2f (m^3/s)" % i["PQ"])
    L.append("        引水分岔角                  F1=%7.2f (度)" % i["F1"])
    L.append("        引水分岔支管直径            DD=%7.2f (m)" % i["DD"])
    L.append("")
    L.append("")
    L.append("           引 水 分 岔 直 通 管 直 径 (m)")
    L.append(_THIN)
    L.append("          I0          DC            DN")
    L.append("         点号      主管内径      支管内径")
    L.append(_THIN)
    for k in range(i["M1"]):
        L.append("%12d%14.2f%14.2f" % (k + 1, i["DC"][k], i["DN"][k]))
    L.append("")
    L.append("")
    L.append("           尾 水 分 岔 直 通 管 直 径 (m)")
    L.append(_THIN)
    L.append("        I2         DC1          DN1         DD1")
    L.append("       点号     主管内径     支管内径    旁支管内径")
    L.append(_THIN)
    for k in range(i["M2"]):
        L.append("%9d%13.2f%13.2f%13.2f"
                 % (k + 1, i["DC1"][k], i["DN1"][k], i["DD1"][k]))
    L.append("")
    L.append("")
    L.append("              其 它 尾 水 分 岔 管 数 据")
    L.append(_THIN)
    L.append(" I2     CQ1 (台)       ANQ1(台)      DQ1 (台)     FF(度)")
    L.append(" 点    直通管主管     直通管支管     旁支管过      分岔")
    L.append(" 号   过流量机组数   过流量机组数   流量机组数     夹角")
    L.append(_THIN)
    for k in range(i["M2"]):
        L.append("%3d%12.2f%15.2f%13.2f%11.2f"
                 % (k + 1, i["CQ1"][k], i["ANQ1"][k], i["DQ1"][k], i["FF"][k]))
    L.append("")
    L.append("")
    L.append("   RESULT OF P HEADRACE SHAPE MANIFOLD ON TURBINE MODE")
    L.append("           引 水 卜 型 岔 管 发 电 工 况 成 果")
    L.append(_STAR)
    L.append("")
    L.append("    I0          Z1             Z2             G1")
    L.append("   分岔     直通管水头     旁支管水头    本岔管总水头")
    L.append("    号      损失表达式     损失表达式     损失表达式")
    L.append(_THIN)
    for r in result["引水发电"]:
        L.append("%6d%16.8f%15.8f%15.8f" % (r["i"], r["Z1"], r["Z2"], r["G1"]))
    L.append("")
    L.append("")
    L.append("           RESULT OF P.S.T.M.  ON PUMPING MODE")
    L.append("           引 水 卜 型 岔 管 抽 水 工 况 成 果")
    L.append(_STAR)
    L.append("")
    L.append("    I0         Z11            Z21            G11")
    L.append("   分岔     直通管水头     旁支管水头    本岔管总水头")
    L.append("    号      损失表达式     损失表达式     损失表达式")
    L.append(_THIN)
    for r in result["引水抽水"]:
        L.append("%6d%16.8f%15.8f%15.8f" % (r["i"], r["Z11"], r["Z21"], r["G11"]))
    L.append("")
    L.append("")
    L.append("         RESULT OF P. SHAPE TAILRACE MANIFOLD")
    L.append("           ON TURBINE MODE & PUMPING MODE")
    L.append("          尾水卜型岔管发电工况及抽水工况成果")
    L.append(_STAR)
    L.append("")
    L.append("  I0    Z31 (发电)   Z41 (发电)    Z3 (抽水)    Z4 (抽水)")
    L.append(" 分岔   直通管水头   旁支管水头   直通管水头   旁支管水头")
    L.append("  号    损失表达式   损失表达式   损失表达式   损失表达式")
    L.append(" " + "-" * 56)
    for r in result["尾水"]:
        L.append("%4d%14.8f%13.8f%13.8f%13.8f"
                 % (r["i"], r["Z31"], r["Z41"], r["Z3"], r["Z4"]))
    L.append("")
    L.append("")
    L.append("   发电工况本尾水分岔管的总水头损失系数      ZA=%9.3f" % S["ZA"])
    L.append("   抽水工况本尾水分岔管的总水头损失系数      ZB=%9.3f" % S["ZB"])
    L.append("   发电工况本电站岔管最大总水头损失系数    ZMAX=%9.3f" % S["ZMAX"])
    L.append("   抽水工况本电站岔管最大总水头损失系数   ZMAX1=%9.3f" % S["ZMAX1"])
    L.append("")
    L.append("")
    names = ["第一个", "第二个", "第三个", "第四个"]
    for k in range(4):
        L.append("   %s岔管发电工况总水头损失值         HT%d=%7.2f (m)"
                 % (names[k], k + 1, S["HT"][k] if k < len(S["HT"]) else 0.0))
        L.append("   %s岔管抽水工况总水头损失值         HP%d=%7.2f (m)"
                 % (names[k], k + 1, S["HP"][k] if k < len(S["HP"]) else 0.0))
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
