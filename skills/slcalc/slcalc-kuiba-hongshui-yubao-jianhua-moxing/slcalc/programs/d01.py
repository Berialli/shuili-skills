# -*- coding: utf-8 -*-
"""
D-1 管、洞水力学计算程序 —— 内核
=================================
复刻《水利程序集》D-1 程序（作者：陈靖齐，水电部天津勘测设计院；
数据来源：《水利水电工程设计计算程序集》公之于众版，
乌鲁木齐正海水利科技有限公司，张校正教授级高工技术总负责）。

一、功能（原著「二、功能」）
--------------------------
  (1) 正问题：按几何尺寸（d、l）等参数校核流量 Q；
  (2) 反问题：按设计流量等参数确定管径 d（或方管边长 a）；
  (3) 按流量与尺寸确定水头 ΔH。
  原著只编了前两个问题的程序；本内核亦只做 (1)(2)。

  管洞分类号 T1：  倒虹吸管=1，涵洞=2，取水涵管=3
  取水涵管形式号 T2：分级卧管式=1，转动圆盘闸门式=2

二、公式（来源：原著 D-1Intro.rtf 的 13 个 OLE 公式对象 +
        D-1-1.OUT「三、算法与公式」正文；两者互校后一致）
------------------------------------------------------------------
本次把 D-1Intro.rtf 的 13 个 ``\\objdata`` 解包为 OLE2 → "Equation Native"
→ MTEF v3 记录，逐字节还原字符流（与 D-3/D-4 同一方法），得到原著公式真形：

  obj01  V  = μ·√(2g·ΔH)                                        (1)
  obj02  μ  = 1/√(∑ξ + λ·ℓ/d)                                   (2)
  obj03  C  = (1/n)·R^y                                         (3)
  obj04  y  = 2.5√n − 0.13 − 0.75·√R·(√n − 0.1)                 (4)
  obj05  Q  = V·ω                                               (5)
  obj06  Q  = μ·ω·√(2g·ΔH)                                      (6)
  obj07  Q  = μ·ω·√(2g)·(√H₁ + √H₂)      （取水口，两孔）        (7)
  obj08  ω  = (π/4)·d²                   （圆管；方管 ω=a²）      (8)
  obj09  Q  = ω·C·√(R·i)                 （无压长管·明渠）        (9)
  obj10  C  = (1/n)·R^(1/6)              （无压长管谢才系数）     (10)
  obj11  h″ = h′/2·(√(1 + 8Q²/(g·b²·h′³)) − 1)   （跃后水深）    (11)
  obj12  V小 = N/8 = 9.8·Q·H/8 (m³)      （消力井最小体积）       (12)
  obj13  Q  = μ·ω·√(2g·ΔH)              （同 obj06，有压涵洞）   (13)

另有两条只在正文给出、无公式对象：
  (14) 无压短管（宽顶堰）：Q = m·b·√(2g)·H₀^1.5
  (15) 消力池：Lk = (4.83~5.52)·(h″−h′)；S = 1.25·(h″−h₀)；b₀ = b + 0.4

三、系数与默认值（原著「三、倒虹吸管」「四、取水涵管」正文）
------------------------------------------------------------
  倒虹吸管粗取 ∑ξ：  弯管式 1.8~2.1；缓坡式 1.5；竖井式 2.5
  进口 ξ进：          修圆 0.05~0.10；稍修圆 0.20~0.25；棱角 0.5
  弯管 ξ弯(α)：       20°→0.2  40°→0.3  50°→0.4  60°→0.55
                      70°→0.70 80°→0.90 90°→1.10（线性内插）
  材料摩擦系数 λ：     混凝土 1/45；砌石 1/26
  有压涵洞：λ = 8g/C²，C=(1/n)R^y（谢才系数，指数见 obj04）
  取水口 μ 取 0.62；卧管过水高度 h = 0.40·d（圆卧管）、h = (1/4~1/3)·a（方卧管）；
  涵管明流 i = 1/100~1/200，h = 0.75·d；混凝土圆管 n=0.017，砌石方管 n=0.025。

四、输入数据字段顺序（由 6 个 .INT 算例文件反推，原著未公开字段表）
-----------------------------------------------------------------
  T1=1 倒虹吸管：T1, S, MA, λ, F, ∑ξ, [F=0 时: IT, ξ进, α, ξ弯], L, ΔH, Q, TAG
      S  管形号 1方管 2圆管      MA 材料号 1混凝土 2砌石
      λ=0 时按 MA 取（1/45 或 1/26）
      F  倒虹吸形式号 0不选(细算) 1弯管式 2缓坡式 3竖井式
      ∑ξ=0 且 F≠0 时按 F 粗取；F=0 时按自定义件细算 ξ合计
  T1=2 涵洞：   T1, S, H/d, CM, LS, n, d(a), λ, Lm, ξ进, α, ξ弯, L, ΔH, X14, X15
      CM 涵洞类型 0无压 1有压；LS 管型 0短管 1长管
  T1=3 取水涵管：T1, Q, H₁, H₂, X4, T2, [T2=1: SP, X7, n, i(卧管底坡), …]
                                          [T2=2: X7]
      SP 卧管管形号 1方管 2圆管
  说明：该字段序未被原著字段表确认，系由 6 个算例的数值逐一反推（见 d1_verify.py
        第⑥段反证），属本内核的「非唯一反演点」U1。

五、验证（2026-09-11）
---------------------
  ① 算例 1 逐位对拍权威 D-1-1.OUT（GBK，3360 B）：
       计算管径 0.8021／校核管径 0.8020／μ 0.5943／V 1.8605／Q 0.9398
       内核给出 0.8021／0.8020／0.5943／1.8605／0.9399
     → 5 个显示值 4 EXACT + 1 末位差 1（Q 相对偏差 8.5e-6）。
       见 U2。
  ② 算例 2~6 对拍说明书参考答案（口径差异见 d1_verify.py）。
  ③ 例 6（有压涵洞）内核 Q=392.34，武水 392／清华 392.5（+0.09%）。

六、未闭合点（如实标注）
-----------------------
  U1. .INT 字段顺序：原著未公开字段表，且倒虹吸 F=0 支与 F≠0 支、涵洞有压/无压支、
      取水涵管 T2=1/T2=2 支的字段数不同（9/13/15/9/6 个），只能用算例反推。
      本内核采用「公共头 + 条件支」的读法，6 个算例全部落在合理物理量上
      （见 d1_verify.py 第⑥段反证：另两种候选读法会给出无物理意义的值）。
  U2. 例 1「校核管径的流量 Q」权威值 0.9398，内核 0.9399（末位差 1）。
      已排除的假设（d1_verify.py 第①段穷举）：
        · 取 π=3.14159（原著常数池 0x1488 确为 3.14159）时 Q(0.802)=0.939857；
          取 π=3.1415 才能得 0.9398，但该常数在 EXE 中不存在；
        · 用 d=0.8021/0.8020 各种组合、V/ω 用不同 d 交叉、g=9.81、
          Single 精度逐步模拟 —— 均无法同时命中 V=1.8605 与 Q=0.9398；
        · 唯一自洽解要求校核管径 d∈[0.801949, 0.801992]，而 0.8020 与
          0.8021 的十进制显示都无法区分该区间 → 判定为原著侧一处不可唯一反演的
          取整/单精度边界（与 D-3「7/10 末位差 1」同性质），不作伪造对齐。
  U3. 校核管径取整规则：权威 OUT 显示 0.8020（计算管径 0.8021），
      本内核取「计算管径四舍五入到 0.001 m」。原著说明书另有一版计算书
      （D-1Intro.txt 附录）显示校核管径 0.8000（即取整到 0.1 m），
      两版规则不同 → 版本差异，非唯一反演点。
  U4. 弯管式 ∑ξ 为区间 1.8~2.1（无唯一值）；进口 ξ进 为区间 0.05~0.10、
      0.20~0.25。本内核分别取区间下限（0.05/0.20）与中值 2.0，
      仅 IT=3（棱角 ξ进=0.5）与缓坡式/竖井式（1.5/2.5）为唯一值。
  U5. 例 4 消力池：说明书记「计算 h″ 时公式有错」，武水书值 Lk=5.225、S=0.99
      与本内核按 obj11 正确公式所得不自洽（反证见 d1_verify.py 第④段），
      故按 DECL 处理。
  U6. 例 5 消力井 V：武水书值 Vs=61.2 m³ 要求 V=(9.8QH)/8 中 H≈21.26 m，
      而算例文件给的两个水位参数为 8.5 与 14.4（组合 22.9/5.9/8.5/14.4 均不给
      21.26）→ 消力井体积口径不可唯一反演，按 DECL 处理；进水口 d 用
      H=|H₂−H₁|=5.9 m 得 d²=0.4488（书 0.453，−0.9%）。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "D-1"
TITLE = "管、洞水力学计算书"
AUTHOR = "陈靖齐（水电部天津勘测设计院）"

# ------------------------------------------------------------
# 原著常数
# ------------------------------------------------------------
G = 9.8            # 重力加速度（原著常数池：以 float32(9.8) 加宽为 double 存放）
PI = 3.14159       # 原著圆周率（EXE 常数池 0x1488 为 3.14159，仅 6 位有效数字）
LAM_MA = {1: 1.0 / 45.0, 2: 1.0 / 26.0}      # 材料摩擦系数：混凝土/砌石
XI_COARSE = {1: 2.0, 2: 1.5, 3: 2.5}         # 倒虹吸粗取 ∑ξ（弯管式取 1.8~2.1 中值 2.0）
XI_IN = {1: 0.05, 2: 0.20, 3: 0.5}           # 进口 ξ进（区间的下限；棱角为唯一值 0.5）
BEND_ALPHA = [20, 40, 50, 60, 70, 80, 90]    # 弯管角度（度）
BEND_XI = [0.2, 0.3, 0.4, 0.55, 0.70, 0.90, 1.10]
MU_ORIFICE = 0.62                            # 取水口流量系数（原著四(一)1）
N_CONCRETE = 0.017                           # 混凝土圆管糙率
N_MASONRY = 0.025                            # 砌石方管糙率
RATIO_CIRC = 0.40                            # 圆卧管内过水高度比 h/d
RATIO_SQUARE = 1.0 / 3.0                     # 方卧管 h/a（取 1/4~1/3 的 1/3）
RATIO_SIPHON = 0.75                          # 涵管内过水高度比 h/d（明流）
LK_COEF = 5.52                               # 消力池长度系数 Lk=(4.83~5.52)(h″−h′)，取上限
S_COEF = 1.25                                # 消力池最小深度系数 S=1.25(h″−h₀)
B0_MARGIN = 0.4                              # 消力池宽度加宽 b₀ = b + 0.4
D_ROUND_DISPLAY = 4                          # 计算管径显示/取整位数
D_ROUND_CHECK = 3                            # 校核管径取整位数（工程取值，1 mm）

T1_NAMES = {1: "倒虹吸管", 2: "涵洞", 3: "取水涵管"}
S_NAMES = {1: "方管", 2: "园管"}
MA_NAMES = {1: "混凝土", 2: "砌石"}
F_NAMES = {1: "弯管式", 2: "缓坡式", 3: "竖井式"}
IT_NAMES = {1: "修圆", 2: "稍修圆", 3: "棱角"}
T2_NAMES = {1: "分级卧管式", 2: "转动闸门式"}


# ------------------------------------------------------------
# 数值小工具
# ------------------------------------------------------------

def _interp(xs, ys, x):
    """一维线性插值（端点外沿用首末值）。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


def bend_xi(alpha_deg):
    """弯管水头损失系数 ξ弯(α)（原著表格，线性内插）。"""
    return _interp(BEND_ALPHA, BEND_XI, alpha_deg)


def area(d, s):
    """过水面积 ω：圆管(3.14159/4·d²) 或 方管(a²)。"""
    return PI / 4.0 * d * d if s == 2 else d * d


def mu_pipe(xi, lm, l, d):
    """有压管流速系数 μ = 1/√(∑ξ + λ·L/d)（obj02）。"""
    return 1.0 / math.sqrt(xi + lm * l / d)


def velocity(mu, dh, g=G):
    """有压管流速 V = μ·√(2g·ΔH)（obj01）。"""
    return mu * math.sqrt(2.0 * g * dh)


def discharge_pipe(d, s, dh, xi, lm, l, g=G):
    """有压管流量 Q = μ·ω·√(2g·ΔH)（obj06/obj13）。返回 (Q, μ, V, ω)。"""
    mu = mu_pipe(xi, lm, l, d)
    v = velocity(mu, dh, g)
    om = area(d, s)
    return v * om, mu, v, om


def chezy_y(n, r):
    """谢才系数指数 y = 2.5√n − 0.13 − 0.75·√R·(√n − 0.1)（obj04）。"""
    return 2.5 * math.sqrt(n) - 0.13 - 0.75 * (math.sqrt(r) * (math.sqrt(n) - 0.1))


def lambda_culvert(n, r):
    """有压涵洞 λ = 8g/C²，C = (1/n)·R^y（说明书「对有压涵洞」）。"""
    c = (1.0 / n) * (r ** chezy_y(n, r))
    return 8.0 * G / (c * c), c


def solve_size(q, dh, xi, lm, l, s, g=G, lo=1e-4, hi=20.0, it=200):
    """
    反问题：解 d（或方管边长 a）使 Q = √(2g·ΔH)·ω(d)/√(∑ξ + λ·L/d)。

    Q(d) 在 (0,∞) 上单调递增（ω∝d²，分母随 d 增大而减小），二分求根。
    """
    def f(x):
        return discharge_pipe(x, s, dh, xi, lm, l, g)[0] - q
    flo, fhi = f(lo), f(hi)
    if flo > 0:
        raise ValueError("反问题：下界即已满足流量（ΔH 过大或 Q 过小）")
    if fhi < 0:
        raise ValueError("反问题：上界仍不满足流量（ΔH 过小或 Q 过大）")
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------
# 无压（明渠 / 明流）水力要素
# ------------------------------------------------------------

def circ_partial(y):
    """
    圆管部分充满（水深比 y=h/d）的 (ω/d², R/d) 系数。

      θ = 2·arccos(1 − 2y)；ω = d²/8·(θ − sinθ)；χ = d·θ/2；R = ω/χ
    """
    if y <= 0:
        return 0.0, 0.0
    if y >= 1:
        return PI / 4.0, 0.25
    th = 2.0 * math.acos(1.0 - 2.0 * y)
    om = (th - math.sin(th)) / 8.0
    chi = th / 2.0
    return om, om / chi


def weir_short(q, m, b, g=G):
    """无压短管（宽顶堰）：Q = m·b·√(2g)·H₀^1.5 → 反求 H₀。"""
    if m <= 0 or b <= 0:
        raise ValueError("宽顶堰需要正的流量系数 m 与相当管宽 b")
    h0 = (q / (m * b * math.sqrt(2.0 * g))) ** (2.0 / 3.0)
    return h0, m * b * math.sqrt(2.0 * g) * h0 ** 1.5


def solve_rect_width(q, depth_ratio, n, i, lo=1e-3, hi=20.0, it=200):
    """
    方卧管：矩形断面宽 a、水深 h=k·a，由 Q = ω·C·√(R·i) 反求 a。

      ω = k·a²；χ = a + 2k·a；R = k·a/(1+2k)；C = (1/n)·R^(1/6)
    """
    def f(a):
        h = depth_ratio * a
        om = a * h
        r = om / (a + 2.0 * h)
        c = (1.0 / n) * r ** (1.0 / 6.0)
        return om * c * math.sqrt(r * i) - q
    flo, fhi = f(lo), f(hi)
    if flo > 0:
        raise ValueError("方卧管：下界已满足流量")
    if fhi < 0:
        raise ValueError("方卧管：上界不满足流量（坡度/糙率过大）")
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve_rect_depth(q, b, n, i, lo=1e-4, hi=20.0, it=200):
    """方卧管：给定底宽 b，由明渠公式反求过水深度 h。"""
    def f(h):
        om = b * h
        r = om / (b + 2.0 * h)
        c = (1.0 / n) * r ** (1.0 / 6.0)
        return om * c * math.sqrt(r * i) - q
    flo, fhi = f(lo), f(hi)
    if flo > 0:
        raise ValueError("方卧管：下界已满足流量")
    if fhi < 0:
        raise ValueError("方卧管：上界不满足流量")
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve_circle_partial(q, y, n, i, lo=1e-3, hi=20.0, it=200):
    """圆管明流（水深比 y=h/d）：由 Q = ω·C·√(R·i) 反求管径 d。"""
    ko, kr = circ_partial(y)

    def f(d):
        om = ko * d * d
        r = kr * d
        c = (1.0 / n) * r ** (1.0 / 6.0)
        return om * c * math.sqrt(r * i) - q
    flo, fhi = f(lo), f(hi)
    if flo > 0:
        raise ValueError("圆管明流：下界已满足流量")
    if fhi < 0:
        raise ValueError("圆管明流：上界不满足流量")
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------
# T1 = 1  倒虹吸管（有压流）
# ------------------------------------------------------------

def resolve_xi(p):
    """按输入解析倒虹吸管的水头损失合计 ∑ξ 与来源说明。"""
    if p.get("F", 0) == 0:
        it = int(p.get("IT", 0) or 0)
        x1 = p.get("X1", 0.0) or 0.0
        alpha = p.get("ALPHA", 0.0) or 0.0
        x2 = p.get("X2", 0.0) or 0.0
        src = []
        if x1 <= 0:
            x1 = XI_IN.get(it, 0.0)
            src.append("ξ进=%s(按进口形式%s)" % (x1, IT_NAMES.get(it, "-")))
        else:
            src.append("ξ进=%.4f(输入)" % x1)
        if x2 <= 0 and alpha > 0:
            x2 = bend_xi(alpha)
            src.append("ξ弯=%.4f(按α=%.1f°)" % (x2, alpha))
        else:
            src.append("ξ弯=%.4f(输入)" % x2)
        xi = x1 + x2 + 1.0
        src.append("ξ出=1.0")
        return xi, "0不选(细算)：" + " + ".join(src)
    f = int(p["F"])
    xi_in = p.get("XI", 0.0) or 0.0
    if xi_in > 0:
        return xi_in, "输入 ∑ξ=%.4f" % xi_in
    return XI_COARSE.get(f, 2.0), "粗取 %s ∑ξ=%.2f" % (F_NAMES.get(f, "?"), XI_COARSE.get(f, 2.0))


def run_inverted_siphon(p):
    """倒虹吸管：反求管径（或边长）+ 按取整管径校核流量。"""
    s = int(p["S"])
    ma = int(p.get("MA", 1) or 1)
    lm = p.get("LM", 0.0) or 0.0
    if lm <= 0:
        lm = LAM_MA.get(ma, 1.0 / 45.0)
    xi, xi_src = resolve_xi(p)
    l = float(p["L"])
    dh = float(p["DH"])
    q = float(p["Q"])

    d_raw = solve_size(q, dh, xi, lm, l, s)
    d_calc = round(d_raw, D_ROUND_DISPLAY)
    d_check = round(d_raw, D_ROUND_CHECK)
    q_chk, mu_chk, v_chk, om_chk = discharge_pipe(d_check, s, dh, xi, lm, l)
    q_cal, mu_cal, v_cal, om_cal = discharge_pipe(d_calc, s, dh, xi, lm, l)

    return {
        "T1": 1, "S": s, "MA": ma, "LM": lm, "F": int(p.get("F", 0) or 0),
        "XI": xi, "XI_SRC": xi_src, "L": l, "DH": dh, "Q": q,
        "D_RAW": d_raw, "D_CALC": d_calc, "D_CHECK": d_check,
        "MU_CALC": mu_cal, "V_CALC": v_cal, "OM_CALC": om_cal, "Q_CALC": q_cal,
        "MU": mu_chk, "V": v_chk, "OM": om_chk, "Q_CHECK": q_chk,
    }


# ------------------------------------------------------------
# T1 = 2  涵洞
# ------------------------------------------------------------

def run_culvert(p):
    """
    涵洞：有压（CM=1）按有压管流；无压按明渠（长管 LS=1）或宽顶堰（短管 LS=0）。
    """
    s = int(p["S"])
    cm = int(p.get("CM", 1) or 0)
    ls = int(p.get("LS", 1) or 0)
    n = p.get("N", 0.0) or 0.0
    d = float(p["D"])
    l = float(p["L"])
    dh = float(p["DH"])
    q_in = float(p.get("Q", 0.0) or 0.0)
    x1 = p.get("X1", 0.0) or 0.0
    x2 = p.get("X2", 0.0) or 0.0
    xi = x1 + x2 + 1.0
    lm_in = p.get("LM", 0.0) or 0.0

    out = {"T1": 2, "S": s, "CM": cm, "LS": ls, "N": n, "D": d, "L": l,
           "DH": dh, "XI": xi, "Q_IN": q_in}

    if cm == 1:                                   # 有压涵洞
        r = d / 4.0 if s == 2 else d / 4.0
        if lm_in > 0:
            lm, c = lm_in, None
        else:
            lm, c = lambda_culvert(n or 0.0125, r)
        om = area(d, s)
        out.update(MODE="有压", R=r, C=c, LM=lm, OM=om)
        if d > 0 and dh > 0:                      # 正问题：校核 Q
            q, mu, v, _ = discharge_pipe(d, s, dh, xi, lm, l)
            out.update(MU=mu, V=v, Q=q)
        else:                                     # 反问题：按设计流量求 d
            qdes = p.get("QDES", 0.0) or 0.0
            if qdes <= 0:
                raise ValueError("有压涵洞需给出设计流量（Q 或 QDES）或几何尺寸 d")
            d_calc = solve_culvert_d(qdes, dh, xi, lm, l, s)
            out.update(Q=qdes, D=d_calc, D_CALC=d_calc,
                       D_CALC_R=round(d_calc, D_ROUND_DISPLAY),
                       D_CHECK=round(d_calc, D_ROUND_CHECK))
        return out

    i = p.get("I", 1.0 / 100.0) or 1.0 / 100.0
    if ls == 1:                                   # 无压长管（明渠）
        ko, kr = circ_partial(RATIO_SIPHON) if s == 2 else (RATIO_SIPHON, None)
        if s == 2:
            om = ko * d * d
            r = kr * d
        else:
            h = RATIO_SIPHON * d
            om = d * h
            r = om / (d + 2.0 * h)
        c = (1.0 / (n or 0.025)) * r ** (1.0 / 6.0)
        q = om * c * math.sqrt(r * i)
        out.update(MODE="无压长管", I=i, OM=om, R=r, C=c, Q=q)
    else:                                         # 无压短管（宽顶堰）
        m = p.get("M", 0.0) or 0.0
        b = p.get("B", 0.0) or 0.0
        if m <= 0:
            m = 0.32
        if b <= 0:
            b = d
        h0_in = p.get("H0", 0.0) or 0.0
        if h0_in <= 0:
            h0_in = (p.get("HD", 0.0) or 0.0) * d        # H/d × d = 作用水头 H
        if h0_in <= 0:
            raise ValueError("无压短管（宽顶堰）需要作用水头 H（或 H/d）")
        _, q = weir_short(q_in if q_in > 0 else 1.0, m, b)
        q = m * b * math.sqrt(2.0 * G) * h0_in ** 1.5
        out.update(MODE="无压短管(宽顶堰)", I=i, M=m, B=b, H0=h0_in, Q=q)
    return out


def solve_culvert_d(q, dh, xi, lm, l, s):
    """有压涵洞反问题：解管径 d（与 solve_size 同式）。"""
    return solve_size(q, dh, xi, lm, l, s)


# ------------------------------------------------------------
# T1 = 3  取水涵管
# ------------------------------------------------------------

def run_intake_grade(p):
    """
    分级卧管式取水涵管：
      (1) 取水口 Q = μ√(2g)(√H₁+√H₂)·ω → 反求孔径 d
      (2) 卧管（方管/圆管）明流，h = (1/4~1/3)a 或 0.40d → 反求断面尺寸
      (3) 涵管圆管明流 i、h = 0.75d → 反求管径
      (4) 消力池 Lk=(4.83~5.52)(h″−h′)、S=1.25(h″−h₀)、b₀=b+0.4
    """
    q = float(p["Q"])
    h1 = float(p.get("H1", 0.0) or 0.0)
    h2 = float(p.get("H2", 0.0) or 0.0)
    sp = int(p.get("SP", 1) or 1)                 # 卧管管形号 1方管 2圆管
    n_wg = p.get("N", 0.0) or 0.0
    if n_wg <= 0:
        n_wg = N_MASONRY if sp == 1 else N_CONCRETE
    i_wg = p.get("I", 0.0) or 0.0                 # 卧管底坡
    if i_wg <= 0:
        raise ValueError("分级卧管式需要卧管底坡 i（算例 4 为 1/3）")
    i_sy = p.get("I2", 0.0) or 0.0                # 涵管底坡
    if i_sy <= 0:
        i_sy = 1.0 / 100.0
    n_sy = p.get("N2", 0.0) or 0.0
    if n_sy <= 0:
        n_sy = N_CONCRETE
    mu = p.get("MU", 0.0) or MU_ORIFICE

    # (1) 取水口
    om_or = q / (mu * math.sqrt(2.0 * G) * (math.sqrt(h1) + math.sqrt(h2))) if (h1 + h2) > 0 else 0.0
    d_or = math.sqrt(4.0 * om_or / PI) if om_or > 0 else 0.0

    # (2) 卧管
    if sp == 1:
        a = solve_rect_width(q, RATIO_SQUARE, n_wg, i_wg)
        h_wg = RATIO_SQUARE * a
        r_wg = a * h_wg / (a + 2.0 * h_wg)
        om_wg = a * h_wg
        ratio_txt = "h/a=%.4f" % RATIO_SQUARE
        dim_a, dim_d = a, None
    else:
        d_wg = solve_circle_partial(q, RATIO_CIRC, n_wg, i_wg)
        om_wg, kr = circ_partial(RATIO_CIRC)
        om_wg = om_wg * d_wg * d_wg
        r_wg = kr * d_wg
        h_wg = RATIO_CIRC * d_wg
        ratio_txt = "h/d=%.4f" % RATIO_CIRC
        dim_a, dim_d = None, d_wg
    c_wg = (1.0 / n_wg) * r_wg ** (1.0 / 6.0)
    q_wg = om_wg * c_wg * math.sqrt(r_wg * i_wg)

    # (3) 涵管（圆管明流，h = 0.75d）
    d_sy = solve_circle_partial(q, RATIO_SIPHON, n_sy, i_sy)
    ko, kr2 = circ_partial(RATIO_SIPHON)
    om_sy = ko * d_sy * d_sy
    r_sy = kr2 * d_sy
    c_sy = (1.0 / n_sy) * r_sy ** (1.0 / 6.0)
    h_sy = RATIO_SIPHON * d_sy
    q_sy = om_sy * c_sy * math.sqrt(r_sy * i_sy)

    # (4) 消力池
    hp = h_wg
    b = dim_a if sp == 1 else dim_d
    hpp = hp / 2.0 * (math.sqrt(1.0 + 8.0 * q * q / (G * b * b * hp ** 3)) - 1.0)
    lk = LK_COEF * (hpp - hp)
    s_depth = S_COEF * (hpp - h_sy)
    b0 = b + B0_MARGIN

    return {
        "T1": 3, "T2": 1, "Q": q, "H1": h1, "H2": h2, "SP": sp,
        "N_WG": n_wg, "I_WG": i_wg, "N_SY": n_sy, "I_SY": i_sy, "MU": mu,
        "D_ORIFICE": d_or, "OM_ORIFICE": om_or,
        "WG_A": dim_a, "WG_D": dim_d, "WG_H": h_wg, "WG_R": r_wg,
        "WG_OM": om_wg, "WG_C": c_wg, "WG_Q": q_wg, "WG_RATIO": ratio_txt,
        "SY_D": d_sy, "SY_H": h_sy, "SY_OM": om_sy, "SY_R": r_sy,
        "SY_C": c_sy, "SY_Q": q_sy,
        "HP": hp, "HPP": hpp, "LK": lk, "S": s_depth, "B0": b0, "B": b,
    }


def run_intake_rotary(p):
    """
    转动圆盘闸门式取水涵管：
      进水口 d 由孔口式 Q = μ·ω·√(2g·H) 反求；
      消力井最小体积 V = N/8 = (9.8·Q·H)/8（obj12）。

    水头取法（非唯一反演点 U6）：
      算例文件给出两个水位参数 H1 与 H2B（D-1-5 为 8.5 与 14.4）。
      孔口有效水头取二者之差 |H2B−H1|（D-1-5 → 5.9 m），
      消力井总水头取二者之和 H1+H2B（D-1-5 → 22.9 m）。
      依据：该取法同时最接近武水书值（d²=0.4488 vs 书 0.453，−0.9%；
      V=65.93 vs 书 61.2，+7.7%），其余任何单一水头组合偏差都更大
      （穷举见 d1_verify.py 第⑤段）。
    """
    q = float(p["Q"])
    h1 = float(p.get("H1", 0.0) or 0.0)
    h2 = float(p.get("H2B", p.get("H2", 0.0)) or 0.0)
    mu = p.get("MU", 0.0) or MU_ORIFICE
    heads = [h for h in (h1, h2) if h > 0]
    if not heads:
        raise ValueError("转动闸门式需要有效作用水头 H")
    if len(heads) >= 2:
        h_eff = abs(h2 - h1)
        h_tot = h1 + h2
    else:
        h_eff = h_tot = heads[0]
    om = q / (mu * math.sqrt(2.0 * G * h_eff))
    d_or = math.sqrt(4.0 * om / PI)
    vs = 9.8 * q * h_tot / 8.0
    return {
        "T1": 3, "T2": 2, "Q": q, "H1": h1, "H2": h2, "MU": mu,
        "H_EFF": h_eff, "H_TOT": h_tot, "D_ORIFICE": d_or, "D2_ORIFICE": d_or * d_or,
        "OM_ORIFICE": om, "VS": vs, "N_KW": 9.8 * q * h_tot,
    }


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """
    解析 D-1 输入。data: dict（直通）或 .INT 文件路径。

    .INT 为纯文本数值流，首值为管洞分类号 T1，其后按分支顺序排列
    （字段序由 6 个算例文件反推，详见模块 docstring 四、及未闭合点 U1）。
    """
    if isinstance(data, dict):
        return dict(data)
    nums = read_numbers(data)
    if len(nums) < 2:
        raise ValueError("INT 数据为空或过短")
    t1 = int(round(nums[0]))
    if t1 not in T1_NAMES:
        raise ValueError("未知管洞分类号 T1=%s（应为 1/2/3）" % nums[0])
    p = {"T1": t1, "RAW": nums[1:]}
    v = nums[1:]

    if t1 == 1:                                    # 倒虹吸管
        need = 9
        if len(v) < need:
            raise ValueError("倒虹吸管需要 ≥9 个数值（实得 %d）" % len(v))
        p.update(S=int(round(v[0])), MA=int(round(v[1])), LM=v[2],
                 F=int(round(v[3])), XI=v[4])
        k = 5
        if p["F"] == 0:
            if len(v) < 13:
                raise ValueError("倒虹吸管 F=0（细算）需要 ≥13 个数值（实得 %d）" % len(v))
            p.update(IT=int(round(v[k])), X1=v[k + 1], ALPHA=v[k + 2], X2=v[k + 3])
            k += 4
        p.update(L=v[k], DH=v[k + 1], Q=v[k + 2], TAG=v[k + 3] if len(v) > k + 3 else 0.0)
    elif t1 == 2:                                  # 涵洞
        if len(v) < 14:
            raise ValueError("涵洞需要 ≥14 个数值（实得 %d）" % len(v))
        p.update(S=int(round(v[0])), HD=v[1], CM=int(round(v[2])), LS=int(round(v[3])),
                 N=v[4], D=v[5], LM=v[6], Lm=v[7], X1=v[8],
                 ALPHA=v[9], X2=v[10], L=v[11], DH=v[12],
                 X14=v[13], X15=v[14] if len(v) > 14 else 0.0)
        p["I"] = v[13] if 0 < v[13] < 1 else 1.0 / 100.0
    else:                                          # 取水涵管
        if len(v) < 6:
            raise ValueError("取水涵管需要 ≥6 个数值（实得 %d）" % len(v))
        p.update(Q=v[0], H1=v[1], H2=v[2], X4=v[3], T2=int(round(v[4])))
        rest = v[5:]
        if p["T2"] == 1:
            p.update(SP=int(round(rest[0])) if len(rest) > 0 else 1,
                     X7=rest[1] if len(rest) > 1 else 0.0,
                     N=rest[2] if len(rest) > 2 else 0.0,
                     I=rest[3] if len(rest) > 3 else 0.0)
        else:
            p.update(H2B=rest[0] if len(rest) > 0 else 0.0)
            p["X7"] = p["H2B"]
    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params):
    """执行计算，返回结构化结果 dict。"""
    p = params
    t1 = int(p["T1"])
    res = {"程序": PROGRAM_ID, "管洞分类": T1_NAMES[t1],
           "基本资料": {k: v for k, v in p.items() if k not in ("RAW",)}}
    if t1 == 1:
        c = run_inverted_siphon(p)
        res["计算"] = c
        res["结果"] = {
            "计算管径": c["D_CALC"], "校核管径": c["D_CHECK"],
            "校核流速系数": c["MU"], "校核流速": c["V"], "校核流量": c["Q_CHECK"],
        }
    elif t1 == 2:
        c = run_culvert(p)
        res["计算"] = c
        res["结果"] = {k: c[k] for k in ("MODE", "Q", "OM", "R", "C", "MU", "V")
                       if k in c}
    else:
        t2 = int(p.get("T2", 1))
        c = run_intake_grade(p) if t2 == 1 else run_intake_rotary(p)
        res["计算"] = c
        if t2 == 1:
            res["结果"] = {
                "取水口孔径": c["D_ORIFICE"], "卧管边长": c["WG_A"],
                "卧管管径": c["WG_D"], "卧管水深": c["WG_H"],
                "涵管管径": c["SY_D"], "跃后水深": c["HPP"],
                "消力池长度": c["LK"], "消力池深度": c["S"], "消力池宽度": c["B0"],
            }
        else:
            res["结果"] = {"取水口孔径": c["D_ORIFICE"], "消力井最小体积": c["VS"]}
    return res


# ------------------------------------------------------------
# 输出（复刻原著 D-1.OUT 计算书版式）
# ------------------------------------------------------------

def _f(x, n=4):
    return ("%%.%df" % n) % x


def render(params, result):
    """生成汉字计算书（复刻原著 D-1.OUT 版式与标签）。"""
    p, c = params, result["计算"]
    t1 = int(p["T1"])
    L = []
    L.append("  __________________ 工程、____________ 阶段、_________ 专业、 _________ 部分")
    L.append("")
    L.append("  ==========================================================================")
    L.append("")
    L.append("  一、基本资料与计算假定")
    L.append("")
    if t1 == 1:
        L.append("      ****** 倒虹吸管 ******")
        L.append("  管型：                            %s" % S_NAMES.get(c["S"], "-"))
        L.append("  材料：                            %s" % MA_NAMES.get(c["MA"], "-"))
        L.append("  摩擦系数  Lambda：                     %.4f" % c["LM"])
        L.append("  局部水头损失总和  Xi：                 %.4f" % c["XI"])
        L.append("        （%s）" % c["XI_SRC"])
        L.append("  管长  L ：                            %.4f(米)" % c["L"])
        L.append("  水位差  DH：                           %.4f(米)" % c["DH"])
        L.append("  流量  Q :                              %.4f(米)" % c["Q"])
        L.append("")
        L.append("  二、计算简图")
        L.append("")
        L.append("      见附图一。")
        L.append("")
        L.append(_FORMULA_BLOCK)
        L.append("")
        L.append("  四、计算结果")
        L.append("")
        L.append("  计算管径  D :                          %.4f(米)" % c["D_CALC"])
        L.append("  校核管径  D ：                         %.4f(米)" % c["D_CHECK"])
        L.append("  校核管径相关的流速系数  Mu：           %.4f" % c["MU"])
        L.append("  校核管径相关的流速  V：                %.4f(米/秒)" % c["V"])
        L.append("  校核管径的流量  Q：                    %.4f(立方米/秒)" % c["Q_CHECK"])
        L.append("  输入流量  Q ：                         %.4f(立方米/秒)" % c["Q"])
    elif t1 == 2:
        L.append("      ****** 涵洞 ******")
        L.append("  受力型：                          %s" % ("有压" if c["CM"] == 1 else "无压"))
        L.append("  长度类型：                        %s" % ("长管" if c["LS"] == 1 else "短管"))
        L.append("  管型：                            %s" % S_NAMES.get(c["S"], "-"))
        L.append("  作用水头与管径之比 H/d：          %.4f" % p.get("HD", 0.0))
        L.append("  材料糙率  n：                     %.4f" % c["N"])
        L.append("  园管直径  d ：                    %.4f" % c["D"])
        L.append("  局部水头损失总和  Xi：            %.4f" % c["XI"])
        L.append("  管长  L ：                        %.4f(米)" % c["L"])
        L.append("  水位差  DH：                      %.4f(米)" % c["DH"])
        L.append("")
        L.append("  二、计算结果")
        L.append("")
        L.append("  管流类型：                        %s" % c["MODE"])
        if "LM" in c and c["LM"]:
            L.append("  LM=8*g/C^2                        %.5f" % c["LM"])
        if c.get("C"):
            L.append("  谢才系数  C：                     %.4f" % c["C"])
        L.append("  水力半径  R：                     %.4f(米)" % c["R"])
        L.append("  过水面积  Om：                    %.4f(平方米)" % c["OM"])
        L.append("  计算流量  Q ：                    %.4f(立方米/秒)" % c["Q"])
    else:
        t2 = int(p.get("T2", 1))
        L.append("      ****** 取水涵管 ******")
        L.append("  取水涵管形式：                    %s" % T2_NAMES.get(t2, "-"))
        L.append("  孔1水头  H1：                     %.4f(米)" % c["H1"])
        L.append("  孔2水头  H2：                     %.4f(米)" % c["H2"])
        L.append("  流量  Q ：                        %.4f(立方米/秒)" % c["Q"])
        L.append("")
        L.append("  ******  试算条件及结果  ******")
        L.append("")
        L.append("  取水口流量  Q：                   %.4f(立方米)" % c["Q"])
        L.append("  计算管径  D ：                    %.4f(米)" % c["D_ORIFICE"])
        if t2 == 1:
            L.append("  卧管管形：                        %s" % S_NAMES.get(c["SP"], "-"))
            L.append("  材料糙率  n ：                    %.4f" % c["N_WG"])
            L.append("  卧管底坡  i ：                    %.4f" % c["I_WG"])
            if c["WG_A"]:
                L.append("  方管边长  B :                     %.4f(米)" % c["WG_A"])
            if c["WG_D"]:
                L.append("  园管管径  d ：                    %.4f(米)" % c["WG_D"])
            L.append("  卧管内过水高度  h0：              %.4f(米)" % c["WG_H"])
            L.append("  过水面积  Om：                    %.4f(平方米)" % c["WG_OM"])
            L.append("  水力半径  R：                     %.4f(米)" % c["WG_R"])
            L.append("  谢才系数  C：                     %.4f" % c["WG_C"])
            L.append("  园管管径  d（涵管）：             %.4f(米)" % c["SY_D"])
            L.append("  涵洞底坡  i ：                    %.4f" % c["I_SY"])
            L.append("  涵管水深  h0：                    %.4f(米)" % c["SY_H"])
            L.append("  跃前水深  h'：                    %.4f(米)" % c["HP"])
            L.append("  跃后水深  h''：                   %.4f(米)" % c["HPP"])
            L.append("  消力池长度  LK：                  %.4f(米)" % c["LK"])
            L.append("  消力池最小深度  S：               %.4f(米)" % c["S"])
            L.append("  消力池宽度  b0：                  %.4f(米)" % c["B0"])
        else:
            L.append("  有效作用水头  H ：                %.4f(米)" % c["H_EFF"])
            L.append("  过水面积  Om:                     %.4f(平方米)" % c["OM_ORIFICE"])
            L.append("  水流能量  N ：                    %.4f(千瓦)" % c["N_KW"])
            L.append("  消力井最小体积  vs：              %.4f(立方米)" % c["VS"])
    L.append("")
    L.append("  ======================================================================")
    L.append("")
    L.append("  校核者 ___________  计算者 ______________  计算单位 __________________")
    return render_text(PROGRAM_ID, TITLE, [("", L)])


_FORMULA_BLOCK = """  三、算法与公式

      管,洞水流分为有压流与无压流,有压管又分长管与短管.前者以沿程
  水头损失为主;後者还要考虑局部水头损失.
  1.有压流:
  (1) 流速公式
          V=MU*SQR(2*g*dH)
              式中 MU--流速系数; SQR--开方符号; dH--水位差;
          MU=1/SQR(Xi+LM*L/d)
              式中 Xi--局部水头损失系数总和(进口损失+弯管损失+1.0)
                   LM--材料摩擦系数; L--管长(m); d--管径(m)
      对有压涵洞
          LM=8*g/C^2
              式中  谢才系数  C=1/n*R^Y
              其中  n--材料糙率;  R--水力半径;
                    y=2.5*SQR(n)-0.13-0.75*(SQR(R)*(SQR(n)-0.1)
  (2) 流量公式
          Q=V*OM
              式中 OM--截面积(m^2); (圆管(3.14159/4*d^2) 方管(a^2))
                   V--流速(m/sec);   Q--流量(m^3/sec)
  2.无压长管,按明渠计算
          Q=OM*C*SQR(R*i)
              式中 OM--过水面积(m^2); i--坡度;
          C=1/n*R^(1/6)
  3.无压短管按宽顶堰计算.
          Q=m*b*SQR(2*g)*Ho(3^2)
              式中  m--流量系数;  b--相当管宽(m);  H0=H+V0^2/2/g
  4.倒虹吸管:按有压流计算."""


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 文件路径 | dict。"""
    params = parse(data)
    result = compute(params)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [], result)
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
