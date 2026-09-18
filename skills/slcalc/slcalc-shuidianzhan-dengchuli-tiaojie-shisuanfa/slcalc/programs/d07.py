# -*- coding: utf-8 -*-
"""
D-7 渠道水力学及水面曲线计算程序 —— 内核
==========================================
复刻《水利水电工程设计计算程序集》D-7 程序（作者：张校正，新疆水利厅）。

功能（原著两种输入模式）
------------------------
  输入 1：计算渠道的正常水深 h0 与临界水深 hk（用于渠道设计）
  输入 2：计算恒定非均匀渐变流水面曲线（棱柱形 / 非棱柱形渠道均可；
          可由上游向下游推算，也可由下游向上游推算）

断面类型（T 组，参数顺序与原著一致）
------------------------------------
  U=1 复式断面：U, n, B, M1, M2, HA, HB, BA, BB, MA, MB, HU
       B   主槽底宽; M1/M2 主槽左/右边坡系数（H:V，0 = 直立）
       HA/HB 左/右岸（滩地）高程（相对主槽底）; 要求 HB<HA，否则左右岸对调
       BA/BB 左/右滩地宽度;  MA/MB 左/右滩地外边坡系数
       HU  主槽带弧形底标志（1=带，0=不带）
  U=2 梯形断面：U, n, B, M
  U=3 圆形断面：U, n, R
  U=4 梯形带圆弧底（M=0 时为 U 形渠槽）：U, n, B, M

核心公式（均经原始 OUT 逐位回归验证，见 d7_verify.py）
-----------------------------------------------------
  过水面积 ω(h)、水面宽 B(h)、湿周 χ(h) 为断面几何函数（见 section_geometry）。
  谢才系数（曼宁）：C = R^(1/6)/n
  流量模数：        K = ω·R^(2/3)/n           R = ω/χ
  正常水深 h0：解  Q = ω·C·√(R·i)  即  f(h) = 1 - Q·n/(ω·R^(2/3)·√i) = 0
                 （即 K(h)·√i - Q = 0；二分法；i≤0 时 h0 无定义，输出 0）
  临界水深 hk：解  αQ²·B/(g·ω³) = 1        （任意断面，h 的隐函数；二分法）
  水面曲线（人工渠槽断面单位能量沿程变化微分方程）：
      dEs/ds = i - Jf
      Es = h + α·v²/(2g)，  v = Q/ω，  Jf = Q²/K²
      差分格式（原著说明书式 17）：
        Es₂ = Es₁ + (i - Jf)·Δs
      二分函数（说明书式 24 解出的 D、G 函数，逐字与 OLE 公式对象核对）：
        F(h₂) = Es₁ - Es₂(h₂) + (i - J̄f)·Δs ,  J̄f = (Jf₁ + Jf₂)/2
      即"E₂ - E₁ = (i - J̄f)·Δs"的隐式离散，二分求 F(h₂)=0。
      原著判据：D、G 同号 → h上限 = h₂；异号 → h下限 = h₂。

曲线类型判据
------------
  i>0： h0 > hk（缓坡 1）/ h0 < hk（陡坡 2）/ h0 = hk（临界坡 3）
        a：h > h0（缓坡）/ h > hk（陡坡）
        b：hk < h < h0（缓坡）/ h0 < h < hk（陡坡）
        c：h < hk（缓坡）/ h < h0（陡坡）
  i≤0：平底 / 倒坡，原著输出代号 "Co"（见 docstring 末尾"未闭合点"）

输出（与原著 .OUT 逐字对齐）
----------------------------
  一.原始数据：流量 / 分段数 / 断面类型数 / 起始断面水深 / 流速分布不均匀系数；
              各断面水力要素；各分段断面类型序号；各分段长度；各分段纵坡
  二.计算结果：分段 / 水深 h / 正常水深 ho / 临界水深 hk / 流速 v / 分段长 dL /
              累计长 L / 曲线类型

未闭合点（详见 d7_verify.py 与 SKILL.md）
-----------------------------------------
  1. 原著二分法"允许误差"未在说明书中给出数值。本内核用严格收敛（1e-10），
     与权威 D-7.OUT 相比 h 最大偏差 0.0053 m（打印精度 0.001 m），
     可解释为原著二分法以 ~0.01 m 为停止阈值（说明书：「直到 |h1-h2| ≤ 允许误差」）。
     ho / hk 同量级（≤0.005 m）。
  2. 复式断面（U=1）的滩地湿周量纲未见于说明书文本（原著为插图）。本内核采用
     水力学通用"分区法"：K = Σ ωᵢ·Rᵢ^(2/3)/n，主槽湿周在主槽岸顶封顶、
     滩地湿周 = 滩宽 + 滩面外坡斜长。对说明书例 2 的三组复式断面，
     hk 误差 ≤0.003 m、h0 误差 ≤0.002 m（sec1/sec3/sec4），
     但 sec2（左岸边坡 M1=1、左滩宽 BA=0）h0 本内核 2.306 对原著 2.330（差 0.024 m），
     由此使例 2 第 1 段曲线类型判为 a2（原著打印 a3，因原著 h0 与 hk 打印值相等）。
     该点为已知不可唯一反演点（缺少对应 .OUT 与断面插图）。
  3. 平底 / 倒坡（i≤0）的水面曲线类型：原著说明书正文写
     「b0(h<hK)、C0(h>hK)」，而权威 OUT 实例中 i=0、h<hk 时打印 "Co"，
     与正文表述相反。本内核按实测 OUT 输出 "Co"。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-7"
TITLE = "恒定非均匀渐变流水面曲线计算书"
AUTHOR = "张校正（新疆水利厅）"
HEAD_NAME = "D-7X"

G = 9.81
DEFAULT_ALPHA = 1.1          # 4.0 版之前的数据文件无 α，原著按 1.1 输出
TOL = 1e-10                  # 内核二分收敛阈值（原著为 ~1e-2，见未闭合点 1）


# ============================================================
# 断面几何与水力要素
# ============================================================

def section_geometry(sec, h):
    """
    断面几何要素。sec 为 dict（见 parse_section）。
    返回 (omega 过水面积, Bw 水面宽, chi 湿周)。
    复式断面（U=1）另返回分块信息 subs（用于分区法流量模数）。
    """
    U = sec["U"]
    if h <= 0:
        return 0.0, 0.0, 0.0

    if U == 2:                                   # 梯形（含矩形 M=0）
        B, M = sec["B"], sec["M"]
        A = B * h + M * h * h
        Bw = B + 2 * M * h
        chi = B + 2 * h * math.sqrt(1 + M * M)
        return A, Bw, chi

    if U == 3:                                   # 圆形
        R = sec["R"]
        if h >= 2 * R:
            h = 2 * R
        t = h - R
        s = math.sqrt(max(R * R - t * t, 0.0))
        A = R * R * math.acos(-t / R) + t * s
        Bw = 2 * s
        chi = 2 * R * math.acos(-t / R)
        return A, Bw, chi

    if U == 4:                                   # 梯形带圆弧底（M=0 → U 形渠槽）
        B, M = sec["B"], sec["M"]
        # 圆弧与边坡相切：R = B / (2(√(1+M²) − M))，M=0 时 R=B/2
        R = B / (2.0 * (math.sqrt(1 + M * M) - M))
        r1m2 = math.sqrt(1 + M * M)
        yt = R * (1 - M / r1m2)                  # 切点高程
        xt = R / r1m2                            # 切点半宽
        if h <= yt:
            # 扇形 − 三角形（弦以下圆缺面积）
            c = (R - h) / R
            A = R * R * math.acos(c) - (R - h) * math.sqrt(max(2 * R * h - h * h, 0.0))
            Bw = 2 * math.sqrt(max(2 * R * h - h * h, 0.0))
            chi = 2 * R * math.acos(c)
            return A, Bw, chi
        c0 = (R - yt) / R
        A0 = R * R * math.acos(c0) - (R - yt) * math.sqrt(max(2 * R * yt - yt * yt, 0.0))
        chi0 = 2 * R * math.acos(c0)
        dh = h - yt
        A = A0 + 2 * xt * dh + M * dh * dh
        Bw = 2 * xt + 2 * M * dh
        chi = chi0 + 2 * dh * r1m2
        return A, Bw, chi

    if U == 1:                                   # 复式断面
        B, M1, M2 = sec["B"], sec["M1"], sec["M2"]
        HA, HB = sec["HA"], sec["HB"]
        BA, BB = sec["BA"], sec["BB"]
        MA, MB = sec["MA"], sec["MB"]
        A = B * h + (M1 + M2) * h * h / 2.0      # 主槽（全水深梯形）
        Bw = B + (M1 + M2) * h
        # 主槽湿周：两岸均封顶于岸顶高程（有滩地的一侧由滩地湿周接管）
        cA = min(h, HA)
        cB = min(h, HB)
        Am = A
        Xm = B + cA * math.sqrt(1 + M1 * M1) + cB * math.sqrt(1 + M2 * M2)
        subs = []
        if BA > 0 and h > HA:
            dh = h - HA
            As = BA * dh + MA * dh * dh / 2.0
            Xs = BA + dh * math.sqrt(1 + MA * MA)
            subs.append((As, Xs))
            A += As
            Bw += BA + MA * dh
        if BB > 0 and h > HB:
            dh = h - HB
            As = BB * dh + MB * dh * dh / 2.0
            Xs = BB + dh * math.sqrt(1 + MB * MB)
            subs.append((As, Xs))
            A += As
            Bw += BB + MB * dh
        chi = Xm + sum(s[1] for s in subs)
        return A, Bw, (Am, Xm, subs)

    raise ValueError("未知断面类型 U=%s" % U)


def omega(sec, h):
    return section_geometry(sec, h)[0]


def conveyance(sec, h, n):
    """流量模数 K。U=1 复式断面按分区法 K = Σ ωᵢ·Rᵢ^(2/3)/n。"""
    if sec["U"] == 1:
        r = section_geometry(sec, h)
        Am, Xm, subs = r[2]
        if Am <= 0 or Xm <= 0:
            return 0.0
        Rm = Am / Xm
        K = Am * Rm ** (2.0 / 3.0) / n
        for As, Xs in subs:
            if As > 0 and Xs > 0:
                K += As * (As / Xs) ** (2.0 / 3.0) / n
        return K
    A, Bw, chi = section_geometry(sec, h)
    if A <= 0 or chi <= 0:
        return 0.0
    return A * (A / chi) ** (2.0 / 3.0) / n


def normal_depth(sec, Q, i):
    """正常水深 h0：解 K(h)·√i = Q。i≤0 返回 0（原著打印 0.000）。"""
    if i <= 0:
        return 0.0
    lo, hi = 1e-9, 1.0
    while conveyance(sec, hi, sec["n"]) * math.sqrt(i) < Q and hi < 1e6:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if conveyance(sec, mid, sec["n"]) * math.sqrt(i) < Q:
            lo = mid
        else:
            hi = mid
        if hi - lo < TOL:
            break
    return 0.5 * (lo + hi)


def critical_depth(sec, Q, alpha):
    """临界水深 hk：解 αQ²·B/(g·ω³) = 1。"""
    def f(h):
        A, Bw, _ = section_geometry(sec, h)[:3]
        if A <= 0:
            return 1.0
        return alpha * Q * Q * Bw / (G * A ** 3) - 1.0
    lo, hi = 1e-9, 1.0
    while f(hi) > 0 and hi < 1e6:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < TOL:
            break
    return 0.5 * (lo + hi)


def curve_type(i, h, h0, hk, tie=0.03):
    """
    水面曲线类型代号（原著打印格式）。
    i ≤ 0：平底 / 倒坡 —— 原著实例打印 "Co"（见未闭合点 3）。
    i > 0：a1/b1/C1/a2/b2/C2/a3/C3（原著对 c 型用大写 C）。
    tie：临界坡度判别容差（|h0-hk| ≤ tie 视为临界坡）。
    """
    if i <= 0:
        return "Co"
    if abs(h0 - hk) <= tie:
        return "a3" if h > h0 else "C3"
    if h0 > hk:                      # 缓坡
        if h > h0:
            return "a1"
        if h > hk:
            return "b1"
        return "C1"
    else:                            # 陡坡
        if h > hk:
            return "a2"
        if h > h0:
            return "b2"
        return "C2"


# ============================================================
# 解析
# ============================================================

def parse_section(nums, pos):
    """按 U 消费断面参数，返回 (sec_dict, new_pos)；数据不足抛 ValueError。"""
    if pos >= len(nums):
        raise ValueError("断面参数不足")
    U = int(round(nums[pos]))
    if U not in (1, 2, 3, 4):
        raise ValueError("未知断面类型 U=%s" % nums[pos])
    pos += 1
    if pos >= len(nums):
        raise ValueError("断面糙率缺失")
    n = nums[pos]
    pos += 1
    if U == 1:
        if pos + 10 > len(nums):
            raise ValueError("U=1 断面参数不足（需 10 个）")
        keys = ["B", "M1", "M2", "HA", "HB", "BA", "BB", "MA", "MB", "HU"]
        vals = nums[pos:pos + 10]
        pos += 10
        sec = {"U": 1, "n": n}
        sec.update(dict(zip(keys, vals)))
        sec["HU"] = int(round(sec["HU"]))
        if sec["HB"] > sec["HA"]:              # 原著：HB<HA，否则左右岸对调
            sec["HA"], sec["HB"] = sec["HB"], sec["HA"]
            sec["BA"], sec["BB"] = sec["BB"], sec["BA"]
            sec["MA"], sec["MB"] = sec["MB"], sec["MA"]
            sec["M1"], sec["M2"] = sec["M2"], sec["M1"]
        return sec, pos
    if U == 2:
        if pos + 2 > len(nums):
            raise ValueError("U=2 断面参数不足（需 2 个）")
        sec = {"U": 2, "n": n, "B": nums[pos], "M": nums[pos + 1]}
        return sec, pos + 2
    if U == 3:
        if pos + 1 > len(nums):
            raise ValueError("U=3 断面参数不足（需 1 个）")
        sec = {"U": 3, "n": n, "R": nums[pos]}
        return sec, pos + 1
    if pos + 2 > len(nums):
        raise ValueError("U=4 断面参数不足（需 2 个）")
    sec = {"U": 4, "n": n, "B": nums[pos], "M": nums[pos + 1]}
    return sec, pos + 2


def _try_mode2(nums, pos, alpha):
    """尝试以给定 alpha 解析模式 2；返回 (结果, 结束位置) 或 None。"""
    try:
        Q = nums[pos]
        K = int(round(nums[pos + 1]))
        T = int(round(nums[pos + 2]))
        H0 = nums[pos + 3]
        pos += 4
        if alpha:
            a = nums[pos]
            pos += 1
        else:
            a = DEFAULT_ALPHA
        sections = []
        for _ in range(T):
            sec, pos = parse_section(nums, pos)
            sections.append(sec)
        if pos + 3 * (K + 1) != len(nums):
            return None
        idx = [int(round(nums[pos + k])) for k in range(K + 1)]
        pos += K + 1
        dL = [nums[pos + k] for k in range(K + 1)]
        pos += K + 1
        slope = [nums[pos + k] for k in range(K + 1)]
        pos += K + 1
    except (IndexError, ValueError):
        return None
    return {"mode": 2, "Q": Q, "K": K, "T": T, "H0": H0, "alpha": a,
            "sections": sections, "idx": idx, "dL": dL, "slope": slope}, pos


def _try_mode1(nums, pos, alpha):
    """尝试以给定 alpha 解析模式 1；返回 (结果, 结束位置) 或 None。"""
    try:
        Q = nums[pos]
        i = nums[pos + 1]
        pos += 2
        if alpha:
            a = nums[pos]
            pos += 1
        else:
            a = DEFAULT_ALPHA
        sec, pos = parse_section(nums, pos)
        if pos != len(nums):
            return None
    except (IndexError, ValueError):
        return None
    return {"mode": 1, "Q": Q, "i": i, "alpha": a, "section": sec}, pos


def parse(data):
    """
    解析输入。data：dict | INT 文件路径。
    返回 dict：
      mode=1：{mode,Q,i,alpha,section}
      mode=2：{mode,Q,K,T,H0,alpha,sections,idx,dL,slope}
    说明：4.0 版起数据文件含 α；旧格式无 α（按原著输出 1.1）。
          本内核以"数值流恰好消费完"为判据自动识别。
    """
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    mode = int(round(nums[0]))
    pos = 1
    tried = []
    for alpha in (True, False):
        r = _try_mode1(nums, pos, alpha) if mode == 1 else _try_mode2(nums, pos, alpha)
        if r is not None:
            tried.append(r)
    if not tried:
        raise ValueError("INT 数据无法解析（模式 %d，数据长度 %d）" % (mode, len(nums)))
    return tried[0][0]


# ============================================================
# 计算
# ============================================================

def _segment_root(sec1, sec2, h1, dL, i, Q, alpha, n):
    """
    二分求渠段末端水深 h2（取与 h1 连续的物理根）。
    F(h) = Es₁ - Es(h) + (i - (Jf₁+Jf₂(h))/2)·Δs
    """
    A1 = omega(sec1, h1)
    Es1 = h1 + alpha * (Q / A1) ** 2 / (2 * G)
    K1 = conveyance(sec1, h1, n)
    J1 = Q * Q / (K1 * K1) if K1 > 0 else 0.0

    def F(h):
        A = omega(sec2, h)
        if A <= 0:
            return -1e18
        Es = h + alpha * (Q / A) ** 2 / (2 * G)
        K = conveyance(sec2, h, n)
        J2 = Q * Q / (K * K) if K > 0 else 0.0
        return Es1 - Es + (i - 0.5 * (J1 + J2)) * dL

    # 粗扫描定位全部变号点，取离 h1 最近的根（物理连续解）
    hmax = min(max(8.0 * max(h1, 1.0), 20.0), 400.0)
    step = max(hmax / 20000.0, 1e-3)
    roots = []
    prev_x, prev_f = 1e-6, F(1e-6)
    x = prev_x + step
    while x <= hmax:
        f = F(x)
        if prev_f * f <= 0:
            a, b, fa = prev_x, x, prev_f
            for _ in range(120):
                m = 0.5 * (a + b)
                fm = F(m)
                if fa * fm <= 0:
                    b = m
                else:
                    a, fa = m, fm
                if b - a < TOL:
                    break
            roots.append(0.5 * (a + b))
        prev_x, prev_f = x, f
        x += step
    if not roots:
        raise ValueError("渠段 F(h)=0 无根（起始水深 %.4f，ΔL=%g，i=%g）" % (h1, dL, i))
    return min(roots, key=lambda r: abs(r - h1))


def compute(params):
    """执行 D-7 计算（模式 1 / 模式 2）。"""
    mode = params.get("mode", 2)
    if mode == 1:
        return _compute_mode1(params)
    return _compute_mode2(params)


def _compute_mode1(params):
    Q = params["Q"]
    i = params["i"]
    alpha = params["alpha"]
    sec = params["section"]
    h0 = normal_depth(sec, Q, i)
    hk = critical_depth(sec, Q, alpha)
    A0 = omega(sec, h0)
    v0 = Q / A0 if A0 > 0 else 0.0
    return {
        "程序": PROGRAM_ID, "模式": 1, "流量Q": Q, "底坡i": i,
        "流速分布不均匀系数": alpha, "断面": _sec_desc(sec),
        "正常水深h0": round(h0, 4), "临界水深hk": round(hk, 4),
        "正常流速v0": round(v0, 4),
        "流态": "急流" if h0 < hk else ("缓流" if h0 > hk else "临界流"),
    }


def _compute_mode2(params):
    Q = params["Q"]
    K = params["K"]
    T = params["T"]
    H0 = params["H0"]
    alpha = params["alpha"]
    sections = params["sections"]
    idx = params["idx"]
    dL = params["dL"]
    slope = params["slope"]

    def sec_at(k):
        return sections[idx[k] - 1]

    rows = []
    h = H0
    L = 0.0
    for k in range(K + 1):
        sec = sec_at(k)
        i = slope[k]
        if k == 0:
            L = dL[0]
        else:
            L = L + dL[k]
        h0 = normal_depth(sec, Q, i)
        hk = critical_depth(sec, Q, alpha)
        A = omega(sec, h)
        v = Q / A if A > 0 else 0.0
        rows.append({
            "k": k, "h": h, "h0": h0, "hk": hk, "v": v,
            "dL": dL[k], "L": L,
            "type": curve_type(i, h, h0, hk),
        })
        if k < K:
            h = _segment_root(sec, sec_at(k + 1), h, dL[k + 1], slope[k + 1],
                              Q, alpha, sec["n"])
    return {
        "程序": PROGRAM_ID, "模式": 2, "流量Q": Q, "分段数K": K,
        "断面类型数T": T, "起始断面水深H0": H0, "流速分布不均匀系数": alpha,
        "断面": [_sec_desc(s) for s in sections],
        "分段断面序号": idx, "分段长": dL, "分段纵坡": slope,
        "行": rows,
    }


def _sec_desc(sec):
    U = sec["U"]
    if U == 1:
        return {"U": 1, "n": sec["n"], "B": sec["B"], "M1": sec["M1"], "M2": sec["M2"],
                "HA": sec["HA"], "HB": sec["HB"], "BA": sec["BA"], "BB": sec["BB"],
                "MA": sec["MA"], "MB": sec["MB"], "HU": sec["HU"]}
    if U == 2:
        return {"U": 2, "n": sec["n"], "B": sec["B"], "M": sec["M"]}
    if U == 3:
        return {"U": 3, "n": sec["n"], "R": sec["R"]}
    return {"U": 4, "n": sec["n"], "B": sec["B"], "M": sec["M"]}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def _str_vb(x):
    """复刻 VB6 Str$()：正数带前导空格；小数不留尾零（如 0.4 → ' .4'）。"""
    neg = x < 0
    ax = abs(float(x))
    if ax == int(ax) and ax < 1e15:
        s = str(int(ax))
    else:
        s = ("%.10f" % ax).rstrip("0").rstrip(".")
    return ("-" if neg else " ") + s


def _sec_line(sec):
    """断面参数回显行（逐字复刻原著，见 SKILL.md 格式说明）。"""
    U = sec["U"]
    line = _str_vb(U) + "  " + ("%.4f" % sec["n"])
    if U == 1:
        rest = [sec["B"], sec["M1"], sec["M2"], sec["HA"], sec["HB"],
                sec["BA"], sec["BB"], sec["MA"], sec["MB"], sec["HU"]]
    elif U == 2:
        rest = [sec["B"], sec["M"]]
    elif U == 3:
        rest = [sec["R"]]
    else:
        rest = [sec["B"], sec["M"]]
    for i, p in enumerate(rest):
        line += ("" if i == 0 else " ") + _str_vb(p)
    return line + " "


def render(params, result):
    """生成原著风格文本计算书（.OUT）。"""
    L = []
    L.append("")
    L.append(" ***********************************************************************")
    L.append(" *****              恒定非均匀渐变流水面曲线计算书 D-7X            *****")
    L.append(" ***********************************************************************")
    L.append("")
    L.append("")
    L.append("")
    if result["模式"] == 1:
        sec = params["section"]
        L.append("          一.原始数据:")
        L.append("")
        L.append(" 流量   渠道纵坡     流速分布不均匀系数")
        L.append("%8.2f%8.4f%16.2f" % (params["Q"], params["i"], params["alpha"]))
        L.append("")
        L.append(" 渠道横断面参数:")
        L.append(_sec_line(sec))
        L.append("")
        L.append("          二.计算结果:")
        L.append("")
        L.append(" 正常水深h0   临界水深hk   正常流速v0   流态")
        L.append("%10.4f%12.4f%12.4f   %s"
                 % (result["正常水深h0"], result["临界水深hk"],
                    result["正常流速v0"], result["流态"]))
        L.append("")
        L.append("")
        return "\n".join(L)

    L.append("          一.原始数据:")
    L.append("")
    L.append(" 流量   分段数     断面类型数   起始断面水深    流速分布不均匀系数")
    L.append("%8.2f%4d%10d%14.2f%16.2f"
             % (params["Q"], params["K"], params["T"], params["H0"], params["alpha"]))
    L.append("")
    L.append(" 断面水力要素:")
    for s in params["sections"]:
        L.append(_sec_line(s))
    L.append("")
    L.append(" 各分段之断面类型:")
    L.append(" " + "  ".join(str(x) for x in params["idx"]) + " ")
    L.append("")
    L.append(" 各分段之长度:")
    L.append(_str_vb(params["dL"][0])
             + "".join(" " + _str_vb(x) for x in params["dL"][1:]) + " ")
    L.append("")
    L.append(" 各分段之纵坡:")
    L.append(" " + " ".join("%.4f" % x for x in params["slope"]) + " ")
    L.append("")
    L.append("")
    L.append("          二.计算结果:")
    L.append("")
    L.append(" 分段    水深h 正常水深ho 临界水深hk  流速v   分段长dL  累计长L 曲线类型")
    for r in result["行"]:
        L.append("%4d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f   %s"
                 % (r["k"], r["h"], r["h0"], r["hk"], r["v"], r["dL"], r["L"], r["type"]))
    L.append("")
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
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(txt)
