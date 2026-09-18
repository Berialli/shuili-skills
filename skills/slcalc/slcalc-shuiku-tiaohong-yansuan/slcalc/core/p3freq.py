# -*- coding: utf-8 -*-
"""
P-Ⅲ（皮尔逊Ⅲ型）频率曲线计算模块
===================================
《水利程序集》A-3 水文频率计算程序的核心数学部分。

P-Ⅲ 密度函数三参数（Qa, Cv, Cs）与伽马分布参数的转换：
  alpha = 4 / Cs^2
  beta  = 2 / (Cv * Cs * Qa)
  delta = Qa * (1 - 2*Cv/Cs)
  x = delta + g / beta,  g ~ Gamma(alpha, 1)

设计值（频率 P，即超过概率 P(X > Qp) = P/100）：
  F = 1 - P/100 = 累积概率 P(X <= Qp)
  g = Gammaincinv(alpha, F)   # 不完全伽马逆函数
  Qp = delta + g / beta

实现完全自包含（级数展开不完全伽马 + 二分求逆），不依赖 scipy，
保证内核可离线单文件分发。scipy 可用时自动优先（精度更高）。
"""
import math

try:
    from scipy.special import gammaincinv as _scipy_gammaincinv
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


# ------------------------------------------------------------
# 不完全伽马函数（自实现，级数展开）
# ------------------------------------------------------------

def _gammainc_lower(a, x, eps=1e-15, max_iter=20000):
    """
    下不完全伽马函数 γ(a,x)（未归一化）。
    γ(a,x) = ∫₀ˣ t^(a-1) e^(-t) dt
    级数展开：γ(a,x) = x^a e^(-x) Σ x^k / (a(a+1)...(a+k))
    """
    if x <= 0:
        return 0.0
    term = 1.0 / a
    total = term
    for k in range(1, max_iter):
        term *= x / (a + k)
        total += term
        # 相对收敛判据（更严格）
        if abs(term) < eps * abs(total):
            break
    return total * math.exp(-x + a * math.log(x))


def _gammainc_upper_series(a, x, eps=1e-15, max_iter=20000):
    """
    上不完全伽马函数 Γ(a,x)（未归一化），级数展开：
    Γ(a,x) = Γ(a) - γ(a,x)
    """
    return math.lgamma(a) - _gammainc_lower(a, x, eps, max_iter)


def gamma_pdf(a, x):
    """伽马分布密度 f(x) = x^(a-1) e^(-x) / Γ(a)"""
    if x < 0:
        return 0.0
    if x == 0:
        return 0.0 if a > 1 else (1.0 / math.gamma(a) if a == 1 else float("inf"))
    return math.exp((a - 1) * math.log(x) - x - math.lgamma(a))


def gammainc_regularized(a, x, eps=1e-15):
    """
    正则化不完全伽马 P(a,x) = γ(a,x)/Γ(a)（下尾概率）。
    """
    if x <= 0:
        return 0.0
    # 级数展开在小 x 稳定，连分数在大 x 稳定；阈值取 a+10 更保守
    if x < a + 10:
        # 级数展开
        return _gammainc_lower(a, x, eps) / math.exp(math.lgamma(a))
    else:
        # 连分数展开（大 x，数值稳定）
        return 1.0 - _gammainc_upper_cf(a, x, eps) / math.exp(math.lgamma(a))


def _gammainc_upper_cf(a, x, eps=1e-15, max_iter=20000):
    """
    上不完全伽马函数 Γ(a,x)（未归一化），连分数展开（Lentz 算法）。
    用于 x 较大时数值稳定。
    """
    # Numerical Recipes 标准实现
    b = x + 1.0 - a
    C = 1.0
    D = 1.0 / b
    f = D
    for i in range(1, max_iter):
        an = -i * (i - a)
        b += 2.0
        D = b + an * D
        if D == 0:
            D = 1e-30
        C = b + an / C
        if C == 0:
            C = 1e-30
        D = 1.0 / D
        delta = C * D
        f *= delta
        if abs(delta - 1.0) < eps:
            break
    return math.exp(-x + a * math.log(x)) * f


def gammaincinv(a, p, tol=1e-14, max_iter=500):
    """
    正则化不完全伽马逆函数：求 x 使 P(a,x) = p（0<p<1）。
    自实现：二分法（不依赖 scipy）。
    注意：极端分位数（p→0 或 p→1）时需高精度，tol 收紧至 1e-14。
    """
    if p <= 0 or p >= 1:
        raise ValueError("p 必须在 (0,1) 之间")
    # 初始上界估计
    lo, hi = 0.0, max(a + 1, 1.0)
    # 倍增找到上界
    while gammainc_regularized(a, hi) < p and hi < 1e8:
        hi *= 2.0
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = gammainc_regularized(a, mid)
        if abs(pm - p) < tol or (hi - lo) < tol * max(1.0, mid):
            return mid
        if pm < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def gammaincinv_std(a, p):
    """逆函数统一入口（scipy 优先）。"""
    if _HAS_SCIPY:
        try:
            from scipy.special import gammaincinv as _sgi
            return float(_sgi(a, p))
        except Exception:
            pass
    return gammaincinv(a, p)


# ------------------------------------------------------------
# P-Ⅲ 频率计算
# ------------------------------------------------------------

def p3_params(Qa, Cv, Cs):
    """三参数 → 伽马分布参数 (alpha, beta, delta)。"""
    if Cv <= 0 or Cs == 0:
        raise ValueError("Cv 必须 > 0，Cs 不能为 0")
    alpha = 4.0 / (Cs * Cs)
    beta = 2.0 / (Cv * Cs * Qa)
    delta = Qa * (1.0 - 2.0 * Cv / Cs)
    return alpha, beta, delta


def p3_quantile(Qa, Cv, Cs, P):
    """
    计算频率 P（%）对应的设计值 Qp 与模比系数 Kp。
    P: 超过概率百分比（0~100，如 0.1 表示 0.1%）。
    """
    alpha, beta, delta = p3_params(Qa, Cv, Cs)
    F = 1.0 - P / 100.0
    if F <= 0 or F >= 1:
        raise ValueError(f"频率 P={P}% 超出可计算范围")
    g = gammaincinv_std(alpha, F)
    Qp = delta + g / beta
    Kp = Qp / Qa if Qa != 0 else 0.0
    return Qp, Kp


def p3_series(Qa, Cv, Cs, freqs):
    """批量计算设计值表。freqs: [P1, P2, ...]（%）。返回 [(P, Kp, Qp), ...]"""
    return [(P, *p3_quantile(Qa, Cv, Cs, P)) for P in freqs]


# ------------------------------------------------------------
# 矩法统计参数（连续系列）
# ------------------------------------------------------------

def moment_params(xs):
    """
    矩法估计均值、Cv、Cs（A-3 原著采用矩法公式）。
      Qa = mean
      Cv = sigma / mean
      Cs = n * Σ(x-Qa)³ / ((n-1)(n-2) * sigma³)   （无偏修正）
    """
    n = len(xs)
    if n < 3:
        raise ValueError("矩法至少需要 3 个数据")
    Qa = sum(xs) / n
    s2 = sum((x - Qa) ** 2 for x in xs) / (n - 1)
    sigma = math.sqrt(s2)
    Cv = sigma / Qa if Qa != 0 else 0.0
    # 偏度（无偏估计）
    m3 = sum((x - Qa) ** 3 for x in xs)
    Cs = (n * m3) / ((n - 1) * (n - 2) * sigma ** 3) if sigma > 0 else 0.0
    return {"Qa": Qa, "Cv": Cv, "Cs": Cs, "sigma": sigma, "n": n}


# ------------------------------------------------------------
# 经验频率（连续/不连续系列）
# ------------------------------------------------------------

def empirical_freqs_continuous(xs, n_years=None):
    """
    连续系列经验频率（数学期望公式 P = m/(n+1)）。
    xs: 从小到大排序后的系列（算法内部会降序处理）。
    返回 [(序号, 值, 频率%), ...]，按值从大到小排列。
    """
    xs_sorted = sorted(xs, reverse=True)
    n = len(xs_sorted)
    out = []
    for i, x in enumerate(xs_sorted, start=1):
        P = i / (n + 1.0) * 100.0
        out.append((i, x, P))
    return out


def empirical_freqs_discontinuous(special_vals, rest_vals, N, n, a, l):
    """
    不连续系列经验频率（含特大值）—— SL44 规范公式，以 A-3X 算例输出验证。
      N: 重现期；n: 实测值个数（含从实测抽出的特大值）；a: 特大值总数；
      l: 实测值中抽出的特大值项数
      special_vals: 全部 a 个特大值（含历史调查与从实测抽出的，乱序）
      rest_vals:    实测中剩余的 n-l 个普通值（已剔除抽出的特大值）
    方法：全部 a+(n-l) 个值合并按从大到小排序，前 a 个（整体第 i 位）用
      P_M = i/(N+1)*100%  (i=1..a)
    其余（整体第 i 位，实测排序 m=i-a）用
      P_m = [a/(N+1) + (1-a/(N+1))*m/(n-l+1)]*100%  (m=1..n-l)
    返回 [(序号, 值, 频率%), ...]。
    验证基准（A-32 算例，N=90,n=26,a=4,l=2）：
      5670/4900/4400/4050 → 1.10/2.20/3.30/4.40%，2860→8.22% ... 67.6→96.20%，
      共 28 点与原著输出吻合（个别点±0.05 为 PC-1500 单精度舍入）。
    """
    if N <= 0:
        raise ValueError("重现期 N 必须 > 0")
    a, l, n = int(a), int(l), int(n)
    n_rest = n - l
    if n_rest <= 0:
        raise ValueError("实测剩余项数 n-l 必须 > 0")

    all_sorted = sorted(list(special_vals) + list(rest_vals), reverse=True)
    out = []
    for i, x in enumerate(all_sorted, start=1):
        if i <= a:
            P = i / (N + 1.0) * 100.0
        else:
            m = i - a  # 实测排序 1..n-l
            P = (a / (N + 1.0) + (1.0 - a / (N + 1.0)) * m / (n_rest + 1.0)) * 100.0
        out.append((i, x, P))
    return out


def moment_params_discontinuous(special_vals, rest_vals, N, n, a, l):
    """
    不连续系列矩法统计参数（修正矩法，SL44 规范）。
      Qa = [Σspecial + (N-a)/(n-l) * Σrest] / N
      Cv² = [Σ(QM-Qa)² + (N-a)/(n-l) * Σ(Qm-Qa)²] / ((N-1) * Qa²)
      Cs = [Σ(QM-Qa)³ + (N-a)/(n-l) * Σ(Qm-Qa)³] / ((N-1) * Qa³ * Cv³)
    验证基准（A-32 算例）：Qa=856.00（理论 855.58），Cv=1.293。
    """
    N = float(N)
    a = int(a)
    l = int(l)
    n = int(n)
    n_rest = n - l
    if N <= 0 or n_rest <= 0:
        raise ValueError("N 与 n-l 必须 > 0")

    factor = (N - a) / n_rest
    sum_special = sum(special_vals)
    sum_rest = sum(rest_vals)
    Qa = (sum_special + factor * sum_rest) / N
    if Qa == 0:
        raise ValueError("均值 Qa 为 0，无法计算")

    s2 = sum((x - Qa) ** 2 for x in special_vals) + factor * sum((x - Qa) ** 2 for x in rest_vals)
    Cv = math.sqrt(s2 / ((N - 1.0) * Qa * Qa))
    sigma = Cv * Qa
    if sigma > 0:
        m3 = sum((x - Qa) ** 3 for x in special_vals) + factor * sum((x - Qa) ** 3 for x in rest_vals)
        Cs = m3 / ((N - 1.0) * sigma ** 3)
    else:
        Cs = 0.0
    return {"Qa": Qa, "Cv": Cv, "Cs": Cs, "sigma": sigma, "n": n, "N": N, "a": a, "l": l}


# ------------------------------------------------------------
# 适线优选目标函数（离差平方和最小）
# ------------------------------------------------------------

def fitting_objective(cv_cs, Qa, emp_series):
    """
    适线优选目标函数：经验点据与理论频率曲线离差平方和。
    cv_cs: [Cv, Cs]；emp_series: [(P%, 实测值), ...]
    返回 S = Σ(Qi - Q̂i)²
    """
    Cv, Cs = cv_cs
    if Cv <= 0 or Cs <= 0:
        return 1e10  # 罚函数
    try:
        S = 0.0
        for P, Qi in emp_series:
            Qp, _ = p3_quantile(Qa, Cv, Cs, P)
            S += (Qi - Qp) ** 2
        return S
    except Exception:
        return 1e10
