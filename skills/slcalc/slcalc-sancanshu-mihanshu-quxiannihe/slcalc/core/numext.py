# -*- coding: utf-8 -*-
"""
数值方法库（numext）
====================
复刻《水利程序集》中反复出现的数值方法：
  - 线性相关（最小二乘回归）与相关系数检验
  - 二元三点插值（P-Ⅲ 离均系数 Φ 查表用）
  - 二分法求根
  - 龙格-库塔法（水库调洪数值解）
  - 单纯形加速法（多维无约束寻优，P-Ⅲ 适线参数优选用）
  - 正态分布积分近似公式（概率格纸坐标换算用，原著常数）

所有方法均以纯 Python + math 实现，无第三方硬依赖，
保证内核可单文件分发、离线运行。
"""
import math

# ============================================================
# 1. 直线相关（最小二乘回归）—— A-1 程序核心
# ============================================================

def linear_regression(xs, ys, x_interp=None):
    """
    两系列直线相关计算（A-1 算法，经原始 OUT 回归验证）。

    参数:
      xs, ys   : 等长的两系列观测值（列表/元组）
      x_interp : 需要插补 Y 的 X 值列表（可省略）

    返回:
      dict: 包含 a,b,r,Kr,Er,Kx,Ky,mean_x,mean_y,y_interp
            （与原始 OUT 输出字段一一对应）

    公式（与原始计算书 A-1.OUT 核对一致）：
      b = Sxy / Sxx ; a = my - b*mx
      r = Sxy / sqrt(Sxx*Syy)
      Kr = (1 - r^2) / sqrt(n)          # 相关系数标准差
      Er = 0.6745 * Kr                  # 相关系数机误
      Kx = sqrt(Sxx/(n-1))  Ky = sqrt(Syy/(n-1))   # 样本标准差
    """
    n = len(xs)
    if n < 2 or len(ys) != n:
        raise ValueError("两系列必须等长且不少于 2 个点据")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0 or syy == 0:
        raise ValueError("系列无离差，无法计算相关")

    b = sxy / sxx
    a = my - b * mx
    r = sxy / math.sqrt(sxx * syy)
    kr = (1 - r * r) / math.sqrt(n)
    er = 0.6745 * kr
    kx = math.sqrt(sxx / (n - 1))
    ky = math.sqrt(syy / (n - 1))

    y_interp = None
    if x_interp:
        y_interp = [a + b * x for x in x_interp]

    return {
        "n": n, "a": a, "b": b, "r": r,
        "Kr": kr, "Er": er, "Kx": kx, "Ky": ky,
        "mean_x": mx, "mean_y": my, "y_interp": y_interp,
    }


# ============================================================
# 2. 二元三点插值 —— A-3 P-Ⅲ 离均系数 Φp 查表
# ============================================================

def bilinear_interp(x0, y0, xs, ys, z):
    """
    二元三点插值（双线性插值）。
    xs: x 方向节点（升序）；ys: y 方向节点（升序）；
    z : 二维表 z[i][j]，对应 (xs[i], ys[j])。
    """
    nx, ny = len(xs), len(ys)
    if nx < 2 or ny < 2:
        raise ValueError("插值表至少需要 2x2 节点")

    # 定位 x 区间
    if x0 <= xs[0]:
        i1, i2 = 0, 1
    elif x0 >= xs[-1]:
        i1, i2 = nx - 2, nx - 1
    else:
        i2 = 1
        while i2 < nx - 1 and xs[i2] < x0:
            i2 += 1
        i1 = i2 - 1

    # 定位 y 区间
    if y0 <= ys[0]:
        j1, j2 = 0, 1
    elif y0 >= ys[-1]:
        j1, j2 = ny - 2, ny - 1
    else:
        j2 = 1
        while j2 < ny - 1 and ys[j2] < y0:
            j2 += 1
        j1 = j2 - 1

    x1, x2 = xs[i1], xs[i2]
    y1, y2 = ys[j1], ys[j2]
    dx = x2 - x1
    dy = y2 - y1
    if dx == 0 or dy == 0:
        raise ValueError("插值节点重复")

    z11, z12 = z[i1][j1], z[i1][j2]
    z21, z22 = z[i2][j1], z[i2][j2]

    # 双线性：先在 x 方向，再在 y 方向
    f1 = z11 * (x2 - x0) / dx + z21 * (x0 - x1) / dx
    f2 = z12 * (x2 - x0) / dx + z22 * (x0 - x1) / dx
    return f1 * (y2 - y0) / dy + f2 * (y0 - y1) / dy


# ============================================================
# 3. 二分法求根 —— D-14A 天然河道水面线推求
# ============================================================

def bisect(func, a, b, tol=1e-8, max_iter=200):
    """
    二分法求 f(x)=0 的根（须保证 f(a)*f(b) < 0）。
    """
    fa, fb = func(a), func(b)
    if fa * fb > 0:
        raise ValueError("二分法要求 f(a) 与 f(b) 异号")
    for _ in range(max_iter):
        c = 0.5 * (a + b)
        fc = func(c)
        if abs(fc) < tol or (b - a) / 2 < tol:
            return c
        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return 0.5 * (a + b)


# ============================================================
# 4. 四阶龙格-库塔法 —— C-2 水库调洪数值解
# ============================================================

def rk4_step(f, t, y, dt, *args):
    """单步四阶 Runge-Kutta，f(t,y,*args) 返回 dy/dt。"""
    k1 = f(t, y, *args)
    k2 = f(t + dt / 2, y + dt / 2 * k1, *args)
    k3 = f(t + dt / 2, y + dt / 2 * k2, *args)
    k4 = f(t + dt, y + dt * k3, *args)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def rk4_solve(f, t0, y0, t_end, dt, *args):
    """
    定步长四阶 RK 积分，返回 (t_list, y_list)。
    """
    ts, ys = [t0], [y0]
    t, y = t0, y0
    while t < t_end - 1e-12:
        dt_use = min(dt, t_end - t)
        y = rk4_step(f, t, y, dt_use, *args)
        t += dt_use
        ts.append(t)
        ys.append(y)
    return ts, ys


# ============================================================
# 5. 单纯形加速法（Nelder-Mead）—— A-3 适线参数优选
# ============================================================

def nelder_mead(func, x0, step=0.1, tol=1e-8, max_iter=2000):
    """
    单纯形加速法（Nelder-Mead）多维无约束寻优（标准实现，Numerical Recipes 版）。
    func: 目标函数（取极小）；x0: 初始点。
    返回 (最优解列表, 最优目标值)。
    """
    n = len(x0)
    # 构造初始单纯形
    simplex = [list(x0)]
    for i in range(n):
        pt = list(x0)
        pt[i] += step
        simplex.append(pt)

    fvals = [func(p) for p in simplex]
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5

    for _ in range(max_iter):
        # 排序（最好 → 最差）
        order = sorted(range(n + 1), key=lambda i: fvals[i])
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]

        # 收敛判断：顶点间最大距离
        spread = max(
            max(abs(simplex[i][j] - simplex[0][j]) for j in range(n))
            for i in range(1, n + 1)
        )
        if spread < tol:
            return simplex[0], fvals[0]

        # 质心（去掉最差点）
        centroid = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]

        # 反射
        xr = [centroid[j] + alpha * (centroid[j] - simplex[n][j]) for j in range(n)]
        fr = func(xr)
        if fr < fvals[0]:
            # 扩张
            xe = [centroid[j] + gamma * (xr[j] - centroid[j]) for j in range(n)]
            fe = func(xe)
            if fe < fr:
                simplex[n], fvals[n] = xe, fe
            else:
                simplex[n], fvals[n] = xr, fr
        elif fr < fvals[n - 1]:
            # 反射点优于次差点，接受
            simplex[n], fvals[n] = xr, fr
        else:
            # 收缩
            xc = [centroid[j] + rho * (simplex[n][j] - centroid[j]) for j in range(n)]
            fc = func(xc)
            if fc < fvals[n]:
                simplex[n], fvals[n] = xc, fc
            else:
                # 整体缩小（向最优点收缩）
                for i in range(1, n + 1):
                    for j in range(n):
                        simplex[i][j] = simplex[0][j] + sigma * (simplex[i][j] - simplex[0][j])
                    fvals[i] = func(simplex[i])

    order = sorted(range(n + 1), key=lambda i: fvals[i])
    return simplex[order[0]], fvals[order[0]]


# ============================================================
# 6. 正态分布积分近似公式 —— 频率格纸坐标换算
# ============================================================

# 原著（A-3 说明文档）给出的积分近似公式常数
_B_COEF = [
    1.570796288, 3.706987906e-2, -8.364353589e-4, -2.250947176e-4,
    6.841218299e-6, 5.824235515e-6, -1.04527497e-6, 8.360937017e-8,
    -3.231081277e-9, 3.657763036e-11, 6.936233982e-13,
]


def norm_integral_approx(x):
    """
    正态分布积分曲线近似公式（原著 A-3 常数）。
    用于频率格纸横坐标与 P=50% 原点的相对距离换算。
    """
    y = 1.0 / (1.0 + x)
    poly = sum(c * (y ** i) for i, c in enumerate(_B_COEF))
    return 1.0 - y * poly


def freq_axis_coord(p):
    """
    将频率 p（0~1，如 0.001=0.1%）换算为概率格纸相对距离（以 P=50% 为原点，
    加 3.090 平移至 P=0.1% 处，与原著一致）。
    """
    if p <= 0 or p >= 1:
        raise ValueError("频率 p 必须在 (0,1) 之间")
    # 用分位数近似：u = Φ^{-1}(p)
    # 原著采用误差函数近似反算，这里用 Python 标准库等价实现
    from math import erf, sqrt

    def _inv_norm(p_):
        # Acklam 近似（精度 1e-9，足够工程配线）
        if p_ <= 0.5:
            return -_inv_norm(1.0 - p_)
        a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2,
             1.383577518672690e2, -3.066479806614716e1, 2.506628277459239]
        b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2,
             6.680131188771972e1, -1.328068155288572e1]
        c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838,
             -2.549732539343734, 4.374664141464968, 2.938163982698783]
        d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996,
             3.754408661907416]
        plow = 0.02425
        q = p_ - 0.5
        if abs(q) <= plow:
            r = q * q
            return q * (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) / \
                   (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
        if p_ < 1.0 - plow:
            r = math.sqrt(-2.0 * math.log(1.0 - p_))
            return -(((((c[0] * r + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r + c[5]) / \
                   ((((d[0] * r + d[1]) * r + d[2]) * r + d[3]) * r + 1.0)
        return float("inf")

    u = _inv_norm(p)
    # 原著：l = u + 3.090（原点移至 0.1%），实际换算由绘图层处理
    return u + 3.090


# ============================================================
# 7. 曼宁公式与渠道水力要素 —— D-6 渠道水力计算
# ============================================================

def manning_q(A, R, S, n):
    """曼宁公式：Q = A * R^(2/3) * S^(1/2) / n"""
    return A * (R ** (2.0 / 3.0)) * math.sqrt(S) / n


def trapezoid_elements(b, h, m, n, S):
    """
    梯形断面渠道水力要素。
    b: 底宽；h: 水深；m: 边坡系数（1:m）；n: 糙率；S: 底坡。
    返回 dict: A(面积) P(湿周) R(水力半径) Q(流量) v(流速)
    """
    A = h * (b + m * h)
    P = b + 2.0 * h * math.sqrt(1.0 + m * m)
    R = A / P if P > 0 else 0.0
    Q = manning_q(A, R, S, n) if R > 0 else 0.0
    v = Q / A if A > 0 else 0.0
    return {"A": A, "P": P, "R": R, "Q": Q, "v": v}


# ============================================================
# 8. 伯努利方程水面线推求 —— D-14A 天然河道
# ============================================================

def bernoulli_step(z_up, z_down, v_up, v_down, L, n, R_up, R_down, k_xi=0.0):
    """
    伯努利方程逐段推求（局部损失 + 沿程损失）：
      z_up = z_down + hf + hj + (v_up^2 - v_down^2)/(2g)
      hf = (n^2 * v_mean^2 * L) / (R_mean^(4/3))  （曼宁形式沿程损失）
      hj = k_xi * (v_up - v_down)^2 / (2g)         （局部损失）
    返回能量方程残差（用于二分法求 z_up）。
    """
    g = 9.81
    v_mean = 0.5 * (v_up + v_down)
    R_mean = 0.5 * (R_up + R_down)
    if R_mean <= 0:
        raise ValueError("水力半径必须大于 0")
    hf = (n * n * v_mean * v_mean * L) / (R_mean ** (4.0 / 3.0))
    dv = v_up - v_down
    hj = k_xi * dv * dv / (2.0 * g)
    return z_up - (z_down + hf + hj + (v_up * v_up - v_down * v_down) / (2.0 * g))
