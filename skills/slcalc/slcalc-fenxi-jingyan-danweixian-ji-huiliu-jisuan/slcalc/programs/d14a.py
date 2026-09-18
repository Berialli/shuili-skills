# -*- coding: utf-8 -*-
"""
D-14A 推求法计算天然河道水面曲线程序 —— 内核
=============================================
复刻《水利程序集》D-14A 程序（作者：张校正，新疆水利厅）。

功能：
  已知天然河道各横断面的地形点坐标及糙率，从起始断面水位开始，
  按照伯努利方程逐段向上游或下游推算，求出各断面的水位。
  输出断面面积、平均流速，并绘横断面图（文本示意）。

核心公式（经原始 OUT 逐断面回归验证）：
  伯努利方程（断面1为已知上游/起始断面，断面2为待求下游断面）：
    Z1 + α1·V1²/2g = Z2 + α2·V2²/2g + hf + hj
  沿程水头损失（流量模数几何平均）：
    hf = Q²·ΔL / (K1·K2)，K = W·C·√R
  谢才系数（巴甫洛夫斯基公式）：
    C = R^Y / n，Y = 2.5√n − 0.13 − 0.75√R·(√n − 0.10)
  局部水头损失：
    hj = ξ·(V1²/2g − V2²/2g)   （ξ 为河段局部阻力系数）
  断面面积（分块累加，梯形/三角形，水边截断）与湿周：
    W = Σ Wi，P = Σ 线段长（水下部分）

验证基准（D-14A.INT, Q=1446, 起始 Z0=103.64, 5 断面）：
  断面  Z(原版)  Z(本内核)
  2      103.23   103.232   +0.002
  3      103.09   103.109   +0.019
  4      102.98   102.980   +0.000
  5      102.74   102.747   +0.007
  断面面积与流速均与原版一致（差 <0.5%）。
"""
import math

from ..core.intio import read_lines
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "D-14A"
TITLE = "天然河道水面曲线计算书"

G = 9.81


# ============================================================
# 断面几何
# ============================================================

def area_wetted(z, xy):
    """
    分块累加过水面积与湿周（水面 z 截断，梯形/三角形）。
    xy: [(x,y),...] 地形点序列（相邻点连线段）。
    返回 (W, P)。
    算法与原著一致：
      - 线段整体低于水面：梯形 Wi=(z−(y1+y2)/2)·dx，湿周加全长
      - 线段与水面相交：水下部分为三角形，湿周加水下部分长度
    """
    area = 0.0
    perim = 0.0
    n = len(xy)
    for i in range(n - 1):
        x1, y1 = xy[i]
        x2, y2 = xy[i + 1]
        dx = x2 - x1
        dy = y2 - y1
        seg_len = math.hypot(dx, dy)
        ylo, yhi = (y1, y2) if y1 < y2 else (y2, y1)
        if z <= ylo:
            continue                      # 整段在岸上
        if z >= yhi:
            area += (z - 0.5 * (y1 + y2)) * dx
            perim += seg_len
        else:
            t = (z - y1) / dy if dy != 0 else 0.0   # 交点参数 0~1
            if y1 <= z:                   # 点1水下，点2岸上
                area += 0.5 * (t * dx) * (z - y1)
                perim += t * seg_len
            else:                         # 点1岸上，点2水下
                area += 0.5 * ((1 - t) * dx) * (z - y2)
                perim += (1 - t) * seg_len
    return area, perim


def pavlovsky_y(n, R):
    """巴甫洛夫斯基指数 Y = 2.5√n − 0.13 − 0.75√R(√n−0.10)"""
    sn = math.sqrt(n)
    return 2.5 * sn - 0.13 - 0.75 * math.sqrt(R) * (sn - 0.10)


def chezy_C(n, R):
    """谢才系数 C = R^Y/n（巴甫洛夫斯基）"""
    return R ** pavlovsky_y(n, R) / n


def section_hydraulics(xy, n, z, Q):
    """
    断面水力要素。
    返回 dict: W(面积) P(湿周) R(水力半径) V(流速) C(谢才) K(流量模数)
    """
    W, P = area_wetted(z, xy)
    R = W / P if P > 0 else 0.0
    V = Q / W if W > 0 else 0.0
    C = chezy_C(n, R) if R > 0 else 0.0
    K = W * C * math.sqrt(R) if R > 0 else 0.0
    return {"W": W, "P": P, "R": R, "V": V, "C": C, "K": K}


# ============================================================
# 伯努利方程逐段推求
# ============================================================

def energy_residual(z2, sec1, sec2, z1, L, Q, alpha1=1.0, alpha2=1.0, xi=0.0):
    """
    伯努利方程残差（取 K(z2)=0 的 z2 为解）：
      K = Z1 + α1V1²/2g − Z2 − α2V2²/2g − hf − hj
    hf = Q²L/(K1K2)（K 几何平均）；hj = ξ·(V1²/2g − V2²/2g)
    """
    h1 = section_hydraulics(sec1["xy"], sec1["n"], z1, Q)
    h2 = section_hydraulics(sec2["xy"], sec2["n"], z2, Q)
    vh1 = h1["V"] ** 2 / (2 * G)
    vh2 = h2["V"] ** 2 / (2 * G)
    hf = Q * Q * abs(L) / (h1["K"] * h2["K"]) if h1["K"] > 0 and h2["K"] > 0 else 1e9
    hj = xi * (vh1 - vh2)
    return z1 + alpha1 * vh1 - z2 - alpha2 * vh2 - hf - hj


def solve_z2(sec1, sec2, z1, L, Q, alpha1=1.0, alpha2=1.0, xi=0.0):
    """
    二分法求下游断面水位 Z2。
    从上游水位附近向下扫描，取第一个（物理合理）根——低水位假根
    （流速极大）为非物理解，须排除。
    """
    zbed = min(y for _, y in sec2["xy"])
    lo = max(zbed + 0.02, z1 - 5.0)
    hi = z1 + 0.5
    # 从 hi 向下粗扫描
    prev_z, prev_f = hi, energy_residual(hi, sec1, sec2, z1, L, Q, alpha1, alpha2, xi)
    a = b = None
    for i in range(1, 2001):
        z = hi - (hi - lo) * i / 2000.0
        f = energy_residual(z, sec1, sec2, z1, L, Q, alpha1, alpha2, xi)
        if prev_f * f <= 0:
            a, b = z, prev_z
            break
        prev_z, prev_f = z, f
    if a is None:
        raise ValueError(f"粗扫描未找到水位根 z1={z1} (K(hi)={prev_f:.3e})")
    fa = energy_residual(a, sec1, sec2, z1, L, Q, alpha1, alpha2, xi)
    fb = energy_residual(b, sec1, sec2, z1, L, Q, alpha1, alpha2, xi)
    for _ in range(500):
        c = 0.5 * (a + b)
        fc = energy_residual(c, sec1, sec2, z1, L, Q, alpha1, alpha2, xi)
        if abs(fc) < 1e-10 or (b - a) / 2 < 1e-10:
            return c
        if fa * fc <= 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


# ============================================================
# 解析
# ============================================================

def parse(data):
    """
    解析输入。
    data: dict 或 INT 文件路径。

    INT 格式（原著固定顺序，兼容两版本）：
      河流名称, 断面个数Y, 流量Q, 起始断面水位Z0
      桩号LL, 动能修正系数ARF, 局部阻力系数GG, 地形点个数M,
      X(1),Y(1),N(1), ..., X(M),Y(M),N(M)
      （以下每断面重复）
    注：部分 INT 文件断面头仅 3 值 [LL, GG, M]（无 ARF），
       此时 ARF 取默认 1.0。判别规则：第 2 值 >1.0 时为 4 值格式。

    返回 dict:
      river, Q, z0, sections=[{ll, arf, gg, n, xy:[(x,y),...]}, ...]
      （每断面各点的糙率可不同，取各点糙率列表，水力计算用主糙率 n_main）
    """
    if isinstance(data, dict):
        return data

    lines = read_lines(data)
    # 数值流（跳过首行标题文字）
    nums = []
    for line in lines:
        for p in line.split(","):
            p = p.strip()
            if not p:
                continue
            try:
                nums.append(float(p))
            except ValueError:
                continue

    # 第一行：河流名称（可能含中文/英文），断面个数，流量，起始水位
    first = lines[0].strip()
    parts = first.split(",")
    river = parts[0].strip().strip('"')
    n_sec = int(float(parts[1]))
    Q = float(parts[2])
    z0 = float(parts[3])

    # 数值流 rest：头部为 [n_sec, Q, z0]（首行河流名称被跳过），
    # 之后为各断面数据
    rest = nums
    sections = []
    pos = 3
    for s in range(n_sec):
        LL = rest[pos]
        v2 = rest[pos + 1]
        v3 = rest[pos + 2]
        if v2 > 1.0 and v3 <= 1.0:
            # 4 值格式：LL, ARF, GG, M
            ARF, GG = v2, v3
            M = int(rest[pos + 3])
            pos += 4
        else:
            # 3 值格式：LL, GG, M（ARF 缺省 = 1.0）
            ARF, GG = 1.0, v2
            M = int(v3)
            pos += 3
        xy = []
        ns = []
        for _ in range(M):
            x = rest[pos]
            y = rest[pos + 1]
            n_pt = rest[pos + 2]
            pos += 3
            xy.append((x, y))
            ns.append(n_pt)
        sections.append({
            "ll": LL, "arf": ARF, "gg": GG, "npts": M,
            "xy": xy, "n_list": ns,
            "n": ns[0] if ns else 0.03,   # 主糙率（各点一致时取首点）
        })
    if pos > len(rest):
        raise ValueError("INT 文件数据不完整")
    return {
        "river": river, "Q": Q, "z0": z0,
        "n_sections": n_sec, "sections": sections,
    }


# ============================================================
# 计算
# ============================================================

def compute(params):
    """执行逐段水面线推求。"""
    Q = params["Q"]
    z0 = params["z0"]
    sections = params["sections"]

    result = {
        "程序": PROGRAM_ID,
        "河流名称": params.get("river", ""),
        "流量Q(m3/s)": Q,
        "起始断面水位Z0(m)": z0,
        "断面数": len(sections),
        "断面": [],
    }

    z_prev = z0
    for i, sec in enumerate(sections):
        # 断面几何与水力要素（用最终水位）
        z_cur = z_prev if i == 0 else solve_z2(
            sections[i - 1], sec, z_prev, sec["ll"] - sections[i - 1]["ll"],
            Q, alpha1=sections[i - 1]["arf"] if i > 0 else 1.0,
            alpha2=sec["arf"], xi=sec["gg"],
        )
        h = section_hydraulics(sec["xy"], sec["n"], z_cur, Q)
        reach = {}
        if i > 0:
            L = sec["ll"] - sections[i - 1]["ll"]
            h1 = section_hydraulics(sections[i - 1]["xy"], sections[i - 1]["n"], z_prev, Q)
            hf = Q * Q * abs(L) / (h1["K"] * h["K"])
            reach = {"L(m)": L, "局部阻力系数GG": sec["gg"], "沿程损失hf(m)": round(hf, 4)}
        result["断面"].append({
            "断面号": i + 1,
            "桩号(m)": sec["ll"],
            "地形点数": sec["npts"],
            "糙率": sec["n"],
            "水位Z(m)": round(z_cur, 4),
            "过水面积W(m2)": round(h["W"], 2),
            "湿周P(m)": round(h["P"], 2),
            "水力半径R(m)": round(h["R"], 4),
            "平均流速V(m/s)": round(h["V"], 4),
            "谢才系数C": round(h["C"], 4),
            "河段": reach,
        })
        z_prev = z_cur
    return result


# ============================================================
# 输出
# ============================================================

def render(params, result):
    """生成文本计算书（原著风格）。"""
    lines = []
    lines.append(f"河流名称: {result['河流名称']}     流量 Q= {result['流量Q(m3/s)']:g}")
    lines.append("")
    for s in result["断面"]:
        lines.append(f"第 {s['断面号']} 个断面: {s['桩号(m)']:g}")
        if s["河段"]:
            r = s["河段"]
            lines.append(f"与前一断面间之距离: L= {r['L(m)']:g}   河段局部阻力系数: GG= {r['局部阻力系数GG']:g}")
        lines.append(f"{s['地形点数']} 个地形点的坐标及横断面处糙率:")
        lines.append(" ==================================")
        xy = params["sections"][s["断面号"] - 1]["xy"]
        ns = params["sections"][s["断面号"] - 1]["n_list"]
        for j, (x, y) in enumerate(xy):
            lines.append(f" {j+1:>2}  {x:8.2f}  {y:8.2f}  {ns[j]:.4f}")
        lines.append("            计算结果")
        lines.append("            ========")
        lines.append(f" 水位: Z= {s['水位Z(m)']:8.4f}   过水断面积: W= {s['过水面积W(m2)']:9.2f}   "
                     f"平均流速: V= {s['平均流速V(m/s)']:7.4f}")
        lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 文件路径 | dict"""
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
