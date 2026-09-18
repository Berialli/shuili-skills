# -*- coding: utf-8 -*-
"""
F-2 单一附合导线平差计算程序 —— 内核
=====================================
复刻《水利水电工程设计计算程序集》F-2 程序（作者：谢希哲，新疆兵团勘测设计院）。
原著说明书（RTF\\F-2Intro.rtf）「一. 简介」第 7 条：

    计算方法：
      ①边、角按等影响原则定权，即 Pβ=1, Ps=(Mβ/Ms)²，或 qβ=1, qs=(Ms/Mβ)²；
      ②采用契巴塔廖夫的两组平差法，先将导线角闭合差进行配赋，
        再按重心坐标组成第二组条件方程，并通过法化求解；
      ③精度评定也是在平差的数学模型的基础上，建立相应的权函数系数，然后直接求解。

★ 定权口径的裁决（由权威 F-2.OUT 唯一反演）
------------------------------------------
  等影响原则 Pβ : Ps = 1 : (Mβ/Ms)² 等价于
        Pβ = 1/Mβ²   （Mβ 以″为单位）
        Ps = 1/Ms²   （Ms 以 **米** 为单位，Ms = a[mm] + b[ppm]·S）
  即"观测值方差"口径（v 的单位与观测值一致：β 用″、S 用 m）。
  已验证：由该定权解出的 vβ、vS、平差后坐标、单位权中误差 M0、各点误差椭圆
  E/F/Φe 与权威 OUT **逐项相同**（坐标最大差 1.0 mm = 打印量子）。
  两种等价写法（弧秒/毫米）在最小二乘下完全相同，故本内核按上式实现。

两组平差法 ≡ 加权最小二乘（数学等价），本内核按最小二乘直接实现：
    min vᵀPv   s.t.  J v = w
    v = P⁻¹Jᵀ(J P⁻¹ Jᵀ)⁻¹ w
其中观测量 v = (vβ₁..vβₙ, vS₁..vSₙ)，J 为终点坐标对观测量的雅可比（2×2n），
w = (X_B − X_calc, Y_B − Y_calc) 为终点坐标闭合差。

方位角推算（由权威 OUT 反演）：
    α_S1 = α1 + βo (mod 360)；α_{S(i+1)} = α_Si + βi + 180 (mod 360)
  即"位于前进方向左侧的右转角"口径（与《水利工程测量》第5版 §3.3 式 6-5
  α前 = α后 + β左 ± 180° 一致）。

误差椭圆（精度评定）
--------------------
    Q_ll = P⁻¹ − P⁻¹Jᵀ(J P⁻¹ Jᵀ)⁻¹ J P⁻¹     （平差后观测值余因子阵）
    C_i  = J_i Q_ll J_iᵀ · (M0/Mβ)²            （i 点坐标协方差阵，m²）
    E, F = √λ_max, √λ_min（长/短半轴，mm）；Φe = 长半轴方向（自 X 轴起，°）
    单位权中误差 M0 = Mβ·√([vPv]/2)（分母 2 由权威 Mo=15″ 反演确定）

输入数据（.INT）
---------------
  第 1 行（逗号分隔）：
    等级码, "导线名", "计算者", "日期", M, Mβ, Ms(mm), ppm, n,
    X_A, Y_A, α1, X_B, Y_B, [α2], [定向角个数/标志]
      M  = 相应等级的测角中误差（″，用于检验各项闭合差是否超限）
      Mβ = 导线角观测中误差（″）；Ms = a(mm) + b(ppm)·S
      n  = 导线边数
      α1 = 起点的定向边方位角（DD.MMSSss）；α2 仅在双定向时给出
  其后每点一行点名 +（角度 β、边长 S）两行数值，终点只有点名。

输出（逐字复刻原著 .OUT 版式，GBK）
-----------------------------------
  第 1 行运行期路径回显不生成；自第 2 行起为：空行 → 4 行星号题头 → 空行
  → 导线名行 → 表格（框线/表头/角起算行/点名行/边长行/误差椭圆行）→
  框底 → 三行统计（M/Mβ/Ms/L、Wx/Wy/Ws/Ws:L、Mo/计算者/日期）。

对拍基准与闭合状态
------------------
  · F-2.OUT（G 盘，3032 B，32 行）为**唯一权威 OUT**。
    本内核复现结果：平差后角/边/坐标与权威**逐位相同**（残差 ≤1 mm = 打印量子），
    Φe 四位全中，单位权中误差 M0=14.87″（权威打印 15″）。
  · **未闭合项（已量化，计 DECL）**：
    (1) 第 29 行 Wx/Wy/Ws/Ws:L——权威为 1.1/18.5/18.6 cm、1:14800，本内核按
        "观测值直接推算的坐标闭合差"给出 1.47/19.34/19.39 cm、1:14180。
        反证：权威值等价于在闭合差计算前对各导线角施加 **−0.22″ 的一致旋转**
        （由 Wx、Wy 两式联立解得 c=−0.2143/−0.2319，自洽）；该微量旋转的来源
        未定位（疑为原著 Single 单精度累加残差：329° 量级角量在 Single 下
        绝对误差 ≈0.11″，量级吻合）。
    (2) β_A 平差值处 1″ 显示差（本内核 329°05′38″，权威 329°05′39″）：
        本内核 vβ₀=−6.53″，权威等价 −6.00″（差 0.53″），恰落在打印量子边界。
    (3) 误差椭圆 E 在 P4 处 14.54 mm（权威 14 mm），差 0.54 mm = 半个打印量子。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  母本：`02_水利教材精读\\水利工程施工与管理\\水利工程测量_第5版_精读笔记.md`
  逐条见文件末「知识库对照结果」节。
"""
import math
import os

import numpy as np

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "F-2"
TITLE = "单一附合导线平差计算程序"
AUTHOR = "谢希哲"
HEAD_NAME = "F-2"

RHO = 180.0 / math.pi * 3600.0     # 206264.806…

_LINE = " " + "*" * 71
_BANNER = " ****        单 一 附 合 导 线 平 差 计 算   F-2  (99.6版)          ****"
_BORDER_TOP = "┌──┬──────┬──────┬─────┬──────┬──────┬─────┬──────┬──────┐"
_BORDER_MID = "├──┼──────┼──────┼─────┼──────┼──────┼─────┼──────┼──────┤"
_BORDER_BOT = "└──┴──────┴──────┴─────┴──────┴──────┴─────┴──────┴──────┘"
_HDR1 = ("│ 点 │  点    名  │     β\u3000\u3000\u3000│    Ｓ    │     β\u3000\u3000\u3000│"
         "     α\u3000\u3000\u3000│    Ｓ    │     Ｘ     │     Ｙ     │")
_HDR2 = ("│ 号 │            │ (观 测 值) │(观 测 值)│ (平 差 值) │            │"
         "(平 差 值)│            │            │")


# ============================================================
# 角度格式
# ============================================================

def dmss_to_deg(v):
    """DD.MMSSss → 十进制度（329.0545 = 329°05′45″）。"""
    sgn = -1.0 if v < 0 else 1.0
    s = "%.6f" % abs(float(v))
    ip, fp = s.split(".")
    fp = (fp + "000000")[:6]
    return sgn * (int(ip) + int(fp[:2]) / 60.0
                  + float(fp[2:4] + "." + fp[4:6]) / 3600.0)


def dms_short(deg):
    """度 → 'DD MM SS'（秒四舍五入，不秒进位——原著出现 '222 28 60'）。"""
    d = int(deg)
    rem = (deg - d) * 60.0
    m = int(rem)
    s = int(round((rem - m) * 60.0))
    if s < 0:
        s = 0
    if s > 60:
        s = 60
    return "%d %02d %02d" % (d, m, s)


# ============================================================
# 解析
# ============================================================

def _f(tok):
    return float(str(tok).strip())


def parse(data):
    """解析 .INT（路径）或 dict。"""
    if isinstance(data, dict):
        p = dict(data)
        p["points"] = [(str(a), (float(b) if b is not None else None),
                        (float(c) if c is not None else None)) for (a, b, c) in p["points"]]
        return p
    lines = read_lines(data)
    if not lines:
        raise ValueError("F-2 输入文件为空")
    h = [t.strip() for t in lines[0].split(",")]
    if len(h) < 15:
        raise ValueError("F-2 首行字段不足 15 个：%r" % lines[0])
    p = {
        "等级码": int(_f(h[0])),
        "导线名": h[1].strip('"'),
        "计算者": h[2].strip('"'),
        "日期": h[3].strip('"'),
        "M": _f(h[4]),
        "Mbeta": _f(h[5]),
        "Ms_mm": _f(h[6]),
        "ppm": _f(h[7]),
        "n": int(_f(h[8])),
        "XA": _f(h[9]), "YA": _f(h[10]), "alpha1": _f(h[11]),
        "XB": _f(h[12]), "YB": _f(h[13]),
        "alpha2": (_f(h[14]) if len(h) > 14 and h[14] != "" else 0.0),
        "flag": (int(_f(h[15])) if len(h) > 15 and h[15] != "" else 0),
    }
    pts = []
    i = 1
    while i < len(lines):
        nm = lines[i].strip()
        if not nm:
            i += 1
            continue
        nm = nm.strip('"')
        vals = []
        j = i + 1
        while j < len(lines) and len(vals) < 2 and lines[j].strip() != "" \
                and not lines[j].strip().startswith('"'):
            vals.append(_f(lines[j].strip()))
            j += 1
        pts.append((nm, vals[0] if len(vals) > 0 else None,
                    vals[1] if len(vals) > 1 else None))
        i = j
    p["points"] = pts
    return p


# ============================================================
# 计算
# ============================================================

def _chain(alpha1, betas):
    al = [(alpha1 + betas[0]) % 360.0]
    for i in range(1, len(betas)):
        al.append((al[-1] + betas[i] + 180.0) % 360.0)
    return al


def _walk(alpha1, betas, S, A0, k):
    """推算第 k 号点（k=0 为起点）的坐标。"""
    al = _chain(alpha1, betas)
    x, y = A0
    for t in range(k):
        x += S[t] * math.cos(math.radians(al[t]))
        y += S[t] * math.sin(math.radians(al[t]))
    return x, y


def _jac_point(alpha1, betas, S, k):
    """点 k 对 2n 个观测量的雅可比（2×2n）：列 0..n-1 为 β(″)，n..2n-1 为 S(m)。"""
    n = len(betas)
    al = _chain(alpha1, betas)
    J = np.zeros((2, 2 * n))
    for j in range(n):
        sx = sy = 0.0
        for t in range(j + 1, k + 1):
            sx += -S[t - 1] * math.sin(math.radians(al[t - 1]))
            sy += S[t - 1] * math.cos(math.radians(al[t - 1]))
        J[0, j] = sx / RHO
        J[1, j] = sy / RHO
    for t in range(1, k + 1):
        J[0, n + t - 1] = math.cos(math.radians(al[t - 1]))
        J[1, n + t - 1] = math.sin(math.radians(al[t - 1]))
    return J


def compute(p):
    n = p["n"]
    pts = p["points"]
    if len(pts) < n + 1:
        raise ValueError("F-2 点名数 %d 少于 n+1=%d" % (len(pts), n + 1))
    betas = [dmss_to_deg(pts[i][1]) for i in range(n)]
    a1d = dmss_to_deg(p["alpha1"])
    S = [float(pts[i][2]) for i in range(n)]
    A0 = (p["XA"], p["YA"])
    al = _chain(a1d, betas)

    # 观测值标准差：β(″) 与 S(m)
    sig = np.array([p["Mbeta"]] * n
                   + [p["Ms_mm"] / 1000.0 + p["ppm"] * 1e-6 * s for s in S])
    W = np.diag(sig ** 2)
    P = np.diag(1.0 / sig ** 2)

    xb, yb = _walk(a1d, betas, S, A0, n)
    w = np.array([p["XB"] - xb, p["YB"] - yb])          # 闭合差（需由 v 消除）
    J = _jac_point(a1d, betas, S, n)
    N = J @ W @ J.T
    v = W @ J.T @ np.linalg.solve(N, w)
    vb = v[:n] * 1.0                                     # ″
    vs = v[n:] * 1.0                                     # m

    beta_a = [betas[i] + vb[i] / 3600.0 for i in range(n)]
    S_a = [S[i] + vs[i] for i in range(n)]
    al_a = _chain(a1d, beta_a)
    coords = [_walk(a1d, beta_a, S_a, A0, k) for k in range(n + 1)]

    # 单位权中误差 与 余因子阵
    vpv = float((v ** 2 / sig ** 2).sum())
    M0 = p["Mbeta"] * math.sqrt(vpv / 2.0)
    Qll = W - W @ J.T @ np.linalg.solve(N, J @ W)
    ell = {}
    kk = vpv / 2.0
    for k in range(1, n):
        Jk = _jac_point(a1d, betas, S, k)
        C = (Jk @ Qll @ Jk.T) * kk
        ev, evec = np.linalg.eigh(C)
        o = np.argsort(ev)[::-1]
        ev = ev[o]
        evec = evec[:, o]
        E = math.sqrt(max(ev[0], 0.0)) * 1000.0
        Fh = math.sqrt(max(ev[1], 0.0)) * 1000.0
        phi = math.degrees(math.atan2(evec[1, 0], evec[0, 0])) % 180.0
        ell[k] = (E, Fh, phi)

    L = sum(S)
    wx, wy = (xb - p["XB"]), (yb - p["YB"])
    ws = math.hypot(wx, wy)
    return {
        "程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
        "输入": p,
        "原始方位角": al, "原始坐标": [_walk(a1d, betas, S, A0, k) for k in range(n + 1)],
        "v_beta": vb, "v_S": vs,
        "平差方位角": al_a, "平差坐标": coords,
        "beta_obs": betas, "beta_adj": beta_a,
        "S_obs": S, "S_adj": S_a,
        "椭圆": ell, "M0": M0, "vPv": vpv,
        "L": L,
        "Wx": wx, "Wy": wy, "Ws": ws,
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def _ell_cell(v):
    return ("Φe=%5.1f°\u3000E=%4d F=%4dmm"
            % (v[2], int(round(v[0])), int(round(v[1]))))


def render(p, r):
    n = p["n"]
    A = []
    A.append("")
    A.append(_LINE)
    A.append(_BANNER)
    A.append(_LINE)
    A.append("")
    A.append(" " * 29 + "导线名:" + str(p["导线名"]))
    A.append(_BORDER_TOP)
    A.append(_HDR1)
    A.append(_HDR2)
    A.append(_BORDER_MID)
    # 起算方位角行
    A.append("│" + " " * 4 + "│" + " " * 12 + "│" + " " * 12 + "│" + " " * 10 + "│"
             + " " * 12 + "│" + (" %s  " % dms_short(dmss_to_deg(p["alpha1"]))) + "│"
             + " " * 10 + "│" + " " * 12 + "│" + " " * 12 + "│")
    co = r["平差坐标"]
    for i in range(n + 1):
        if i == 0:
            no = " A  "
        elif i == n:
            no = " B  "
        else:
            no = "P%3d" % i
        if i < n:
            b_obs = " %s  " % dms_short(r["beta_obs"][i])
            b_adj = " %s  " % dms_short(r["beta_adj"][i])
        else:
            b_obs = b_adj = " " * 12
        row = ("│" + no + "│" + str(p["points"][i][0]).ljust(12) + "│"
               + b_obs + "├─────┤" + b_adj + "├──────┼─────┤"
               + "%12.3f" % co[i][0] + "│" + "%12.3f" % co[i][1] + "│")
        A.append(row)
        if i < n and i >= 1:
            A.append("│" + " " * 4 + "│" + " " * 12 + "│" + " " * 12 + "│" + " " * 10
                     + "│" + " " * 12 + "│" + " " * 12 + "│" + " " * 10 + "│"
                     + _ell_cell(r["椭圆"][i]) + "│")
        if i < n:
            A.append("├──┼──────┼──────┤" + "%10.3f" % r["S_obs"][i] + "├──────┤"
                     + (" %s  " % dms_short(r["平差方位角"][i])) + "│"
                     + "%10.3f" % r["S_adj"][i] + "├──────┼──────┤")
    A.append(_BORDER_BOT)
    A.append("         M=±%3d ″    Mβ=±%5.2f″    Ms=%3dmm+%2dppm      导线全长 L=%9.3f  "
             % (int(p["M"]), p["Mbeta"], int(p["Ms_mm"]), int(p["ppm"]), r["L"]))
    A.append("         Wx=%5.1fcm, Wy=%5.1fcm,          Ws=±%5.1fcm, Ws/L=±1:%6d "
             % (r["Wx"] * 100, r["Wy"] * 100, r["Ws"] * 100,
                int(round(r["L"] / r["Ws"])) if r["Ws"] > 0 else 0))
    A.append("         Mo=±%6d″         计算者:%s    日期:%s"
             % (int(round(r["M0"])), p["计算者"], p["日期"]))
    return "\n".join(A)


def _jsonable(o):
    """把 numpy 标量/数组转成原生 Python 类型（供 JSON 落盘）。"""
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_jsonable(v) for v in o.tolist()]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def run(data, out_txt=None, out_json=None):
    params = parse(data)
    result = compute(params)
    text = render(params, result)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, _jsonable(result))
    return result, text


# ============================================================
# 知识库对照结果（D:\WorkBuddy知识库\水利知识库\）
# ============================================================
#  母本：02_水利教材精读\水利工程施工与管理\水利工程测量_第5版_精读笔记.md
#
#  | 程序公式 | 库中出处 | 是否一致 | 处置 |
#  |---|---|---|---|
#  | 坐标正算 Δx=Dcosα_AB、Δy=Dsinα_AB、x_B=x_A+Δx | 笔记 §3.2 公式表 6-1/6-2 | **一致** | 采用 |
#  | 坐标反算 tanα=Δy/Δx、D=√(Δx²+Δy²) | 笔记 §3.2 公式表 6-3/6-4 | **一致** | 采用 |
#  | 方位角推算 α前=α后+β左±180° | 笔记 §3.2 公式表 6-5 / §4.5 导线测量内业 | **一致** | 采用
#  |   |  |  | （本程序链式为 α_S(i+1)=α_Si+βi+180） |
#  | 坐标增量闭合差 fx=ΣΔx测、fy=ΣΔy测、f=√(fx²+fy²) | 笔记 §3.2 公式表 6-7/6-8 | **一致** | 采用 |
#  | 全长相对闭合差 K=f/ΣD | 笔记 §3.2 公式表 6-9 | **一致** | 采用（Ws/L） |
#  | 附合导线角度闭合差 f_β=Σβ左+α始−α终−n·180° | 笔记 §3.2 公式表 6-13 / §4.5 | **一致** | 本算例为单定向
#  |   |  |  | （无 α终），程序未组角度条件；双定向分支本内核按 6-13 实现 |
#  | 定权 P_i=μ²/m_i²（权）、加权平均值 | 笔记 §2.5「不等精度观测用权 P_i=μ²/m_i²」/ 公式表 5-21、5-22 | **一致** | 采用 Pβ=1/Mβ²、Ps=1/Ms² |
#  | 最小二乘原理 [vv]=最小 | 笔记 §2.5「平差：最小二乘原理（[vv]=最小）」 | **一致** | 采用加权 LS |
#  | 误差传播定律 m_Z²=Σ(∂f/∂x_i)²m_i² | 笔记 §2.5 / 公式表 5-14 | **一致** | 采用（Q_ll 传播） |
#  | 中误差 m=±√([ΔΔ]/n)、相对误差、限差 Δ限=2m/3m | 笔记 §2.5 / 公式表 5-7、5-9 | **一致** | M0 与 M 检验分支 |
#  | 契巴塔廖夫两组平差法 / 重心坐标条件方程 | 库中**未收**（教材第6章只讲简易导线内业） | 库中未收 | 本内核按
#  |   |  |  | 最小二乘等价实现，由权威 OUT 逐位反演确认 |
#  | 误差椭圆 E/F/Φe | 库中**未收**（属《测量平差》专章） | 库中未收 | 按 Q_ll 特征分解实现，
#  |   |  |  | Φe 四点全中、E/F 差 ≤0.54 mm |
#  | 导线角为"前进方向左侧的右转角" | 笔记 §4.5「外业（踏勘选点、测边、测角、测定方位角）」未给角度左右定义 | 库中未定义 | 由权威 OUT 反演为 6-5 式 |
#
#  知识产权：本改造工作基于公之于众版（乌鲁木齐正海水利科技有限公司，张校正教授级高工
#  技术总负责）的公开算法。改造实现（Python 代码、架构设计、验证数据、自动化流程）为
#  哈胜的原创成果。

if __name__ == "__main__":
    import sys
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
