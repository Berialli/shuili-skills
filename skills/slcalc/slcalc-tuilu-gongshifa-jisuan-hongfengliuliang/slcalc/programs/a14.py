# -*- coding: utf-8 -*-
"""
A-14 马斯京根模型最优参数估计程序 —— 内核
===========================================
复刻《水利程序集》A-14 程序（作者：叶泽纲，湖南省水文总站）。

功能：
  已知河段一场洪水的实测入流 I(1..N)、出流 Q(1..N) 与计算时段 T(h)，
  用"槽蓄量排队 + 循环差统计"的新方法直接估计马斯京根模型参数 X、K，
  并按所得参数演算出流过程 O(1..N)。

算法（叶泽纲引文新法；原著公式为 Equation.3 对象，经 MTEF 解码 +
     算例黑盒反推验证，2026-09-03 存档 a14_research*.py）：
  1. 河段槽蓄量 W(1..N)（梯形累计，起点 W(1)=0，时段长 T 小时）：
       W(i) = W(i−1) + T·[ (I(i−1)+I(i))/2 − (Q(i−1)+Q(i))/2 ]
     I/Q 单位 m3/s、T 单位 h → W 单位 m3/s·h（形状用于回归，量纲由 K 吸收）。
  2. W 从小到大排队，I、Q 相应重排（说明书步骤③；循环差估计依赖此序）。
  3. 参数 X（说明书步骤④，eq_01 公式对象解码 + 数值枚举）：
     以重排后序列定义三时刻循环差（下标 k=1..N−2，差基 W 循环）
       V(k) = (Q_k−I_k)(W_{k+1}−W_{k+2}) + (Q_{k+1}−I_{k+1})(W_{k+2}−W_k)
              + (Q_{k+2}−I_{k+2})(W_k−W_{k+1})
       Y(k) = Q_k(W_{k+1}−W_{k+2}) + Q_{k+1}(W_{k+2}−W_k)
              + Q_{k+2}(W_k−W_{k+1})
       X = ΣV(k)Y(k) / ΣV(k)²
     （黑盒枚举：V=Q−I、Y=Q 组合唯一命中 X=0.101237→.101；其余候选
       V=I·Q/I/Q 等误差 >0.005。eq_01 MTEF 视觉结构中 V 项首因子
       究竟是 (Q−I) 差还是含乘号的 I·Q，仍有 worker 解码确认中；
       数值结果以 (Q−I) 为准。）
     输出 X 时取 3 位小数（F6.3），后续 K、演算均用显示 X。
  4. 示蓄流量（说明书步骤⑤，eq_02）：φ(i) = X·I(i) + (1−X)·Q(i)。
  5. 传播时间 K（说明书步骤⑥，eq_03，W 对 φ 的带截距回归斜率）：
       K = [NΣφW − Σφ·ΣW] / [NΣφ² − (Σφ)²]
     （排序与否 LS 不变；数值 K=18.2120，输出 2 位小数 18.21。）
  6. 演算（马斯京根正演，A-8/A-9 同式）：
       den = 2K(1−X)+T
       C0 = (T−2KX)/den，C1 = (T+2KX)/den，C2 = (2K(1−X)−T)/den
       O(1)=Q(1)、O(2)=Q(2)（前两点取实测出流——A-14.OUT 首两行
       演算值=出流值，与 A-8 的 O(1)=I(1) 不同）
       O(i) = C0·I(i) + C1·I(i−1) + C2·O(i−1)，i=3..N
     O 四舍五入取整打印。

验证基准（A-14.INT 算例，原著 A-14.OUT 全文核对）：
  N=14、T=6h；I=(56,66,250,550,595,420,295,210,147,100,74,60,51,46)
  Q=(70,66,102,185,265,335,370,368,310,245,200,165,132,100)
  → X=0.1012373（内部）→ 输出 X= .101；K=18.2120 → 输出 18.21；
  演算 O=(70,66,77,149,276,364,374,344,299,249,201,161,129,105)，
  A-14.OUT 14 行演算值逐位一致（round 全中）。

实现说明（黑盒反推的取舍，标注如下）：
  * W 起点 W(1)=0：说明书注明"Wi 值可自行假定"，回归斜率 K 与 X 的
    循环差结构对 W 平移不变（循环差消常数），故起点任意。
  * 排序后的重排序列仅用于步骤③的 V/Y/X；K 的回归用原始序列或排序
    序列结果相同（LS 对行置换不变），本内核在排序序列上统一计算。
  * VB 内部精度：X 先舍入到 3 位小数（0.101）参与 φ/K/演算，使演算
    第 4 点 148.6→149 命中（全精度 X=0.101237 会给 148.5→148）。
    K 输出 2 位小数。演算结果取整打印。
"""
import os

from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-14"
TITLE = "马斯京根模型最优参数估计计算书"


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def _num(tok):
    tok = tok.strip()
    if not tok:
        raise ValueError("空数值字段")
    low = tok.lower()
    if "." in low or "e" in low:
        return float(tok)
    return int(tok)


def _tokenize(line):
    return [_num(t) for t in line.split(",") if t.strip()]


def parse(data):
    """解析 A-14.INT。data: INT 路径 | 文本 | dict。

    数据文件顺序（原著，A-14.INT 见模块 docstring 示例）：
      14,6,56,70,66,66,250,102,550,185,595,265,420,335,295,370,210
      368,147,310,100,245,74,200,60,165,51,132,46,100
    即：N, T, I(1),Q(1),I(2),Q(2),...,I(N),Q(N)  （I/Q 交替！）
    注意与 A-8 的"先全部 I 再全部 Q"不同。
    """
    if isinstance(data, dict):
        return _parse_dict(data)
    if isinstance(data, str) and os.path.isfile(data):
        with open(data, encoding="utf-8", errors="replace") as f:
            text = f.read()
    else:
        text = str(data)
    text = text.replace("\r", "\n").replace("\x1a", "").strip()
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("A-14 无数据内容")
    nums = [x for ln in lines for x in _tokenize(ln)]
    if len(nums) < 4:
        raise ValueError("A-14 数据不足（需 N, T 与至少一组 I/Q）")

    N = int(nums[0])
    T = float(nums[1])
    if N <= 0:
        raise ValueError(f"A-14 时段数 N={N} 非法")
    need = 2 + 2 * N
    if len(nums) < need:
        raise ValueError(f"A-14 数据不足：需 N,T + {N} 对 I/Q，收到 {len(nums)-2} 个流量值")
    I = [float(nums[2 + 2 * j]) for j in range(N)]
    Q = [float(nums[3 + 2 * j]) for j in range(N)]
    return {"T": T, "N": N, "I": I, "Q": Q}


def _parse_dict(d):
    return {"T": float(d["T"]), "N": int(d["N"]),
            "I": [float(x) for x in d["I"]],
            "Q": [float(x) for x in d["Q"]]}


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _storage_series(I, Q, T):
    """槽蓄量梯形累加，W(1)=0。W(i)=Σ_{k<=i} T·[(I(k−1)+I(k))/2−(Q(k−1)+Q(k))/2]。"""
    W = []
    acc = 0.0
    for i in range(len(I)):
        if i == 0:
            acc += 0.0  # W(1)=0，起点前值不参与
        else:
            acc += T * ((I[i - 1] + I[i]) / 2.0 - (Q[i - 1] + Q[i]) / 2.0)
        W.append(acc)
    return W


def _cyclic_xy(I, Q, W):
    """eq_01：三时刻循环差 V(k)/Y(k)（k=0..N−3，对应 i=1..N−2）。

    差基为 W 的循环差分：对任意序列 A，
      cyc(A; k) = A_k(W_{k+1}−W_{k+2}) + A_{k+1}(W_{k+2}−W_k) + A_{k+2}(W_k−W_{k+1})
    V(k) = cyc(Q−I; k)，Y(k) = cyc(Q; k)。
    返回 X = ΣV(k)Y(k)/ΣV(k)²。
    """
    m = len(I)
    V = []
    Y = []
    for k in range(m - 2):
        d1 = W[k + 1] - W[k + 2]
        d2 = W[k + 2] - W[k]
        d3 = W[k] - W[k + 1]
        V.append((Q[k] - I[k]) * d1 + (Q[k + 1] - I[k + 1]) * d2
                 + (Q[k + 2] - I[k + 2]) * d3)
        Y.append(Q[k] * d1 + Q[k + 1] * d2 + Q[k + 2] * d3)
    num = sum(V[i] * Y[i] for i in range(len(V)))
    den = sum(v * v for v in V)
    if abs(den) < 1e-12:
        raise ValueError("A-14 循环差 V 平方和为 0，无法估计 X")
    return num / den


def _lsq_slope(W, Y):
    """带截距最小二乘：W = K·Y + b → K。K=[NΣYW−ΣYΣW]/[NΣY²−(ΣY)²]。"""
    n = len(W)
    ym = sum(Y) / n
    wm = sum(W) / n
    sxy = sum((Y[i] - ym) * (W[i] - wm) for i in range(n))
    sxx = sum((Y[i] - ym) ** 2 for i in range(n))
    if sxx <= 0:
        raise ValueError("A-14 示蓄流量无变幅，无法回归 K")
    return sxy / sxx


def _muskingum_coeffs(X, K, T):
    """马斯京根演算系数：den=2K(1−X)+T；C0/C1/C2 和为 1。"""
    den = 2.0 * K * (1.0 - X) + T
    if abs(den) < 1e-12:
        raise ValueError("A-14 演算系数分母为 0（K、X、T 组合非法）")
    C0 = (T - 2.0 * K * X) / den
    C1 = (T + 2.0 * K * X) / den
    C2 = (2.0 * K * (1.0 - X) - T) / den
    return C0, C1, C2


def _route(I, Q, X, K, T):
    """马斯京根正演。O(1)=Q(1)、O(2)=Q(2)（前两点取实测出流），
    O(i)=C0·I(i)+C1·I(i−1)+C2·O(i−1)，i=3..N。"""
    C0, C1, C2 = _muskingum_coeffs(X, K, T)
    O = [0.0] * len(I)
    O[0] = Q[0]
    O[1] = Q[1]
    for i in range(2, len(I)):
        O[i] = C0 * I[i] + C1 * I[i - 1] + C2 * O[i - 1]
    return O, (C0, C1, C2)


def compute(params, mode="original"):
    """执行 A-14。返回结构化结果 dict。

    mode：
      "original"（默认）：按原著算法完整估计（W排队→X→φ→K→演算）。
      "theory"：同 original（A-14 无反演悬案，两模式等价，保留参数以
        兼容统一调用约定）。
    """
    T = params["T"]
    N = params["N"]
    I = params["I"]
    Q = params["Q"]
    if len(I) != N or len(Q) != N:
        raise ValueError(f"A-14 I/Q 应各有 N={N} 个值，收到 {len(I)}/{len(Q)}")
    if T <= 0:
        raise ValueError(f"A-14 计算时段 T={T} 非法")
    if N < 3:
        raise ValueError(f"A-14 时段数 N={N} 过少（循环差需 ≥3 点）")

    # ① 槽蓄量
    W_raw = _storage_series(I, Q, T)

    # ② W 排队，I/Q 相应重排（stable argsort，保持等值原序）
    order = sorted(range(N), key=lambda i: W_raw[i])
    Is = [I[i] for i in order]
    Qs = [Q[i] for i in order]
    Ws = [W_raw[i] for i in order]

    # ③ X（eq_01 循环差，用排序后序列）
    X_full = _cyclic_xy(Is, Qs, Ws)
    X = round(X_full, 3)          # VB 输出 F6.3，且后续 φ/K/演算用显示值
    X_disp = X

    # ④ 示蓄流量 φ = X·I+(1−X)·Q（排序序列）
    Phi_s = [X * Is[i] + (1.0 - X) * Qs[i] for i in range(N)]

    # ⑤ K = 带截距回归斜率（排序序列；与原始序列同解）
    K_full = _lsq_slope(Ws, Phi_s)
    K = round(K_full, 2)          # VB 输出 F8.2

    # ⑥ 演算（原始时序，用显示 X/K）
    O_full, (C0, C1, C2) = _route(I, Q, X, K, T)
    O_int = [round(o) for o in O_full]

    return {
        "程序": PROGRAM_ID,
        "模式": "槽蓄量排队+循环差统计估计 X/K 并演算",
        "输入": {
            "计算时段T(h)": T,
            "时段数N": N,
            "实测入流I": [round(x, 2) for x in I],
            "实测出流Q": [round(x, 2) for x in Q],
        },
        "结果": {
            "流量比重因素X": X,
            "传播时间K(h)": K,
            "C0": C0,
            "C1": C1,
            "C2": C2,
            "演算出流O(取整)": O_int,
            "演算出流O(全精度)": [round(o, 4) for o in O_full],
        },
        "中间量": {
            "槽蓄量W": [round(w, 2) for w in W_raw],
            "W排队序": [i + 1 for i in order],
            "排序后入流Is": [round(x, 2) for x in Is],
            "排序后出流Qs": [round(x, 2) for x in Qs],
            "排序后槽蓄Ws": [round(w, 2) for w in Ws],
            "X全精度": X_full,
            "K全精度": K_full,
            "示蓄流量φ": [round(p, 4) for p in Phi_s],
        },
    }


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def render(params, result):
    """生成文本计算书（原著 .OUT 风格）。"""
    inp = result["输入"]
    res = result["结果"]
    mid = result["中间量"]
    lines = []
    lines.append("（一）基本数据")
    lines.append("")
    lines.append(f"实测入流/出流点数 N=   {inp['时段数N']}")
    lines.append(f"计算时距 T=   {inp['计算时段T(h)']:8.2f} 小时")
    lines.append("")
    lines.append("（二）参数估计")
    lines.append("")
    lines.append(f"流量比重因素 X=   {res['流量比重因素X']:.3f}")
    lines.append(f"传播时间 K=   {res['传播时间K(h)']:.2f} 小时")
    lines.append("")
    lines.append("（三）演算结果")
    lines.append("")
    lines.append("  实测入流值    出  流  值    演算出流值")
    for i in range(inp["时段数N"]):
        lines.append(f"{inp['实测入流I'][i]:10.0f}{inp['实测出流Q'][i]:13.0f}"
                     f"{res['演算出流O(取整)'][i]:13.0f}")
    lines.append("")
    lines.append("（四）演算系数")
    lines.append("")
    lines.append(f"C0= {res['C0']:.4f}    C1= {res['C1']:.4f}    C2= {res['C2']:.4f}")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text", mode="original"):
    """统一入口。data: INT 路径 | dict"""
    params = parse(data)
    result = compute(params, mode=mode)
    text = render(params, result)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None,
                   mode="auto" if "--auto" in sys.argv else "original")
    print(txt)
