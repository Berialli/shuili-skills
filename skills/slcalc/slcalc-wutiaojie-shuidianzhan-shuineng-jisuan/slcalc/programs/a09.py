# -*- coding: utf-8 -*-
"""
A-9 马司京根法分段连续演算程序 —— 内核
========================================
复刻《水利程序集》A-9 程序（作者：郝福良，水电部天津勘测设计院）。

功能：
  已知马斯京根演算参数 X、K、演算时段 T、连续演算河段数 M 及上断面入流
  I(1..N)，分段连续演算，计算下断面出流 O(1..N)。

算法（马斯京根分段连续演算；原著公式为图像，数值结构经算例黑盒反推
     验证，2026-09-03 存档 a9_research*.py）：
  1. 将整段河段分成 M 个连续子河段，各子段传播时间（分段传播时间）：
       KL = K / M
  2. 子段流量比重因素（分段流量比重因素）按下式折算（华东水利学院/
     赵人俊经典分段公式，M=2 时等价于 2X−0.5）：
       XL = 0.5 − (M/2)·(1 − 2X)   （即 XL = 0.5 − 0.5·M·(1−2X)）
     当 M=1 时 XL=X、KL=K（与 A-8 单段演算一致）。
  3. 每个子段按马斯京根演算系数逐时段演算：
       den = 2·KL·(1−XL) + T
       C0 = (T − 2·KL·XL) / den
       C1 = (T + 2·KL·XL) / den
       C2 = (2·KL·(1−XL) − T) / den     （C0+C1+C2=1）
       O(1) = I(1)                       （子段首点沿用其入流首值）
       O(i) = C0·I(i) + C1·I(i−1) + C2·O(i−1)，i=2..N
    前一子段的出流过程即作为后一子段的入流过程（分段连续演算）。
  4. 段间结转：上一子段输出的序列整体传入下一子段。原程序 O 为整型打印
     （F 格式无小数位，等价 FORTRAN ANINT/print 圆整），故在 M>1 且非
     末段时按 "打印取整" 结转——该规则使本算例 35 点中 34 点与 A-9.OUT
     逐位一致（±0），仅第 14 点因源程序内部舍入细节差 1（见"实现说明"）。
  5. 末段输出序列 O(1..N) 打印取整。

验证基准（A-9.INT 算例，原著 A-9.OUT 核对）：
  M=2、X=0.2609、K=8.9095、T=4、N=35，
  XL=0.5−(2/2)(1−2×0.2609)=2×0.2609−0.5=0.0218、
  KL=8.9095/2=4.45475。
  C0=(4−2×4.45475×0.0218)/den≈0.4334、
  C1=(4+2×4.45475×0.0218)/den≈0.5036、
  C2=(2×4.45475×(1−0.0218)−4)/den≈0.0630。
  O(1..35) 与 A-9.OUT 相比：34 点逐位一致、O(14) 差 +1
  （拟合命中 34/35 点 ±0、35/35 点 ±1），见 a09_verify.py。

实现说明（黑盒反推的取舍，标注如下）：
  * 原著公式区为图像（OCR 缺失）。分段折算规则经多轮黑盒反推确认：
    M 等分传播时间 K→KL=K/M（research2 证明 K 不折分时误差 225）；
    分段 X 折算必为 XL=0.5−(M/2)(1−2X)（M=2 得 0.0218，research1/4 在
    ±1 内；赵人俊/华东水利学院分段演算公式，WebSearch 交叉验证）。
  * 段间结转取整规则（逐段把序列打印取整后再进下一段）由反推确定：
    M=2 时"末段 34/35 点 ±0"优于"全精度段间传播"（后者点 13/30 差 1）。
  * 剩余 1 点（O14=1156 vs 理论 1157）为源程序单精度/舍入细节不可完全
    复现（自由 2D 扫描 X'∈[0.018,0.03]、K'∈[K/2±0.15] 均无全命中，
    research13 存档）。允许 ±1，视为 reverse 未完全闭合的悬案。
  * compute 提供 mode="theory"（默认）与 mode="direct"（按给定单段
    X,K 与 M 直接演算，不做分段折算——供校验对比，一般不用）。
"""
import os

from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-9"
TITLE = "马司京根法分段连续演算计算书"


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
    """解析 A-9.INT。data: INT 路径 | 文本 | dict。

    数据文件顺序（原著）：
      M, N, X, K, T
      I(1), ..., I(N)
    各数值以逗号分隔，可一行多值 / 多行；兼容含 ^Z(0x1A) 的旧式文本。
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
        raise ValueError("A-9 无数据内容")
    nums = [x for ln in lines for x in _tokenize(ln)]
    if len(nums) < 5:
        raise ValueError("A-9 数据不足（需 M,N,X,K,T 与至少 1 个入流）")

    M = int(nums[0])
    N = int(nums[1])
    X = float(nums[2])
    K = float(nums[3])
    T = float(nums[4])
    if M <= 0:
        raise ValueError(f"A-9 演算河段数 M={M} 非法（需 ≥1）")
    if N <= 0:
        raise ValueError(f"A-9 时段数 N={N} 非法")
    need = 5 + N
    if len(nums) < need:
        raise ValueError(f"A-9 数据不足：需 M,N,X,K,T + {N} 个 I，"
                         f"收到 {len(nums) - 5} 个流量值")
    I = [float(nums[5 + j]) for j in range(N)]
    return {"M": M, "N": N, "X": X, "K": K, "T": T, "I": I}


def _parse_dict(d):
    return {"M": int(d["M"]), "N": int(d["N"]), "X": float(d["X"]),
            "K": float(d["K"]), "T": float(d["T"]),
            "I": [float(x) for x in d["I"]]}


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _seg_round(v):
    """FORTRAN 打印取整：ANINT(v)=floor(v+0.5)（对正值等价四舍五入，
    避开 Python round 的 banker's rounding）。"""
    import math
    return math.floor(v + 0.5)


def _seg_coeffs(XL, KL, T):
    """子段马斯京根演算系数。den=2KL(1−XL)+T；C0/C1/C2 和为 1。"""
    den = 2.0 * KL * (1.0 - XL) + T
    if abs(den) < 1e-12:
        raise ValueError("A-9 演算系数分母为 0（KL/XL/T 组合非法）")
    C0 = (T - 2.0 * KL * XL) / den
    C1 = (T + 2.0 * KL * XL) / den
    C2 = (2.0 * KL * (1.0 - XL) - T) / den
    return C0, C1, C2


def _route_seg(I, XL, KL, T, keep_round=False):
    """单子段逐时段演算。O(1)=I(1)；
    O(i)=C0·I(i)+C1·I(i−1)+C2·O(i−1)，i=2..N。
    keep_round=True 时每一步都取打印整数值（验证用，一般内核不启用）。
    """
    C0, C1, C2 = _seg_coeffs(XL, KL, T)
    O = [0.0] * len(I)
    O[0] = I[0]
    for i in range(1, len(I)):
        v = C0 * I[i] + C1 * I[i - 1] + C2 * O[i - 1]
        O[i] = _seg_round(v) if keep_round else v
    return O, (C0, C1, C2)


def compute(params, mode="theory"):
    """执行 A-9。返回结构化结果 dict。

    参数：
      M  —— 连续演算河段数（≥1）
      N  —— 入流时段数
      X  —— 流量比重因素（整河段）
      K  —— 传播时间（整河段，小时）
      T  —— 演算时段（小时）
      I  —— 上断面入流过程（长度 N）

    mode：
      "theory"（默认）：分段折算 XL=0.5−(M/2)(1−2X)、KL=K/M 后分段
        连续演算；M=1 时退化单段。
      "direct"：跳过折算，直接以 (X, K) 作为子段参数连演 M 段
        （对照/调试用，非原著流程）。
    """
    M = params["M"]
    N = params["N"]
    X = params["X"]
    K = params["K"]
    T = params["T"]
    I = params["I"]
    if len(I) != N:
        raise ValueError(f"A-9 I 应有 N={N} 个值，收到 {len(I)}")
    if T <= 0:
        raise ValueError(f"A-9 演算时段 T={T} 非法")
    if K <= 0:
        raise ValueError(f"A-9 传播时间 K={K} 非法")

    # 分段折算
    if mode == "theory":
        # XL = 0.5 − (M/2)(1 − 2X)；M=1 时退化为 X
        XL = 0.5 - (M / 2.0) * (1.0 - 2.0 * X)
        KL = K / float(M)
    else:  # direct（调试对照）
        XL = X
        KL = K

    # 逐子段连续演算：上段出流作为下段入流；M>1 时前段序列打印取整后结转
    cur = list(I)
    C_all = []
    for seg in range(M):
        O, Cs = _route_seg(cur, XL, KL, T)
        C_all.append(Cs)
        if seg < M - 1:
            cur = [_seg_round(v) for v in O]   # 段间取整结转
        else:
            cur = O                            # 末段保留全精度待打印

    O_full = cur
    O_int = [_seg_round(v) for v in O_full]

    den = 2.0 * KL * (1.0 - XL) + T
    C0, C1, C2 = C_all[-1]

    # 输出段间结转的级联表（逐段首点与末点，诊断用）
    return {
        "程序": PROGRAM_ID,
        "模式": ("M 段连续演算（分段折算）" if mode == "theory"
                 else "M 段连续演算（直接参数，对照）"),
        "输入": {
            "演算河段数M": M,
            "时段数N": N,
            "流量比重因素X": X,
            "传播时间K(h)": K,
            "演算时段T(h)": T,
            "上断面入流I": [round(x, 2) for x in I],
        },
        "结果": {
            "分段流量比重因素XL": XL,
            "分段传播时间KL(h)": KL,
            "演算系数C0": C0,
            "演算系数C1": C1,
            "演算系数C2": C2,
            "演算出流O(取整)": O_int,
            "演算出流O(全精度)": [round(o, 4) for o in O_full],
        },
    }


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def _col_lines(label, vals, N):
    """原著 4 列布局（列优先）：每列长 S=(N+3)//4，共 S 行，
    第 r 行（r=0..S-1）含 r, r+S, r+2S, r+3S 号。"""
    S = (N + 3) // 4
    lines = []
    for r in range(S):
        cells = []
        for c in range(4):
            idx = r + c * S
            if idx < N:
                cells.append(f"{label}({idx + 1:2d})= {vals[idx]:8.1f}")
        lines.append("    ".join(cells))
    return lines


def render(params, result):
    """生成文本计算书（原著 .OUT 风格，4 列布局）。"""
    inp = result["输入"]
    res = result["结果"]
    N = inp["时段数N"]
    lines = []

    lines.append("（一）基本数据")
    lines.append("")
    lines.append("演算河段数   M= %2d            N= %2d"
                 % (inp["演算河段数M"], N))
    lines.append("流量比重因素 X= %.4f         传播时间 K=   %.4f小时"
                 % (inp["流量比重因素X"], inp["传播时间K(h)"]))
    lines.append("演算时段     T=   %.4f小时" % inp["演算时段T(h)"])
    lines.append("")
    lines.append(" 入流(秒立米)")
    lines.extend(_col_lines("I", inp["上断面入流I"], N))
    lines.append("")
    lines.append("（二）输出演算结果")
    lines.append("")
    lines.append(" 演算出流(秒立米)")
    lines.extend(_col_lines("O", res["演算出流O(取整)"], N))
    lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text", mode="theory"):
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
                   mode="direct" if "--direct" in sys.argv else "theory")
    print(txt)
