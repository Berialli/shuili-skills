# -*- coding: utf-8 -*-
"""
A-15 由直方图限定的曲线拟合程序 —— 内核
===========================================
复刻《水利程序集》A-15 程序（作者：潘东海、朱凤娟，水利部天津勘测设计研究院；
算法见潘东海、刘世芬《由直方图限定的曲线拟合》，中国民航学院学报 1992 年第 4 期）。

功能：
  将直方图（分段常值的矩形序列）还原成一条连续曲线，满足：
    ① 每个子区间 [X_{j−1}, X_j] 上曲线与 X 轴所围面积 = 该区间矩形面积 Q_j；
     ② 两端满足边界条件（初值/末值标高及其导数）；
     ③ 各段间曲线的衔接具有不高于一阶的连续性（C1，可控）。
  当边界条件（端导）与内部衔接互相冲突（超定）时，以最小二乘意义求解，
  权值 α、β（α+β=1，键盘 ALF/BT）分别控制边界条件与内部衔接的权重。

算法（原著公式为 Equation.3 OLE 对象，经 MTEF token 解码 + 例1/例2 101 点
     黑盒反推全链验证，2026-09-03 存档 a15_*.py 研究脚本）：

  ⚠️ 黑盒反推发现的权重修正（关键）：
     原著说明书给出目标 min α(M_0²+M_N²) + β·Σ_j L_j²，但以例2
     （ALF=0.3、BT=0.7）权威输出反推，实际求解的等效目标为
       min 4α(M_0²+M_N²) + β·Σ_j L_j²
     （等效权重比 4α:β = 1.2:0.7；数值反推 r*=1.714283 ≈ 12/7 精确吻合）。
     可能因 L_j（eq_07）实际定义为 C1 差之半（半差平滑），平方后产生 1/4，
     使 β 项等效缩为 β/4 —— 对解的影响等价于端导加权翻倍。α=0 时不受影响
     （例1 仍为纯 C1 系统，101/101 逐位命中），两例自洽。
  1. 每子区间 [X_{j−1},X_j] 内取二次曲线（归一化坐标 t=(X−X_{j−1})/H_j，
     H_j = X_j − X_{j−1}）：
       S_j(X) = Y_{j−1}(1−t) + Y_j·t + C_j·t(1−t)
     其中 Y_{j−1}=S_j(X_{j−1})、Y_j=S_j(X_j) 为节点曲线标高（待求内部 Y），
     二次项 C_j 由面积守恒条件（∫S_j dX = Q_j = 柱高 y_j·H_j）解出：
       C_j = 6Q_j/H_j − 3Y_{j−1} − 3Y_j = 6·y_j − 3Y_{j−1} − 3Y_j
  2. 未知量为内部节点 Y_1..Y_{N−1}（共 N−1 个；Y_0、Y_N 为已知边界标高）。
     段间 C1 连续（L_j = 0，j=1..N−1）给出 N−1 个方程：
       (Y_j−Y_{j−1}−C_j)/H_j = (Y_{j+1}−Y_j+C_{j+1})/H_{j+1}   （左导=右导）
     两端导数条件（M_0、M_N）为额外 2 个方程：
       M_0: 段 1 左端导 − Y'_0 = 0
       M_N: 段 N 右端导 − Y'_N = 0
     合共 N+1 个方程、N−1 个未知数 → 超定 → 最小二乘：
       min α(M_0² + M_N²) + β Σ_j L_j²
  3. α=ALF、β=BT 为权值（键盘输入，α+β=1）。α=0 时纯内部衔接（边界导
     自由）；β=0 时纯边界条件。例 1：ALF=0、BT=1；例 2：ALF=0.3、BT=0.7。
  4. 输出：从 X_0 到 X_N 每段 10 等分共 10N+1 个采样点（格式 F10.2/F8.2，
     曲线标高保留 2 位小数），以及"直方图面积"（ΣQ_j）与
     "直方图所限定曲线面积"（对输出 2 位小数的采样点逐点梯形积分）。

验证基准（A-15-1.INT 算例，原著 A-15-1.OUT 逐点核对）：
  N=10；X=(2,4,5,7.5,10,12,14.5,17,18,20,22)；y=(2,4,6,5,3,2,3,4,6,5)
  Y0=1.1、YM=3.4、Y1=1、Y2=−1（Y1/Y2 为端导，例1 α=0 不使用）
  KDX=2、KDY=1、YMIN=0、ALF=0、BT=1
  → 内部节点 LSQ 解 Y=(3.20128,4.84615,6.07376,3.85881,2.24928,2.30472,
     3.53185,4.78892,6.20277)，输出 2 位 (3.20,4.85,6.07,3.86,2.25,2.30,
     3.53,4.79,6.20)；101 个采样点（每段 10 等分）逐点 2 位小数
     与 A-15-1.OUT 完全一致（101/101）。
  曲线面积：直方图面积 80.000；曲线面积 79.976 = 输出采样点梯形积分。

  验证基准 2（A-15-2.INT 算例，权威 OUT 为 RTF 内嵌文本提取）：
  N=10；X=(20,40,50,75,100,120,145,170,180,200,220)；y=(200,400,600,500,
  300,200,300,400,600,500)；Y0=110、YM=340、Y1=1、Y2=−1；KDX=20、KDY=100、
  YMIN=0、ALF=0.3、BT=0.7。
  → 按等效 4α:β 加权解：9 节点与权威显示值全中（max|Δ|=0.0044<0.005），
    101 采样点 100/101 逐位一致（仅 X=90 全精度 479.015 处于 2 位舍入边界
    差 0.01），曲线面积 79986.188 vs 权威 79986.148（相对差 5e-7）。
  → 若用输入 α:β=0.3:0.7 直解，末节点 Y9=625.94 vs 权威 641.06（差 15）——
    证明 4 倍修正必要。两例自洽，算法链闭合。

实现说明：
  * 内部节点用带截距线性最小二乘一次解出（残差对节点线性）。
  * 段内二次曲线、C_j 公式、C1 衔接方程、加权 LSQ 均严格按原著公式。
  * 采样点标高按原著 2 位小数输出（VB Format "0.00"）。
"""
import os

import numpy as np

from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-15"
TITLE = "由直方图所限定的曲线拟合计算书"


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
    return [_num(t) for t in line.replace(";", ",").split(",") if t.strip()]


def parse(data):
    """解析 A-15.INT。data: INT 路径 | 文本 | dict。

    数据文件顺序（原著）：
      第 1 行 N（时段个数）
      第 2 行 X(1..N+1)  各节点水平坐标
      第 3 行 y(1..N)    各直方块高度（垂直坐标）
      第 4 行 Y0, YM, Y1, Y2 （初值标高、终值标高、初值处导数、终值处导数）
      第 5 行 KDX, KDY, YMIN, ALF, BT （键盘输入：水平刻度、垂直刻度、
              垂直坐标原点、权值 α、权值 β；α+β=1）
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
        raise ValueError("A-15 无数据内容")
    nums = [x for ln in lines for x in _tokenize(ln)]
    if len(nums) < 7:
        raise ValueError("A-15 数据不足（需 N, X(≥2), y, Y0/YM/Y1/Y2）")

    N = int(nums[0])
    if N <= 0:
        raise ValueError(f"A-15 时段个数 N={N} 非法")
    need = 1 + (N + 1) + N + 4
    if len(nums) < need:
        raise ValueError(f"A-15 数据不足：需 1+{N+1}+{N}+4 个数，收到 {len(nums)}")
    X = [float(nums[1 + i]) for i in range(N + 1)]
    y = [float(nums[1 + (N + 1) + i]) for i in range(N)]
    off = 1 + (N + 1) + N
    Y0, YM, Y1, Y2 = (float(nums[off]), float(nums[off + 1]),
                      float(nums[off + 2]), float(nums[off + 3]))
    # 可选第 5 行键盘参数
    KDX, KDY, YMIN, ALF, BT = 1.0, 1.0, 0.0, 0.0, 1.0
    if len(nums) >= off + 4 + 1:
        rest = nums[off + 4:]
        # 容错：可能缺 YMIN 只有 KDX,KDY,ALF,BT 等
        if len(rest) >= 5:
            KDX, KDY, YMIN, ALF, BT = (float(rest[0]), float(rest[1]),
                                       float(rest[2]), float(rest[3]),
                                       float(rest[4]))
        elif len(rest) == 4:
            KDX, KDY, ALF, BT = (float(rest[0]), float(rest[1]),
                                 float(rest[2]), float(rest[3]))
    for i in range(N):
        if X[i + 1] <= X[i]:
            raise ValueError(f"A-15 X 坐标须严格递增，X[{i+1}]={X[i+1]} <= X[{i}]={X[i]}")
    return {"N": N, "X": X, "y": y,
            "Y0": Y0, "YM": YM, "Y1": Y1, "Y2": Y2,
            "KDX": KDX, "KDY": KDY, "YMIN": YMIN,
            "ALF": ALF, "BT": BT}


def _parse_dict(d):
    return {"N": int(d["N"]),
            "X": [float(x) for x in d["X"]],
            "y": [float(x) for x in d["y"]],
            "Y0": float(d["Y0"]), "YM": float(d["YM"]),
            "Y1": float(d["Y1"]), "Y2": float(d["Y2"]),
            "KDX": float(d.get("KDX", 1.0)), "KDY": float(d.get("KDY", 1.0)),
            "YMIN": float(d.get("YMIN", 0.0)),
            "ALF": float(d.get("ALF", 0.0)), "BT": float(d.get("BT", 1.0))}


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _residuals(Yv, H, yq, Yp0, Ypn, ALF, BT):
    """残差向量：内部 C1 衔接 L_j 与（加权）端导条件。

    Yv: 全节点标高 [Y0, Y1..Y_{N-1}, YN]
    H : 各段宽
    yq: 各段柱高
    返回 sqrt(BT_w)*L (N-1), sqrt(ALF_w)*M0, sqrt(ALF_w)*MN

    权重修正（黑盒反推）：实际目标 = 4α(M0²+MN²) + β·ΣL²，
    故 ALF_w:BT_w = 4ALF:BT（归一化）。
    """
    ALF_e, BT_e = 4.0 * ALF, BT          # 等效权重
    s = ALF_e + BT_e
    if s <= 0:
        ALF_w, BT_w = 0.0, 1.0
    else:
        ALF_w, BT_w = ALF_e / s, BT_e / s
    N = len(H)
    C = 6.0 * np.asarray(yq) - 3.0 * Yv[:-1] - 3.0 * Yv[1:]   # 面积守恒二次项
    L = np.zeros(N - 1)
    for j in range(N - 1):
        L[j] = ((Yv[j + 1] - Yv[j] - C[j]) / H[j]
                - (Yv[j + 2] - Yv[j + 1] + C[j + 1]) / H[j + 1])
    M0 = (Yv[1] - Yv[0] + C[0]) / H[0] - Yp0        # 段1 左端导 − Y'0
    MN = (Yv[N] - Yv[N - 1] - C[N - 1]) / H[N - 1] - Ypn  # 段N 右端导 − Y'N
    return np.sqrt(BT_w) * L, np.sqrt(ALF_w) * M0, np.sqrt(ALF_w) * MN


def _solve_nodes(X, yq, Y0, YN, Yp0, Ypn, ALF, BT):
    """解内部节点 Y_1..Y_{N-1}（残差对节点线性 → 一次最小二乘）。"""
    N = len(yq)
    n = N - 1
    H = np.diff(np.asarray(X, float))
    yq = np.asarray(yq, float)

    def r_of(Yin):
        Yv = np.concatenate(([Y0], Yin, [YN]))
        L, M0, MN = _residuals(Yv, H, yq, Yp0, Ypn, ALF, BT)
        return np.concatenate((L, [M0], [MN]))

    y0v = np.zeros(n)
    r0 = r_of(y0v)
    A = np.array([r_of(np.eye(1, n, k)[0]) - r0 for k in range(n)]).T
    sol, *_ = np.linalg.lstsq(A, -r0, rcond=None)
    return sol, H, yq


def _sample_curve(X, Yv, C, per_seg=10):
    """每段 per_seg 等分采样（含首尾共享端点）→ (xs, ys_raw 全精度)。"""
    N = len(Yv) - 1
    H = np.diff(np.asarray(X, float))
    xs, ys = [], []
    for j in range(N):
        a, b = X[j], X[j + 1]
        for i in range(per_seg):
            t = i / per_seg
            x = a + t * (b - a)
            xs.append(x)
            ys.append(Yv[j] * (1 - t) + Yv[j + 1] * t + C[j] * t * (1 - t))
    xs.append(X[-1])
    ys.append(Yv[-1])
    return np.array(xs), np.array(ys)


def compute(params, mode="original"):
    """执行 A-15。返回结构化结果 dict。

    mode "original" / "theory" 等价（A-15 无反演悬案）。
    """
    N = params["N"]
    X = params["X"]
    y = params["y"]
    Y0, YM = params["Y0"], params["YM"]
    Y1, Y2 = params["Y1"], params["Y2"]
    ALF, BT = params["ALF"], params["BT"]
    if len(X) != N + 1:
        raise ValueError(f"A-15 X 应有 N+1={N+1} 个，收到 {len(X)}")
    if len(y) != N:
        raise ValueError(f"A-15 y 应有 N={N} 个，收到 {len(y)}")
    if ALF < 0 or BT < 0:
        raise ValueError(f"A-15 权值 ALF/BT 应非负：ALF={ALF}, BT={BT}")

    H = np.diff(np.asarray(X, float))
    yq = np.asarray(y, float)
    Q = yq * H                       # 各矩形块面积

    # ① 解内部节点
    Yin, Hh, _ = _solve_nodes(X, yq, Y0, YM, Y1, Y2, ALF, BT)
    Yv = np.concatenate(([Y0], Yin, [YM]))

    # ② 面积守恒二次项 C_j
    C = 6.0 * yq - 3.0 * Yv[:-1] - 3.0 * Yv[1:]

    # ③ 每段 10 等分采样
    xs, ys = _sample_curve(X, Yv, C, per_seg=10)
    ys_r2 = [round(v, 2) for v in ys]          # 原著 2 位小数输出

    # ④ 面积
    hist_area = float(np.sum(Q))
    # 曲线面积 = 输出采样点（2位）逐点梯形积分（原著）
    curve_area = float(np.trapezoid(ys_r2, xs)) if hasattr(np, "trapezoid") \
        else float(np.trapz(ys_r2, xs))

    # C1 检查（L 残差）
    L, M0, MN = _residuals(Yv, Hh, yq, Y1, Y2, ALF, BT)

    # 每段面积校验（解析）
    seg_area = Hh * ((Yv[:-1] + Yv[1:]) / 2.0 + C / 6.0)

    return {
        "程序": PROGRAM_ID,
        "模式": f"加权最小二乘解节点（输入 α={ALF:g}、β={BT:g}；"
                f"等效 4α:β={4*ALF:g}:{BT:g}）",
        "输入": {
            "时段个数N": N,
            "节点坐标X": [round(x, 4) for x in X],
            "直方块高y": [round(v, 4) for v in y],
            "初值标高Y0": Y0, "终值标高YM": YM,
            "初值处导数Y1": Y1, "终值处导数Y2": Y2,
            "权值ALF(α)": ALF, "权值BT(β)": BT,
            "水平刻度KDX": params["KDX"], "垂直刻度KDY": params["KDY"],
            "垂直原点YMIN": params["YMIN"],
        },
        "结果": {
            "节点标高Y(0..N)": [round(v, 6) for v in Yv],
            "节点标高(2位)": [round(v, 2) for v in Yv],
            "二次项C_j": [round(v, 6) for v in C],
            "采样点数": len(xs),
            "采样点X": [round(v, 4) for v in xs],
            "采样点标高Y(2位)": ys_r2,
            "直方图面积": round(hist_area, 3),
            "直方图所限定曲线面积": round(curve_area, 3),
        },
        "中间量": {
            "段宽H_j": [round(v, 6) for v in Hh],
            "矩形面积Q_j": [round(v, 6) for v in Q],
            "段面积(解析校验)": [round(v, 6) for v in seg_area],
            "内部衔接残差L_j": [0.0 if abs(v) < 1e-6 else round(float(v), 6) for v in L],
            "端导残差M0": float(M0), "端导残差MN": float(MN),
            "内部节点(全精度)": [round(v, 8) for v in Yin],
        },
    }


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def render(params, result):
    """生成文本计算书（原著 .OUT 风格：逐采样点两列 10 段）。"""
    inp = result["输入"]
    res = result["结果"]
    N = inp["时段个数N"]
    Xn = inp["节点坐标X"]
    Yn = inp["直方块高y"]
    Xs = res["采样点X"]
    Ys = res["采样点标高Y(2位)"]

    lines = []
    lines.append("一、基本数据")
    lines.append(f"  时段个数 N = {N}")
    lines.append("  节点水平坐标 X(I)：")
    lines.append("  " + "  ".join(f"{x:9.3f}" for x in Xn))
    lines.append("  直方块垂直坐标 y(I)：")
    lines.append("  " + "  ".join(f"{v:9.3f}" for v in Yn))
    lines.append("")
    lines.append(f"  初值标高 Y0={inp['初值标高Y0']:g}   终值标高 YM={inp['终值标高YM']:g}")
    lines.append(f"  初值处导数 Y1={inp['初值处导数Y1']:g}   终值处导数 Y2={inp['终值处导数Y2']:g}")
    lines.append(f"  权值 α(ALF)={inp['权值ALF(α)']:g}   β(BT)={inp['权值BT(β)']:g}")
    lines.append("")
    lines.append("二、曲线采样点（每段 10 等分）")
    lines.append("")
    hdr = f"{'I':>4} {'X(I)':>10} {'Y(I)':>10}   {'I':>4} {'X(I)':>10} {'Y(I)':>10}"
    lines.append(hdr)
    lines.append("-" * 62)
    # 双列输出（原著 2 列，I 从 1 到 101）
    n = len(Xs)
    col = (n + 1) // 2
    for r in range(col):
        i1 = r
        i2 = col + r
        if i2 < n:
            lines.append(f"{i1+1:4d} {Xs[i1]:10.3f} {Ys[i1]:10.2f}"
                         f"   {i2+1:4d} {Xs[i2]:10.3f} {Ys[i2]:10.2f}")
        else:
            lines.append(f"{i1+1:4d} {Xs[i1]:10.3f} {Ys[i1]:10.2f}")
    lines.append("-" * 62)
    lines.append("")
    lines.append(f" 直方图面积                {res['直方图面积']:8.3f}")
    lines.append(f" 直方图所限定曲线面积      {res['直方图所限定曲线面积']:8.3f}")
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
