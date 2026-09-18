# -*- coding: utf-8 -*-
"""
A-11 最大24小时洪量计算程序 —— 内核
=====================================
复刻《水利程序集》A-11 程序（作者：陈沂，水电部天津勘测设计院）。

功能：
  直接采用洪水水文要素摘录资料（时刻—流量过程），使用面积包围法计算
  最大 24 小时洪量。

计算方法（面积包围法；原著公式为图像，数值结构经 A-11.OUT 基准黑盒
     反推验证，2026-09-03 存档 a11_research1..4.py）：
  1. 摘录时刻 T(I) 采用"时.分"编码：整数部分为时，小数×100 为分
     （如 17.5 = 17:50、0.42 = 0:42、1.18 = 1:18、23.21 = 23:21）。
     时刻从 DO 日 UO 时起摘录；数值发生"回跳"（后点小于前点）表示跨入
     次日，绝对时刻 +24h（+1440min）。UO 即摘录起点钟点（=第一个 T）。
  2. 逐摘录段梯形累计水量，得摘录点处累积水量曲线 C(t)：
       C(T(i)) = Σ 段梯形  = C(T(i−1)) + (Y(i−1)+Y(i))/2 · ΔT
     C 的单位取 m3/s·min（后续换算按量纲齐性折算万 m3）。
  3. "求累积水量过程，用直线内插法"：C 在相邻摘录点间视为直线
     （内插值按线性比例）；即每小时洪量由内插后的累积曲线差分得出，
     而非对流量折线直接积分——当整点落在某两个摘录点之间时，该区间
     的水量被线性均摊到各小时（这是面积包围法的关键，已黑盒校准）。
  4. 逐时段（每小时）洪量 W(I)（I=1..W，W 为需内插个数）：
       起讫整点为 H(I−1)=UO+(I−1)h 与 H(I)=UO+I·h
       W(I) = round( [C(H(I)) − C(H(I−1))] · 0.006 )   （单位：万 m3）
     换算系数：C 差单位 m3/s·min，除以 60 得 m3/s·h 当量，再 ×0.36
     得万 m3，即 ×0.006；舍入与原著逐时洪量表一致（half-even 由反推
     确认，598.5→598、805.5→806；对 .5 依赖二进制值，实测命中）。
  5. 24h 累计 WW = 连续 24 个逐时洪量之和，滑动取值，取最大为 W24：
       W24 = max Σ W(s..s+23)，s=1..(W−23)
     （原著输出"序号 No=1-24"，基准窗口即起于第一个整点 UO。）
  6. 结果打印：W24 取整（万 m3）。

验证基准（A-11.INT 算例，原著 A-11.OUT）：
  N=38、W=35、DO=23、UO=17（23日17时起摘录）；
  摘录覆盖绝对 [17h,52h)（至 25 日 4 时），需内插 35 个整点小时。
  逐时洪量表 35 行与基准逐位一致（35/35 ±0），
  最大 24h 洪量 W24 = 10142 万 m3（±0，窗口 No=1..24）。
  黑盒校准点：逐时 5–6 与 6–7 时同为 509（摘录 [5,7) 内无整点摘录点，
  区间水量经 C 内插均分）；17–18(24日) 等 3 段同为 169，同理由
  [17,20) 段均摊；两处均与流量梯形积分的直接分划（596/421、180）
  不同，从而唯一确定面积包围法的"累积曲线直线内插"语义。
"""
import os

from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-11"
TITLE = "最大24小时洪量计算书"

# 换算系数：C 差(m3/s·min) × 0.006 = 万 m3
#   C差/60 = 平均流量(m3/s)相当小时水量(m3/s·h)，×0.36 万 m3·h/(m3/s)
CF = 0.006


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
    """解析 A-11.INT。data: INT 路径 | 文本 | dict。

    数据文件顺序（原著）：
      N, W
      T(1), ..., T(N)      （时.分 编码；跨日回跳）
      Y(1), ..., Y(N)      （对应时刻流量 m3/s）
      DO, UO               （初始日期、初始内插时刻）
    各数值以逗号分隔，可一行多值/多行；兼容含 ^Z(0x1A) 的旧式文本。
    dict 形式键：N, W, T(list), Y(list), DO, UO。
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
        raise ValueError("A-11 无数据内容")

    nums = [x for ln in lines for x in _tokenize(ln)]
    if len(nums) < 4:
        raise ValueError("A-11 数据不足（需 N,W 与至少 1 个 T/Y）")

    N = int(nums[0])
    W = int(nums[1])
    if N <= 0 or W <= 0:
        raise ValueError(f"A-11 N={N}、W={W} 非法（需 ≥1）")
    need = 2 + 2 * N + 2
    if len(nums) < need:
        raise ValueError(f"A-11 数据不足：需 N,W + {N} 个 T + {N} 个 Y + DO,UO，"
                         f"收到数值 {len(nums)} 个（需 ≥{need}）")

    T = [float(nums[2 + j]) for j in range(N)]
    Y = [float(nums[2 + N + j]) for j in range(N)]
    DO = int(nums[2 + 2 * N])
    UO = float(nums[2 + 2 * N + 1])

    # 摘录时刻必须单调（允许跨日回跳）
    if any(Y[j] < 0 for j in range(N)):
        raise ValueError("A-11 流量 Y 含负值")

    return {"N": N, "W": W, "T": T, "Y": Y, "DO": DO, "UO": UO}


def _parse_dict(d):
    return {
        "N": int(d["N"]), "W": int(d["W"]),
        "T": [float(x) for x in d["T"]],
        "Y": [float(x) for x in d["Y"]],
        "DO": int(d.get("DO", d.get("d0", 1))),
        "UO": float(d.get("UO", d.get("uo", d["T"][0]))),
    }


# ------------------------------------------------------------
# 时刻解码
# ------------------------------------------------------------

def hhmm(t):
    """'时.分'编码 → 绝对分钟（小数×100 为分钟）。"""
    h = int(t)
    m = round((t - h) * 100)
    if m >= 60:                       # 容错：分=60 进位
        h += 1
        m -= 60
    return h * 60 + m


def norm_times(T):
    """摘录时刻序列 → 绝对分钟序列（回跳视为跨日 +1440min）。"""
    out, base, prev = [], 0, None
    for v in T:
        if prev is not None and v < prev - 1e-9:
            base += 1440
        out.append(base + hhmm(v))
        prev = v
    return out


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _bank_round(v):
    """与原著逐时打印一致的舍入。用 round()（half-even）；
    本算例 598.5→598、805.5→806，实测与 Python round 一致。"""
    return round(v)


def compute(params):
    """执行 A-11。params: parse 返回 dict。返回结构化结果。"""
    N = params["N"]
    W = params["W"]
    T = params["T"]
    Y = params["Y"]
    UO = params["UO"]
    if len(T) != N or len(Y) != N:
        raise ValueError(f"A-11 T/Y 应各有 N={N} 个值，收到 {len(T)}/{len(Y)}")
    if W > 24 * 10:      # 防御：小时数过大无意义
        raise ValueError(f"A-11 W={W} 过大")

    ta = norm_times(T)

    # 摘录点累计水量 C（m3/s·min，梯形累加）
    cum = [0.0]
    for i in range(N - 1):
        if ta[i + 1] <= ta[i]:
            raise ValueError(f"A-11 摘录时刻非增（第 {i + 1}/{i + 2} 点）")
        cum.append(cum[-1] + (Y[i] + Y[i + 1]) / 2.0 * (ta[i + 1] - ta[i]))

    def C_at(t):
        """累积水量曲线直线内插（绝对分钟 t → m3/s·min）。"""
        if t <= ta[0]:
            if N == 1:
                return cum[0]
            return cum[0] + (cum[1] - cum[0]) * (t - ta[0]) / (ta[1] - ta[0])
        if t >= ta[-1]:
            return cum[-1]
        lo, hi = 0, N - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if ta[mid] <= t:
                lo = mid
            else:
                hi = mid
        x1, x2 = ta[lo], ta[hi]
        if x2 == x1:
            return cum[lo]
        return cum[lo] + (cum[hi] - cum[lo]) * (t - x1) / (x2 - x1)

    # 逐整点小时（起于 UO，共 W 个内插小时）
    h0 = hhmm(UO)
    Wi_full = []
    Wi_int = []
    for k in range(W):
        a = h0 + k * 60
        b = a + 60
        v = (C_at(b) - C_at(a)) * CF          # 万 m3
        Wi_full.append(v)
        Wi_int.append(_bank_round(v))

    # 24h 滑动累计（对取整后的逐时洪量求和，与原著一致）
    windows = []
    for s in range(W - 24 + 1):
        windows.append(sum(Wi_int[s:s + 24]))
    best = max(range(len(windows)), key=lambda i: windows[i])
    W24 = windows[best]

    # 逐时行附日期标签：绝对小时 A = UO_whole + k；钟点 = A mod 24；
    # 日期 = DO + A//24（UO 属 DO 日；A 越过 0 点即跨自然日）
    base_day = int(params["DO"])
    labels = []
    for k in range(W):
        A = h0 // 60 + k
        labels.append((base_day + A // 24, A))

    return {
        "程序": PROGRAM_ID,
        "输入": {
            "数据组数N": N,
            "需内插个数W": W,
            "时刻T(时.分)": T,
            "流量Y(m3/s)": Y,
            "初始日期DO": params["DO"],
            "初始内插时刻UO": UO,
        },
        "中间量": {
            "摘录绝对分钟ta": ta,
            "摘录点累计水量C(m3/s·min)": [round(c, 3) for c in cum],
        },
        "结果": {
            "逐时洪量表(Wi)": Wi_int,          # 长度 W，单位万 m3
            "逐时洪量全精度": [round(v, 4) for v in Wi_full],
            "日期标签": labels,
            "24h累计候选": windows,
            "最大窗口序号No": best + 1,
            "最大24小时洪量W24(万m3)": W24,
        },
    }


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def render(params, result):
    """生成文本计算书（原著 .OUT 风格）。"""
    inp = result["输入"]
    res = result["结果"]
    N = inp["数据组数N"]
    T = inp["时刻T(时.分)"]
    Y = inp["流量Y(m3/s)"]
    lines = []

    lines.append("（一）基本数据：")
    lines.append("")
    lines.append("序号  时刻系列数   流量系列数    序号   时刻系列数   流量系列数")
    lines.append("")
    half = (N + 1) // 2
    for r in range(half):
        cells = []
        for c in range(2):
            i = r + c * half
            if i < N:
                cells.append(f"{i + 1:2d}      {T[i]:9.2f}    {Y[i]:9.2f}")
        lines.append("   ".join(cells))
    lines.append("")
    lines.append("日期  序号     时间        每小时洪量")
    lines.append("---------------------------------------")

    labels = res["日期标签"]
    Wi = res["逐时洪量表(Wi)"]
    # 原著按自然日分组打印：每行 <时刻起> - <时刻讫>:00，逢新日期先空一行，
    # 新日期首行行首印该日期号（其后同日期行行首留空）。讫点为 0 时显示 24。
    def fmt_clock(h):
        return "24" if h % 24 == 0 else f"{h % 24:2d}"

    last_day = None
    for k in range(len(Wi)):
        day, A = labels[k]
        hh = A % 24
        hh2 = (A + 1) % 24
        if k > 0 and day != last_day:
            lines.append("")
        if day != last_day:
            lines.append(f"{day:3d}  {k + 1:3d}    {hh:2d} - {fmt_clock(hh2)}:00    {Wi[k]:5d}")
        else:
            lines.append(f"      {k + 1:3d}    {hh:2d} - {fmt_clock(hh2)}:00    {Wi[k]:5d}")
        last_day = day

    lines.append("")
    lines.append(f"  序号 No = {res['最大窗口序号No']:2d}  - {res['最大窗口序号No'] + 23:2d}")
    lines.append(f"  最大24小时洪量 W 24 =  {res['最大24小时洪量W24(万m3)']}万立米")
    lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 路径 | dict"""
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
