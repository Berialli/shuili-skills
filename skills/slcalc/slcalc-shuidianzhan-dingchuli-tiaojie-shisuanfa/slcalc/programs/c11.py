# -*- coding: utf-8 -*-
"""
C-11 无调节(日)水电站常规水能计算程序 —— 内核
================================================
复刻《水利程序集》C-11 程序（作者：黄璜，新疆水利水电勘测设计院）。

功能：将全部径流系列的日平均流量资料分成若干个流量阶梯进行统计计算，
并据此进行能量计算。输出：
  ① 流量~保证率~出力~发电时间~电量 阶梯统计表
  ② 不同装机容量(Nz)方案比较（年电量 / 发电流量 / 发电持续时间）
  ③ 电站保证出力与设计流量汇总

【2026-09-05 破译结论】（对照权威 C-11.OUT 库车河兰干水电站算例）
  1. 分档几何 = 五段等分（用户给的分档参数 D1,I1,D2,I2,...,I5）：
       段1 [Sjmin,D1] 等分 I1 档、…、段5 [D4,Sjmax] 等分 I5 档。
  2. 逐日流量归类：先做 VB6 单精度变换  x -> single(x * single(0.95))，
     再按 double 档界 [a,b) 半开归档（0 值流量天然落入首档下限之外，
     因 single(0*0.95)=0 < 1.0 而无档，恰与权威 n 合计 364/空档吻合）。
     —— 这是 double 运算下 2.0*0.95==1.9 会越界、而 VB single 下略小于
        1.9 归入首档这一“第 55 个流量”悬案的根因。
  3. 平均流量 QP = 档内 double(x*0.95) 均值，打印 2 位小数（half-away）。
     注意：归属用 single，均值用 double；两者仅影响尾数 5 舍入方向。
  4. 出力 Ni = A·H·QP/10000（万千瓦），打印 3 位。
  5. 频率 P = 已排除“低于本档下界”的天后剩余占比 = 本行门槛对应的
     保证率；发电时间 T = P×8760（T=365×24 起、逐档递减 24n）。
     即行 k 门槛为“第 k 档下界 a_k”：P_k = (365 - Σ_{m<k} n_m)/365。
  6. 电量：E_k = Σ_{m≤k} ΔNi_m·T_m（ΔNi_0 = Ni_0，ΔNi_m = Ni_m−Ni_{m-1}），
     打印单位亿 kWh；ΔE 打印行 0 = 0，其后 = ΔNi·T（亿）。
     —— 语义为“装机自 0 逐步抬高至 Ni_k 的累计电量”，权威 36 行全列逐值命中。
  7. 装机方案：给定 Nz（万千瓦），
      发电流量 Q  = round(Nz×10⁴/(A·H), 2)            （raw，不含 F）
      年电量 E    = E_j + (Nz − Ni_j)·T_{j+1}/10⁴（亿）  j = max{行: Ni_j≤Nz}
      —— 与权威 5 方案 Q/E 全部逐位命中。
     发电持续时间 T：权威为非 24 倍数小数小时，与任何离散主表/线性内插链
     不一致，未能闭式还原（详见报告），内核采用“主表 T 链两相邻行线性内插”。

数据文件（C-11.INT）顺序：
  P, A, H, F, Sjall, Sjmin, Sjmax, D1,I1,D2,I2,D3,I3,D4,I4,I5, 然后 365 个流量

验证基准：权威 C-11.OUT —— 主表 38 行（2 个空档）n 列 36/36 硬性全中、
非空档 36 行 QP/Ni/P/ΔNi/T/ΔE/E 逐值命中；装机比较 5 方案 Q、E 全中。
"""
import struct
import math
import os
from ..core.intio import read_numbers, smart_read_text
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "C-11"
TITLE = "无调节(日)水电站常规水能计算书 C-11"


def _f32(x):
    """VB6 Single 语义：double -> float32 -> double。"""
    return struct.unpack('f', struct.pack('f', x))[0]


def _fmt(v, nd):
    """VB FormatNumber half-away（正数）：floor(v*10^nd+0.5)/10^nd。"""
    if v is None:
        return None
    s = 10.0 ** nd
    return math.floor(v * s + 0.5) / s


def _build_bounds(Sjmin, Sjmax, D, Iseg):
    """五段等分构造 38 档 (lo, hi) 列表（double）。"""
    bounds = []
    lo = Sjmin
    for s in range(5):
        hi = D[s + 1] if s < 4 else Sjmax
        nsub = Iseg[s]
        w = (hi - lo) / nsub
        for k in range(nsub):
            bounds.append((lo + k * w, lo + (k + 1) * w))
        lo = hi
    return bounds


def _classify_single(x, bounds, F32):
    """VB6 single 语义归类：single(single(x)*single(F)) 落在 [lo,hi)。"""
    v = _f32(_f32(x) * F32)
    for i, (a, b) in enumerate(bounds):
        if a <= v < b:
            return i
    return None


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    if isinstance(data, dict):
        return dict(data)
    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}
    idx = 0
    p["P"] = nums[idx]; idx += 1          # 设计保证率
    p["A"] = nums[idx]; idx += 1          # 出力系数
    p["H"] = nums[idx]; idx += 1          # 设计净水头
    p["F"] = nums[idx]; idx += 1          # 水量利用系数
    p["Sjall"] = int(nums[idx]); idx += 1
    p["Sjmin"] = nums[idx]; idx += 1
    p["Sjmax"] = nums[idx]; idx += 1
    D = [None]
    D.append(nums[idx]); idx += 1       # D1
    I1 = int(nums[idx]); idx += 1       # I1
    D.append(nums[idx]); idx += 1       # D2
    I2 = int(nums[idx]); idx += 1       # I2
    D.append(nums[idx]); idx += 1       # D3
    I3 = int(nums[idx]); idx += 1       # I3
    D.append(nums[idx]); idx += 1       # D4
    I4 = int(nums[idx]); idx += 1       # I4
    I5 = int(nums[idx]); idx += 1       # I5
    Iseg = [I1, I2, I3, I4, I5]
    p["D"] = D
    p["Iseg"] = Iseg
    n = p["Sjall"]
    p["Q"] = nums[idx:idx + n]
    if len(p["Q"]) < n:
        raise ValueError(f"C-11 流量个数不足: 需 {n} 实 {len(p['Q'])}")
    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params, cfg=None):
    P = float(params["P"])
    A = float(params["A"])
    H = float(params["H"])
    F = float(params["F"])
    Sjmin, Sjmax = float(params["Sjmin"]), float(params["Sjmax"])
    Q = params["Q"]
    bounds = _build_bounds(Sjmin, Sjmax, params["D"], params["Iseg"])
    F32 = _f32(F)

    # --- 逐日归类（VB6 single 语义） ---
    members = [[] for _ in range(len(bounds))]
    for x in Q:
        c = _classify_single(x, bounds, F32)
        if c is not None:
            members[c].append(x)

    # --- 主表阶梯行（仅非空档打印） ---
    rows = []                 # 打印行序
    done_prev = 0             # 本行门槛(=档下界)前已排除天数
    Eacc = 0.0                # 万 kWh
    ni_prev = None
    for i, (a, b) in enumerate(bounds):
        if not members[i]:
            continue
        n = len(members[i])
        xf = [x * F for x in members[i]]            # double 均值语义
        qbar = sum(xf) / n
        ni = A * H * qbar / 10000.0                 # 万 kW
        Pct = (365 - done_prev) / 365.0             # 门槛=a 的保证率
        T = Pct * 8760                              # h（未含本档 n）
        dni_eff = ni if ni_prev is None else ni - ni_prev
        Eacc += dni_eff * T
        dE_show = dni_eff * T
        rows.append({
            "idx": i,
            "下界": a, "上界": b,
            "下界_raw": a, "上界_raw": b,
            "个数": n,
            "QP_raw": qbar, "QP": _fmt(qbar, 2),
            "Ni_raw": ni, "Ni": _fmt(ni, 3),
            "P_raw": Pct * 100.0, "P": _fmt(Pct * 100.0, 2),
            "T": int(T),
            "dNi_raw": 0.0 if ni_prev is None else ni - ni_prev,
            "dNi": 0.0 if ni_prev is None else _fmt(ni - ni_prev, 2),
            "dE_raw": dE_show, "dE": 0.0 if ni_prev is None else _fmt(dE_show / 10000.0, 3),
            "E_raw": Eacc, "E": _fmt(Eacc / 10000.0, 3),
            "成员": members[i],
        })
        done_prev += n
        ni_prev = ni

    # --- 装机容量方案比较 ---
    Nz_list = None
    if cfg and cfg.get("Nz") is not None:
        Nz_list = cfg["Nz"]
    if Nz_list is None:
        # 默认：无调节电站典型方案（可用 --cfg 覆盖）
        Nz_list = [0.08, 0.12, 0.16, 0.24, 0.36]
    schemes = []
    if Nz_list:
        Ni = [r["Ni_raw"] for r in rows]
        Eacc_list = [r["E_raw"] for r in rows]
        T_list = [r["T"] for r in rows]
        for Nz in Nz_list:
            # j = 最后一个 Ni<=Nz
            j = -1
            for k, niv in enumerate(Ni):
                if niv <= Nz:
                    j = k
            Q_disp = _fmt(Nz * 10000.0 / (A * H), 2)      # 发电流量（raw）
            if j < 0:
                E = Nz * T_list[0] / 10000.0              # 亿
                T_disp = T_list[0]
            else:
                E = Eacc_list[j] / 10000.0
                extra = (Nz - Ni[j]) * (T_list[j + 1] if j + 1 < len(T_list) else 0)
                E += extra / 10000.0
                if j + 1 < len(T_list):
                    # T 近似：主表 T 链相邻行线性内插（权威为非24倍数，未能闭式还原）
                    n0, n1 = Ni[j], Ni[j + 1]
                    if n1 > n0:
                        f = (Nz - n0) / (n1 - n0)
                        T_disp = T_list[j] + f * (T_list[j + 1] - T_list[j])
                    else:
                        T_disp = T_list[j + 1]
                else:
                    T_disp = T_list[j]
            schemes.append({
                "Nz": Nz,
                "Q": Q_disp,
                "E": _fmt(E, 3),
                "E_raw": E,
                "T": T_disp,
            })

    # --- 汇总：保证出力与设计流量（75% 例：Ni=0.082, Q=2.65, E=0.065） ---
    summary = None
    if P > 0:
        # 频率链 P_k(门槛=行 k 档下界 a_k, 单位%) 与行出力 Ni_k。
        # 保证出力 Ni_guar = 在 P=P% 处对 (Ni_k, P_k) 相邻行线性内插。
        Ni_arr = [r["Ni_raw"] for r in rows]
        P_arr = [r["P_raw"] for r in rows]          # %（门槛保证率）
        f_guar = None
        if P_arr[0] >= P * 100.0:
            for k in range(len(rows) - 1):
                if P_arr[k + 1] <= P * 100.0 <= P_arr[k]:
                    f_guar = (P_arr[k] - P * 100.0) / (P_arr[k] - P_arr[k + 1])
                    Ni_guar = Ni_arr[k] + f_guar * (Ni_arr[k + 1] - Ni_arr[k])
                    break
        if f_guar is None:
            Ni_guar = Ni_arr[-1]
            k = len(rows) - 2
        else:
            k = next(i for i in range(len(rows) - 1)
                     if P_arr[i + 1] <= P * 100.0 <= P_arr[i])
        # 设计流量 = 保证出力折算（raw，不含 F），与装机方案同式
        q_guar = Ni_guar * 10000.0 / (A * H)
        T_guar = P * 100.0 / 100.0 * 8760.0
        # 年电量：同装机 E 模型
        j = -1
        for ki, niv in enumerate(Ni_arr):
            if niv <= Ni_guar:
                j = ki
        E_arr = [r["E_raw"] for r in rows]
        if j < 0:
            E_guar = Ni_guar * T_guar / 10000.0
        else:
            E_guar = E_arr[j] / 10000.0
            if j + 1 < len(rows):
                E_guar += (Ni_guar - Ni_arr[j]) * rows[j + 1]["T"] / 10000.0
            else:
                E_guar += (Ni_guar - Ni_arr[j]) * rows[j]["T"] / 10000.0
        summary = {
            "P": P * 100.0,
            "Q": _fmt(q_guar, 2),
            "Ni": _fmt(Ni_guar, 3),
            "E": _fmt(E_guar, 3),
            "T": T_guar,
        }

    return {
        "程序": PROGRAM_ID,
        "P": P, "A": A, "H": H, "F": F,
        "Sjmin": Sjmin, "Sjmax": Sjmax,
        "D": params["D"], "Iseg": params["Iseg"],
        "bounds": bounds,
        "rows": rows,
        "schemes": schemes,
        "summary": summary,
    }


# ------------------------------------------------------------
# 输出
# ------------------------------------------------------------

def render(params, result, table=None):
    L = []
    L.append("一、工程名称:suanli")
    L.append("")
    L.append("二、计算参数")
    L.append(f"  1.电站设计保证率 P= {params['P']*100:.0f} %")
    L.append(f"  2.发电设备综合利用系数 A= {params['A']:.1f}")
    L.append(f"  3.电站设计净水头 H= {params['H']:.1f} 米")
    L.append(f"  4.考虑发电损失后水的利用系数 F= {params['F']:.2f}")
    L.append("")
    L.append("")
    L.append("                              一. 统计计算成果")
    L.append("")
    L.append("   流量阶梯   平均流量   出力  流量个数  频率   出力差值  发电时间  电量差值  累计电量")
    L.append("   (m^3/s)    (m^3/s)  (万kW)   (个)     (%)     (万kW)     (t)     (亿kWh)   (亿kWh)")
    for r in result["rows"]:
        a, b = r["下界"], r["上界"]
        L.append(
            f"{a:>6.1f}-{b:>7.1f}  {r['QP']:>7.2f}  {r['Ni']:>6.3f}  {r['个数']:>5d}  "
            f"{r['P']:>6.2f}  {r['dNi']:>6.2f}  {r['T']:>6.0f}  {r['dE']:>8.3f}  {r['E']:>8.3f}")
    L.append("")
    L.append("")
    L.append("                              二. 装机容量比较成果")
    L.append("")
    L.append("           方案    装机容量    年电量     发电流量   发电持续时间")
    L.append("                     (万kW)    (亿kW.h)   (m^3/s)       ( t )")
    for k, s in enumerate(result["schemes"]):
        L.append(f"             {k+1:>2d}     {s['Nz']:>7.3f}   {s['E']:>6.3f}   "
                 f"{s['Q']:>7.2f}   {s['T']:>8.0f}")
    L.append("")
    L.append("")
    L.append("                              三. 汇 总 成 果")
    L.append("")
    L.append("      设计保证率   设计流量   保证出力   装机容量   年发电量   发电持续时间")
    L.append("        ( % )       (m^3/s)    (万kW)     (万kW)    (亿kW.h)      ( t )")
    if result.get("summary"):
        s = result["summary"]
        L.append(f"         {s['P']:>7.2f}     {s['Q']:>5.2f}    {s['Ni']:>5.3f}   "
                 f"{s['Ni']:>6.3f}    {s['E']:>6.4f}     {s['T']:>6.0f}")
    L.append("")
    L.append("                        计算全部结束.谢谢 !")
    return render_text(PROGRAM_ID, TITLE, [("", L)])


def run(data, out_txt=None, out_json=None, fmt="text", cfg=None):
    params = parse(data)
    result = compute(params, cfg)
    text = render(params, result, None)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    cfg = None
    if len(sys.argv) > 2:
        nz = [float(v) for v in sys.argv[2].split(",")]
        cfg = {"Nz": nz}
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None, cfg=cfg)
    print(txt)
