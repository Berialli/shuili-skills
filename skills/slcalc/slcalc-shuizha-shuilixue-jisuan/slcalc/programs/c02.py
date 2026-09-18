# -*- coding: utf-8 -*-
"""
C-2 水库调洪演算的数值解程序 —— 内核
======================================
复刻《水利程序集》C-2 程序（作者：张校正，新疆水利厅）。

原理：水库水量平衡微分方程
    F(z) * dz/dt = Q(t) - q(z)
  F(z) —— 水位~水面面积关系（线性插值）
  Q(t) —— 入库洪水过程（给定）
  q(z) —— 出库流量 = 泄洪洞 + 溢洪道 + 电站/灌溉常流量
用定步长四阶龙格-库塔法求解水位过程 z(t)。

泄流公式（经原始 OUT 回归验证，误差 <0.1%）：
  泄洪洞（深孔，孔口出流）:
    Q = M1 * A * sqrt(2g * (Z - 孔口顶高程))     A = B1*A1（全开）
  溢洪道（堰流）:
    Q = M2 * B * sqrt(2g) * (Z - C2)^1.5

数据文件顺序（C-2.INT / C-2X-1.INT）：
  工程名, 频率
  K                    （水位~面积曲线节点数）
  ZP,F, ZP,F, ...      （K 对）
  J, T                 （洪水时段数，时段间隔秒）
  QL(0..J)             （J+1 个洪水流量）
  ZM, Z0, Q0, W        （防洪下限水位，起始水位，泄洪起始流量，电站常流量）
  QAAA                 （下游安全限泄对个数）
  QAZ,QAQ, ...         （QAAA 对：水位,限泄流量）
  G0                   （泄洪洞个数）
  B1,A1,C1,M11,B11,A11,C11,QX, ...  （G0 组）
  H0                   （溢洪道个数）
  M2,B2,C2,GMH, ...    （H0 组）
  KK1, ...             （变宽变高泄流孔，暂支持 0）
  KK0, ...             （变宽溢洪道，暂支持 0）
  KK2, ...             （其它泄流方式曲线，暂支持 0）

验证基准：C-2X-1.INT（算例1，1/1000 频率）
  时段12: Z=45.77  洞=592.94  溢=925.17  总=1578.11
  时段43: Z=50.88  洞=656.00  溢=3126.73  总=3842.73
"""
import math

from ..core.numext import rk4_step
from ..core.intio import read_numbers
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "C-2"
TITLE = "水库调洪演算数值解计算书 C-2"

G = 9.81


# ------------------------------------------------------------
# 插值
# ------------------------------------------------------------

def interp1d(xs, ys, x):
    """线性插值（xs 升序）。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    i = 1
    while i < len(xs) - 1 and xs[i] < x:
        i += 1
    t = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
    return ys[i - 1] + t * (ys[i] - ys[i - 1])


# ------------------------------------------------------------
# 泄流设施
# ------------------------------------------------------------

def make_outflow_func(params):
    """
    构造出库流量函数 q(z)。
    params 含 holes（泄洪洞列表）、weirs（溢洪道列表）、W（常流量）、
    limits（安全限泄 [(水位,限泄)]）。
    """
    holes = params.get("holes", [])
    weirs = params.get("weirs", [])
    W = params.get("W", 0.0)

    def q_holes(z):
        q = 0.0
        for h in holes:
            # 孔口顶高程 = 闸门底高程 C11 + 闸门高 A11
            crest = h["C11"] + h["A11"]
            if z > crest:
                q += h["M11"] * (h["B1"] * h["A1"]) * math.sqrt(2 * G * (z - crest))
        return q

    def q_weirs(z):
        q = 0.0
        for w in weirs:
            if z > w["C2"]:
                q += w["M2"] * w["B2"] * math.sqrt(2 * G) * (z - w["C2"]) ** 1.5
        return q

    def qout(z):
        q = q_holes(z) + q_weirs(z) + W
        # 安全限泄：水位超过限泄水位后，总泄量不超过限泄值
        for z_limit, q_limit in params.get("limits", []):
            if z >= z_limit:
                q = min(q, q_limit)
        return q

    return qout, q_holes, q_weirs


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """
    解析输入。
    data: dict 或 INT 文件路径。
    返回标准参数字典。
    """
    if isinstance(data, dict):
        return data

    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}

    idx = 0
    # 工程名与频率在数值流之前，需从原始文本读
    from ..core.intio import smart_read_text
    text = smart_read_text(data)
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    # 首行: 工程名,频率
    if lines and "," in lines[0]:
        head = lines[0].split(",")
        p["工程名"] = head[0].strip()
        p["频率"] = head[1].strip() if len(head) > 1 else ""

    p["K"] = int(nums[idx]); idx += 1
    zps, fs = [], []
    for _ in range(p["K"]):
        zps.append(nums[idx]); fs.append(nums[idx + 1]); idx += 2
    p["zps"] = zps
    p["fs"] = fs  # 单位：万平方米 → 计算时 ×1e4

    p["J"] = int(nums[idx]); p["T"] = nums[idx + 1]; idx += 2
    p["QL"] = nums[idx:idx + p["J"] + 1]; idx += p["J"] + 1

    p["ZM"] = nums[idx]; p["Z0"] = nums[idx + 1]
    p["Q0"] = nums[idx + 2]; p["W"] = nums[idx + 3]; idx += 4

    qaaa = int(nums[idx]); idx += 1
    p["limits"] = []
    for _ in range(qaaa):
        p["limits"].append((nums[idx], nums[idx + 1])); idx += 2

    g0 = int(nums[idx]); idx += 1
    p["holes"] = []
    for _ in range(g0):
        p["holes"].append({
            "B1": nums[idx], "A1": nums[idx + 1], "C1": nums[idx + 2],
            "M11": nums[idx + 3], "B11": nums[idx + 4], "A11": nums[idx + 5],
            "C11": nums[idx + 6], "QX": nums[idx + 7],
        })
        idx += 8

    h0 = int(nums[idx]); idx += 1
    p["weirs"] = []
    for _ in range(h0):
        p["weirs"].append({
            "M2": nums[idx], "B2": nums[idx + 1],
            "C2": nums[idx + 2], "GMH": nums[idx + 3],
        })
        idx += 4

    # KK1, KK0, KK2 本版暂只支持 0（变宽变高/变宽/其它曲线）
    p["KK1"] = int(nums[idx]) if idx < len(nums) else 0; idx += 1
    if p["KK1"] > 0:
        idx += p["KK1"] * 6
    p["KK0"] = int(nums[idx]) if idx < len(nums) else 0; idx += 1
    if p["KK0"] > 0:
        idx += p["KK0"] * 5
    p["KK2"] = int(nums[idx]) if idx < len(nums) else 0; idx += 1
    if p["KK2"] > 0:
        idx += p["KK2"] * 2

    if p["KK1"] > 0 or p["KK0"] > 0 or p["KK2"] > 0:
        print("⚠ 变宽变高泄流孔/变宽溢洪道/其它泄流曲线本版未实现，按 0 处理")

    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params):
    """
    龙格-库塔法调洪演算。
    返回 (result_dict, 过程表)。
    """
    zps, fs = params["zps"], params["fs"]
    J, T = int(params["J"]), float(params["T"])
    QL = params["QL"]
    z0 = params["Z0"]
    # 面积单位：万平方米 → m²
    fs_m2 = [f * 1e4 for f in fs]

    def area(z):
        return interp1d(zps, fs_m2, z)

    qout, q_holes, q_weirs = make_outflow_func(params)

    def dzdz(t, z):
        """dz/dt = (Qin(t) - qout(z)) / F(z)
        水位下限约束：z <= Z0 且来水 <= 泄流能力时维持起调水位（出库=来水）。"""
        # Qin 线性插值（洪水过程按时段）
        i = min(int(t / T), J - 1)
        frac = (t / T) - i
        qin = QL[i] + frac * (QL[i + 1] - QL[i])
        fz = area(z)
        if fz <= 0:
            raise ValueError("水面面积必须大于 0")
        qout_cur = qout(z)
        # 水位下限约束：水位不降到起调水位以下
        if z <= z0 and qin <= qout_cur:
            return 0.0
        return (qin - qout_cur) / fz

    # RK4 积分
    ts, zs = [], []
    z = z0
    n_steps = J
    # 每时段 4 个子步（精度与速度平衡）
    sub = 4
    for n in range(n_steps):
        t0 = n * T
        for k in range(sub):
            z = rk4_step(dzdz, t0 + k * T / sub, z, T / sub)
            # 水位下限约束
            z = max(z, z0)
        ts.append((n + 1) * T)
        zs.append(z)

    # 汇总过程表
    table = []
    for n in range(J):
        z_cur = zs[n]
        qh = q_holes(z_cur)
        qw = q_weirs(z_cur)
        qt = qh + qw + params["W"]
        # 限泄
        for z_limit, q_limit in params.get("limits", []):
            if z_cur >= z_limit:
                qt = min(qt, q_limit)
        qin = QL[n + 1]  # 时段末来水
        # 水位下限约束：起调水位时来水即泄（出库=来水），洞/溢显示 -------
        if z_cur <= params["Z0"] + 1e-6 and qin <= qt:
            table.append({
                "时段": n + 1, "水位": round(z_cur, 2),
                "来水": round(qin, 2), "泄洪洞": None,
                "溢洪道": None, "下泄总量": round(qin, 2),
            })
            continue
        table.append({
            "时段": n + 1, "水位": round(z_cur, 2),
            "来水": round(qin, 2), "泄洪洞": round(qh, 2),
            "溢洪道": round(qw, 2), "下泄总量": round(qt, 2),
        })

    z_max = max(zs)
    z_max_t = (zs.index(z_max) + 1) * T
    return {
        "程序": PROGRAM_ID,
        "工程名": params.get("工程名", ""),
        "频率": params.get("频率", ""),
        "时段数": J,
        "时段间隔": T,
        "起始水位": z0,
        "最高水位": round(z_max, 2),
        "最高水位时刻": round(z_max_t, 0),
        "泄洪洞个数": len(params.get("holes", [])),
        "溢洪道个数": len(params.get("weirs", [])),
        "过程表": table,
    }, table


def render(params, result, table):
    """生成文本计算书（原著风格）。"""
    lines = []
    lines.append(f"工程名:{params.get('工程名', '')}  频率:{params.get('频率', '')}")
    lines.append("")
    lines.append("一. 原始数据:")
    lines.append(f"水位~水面面积关系曲线结点数 K= {params['K']}")
    zps, fs = params["zps"], params["fs"]
    half = (len(zps) + 1) // 2
    lines.append("水位(米) 水面面积(万平方米)   水位(米) 水面面积(万平方米)")
    for i in range(half):
        j = i + half
        if j < len(zps):
            lines.append(f"{zps[i]:>8.2f}  {fs[i]:>12.2f}   {zps[j]:>8.2f}  {fs[j]:>12.2f}")
        else:
            lines.append(f"{zps[i]:>8.2f}  {fs[i]:>12.2f}")
    lines.append(f"洪水过程时段数 J= {params['J']}   时段间隔 T= {params['T']}")
    lines.append("洪水过程(立方米/秒):")
    ql = params["QL"]
    for i in range(0, len(ql), 8):
        lines.append("  " + "  ".join(f"{v:>8.1f}" for v in ql[i:i + 8]))
    lines.append(f"防洪下限水位 ZM= {params['ZM']:.2f}  调洪起始水位 Z0= {params['Z0']:.2f}")
    lines.append(f"泄洪起始流量 Q0= {params['Q0']:.1f}  电站常流量 W= {params['W']:.1f}")
    lines.append("")

    lines.append("二. 计算结果:")
    lines.append("时段  水位     河道    泄洪洞   溢洪道     下泄")
    lines.append("            来水量    流量     流量     总流量")
    for row in table:
        qh_s = f"{row['泄洪洞']:>9.2f}" if row['泄洪洞'] is not None else "  -------"
        qw_s = f"{row['溢洪道']:>9.2f}" if row['溢洪道'] is not None else "  -------"
        lines.append(f"{row['时段']:>4} {row['水位']:>7.2f} {row['来水']:>9.2f} "
                     f"{qh_s} {qw_s} {row['下泄总量']:>9.2f}")
    lines.append("")
    lines.append(f"最高水位: {result['最高水位']:.2f} 米（时刻 {result['最高水位时刻']:.0f} 秒）")
    lines.append("")

    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """
    统一入口。
    data: INT 文件路径 | dict
    """
    params = parse(data)
    result, table = compute(params)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [("一", ["原始数据"]), ("二", ["结果见 JSON"])], result)
    else:
        text = render(params, result, table)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    res, txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(txt)
