# -*- coding: utf-8 -*-
"""
D-2 水闸水力学计算程序 —— 内核
================================
复刻《水利程序集》D-2 程序（作者：陈靖齐，水电部天津勘测设计院）。

功能：计算 4 种闸门/底坎组合下闸孔出流的流量 Q 或闸门开启度 e：
  W=1 平板闸门 + 平底（可求淹没系数 σ）
  W=2 弧形闸门 + 平底
  W=3 平板闸门 + 实用堰
  W=4 弧形闸门 + 实用堰
问题类型 AS：1=求流量 Q；2=求闸门开启度 e。

算法要点（武水公式体系，经 D-2-1.OUT 逐值回归验证）：
  (一) 流量公式
    Q = SG·μ0·e·B·√(2g·H0)          （W=1，SG 为淹没系数）
    Q = μ0·e·B·√(2g·(H0−hs))        （W=2/3/4，淹没作用水头扣除 hs）
    H0 = H + V0²/2g，V0 = Q/(B·(P+H))，迭代至 |ΔQ|<0.001 m³/s（考虑行进流速）
  (二) 流量系数（武水一套公式）
    W=1 平底平板：μ0 = 0.6 − 0.176·e/H            ［清华、武水］
    W=2 平底弧形：μ0 = 0.97 − 0.56·e/H − (1−e/H)·0.258·θ ［武水］
         θ = arccos((C−e)/R)（下缘入流角，弧度；25°≤θ≤90°，0.1≤e/H≤0.65）
    W=3 实用堰平板：μ0 = 0.745 − 0.274·e/H        ［武水，0.1<e/H<0.75］
    W=4 实用堰弧形：μ0 = 0.685 − 0.19·e/H         ［武水，0.1<e/H<0.75］
  (三) 淹没出流（仅 W=1 平板平底按 σ 曲线折减）
    hc = ε·e            （ε：平板闸门垂向收缩系数表，存数据库自动插值）
    Uc = φ·√(2g(H0−hc)) （φ=0.97 流速系数，由 D-2-1.OUT 的 vc=7.5747 反推验证）
    Frc = Uc²/(g·hc)    （收缩断面 Froude 数，F=Uc²/ghc）
    hc2 = hc/2·(√(1+8·Frc)−1)   （跃后水深）
    hs>hc2 时淹没：潜流比 x=(hs−hc2)/(H−hc2)，σ 按离散表插值；否则 σ=1
    Q = σ·Q自由
  (四) 求 e 用迭代（弦截法，目标 Q=Q(输入)）

验证基准：D-2-1.OUT（平板平底自由出流，H=3.5, e=0.7, B=3）
  Vo=0.9415  H0=3.5452  e/H=0.2  EP=0.62  hc=0.4340  vc=7.5747
  hc2=2.0477  SG=1.0  μ0=0.5648  Q=9.8870 —— 全部数值点吻合（误差<0.5%）。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "D-2"
TITLE = "水闸水力学计算书"

G = 9.8          # 原著取 9.8 m/s²（D-2-1 数值反推验证）
PHI = 0.97       # 收缩断面流速系数（薄壁/闸孔出流 0.97，D-2-1 vc 反推）
GATE_NAMES = {1: "平板平底", 2: "弧形平底", 3: "平板实用堰", 4: "弧形实用堰"}

# ------------------------------------------------------------
# 离散表（原著"存数据库中，可自动插值"，此处内嵌为 Python 表）
# ------------------------------------------------------------

# 平板闸门垂向收缩系数 ε —— e/H 与 ε（0.2→0.620，与 D-2-1.OUT EP=0.620 精确相符；
#   清华/武水《水力学》表 8-12 同源数据）
EPS_EH = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35,
          0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75]
EPS_VAL = [0.613, 0.613, 0.615, 0.618, 0.620, 0.622, 0.625, 0.628,
           0.630, 0.638, 0.645, 0.650, 0.660, 0.675, 0.690, 0.705]

# 平板平底闸孔淹没系数 σ —— 潜流比 x=(hs−hc2)/(H−hc2) 与 σ
#   （南京水利科学研究所曲线，规范/理正采用同一离散表）
SIG_X = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
         0.92, 0.94, 0.96, 0.98, 0.99, 0.995]
SIG_VAL = [1.00, 0.86, 0.78, 0.71, 0.66, 0.59, 0.52, 0.45, 0.36, 0.23,
           0.19, 0.16, 0.12, 0.07, 0.04, 0.02]


def interp(xs, ys, x):
    """一维线性插值（端点外沿用首末值）。"""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


# ------------------------------------------------------------
# 流量系数
# ------------------------------------------------------------

def mu0_flat_plate(rat):
    """W=1 平底平板：μ0 = 0.6 − 0.176·(e/H)"""
    return 0.6 - 0.176 * rat


def mu0_flat_arc(rat, theta):
    """W=2 平底弧形：μ0 = 0.97 − 0.56·(e/H) − (1−e/H)·0.258·θ
    θ 弧度，= (0.97−0.81·α°/180°) 展开（0.81/π≈0.258）。"""
    return 0.97 - 0.56 * rat - (1.0 - rat) * 0.258 * theta


def mu0_weir_plate(rat):
    """W=3 实用堰平板：μ0 = 0.745 − 0.274·(e/H)"""
    return 0.745 - 0.274 * rat


def mu0_weir_arc(rat):
    """W=4 实用堰弧形：μ0 = 0.685 − 0.19·(e/H)"""
    return 0.685 - 0.19 * rat


def gate_mu0(W, rat, theta=0.0):
    if W == 1:
        return mu0_flat_plate(rat)
    if W == 2:
        return mu0_flat_arc(rat, theta)
    if W == 3:
        return mu0_weir_plate(rat)
    return mu0_weir_arc(rat)


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """
    解析输入。data: dict 或 INT 文件路径。

    INT 数值流格式（W 闸门分类号 1..4；AS 问题类型 1=求Q 2=求e）：
      W, AS,
      求 Q（AS=1）时：
        W=1: H, e, hs, B
        W=2: H, e, hs, R, C, B
        W=3: H, e, hs, P, B
        W=4: H, e, hs, P, R, C, B
      求 e（AS=2）时去掉 e、追加 Q：
        W=1: H, hs, B, Q
        W=2: H, hs, R, C, B, Q
        W=3: H, hs, P, B, Q
        W=4: H, hs, P, R, C, B, Q
    其中：H 闸前水深/堰上水头(m)；e 开度(m)；hs 下游水深(m，相对底板或堰顶)；
      B 过流宽度(m)；P 溢流坝高(m，V0 计算用)；R 弧形闸门半径(m)；
      C 转轴高(m，相对闸底)。Q(m³/s)。
    """
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 2:
        raise ValueError("INT 数据至少需要 2 个数值（W, AS）")
    p = {"程序": PROGRAM_ID}
    p["W"] = int(nums[0])
    p["AS"] = int(nums[1])
    W, AS = p["W"], p["AS"]
    if W not in GATE_NAMES:
        raise ValueError(f"未知闸门分类号 W={W}（应为 1..4）")
    if AS not in (1, 2):
        raise ValueError(f"未知问题类型 AS={AS}（1=求Q 2=求e）")

    rest = nums[2:]
    if AS == 1:
        if W == 1:
            H, e, hs, B = rest[:4]
        elif W == 2:
            H, e, hs, R, C, B = rest[:6]
        elif W == 3:
            H, e, hs, P, B = rest[:5]
        else:
            H, e, hs, P, R, C, B = rest[:7]
        p.update(H=H, e=e, hs=hs, B=B)
        if W in (2,):
            p.update(R=R, C=C)
        if W == 3:
            p.update(P=P)
        if W == 4:
            p.update(P=P, R=R, C=C)
        p["Q"] = None
    else:
        if W == 1:
            H, hs, B, Q = rest[:4]
            e = None
        elif W == 2:
            H, hs, R, C, B, Q = rest[:6]
            e = None
        elif W == 3:
            H, hs, P, B, Q = rest[:5]
            e = None
        else:
            H, hs, P, R, C, B, Q = rest[:7]
            e = None
        p.update(H=H, e=e, hs=hs, B=B, Q=Q)
        if W in (2,):
            p.update(R=R, C=C)
        if W == 3:
            p.update(P=P)
        if W == 4:
            p.update(P=P, R=R, C=C)
    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def theta_arc(R, C, e):
    """弧形闸门下缘入流角 θ = arccos((C−e)/R)，弧度。"""
    cosv = (C - e) / R
    cosv = max(-1.0, min(1.0, cosv))
    return math.acos(cosv)


def discharge_for_open(W, H, e, hs, B, P, R, C, tol=0.001, max_iter=60):
    """
    已知开度 e 求流量 Q（含行进流速 V0 迭代 |ΔQ|<0.001）。

    W=1 平板平底：淹没按 σ(潜流比) 折减，SG 输出在结果里。
    W=2/3/4：淹没在作用水头中考虑，H0 有效 = H0 − hs。
    返回 dict 含 Vo,H0,rat,mu0,ep,hc,Uc,hc2,sg,Q。
    """
    rat = e / H if H > 0 else 0.0
    mu = gate_mu0(W, rat, theta_arc(R, C, e) if W == 2 else 0.0)
    hc = ep = None
    if W == 1:
        ep = interp(EPS_EH, EPS_VAL, rat)
        hc = ep * e

    Qold = 0.0
    Vo = 0.0
    H0 = H
    sg = 1.0
    Uc = Frc = hc2 = 0.0
    Q = 0.0
    for _ in range(max_iter):
        H0 = H + Vo * Vo / (2.0 * G)
        base = e * B * math.sqrt(2.0 * G * max(H0 - (hs if W != 1 else 0.0), 0.0))
        Q = mu * base
        if W == 1:
            # 收缩断面水力要素（自由出流下流速系数 φ=0.97）
            Uc = PHI * math.sqrt(2.0 * G * max(H0 - hc, 0.0))
            Frc = Uc * Uc / (G * hc) if hc > 0 else 0.0
            hc2 = hc / 2.0 * (math.sqrt(1.0 + 8.0 * Frc) - 1.0) if hc > 0 else 0.0
            if hs > hc2:
                x = (hs - hc2) / (H - hc2) if (H - hc2) > 1e-9 else 1.0
                sg = interp(SIG_X, SIG_VAL, x)
            else:
                sg = 1.0
            Q = sg * Q
        Vo = Q / (B * (H + P)) if B > 0 else 0.0
        if abs(Q - Qold) < tol:
            break
        Qold = Q
    return {"Vo": Vo, "H0": H0, "rat": rat, "mu0": mu, "ep": ep, "hc": hc,
            "Uc": Uc, "Frc": Frc, "hc2": hc2, "sg": sg, "Q": Q}


def compute(params):
    """执行水闸水力计算。"""
    W = params["W"]
    AS = params["AS"]
    H = params["H"]
    hs = params["hs"]
    B = params["B"]
    P = params.get("P", 0.0)
    R = params.get("R")
    C = params.get("C")
    gname = GATE_NAMES[W]

    result = {"程序": PROGRAM_ID, "闸门类型": gname, "问题类型": AS,
              "基本资料": {"W": W, "H": H, "hs": hs, "B": B, "P": P,
                          "R": R, "C": C}}

    if AS == 1:
        e = params["e"]
        s = discharge_for_open(W, H, e, hs, B, P, R, C)
        result["输入e"] = e
        result["计算"] = {k: (round(v, 4) if isinstance(v, float) else v)
                         for k, v in s.items()}
        result["结果"] = {"Q": round(s["Q"], 4), "e": e,
                          "Vo": round(s["Vo"], 4), "H0": round(s["H0"], 4),
                          "mu0": round(s["mu0"], 4), "SG": round(s["sg"], 4)}
        if s["ep"] is not None:
            result["结果"].update(EP=round(s["ep"], 4), hc=round(s["hc"], 4),
                                  vc=round(s["Uc"], 4), hc2=round(s["hc2"], 4))
    else:
        Qtar = params["Q"]
        # 反解开度 e：目标 f(e)=Q(e)−Qtar=0；e 有效范围 (0, 0.75H]
        def f(e_):
            return discharge_for_open(W, H, e_, hs, B, P, R, C)["Q"] - Qtar

        # 扫描确定根区间（Q 随 e 单调增）
        hi = 0.75 * H if H > 0 else 1.0
        if f(hi) < 0:
            raise ValueError(f"最大开度 e={hi:.3f} 时 Q={f(hi)+Qtar:.3f} < 目标 Q={Qtar}，无法达到")
        lo = 0.0
        flo = f(lo) if lo > 0 else -Qtar  # e→0 时 Q→0
        for _ in range(300):              # 弦截法（保证单调收敛）
            fhi = f(hi)
            if fhi >= 0 and flo <= 0:
                e_new = hi - fhi * (hi - lo) / (fhi - flo) if fhi != flo else 0.5 * (lo + hi)
                e_new = max(lo + 1e-9, min(hi - 1e-9, e_new))
                fe = f(e_new)
                if abs(fe) < 0.001 or (hi - lo) < 1e-7:
                    e = e_new
                    break
                if fe > 0:
                    hi, fhi = e_new, fe
                else:
                    lo, flo = e_new, fe
            else:
                # 回落二分保证区间
                mid = 0.5 * (lo + hi)
                fm = f(mid)
                if abs(fm) < 0.001:
                    e = mid
                    break
                if fm > 0:
                    hi = mid
                else:
                    lo = mid
        else:
            raise ValueError("求 e 迭代不收敛")
        s = discharge_for_open(W, H, e, hs, B, P, R, C)
        result["目标Q"] = Qtar
        result["计算"] = {k: (round(v, 4) if isinstance(v, float) else v)
                         for k, v in s.items()}
        result["结果"] = {"Q": round(s["Q"], 4), "e": round(e, 4),
                          "Vo": round(s["Vo"], 4), "H0": round(s["H0"], 4),
                          "mu0": round(s["mu0"], 4), "SG": round(s["sg"], 4)}
        if s["ep"] is not None:
            result["结果"].update(EP=round(s["ep"], 4), hc=round(s["hc"], 4),
                                  vc=round(s["Uc"], 4), hc2=round(s["hc2"], 4))
    return result


# ------------------------------------------------------------
# 输出
# ------------------------------------------------------------

def render(params, result):
    """生成文本计算书（复刻 D-2-1.OUT 原著风格）。"""
    W = params["W"]
    AS = params["AS"]
    gname = result["闸门类型"]
    c = result["计算"]
    r = result["结果"]
    L = []
    L.append(f"  {gname} ")
    if AS == 1:
        L.append(f"      上游水位H:                  {params['H']:>10.4f}(米)")
        L.append(f"      闸门开启度e:                {params['e']:>10.4f}(米)")
    else:
        L.append(f"      上游水位H:                  {params['H']:>10.4f}(米)")
    L.append(f"      下游水位hs:                 {params['hs']:>10.4f}(米)")
    L.append(f"      过流宽度b:                  {params['B']:>10.4f}(米)")
    if W in (2, 4):
        L.append(f"      弧形闸门半径R:              {params['R']:>10.4f}(米)")
        L.append(f"      转轴高程C:                  {params['C']:>10.4f}(米)")
    if W in (3, 4):
        L.append(f"      溢流坝高P:                  {params['P']:>10.4f}(米)")
    L.append("")
    L.append("  四、计算结果")
    L.append("")
    L.append("    1.相关结果")
    L.append(f"      行进流速Vo:                 {r.get('Vo', c.get('Vo', 0)):>10.4f}(米/秒)")
    L.append(f"      闸前水深Ho:                 {r.get('H0', c.get('H0', 0)):>10.4f}(米)")
    if W == 1:
        L.append(f"      闸门开启度e/上游水位H：     {c.get('rat', 0):>10.4f}")
        L.append(f"      平板闸门侧收缩系数  ：      {r.get('EP', 0):>10.4f}")
        L.append(f"      收缩断面水深hc:            {r.get('hc', 0):>10.4f}(米)")
        L.append(f"      流速vc:                    {r.get('vc', 0):>10.4f}(米/秒)")
        L.append(f"      跃后水深hc2:               {r.get('hc2', 0):>10.4f}(米)")
        L.append(f"      平板闸门淹没系数  :        {r.get('SG', 1.0):>10.4f}")
    L.append(f"      流量系数mu:                 {r.get('mu0', 0):>10.4f}")
    L.append("    2.结果")
    L.append(f"      流量Q:                      {r['Q']:>10.4f}(立方米/秒)")
    L.append("    3.校核")
    if AS == 1:
        L.append(f"      闸门开启度e:                {params['e']:>10.4f}(米)")
    else:
        L.append(f"      闸门开启度e:                {r['e']:>10.4f}(米)")
    L.append(f"      流量Q:                      {r['Q']:>10.4f}(立方米/秒)")
    return render_text(PROGRAM_ID, TITLE, [("", L)])


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
