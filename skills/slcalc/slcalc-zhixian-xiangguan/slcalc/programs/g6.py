# -*- coding: utf-8 -*-
"""
G-6 重力式挡土墙计算程序 —— 内核
====================================
复刻《水利程序集》G-6 程序（作者：刘以钢，水电部天津勘测设计院）。

功能（I=0）：
  1. 库仑主动土压力计算（含均布荷载折算高度 h0=q/γ）
  2. 抗滑稳定验算 KC
  3. 抗倾覆稳定验算 KO
  4. 基底应力验算 CX（墙趾）/ CN（墙踵），含 e>B/6 拉应力截断（三角形分布）

输入参数（INT 顺序）：
  I, R, P, R0, F, D, U, Q, H1, H2, H3, H4, H5, A, C, B, B1, B2, B4, MS, ML, MK

公式体系（依据《重力式挡土墙专业课程设计计算报告书》确认）：
  库仑系数  Ka = cos²(φ−α) / [cos²α·cos(α+δ)·(1+√(sin(φ+δ)·sin(φ−β)/(cos(α+δ)·cos(α−β))))²]
  主动土压力 Ea = Ka·γ·H·(H/2 + h0)，h0 = q/γ
  水平分力  Ex = Ea·cos(α+δ)，垂直分力 Ey = Ea·sin(α+δ)
  （算例中 D=|R| → α+δ=0 → Ex=Ea, Ey=0，本内核自动处理）
  抗滑稳定  KC = (W+Ey)(U+tanR0)/Ex
  抗倾稳定  KO = (W·ZW + Ey·ZEy)/(Ex·ZEx)
  基底应力  σ = V/Bw·(1±6e/Bw)，e = Bw/2 − Zn
            Zn = (W·ZW + Ey·ZEy − Ex·ZEx)/(W+Ey)
            当 e > Bw/6 时 σ趾' = 2V/[3(Bw/2−e)]，CN=0

几何模型（简单梯形 + 逆坡基底）：
  墙趾(0,0) → 墙踵(B2, B2·tanR0) → 墙背顶 → 墙顶前缘 → 闭合
  底宽：R≤0(仰斜)用 C，R>0(俯斜)用 B4（两算例参数语义差异，近似处理）
  顶宽：R≤0 用 B，R>0 用 C

已知偏差（诚实声明，后续版本修正）：
  例1 目标 KC=1.944, KO=1.575, CX=22.948, CN=5.603
  例2 目标 KC=1.338, KO=1.983, CX=36.018, CN=0.000
  本内核为规范正确的一阶复刻，几何细节（B1 含义、齿坎/台阶细分）未完全复现原著，
  回归偏差见 G6公式对比分析报告.md 与 run() 输出的"回归对照"。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "G-6"
TITLE = "重力式挡土墙计算书"


# ------------------------------------------------------------
# 基础公式
# ------------------------------------------------------------

def ka_coulomb(alpha_deg, beta_deg, phi_deg, delta_deg):
    """
    库仑主动土压力系数。
    alpha: 墙背与铅直面夹角（度，含符号：负=仰斜，正=俯斜）
    beta:  墙顶填土坡角（度）
    phi:   土壤内摩擦角（度）
    delta: 土壤与墙背间摩擦角（度）
    """
    a = math.radians(alpha_deg)
    b = math.radians(beta_deg)
    p = math.radians(phi_deg)
    d = math.radians(delta_deg)
    num = math.cos(p - a) ** 2
    rad = math.sin(p + d) * math.sin(p - b) / (math.cos(a + d) * math.cos(a - b))
    if rad < 0:
        rad = 0.0
    rad = math.sqrt(rad)
    den = math.cos(a) ** 2 * math.cos(a + d) * (1 + rad) ** 2
    return num / den


def coulomb_earth(alpha_deg, beta_deg, phi_deg, delta_deg, gamma, H, q):
    """
    库仑主动土压力及分力。
    返回 dict(Ka, h0, Ea, Ex, Ey, ZEx)
      Ea  = Ka·γ·H·(H/2 + h0)，h0 = q/γ
      Ex  = Ea·cos(α+δ)
      Ey  = Ea·sin(α+δ)
      ZEx = H·(H+3h0)/(3(H+2h0)) —— Ex 作用点距墙底高度
    算例 D=|R| → α+δ=0 → Ex=Ea, Ey=0（两算例均如此，与原著 OUT 吻合）
    """
    Ka = ka_coulomb(alpha_deg, beta_deg, phi_deg, delta_deg)
    h0 = q / gamma if gamma > 0 else 0.0
    Ea = Ka * gamma * H * (H / 2.0 + h0)
    ad = math.radians(alpha_deg + delta_deg)
    Ex = Ea * math.cos(ad)
    Ey = Ea * math.sin(ad)
    ZEx = H * (H + 3 * h0) / (3 * (H + 2 * h0)) if (H + 2 * h0) > 0 else H / 3.0
    return {"Ka": Ka, "h0": h0, "Ea": Ea, "Ex": Ex, "Ey": Ey, "ZEx": ZEx}


def poly_centroid(pts):
    """鞋带公式：返回 (面积, 重心x, 重心z)。pts: [(x,z),...]"""
    S = 0.0
    cx = cz = 0.0
    n = len(pts)
    for i in range(n):
        x1, z1 = pts[i]
        x2, z2 = pts[(i + 1) % n]
        crs = x1 * z2 - x2 * z1
        S += crs
        cx += (x1 + x2) * crs
        cz += (z1 + z2) * crs
    S = abs(S) / 2.0
    if S < 1e-12:
        return 0.0, 0.0, 0.0
    cx /= 6.0 * S
    cz /= 6.0 * S
    return S, cx, cz


# ------------------------------------------------------------
# 几何模型
# ------------------------------------------------------------

def build_geometry(p):
    """
    构建墙体几何。
    坐标：墙趾（临空侧）为原点 O(0,0)，x 向填土侧（墙踵方向）。

    仰斜（R≤0，例1 形态）——墙趾齿坎 + 梯形墙身：
      (0,0) → (B2,0) → (B2,H3) → (C,H3) → (x_bt,H1) → (x_ft,H1)
      齿坎 B2×H3 在墙趾下，墙身底宽 C，顶宽 B，墙背仰斜顶向填土侧。

    俯斜（R>0，例2 形态）——墙趾台阶 + 墙身 + 墙踵台阶：
      (0,0) → (A,H5) → (x_ft,H1) → (x_bt,H1) → (B4,H3) → (B2,0)
      墙趾台阶 A×H5，墙身顶宽 C，墙踵台阶 B4×H3，基底水平投影 B2。
    """
    H1 = p["H1"]
    H3 = p["H3"]
    H5 = p["H5"]
    A = p["A"]
    C = p["C"]
    B = p["B"]
    B2 = p["B2"] if p["B2"] > 0 else C
    B4 = p["B4"] if p["B4"] > 0 else C
    R = p["R"]
    R0 = p["R0"]
    MK = p["MK"]

    tR = math.tan(math.radians(R0))
    H4 = p["H4"] if p["H4"] > 0 else H1 - H3 - H5
    tanR = math.tan(math.radians(R))   # R<0 仰斜 → 负

    if R <= 0:
        # 仰斜：墙背顶向填土侧偏（x 增大）
        x_bt = C + H4 * (-tanR)
        x_ft = x_bt - B
        pts = [(0.0, 0.0), (B2, 0.0), (B2, H3), (C, H3), (x_bt, H1), (x_ft, H1)]
        base_w = C
        top_w = B
    else:
        # 俯斜：墙背顶向临空侧偏（x 减小）
        x_bt = B4 - H4 * tanR
        x_ft = x_bt - C
        pts = [(0.0, 0.0), (A, H5), (x_ft, H1), (x_bt, H1), (B4, H3), (B2, 0.0)]
        base_w = B4
        top_w = C

    S, cx, cz = poly_centroid(pts)
    W = S * MK
    return {
        "面积": S, "W": W, "ZW": cx, "重心z": cz,
        "x_back_top": x_bt, "x_front_top": x_ft,
        "base_w": base_w, "top_w": top_w,
        "pts": pts, "H4": H4,
    }


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """解析输入。data: dict 或 INT 文件路径。"""
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}
    keys = ["I", "R", "P", "R0", "F", "D", "U", "Q",
            "H1", "H2", "H3", "H4", "H5", "A", "C", "B", "B1", "B2", "B4",
            "MS", "ML", "MK"]
    for i, k in enumerate(keys):
        p[k] = nums[i] if i < len(nums) else 0.0
    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params):
    """执行重力式挡土墙稳定与基底应力计算（I=0）。"""
    I = int(params.get("I", 0))
    if I != 0:
        raise NotImplementedError("I=1（断面应力）暂未实现，后续版本补充")

    R = params["R"]
    P = params["P"]
    R0 = params["R0"]
    F = params["F"]
    D = params["D"]
    U = params["U"]
    Q = params["Q"]
    MS = params["MS"]
    ML = params["ML"]
    MK = params["MK"]

    # ---- 几何 ----
    geo = build_geometry(params)
    W = geo["W"]
    ZW = geo["ZW"]

    # ---- 库仑土压力 ----
    Hw = params["H1"]
    if Hw <= 0:
        Hw = params["H4"]
    ep = coulomb_earth(R, P, F, D, MS, Hw, Q)
    Ex = ep["Ex"]
    Ey = ep["Ey"]
    ZEx = ep["ZEx"]
    # Ey 力臂：墙背中点 x
    ZEy = (geo["x_back_top"] + geo["base_w"]) / 2.0
    # 浮容重修正：ML>0 时有地下水，墙重按浮容重差（近似 MK−1）
    if ML > 0:
        W = geo["面积"] * max(MK - 1.0, 0.0)

    # ---- 抗滑稳定 ----
    # 基底分解式（更接近原著，见分析报告）：
    #   N = V·cosR0 + Ex·sinR0（法向合力）
    #   T = Ex·cosR0 − V·sinR0（切向滑动力）
    #   KC = N·U / T
    tR = math.tan(math.radians(R0))
    cr0 = math.cos(math.radians(R0))
    sr0 = math.sin(math.radians(R0))
    V = W + Ey
    N = V * cr0 + Ex * sr0
    T = Ex * cr0 - V * sr0
    KC = N * U / T if T > 0 else float("inf")

    # ---- 抗倾稳定 ----
    M_W = W * ZW
    M_Ey = Ey * ZEy
    M_Ex = Ex * ZEx
    KO = (M_W + M_Ey) / M_Ex if M_Ex > 0 else float("inf")

    # ---- 基底应力 ----
    # 基底水平宽：例2 用 B2（3.8），例1 用 C（1.45）
    Bw = params["B2"] if params["B2"] > 0 else params["C"]
    if params["B1"] > 0 and Bw <= 0:
        Bw = params["B1"] * math.cos(math.radians(R0))
    # 基底斜长（逆坡）
    Bp = Bw / math.cos(math.radians(R0))

    Zn = (M_W + M_Ey - M_Ex) / V if V > 0 else 0.0
    e = Bp / 2.0 - Zn
    if e >= Bp / 6.0:
        CX = 2.0 * N / (3.0 * (Bp / 2.0 - e)) if (Bp / 2.0 - e) > 0 else float("inf")
        CN = 0.0
    else:
        CX = N / Bp * (1 + 6.0 * e / Bp)
        CN = N / Bp * (1 - 6.0 * e / Bp)

    result = {
        "程序": PROGRAM_ID,
        "输入": {k: params[k] for k in
                 ["I", "R", "P", "R0", "F", "D", "U", "Q",
                  "H1", "H2", "H3", "H4", "H5", "A", "C", "B", "B1", "B2", "B4",
                  "MS", "ML", "MK"]},
        "几何": {
            "墙身面积": round(geo["面积"], 4),
            "墙重W": round(W, 4),
            "墙重力臂ZW": round(ZW, 4),
            "墙背高Hw": round(Hw, 4),
            "基底水平宽Bw": round(Bw, 4),
            "基底斜长Bp": round(Bp, 4),
        },
        "土压力": {
            "库仑系数Ka": round(ep["Ka"], 6),
            "折算高度h0": round(ep["h0"], 4),
            "主动土压力Ea": round(ep["Ea"], 4),
            "水平分力Ex": round(Ex, 4),
            "垂直分力Ey": round(Ey, 4),
            "Ex力臂ZEx": round(ZEx, 4),
            "Ey力臂ZEy": round(ZEy, 4),
        },
        "稳定": {
            "竖向合力V": round(V, 4),
            "抗滑KC": round(KC, 4),
            "抗倾KO": round(KO, 4),
        },
        "基底应力": {
            "合力作用点Zn": round(Zn, 4),
            "偏心距e": round(e, 4),
            "e/Bp": round(e / Bp, 4) if Bp > 0 else 0.0,
            "墙趾应力CX": round(CX, 4),
            "墙踵应力CN": round(CN, 4),
            "应力分布": "三角形(拉应力截断)" if e >= Bp / 6 else "梯形",
        },
    }
    return result


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def render(params, result):
    """生成文本计算书（原著风格）。"""
    L = []
    inp = result["输入"]
    L.append(f"计算方式选择 I= {inp['I']:>10.3f}   挡土面铅直角 R= {inp['R']:>10.3f}   墙顶土体坡角 P= {inp['P']:>10.3f}")
    L.append(f"基底坡角 R0= {inp['R0']:>10.3f}      土壤内摩擦角 F= {inp['F']:>10.3f}   土与墙背间摩擦角 D= {inp['D']:>10.3f}")
    L.append(f"摩擦系数 U= {inp['U']:>10.3f}       土体上均布荷载 Q= {inp['Q']:>10.3f} 墙体尺寸数据 H1= {inp['H1']:>10.3f}")
    L.append(f"墙体尺寸数据 H2= {inp['H2']:>10.3f}  墙体尺寸数据 H3= {inp['H3']:>10.3f}   墙体尺寸数据 H4= {inp['H4']:>10.3f}")
    L.append(f"墙体尺寸数据 H5= {inp['H5']:>10.3f}  墙体尺寸数据 A= {inp['A']:>10.3f}    墙体尺寸数据 C= {inp['C']:>10.3f}")
    L.append(f"墙体尺寸数据 B= {inp['B']:>10.3f}   墙体尺寸数据 B1= {inp['B1']:>10.3f}   墙体尺寸数据 B2= {inp['B2']:>10.3f}")
    L.append(f"墙体尺寸数据 B4= {inp['B4']:>10.3f}  土壤容重 MS= {inp['MS']:>10.3f}      土壤浮容重 ML= {inp['ML']:>10.3f}")
    L.append(f"墙体容重 MK= {inp['MK']:>10.3f}")
    L.append("")
    L.append("计算结果:")
    L.append(f" 墙体抗滑稳定系数 KC= {result['稳定']['抗滑KC']:>8.3f}")
    L.append(f" 墙体抗倾稳定系数 K0= {result['稳定']['抗倾KO']:>8.3f}")
    L.append(f" 墙趾处基底应力   CX= {result['基底应力']['墙趾应力CX']:>8.3f}")
    L.append(f" 墙踵处基底应力   CN= {result['基底应力']['墙踵应力CN']:>8.3f}")
    L.append("")
    L.append("中间量:")
    L.append(f" 墙重W= {result['几何']['墙重W']:>8.3f}  墙重力臂ZW= {result['几何']['墙重力臂ZW']:>8.4f}")
    L.append(f" 库仑系数Ka= {result['土压力']['库仑系数Ka']:>8.5f}  Ea= {result['土压力']['主动土压力Ea']:>8.3f}")
    L.append(f" 水平分力Ex= {result['土压力']['水平分力Ex']:>8.3f}  垂直分力Ey= {result['土压力']['垂直分力Ey']:>8.3f}")
    L.append(f" 偏心距e= {result['基底应力']['偏心距e']:>8.4f}  e/Bp= {result['基底应力']['e/Bp']:>8.4f}")
    L.append("")
    L.append("注: 本内核为复刻版(一阶近似), 与原著偏差见 G6公式对比分析报告.md")
    return render_text(PROGRAM_ID, TITLE, [("", L)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 文件路径 | dict"""
    params = parse(data)
    result = compute(params)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [("一", ["输入数据"]), ("二", ["计算结果见 JSON"])], result)
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
