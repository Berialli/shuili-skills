# -*- coding: utf-8 -*-
"""
D-6 常用断面渠道水力学计算程序 —— 内核
========================================
复刻《水利程序集》D-6 程序（作者：杨志河，水电部天津勘测设计院）。

支持梯形断面（含矩形/三角形）与复式断面渠道水力计算，5 种计算功能：
  1. 过水能力校核（已知 B,H,M1,M2,N,I → Q,V,Bo）
  2. 求均匀流水深（已知 B,Q,M1,M2,N,I → H）
  3. 求底坡（已知 B,H,Q,M1,M2,N → I）
  4. 求底宽（已知 H,Q,M1,M2,N,I → B）
  5. 水力最优断面（已知 Q,M,N,I → Bm,Hm）

核心公式（经原始 OUT 回归验证）：
  面积 W = B*H + (M1+M2)*H²/2
  湿周 X = B + H*(√(1+M1²) + √(1+M2²))
  水力半径 R = W/X
  谢才系数 C = (1/N)*R^Y
    Y=0（INT 默认）→ 巴甫洛夫斯基：Y = 2.5√N - 0.13 - 0.75√R(√N-0.10)
    Y=1/6 → 曼宁公式
  流量 Q = W*C*√(R*I)  流速 V = C*√(R*I)
  水面宽 Bo = B + (M1+M2)*H
  水力最优宽深比 Pm = Bm/Hm = 2*(√(1+M²)-M)

验证基准：D-61-1.INT（B=7,H=1.4403,M=1.5,N=0.025,I=0.0003）
  Q=9.6799  V=0.7337  Bo=11.3209
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "D-6"
TITLE = "常用断面灌溉渠道水力学计算书"

G = 9.81


def pavlovsky_y(N, R):
    """巴甫洛夫斯基公式指数 Y = 2.5√N - 0.13 - 0.75√R(√N-0.10)"""
    return 2.5 * math.sqrt(N) - 0.13 - 0.75 * math.sqrt(R) * (math.sqrt(N) - 0.10)


def chezy_C(N, R, Y):
    """谢才系数 C = (1/N)*R^Y"""
    return (1.0 / N) * (R ** Y)


def trapezoid_elements(B, H, M1, M2, N, I, Y=None):
    """
    梯形断面水力要素。
    Y=None → 巴甫洛夫斯基自动；Y=1/6 → 曼宁。
    返回 dict(W,X,R,C,Q,V,Bo,Y)
    """
    W = B * H + (M1 + M2) * H * H / 2.0
    X = B + H * (math.sqrt(1 + M1 * M1) + math.sqrt(1 + M2 * M2))
    R = W / X if X > 0 else 0.0
    if Y is None:
        Y = pavlovsky_y(N, R)
    C = chezy_C(N, R, Y)
    V = C * math.sqrt(R * I) if R > 0 else 0.0
    Q = W * V
    Bo = B + (M1 + M2) * H
    return {"W": W, "X": X, "R": R, "C": C, "Q": Q, "V": V, "Bo": Bo, "Y": Y}


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def parse(data):
    """
    解析输入。
    data: dict 或 INT 文件路径。
    INT 格式：断面类型, 功能号, B, Q, M1, M2, N, Y, I  （梯形）
    复式断面另行处理。
    """
    if isinstance(data, dict):
        return data

    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}
    p["断面类型"] = int(nums[0])
    p["功能"] = int(nums[1])
    if p["断面类型"] == 1:  # 梯形
        # 字段顺序按功能不同：见说明文档
        # 功能 1: 求底宽  输入 H, Q
        # 功能 2: 求流量  输入 B, H
        # 功能 3: 求水深  输入 B, Q
        # 功能 4: 求底坡  输入 B, H, Q（无 I，I 为输出）
        func = p["功能"]
        # 尾部字段数：功能4无I，其余功能有I
        tail_len = 4 if func == 4 else 5  # M1,M2,N,Y[,I]
        tail = nums[-tail_len:]
        p["M1"] = tail[0]
        p["M2"] = tail[1]
        p["N"] = tail[2]
        p["Y"] = tail[3]  # 0=自动（巴甫洛夫斯基）
        p["I"] = tail[4] if tail_len == 5 else 0.0  # 功能4时I待求，暂置0
        head = nums[2:-tail_len]  # 前面的参数
        if func == 1:
            p["H"] = head[0]
            p["Q"] = head[1]
            p["B"] = None  # 待求
        elif func == 2:
            p["B"] = head[0]
            p["H"] = head[1]
            p["Q"] = None  # 待求
        elif func == 3:
            p["B"] = head[0]
            p["Q"] = head[1]
            p["H"] = None  # 待求
        elif func == 4:
            p["B"] = head[0]
            p["H"] = head[1]
            p["Q"] = head[2]
        else:
            raise ValueError(f"未知功能号 {func}")
    else:
        raise NotImplementedError("复式断面（类型 2）将在后续版本实现")
    return p


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def compute(params):
    """执行渠道水力计算（按功能号）。"""
    stype = params.get("断面类型", 1)
    func = params.get("功能", 1)
    if stype != 1:
        raise NotImplementedError("复式断面暂未实现")

    B = params.get("B")
    Q = params.get("Q")
    M1, M2, N, I = params["M1"], params["M2"], params["N"], params["I"]
    Y_raw = params.get("Y", 0.0)
    Y_fixed = None if Y_raw == 0 else float(Y_raw)  # 0=自动

    result = {"程序": PROGRAM_ID, "断面类型": "梯形", "功能": func}
    result["基本资料"] = {
        "B": B, "Q": Q, "M1": M1, "M2": M2, "N": N, "I": I,
        "Y": "自动(巴甫洛夫斯基)" if Y_raw == 0 else Y_raw,
    }

    if func == 1:
        # 求底宽 B：给定 H, Q
        H = params["H"]
        def err(b):
            el = trapezoid_elements(b, H, M1, M2, N, I, Y_fixed)
            return el["Q"] - Q
        B_calc = _bisect_solve(err, 1e-3, 100.0)
        el = trapezoid_elements(B_calc, H, M1, M2, N, I, Y_fixed)
        result["计算"] = {"H": H, "B(试算)": round(B_calc, 4)}
        result["结果"] = {"B": round(B_calc, 4), "Q(校核)": round(el["Q"], 4),
                         "V": round(el["V"], 4), "Bo": round(el["Bo"], 4), "H": H}
    elif func == 2:
        # 求流量（过水能力校核）：给定 B, H
        B, H = params["B"], params["H"]
        el = trapezoid_elements(B, H, M1, M2, N, I, Y_fixed)
        result["计算"] = {"H": H, **{k: round(v, 4) for k, v in el.items()}}
        result["结果"] = {"Q": round(el["Q"], 4), "V": round(el["V"], 4),
                         "Bo": round(el["Bo"], 4), "R": round(el["R"], 4), "H": H}
    elif func == 3:
        # 求均匀流水深 H：给定 B, Q
        B, Q = params["B"], params["Q"]
        def err(h):
            el = trapezoid_elements(B, h, M1, M2, N, I, Y_fixed)
            return el["Q"] - Q
        H = _bisect_solve(err, 1e-4, 100.0)
        el = trapezoid_elements(B, H, M1, M2, N, I, Y_fixed)
        result["计算"] = {"H(试算)": round(H, 4), **{k: round(v, 4) for k, v in el.items()}}
        result["结果"] = {"H": round(H, 4), "Q(校核)": round(el["Q"], 4),
                         "V": round(el["V"], 4), "Bo": round(el["Bo"], 4)}
    elif func == 4:
        # 求底坡 I：给定 B, H, Q
        B, H = params["B"], params["H"]
        el = trapezoid_elements(B, H, M1, M2, N, I, Y_fixed)
        I_calc = (Q / (el["W"] * el["C"])) ** 2 / el["R"]
        el2 = trapezoid_elements(B, H, M1, M2, N, I_calc, Y_fixed)
        result["计算"] = {"H": H, "W": round(el["W"], 4), "C": round(el["C"], 4), "R": round(el["R"], 4)}
        result["结果"] = {"I": I_calc, "Q(校核)": round(el2["Q"], 4), "V": round(el2["V"], 4), "H": H}
    elif func == 5:
        # 水力最优断面
        m = M1  # 对称
        Pm = 2 * (math.sqrt(1 + m * m) - m)
        def err(h):
            b = Pm * h
            el = trapezoid_elements(b, h, m, m, N, I, Y_fixed)
            return el["Q"] - Q
        Hm = _bisect_solve(err, 1e-4, 100.0)
        Bm = Pm * Hm
        el = trapezoid_elements(Bm, Hm, m, m, N, I, Y_fixed)
        result["计算"] = {"宽深比 Pm": round(Pm, 4)}
        result["结果"] = {"Bm": round(Bm, 4), "Hm": round(Hm, 4),
                         "Q(校核)": round(el["Q"], 4), "V": round(el["V"], 4)}
    return result


def _bisect_solve(f, lo, hi, tol=1e-8, max_iter=200):
    """二分法求 f(x)=0（f 单调）。"""
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise ValueError(f"二分法初始区间 f(lo)={flo:.3g}, f(hi)={fhi:.3g} 同号")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol or (hi - lo) < tol:
            return mid
        if flo * fm < 0:
            hi = mid
        else:
            lo = mid
            flo = fm
    return 0.5 * (lo + hi)


def render(params, result):
    """生成文本计算书（原著风格）。"""
    lines = []
    lines.append("断面形式:                      梯形")
    r = result["结果"]
    lines.append(f"梯形断面渠道的底宽B:       {r.get('B', params['B']):>10.4f}(米)")
    lines.append(f"均匀流水深H:               {r.get('H', params.get('H', 0)):>10.4f}(米)")
    lines.append(f"边坡系数 M1:               {params['M1']:>10.4f}")
    lines.append(f"边坡系数 M2:               {params['M2']:>10.4f}")
    lines.append(f"糙率系数N:                 {params['N']:>10.4f}")
    lines.append(f"渠道底坡I:                 {params['I']:>10.6f}")
    lines.append("")
    lines.append("四、计算结果")
    if "H" in r:
        lines.append(f"断面水深(试算值)H:       {r['H']:>10.4f}(米)")
    if "Q" in r:
        lines.append(f"过水断面流量Q:           {r['Q']:>10.4f}(立方米/秒)")
    if "V" in r:
        lines.append(f"断面平均流速V:           {r['V']:>10.4f}(米/秒)")
    if "Bo" in r:
        lines.append(f"水面宽度B0:              {r['Bo']:>10.4f}(米)")
    if "I" in r:
        lines.append(f"渠道底坡I:               {r['I']:>10.6f}")
    if "B" in r:
        lines.append(f"渠道底宽B:               {r['B']:>10.4f}(米)")
    if "Bm" in r:
        lines.append(f"水力最优底宽Bm:          {r['Bm']:>10.4f}(米)")
        lines.append(f"水力最优水深Hm:          {r['Hm']:>10.4f}(米)")
    lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """
    统一入口。
    data: INT 文件路径 | dict
    """
    params = parse(data)
    # 功能 1 需要 H，从 params 中取（dict 输入时用户提供）
    result = compute(params)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE, [("一", ["基本资料"]), ("二", ["结果见 JSON"])], result)
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
