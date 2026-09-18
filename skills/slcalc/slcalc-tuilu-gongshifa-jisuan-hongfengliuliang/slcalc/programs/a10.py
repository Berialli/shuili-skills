# -*- coding: utf-8 -*-
"""
A-10 推理公式法计算洪峰流量程序 —— 内核
=========================================
复刻《水利程序集》A-10 程序（作者：郝福良，水电部天津勘测设计院）。

功能：
  已知流域特征与暴雨特性（平均比降 J、最长距离 L、流域面积 F、
  平均入渗强度 U、汇流参数 M、暴雨递减指数 N、控制精度 E）及
  设计频率最大 24h 雨量 HP，用推理公式法计算最大洪峰流量 QM。

算法（推理公式法，水利院通用体系；原著公式为图像，数值结构经
     A-10.INT 双基准黑盒验证，2026-09-03 存档 a10_research*.py）：
  1. 雨力（设计频率最大时雨量）由 24h 雨量 HP 折算：
       Sp = HP · 24^(N−1)        （HP = Sp·24^(1−N)）
  2. 汇流历时 τ 与洪峰流量 QM 互为隐函数，用控制精度 E 迭代：
       τ = 0.278·L / (M·J^(1/3)·QM^(1/4))
       （L 单位 km，J 平均比降，M 汇流参数，QM 单位 m3/s）
  3. 产流（出流）历时：
       Tc = ((1−N)·Sp / U)^(1/N)
  4. 汇流状态判别与 QM 计算：
       Tc ≥ τ（全面汇流）：
         QM = 0.278·(Sp/τ^N − U)·F
       Tc < τ（部分汇流）：
         QM = 0.278·N·Sp·Tc^(1−N)·F / τ
     （系数 0.278 为 m3/s·km2·mm/h 单位换算；损失取常强度 U。）
  5. 迭代收敛：|QM_new − QM_old| ≤ E（控制精度，E 通常 1e-5，相对量级）。
     QM 初值 1.0 即可；收敛点与初值/判据形态无关（研究2 已证）。
  6. 打印时 QM 取整。

验证基准（A-10.INT 算例，原著 A-10.OUT）：
  J=0.0071, L=39.6 km, F=198 km², U=10 mm/h, M=1.2, N=0.7, E=1e-5；
  HP=482 → Sp=185.7735, τ=7.100, Tc=11.638（全面汇流）→ QM=2042.6→2043；
  HP=339 → Sp=130.6581, τ=8.311, Tc=7.039（部分汇流）→ QM=1087.7→1088。
  两值同时命中（±0）全 PASS。
"""
import os

from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-10"
TITLE = "推理公式法计算洪峰流量计算书"

# 输入参数顺序（原著 A-10.INT）：J, L, F, U, M, N, E
KEYS = ["J", "L", "F", "U", "M", "N", "E"]


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
    """解析 A-10.INT。data: INT 路径 | 文本 | dict。

    数据文件组织（原著）：
      首行（基础参数，逗号分隔）：
        J, L, F, U, M, N, E
      后续行（每个设计频率一行，可多个；缺省则必须由 dict 提供 hp）：
        HP1[, HP2, ...]
    兼容含 ^Z(0x1A) 的旧式文本；HP 也可经 dict 的 hp/hp_list/HP 提供。
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
        raise ValueError("A-10 无数据内容")

    head = [x for ln in lines[:1] for x in _tokenize(ln)]
    if len(head) < 7:
        raise ValueError("A-10 基础参数不足（需 J,L,F,U,M,N,E 7 个值）")
    base = {k: float(v) for k, v in zip(KEYS, head[:7])}
    if base["F"] <= 0 or base["L"] <= 0 or base["J"] <= 0:
        raise ValueError("A-10 流域参数非法（需 J>0、L>0、F>0）")
    if base["N"] <= 0 or base["N"] >= 1:
        raise ValueError(f"A-10 暴雨递减指数 N={base['N']} 非法（需 0<N<1）")
    if base["E"] <= 0:
        raise ValueError(f"A-10 控制精度 E={base['E']} 非法（需 >0）")

    hp = []
    for ln in lines[1:]:
        if not ln or ln.startswith("#") or ln.startswith("//"):
            continue
        hp.extend(float(x) for x in _tokenize(ln) if x > 0)
    return {**base, "HP": hp}


def _parse_dict(d):
    base = {k: float(d[k]) for k in KEYS}
    hp = d.get("HP") or d.get("hp") or d.get("hp_list") or []
    if not isinstance(hp, (list, tuple)):
        hp = [hp]
    if isinstance(hp, (list, tuple)):
        hp = [float(h) for h in hp]
    if not hp:
        raise ValueError("A-10 需要至少一个设计频率 HP")
    return {**base, "HP": hp}


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _solve_hp(J, L, F, U, M, N, E, HP):
    """对单一设计频率 HP 解推理公式，返回中间量与 QM。

    Sp=HP·24^(N−1)；τ=0.278L/(M·J^(1/3)·Q^(1/4))；
    Tc=((1−N)·Sp/U)^(1/N)；Tc≥τ 全面汇流 Q=0.278·(Sp/τ^N−U)·F；
    Tc<τ 部分汇流 Q=0.278·N·Sp·Tc^(1−N)·F/τ。收敛 |ΔQ|≤E。
    """
    Sp = HP * 24.0 ** (N - 1.0)
    Q = 1.0  # 洪峰流量初值（收敛点与初值无关）
    iters = 0
    tau = tc = 0.0
    for _ in range(500):
        tau = 0.278 * L / (M * J ** (1.0 / 3.0) * Q ** 0.25)
        tc = ((1.0 - N) * Sp / U) ** (1.0 / N)
        if tc >= tau:
            Qn = 0.278 * (Sp / tau ** N - U) * F      # 全面汇流
        else:
            Qn = 0.278 * N * Sp * tc ** (1.0 - N) * F / tau   # 部分汇流
        iters += 1
        if abs(Qn - Q) <= max(E, E * abs(Qn)):
            Q = Qn
            break
        Q = Qn
    if Q <= 0:
        raise ValueError(f"A-10 HP={HP} 计算 QM 非法（{Q}）")

    i_tau = Sp / tau ** N if tau > 0 else 0.0
    mode = "全面汇流(Tc≥τ)" if tc >= tau else "部分汇流(Tc<τ)"
    return {
        "设计频率24h雨量HP(mm)": HP,
        "雨力Sp(mm/h)": Sp,
        "汇流历时τ(h)": tau,
        "出流历时Tc(h)": tc,
        "汇流状态": mode,
        "τ历时最大雨强iτ(mm/h)": i_tau,
        "洪峰径流系数ψ": max(0.0, 1.0 - U / i_tau) if i_tau > 0 else 0.0,
        "洪峰流量QM(取整)": round(Q),
        "洪峰流量QM(全精度)": Q,
        "迭代次数": iters,
    }


def compute(params):
    """执行 A-10。params: parse 返回 dict（含 HP 列表）。"""
    J = params["J"]
    L = params["L"]
    F = params["F"]
    U = params["U"]
    M = params["M"]
    N = params["N"]
    E = params["E"]
    hp_list = params["HP"]
    if not hp_list:
        raise ValueError("A-10 无 HP（设计频率 24h 雨量）")

    cases = [_solve_hp(J, L, F, U, M, N, E, hp) for hp in hp_list]
    return {
        "程序": PROGRAM_ID,
        "输入": {
            "平均比降J": J,
            "最长距离L(km)": L,
            "流域面积F(km²)": F,
            "平均入渗强度U(mm/h)": U,
            "经验性汇流参数m": M,
            "暴雨递减指数N": N,
            "控制精度E": E,
        },
        "结果": {
            "计算频率数": len(cases),
            "各频率成果": cases,
        },
    }


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def render(params, result):
    """生成文本计算书（原著 .OUT 风格）。"""
    inp = result["输入"]
    cases = result["结果"]["各频率成果"]
    lines = []

    lines.append("输出基本数据:")
    lines.append(f"平均比降 J=  {inp['平均比降J']:.4f}")
    lines.append(f"最长距离 L=  {inp['最长距离L(km)']:.2f}公里")
    lines.append(f"流域面积 F=  {inp['流域面积F(km²)']:.2f}平方公里")
    lines.append(f"平均入渗强度 U=   {inp['平均入渗强度U(mm/h)']:.2f}毫米/小时")
    lines.append(f"经验性汇流参数 m=  {inp['经验性汇流参数m']:.2f}")
    lines.append(f"暴雨递减指数 N=  {inp['暴雨递减指数N']:.2f}")
    lines.append(f"控制精度循环值 E={inp['控制精度E']:.6f}")
    lines.append("")
    for case in cases:
        lines.append(f"设计频率为P的最大24小时雨量 HP=   {case['设计频率24h雨量HP(mm)']:.2f}")
        lines.append("计算结果:")
        lines.append(f" 最大洪峰流量 QM=  {case['洪峰流量QM(取整)']:8.1f}秒立米")
        lines.append("  汇流历时 τ= %.4f 小时  出流历时 Tc= %.4f 小时   %s"
                     % (case['汇流历时τ(h)'], case['出流历时Tc(h)'], case['汇流状态']))
        lines.append("  雨力 Sp= %.4f 毫米/小时   洪峰径流系数 ψ= %.4f"
                     % (case['雨力Sp(mm/h)'], case['洪峰径流系数ψ']))
        lines.append("*********************")
        lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text"):
    """统一入口。data: INT 路径 | dict（dict 须含 HP）"""
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
