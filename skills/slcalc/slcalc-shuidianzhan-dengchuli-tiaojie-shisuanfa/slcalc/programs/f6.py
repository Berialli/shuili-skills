# -*- coding: utf-8 -*-
"""
F-6 任意高程网平差计算程序 —— 内核
==================================
复刻《水利水电工程设计计算程序集》F-6（作者：谢希哲，新疆兵团勘测设计院）。

功能
----
对水准网（L）或三角高程网（T）作间接观测平差：解算各结点高程、结点高程中误差、
各边高差改正数与平差后高差，并给出平差后单位权中误差（每公里中误差）。

公式（F-6Intro「一. 简介」第 4~7 条 + 权威 F-6-1.OUT 逐位反演）
--------------------------------------------------------------
  · 未知数        ：各结点（未知点）高程 H_k
  · 观测方程      ：H_j − H_i = Δh_ij（观测）＋ v
  · 定权          ：水准网   P = 1/S      （S = 边长，km）
                    三角高程  P = 1/S²    （三角高程网多条边串联时 S = √ΣS_i²）
                    依据：本程序以 1km 边长的权为单位权（P=1），
                    故表下方 m₀ 即「每公里中误差」。
  · 单位权中误差  ：m₀ = ±√( [p·v·v] / (n − u) )      n=边数，u=结点数
  · 结点高程中误差：Mh_k = m₀·√(Q_kk)，  Q = N⁻¹，N = Aᵗ·P·A
  · 平差后高差    ：Δĥ = Δh + v；高程逐边推算，顺序由边号保证。
  逐位反演证据（F-6-1，L 网，n=7、u=6、r=1）：
    m₀ = √([vv/S]/1) = 0.025862 m/√km → 打印 25.9 mm ✓
    Q_kk = S_前·S_后/S_总（单链闭式），Mh = m₀√Q_kk = 9.03/12.51/13.50/12.61/9.30/5.65 mm
    → 打印 9.0 / 12.5 / 13.5 / 12.6 / 9.3 / 5.6 ✓（逐位命中）
    V = −f_h·S_i/ΣS（f_h = ΣΔh − ΔH_已知）→ 3.467/4.954/4.954/4.954/4.954/2.477/1.239 mm
    → 打印 3.5 / 5.0 / 5.0 / 5.0 / 5.0 / 2.5 / 1.2 ✓

输入数据（.INT）
----------------
  第 1 行： "网名","计算者","日期",网型 或等级 等
      "网名","计算者","日期","L|T",等级,结点数,边数
  随后为点名表（共 结点数 + 已知点数 个）：
      已知点： 点名行（引号）＋ 高程行（裸数）
      未知点： 点名行（引号）
  再随 边数 行高差观测： 起点号,终点号,高差Δh(m),边长S(km)
  边号顺序须能逐边推算高程（与手算一致）。

输出（逐字复刻原著 .OUT 版式，GBK；左右双表）
--------------------------------------------
  首行「文件：…」为运行期路径回显，本内核不生成。
  左表：No./起点/终点/Δh(观)/S(km)/V(mm)/Δh(平)；
  右表：No./点名/高程(H)/Mh(mm)。左右两栏按显示宽度 62 列拼接
  （全角制表符按 2 列计）。

基准与闭合状态
--------------
  · 权威 OUT：F-6-1.OUT（2331 B，27 行，L 网 支线算例）—— 逐行对拍。
  · F-6-2（L 网）、F-6-3（T 网）G 盘无 OUT；T 网取位口径按说明书「三角高程网
    高程给到 0.01m」实现，计 DECL（已量化）。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
----------------------------------------------
  母本：`02_水利教材精读\\水利工程施工与管理\\水利工程测量_第5版_精读笔记.md`
  逐条对照：
  | 程序公式 | 库中出处 | 是否一致 | 处置 |
  |---|---|---|---|
  | 最小二乘平差 [pvv]=最小 | §2.5「平差：最小二乘原理（[vv]=最小）」；§4.3；§5.1 关联表 | **一致** | 采用 |
  | 定权 P=μ²/m_i²（本程序 1km 单位权） | §2.5 + 公式 5-21「P_i = μ²/m_i²」；§4.1 不等精度平差 | **一致** | 采用（水准 P=1/S、三角高程 P=1/S²） |
  | 闭合差按路线长成比例分配 Δh_i=−f_h·L_i/ΣL | 公式 2-13「Δh_i = −f_h·(L_i/ΣL) 或 −f_h·(n_i/Σn)」；§4.1 第 5 步 | **一致** | 采用（V=−f_h·S_i/ΣS，由 P=1/S 导出，与最小二乘严格等价） |
  | 精度评定 Mh=m₀·√Q_kk | 公式 5-12（线性函数误差传播）、5-14（一般函数误差传播）；§2.5 | 原则**一致**（库给出传播通式） | 采用 Q=N⁻¹ 协因数阵 |
  | 单位权中误差 m₀=±√([pvv]/(n−u)) | §2.5「单位权中误差」；§4.1 | **一致** | 采用 |
  | 高程网分级 / 水准点 BM 命名 | §2.6「国家一~四等水准网（水准点 BM）；水利水电分三级」 | **一致** | 采用（高差表以 BM/P + 序号显示） |
  | 三角高程网串联边长 S=√ΣS_i² | F-6Intro 第 5 条；库中公式 7-1、7-4 未给出串联定权式 | 库中无对应条款（说明书给出） | 按说明书实现，计 DECL（无基准） |

  分歧说明：库中《水利工程测量》第 5 版以「闭合/附合水准路线」手算分配为主，
  未展开任意高程网的间接观测平差；两者在「[pvv]=最小 + P=μ²/m²」层面同源，
  差异仅在线路 vs 网形的解算规模 —— 非分歧，属教材深度差异。
"""
import math
import os

from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "F-6"
TITLE = "任意高程网平差计算程序"
AUTHOR = "谢希哲(新疆兵团勘测设计院)"
HEAD_NAME = "F-6"

# ---------------------------------------------------------------- 版式常量
STAR = "*" * 71
TITLE_LINE = " ****          任 意 高 程 网 平 差 计 算   F-6  (99.6版)           ****"

L_TOP = "┌──┬───┬───┬────┬────┬───┬────┐"
L_HDR = "│ No.│ 起点 │ 终点 │ Δh(观)│  S(km) │ V(mm)│ Δh(平)│"
L_SEP = "├──┼───┼───┼────┼────┼───┼────┤"
L_BOT = "└──┴───┴───┴────┴────┴───┴────┘"
R_TOP = "┌──┬──────┬────┬───┐"
R_HDR = "│ No.│  点    名  │ 高程(H)│Mh(mm)│"
R_SEP = "├──┼──────┼────┼───┤"
R_BOT = "└──┴──────┴────┴───┘"

COL_W = 62          # 左栏显示宽度（全角按 2 列）


def disp_width(s):
    return sum(2 if ord(c) > 127 else 1 for c in s)


def pad_disp(s, w=COL_W):
    n = w - disp_width(s)
    return s + (" " * n if n > 0 else "")


def split_csv(line):
    """按逗号切分一行，保留双引号内的逗号并去引号。"""
    out, cur, inq = [], [], False
    for c in line:
        if c == '"':
            inq = not inq
            continue
        if c == "," and not inq:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(c)
    out.append("".join(cur).strip())
    return out


# ============================================================
# 解析
# ============================================================

def parse(data):
    if isinstance(data, dict):
        p = dict(data)
        p.setdefault("点", [])
        p.setdefault("边", [])
        p.setdefault("网型", "L")
        p.setdefault("网名", "")
        p.setdefault("计算者", "")
        p.setdefault("日期", "")
        p.setdefault("等级", None)
        return p
    lines = read_lines(data)
    head = split_csv(lines[0])
    if len(head) < 7:
        raise ValueError("F-6 首行字段不足：%s" % lines[0])
    netname, author, date = head[0], head[1], head[2]
    wtype = head[3].strip().upper() or "L"
    grade = float(head[4])
    n_node = int(float(head[5]))
    n_edge = int(float(head[6]))

    pts, edges = [], []
    i = 1
    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            i += 1
            continue
        if ln.startswith('"'):
            name = ln.strip().strip('"')
            h = None
            if i + 1 < len(lines):
                nxt = lines[i + 1].strip()
                if nxt and ("," not in nxt):
                    try:
                        h = float(nxt)
                        i += 1
                    except ValueError:
                        h = None
            pts.append({"名": name, "H": h, "已知": h is not None})
            i += 1
            continue
        break
    while i < len(lines):
        ln = lines[i].strip()
        i += 1
        if not ln or "," not in ln:
            continue
        f = [x for x in ln.split(",") if x.strip() != ""]
        if len(f) < 4:
            continue
        edges.append((int(float(f[0])), int(float(f[1])),
                      float(f[2]), float(f[3])))
    return {"网名": netname, "计算者": author, "日期": date, "网型": wtype,
            "等级": grade, "结点数": n_node, "边数": n_edge,
            "点": pts, "边": edges, "源": str(data)}


# ============================================================
# 计算
# ============================================================

def _solve(A, P, l):
    """最小二乘 A x = l（权 P）→ x，N⁻¹，v。纯 Python（无 numpy 依赖）。"""
    nrow = len(A)
    ncol = len(A[0]) if nrow else 0
    N = [[0.0] * ncol for _ in range(ncol)]
    w = [0.0] * ncol
    for r in range(nrow):
        pr = P[r]
        for a in range(ncol):
            if A[r][a] == 0.0:
                continue
            w[a] += pr * A[r][a] * l[r]
            for b in range(a, ncol):
                if A[r][b] == 0.0:
                    continue
                N[a][b] += pr * A[r][a] * A[r][b]
    for a in range(ncol):
        for b in range(a):
            N[a][b] = N[b][a]
    Q, x = _inv_solve(N, w)
    v = [0.0] * nrow
    for r in range(nrow):
        s = 0.0
        for a in range(ncol):
            if A[r][a] != 0.0:
                s += A[r][a] * x[a]
        v[r] = s - l[r]
    return x, Q, v


def _inv_solve(N, w):
    """高斯—约当法解 N x = w，并同时求 N⁻¹。"""
    n = len(N)
    M = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] + [w[i]]
         for i, row in enumerate(N)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-300:
            raise ValueError("F-6 法方程奇异：网形不足以解算（第 %d 列）" % (c + 1))
        M[c], M[piv] = M[piv], M[c]
        d = M[c][c]
        M[c] = [v / d for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0.0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    # [N | I | w] → [I | N⁻¹ | x]
    Q = [[M[i][n + j] for j in range(n)] for i in range(n)]
    x = [M[i][2 * n] for i in range(n)]
    return Q, x


def compute(p):
    pts = p["点"]
    edges = p["边"]
    wtype = p["网型"].upper()
    n_pt = len(pts)
    idx_unk = [k for k, q in enumerate(pts) if not q["已知"]]
    u = len(idx_unk)
    n = len(edges)
    pos = {k: j for j, k in enumerate(idx_unk)}

    A, P, l = [], [], []
    for (i, j, dh, S) in edges:
        if not (1 <= i <= n_pt and 1 <= j <= n_pt):
            raise ValueError("F-6 边端点越界：%d→%d（共 %d 点）" % (i, j, n_pt))
        if S <= 0:
            raise ValueError("F-6 边长须为正：边 %d→%d S=%g" % (i, j, S))
        row = [0.0] * u
        if (j - 1) in pos:
            row[pos[j - 1]] += 1.0
        if (i - 1) in pos:
            row[pos[i - 1]] -= 1.0
        h0i = pts[i - 1]["H"] or 0.0
        h0j = pts[j - 1]["H"] or 0.0
        A.append(row)
        P.append(1.0 / (S * S) if wtype == "T" else 1.0 / S)
        l.append(dh - (h0j - h0i))

    if u == 0:
        raise ValueError("F-6 无未知结点（全部为已知点）")
    x, Q, v = _solve(A, P, l)
    r = n - u
    pvv = sum(P[k] * v[k] * v[k] for k in range(n))
    if r > 0:
        m0 = math.sqrt(pvv / r)
    else:
        m0 = 0.0          # 无多余观测：无精度评定

    H = [q["H"] for q in pts]
    for j, k in enumerate(idx_unk):
        H[k] = x[j]
    # 逐边推算（边号顺序保证起点高程已知/已推出），作为独立校验
    Hseq = [q["H"] for q in pts]
    ok_seq = True
    for k, (i, j, dh, S) in enumerate(edges):
        if Hseq[i - 1] is None:
            ok_seq = False
            break
        Hseq[j - 1] = Hseq[i - 1] + dh + v[k]
    seq_res = None
    if ok_seq:
        d = [abs(Hseq[j] - H[j]) for j in range(n_pt) if H[j] is not None]
        seq_res = max(d) if d else None

    Mh = [None] * n_pt
    for j, k in enumerate(idx_unk):
        Mh[k] = m0 * math.sqrt(max(Q[j][j], 0.0))

    rows = []
    for k, (i, j, dh, S) in enumerate(edges):
        rows.append({"序号": k + 1, "起": i, "终": j, "dh": dh, "S": S,
                     "V": v[k] * 1000.0, "dh平": dh + v[k]})
    prows = []
    for k in range(n_pt):
        prows.append({"序号": k + 1, "名": pts[k]["名"], "H": H[k],
                      "Mh": None if Mh[k] is None else Mh[k] * 1000.0,
                      "已知": pts[k]["已知"]})
    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR,
            "输入": {"网名": p["网名"], "计算者": p["计算者"], "日期": p["日期"],
                     "网型": wtype, "等级": p.get("等级"), "结点数": u,
                     "边数": n, "点数": n_pt},
            "高差表": rows, "高程表": prows, "m0": m0,
            "Mo": m0 * 1000.0, "pvv": pvv, "r": r, "逐边闭合残差": seq_res,
            "H": H}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

def render(p, res):
    wtype = res["输入"]["网型"]
    pre = "BM" if wtype == "L" else "P"
    hdec = 3 if wtype == "L" else 2
    B = []
    B.append("")
    B.append(" " + STAR)
    B.append(TITLE_LINE)
    B.append(" " + STAR)
    B.append("")
    B.append(" " * 35 + "网名:" + res["输入"]["网名"])

    L = [L_TOP, L_HDR, L_SEP]
    for row in res["高差表"]:
        L.append("│%3d │%-6s│%-6s│%8.3f│%8.3f│%6.1f│%8.3f│" % (
            row["序号"], pre + "%3d" % row["起"], pre + "%3d" % row["终"],
            row["dh"], row["S"], row["V"], row["dh平"]))
        L.append(L_SEP)
    L[-1] = L_BOT
    L.append("        计算者:%s      日期:%s" % (res["输入"]["计算者"],
                                                res["输入"]["日期"]))
    L.append("        Mo=± %4.1fmm" % res["Mo"])

    R = [R_TOP, R_HDR, R_SEP]
    for pr in res["高程表"]:
        mhs = "      " if pr["Mh"] is None else "%6.1f" % pr["Mh"]
        R.append("│%3d │%-12s│%8.*f│%s│" % (pr["序号"], pr["名"], hdec,
                                            pr["H"], mhs))
        R.append(R_SEP)
    R[-1] = R_BOT

    for k in range(max(len(L), len(R))):
        lft = L[k] if k < len(L) else ""
        rgt = R[k] if k < len(R) else ""
        B.append(pad_disp(lft) + ((" " + rgt) if rgt else ""))
    return "\n".join(B) + "\n"


def run(data, out_txt=None, out_json=None, fmt="text", **kw):
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
    _r, _t = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_t)
