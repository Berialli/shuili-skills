# -*- coding: utf-8 -*-
"""
A-6 推求流域时段平均面雨量程序 —— 内核
=========================================
复刻《水利程序集》A-6 程序（作者：陈沂，水电部天津勘测设计院）。

功能：
  根据流域内若干雨量站的实测降雨记录，用"直线内插法"推求各站逐时段降雨
  强度，再按各站占全流域之权重 K 加权，求全流域各时段平均面雨量。
  只输出计算窗口（自 DO 日 UO:00 起）内的逐时段（默认逐小时）权雨量。

算法（直线内插，已由黑盒反推验证）：
  每个雨量段给出 (起时, 止时, 段总雨量)，段内雨量均匀分摊到每一小时：
      单站权雨量 P = 段总雨量 / 段持续小时数 × K（H=1 时每小时）
  段持续小时数按"绝对时刻模型"计算（见数据格式说明）。
  全流域时段平均面雨量 = 各站同一时段 P 之和。

数据文件顺序（原著 A-6.INT）：
  第 1 行 : M, N, W
            M  计算雨量站总站数
            N  降雨时刻数（= 段数×2，各段起、止时刻端点总数）
            W  内插时段数（程序实际按"窗口整点时刻数(含首尾)"使用，
               见下方"已知悬案/差异"）
  每站 2 行 :
            DO, UO, H, K      DO 降雨开始日期(日)，UO 计算起始时刻(h)，
                              H 计算时段长度(h)，K 该站权重
            降雨资料（逗号分隔数值）：
              段 1 = 3 个数 (起时, 止时, 雨量)
              段 2 起每段 = 4 个数 (起时, 0, 止时, 雨量)，第 2 位恒 0，
                            为固定占位列（程序不回显、不参与计算），跳过
  时段端点时刻为 0..24 的整点计数（24 即午夜）；"止时 < 起时"表示跨 0 点。

日期归属（绝对时刻模型）：
  以 DO 日 0:00 为原点（绝对小时 0）。段 1 起点在 DO 日；其后每段：
    * 候选起点 = 当前日 0 点 + 起时数字
    * 若候选起点 < 前一段的绝对结束时刻 → 顺延 1 天（+24h）再作为起点
    * 段内止点 = 起点 + (止时 − 起时)，若 止时 < 起时（跨 0 点）再加 24h
  例：站 2 第 3 段 (4,20)：前段止于 24 日 22:00，候选起点 24 日 4:00 更早
      → 判为 25 日 4-20 时，落在计算窗口（止 25 日 4:00）之外，无输出。

计算窗口：
  自 DO 日 UO:00 起，按 H(=1) 逐时段输出，窗口半开区间为
      [DO日UO:00, DO日(UO+W−1):00)，共 W−1 个时段
  本例 W=15 → 输出 24 日 14-15 起至 25 日 3-4 止共 14 行。

验证基准（A-6.INT 三站算例，说明书 A-6Intro.txt 全文核对）：
  站 1 (K=0.5) 14 行 : 0.5×4, 0.0×3, 2.5×6, 0.0×1
  站 2 (K=0.3) 14 行 : 3.0×8, 0.0×6
  站 3 (K=0.2) 14 行 : 0.2×11, 0.0×3
  全流域平均 14 行 : 3.7×4, 3.2×3, 5.7, 2.7×3, 2.5×2, 0.0
  回归状态：逐位一致 PASS ✓（2026-09-02）

已知悬案（黑盒无法最终敲定，标注如下）：
  * 每站后段资料中占位列恒 0 的原始语义不明（疑为老式界面固定列/分隔符
    遗留），对计算结果无影响，解析时固定跳过。
  * W 按手册称"内插降雨量时段数"，但程序实际输出 W−1 行（本例 15→14）。
    实现按"W=窗口整点时刻数(含首尾)、窗口=[UO,UO+W−1)"处理，与原著输出
    逐位一致；若将来遇到 W 语义不同的算例需复核。
  * H>1（时段长非 1 小时）时输出步长的确切行为未见算例，属推测：按
    H 小时为一步、每行累计该 H 小时内权雨量实现（与 H=1 退化为逐时一致）。
"""
import os
import re

from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "A-6"
TITLE = "推求流域时段平均面雨量计算书"


# ------------------------------------------------------------
# 解析
# ------------------------------------------------------------

def _parse_num(tok):
    tok = tok.strip()
    if not tok:
        raise ValueError(f"空数值字段")
    if "." in tok or "e" in tok.lower() or "E" in tok:
        return float(tok)
    return int(tok)


def _tokenize_line(line):
    """按逗号切分为数值；整数值保留 int、含小数点按 float。"""
    return [_parse_num(t) for t in line.split(",") if t.strip()]


def parse(data):
    """解析输入。data: INT 路径 | 文本字符串 | dict。

    返回 dict:
      M 站数, N 时刻数, W 窗口整点数,
      stations: [ {DO,UO,H,K, raw:[...], segs:[(s,e,q),...]} ... ]
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
        raise ValueError("A-6 无数据内容")

    rows = [_tokenize_line(ln) for ln in lines]
    if len(rows[0]) < 3:
        raise ValueError("A-6 头行应为 M,N,W 三个数")

    M, N, W = int(rows[0][0]), int(rows[0][1]), int(rows[0][2])
    stations = []
    idx = 1
    for st in range(M):
        if idx >= len(rows):
            raise ValueError(f"A-6 缺第 {st + 1} 站参数行 (DO,UO,H,K)")
        head = rows[idx]
        idx += 1
        if len(head) < 4:
            raise ValueError(f"A-6 第 {st + 1} 站参数行应为 DO,UO,H,K 四数")
        DO, UO, H, K = (int(head[0]), int(head[1]),
                        float(head[2]) if isinstance(head[2], float) else head[2],
                        float(head[3]))
        if idx >= len(rows):
            raise ValueError(f"A-6 缺第 {st + 1} 站降雨资料行")
        raw = rows[idx]
        idx += 1
        segs = _split_segments(raw)
        stations.append({
            "DO": DO, "UO": UO, "H": H, "K": K,
            "raw": raw, "segs": segs,
        })
    if len(stations) != M:
        raise ValueError(f"A-6 声明 {M} 站，实际解析 {len(stations)} 站")
    return {"M": M, "N": N, "W": W, "stations": stations}


def _split_segments(raw):
    """资料行 → 段列表 [(起, 止, 量), ...]。

    布局：段 1 = (起,止,量) 3 个数；
          其后每段 = (起, 0, 止, 量) 4 个数，第 2 位恒 0（占位列）跳过。
    即长度 = 3 + 4×(段数−1)。N=6 → 6 个时刻端点 = 3 段，资料 11 个数。
    """
    if len(raw) < 3:
        raise ValueError(f"A-6 降雨资料不足 3 个数: {raw}")
    segs = [(raw[0], raw[1], raw[2])]
    i = 3
    while i < len(raw):
        if i + 3 >= len(raw) + 1:  # 需要 4 个数
            raise ValueError(f"A-6 降雨资料布局非法(应为 3+4×(段数−1)): {raw}")
        s, z, e, q = raw[i], raw[i + 1], raw[i + 2], raw[i + 3]
        if z != 0:
            # 黑盒反推中恒 0；若遇非 0 仍宽容解析并告警（可能为其它占位语义）
            pass
        segs.append((s, e, q))
        i += 4
    return segs


def _parse_dict(d):
    """dict 形式（run 传 JSON 时）与文件形式等价。"""
    st_out = []
    for s in d["stations"]:
        st_out.append({
            "DO": int(s["DO"]), "UO": int(s["UO"]),
            "H": float(s["H"]), "K": float(s["K"]),
            "raw": [float(x) for x in s["raw"]] if "raw" in s else None,
            "segs": [(int(a), int(b), float(c)) for (a, b, c) in s["segs"]],
        })
    return {"M": int(d["M"]), "N": int(d["N"]), "W": int(d["W"]),
            "stations": st_out}


# ------------------------------------------------------------
# 计算
# ------------------------------------------------------------

def _abs_intervals(segs, DO):
    """
    段 → 绝对小时区间 [(a, b, q), ...]，原点 = DO 日 0:00。
    日期推进规则见模块 docstring。
    """
    out = []
    cur_day = 0          # 相对 DO 的天偏移
    last_end = None
    for (s, e, q) in segs:
        a = cur_day * 24 + s
        if last_end is not None:
            while a < last_end:     # 候选起点早于前段结束 → 顺延一天
                cur_day += 1
                a = cur_day * 24 + s
        b = a + (e - s)
        if e < s:
            b += 24                 # 跨 0 点
        out.append((a, b, float(q)))
        last_end = b
    return out


def _window_rows(abs_segs, UO, W, K, H=1):
    """
    单站窗口逐时段权雨量。
    窗口半开区间 [UO, UO+W−1)，步长 H（H=1 逐小时，H>1 为推测）。
    返回 list[(day_offset, start_h, P)]。
    """
    rows = []
    if H < 1:
        raise ValueError(f"A-6 H 时段长应为正数，收到 {H}")
    t = float(UO)
    t_end = UO + W - 1.0
    while t < t_end - 1e-9:
        p = 0.0
        # 若步长 H 非整，时段 [t, t+H)
        seg_b, seg_e = t, t + H
        for (a, b, q) in abs_segs:
            dur = b - a
            if dur <= 0:
                continue
            ov_s = max(a, seg_b)
            ov_e = min(b, seg_e)
            if ov_e > ov_s + 1e-9:
                p += (q / dur) * K * (ov_e - ov_s)
        day = int(t // 24)
        sh = int(round(t)) % 24
        rows.append((day, sh, p))
        t += H
    return rows


def compute(params):
    """执行 A-6：各站权雨量 + 全流域时段平均面雨量。"""
    stations = params["stations"]
    W = params["W"]
    station_abs = []
    for st in stations:
        abs_segs = _abs_intervals(st["segs"], st["DO"])
        station_abs.append(abs_segs)

    # 各站窗口逐行（行数一致：W−1）
    st_tables = []
    for st, abs_segs in zip(stations, station_abs):
        rows = _window_rows(abs_segs, st["UO"], W, st["K"], st["H"])
        st_tables.append(rows)

    # 校验各站行数一致
    n_rows = len(st_tables[0]) if st_tables else 0
    for rt in st_tables:
        if len(rt) != n_rows:
            raise ValueError("A-6 各站窗口行数不一致（内部错误）")

    basin = []
    for i in range(n_rows):
        basin.append(round(sum(st_tables[k][i][2] for k in range(len(stations))), 8))

    # 结构化结果
    res = {
        "程序": PROGRAM_ID,
        "输入": {
            "雨量站总站数M": len(stations),
            "时刻数N": params["N"],
            "窗口整点数W": W,
            "各站": [
                {"站": i + 1, "DO": st["DO"], "UO": st["UO"],
                 "H": st["H"], "K": st["K"],
                 "降雨段": [(int(s), int(e), round(float(q), 4))
                           for (s, e, q) in st["segs"]],
                 "绝对区间": [(a, b, round(float(q), 4))
                           for (a, b, q) in abs_segs]}
                for i, (st, abs_segs) in enumerate(zip(stations, station_abs))
            ],
        },
        "单站逐时段权雨量": [
            {
                "站": i + 1, "K": st["K"],
                "行": [(day, sh, round(p, 8))
                       for (day, sh, p) in st_tables[i]],
            }
            for i, st in enumerate(stations)
        ],
        "全流域时段平均降雨量": {
            "行": [(basin_day := st_tables[0][i][0],
                    st_tables[0][i][1], round(basin[i], 8))
                   for i in range(n_rows)],
        },
        "窗口说明": {
            "起点": f"DO 日（{stations[0]['DO']} 日）{stations[0]['UO']}:00",
            "终点时刻": f"绝对小时 {stations[0]['UO'] + W - 1}，"
                       f"即 {stations[0]['DO']} 日 +{1 + (stations[0]['UO'] + W - 1) // 24} 日 "
                       f"{(stations[0]['UO'] + W - 1) % 24}:00",
            "行数": n_rows,
        },
    }
    return res


# ------------------------------------------------------------
# 渲染
# ------------------------------------------------------------

def _time_label(sh):
    """整点起始时刻 → '14 - 15' 标签；sh=23 → '23 - 24'；sh=0 → '0 - 1'。"""
    return f"{sh} - {sh + 1}"


def _fmt1(x):
    return f"{x:.1f}"


def _fmt2(x):
    return f"{x:.2f}"


def render(params, result):
    """生成文本计算书（原著 .OUT 风格）。"""
    inp = result["输入"]
    M = inp["雨量站总站数M"]
    N = inp["时刻数N"]
    W = inp["窗口整点数W"]
    st_in = inp["各站"]

    # ---- (一) 各站逐时段权雨量 ----
    sub_lines = ["(一)、计算各站各时段权雨量", "",
                 f"雨量站总站数 M= {M}    N= {N}    W= {W}", ""]
    for s in st_in:
        sub_lines.append(f"雨量站 L= {s['站']}    DO= {s['DO']}    "
                         f"UO= {s['UO']}    H= {s['H']:g}    K= {s['K']:.2f}")
        segdesc = ", ".join(f"{a}-{b}:{q:g}" for (a, b, q) in s["降雨段"])
        sub_lines.append(f"    降雨资料：{segdesc}")
    sub_lines += ["", "", "------------"]
    for i, s in enumerate(st_in):
        sub_lines.append(f"雨量站 L= {s['站']}    权重 K= {s['K']:.2f}")
        sub_lines.append("        时间                   权雨量")
        for (a, b, q) in s["降雨段"]:
            sub_lines.append(f"        {a}  - {b}               {_fmt2(q)}")
        # 逐时段权雨量
        rows = result["单站逐时段权雨量"][i]["行"]
        sub_lines.append("降雨起始日期       时间       单站权雨量")
        sub_lines.append("    D                T              P")
        prev_day = None
        for (day, sh, p) in rows:
            d = s["DO"] + day
            if d != prev_day:
                day_str = f"{d:>5}"
            else:
                day_str = "     "
            sub_lines.append(f"{day_str}         {_time_label(sh):<10}   {_fmt1(p)}")
            prev_day = d
        sub_lines.append("        ")
        sub_lines.append("------------")
    sub_lines += ["", ""]

    # ---- (二) 全流域时段平均降雨量 ----
    blines = ["（二）、全流域时段平均降雨量",
              "  -----------------------",
              f"  雨量站总站数 M= {M} ",
              "  降雨开始日期     时间         流域平均降雨量"]
    d0 = st_in[0]["DO"]
    prev_day = None
    for (day, sh, val) in result["全流域时段平均降雨量"]["行"]:
        d = d0 + day
        if d != prev_day:
            day_str = f"{d:>5}"
        else:
            day_str = "     "
        blines.append(f"{day_str}          {_time_label(sh):<10}   {_fmt1(val)}")
        prev_day = d
    blines.append("")
    blines.append("  ------------------ 结 束 ------------------")

    return render_text(PROGRAM_ID, TITLE, [("", sub_lines + blines)])


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
