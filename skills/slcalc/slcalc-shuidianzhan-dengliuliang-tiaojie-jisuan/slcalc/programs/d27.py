# -*- coding: utf-8 -*-
"""
D-27 压力供水管道管径选择程序 —— 内核
======================================
复刻《水利水电工程设计计算程序集》公之于众版 D-27（作者：张校正，新疆水利厅）。

功能
----
树枝状管网（重力流压力管道）的**管径选择**：已知各节点名称、桩号、设计流量、
管壁类型、局部水头损失百分比、地面高程与要求的水头压力，逐段配置管径并给出
每段两种相邻标准管径的长度（变径点），免去试算。

原著原理（D-27Intro「二、选择管径的计算原理」）
------------------------------------------------
已知设计流量 Q、管道总长 L、管道两端水位差 Hw；选定 D1、L1、D2、L2：
全线用较小的 D1 时水头损失偏大（不够），用较大的 D2 时偏小（有富裕）；
两种**相邻标准管径**配合可使 Σi·L 恰等于 Hw —— 即 2×2 线性方程组
      i1·L1 + i2·L2 = Hw ，  L1 + L2 = L
（i 为每米水头损失）。解出：
      L2 = (Hw − L·i1)/(i2 − i1)  （本内核按「大管径在前」排布，与权威 OUT 一致）
程序自起始点逐段配置，算出变径点桩号；变径点名称 = 其**上游**整桩号节点名前加 'J'；
最后列出各节点的水位 Z、地面高程 GD 与水头压力 GY = Z − GD。

9 种管壁类型（D-27Intro「二」）与其单米水头损失 i = C·Q^a·d^(−b)
--------------------------------------------------------------------
  类型           名称                                     a        b        C            来源
  1  混凝土及钢筋混凝土管 n=0.015  曼宁     2        16/3     K_M·n²       卷8 第4章「算法溯源」
  2  混凝土及钢筋混凝土管 n=0.014  曼宁     2        16/3     K_M·n²       （λ=8g·n²/R^(1/3)，R=D/4）
  3  混凝土及钢筋混凝土管 n=0.013  曼宁     2        16/3     K_M·n²
  4  旧钢管及旧铸铁管                        1.90     5.10     6.25e5※      卷9 §5.4 表5.4-9
  5  石棉水泥管                              —— 未标定（知识库未收录、无权威算例）
  6  硬聚氯乙烯 UPVC (CECS 82:96)            1.761    4.761    0.000918437  ★权威 D-27-1.OUT 反演
  7  聚乙烯 PE / 聚丙烯 PP (CECS 82:96)      1.774    4.774    0.000912      D-27vb.EXE 内嵌常量（未标定）
  8  塑料硬管 (《喷灌系统技术规程》)          1.77     4.77     0.948e5※     卷9 §5.4 表5.4-9
  9  铝管及铝合金管                          —— 未标定（知识库未收录、无权威算例）
  ※ 类型 4/8 的原表 f 以 Q(m³/h)、d(mm) 给出，本内核按 C_SI = f·3600^a/1000^b 折算为 SI。

9 种管壁类型（D-27Intro「二」；f/m/b 取自 D-28Intro「二、计算原理」原表）
--------------------------------------------------------------------------------
  类型           名称                                    a        b        C            来源
  1  混凝土及钢筋混凝土管 n=0.013（谢才公式）   2        5.33     0.0131        D-28Intro 原表
  2  混凝土及钢筋混凝土管 n=0.014（谢才公式）   2        5.33     0.01553       D-28Intro 原表
  3  混凝土及钢筋混凝土管 n=0.014（谢才公式）   2        5.33     0.01782       D-28Intro 原表
  4  旧钢管、旧铸铁管（安德列维舍夫公式）       1.9      5.1      0.00625       D-28Intro 原表
  5  石棉水泥管（阿勃拉莫夫公式）               1.85     4.89     0.001455      D-28Intro 原表
  6  硬聚氯乙烯 UPVC（《农村给水设计规范》公式）1.761    4.761    0.000918437   ★权威 OUT 反演（表值 0.000875）
  7  聚乙烯 PE、聚丙烯 PP（《农村给水设计规范》）1.774    4.774    0.000951      D-28Intro 原表
  8  塑料硬管（谢维列夫公式）                   1.77     4.77     0.000946      D-28Intro 原表
  9  铝管及铝合金管（武汉水利电力学院公式）     1.74     4.74     0.000861      D-28Intro 原表
 10  旧钢管 / 铸铁管按流速分段（《农村给水设计规范》）—— Intro 未给出公式 → 未实现

※ D-28Intro 原表把单位写作「Q（L/s）、d（m）」，但该口径**与 Intro 自身算例及本程序行为
   均不符**：按 L/s 计，材料 6 的 i 将达 1895 m/m（算例同点仅 0.0104 m/m），偏 5 个量级。
   由 Intro 例题（站 0/1200/12000，Q=0.043238 m³/s，D=0.188/0.150，损失 12.47/328.75 m）
   反解，f 必须配以 **Q（m³/s）、d（m）**，即原表单位注记之「L/s」系文档笔误
   （同一原表内各行 f 亦无法用任何单一单位折算互相自洽）。本内核统一按 Q m³/s、d m 实现。

常量反演与裁决（详见模块尾「反演与未闭合」）
-------------------------------------------
  π = 3.14159（仅用于 V = Q/A 显示；与 D-17 同口径）
  g = 9.81    （仅用于 V²/2g；本程序 D-27 不打印 V²/2g，D-28 打印且 2 位小数不敏感）
  类型 6：由权威 D-27-1.OUT 的 8 个跨度（每跨度两段、水头恰好用尽）与 D-28.OUT 的
  总水头约束联合反演得 a = 1.761、b = 4.761（与 D-28Intro 原表 m/b **逐位相同**，
  且 b − a ≡ 3，即 λ ∝ Re^−0.239 的达西族结构）；C 由「每跨度恰好用尽可用水头」唯一确定。
  D-27vb.EXE 内嵌 double 常量恰为 1.761000633 / 4.761001587（= float32(1.761)/float32(4.761)）
  与 0.0009120004252，指数与原表、与 OUT 反演三方吻合；系数见表值 0.000875（详见「未闭合 ②」）。

反演与未闭合（量化 + 反证 + 出处）
----------------------------------
① 指数 (a, b) = (1.761, 4.761)：**三路互证**。
   出处 A：D-28Intro「二、计算原理」原表第 6 行（硬聚氯乙烯 UPVC 管）明列 m=1.761、b=4.761；
   出处 B：D-27vb.EXE / D-28vb.exe / D-28vb0.EXE 三个 EXE 的常量池内均含相邻存放的
           double 1.761000633 与 4.761001587（= float32(1.761)/float32(4.761)）；
   出处 C：以权威 D-27-1.OUT 8 个跨度的「显示 L」为约束做 (a,b) 二维扫描
           （a∈[1.70,1.82]、b∈[4.70,4.85]，步长 1e-3，另 1e-4 细扫），
           使 8 跨度 C 相对离散最小者恒为 (1.761, 4.761)（离散 2.3e-5）；
           相邻候选 (1.760,4.760) 离散 3.5e-3、(1.774,4.774) 离散 4.5e-2，均劣 2~3 个量级。
   结构旁证：b − a = 3.000 精确成立，与 i = λ·v²/(2gd)、λ ∝ Re^−p（a = 2−p, b = 5−p）
           的达西形式自洽（p = 0.239）。
② 系数 C6 = 0.000918437（**权威 OUT 反演值**）vs 0.000875（**D-28Intro 原表值**），差 +4.96%。
   反证 1：以 C = 0.000875 计，D-28 算例首段（L=200、d=0.188）损失为 1.98 m，权威 2.08 m，
           相对差 −4.96%——8 个跨度显示 L 将全部偏移，**0/8 命中**；
           以 C = 0.000918437 计，8 跨度显示 L 有 6/8 逐位命中（余 2 项差 0.01 m），
           且 D-28.OUT 全 41 行逐字命中。
   反证 2：C 与指数无关地由「每跨度恰好用尽 H」唯一确定（8 个跨度的相容区间宽仅 2.3e-5
           相对），故该值不是拟合自由度，而是被数据钉死的量。
   反证 3：三值递进 0.000875（原表）< 0.0009120004252（EXE 内嵌常量，+4.23%）
           < 0.000918437（OUT 反演，+4.96%），且 D-27-1.OUT 存档于 2010-09-25，
           早于三个 EXE（2014-10-08 / 2018-02-13）——即原表、EXE、OUT 分属三次不同构建。
   处置：按任务「不确定时以 OUT 为准」，内核取 C6 = 0.000918437；0.000875 记原表值，
        0.000912 记 EXE 内嵌值，三值并列备查。
③ 局部水头损失百分比**未参与计算**（D-28Intro「2、局部水头损失可以按沿程水头损失的比例
   （5%～10%）考虑」有述，但两程序均未消费；由 OUT 反演裁决）。
   反证：8 个跨度的「大径段 + 小径段」打印损失之和恰等于该跨度可用水头 H
        （C→10：0.12+0.38=0.50=H，逐跨度成立），即 Σi·L ≡ H 已由纯沿程损失用尽；
        若再叠加 loc=0.1 的局部损失（1.1 倍），则 H 需放大 10%，C 将为 0.0008349，
        8 跨度显示 L 全部偏差 3~25 m。故原著该输入项在本程序内**未被消费**。
        （D-28 同：24 段损失之和 = 341.21 = Zo，未叠加 10%；若叠加则解出 Q=0.040940，
          与权威 0.043238 不符。）
④ 未闭合：跨度 1、2 的变径点桩号/距离 4 项差 0.01 m（= 打印末位 1 个单位）。
   量化：跨度1 模型 L大=144.5166（权威 144.51，Δ=+0.0066 m）；跨度2 模型 350.7105
        （权威 350.70，Δ=+0.0105 m）；其余 6 跨度 |Δ| ≤ 0.0046 m。
   反证：① 8 跨度联立的 C 可行区间交集为空——跨度 1,2 要求 C ≤ 0.9184299e-3，
        跨度 3,6 要求 C ≥ 0.9184356e-3，缺口 6e-9（相对 6.5e-6）；
        ② 在 (a,b) ∈ [1.755,1.770]×[4.740,4.770]（步长 1e-4）与
          i=C·Q^a·(d+δ)^−b、i=C·Q^a·d^−b·(1+k/d)（δ,k∈[−5,5]e-3，步长 1e-5）
          两类修正模型下，max|ΔL| 的下确界恒为 0.0105 m（最优 δ=k=0）；
        ③ float32 单精度模拟（4 种 f32 组合）与 double 结果逐位相同，排除精度解释。
        结论：纯幂律模型已达该数据集的**精度下界**（打印分辨率 0.01 m），残差属原著
        舍入/刷值边界，非公式差异。verify 记 DECL（容差 0.011 m），非 FAIL。
⑤ 版式：权威 D-27-1.OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-27-1.out」为
   原著运行期路径回显（机器相关），本内核不生成；比对自第 2 行起算（与 D-1/D-3/D-4/
   D-7/D-8/D-11/D-14/D-16/D-17/D-19/D-23~D-28 同一口径）。
"""
import os
from ..core.intio import read_lines
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-27"
TITLE = "压力供水管道管径选择"
AUTHOR = "张校正（新疆水利厅）"
HEAD_NAME = "D-27"

#: 圆周率（仅用于流速显示；见 docstring）
PI = 3.14159
#: 重力加速度 (m/s^2)
G = 9.81
#: 类型 6（硬聚氯乙烯 UPVC，《农村给水设计规范》公式）：由权威 OUT 反演的 (a,b,C)
#: （Q m³/s, d m）；与 D-28Intro 原表的 m/b 逐位相同
UPVC = (1.761, 4.761, 0.000918437)
#: 类型 6 在 D-28Intro 原表中登记的 f（名义值，见「未闭合 ②」）
UPVC_F_TABLE = 0.000875
#: 类型 6 在 D-27vb.EXE / D-28vb.exe 常量池中内嵌的源常量（名义值，见「未闭合 ②」）
UPVC_C_EXE = 0.0009120004252

#: 管壁类型表：code -> (名称, a, b, C, 标定状态)
#:   f/m/b 取自 D-28Intro「二、计算原理」原表（Q m³/s、d m；原表注记的 L/s 系文档笔误）；
#:   类型 6 的 C 取权威 OUT 反演值（原表值 0.000875，见「未闭合 ②」）。
#:   标定状态 "★OUT"=权威算例反演；"名义"=有原表出处但无权威算例；None=未实现
MATERIALS = {
    1: ("混凝土及钢筋混凝土管 n=0.013（谢才公式）", 2.0, 5.33, 0.0131, "名义"),
    2: ("混凝土及钢筋混凝土管 n=0.014（谢才公式）", 2.0, 5.33, 0.01553, "名义"),
    3: ("混凝土及钢筋混凝土管 n=0.014（谢才公式）", 2.0, 5.33, 0.01782, "名义"),
    4: ("旧钢管、旧铸铁管（安德列维舍夫公式）", 1.9, 5.1, 0.00625, "名义"),
    5: ("石棉水泥管（阿勃拉莫夫公式）", 1.85, 4.89, 0.001455, "名义"),
    6: ("硬聚氯乙烯UPVC(《农村给水设计规范》公式)", UPVC[0], UPVC[1], UPVC[2], "★OUT"),
    7: ("聚乙烯PE聚丙烯PP(《农村给水设计规范》公式)", 1.774, 4.774, 0.000951, "名义"),
    8: ("塑料硬管（谢维列夫公式）", 1.77, 4.77, 0.000946, "名义"),
    9: ("铝管及铝合金管（武汉水利电力学院公式）", 1.74, 4.74, 0.000861, "名义"),
    10: ("旧钢管/铸铁管按流速分段(《农村给水设计规范》)", None, None, None, None),
}


class MaterialError(ValueError):
    """管壁类型公式未标定 / 未实现。"""


def material_params(code):
    """返回 (a, b, C)；类型 10 抛 MaterialError。"""
    code = int(code)
    if code not in MATERIALS:
        raise MaterialError("管壁类型 %d 不在 1~10 之列" % code)
    name, a, b, c, tag = MATERIALS[code]
    if a is None:
        raise MaterialError(
            "管壁类型 %d（%s）未实现：D-28Intro 仅说明其存在（《农村给水设计规范》中"
            "旧钢管与铸铁管按流速分段计算每米水头损失的公式），未在说明书中给出 f/m/b；"
            "本内核实现 1~9 号管材。" % (code, name))
    return a, b, c


def loss_rate(code, q, d):
    """单米水头损失 i (m/m)。q：m³/s；d：m。"""
    a, b, c = material_params(code)
    return c * (q ** a) * (d ** (-b))


def area(d):
    """圆管断面积 A = π·d²/4 (m²)。"""
    return PI * d * d / 4.0


def velocity(q, d):
    """流速 V = Q/A (m/s)。"""
    return q / area(d)


def fmt_station(s):
    """桩号格式化：15+200.00（km+m，两位小数）。"""
    return "%d+%06.2f" % (int(s // 1000), s - 1000 * int(s // 1000))


# ============================================================
# 解析
# ============================================================

def parse(data):
    """
    解析输入。data：dict（原样返回）| .INT 文件路径。

    .INT 结构（D-27Intro「三、输入数据顺序」）：
      工程名,节点数
      (以下按节点输入，共「节点数」行)
        节点名称,桩号,设计流量(L/s),管壁类型,局部水头损失百分比,地面高程(m),地面处水头(m)
      数据个数,内径值(mm),内径值,...
    """
    if isinstance(data, dict):
        return data
    lines = [ln for ln in read_lines(data) if ln.strip()]
    if not lines:
        raise ValueError("D-27 输入为空")
    head = lines[0].split(",")
    name = head[0].strip()
    nnode = int(round(float(head[1])))
    nodes = []
    for ln in lines[1:1 + nnode]:
        f = [x.strip() for x in ln.split(",")]
        nodes.append({
            "name": f[0],
            "station": float(f[1]),
            "q": float(f[2]),          # L/s
            "mat": int(round(float(f[3]))),
            "loc": float(f[4]),        # 局部水头损失百分比（原著未参与计算，见 docstring）
            "gd": float(f[5]),
            "head": float(f[6]),
        })
    dia_line = [x.strip() for x in lines[1 + nnode].split(",") if x.strip()]
    ndia = int(round(float(dia_line[0])))
    dias = [float(x) for x in dia_line[1:1 + ndia]]
    return {"name": name, "nnode": nnode, "nodes": nodes, "diameters_mm": dias}


# ============================================================
# 计算
# ============================================================

def _pick_pair(dias_m, code, q, h, length):
    """在标准内径表中选出相邻管径对 (d_large, d_small, i_large, i_small)。

    规则：按内径升序扫描，找到第一个使 i(d)·L < H 的 d（该管径过大、损失不足），
    则 d_large = d（损失略小于 H）、d_small = 前一号（损失恰大于 H）。
    """
    loss_full = [loss_rate(code, q, d / 1000.0) * length for d in dias_m]
    for j in range(len(dias_m)):
        if loss_full[j] < h:
            if j == 0:
                return None
            dl, ds = dias_m[j], dias_m[j - 1]
            il = loss_rate(code, q, dl / 1000.0)
            is_ = loss_rate(code, q, ds / 1000.0)
            return dl, ds, il, is_
    return None


def compute(p):
    """执行 D-27 计算，返回节点成果表。"""
    nodes = p["nodes"]
    dias_m = p["diameters_mm"]

    def z(n):
        return n["gd"] + n["head"]

    rows = []          # {name, station, L, D, N, Z, GD, Q}
    spans = []         # 每跨度明细（供 JSON / 校核）

    # 起始行
    first = None
    for k in range(len(nodes) - 1):
        n0, n1 = nodes[k], nodes[k + 1]
        q = n1["q"] / 1000.0
        code = n1["mat"]
        h = z(n0) - z(n1)
        lt = n1["station"] - n0["station"]
        if h <= 0:
            raise ValueError("水压曲线有升高要求，请重新调整节点地面水头。")
        if lt <= 0:
            raise ValueError("桩号未递增（节点 %s → %s）" % (n0["name"], n1["name"]))
        picked = _pick_pair(dias_m, code, q, h, lt)
        if picked is None:
            raise ValueError("节点 %s→%s：内径资料无法覆盖所需水头（H=%.4f m, L=%.2f m）"
                             % (n0["name"], n1["name"], h, lt))
        dl, ds, il, is_ = picked
        L_large = (h - lt * is_) / (il - is_)
        if L_large < 0.0:
            L_large = 0.0
        if L_large > lt:
            L_large = lt
        L_small = lt - L_large
        spans.append({
            "i": k + 1, "从": n0["name"], "到": n1["name"], "Q": q, "管壁类型": code,
            "H": h, "L总": lt, "D大": dl / 1000.0, "D小": ds / 1000.0,
            "L大": L_large, "L小": L_small,
            "i大": il, "i小": is_,
        })
        if first is None:
            first = dl / 1000.0

    # 组装成果行
    z0 = z(nodes[0])
    rows.append({"name": nodes[0]["name"], "station": nodes[0]["station"],
                 "L": 0.0, "D": first, "N": nodes[0]["mat"], "Z": z0,
                 "GD": nodes[0]["gd"], "Q": nodes[0]["q"] / 1000.0})
    for k, sp in enumerate(spans):
        n0, n1 = nodes[k], nodes[k + 1]
        zk = z0 - sum(s["H"] for s in spans[:k])
        have_j = sp["L大"] > 1e-9
        if have_j:
            zj = zk - sp["i大"] * sp["L大"]
            gd_j = n0["gd"] + (n1["gd"] - n0["gd"]) * (sp["L大"] / sp["L总"])
            rows.append({"name": "J" + n0["name"], "station": n0["station"] + sp["L大"],
                         "L": sp["L大"], "D": sp["D大"], "N": n1["mat"], "Z": zj,
                         "GD": gd_j, "Q": sp["Q"]})
        rows.append({"name": n1["name"], "station": n1["station"],
                     "L": sp["L小"] if have_j else sp["L总"],
                     "D": sp["D小"] if have_j else sp["D大"],
                     "N": n1["mat"], "Z": z(n1),
                     "GD": n1["gd"], "Q": sp["Q"]})

    return {"程序": PROGRAM_ID, "标题": TITLE, "作者": AUTHOR, "输入": p,
            "起始桩号": nodes[0]["station"], "起始水位": z0,
            "分段数": len(rows) - 1, "跨度": spans, "成果": rows}


# ============================================================
# 输出（逐字复刻原著 .OUT 版式，自第 2 行起）
# ============================================================

_LINE = " " + "*" * 71
_TITLE = " ******             压力供水管道管径选择计算书 D-27               ******"
_SEP = " " + "-" * 85
_HDR1 = "         终止                           调整  管壁       水头    水压    地面    地面"
_HDR2 = " 节点名           流量    距离    内径  内径  类型 流速  损失    高程    高程    水头"
_HDR3 = "  I      桩号       Q       L       D           N    V    和       Z      GD      GY "


def _r2(x):
    """按两位小数四舍五入（与原著打印口径一致）。"""
    return float("%.2f" % x)


def _row_fmt(r):
    """数据行格式（D-27 成果表，85 字符）。

    GY（水头）按原著口径取**打印后**的 Z、GD 之差（Z−GD），
    与权威 OUT 的 J12/J14 行一致（若按全精度相减将差 0.01）。
    """
    return (" %-4s%9s %8.6f%8.2f %6.3f%8s%2d %5.2f%7.2f%8.2f%8.2f%8.2f"
            % (r["name"], fmt_station(r["station"]), r["Q"], r["L"], r["D"],
               "", r["N"], velocity(r["Q"], r["D"]),
               loss_rate(r["N"], r["Q"], r["D"]) * r["L"],
               r["Z"], r["GD"], _r2(r["Z"]) - _r2(r["GD"])))


def render(p, r):
    """生成原著风格文本计算书（.OUT）。返回自第 2 行（空行）起的文本。"""
    L = []
    L.append("")
    L.append(_LINE)
    L.append(_TITLE)
    L.append(_LINE)
    L.append("")
    L.append("                        工程名  %s" % p["name"])
    L.append("")
    L.append("                        管道分段数 M= %d " % r["分段数"])
    L.append("    起始桩号: %s                         起始水位  Zo=%.2f"
             % (fmt_station(r["起始桩号"]), r["起始水位"]))
    L.append("")
    L.append("                             分段计算结果")
    L.append("                             ============")
    L.append("")
    L.append(_HDR1)
    L.append(_HDR2)
    L.append(_HDR3)
    L.append("")
    L.append(_row_fmt(r["成果"][0]))
    for k, row in enumerate(r["成果"][1:]):
        if k % 2 == 0:
            L.append(_SEP)
        L.append(_row_fmt(row))
    L.append(_SEP)
    L.append("")
    return "\n".join(L)


def run(data, out_txt=None, out_json=None):
    """统一入口。data：INT 文件路径 | dict。"""
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
    _res, _txt = run(sys.argv[1] if len(sys.argv) > 1 else None)
    print(_txt)
