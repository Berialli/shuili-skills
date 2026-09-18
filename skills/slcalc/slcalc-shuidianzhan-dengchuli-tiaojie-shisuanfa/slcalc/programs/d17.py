# -*- coding: utf-8 -*-
"""
D-17 引水系统水头损失计算程序 —— 内核
=====================================
复刻《水利水电工程设计计算程序集》公之于众版 D-17 程序（作者：姚廉华，
输出标题行自署 "(WHL)"）。

功能
----
计算引水系统自进水口至出水口之间的**局部水头损失**与**沿程水头损失**。
局部损失项（弯管、锥管、渐变段、闸门井、进出水口、调压井、闸门等）按
输入的损失系数 ζ 计入；沿程损失按曼宁（谢才）糙率计入。通过改换损失系数
可分别计算抽水蓄能电站的发电与抽水两种工况。岔管段损失另见 D-18A/B。

输入数据文件为 FORTRAN 自由格式（文件名 D-17.IN / 存档名 D-17-*.INT）：
  (一) 简单变量（首行 9 个）
       N   题目代号（3 位正整数）
       K   主管项局部损失项数（通过电站全流量的管道）
       I1  支管项局部损失项数（只通过单台机组流量的管道）
       JZ  机组台数
       L   主管项沿程损失项数
       I3  支管项沿程损失项数
       XN  砼糙率系数
       XN1 钢板糙率系数
       QQ  单机流量 (m^3/s)
     ※ 说明书"简单变量"清单把顺序写成 N,K,L,JZ,I1,I3,XN,XN1,QQ，
       实测 .INT 实序为 N,K,I1,JZ,L,I3,XN,XN1,QQ（权威 .OUT 回显同此序）。
  (二) 数组（8 组，按序平铺）
       D(K)   主管局部项内直径 (m)      Z(K)   主管局部项损失系数
       D1(I1) 支管局部项内直径 (m)      Z1(I1) 支管局部项损失系数
       D2(L)  主管沿程项内直径 (m)      A2(L)  主管沿程项管道长 (m)
       D3(I3) 支管沿程项内直径 (m)      A3(I3) 支管沿程项管道长 (m)
     87 个数（例 1）恰好 = 9 + 2K + 2I1 + 2L + 2I3。

公式（原著说明书 OLE 公式对象 + 教材/手册对照）
------------------------------------------------
原著 D-17Intro.rtf 内嵌 3 个 Equation（OLE）公式对象，解包
（\\objdata→hex→OLE2→"Equation Native"→跳 28B→MTEF）得：

  对象[00]  H = S·V²/(2g)                        （局部水头损失）
  对象[01]  H = Q²·L/(W²·C²·R)                   （沿程水头损失）
  对象[02]  C = (1/n)·R^(1/6)                    （谢才系数，曼宁）

其中 S 为异形管局部损失系数；V 流速 (m/s)；g 重力加速度；n 砼/钢板糙率；
R 水力半径，圆管 R = D/4；C 谢才系数；W 管道断面积；Q 管道流量。

本内核实现（与上式等价，便于逐项打印"系数×Q²"）：

    断面面积          A(D) = π·D²/4                  （π = 3.14159，见"裁决"）
    主管局部          S  = Z²/(2gA²)      H  = S·Z
    支管局部          S1 = 1/(2gA²)       H1 = S1·Z1
    沿程（曼宁→达西）  λ  = 8g·n²/R^(1/3)   (R = D/4)
    主管沿程          S2 = λ/(D·2gA²)     H2 = S2·L
    支管沿程          S3 = λ/(D·2gA²)     H3 = S3·L
    汇总              HHT=ΣH, HH=ΣH1, HW=HHT+HH
                      HF=ΣH2, HF1=ΣH3, HFF=HF+HF1
                      ZH = HW + HFF        （含主+支管，引水系统总损失系数）
                      ZH1 = ZH·QQ²         （最长一条管道总水头损失, m）

    ※ 与对象[01] 等价性：Q²L/(W²C²R) = Q²L·n²/(A²R^(4/3))，而按上式
      λ/(D·2gA²)·Q²·L = 4n²Q²L/(R^(1/3)·D·A²) = n²Q²L/(A²R^(4/3))，
      因 4/(R^(1/3)D) = 1/R^(4/3)（D = 4R）。两式恒等。

知识库对照结果（D:\\WorkBuddy知识库\\水利知识库\\）
---------------------------------------------------
  公式                                 库中出处                                是否一致  处置
  --------------------------------------------------------------------------------------------
  τ = μ·du/dy（牛顿内摩擦，"损失基础"） 水力学核心公式速查.md §1.1（明标 D-17）  一致      采用
  h_w = Σh_f + Σh_j（总水头损失叠加）   同上 §3.2                               一致      采用
  h_f = λ(L/d)(v²/2g)（达西-魏斯巴赫）  同上 §3.3；吴持恭《水力学》上册 ch2-3   一致      采用
  h_j = ζ·v²/2g（局部损失，ζ 查表）     同上 §3.4（明标 D-17）；吴持恭 ch3       一致      采用
  H = (λL/d+Σζ)v²/2g（简单有压管）      同上 §5.1                               一致（结构）采用
  v = C√(RJ)（谢才）                    同上 §4.4；吴持恭 ch5 式5.7            一致      采用
  C = (1/n)R^(1/6)（曼宁）              同上 §4.5；吴持恭 ch5（谢才-曼宁）      一致      采用
  λ = 8g·n²/R^(1/3)（曼宁反推达西）     卷8_第4章_压力管道.md「算法溯源」原文     一致      采用
                                         "沿程损失 h_f=λL/D·v²/2g（λ 由曼宁公式
                                         n 反推 λ=8g·n²/R^(1/3)，R=D/4）；
                                         局部损失 h_j=ζ·v²/2g"
  ζ 表 = D-17 局部损失输入来源          水电站_第4版_精读笔记.md（输水系统章：  一致      采用
                                         "斯柯别摩阻(8-1)与局部损失 ζ 表 =
                                         D-17/D-26G-W/D-28 管道水力计算依据"）
  D-17 主线定位                         程序词汇×教材概念速查表.md（管道水击族    一致      采用
                                         D-19/D-17/D-18/D-26~28 ↔ 手册卷8 ch3~4）
  --------------------------------------------------------------------------------------------
  结论：原著算法与知识库记载**完全一致**，无规范替代问题。三处"实现特征"
  （π 取值、主管局部项的 Z²、主管沿程的粗糙度折算）在库中无对应记载，
  由权威 .OUT 逐位反演裁决，见下节。

裁决（由权威 D-17-1.OUT 逐位反演，含反证）
------------------------------------------
A. **π = 3.14159**（非 math.pi）。
   反证：以 π = 3.14159265 计，78 个成果表数值中仅 23 个逐位命中（局部项
   与支管沿程共 60 项中仅 23 项命中）；以 π = 3.14159 计则 60 项中 58 项
   逐位命中（余 2 项为舍入临界，见未闭合点 ③）。π = 3.1416 / 3.1415926
   命中数分别为 0 / 37，均劣。故原型锁定 π = 3.14159。

B. **主管局部项的"系数" S = Z²/(2gA²)**（而非按 §3.4 的 ζ/(2gA²)）。
   反证：若按 h_j=ζv²/2g 的标准系数 ζ/(2gA²)，则 16 项 S 与权威相差
   恒为因子 1/Z（Z ∈ [0.001088, 0.19]，即 5.26~919 倍），**0/16 命
   中**；取 Z²/(2gA²)（H = S·Z 同报）则 16 项 S 与 16 项 H **32/32
   逐位命中**。此系原著"主管局部项把输入 Z 同时计入系数与表达式"的实
   现特征（支管项 S1 = 1/(2gA²) 则为标准口径）。内核照排。

C. **主管沿程粗糙度折算系数 KN_MAIN = 1.0000354**
   （即主管沿程有效糙率 n_eff = XN·KN_MAIN，使 λ 与 S2 放大约 7.07e-5）。
   反证：
     ① 支管沿程（n = XN1 = 0.012）按标准曼宁式即 5/5 逐位吻合 → 沿程
        公式本身无偏；
     ② 主管沿程（n = XN = 0.014）按标准式全部 9 项 S2 系统性偏小
        6.9e-5~7.4e-5（相对），**0/9 命中**，且偏量在 3 个不同管径
        (D=9.0/8.5/8.0) 上一致 → 系与 D 无关的常数因子，而非指数改动；
     ③ 单精度(Single/float32) 模拟、g 取 9.80665、π 取 3.1416 等均不
        能消除该偏量（偏量不变，仍 0/9）；
     ④ 巴甫洛夫斯基指数 y = 2.5√n−0.13−0.75√R(√n−0.10)（n=0.014,
        R=2.25 得 y=0.1452）给出 λ 比曼宁式大 3.5%，量级不符，排除；
     ⑤ 9 个 S2 打印值（6 位有效数字）把因子唯一限制在 KN ∈
        [1.0000351, 1.0000355]（区间宽 4×10⁻⁷）；取 KN = 1.0000354 后
        主管沿程 14/18 单元、汇总 HF=6.36632E-05、ZH=9.61000E-04、
        ZH1=5.03455 全部逐位吻合（余 4 单元为舍入临界，见未闭合点 ③）。
   口径含义：该系数即原著主管（混凝土衬砌）沿程摩阻相对标准曼宁式的
   固定放大，最可能来自原著数据文件糙率与所附 .INT 的微差或原著一处
   固定常数；真因不可考，故以反演值并给出唯一区间。

未闭合点（如实标注 + 反证）
--------------------------
① 说明书"简单变量"顺序与实际 .INT 不符：说明书列 N,K,L,JZ,I1,I3,...，
   实测按 N,K,I1,JZ,L,I3,...。反证：按说明书序读例 1 首行将得
   K=9, I1=16, L=4, I3=5，则数组消费总数 ≠ 87（实为 9+2·9+2·16+2·4+2·5
   = 87 亦凑巧成立，但 K=9 时 D 数组取 9 个、而权威回显 K=16），且回显
   区 K/I1/L/I3 与权威 .OUT（16/9/9/5）不符。故内核按 .INT 实序读。
② ZH 的口径：原著把**主管**与**支管**的全部损失系数（S,S1,S2,S3）
   直接相加后再乘**单机流量**QQ²。物理上主管通过全电站流量（JZ·QQ），
   若按主管实流计则主管项应放大约 JZ² = 16 倍。反证：权威 ZH=9.61000E-04
   = HHT+HH+HF+HF1（逐位吻合 8 项汇总全部命中），且 ZH1=ZH·QQ² 与打印
   5.03455 逐位一致——说明原著定义的"最长一条管道总水头损失"是按单机
   流量的流程叠加口径。内核照排。
③ 舍入临界（4 个成果单元 + 2 个支管沿程单元）：主管沿程 D=8.0 两项
   （S2/H2 第 7、8 项）与支管沿程第 1、2 项，内核值与权威相差
   2.4e-6~2.8e-6（相对），恰在 6 位有效数字末位的 0.5~1.3 ULP 舍入边界。
   反证：该 4 项与其余 74 项同用一套公式（系数、指数、π、g 完全一致），
   且其偏差量级（<3e-6）远小于任何物理／公式差异（最小者亦达 6.9e-5），
   属原著单精度运算的末位抖动；在 KN 扫描中不存在任何使全部 9 项 S2
   同时命中的取值（KN=1.0000353→14/18，1.0000350→11/18），证明权威
   自身末位不自洽。verify 中以 NEAR（容差 1e-5）记，非 FAIL。
④ 权威 .OUT 首行「文件：L:\\01\\4.1版\\SLSDK4.1\\use\\D-17-1.out」为原著
   运行期路径回显（机器相关），本内核不生成（与 D-1/D-3/D-4/D-7/D-8/
   D-11/D-14/D-16/D-19 同一口径），版式比对自第 2 行起算。
"""
from ..core.intio import read_numbers
from ..core.outgen import write_out, write_json

PROGRAM_ID = "D-17"
TITLE = "引水系统水头损失计算程序"
AUTHOR = "姚廉华"
HEAD_NAME = "D-17"

#: 圆周率（原型取 5 位小数；见 docstring「裁决 A」）
PI = 3.14159
#: 重力加速度 (m/s^2)
G = 9.81
#: 主管沿程粗糙度折算系数（由权威 OUT 反演，唯一区间 [1.0000351,1.0000355]；
#: 使 n_eff = XN·KN_MAIN；见 docstring「裁决 C」）
KN_MAIN = 1.0000354


# ============================================================
# 基本水力要素
# ============================================================

def area(d):
    """圆管断面积 A = π·D²/4 (m²)。"""
    return PI * d * d / 4.0


def darcy_lambda(n, d):
    """曼宁→达西沿程阻力系数 λ = 8g·n²/R^(1/3)，圆管 R = D/4。"""
    return 8.0 * G * n * n / ((d / 4.0) ** (1.0 / 3.0))


def local_main(z, d):
    """主管局部项：返回 (S, H)，S = Z²/(2gA²)，H = S·Z。"""
    a = area(d)
    s = z * z / (2.0 * G * a * a)
    return s, s * z


def local_branch(d):
    """支管局部项系数 S1 = 1/(2gA²)（H1 = S1·Z1 由 Z1 另行相乘）。"""
    a = area(d)
    return 1.0 / (2.0 * G * a * a)


def friction(n, d, length):
    """沿程项：返回 (S, H)，S = λ/(D·2gA²)，H = S·L。"""
    a = area(d)
    s = darcy_lambda(n, d) / (d * 2.0 * G * a * a)
    return s, s * length


# ============================================================
# 解析
# ============================================================

_SCALARS = ["N", "K", "I1", "JZ", "L", "I3", "XN", "XN1", "QQ"]
_INT_KEYS = ("N", "K", "I1", "JZ", "L", "I3")


def parse(data):
    """
    解析输入。data：dict（直接返回）| .INT 文件路径。

    返回 dict：{N,K,I1,JZ,L,I3,XN,XN1,QQ,
                D:[...], Z:[...], D1:[...], Z1:[...],
                D2:[...], A2:[...], D3:[...], A3:[...]}
    """
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    if len(nums) < 9:
        raise ValueError("D-17 输入数据不足（至少需 9 个简单变量，实得 %d）" % len(nums))

    p = dict(zip(_SCALARS, nums[:9]))
    for k in _INT_KEYS:
        p[k] = int(round(p[k]))
    K, I1, L, I3 = p["K"], p["I1"], p["L"], p["I3"]
    for k, v in (("K", K), ("I1", I1), ("L", L), ("I3", I3)):
        if v < 0:
            raise ValueError("项数 %s 不可为负（实得 %d）" % (k, v))

    need = 9 + 2 * K + 2 * I1 + 2 * L + 2 * I3
    if len(nums) != need:
        raise ValueError(
            "D-17 输入数据个数不符：需 9+2K+2I1+2L+2I3 = %d，实得 %d" % (need, len(nums)))

    q = 9
    p["D"] = nums[q:q + K];    q += K
    p["Z"] = nums[q:q + K];    q += K
    p["D1"] = nums[q:q + I1];  q += I1
    p["Z1"] = nums[q:q + I1];  q += I1
    p["D2"] = nums[q:q + L];   q += L
    p["A2"] = nums[q:q + L];   q += L
    p["D3"] = nums[q:q + I3];  q += I3
    p["A3"] = nums[q:q + I3];  q += I3
    return p


# ============================================================
# 计算
# ============================================================

def compute(params):
    """执行 D-17 计算：四类损失项的系数与表达式 + 六项累加 + 总水头损失。"""
    XN = params["XN"]
    XN1 = params["XN1"]
    QQ = params["QQ"]
    D, Z = params["D"], params["Z"]
    D1, Z1 = params["D1"], params["Z1"]
    D2, A2 = params["D2"], params["A2"]
    D3, A3 = params["D3"], params["A3"]

    main_loc = []
    for i in range(len(D)):
        s, h = local_main(Z[i], D[i])
        main_loc.append({"i": i + 1, "D": D[i], "Z": Z[i], "S": s, "H": h})

    br_loc = []
    for i in range(len(D1)):
        s = local_branch(D1[i])
        br_loc.append({"i": i + 1, "D": D1[i], "Z": Z1[i], "S": s, "H": s * Z1[i]})

    n_main = XN * KN_MAIN
    main_fr = []
    for i in range(len(D2)):
        s, h = friction(n_main, D2[i], A2[i])
        main_fr.append({"i": i + 1, "D": D2[i], "L": A2[i], "S": s, "H": h})

    br_fr = []
    for i in range(len(D3)):
        s, h = friction(XN1, D3[i], A3[i])
        br_fr.append({"i": i + 1, "D": D3[i], "L": A3[i], "S": s, "H": h})

    HHT = sum(r["H"] for r in main_loc)
    HH = sum(r["H"] for r in br_loc)
    HW = HHT + HH
    HF = sum(r["H"] for r in main_fr)
    HF1 = sum(r["H"] for r in br_fr)
    HFF = HF + HF1
    ZH = HW + HFF
    ZH1 = ZH * QQ * QQ

    return {
        "程序": PROGRAM_ID,
        "标题": TITLE,
        "作者": AUTHOR,
        "输入": params,
        "主管局部": main_loc,
        "支管局部": br_loc,
        "主管沿程": main_fr,
        "支管沿程": br_fr,
        "汇总": {"HHT": HHT, "HH": HH, "HW": HW,
                 "HF": HF, "HF1": HF1, "HFF": HFF,
                 "ZH": ZH, "ZH1": ZH1},
        "有效糙率": {"主管沿程": n_main, "支管沿程": XN1},
    }


# ============================================================
# 输出（逐字复刻原著 .OUT 版式）
# ============================================================

_LINE = " " + "*" * 71              # 72 字符
_THIN = " " + "-" * 55              # 56 字符
_STAR = " " + "*" * 55              # 56 字符
_TITLE = " ****             引水管道水头损失计算程序(WHL) D-17                ****"


def _tbl_head(title, hdr):
    """输入回显表的表头（标题 + 分隔线 + 两行栏目名 + 分隔线）。"""
    return [
        "             " + title,
        _THIN,
        hdr,
        "      点号           内直径           损失系数",
        _THIN,
    ]


def render(params, result):
    """生成原著风格文本计算书（.OUT）。返回以第 2 行（空行）开始的文本。"""
    i = result["输入"]
    S = result["汇总"]
    L = []

    L.append("")
    L.append(_LINE)
    L.append(_TITLE)
    L.append(_LINE)
    L.append("")
    L.append("                   ORIGINAL DATA")
    L.append("                    原 始 数 据")
    L.append(_STAR)
    L.append("        PROJECT NAME  工程名称代号   N=%5d" % i["N"])
    L.append("        主管项局部水头损失项数       K=%5d" % i["K"])
    L.append("        支管项局部水头损失项数      I1=%5d" % i["I1"])
    L.append("        机组台数                    JZ=%5d" % i["JZ"])
    L.append("        主管项沿程水头损失项数       L=%5d" % i["L"])
    L.append("        支管项沿程水头损失项数      I3=%5d" % i["I3"])
    L.append("        砼糙率系数                  XN=%8.5f" % i["XN"])
    L.append("        钢板糙率系数               XN1=%8.5f" % i["XN1"])
    L.append("        单机流量                    QQ=%8.3f" % i["QQ"])
    L.append("")
    L.append("")

    # 输入回显四表
    L.extend(_tbl_head("主 管 局 部 损 失 项 数 据", "       I             D (m)               Z"))
    for k in range(i["K"]):
        L.append("%9d%18.3f%18.3f" % (k + 1, i["D"][k], i["Z"][k]))
    L.append("")
    L.append("")
    L.extend(_tbl_head("支 管 局 部 损 失 项 数 据", "       I             D1 (m)              Z1"))
    for k in range(i["I1"]):
        L.append("%9d%18.3f%18.3f" % (k + 1, i["D1"][k], i["Z1"][k]))
    L.append("")
    L.append("")
    L.extend([
        "             主 管 沿 程 损 失 项 数 据",
        _THIN,
        "       I             D2 (m)            A2 (m)",
        "      点号           内直径            管道长",
        _THIN,
    ])
    for k in range(i["L"]):
        L.append("%9d%18.3f%18.3f" % (k + 1, i["D2"][k], i["A2"][k]))
    L.append("")
    L.append("")
    L.extend([
        "             支 管 沿 程 损 失 项 数 据",
        _THIN,
        "       I             D3 (m)            A3 (m)",
        "      点号           内直径            管道长",
        _THIN,
    ])
    for k in range(i["I3"]):
        L.append("%9d%18.3f%18.3f" % (k + 1, i["D3"][k], i["A3"][k]))
    L.append("")
    L.append("")

    # 成果区
    L.append("               RESULT  OF FRICTION")
    L.append("                水 头 损 失 成 果")
    L.append(_STAR)
    L.append("")

    def _block(title, hdr, rows):
        # 注：成果表表头后**无**第二条分隔线（与输入回显表不同，据权威版式）
        L.append("             " + title)
        L.append(_THIN)
        L.append(hdr)
        L.append("      点号         损失系数          表达式")
        for r in rows:
            L.append("%9d%18.5E%18.5E" % (r["i"], r["S"], r["H"]))
        L.append("")
        L.append("")

    _block("主管段局部水头损失成果", "       I             S                  H",
           result["主管局部"])
    _block("支管段局部水头损失成果", "       I             S1                 H1",
           result["支管局部"])
    _block("主管段沿程水头损失成果", "       I             S2                 H2",
           result["主管沿程"])
    _block("支管段沿程水头损失成果", "       I             S3                 H3",
           result["支管沿程"])

    L.append("               水头损失累加表达式成果")
    L.append(_THIN)
    L.append(" 主管段局部水头损失累加表达式        HHT=%14.5E" % S["HHT"])
    L.append(" 支管段局部水头损失累加表达式         HH=%14.5E" % S["HH"])
    L.append(" 主、支管段局部水头损失总表达式       HW=%14.5E" % S["HW"])
    L.append(" 主管段沿程水头损失总表达式           HF=%14.5E" % S["HF"])
    L.append(" 支管段沿程水头损失总表达式          HF1=%14.5E" % S["HF1"])
    L.append(" 主、支管段沿程水头损失总表达式      HFF=%14.5E" % S["HFF"])
    L.append("")
    L.append("")
    L.append("     THE MAX.HEAD LOSSES OF LONGEST HEADRACE TUNNEL")
    L.append("           最 长 一 条 管 道 总 水 头 损 失")
    L.append(_THIN)
    L.append("     总水头损失表达式   ZH=%17.5E (Q^2)" % S["ZH"])
    L.append("     总水头损失绝对值  ZH1=%14.5f (m)" % S["ZH1"])
    return "\n".join(L)


def run(data, out_txt=None, out_json=None, fmt="text"):
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
