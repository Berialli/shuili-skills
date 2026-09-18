# -*- coding: utf-8 -*-
"""
D-4 各种堰流水力学计算程序 —— 内核
===================================
复刻《水利程序集》D-4 程序（作者：陈靖齐，水电部天津勘测设计院；
数据来源：《水利水电工程设计计算程序集》公之于众版，
乌鲁木齐正海水利科技有限公司，张校正教授级高工技术总负责）。

功能（原著五、程序功能）：
  校核流量 Q；设计水头 Hd；计算溢流宽度 B。用字符变量 A$="Q"/"H"/"B" 控制。
  计算 Hd 用逐次迭代法，|ΔH| < 0.001。

三类堰（TY，判据 δ/H）：
  TY=1 薄壁堰  δ/H < 0.67   ……巴赞(Bazin)流量系数
  TY=2 实用堰  0.67 < δ/H < 2.5 ……WES 剖面
  TY=3 宽顶堰  2.5  < δ/H < 10

────────────────────────────────────────────────────────────────────────
算法与数据来源（本次全部从原著 EXE 二进制还原，非教科书替换）
────────────────────────────────────────────────────────────────────────
原著说明书 D-4Intro.rtf / D-4GIntro.rtf 中 8 个公式为 OLE(Equation.3) 对象，
文本层缺失。本次把 8 个 \\objdata 解包为 OLE2 →“Equation Native”→ 跳过 28 字节
EQNOLEFILEHDR → MTEF v3 记录流，逐记录还原字符（见 _d4_mtef_probe.py），得到
公式的**真实字面形式**；同时 D-4Gvb.EXE 为原生 x86，其公式文本以 UTF-16 常量
存放（见 _d4_exe_probe.py），数据库数组以

    c7 85 2c fa ff ff <idx:4>      mov dword [ebp-0x5d4], idx
    c7 04 (81|90|8a) <val:4>       mov dword [reg+reg*4], Single(val)

循环初始化，本次把全部 264 个数据库条目逐条提取（见 _d4_db_extract.py），
得到原著的 σ 表 / m 表 / ξo 表（下表中括号内为 DB 下标，可溯源）。

一、薄壁堰（TY=1）（EXE 常量字面）
  流量公式      Q = m₀·b·√(2g)·H^(3/2)                     （打印字面 H）
  流量系数(巴赞) m₀ = (0.405 + 0.0027/H)·(1 + 0.55·(H/(H+P))²)
  考虑侧收缩时   m  = (0.405 + 0.0027/H − 0.03·(1 − b/B))·
                       (1 + 0.55·(H/(H+P))²·(b/B)²)
  适用范围 H=0.1~0.6 m；q=0.2~2.0；H ≤ 2P。不包含淹没问题。
  TZ：矩形堰 1、三角堰 2（三角堰输入内角）。

二、宽顶堰（TY=3）
  流量公式  Q = σ·ε·m·B·√(2g)·Ho^(3/2)，Ho = H + Vo²/2g
  淹没系数  σ = f(hs/Ho)，19 点表 hs/Ho=0.80~0.98（DB 153~171 横坐标、173~191 纵坐标）：
     0.80→1.00 0.81→0.995 0.82→0.99 0.83→0.98 0.84→0.97 0.85→0.96 0.86→0.95
     0.87→0.93 0.88→0.90 0.89→0.87 0.90→0.84 0.91→0.82 0.92→0.78 0.93→0.74
     0.94→0.70 0.95→0.65 0.96→0.59 0.97→0.50 0.98→0.40；hs/Ho<0.80 时 σ=1。
  侧收缩系数 ε = 1 − 0.2·(ξk + (n−1)·ξo)·Ho/(n·b)，b = B/n 为单孔宽
     边墩 ξk：直角1→1.0、八字(折角)2→0.7、圆弧3→0.7、流线型4→0.4
     闸墩 ξo = f(SP, hs/Ho)（DB 235~264）：
       横坐标 hs/Ho = 0.75, 0.80, 0.85, 0.90, 0.95（低于 0.75 取首列）
       SP=1 矩形   [0.80, 0.86, 0.92, 0.98, 1.00]
       SP=2 尖角形 [0.45, 0.51, 0.57, 0.63, 0.69]
       SP=3 半圆形 [0.45, 0.51, 0.57, 0.63, 0.69]
       SP=4 尖圆形 [0.25, 0.32, 0.39, 0.46, 0.53]
       SP=5 流线型 [0, 0, 0, 0, 0]
       （n=1 时 ξo 不参与，程序按 0 显示）
  流量系数 m：
     MI=1 直坎  m = 0.32 + 0.01·(3 − P/H)/(0.46 + 0.75·P/H)
     MI=2 圆坎  m = 0.36 + 0.01·(3 − P/H)/(1.2 + 1.5·P/H)
     MI=3 无坎宽顶堰（m 已含翼墙影响，此时 ε 不计 ξk）：
        m = f(B/Bo)（DB 1~11 横坐标 0~1.0 步长 0.1）：
          直角翼墙(以及 r/b=0)  [0.320 0.322 0.324 0.327 0.330 0.334 0.340 0.346 0.355 0.367 0.385]
          八字 ctgθ=0          [0.320 0.322 0.324 0.327 0.330 0.334 0.340 0.346 0.355 0.367 0.385]
          八字 ctgθ=0.5        [0.343 0.344 0.346 0.348 0.350 0.352 0.356 0.360 0.365 0.373 0.385]
          八字 ctgθ=1.0        [0.350 0.351 0.352 0.354 0.356 0.358 0.361 0.364 0.369 0.375 0.385]
          八字 ctgθ=2.0        [0.353 0.354 0.355 0.357 0.358 0.360 0.363 0.366 0.370 0.376 0.385]
          圆角 r/b=0.2         [0.349 0.350 0.351 0.353 0.355 0.357 0.360 0.363 0.368 0.375 0.385]
          圆角 r/b=0.3         [0.354 0.355 0.356 0.357 0.359 0.361 0.363 0.366 0.371 0.376 0.385]
          圆角 r/b=0.5         [0.360 0.361 0.362 0.363 0.364 0.366 0.368 0.370 0.373 0.378 0.385]
  求 B 时由流量公式直接反解（原著五(五)）：B = Q/(σ·m·√(2g)·Ho^1.5) + 0.2·(ξk+(n−1)ξo)·Ho

三、实用堰（TY=2，WES 剖面）
  堰面 y/Hd = 0.5·(x/Hd)^1.85 (x≥0)；上游三圆弧 R1=0.5Hd、R2=0.2Hd、R3=0.04Hd，
     b1=0.175Hd、b2=0.276Hd、b3=0.2818Hd。
  流量公式 Q = σ·ε·m·B·√(2g)·Ho^(3/2)   （EXE 打印字面误作 Ho^(2/3)，见下“裁决”）
    m = f(Ho/Hd)（DB 119~135，横坐标 DB 136~152 = 0~1.6 步长 0.1）：
      0.0→0.385 0.1→0.401 0.2→0.416 0.3→0.430 0.4→0.442 0.5→0.455 0.6→0.466
      0.7→0.476 0.8→0.485 0.9→0.494 1.0→0.502 1.1→0.510 1.2→0.516 1.3→0.522
      1.4→0.530 1.5→0.535 1.6→0.542
    σ = f(hs/Ho)（DB 214~233 纵坐标，横坐标 0.05~1.00 步长 0.05）：
      0.05→1.000 0.10→0.996 0.15→0.991 0.20→0.986 0.25→0.981 0.30→0.976
      0.35→0.970 0.40→0.963 0.45→0.956 0.50→0.948 0.55→0.937 0.60→0.923
      0.65→0.907 0.70→0.886 0.75→0.856 0.80→0.821 0.85→0.778 0.90→0.709
      0.95→0.621 1.00→0.438；hs/Ho<0.05 时 σ=1。
    ε 同宽顶堰（原著“计算公式同宽顶堰，并用同一段程序”）。
  设计水头 Hd：A$="H" 时逐次迭代 |ΔH|<0.001。

────────────────────────────────────────────────────────────────────────
裁决项（正文印刷与实现不一致，已实测锁定）
────────────────────────────────────────────────────────────────────────
  (A) 直坎流量系数基数：EXE/OUT 打印字面为「m=0.33+0.01*(3-P/H)/(0.46+0.75*P/H)」，
      但用权威 D-4G-1.OUT 的 m=0.3625 反算，P/H=0.25 时
        0.33 基数 → 0.372471（≠0.3625，差 +0.0100）
        0.32 基数 → 0.362471（=0.3625，逐位命中）
      故实现取 0.32，打印字面 0.33 属原著文字笔误。
  (B) 实用堰流量公式指数：EXE 打印字面为「Q=…*B*SQR(2*g)*Ho^(2/3)」，
      宽顶堰同一位置为 Ho^(3/2)。堰流为 3/2 次方，2/3 系印刷笔误；
      本内核按 Ho^(3/2) 实现（另见 d4_verify.py 反证③）。
  (C) D-4G-1.OUT 的“四、计算结果”段被程序打印了两遍（同一块重复），
      本内核只输出一次（不影响数值对拍）。

回补 G1~G3（2026-09-11，D-4G 身份判定后对实现颗粒度的修正）
────────────────────────────────────────────────────────────────────────
  G1  INT 第 3 字段语义：由「MJ（翼墙形式号）」更正为「**计算内容模式码**」
      （1=求Q / 2=求H / 3=求B，即原著§五的字符变量 A$）。
      · INT 实序反证：D-4G-10~13 的第 3 字段依次 2/1/3/1，与说明书
        「求水头 / 求流量 / 求堰宽 / 求流量」逐例一一对应；
        D-4G-1~9 同字段依次 1/3/3/3/3，与「例1 求流量、例2~5 求堰宽」亦逐例吻合。
      · 双路反证：§五 程序功能原文「用字符变量 A＄="Q"，"H"，"B"控制」，
        且 D-4Gvb.EXE UI 串表中「计算内容」与「至流量/至水头/至堰宽」并列出现。
      · 向后兼容：D-4 原 9 例中该字段恒与「例1 求Q」一致（值=1，与旧 MJ=1 数值巧合），
        故例 1 对权威 D-4G-1.OUT 的对拍结果不变（EXACT 9/10 + NEAR 1/10）。
      · 副作用：旧实现把该字段当 MJ，仅 MI=3（无坎）时才被使用，而 MI=1/2 时 MJ 本就被忽略，
        故对例 1 无影响；MI=3 的 INT 直读一支本属 U1 不可唯一反演（详见下），改后由
        「紧凑格式」分支（D-4G-10~13）提供可靠通路。

  G2  补 calc_wide() 的 MODE="H"（求水头）分支。旧实现只支持 "Q"/"B"，
      显式传 MODE="H" 会落入 "B" 分支并抛 KeyError('H')（D-4G 例10 即此支）。
      现以 Q=σ·ε·m·B·√(2g)·Ho^1.5 对 H 二分反解，收敛口径沿用原著§五
      「逐次迭代法，│△H│＜0.001」（tol_H=0.001）。
      对拍：例10（B0=6,n=1,B=2,Q=3.39,hs=0.82,X1=3,r/b=0.1,B/Bo=0.333333）
      说明书参考解 H=1.06 → 内核 1.0663（+0.59%，见 G3）。

  G3  无坎宽顶堰 m 表 r/b=0.1 的插值口径 **定案为「r/b 维线性内插」**。
      原著数据库该表横坐标为 DB71~74 = [0.0, 0.2, 0.3, 0.5]（EXE 逐条复抽确认，
      非转抄误差），而 D-4G-10~13 均给 r/B=0.1（落在 0.0 与 0.2 之间）。
      旧实现 min(RB_GRID, key=|v-rb|) 在 0.1 处**并列**、取到 0.0 行（=直角行）。
      定案依据（说明书参考答案为唯一基准；全量穷举见 slcalc/_d4g_g3_search.txt）：

        m 口径（均不计 ξk，与原著「无坎…不计ξk」一致）    例10    例11    例12    例13    最大
        ─────────────────────────────────────────────────────────────────────────────
        ① r/b 线性内插 (0.0,0.2) 中值  ★采用             +0.59%  +0.88%  +0.83%  +4.11%  4.11%
        ② 八字翼墙 ctgθ=0.5 行                           +0.98%  +1.49%  +1.39%  −2.01%  2.01%
        ③ 最近行取 0.2                                    +1.96%  +3.01%  +2.76%  −0.49%  3.01%
        ④ 最近行取 0.0（改前现状）                        +3.30%  +4.76%  +4.65%  −7.69%  7.69%
        ─────────────────────────────────────────────────────────────────────────────
        另：若把 ξk 计入侧收缩（X1=3→ξk=0.7），例10/例12 无解（ε 过小），
            穷举全部 11 行均发散 → **反证** d04.py 原有的 include_xk=False 正确。

      · 采用①的理由：(a) 语义正确——原著正文「本程序已把他们存入数据库中，可供插值用」，
        0.1 位于网格 0.0 与 0.2 之间，内插是唯一非外推口径；(b) 残差最优——① 对
        **同一物理工况的三个互相独立未知量**（例10 求H / 例11 求Q / 例12 求B）同时给出
        +0.59%/+0.88%/+0.83% 的自洽三角，而③④对同一三角分别偏 +2~3%/+3~5%，
        呈系统性偏差；(c) ②③ 在语义上均不成立（说明书明写「无坎宽顶堰园角形翼墙之r/B=0.1」，
        非八字翼墙）。
      · **未闭合项（如实标注）**：例13（B0=8.8,n=2,B=4,hs=0.82,SP=3,r/b=0.1）在①下得
        Q=6.7316，说明书 7.02，偏 **−4.11%**。反证：要让例13 在①下命中 7.02，需 ε=1.0140>1
        （物理不可能），或 m 需取 r/b≈0.25 行——但 r/b=0.25 会使例11 偏 +3.5%，与①的
        +0.88% 冲突。**该 4 例无法由任一行选择同时满足**（取②行时例13 命中 −2.0% 而例11 偏 +1.5%），
        故判定为**原著说明书例13 自身（或其数据录入）与例10~12 不属同一 wing-wall 取值**，
        非本内核公式错误。列 DECL，量化残差 −4.11%。

未闭合点（如实标注）
────────────────────
  U1. INT 数据文件的字段顺序随 A$（求 Q / 求 H / 求 B）与 MI（直坎/圆坎/无坎）
      而变。**第 3 字段已定案为「计算内容模式码」（见 G1）**，其余槽位：
        · 紧凑格式（逗号单行，D-4G-10~13，TY=3）——**已唯一锁定，见 _parse_compact_wide()**：
            求Q: [TY, MI, 1, B0, n, B, H, Hd, hs, X1, (SP|n>1), r/b, B/Bo]
            求H: [TY, MI, 2, B0, n, B, Q, hs, X1, (SP|n>1), r/b, B/Bo]
            求B: [TY, MI, 3, B0, n, Q, H, Hd, hs, X1, (SP|n>1), r/b, B/Bo]
          四例代入后字段数 11/12/12/13 与文件字节数逐例吻合（例13 因 n=2 多一 SP 槽）。
        · 旧式（CRLF 逐行，D-4-1~9 ≡ D-4G-1~9）——仅「求 Q、MI=1/2」一支可用权威
          D-4G-1.OUT 完整锁定：[TY, MI, MODE, Bo, n, B, P, H, Hd, hs, X1, SP, …]
          例2/4/5（求 B）、例3（MI=3）末段字段存在 ±1~2 位漂移
          （例3 的 hs=5.42 落在第 9 位而例4 的 hs=3.7 落在第 10 位），
          该支仍**不可唯一反演**；内核提供 dict 规范化入口规避，详见 parse()。
  U2. 实用堰 σ 表横坐标起点：DB 中 σ 值 20 个（0.05~1.00），而其前另有
      2 个 0.0 常量（DB 192、193），无法唯一判定横坐标是 0.00~1.00(21 点)
      还是 0.05~1.00(20 点)。本内核取后者并令 hs/Ho<0.05 时 σ=1；
      该歧义不影响例 6/7（二者 hs/Ho<0，σ≡1），属不可唯一反演点。
  U3. 圆坎/直坎 m 公式的代码常数（0.32/0.36/0.01/0.46/0.75/1.2/1.5）在 EXE 中
      未以 float/double 明文出现（疑由 VB6 编译期常量折叠或运行时按 Double 计算），
      其中 0.32 由例 1 权威值锁定（见裁决 A），0.36 由例 2/4 书值支持（偏差 +0.3~0.4%）。
  U4. 薄壁堰三角堰（TZ=2）无权威算例（D-4-9.INT 仅 4 个字段且无 OUT），
      本内核按常用公式 Q=(8/15)·μ·√(2g)·tan(θ/2)·H^2.5（μ=0.6）实现并标注。

验证
────
  · 例 1（TY=3、直坎、求 Q）对权威 D-4G-1.OUT 逐位对拍：Vo/H0/hs·Ho⁻¹/σ/ξk/ξo/m/B/ε/Q
    共 10 项全部 4 位小数相同（EXACT 10/10 → 实测 EXACT 9 + NEAR 1）。
  · 例 2/5/8 对说明书参考答案 PASS（偏差 ≤0.4%）；例 3/4/6/7 见 d4_verify.py 逐条判定。
  · 例 10~13（D-4G 增补的无坎宽顶堰 求H/求Q/求B/求Q，基准=说明书参考答案）：
    经 G1~G3 回补后由 INT 紧凑格式直读并计算，统计见 d4_verify.py ③ 节。
"""
import math

from ..core.intio import read_numbers
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "D-4"
TITLE = "堰流水力学计算书"
AUTHOR = "陈靖齐（水电部天津勘测设计院）"

G = 9.8                     # 原著重力加速度（D-2/D-3 同源，取 9.8 m/s²）
SIGMA_FLOOR = 1.0           # 未淹没时 σ=1

# ── 堰型 / 型号名 ────────────────────────────────────────────
TY_NAMES = {1: "薄壁堰", 2: "实用堰", 3: "宽顶堰"}
MI_NAMES = {1: "直坎", 2: "圆坎", 3: "无坎宽顶堰"}
TZ_NAMES = {1: "矩形堰", 2: "三角堰"}
MJ_NAMES = {1: "直角翼墙", 2: "八字翼墙", 3: "圆弧翼墙"}
X1_NAMES = {1: "直角形", 2: "八字形(折角)", 3: "圆弧形", 4: "流线形"}
SP_NAMES = {1: "矩形", 2: "尖角形", 3: "半圆形", 4: "尖圆形", 5: "流线型"}

# ── 原著数据库（逐条抄自 D-4Gvb.EXE，括号为 DB 下标）─────────
BB0_GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]          # DB 1~11

M_ZHI_WING = [0.320, 0.322, 0.324, 0.327, 0.330, 0.334,
              0.340, 0.346, 0.355, 0.367, 0.385]                            # DB 12~22 直角
CTX_GRID = [0.0, 0.5, 1.0, 2.0]                                            # DB 23~26
M_CTG = {                                                                  # DB 27~70
    0.0: [0.320, 0.322, 0.324, 0.327, 0.330, 0.334, 0.340, 0.346, 0.355, 0.367, 0.385],
    0.5: [0.343, 0.344, 0.346, 0.348, 0.350, 0.352, 0.356, 0.360, 0.365, 0.373, 0.385],
    1.0: [0.350, 0.351, 0.352, 0.354, 0.356, 0.358, 0.361, 0.364, 0.369, 0.375, 0.385],
    2.0: [0.353, 0.354, 0.355, 0.357, 0.358, 0.360, 0.363, 0.366, 0.370, 0.376, 0.385],
}
RB_GRID = [0.0, 0.2, 0.3, 0.5]                                             # DB 71~74
M_RB = {                                                                   # DB 75~118
    0.0: [0.320, 0.322, 0.324, 0.327, 0.330, 0.334, 0.340, 0.346, 0.355, 0.367, 0.385],
    0.2: [0.349, 0.350, 0.351, 0.353, 0.355, 0.357, 0.360, 0.363, 0.368, 0.375, 0.385],
    0.3: [0.354, 0.355, 0.356, 0.357, 0.359, 0.361, 0.363, 0.366, 0.371, 0.376, 0.385],
    0.5: [0.360, 0.361, 0.362, 0.363, 0.364, 0.366, 0.368, 0.370, 0.373, 0.378, 0.385],
}

WES_HOHD = [i / 10.0 for i in range(17)]                                   # DB 136~152
WES_M = [0.385, 0.401, 0.416, 0.430, 0.442, 0.455, 0.466, 0.476, 0.485,    # DB 119~135
         0.494, 0.502, 0.510, 0.516, 0.522, 0.530, 0.535, 0.542]

KDR_SG_X = [0.80 + 0.01 * i for i in range(19)]                            # DB 153~171
KDR_SG = [1.00, 0.995, 0.990, 0.980, 0.970, 0.960, 0.950, 0.930, 0.900,    # DB 173~191
          0.870, 0.840, 0.820, 0.780, 0.740, 0.700, 0.650, 0.590, 0.500, 0.400]

WES_SG_X = [0.05 + 0.05 * i for i in range(20)]                            # DB 194~213
WES_SG = [1.000, 0.996, 0.991, 0.986, 0.981, 0.976, 0.970, 0.963, 0.956,   # DB 214~233
          0.948, 0.937, 0.923, 0.907, 0.886, 0.856, 0.821, 0.778, 0.709,
          0.621, 0.438]

XK_TABLE = {1: 1.0, 2: 0.7, 3: 0.7, 4: 0.4}                                # ξk 边墩
XO_X = [0.75, 0.80, 0.85, 0.90, 0.95]                                      # DB 235~239
XO_TABLE = {                                                               # DB 240~264
    1: [0.80, 0.86, 0.92, 0.98, 1.00],     # 矩形
    2: [0.45, 0.51, 0.57, 0.63, 0.69],     # 尖角形
    3: [0.45, 0.51, 0.57, 0.63, 0.69],     # 半圆形
    4: [0.25, 0.32, 0.39, 0.46, 0.53],     # 尖圆形
    5: [0.0, 0.0, 0.0, 0.0, 0.0],          # 流线型
}

TZ_VNOTCH_MU = 0.6        # 三角堰流量系数（无权威算例，见 U4）


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


# ── 系数 ─────────────────────────────────────────────────────

def m_zhikan(P, H):
    """直坎流量系数（基数 0.32，见裁决 A）。"""
    r = P / H if H else 0.0
    return 0.32 + 0.01 * (3.0 - r) / (0.46 + 0.75 * r)


def m_yuankan(P, H):
    """圆坎流量系数（基数 0.36）。"""
    r = P / H if H else 0.0
    return 0.36 + 0.01 * (3.0 - r) / (1.2 + 1.5 * r)


def _rb_bracket(rb):
    """
    G3 定案：r/b 维线性内插的括入行与权重（端点外沿用首/末行）。

    原著 DB71~74 = [0.0, 0.2, 0.3, 0.5]（EXE 逐条复抽确认）。D-4G-10~13 的
    r/B=0.1 落在 0.0 与 0.2 之间，旧实现按「最近行」在 0.1 处并列取到 0.0 行
    （=直角行），对例10/11/12 呈 +3.3~4.8% 的系统性偏差；改为线性内插后
    同一工况的三个未知量（求H/求Q/求B）残差收敛到 +0.6~0.9% 的自洽三角。
    定案穷举见 slcalc/_d4g_g3_search.txt（①行），未闭合的例13 见 docstring G3。
    """
    if rb <= RB_GRID[0]:
        return RB_GRID[0], RB_GRID[0], 0.0
    if rb >= RB_GRID[-1]:
        return RB_GRID[-1], RB_GRID[-1], 0.0
    for i in range(len(RB_GRID) - 1):
        lo, hi = RB_GRID[i], RB_GRID[i + 1]
        if lo <= rb <= hi:
            return lo, hi, (rb - lo) / (hi - lo)
    return RB_GRID[-1], RB_GRID[-1], 0.0


def m_wukan(MJ, bb0, ctg, rb):
    """无坎宽顶堰流量系数（m 已含翼墙影响）：MJ=1 直角 / 2 八字(ctgθ) / 3 圆弧(圆角形, r/b)。"""
    if MJ == 1:
        return interp(BB0_GRID, M_ZHI_WING, bb0)
    if MJ == 2:
        c = min(CTX_GRID, key=lambda v: abs(v - ctg))
        return interp(BB0_GRID, M_CTG[c], bb0)
    lo, hi, w = _rb_bracket(rb)                     # G3：r/b 维线性内插
    if w == 0.0:
        return interp(BB0_GRID, M_RB[lo], bb0)
    row = [M_RB[lo][i] * (1.0 - w) + M_RB[hi][i] * w for i in range(len(BB0_GRID))]
    return interp(BB0_GRID, row, bb0)


def sigma_kdr(hs, H0):
    """宽顶堰淹没系数 σ = f(hs/Ho)，hs/Ho<0.80 取 1。"""
    if H0 <= 0:
        return 1.0
    x = hs / H0
    if x < KDR_SG_X[0]:
        return 1.0
    return interp(KDR_SG_X, KDR_SG, x)


def sigma_wes(hs, H0):
    """实用堰淹没系数 σ = f(hs/Ho)，hs/Ho<0.05 取 1。"""
    if H0 <= 0:
        return 1.0
    x = hs / H0
    if x < WES_SG_X[0]:
        return 1.0
    return interp(WES_SG_X, WES_SG, x)


def xk_of(X1):
    return XK_TABLE.get(int(X1 or 0), 0.7)


def xo_of(SP, hs, H0, n):
    """闸墩形状系数 ξo；n<=1 时不参与，程序按 0 显示。"""
    if n is None or n <= 1:
        return 0.0
    row = XO_TABLE.get(int(SP or 0), XO_TABLE[1])
    x = hs / H0 if H0 > 0 else 0.0
    return interp(XO_X, row, x)


def eps_of(n, b, H0, xk, xo, include_xk=True):
    """侧收缩系数 ε = 1 − 0.2·(ξk+(n−1)ξo)·Ho/(n·b)。无坎时不计 ξk。"""
    if n <= 0 or b <= 0:
        return 1.0
    c = ((xk if include_xk else 0.0) + (n - 1) * xo)
    return 1.0 - 0.2 * c * H0 / (n * b)


def velocity_head(Q, B0, P, H, g=G):
    """行进流速水头 Vo²/2g（B0<=0 时程序取 Vo=0）。"""
    if not B0 or B0 <= 0:
        return 0.0, 0.0
    A = B0 * (P + H)
    if A <= 0:
        return 0.0, 0.0
    v0 = Q / A
    return v0, v0 * v0 / (2.0 * g)


# ── 各堰型计算 ───────────────────────────────────────────────

def _wide_m(mi, P, H, MJ, bb0, ctg, rb):
    if mi == 1:
        return m_zhikan(P, H)
    if mi == 2:
        return m_yuankan(P, H)
    return m_wukan(MJ, bb0, ctg, rb)


def calc_wide(params, mode="Q", g=G, tol=1e-9, max_iter=200, tol_H=0.001):
    """宽顶堰（TY=3）：求 Q / 求 H（G2 新增）/ 求 B。返回结果 dict。"""
    mi = int(params["MI"])
    P = params["P"]
    H = params.get("H", 0.0)          # MODE="H"（求水头）时未知，由 H 分支反解后回填
    hs = params.get("HS", 0.0)
    B0 = params.get("B0", 0.0)
    n = int(params.get("N", 1) or 1)
    X1 = params.get("X1", 0)
    SP = params.get("SP", 1)
    MJ = int(params.get("MJ", 1) or 1)
    ctg = params.get("CTG", 0.0)
    rb = params.get("RB", 0.0)
    xk = xk_of(X1)
    include_xk = (mi != 3)

    if mode == "Q":
        B = params["B"]
        b = B / n if n else B
        xo = xo_of(SP, hs, H + 1e-12, n)
        Q, v0, H0 = 0.0, 0.0, H
        for _ in range(max_iter):
            v0, hv = velocity_head(Q, B0, P, H, g)
            H0 = H + hv
            sig = sigma_kdr(hs, H0)
            xo = xo_of(SP, hs, H0, n)
            ep = eps_of(n, b, H0, xk, xo, include_xk)
            bb0 = (B / B0) if B0 and B0 > 0 else 0.0
            m = _wide_m(mi, P, H, MJ, bb0, ctg, rb)
            Qn = sig * ep * m * B * math.sqrt(2 * g) * H0 ** 1.5
            if abs(Qn - Q) < tol:
                Q = Qn
                break
            Q = Qn
        return dict(TY=3, MI=mi, mode=mode, Vo=v0, H0=H0, hsH0=(hs / H0 if H0 else 0),
                    SG=sigma_kdr(hs, H0), XK=xk, X0=xo, M=m,
                    EP=ep, B=B, Q=Q, b=b, n=n)

    if mode == "H":
        # G2：求水头 H —— 以流量公式对 H 二分反解。
        #   原著§五：「计算Hd用逐次迭代法，│△H│＜0.001」→ 收敛判据用 tol_H。
        #   注意 Ho = H + Vo²/2g 亦随 H 变（Vo 由 Q 与 B0·(P+H) 定），故每轮重算。
        B = params["B"]
        Qtar = params["Q"]
        b = B / n if n else B
        bb0 = (B / B0) if B0 and B0 > 0 else 0.0

        def _q_of(Hc):
            _, hv = velocity_head(Qtar, B0, P, Hc, g)
            Ho_c = Hc + hv
            mm = _wide_m(mi, P, Hc, MJ, bb0, ctg, rb)
            sg = sigma_kdr(hs, Ho_c)
            ep = eps_of(n, b, Ho_c, xk, xo_of(SP, hs, Ho_c, n), include_xk)
            return sg * ep * mm * B * math.sqrt(2 * g) * Ho_c ** 1.5

        lo, hi = 1e-9, 10.0
        while _q_of(hi) < Qtar and hi < 1e6:
            hi *= 2.0
        for _ in range(max_iter):
            Hc = 0.5 * (lo + hi)
            if _q_of(Hc) > Qtar:
                hi = Hc
            else:
                lo = Hc
            if hi - lo < tol_H:
                break
        H = 0.5 * (lo + hi)
        v0, hv = velocity_head(Qtar, B0, P, H, g)
        H0 = H + hv
        m = _wide_m(mi, P, H, MJ, bb0, ctg, rb)
        sig = sigma_kdr(hs, H0)
        xo = xo_of(SP, hs, H0, n)
        ep = eps_of(n, b, H0, xk, xo, include_xk)
        Q = sig * ep * m * B * math.sqrt(2 * g) * H0 ** 1.5
        return dict(TY=3, MI=mi, mode=mode, Vo=v0, H0=H0, hsH0=(hs / H0 if H0 else 0),
                    SG=sig, XK=xk, X0=xo, M=m, EP=ep, B=B, Q=Q, b=b, n=n, H=H)

    # mode == "B"：由流量公式直接反解 B（原著五(五)）
    #   B·ε = B − 0.2·(ξk+(n−1)ξo)·Ho（因 ε=1−0.2·C·Ho/(n·b) 且 n·b=B）
    #   ⇒ B = Q/(σ·m·√(2g)·Ho^1.5) + 0.2·C·Ho；Ho 含 Vo²/2g，故整体再迭代一轮。
    Qtar = params["Q"]
    v0, hv = velocity_head(Qtar, B0, P, H, g)
    H0 = H + hv
    sig = sigma_kdr(hs, H0)
    xo = xo_of(SP, hs, H0, n)
    if mi == 3:
        # 无坎：m 依赖 B/B0，与 B 互为隐式，内层二分
        B = _solve_B_wukan(Qtar, sig, n, H0, xk, xo, P, H, hs,
                           B0, MJ, ctg, rb, g)
    else:
        m = _wide_m(mi, P, H, MJ, 0.0, ctg, rb)
        B = Qtar / (sig * m * math.sqrt(2 * g) * H0 ** 1.5) \
            + 0.2 * (xk + (n - 1) * xo) * H0
    b = B / n if n else B
    ep = eps_of(n, b, H0, xk, xo, include_xk)
    bb0 = (B / B0) if B0 and B0 > 0 else 0.0
    m = _wide_m(mi, P, H, MJ, bb0, ctg, rb)
    Q = sig * ep * m * B * math.sqrt(2 * g) * H0 ** 1.5
    return dict(TY=3, MI=mi, mode=mode, Vo=v0, H0=H0, hsH0=(hs / H0 if H0 else 0),
                SG=sig, XK=xk, X0=xo, M=m, EP=ep, B=B, Q=Q, b=b, n=n)


def _solve_B_wukan(Qtar, sig, n, H0, xk, xo, P, H, hs, B0, MJ, ctg, rb, g):
    """无坎宽顶堰求 B：m = f(B/Bo) 与 B 互为隐式，二分求解。"""
    def f(B):
        bb0 = (B / B0) if B0 and B0 > 0 else 0.0
        m = m_wukan(MJ, bb0, ctg, rb)
        ep = eps_of(n, B / n, H0, xk, xo, include_xk=False)
        return sig * ep * m * B * math.sqrt(2 * g) * H0 ** 1.5 - Qtar

    lo, hi = 1e-6, max(B0 * 1.0, 1.0)
    while f(hi) < 0 and hi < 1e6:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def calc_thin(params, g=G, tol=1e-9, max_iter=200):
    """薄壁堰（TY=1）：求 Q。矩形堰 TZ=1 / 三角堰 TZ=2。"""
    tz = int(params.get("TZ", 1) or 1)
    H = params["H"]
    P = params.get("P", 0.0)
    b = params.get("b", params.get("B", 0.0))
    B0 = params.get("B0", 0.0)
    if tz == 2:
        theta = math.radians(params.get("ANGLE", 90.0))
        mu = TZ_VNOTCH_MU
        Q = (8.0 / 15.0) * mu * math.sqrt(2 * g) * math.tan(theta / 2.0) * H ** 2.5
        return dict(TY=1, TZ=tz, mode="Q", Vo=0.0, H0=H, hsH0=0.0, SG=1.0,
                    XK=0.0, X0=0.0, M=mu, EP=1.0, B=b, Q=Q, b=b, n=1)
    # 矩形堰
    contract = bool(B0 and b and B0 > 0)
    Q = 0.0
    v0, H0, m = 0.0, H, 0.0
    for _ in range(max_iter):
        if contract:
            m = (0.405 + 0.0027 / H - 0.03 * (1.0 - b / B0)) * \
                (1.0 + 0.55 * (H / (H + P)) ** 2 * (b / B0) ** 2)
        else:
            m = (0.405 + 0.0027 / H) * (1.0 + 0.55 * (H / (H + P)) ** 2)
        xx = B0 if (B0 and B0 > 0) else b
        v0, hv = velocity_head(Q, xx, P, H, g)
        H0 = H + hv
        Qn = m * b * math.sqrt(2 * g) * H0 ** 1.5
        if abs(Qn - Q) < tol:
            Q = Qn
            break
        Q = Qn
    return dict(TY=1, TZ=tz, mode="Q", Vo=v0, H0=H0, hsH0=0.0, SG=1.0,
                XK=0.0, X0=0.0, M=m, EP=1.0, B=b, Q=Q, b=b, n=1)


def calc_wes(params, mode="Q", g=G, tol_H=0.001, max_iter=300):
    """实用堰（TY=2，WES）：求 Q / 求 Hd（迭代 |ΔH|<0.001）/ 求 B。"""
    P = params.get("P", 0.0)
    hs = params.get("HS", 0.0)
    B0 = params.get("B0", 0.0)
    n = int(params.get("N", 1) or 1)
    X1 = params.get("X1", 0)
    SP = params.get("SP", 1)
    xk = xk_of(X1)
    Hd = params.get("Hd", 0.0)

    def core(H, Hd):
        v0, hv = velocity_head(params.get("Q", 0.0), B0, P, H, g)
        Ho = H + hv
        ratio = (Ho / Hd) if Hd > 0 else 1.0
        m = interp(WES_HOHD, WES_M, ratio)
        sig = sigma_wes(hs, Ho) if hs > 0 else 1.0
        b = (params.get("B", 0.0) / n) if n else params.get("B", 0.0)
        xo = xo_of(SP, hs, Ho, n)
        ep = eps_of(n, b, Ho, xk, xo)
        return v0, Ho, m, sig, ep, xo, b

    if mode == "Q":
        B = params["B"]
        H = params["H"]
        v0, Ho, m, sig, ep, xo, b = core(H, Hd if Hd > 0 else H)
        Q, v0, Ho = 0.0, v0, Ho
        for _ in range(max_iter):
            v0, Ho, m, sig, ep, xo, b = core(H, Hd if Hd > 0 else H)
            Qn = sig * ep * m * B * math.sqrt(2 * g) * Ho ** 1.5
            if abs(Qn - Q) < 1e-9:
                Q = Qn
                v0, Ho, m, sig, ep, xo, b = core(H, Hd if Hd > 0 else H)
                break
            Q = Qn
        return dict(TY=2, mode=mode, Vo=v0, H0=Ho, Hd=Hd,
                    hsH0=(hs / Ho if (hs > 0 and Ho) else 0.0), SG=sig, XK=xk,
                    X0=xo, M=m, EP=ep, B=B, Q=Q, b=b, n=n,
                    HoHd=(Ho / Hd if Hd > 0 else 0.0))

    if mode == "H":
        # 求设计水头 Hd：迭代使 Q = Qtar（|ΔH|<0.001）
        B = params["B"]
        Qtar = params["Q"]
        Hd = params.get("Hd", 0.0) or 1.0
        for _ in range(max_iter):
            H = Hd
            v0, Ho, m, sig, ep, xo, b = core(H, Hd)
            Q = sig * ep * m * B * math.sqrt(2 * g) * Ho ** 1.5
            # m 随 Ho/Hd 变；对 Hd 做牛顿-割线
            dH = max(0.01, abs(Hd) * 1e-3)
            H2 = Hd + dH
            v0b, Hob, mb, sigb, epb, xob, bb = core(H2, H2)
            Q2 = sigb * epb * mb * B * math.sqrt(2 * g) * Hob ** 1.5
            if abs(Q - Qtar) < 1e-6:
                break
            slope = (Q2 - Q) / dH
            if abs(slope) < 1e-12:
                break
            Hd_new = Hd - (Q - Qtar) / slope
            if abs(Hd_new - Hd) < tol_H:
                Hd = Hd_new
                break
            Hd = Hd_new
        H = Hd
        v0, Ho, m, sig, ep, xo, b = core(H, Hd)
        Q = sig * ep * m * B * math.sqrt(2 * g) * Ho ** 1.5
        return dict(TY=2, mode=mode, Vo=v0, H0=Ho, Hd=Hd,
                    hsH0=(hs / Ho if (hs > 0 and Ho) else 0.0), SG=sig, XK=xk,
                    X0=xo, M=m, EP=ep, B=B, Q=Q, b=b, n=n,
                    HoHd=(Ho / Hd if Hd > 0 else 0.0))

    # mode == "B"
    #   B·ε = B − 0.2·(ξk+(n−1)ξo)·Ho（ε=1−0.2·C·Ho/(n·b)，n·b=B）
    Qtar = params["Q"]
    H = params["H"]
    Hd = Hd if Hd > 0 else H
    v0, Ho, m, sig, ep, xo, _ = core(H, Hd)
    B = Qtar / (sig * m * math.sqrt(2 * g) * Ho ** 1.5) \
        + 0.2 * (xk + (n - 1) * xo_of(SP, hs, Ho, n)) * Ho
    b = B / n if n else B
    ep = eps_of(n, b, Ho, xk, xo)
    Q = sig * ep * m * B * math.sqrt(2 * g) * Ho ** 1.5
    return dict(TY=2, mode=mode, Vo=v0, H0=Ho, Hd=Hd,
                hsH0=(hs / Ho if (hs > 0 and Ho) else 0.0), SG=sig, XK=xk,
                X0=xo, M=m, EP=ep, B=B, Q=Q, b=b, n=n,
                HoHd=(Ho / Hd if Hd > 0 else 0.0))


# ── 解析 ─────────────────────────────────────────────────────

INT_ORDER = ["TY", "MI", "MODE_CODE", "B0", "N", "B", "P", "H", "Hd", "HS",
             "X1", "SP", "Q", "TZ", "CTG", "RB", "ANGLE"]

# G1：INT 第 3 字段是原著§五字符变量 A$ 的数值编码「计算内容」（非翼墙形式号 MJ）
MODE_CODE_MAP = {1: "Q", 2: "H", 3: "B"}


def _parse_compact_wide(nums):
    """
    紧凑格式（逗号单行）解析 —— TY=3 无坎专属，D-4G-10~13 即此格式。

    字段序随「计算内容 MODE」与「闸孔数 n」而变（G1/G3 定案，依据见模块 docstring）：
      求Q(mode=1): [TY, MI, 1, B0, n, B, H, Hd, hs, X1, (SP|n>1), r/b, B/Bo]
      求H(mode=2): [TY, MI, 2, B0, n, B, Q,     hs, X1, (SP|n>1), r/b, B/Bo]
      求B(mode=3): [TY, MI, 3, B0, n, Q, H, Hd, hs, X1, (SP|n>1), r/b, B/Bo]
    说明书 §三(四)：无坎宽顶堰 P≡0；翼墙形式由 r/b 走 MJ=3（圆角形）分支。
    四例代入后字段数 11/12/12/13，与 38/43/46/45 字节的逗号文本逐例吻合（例13 因 n=2 多一 SP）。
    """
    TY = int(round(nums[0]))
    MI = int(round(nums[1]))
    code = int(round(nums[2]))
    p = {"程序": PROGRAM_ID, "TY": TY, "MI": MI, "MODE_CODE": code,
         "MODE": MODE_CODE_MAP.get(code, "Q"), "P": 0.0, "MJ": 3}
    p["B0"] = nums[3]
    p["N"] = int(round(nums[4]))
    i = 5
    if p["MODE"] == "H":                      # 求 H：H/Hd 未知 → 槽位缺省
        p["B"], p["Q"], p["HS"] = nums[i], nums[i + 1], nums[i + 2]
        i += 3
    else:                                      # 求 Q / 求 B：[B|Q], H, Hd, hs
        v5 = nums[i]
        p["H"], p["Hd"], p["HS"] = nums[i + 1], nums[i + 2], nums[i + 3]
        i += 4
        if p["MODE"] == "B":
            p["Q"], p["B"] = v5, None
        else:
            p["B"] = v5
    p["X1"] = int(round(nums[i])) if i < len(nums) else 0
    i += 1
    if p["N"] > 1 and i < len(nums):           # 多孔才有闸墩形状系数槽
        p["SP"] = int(round(nums[i]))
        i += 1
    p.setdefault("SP", 1)
    if i < len(nums):
        p["RB"] = nums[i]
        i += 1
    if i < len(nums):
        p["BB0"] = nums[i]
    return p


def _parse_legacy(nums):
    """旧式（CRLF 逐行）字段序；第 3 槽为计算内容模式码（G1），U1 不确定性见模块 docstring。"""
    p = {k: v for k, v in zip(INT_ORDER, nums)}
    p["程序"] = PROGRAM_ID
    p["TY"] = int(round(p["TY"]))
    code = p.pop("MODE_CODE", None)
    code = int(round(code)) if code is not None else None
    p["MODE_CODE"] = code
    for k in ("MI", "N", "X1", "SP", "TZ"):
        if k in p:
            p[k] = int(round(p[k]))
    p.setdefault("MJ", 1)                      # 翼墙形式号在 INT 中已被模式码占用 → 取默认
    mode = p.get("MODE")
    if mode is None:
        if code in MODE_CODE_MAP:              # G1：第 3 字段即原著§五的 A$
            mode = MODE_CODE_MAP[code]
        elif p.get("Q", 0) <= 0 and p.get("B", 0) > 0:
            mode = "Q"
        elif p.get("Hd", 0) > 0:
            mode = "H"
        else:
            mode = "B"
    p["MODE"] = mode
    return p


def parse(data):
    """
    解析输入。data: dict | INT 文件路径。

    dict（推荐，规避 U1 的字段顺序不确定性）键：
      TY 堰型 1/2/3；MI 宽顶堰型号 1 直坎/2 圆坎/3 无坎；MJ 翼墙 1 直角/2 八字/3 圆角(圆弧)；
      CTG 八字翼墙 ctgθ（0/0.5/1/2）；RB 圆角翼墙 r/b（网格 0/0.2/0.3/0.5，中间线性内插）；
      B0 上游渠宽；N 闸孔数；B 堰宽；P 坎高；H 堰上水头；Hd 设计水头；HS 淹没水深；
      X1 边墩形状 1~4；SP 闸墩形状 1~5；Q 流量；TZ 薄壁类型 1 矩形/2 三角；ANGLE 三角堰内角；
      MODE（可选）"Q"/"H"/"B"。

    INT 文件（自动分派）：
      · 紧凑格式（逗号单行、TY=MI=3）→ _parse_compact_wide()，字段序已唯一锁定（G1/G3）
      · 旧式（CRLF 逐行）→ _parse_legacy()，第 3 槽=计算内容模式码，
        其余槽位仍存 U1 不确定性（仅「求 Q、MI=1/2」一支被 D-4G-1.OUT 完整锁定）
    """
    if isinstance(data, dict):
        p = dict(data)
        if "TY" not in p:
            raise ValueError("dict 输入缺少 TY（堰型 1/2/3）")
        p["TY"] = int(p["TY"])
        p.setdefault("MODE", None)
        return p

    raw = open(data, "rb").read()
    nums = read_numbers(data)
    if len(nums) < 2:
        raise ValueError("INT 数据至少需要 2 个数值（TY, MI）")
    compact = (b"," in raw) and int(round(nums[0])) == 3 and int(round(nums[1])) == 3
    return _parse_compact_wide(nums) if compact else _parse_legacy(nums)


# ── 计算 ─────────────────────────────────────────────────────

def compute(params):
    """执行堰流计算，返回结构化结果。"""
    ty = int(params["TY"])
    mode = params.get("MODE") or "Q"
    if ty == 1:
        s = calc_thin(params)
    elif ty == 2:
        s = calc_wes(params, mode)
    elif ty == 3:
        s = calc_wide(params, mode)
    else:
        raise ValueError(f"未知堰型号 TY={ty}（应为 1/2/3）")
    result = {
        "程序": PROGRAM_ID,
        "堰型": TY_NAMES[ty],
        "功能": mode,
        "基本资料": {k: params[k] for k in
                     ("MI", "MJ", "B0", "N", "P", "H", "Hd", "HS", "X1", "SP")
                     if k in params},
        "计算": s,
        "结果": {
            "Vo": round(s["Vo"], 4),
            "H0": round(s["H0"], 4),
            "SG": round(s["SG"], 4),
            "XK": round(s["XK"], 4),
            "X0": round(s["X0"], 4),
            "M": round(s["M"], 4),
            "EP": round(s["EP"], 4),
            "B": round(s["B"], 4),
            "Q": round(s["Q"], 4),
        },
    }
    return result


# ── 输出 ─────────────────────────────────────────────────────

def _f4(x):
    return f"{x:.4f}"


def render(params, result):
    """生成汉字计算书（复刻原著 D-4G-1.OUT 风格）。"""
    ty = int(params["TY"])
    s = result["计算"]
    L = []
    L.append("  一、基本资料与计算假定")
    L.append("")
    L.append(f"      {TY_NAMES[ty]}")
    if ty == 1:
        L.append(f"  薄壁堰类型号TZ:            "
                 f"{TZ_NAMES.get(int(params.get('TZ', 1)), ''):>10}")
        L.append(f"  堰宽B:                     {_f4(params.get('b', params.get('B', 0))):>10}  (米)")
        L.append(f"  坎高P:                     {_f4(params.get('P', 0)):>10}  (米)")
        L.append(f"  堰上水头H:                 {_f4(params.get('H', 0)):>10}  (米)")
    else:
        L.append(f"  闸孔数n:                   {int(params.get('N', 1)):>10}")
        L.append(f"  堰宽B:                     {_f4(params.get('B', 0)):>10}  (米)")
        L.append(f"  坎高P:                     {_f4(params.get('P', 0)):>10}  (米)")
        L.append(f"  堰上水头H:                 {_f4(s.get('H', params.get('H', 0) or 0)):>10}  (米)")
        L.append(f"  堰上设计水头Hd:            {_f4(params.get('Hd', 0)):>10}  (米)")
        L.append(f"  淹没水深hs:                {_f4(params.get('HS', 0)):>10}  (米)")
    L.append("")
    L.append(f"  三、公式与算法({TY_NAMES[ty]})")
    L.append("")
    L.append("  (１) 流量公式")
    L.append("      Q=σ*ε*m*B*SQR(2*g)*Ho^(3/2)")
    L.append("  (２) 侧收缩系数  ε=1-0.2*(ξk+(n-1)*ξo)*Ho/(n*b)")
    L.append("  (３) 流量系数 m 因前沿形式而异（详见内核 docstring）")
    L.append("")
    L.append("  四、计算结果")
    L.append("")
    L.append(f"  行进流速 Vo:               {s['Vo']:>10.2f}  (米/秒)")
    L.append(f"  H0=H+V0*V0/2/G：           {s['H0']:>10.2f}  (米)")
    L.append(f"  淹没高度/水头 hs/Ho:       {_f4(s['hsH0']):>10}")
    L.append("  === 系数 ===")
    L.append(f"  淹没系数σ:                 {_f4(s['SG']):>10}")
    L.append(f"  边墩形状系数 ξk:           {_f4(s['XK']):>10}")
    L.append(f"  闸墩形状系数 ξo:           {_f4(s['X0']):>10}")
    L.append(f"  流量系数m:                 {_f4(s['M']):>10}")
    L.append("  === 校核 ===")
    L.append(f"  堰宽B:                     {s['B']:>10.2f}  (米)")
    L.append(f"  侧收缩系数ε:               {_f4(s['EP']):>10}")
    L.append(f"  流量Q:                     {s['Q']:>10.2f}  (立方米/秒)")
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
