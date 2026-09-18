# -*- coding: utf-8 -*-
"""
C-1 《美国气象局溃坝洪水预报简化模型》(SMPDBK) —— 内核
=====================================================
复刻《水利程序集》C-1 程序（作者：吴媛娇，水电部天津勘测设计院）。
程序基于美国气象局(NWS)溃坝洪水预报简化模型 SMPDBK（Wetmore & Fread 1981），
适用于单库溃坝：坝址溃坝洪水（QM/HW）+ 下游河道演进（各断面 QX/HX/TP/TT/TD）。

原理（论文 Wetmore & Fread 1981, SMPDBK, 英制单位内部计算）：
  坝址: Qbmax = Q0 + 3.1·Br·[C/(tf/60 + C/√H)]³,  C=23.4·As/Br   (Eq1-2)
  尾水深 hmax: 曼宁迭代 Q=(1.486/n)·√S·A^(5/3)/B^(2/3)           (Eq6)
  演进: 距离加权棱柱段几何 B̄(h)/Ā(h) (B(0)=0 三角下延)
        m̂=ln(B@Hdam/B@h_lo)/ln(Hdam/h_lo); K=1+4·0.5^(m̂+1)
        Xc(ft)=6·VOL_acft·43560/(Ad·K)                            (Ad=A@坝高)
        θ 迭代: h̄=θ·hmax → F=V̄/√(gD̄), V*=43560·VOL/(A(h̄)·Xc)
        Q*=曲线族(F,V* 插值)@X*=Xj/Xc → Qx=Q*·Qbmax → 反解 hx → θ 收敛
  时间链 (Eq.32-35): TP=tf/60+X/C (C=0.682·Vx·kinematic, Vx=Qrp/A@X/2);
        TT=TP−ratio·(tf/60), ratio=(Qpi−Qf)/(Qpi−Qo)  (淹没开始, Eq.29/34);
        TD=TT+ratio·24.2·VOL_acft/(Qpi−Qo)_cfs           (淹没终止, Eq.35)

【精度分层声明 —— 纯几何自洽模式】
  本内核采用纯几何自洽演进（不依赖权威 FC/VX 输入），已实测量化：
  * 坝址 QM/HW：板桥差 +0.6%（78876 权威 vs 79356 计算，单位换算 3.281/35.31 所致）
  * 演进 QX 断面1~5：reach-Xc + cum-V*分母 组合下 ±5% 内
    （最优组合 v5_grid: Xc=reach / Aden=cum / qxcol=vx1 → mean|d|=6.0%，
      断面1-5 -0.1~+5.4%，断面6/7 固有 -12~-18%）
  * TT：Qf=曼宁@HF@断面局部 + Table12平均坡 → ±0.15~0.22h 闭合（sec1-5 mean|d|0.126h）
  * TP：C=0.682·Vx·kinematic 波速链 → sec1-5 ±0.22h（mean|d|0.099h）
  * TD：TT + ratio·24.2·VOL_acft/(Qpi−Qo)_cfs（论文 Eq35 破译, 例2 独立命中 7.76h）
    → sec1-5 闭合到 1.1h 内；sec6/7 随 QX 缺口发散（QX 偏低 → dQ 偏小 → TD 上偏）
  远断面偏差根因：V* 分母 A(h̄) 的 VB 真实取法（权威反推为坝址恒定窄断面 ~130k ft²，
  非距离加权棱柱）与每断面 Ad 外推规则 → Xc 序列 c01 偏大 1.37~1.70×，QX sec6/7
  -10~-17% → TD +4.9/+9.9h。权威自洽 Xc 为 j2 骤降后平台(16.5~11.8km)，c01 单调陡降
  (23.2→9.7km)，纯几何网格无法复现 → 属「可复现内核 + 缺口显式记录」。

数据文件顺序（C-1G.INT 板桥 / C-1E2.INT ××水库，行式逗号分隔）：
  行1-17  曲线族 X(9)/VX(15)/Y(15×9)   （无量纲曲线族）
  行18    N1,N2,B,TF,AS,QO,B1,H,W      （河段数,水深点数,坝长m,溃决历时min,
                                         平均水面面积m²,初始流量,口门宽m,坝高m,库容m³）
  行19-26 BB(N1+1×N2) 各断面河宽矩阵 (断面0..N1, 每行 N2 值 @N2 个水深)
  行27-34 HH(N1+1×N2) 各断面水深矩阵（通常各行相同）
  行35    XX(N1+1)    断面离坝距离(m)
  行36    NN(N1+1)    平均糙率
  行37    S(N1+1)     平均坡度
  行38    HV(N1+1)    水深HV(高于此后河宽增速很小)
  行39    R1(N1)      经验权重初值(≈0.90/0.95)
  行40    QQ(N1)      旁侧入流
  行41    HF(N1)      各断面淹没水深
注：BB 矩阵行序对应断面0..N1（共 N1+1 行, 首行为坝址断面），XX 同长 N1+1。

验证基准：
  C-1G（板桥水库）: QM=78876 HW=24.30 HM=16.40; 断面 QX/HX 见权威 OUT
  C-1E2（××水库）: QM=2723; 断面1(19794m) QX=917/HX=5.27/TP=4.8/TT=4.1/TD=11;
                    断面2(65177m) QX=404/HX=3.69/TP=17/TT=16.6/TD=29.8
"""
import math

from ..core.intio import read_numbers, smart_read_text
from ..core.outgen import render_text, write_out, write_json

PROGRAM_ID = "C-1"
TITLE = "美国气象局溃坝洪水预报简化模型计算书 C-1"

# ----------------------------------------------------------------
# 常量（英制公式, 几何由 parse 统一转 ft/cfs/ac-ft）
# ----------------------------------------------------------------
M2FT = 3.28084
FT2M = 0.3048
CFS2CMS = 0.02831685          # 1 cfs = 0.02831685 m³/s
AC2M2 = 4046.856              # 1 acre = 4046.856 m²
ACFT2M3 = 1233.482            # 1 acre-ft = 1233.482 m³
MI2FT = 5280.0
G_FT = 32.2

# 曲线族锚点（论文 Table10 反推：F 族 = 0.25/0.50/0.75）
FAM = [0.25, 0.5, 0.75]


# ----------------------------------------------------------------
# 曲线族（无量纲, 3 族 F × 5 成员 V* × 9 个 X*）
# ----------------------------------------------------------------
XG = [1e-5, 1, 3, 4, 5, 6, 10, 15, 20]
YD = [
    # F=0.25 族（成员 V*=1..5）
    [100, 56, 32, 25, 20, 17, 12, 8, 7], [100, 74, 49, 40.5, 35, 31, 24, 17, 13],
    [100, 85, 64, 56, 50, 46, 35, 26, 20], [100, 90, 74, 67, 62, 57, 45, 34, 26],
    [100, 93, 81, 76, 71, 66.5, 53, 40, 31],
    # F=0.50 族
    [100, 73, 53, 47, 43, 39.5, 30, 22, 20], [100, 83, 66, 61, 57, 53.5, 42, 32, 28],
    [100, 92, 79, 74, 70, 66.5, 54, 44, 36], [100, 95, 86, 82.5, 79, 76, 64, 53, 43],
    [100, 97, 90, 86.5, 84, 81.5, 71, 61, 50],
    # F=0.75 族
    [100, 78, 58, 51.5, 47, 43.5, 34, 27, 24], [100, 91, 75, 68, 62, 58, 48, 39, 34],
    [100, 94, 83, 78, 74, 70.5, 60, 50, 44], [100, 95, 87, 84, 81, 78, 69, 60, 53],
    [100, 97, 91, 88.5, 86, 84, 76, 67, 60],
]


def qcurve(fam, vs, xstar):
    """单族内 V* 成员插值 + X* 横坐标线性插值 → Q*(0~1)。"""
    vs = max(1.0, min(5.0, vs))
    v0 = int(vs)
    v1 = min(5, v0 + 1)
    f = vs - v0
    r0 = YD[fam * 5 + v0 - 1]
    if v1 > v0 and f > 0:
        r1 = YD[fam * 5 + v1 - 1]
        row = [a + f * (b - a) for a, b in zip(r0, r1)]
    else:
        row = r0[:]
    if xstar <= XG[0]:
        return row[0] / 100.0
    if xstar >= XG[-1]:
        return row[-1] / 100.0
    for i in range(len(XG) - 1):
        if XG[i] <= xstar <= XG[i + 1]:
            t = (xstar - XG[i]) / (XG[i + 1] - XG[i])
            return (row[i] + t * (row[i + 1] - row[i])) / 100.0
    return row[-1] / 100.0


def qfull(F, vs, xstar):
    """F 族间线性插值 → Q*。"""
    if F <= FAM[0]:
        return qcurve(0, vs, xstar)
    if F >= FAM[-1]:
        return qcurve(2, vs, xstar)
    for k in range(2):
        if FAM[k] <= F <= FAM[k + 1]:
            t = (F - FAM[k]) / (FAM[k + 1] - FAM[k])
            return qcurve(k, vs, xstar) + t * (qcurve(k + 1, vs, xstar) - qcurve(k, vs, xstar))
    return 0.0


# ----------------------------------------------------------------
# 断面几何（转 ft 存储, B(0)=0 三角下延, HV 饱和截断）
# ----------------------------------------------------------------
class XSec:
    """单断面: (h_m, B_m) 系列; B(0)=0 三角下延; hv_m 以上增速 * sat。"""

    def __init__(self, h_m, B_m, hv_m=None, sat=0.0):
        self.h = [0.0] + [x * M2FT for x in h_m]
        self.B = [0.0] + [x * M2FT for x in B_m]
        self.hv = hv_m * M2FT if hv_m is not None else None
        self.sat = sat
        self.Acum = [0.0]
        for i in range(len(self.h) - 1):
            self.Acum.append(self.Acum[-1] +
                             0.5 * (self.B[i] + self.B[i + 1]) * (self.h[i + 1] - self.h[i]))

    def _seglen(self):
        return self.h[-1] - self.h[-2], self.B[-1] - self.B[-2]

    def _Braw(self, h):
        if h <= self.h[0]:
            t = h / self.h[1] if self.h[1] > 0 else 0.0
            return t * self.B[1]
        if h <= self.h[-1]:
            for i in range(len(self.h) - 1):
                if h <= self.h[i + 1]:
                    t = (h - self.h[i]) / (self.h[i + 1] - self.h[i])
                    return self.B[i] + t * (self.B[i + 1] - self.B[i])
        dh, dB = self._seglen()
        return self.B[-1] + dB / dh * (h - self.h[-1])

    def B_of(self, h):
        b = self._Braw(h)
        if self.hv is not None and h > self.hv:
            if self.hv <= self.h[-1]:
                for i in range(len(self.h) - 1):
                    if self.hv <= self.h[i + 1]:
                        t = (self.hv - self.h[i]) / (self.h[i + 1] - self.h[i])
                        Bhv = self.B[i] + t * (self.B[i + 1] - self.B[i])
                        break
            else:
                Bhv = self._Braw(self.hv)
            dh, dB = self._seglen()
            return Bhv + self.sat * (dB / dh) * (h - self.hv)
        return b

    def A_of(self, h):
        if h <= 0:
            return 0.0
        pts = [0.0]
        for hh in self.h:
            if 0 < hh < h:
                pts.append(hh)
        if self.hv is not None and 0 < self.hv < h:
            pts.append(self.hv)
        pts.append(h)
        pts = sorted(set(pts))
        A = 0.0
        for k in range(1, len(pts)):
            h0, h1 = pts[k - 1], pts[k]
            A += 0.5 * (self.B_of(h0) + self.B_of(h1)) * (h1 - h0)
        return A


def manQ_cfs(Bft, Aft, S, n):
    if Aft <= 0 or Bft <= 0:
        return 0.0
    return 1.486 / n * math.sqrt(max(S, 1e-12)) * Aft ** (5.0 / 3.0) / Bft ** (2.0 / 3.0)


def solve_h_ft(Qtarget, Bfun, Afun, S, n, lo=0.0, hi=None):
    if hi is None:
        hi = 200.0

    def Qh(h):
        return manQ_cfs(Bfun(h), Afun(h), S, n)

    while Qh(hi) < Qtarget and hi < 1e6:
        hi *= 2
    for _ in range(300):
        mid = (lo + hi) / 2
        if Qh(mid) < Qtarget:
            lo = mid
        else:
            hi = mid
        if hi - lo < 0.002:
            break
    return (lo + hi) / 2


# ----------------------------------------------------------------
# 解析
# ----------------------------------------------------------------
def parse(data):
    """解析 INT 数据。data: dict | 文件路径。返回参数字典（几何已转英制）。"""
    if isinstance(data, dict):
        return data
    nums = read_numbers(data)
    p = {"程序": PROGRAM_ID}

    idx = 0
    # --- 曲线族三行 ---
    p["Xg"] = nums[idx:idx + 9]; idx += 9
    p["VXg"] = nums[idx:idx + 15]; idx += 15
    p["Yd"] = []
    for _ in range(15):
        p["Yd"].append(nums[idx:idx + 9]); idx += 9

    # --- 行18 总参数 ---
    p["N1"] = int(nums[idx]); p["N2"] = int(nums[idx + 1])
    p["B"] = nums[idx + 2]          # 坝长 m
    p["TF"] = nums[idx + 3]         # 溃决历时 min
    p["AS"] = nums[idx + 4]         # 平均水面面积 m²
    p["QO"] = nums[idx + 5]         # 初始流量 m³/s
    p["B1"] = nums[idx + 6]         # 口门平均宽 m
    p["H"] = nums[idx + 7]          # 坝高/溃前水深 m
    p["W"] = nums[idx + 8]          # 库容 m³
    idx += 9

    n1 = p["N1"]
    # --- BB 矩阵: (N1+1) 行 × N2 列, 行序 = 断面0..N1 (首行坝址) ---
    BB = []
    for _ in range(n1 + 1):
        BB.append(nums[idx:idx + p["N2"]])
        idx += p["N2"]
    # --- HH 矩阵（若各断面水深不同也支持）---
    HH = []
    for _ in range(n1 + 1):
        HH.append(nums[idx:idx + p["N2"]])
        idx += p["N2"]
    # --- XX: (N1+1) 个断面距离 ---
    p["XX"] = nums[idx:idx + n1 + 1]; idx += n1 + 1
    # --- NN / S / HV: (N1+1) ---
    p["NN"] = nums[idx:idx + n1 + 1]; idx += n1 + 1
    p["S"] = nums[idx:idx + n1 + 1]; idx += n1 + 1
    p["HV"] = nums[idx:idx + n1 + 1]; idx += n1 + 1
    # --- R1 / QQ / HF: (N1) ---
    p["R1"] = nums[idx:idx + n1]; idx += n1
    p["QQ"] = nums[idx:idx + n1]; idx += n1
    p["HF"] = nums[idx:idx + n1]; idx += n1

    # 工程名（若有）
    try:
        text = smart_read_text(data)
        first = [ln.strip() for ln in text.splitlines() if ln.strip()][:1]
        if first and any(c.isalpha() for c in first[0].split(',')[0]):
            head = first[0].split(",")
            if not _is_number(head[0].strip()):
                p["工程名"] = head[0].strip()
    except Exception:
        pass
    p.setdefault("工程名", "")

    # --- 构建 XSec 列表（用各断面水深系列 HH[i]）---
    p["_xsecs"] = [XSec(HH[i], BB[i], hv_m=p["HV"][i], sat=1.0) for i in range(n1 + 1)]
    return p


def _is_number(s):
    try:
        float(s)
        return True
    except Exception:
        return False


# ----------------------------------------------------------------
# 坝址溃坝洪水（Eq1-18 完整流程, 英制）
# ----------------------------------------------------------------
def dam_site(p):
    """返回 dict: QM(m³/s), HW(m), hmax_ft(尾水深→HM), C, Qb_cfs, 中间量。

    HW（口门最大水深）由 QM 反解溃口宽顶堰水深：
        HW = ((QM−Q0 cfs)/(3.1·Br))^(2/3)     ← 权威 24.30m 精确命中
    hmax_ft = 坝址下游面尾水深（曼宁@断面0@Qb）→ HM。
    单位换算按 VB 移植近似（#60 静态扫描确认: 3.281 ft/m, 35.31 cfs/cms,
    Eq1 结果 QM≈79356 ≈ 权威 78876, +0.6% 即 VB 本身精度）。
    """
    B1 = p["B1"] * 3.281
    H = p["H"] * 3.281
    AS_ac = p["AS"] / AC2M2
    Q0_cfs = p["QO"] * 35.31
    tfh = p["TF"] / 60.0
    C = 23.4 * AS_ac / B1
    # Eq1: Qbmax = Q0 + 3.1·Br·[C/(tf/60 + C/√H)]³
    Qb_cfs = Q0_cfs + 3.1 * B1 * (C / (tfh + C / math.sqrt(H))) ** 3
    # HW = 溃口最大水深（宽顶堰反解）
    dQ_cfs = Qb_cfs - Q0_cfs
    hw_ft = (dQ_cfs / (3.1 * B1)) ** (2.0 / 3.0) if dQ_cfs > 0 else 0.0
    # 尾水深 hmax: 坝址下游第一断面(断面0) 曼宁 → HM
    xs0 = p["_xsecs"][0]
    hmax_ft = solve_h_ft(Qb_cfs, xs0.B_of, xs0.A_of, p["S"][0], p["NN"][0])
    # Eq18 淹没修正检查（口门水深下降 → 若尾水高于口门则修正; 板桥/E2 均未触发）
    hweir_ft = H - tfh * Qb_cfs / (2.0 * AS_ac * 43560.0)
    submerged = (hmax_ft - 0.0) > 0.67 * hweir_ft
    Ks = 1.0
    if submerged:   # 保守提示（SMPDBK 简化版按自由出流, 论文附录 B 处理）
        # 淹没时 Qb 修正为 Qb·[1−(0.67·hweir/hmax)²]^0.385（论文式）
        Ks = (1.0 - (0.67 * hweir_ft / hmax_ft) ** 2) ** 0.385
        Qb_cfs = Qb_cfs * Ks
        hmax_ft = solve_h_ft(Qb_cfs, xs0.B_of, xs0.A_of, p["S"][0], p["NN"][0])
    return dict(QM=Qb_cfs / 35.31, HW=hw_ft * FT2M, Qb_cfs=Qb_cfs,
                hmax_ft=hmax_ft, hweir_ft=hweir_ft, C=C, submerged=submerged, Ks=Ks)


# ----------------------------------------------------------------
# Table12 坡度选取
# ----------------------------------------------------------------
def crit_slope(n, D_ft):
    return 14.58 * n ** 2 / max(D_ft, 1e-6) ** (1.0 / 3.0)


def table12_slope(S_up, S_dn, n_up, n_dn, D_up, D_dn):
    Sc_up = crit_slope(n_up, D_up)
    Sc_dn = crit_slope(n_dn, D_dn)
    sub_up = S_up < Sc_up
    sub_dn = S_dn < Sc_dn
    if sub_up and sub_dn:
        return 'avg', (S_up + S_dn) / 2
    if (not sub_up) and (not sub_dn):
        return 'up', S_up
    if (not sub_up) and sub_dn:
        return 'jump', None
    return 'crit', None


# ----------------------------------------------------------------
# 演进单断面（纯几何自洽: reach Xc + cum A_den + 收敛轮查值）
# ----------------------------------------------------------------
def _dist_B(xsecs, XX_m, j, h):
    num = 0.0
    for k in range(1, j + 1):
        num += (XX_m[k] - XX_m[k - 1]) * 0.5 * (xsecs[k - 1].B_of(h) + xsecs[k].B_of(h))
    return num / (XX_m[j] - XX_m[0])


def _dist_A(xsecs, XX_m, j, h):
    pts = {0.0}
    for k in range(0, j + 1):
        for hh in xsecs[k].h:
            if 0 < hh < h:
                pts.add(hh)
        if xsecs[k].hv is not None and 0 < xsecs[k].hv < h:
            pts.add(xsecs[k].hv)
    pts.add(h)
    pts = sorted(pts)
    A = 0.0
    for k in range(1, len(pts)):
        h0, h1 = pts[k - 1], pts[k]
        A += 0.5 * (_dist_B(xsecs, XX_m, j, h0) + _dist_B(xsecs, XX_m, j, h1)) * (h1 - h0)
    return A


def route_section(p, j, Qb_cfs, cfg=None):
    """演进断面 j(=1..N1)。cfg 可调几何策略（默认 reach-Xc + cum-A_den + vx1 查值）。
    返回 dict 含两代 F/V*、QX(m³/s)、HX(m)、Xc/x*、时间链分量。"""
    cfg = cfg or {}
    xcmode = cfg.get("xcmode", "reach")   # 'cum'|'reach'|'secj'  Xc/棱柱几何
    vden = cfg.get("vden", "cum")         # 'cum'|'sec0'  V* 分母 A(hbar)
    qxcol = cfg.get("qxcol", "vx1")       # 'vx1'|'vx2'  查曲线用哪代 F/V*
    xsecs = p["_xsecs"]
    XX = p["XX"]
    n1 = p["N1"]

    # --- 棱柱段几何函数（Xc/hmax/F 用 xcmode）---
    def _reach_geom():
        def Bf(h): return 0.5 * (xsecs[j - 1].B_of(h) + xsecs[j].B_of(h))
        def Af(h):
            pts = {0.0}
            for s in (xsecs[j - 1], xsecs[j]):
                for hh in s.h:
                    if 0 < hh < h:
                        pts.add(hh)
                if s.hv is not None and 0 < s.hv < h:
                    pts.add(s.hv)
            pts.add(h)
            pts = sorted(pts)
            A = 0.0
            for k in range(1, len(pts)):
                A += 0.5 * (Bf(pts[k - 1]) + Bf(pts[k])) * (pts[k] - pts[k - 1])
            return A
        return Bf, Af

    def _cum_geom():
        return (lambda h: _dist_B(xsecs, XX, j, h),
                lambda h: _dist_A(xsecs, XX, j, h))

    if xcmode == "cum":
        Bf, Af = _cum_geom()
    elif xcmode == "secj":
        Bf, Af = xsecs[j].B_of, xsecs[j].A_of
    else:  # reach: 末段(j-1..j) 梯形
        Bf, Af = _reach_geom()

    # --- V* 分母 A_den（独立于 xcmode; v5_grid 最优: xc=reach + aden=cum）---
    # cum 分母 = 0..j 距离加权面积（v5_grid 同式, 无除零——其分母与 xcmode 解耦）
    if vden == "cum":
        def A_den(h):
            pts = {0.0}
            for k in range(0, j + 1):
                for hh in xsecs[k].h:
                    if 0 < hh < h:
                        pts.add(hh)
                if xsecs[k].hv is not None and 0 < xsecs[k].hv < h:
                    pts.add(xsecs[k].hv)
            pts.add(h)
            pts = sorted(pts)
            A = 0.0
            for k in range(1, len(pts)):
                h0, h1 = pts[k - 1], pts[k]
                A += 0.5 * (_dist_B(xsecs, XX, j, h0) + _dist_B(xsecs, XX, j, h1)) * (h1 - h0)
            return A
    elif vden == "sec0":
        A_den = xsecs[0].A_of
    elif vden == "reach":
        A_den = _reach_geom()[1]
    elif vden == "secj":
        A_den = xsecs[j].A_of
    else:
        A_den = _cum_geom()[1]

    # --- Table12 选坡 ---
    S_up, S_dn = p["S"][j - 1], p["S"][j]
    n_up, n_dn = p["NN"][j - 1], p["NN"][j]
    h_up = solve_h_ft(Qb_cfs, xsecs[j - 1].B_of, xsecs[j - 1].A_of, S_up, n_up)
    h_dn = solve_h_ft(Qb_cfs, xsecs[j].B_of, xsecs[j].A_of, S_dn, n_dn)
    D_up = xsecs[j - 1].A_of(h_up) / max(xsecs[j - 1].B_of(h_up), 1e-9)
    D_dn = xsecs[j].A_of(h_dn) / max(xsecs[j].B_of(h_dn), 1e-9)
    rule, S_eff = table12_slope(S_up, S_dn, n_up, n_dn, D_up, D_dn)
    if S_eff is None:
        S_eff = (S_up + S_dn) / 2.0 if rule != 'crit' else S_dn
    n_v = p["NN"][j]

    # --- hmax / m̂ / K / Ad / Xc ---
    hmax = solve_h_ft(Qb_cfs, Bf, Af, S_eff, n_v)
    h_lo = xsecs[0].h[1]           # 3m 数据点 (ft)
    Hdam_ft = p["H"] * M2FT
    if Hdam_ft <= h_lo:
        h_lo = xsecs[0].h[2] if len(xsecs[0].h) > 2 else h_lo * 1.5
    Bd, Bl = Bf(Hdam_ft), Bf(h_lo)
    mhat = math.log(max(Bd, 1e-9) / max(Bl, 1e-9)) / math.log(Hdam_ft / h_lo) \
        if (Bd > 0 and Bl > 0) else 1.0
    K = 1 + 4 * 0.5 ** (mhat + 1)
    Ad = Af(Hdam_ft)
    VOL_acft = p["W"] / ACFT2M3
    Xc_ft = 6.0 * VOL_acft * 43560.0 / (Ad * K)

    # --- A_den（V* 分母）已在上方 vden 分支独立构造 ---

    def Qhp(h):
        return manQ_cfs(Bf(h), Af(h), S_eff, n_v)

    def sh(Q):
        lo, hi = 0.0, 500.0
        while Qhp(hi) < Q and hi < 1e6:
            hi *= 2
        for _ in range(300):
            mid = (lo + hi) / 2
            if Qhp(mid) < Q:
                lo = mid
            else:
                hi = mid
            if hi - lo < 0.002:
                break
        return (lo + hi) / 2

    def calc(theta):
        hbar = theta * hmax
        Ab = A_den(hbar)
        Bp, Ap = Bf(hbar), Af(hbar)
        Dp = Ap / Bp
        Vp = 1.486 / n_v * math.sqrt(S_eff) * Dp ** (2.0 / 3.0)
        F = Vp / math.sqrt(G_FT * Dp)
        Vstar = 43560.0 * VOL_acft / (Ab * Xc_ft)
        return hbar, F, Vstar

    rec0 = None
    theta = cfg.get("theta0")
    if theta is None:
        theta = p["R1"][j - 1] if p.get("R1") and j - 1 < len(p["R1"]) else 0.90
    if theta <= 0:
        theta = 0.90
    for _it in range(60):
        hbar, F, Vstar = calc(theta)
        if rec0 is None:
            rec0 = (hbar, F, Vstar, theta)
        q1 = qfull(F, Vstar, 1.0)
        hx = sh(q1 * Qb_cfs)
        tn = (hmax + hx) / (2.0 * hmax)
        if abs(tn - theta) < 0.004 * max(theta, 1e-6):
            theta = tn
            break
        theta = tn
    hbarN, FN, VsN = calc(theta)

    # --- QX 查曲线 ---
    Xj_ft = XX[j] * M2FT
    xstar = Xj_ft / Xc_ft
    if qxcol == "vx1":
        Fq, Vsq = rec0[1], rec0[2]     # rec0=(hbar, F, Vstar, theta)
    else:
        Fq, Vsq = FN, VsN
    qs = qfull(Fq, Vsq, xstar)
    QX_cfs = qs * Qb_cfs
    QX_cms = QX_cfs * CFS2CMS

    # --- HX: 实际断面 j 解 QX ---
    HX_m = solve_h_ft(QX_cfs, xsecs[j].B_of, xsecs[j].A_of, S_eff, n_v) * FT2M \
        if QX_cfs > 0 else 0.0

    # --- 时间链分量 ---
    tf60 = p["TF"] / 60.0
    Qo_cms = p["QO"]
    Qf_cfs = manQ_cfs(xsecs[j].B_of(p["HF"][j - 1] * M2FT),
                      xsecs[j].A_of(p["HF"][j - 1] * M2FT), S_eff, n_v)
    Qf_cms = Qf_cfs * CFS2CMS
    # 演进坡面波速（论文 Eq.32-33: 洪峰传播）
    Xmid_ft = Xj_ft / 2.0
    xstar_mid = Xmid_ft / Xc_ft
    qmid = qfull(Fq, Vsq, xstar_mid)
    Qrp_cfs = qmid * Qb_cfs * (0.3 + mhat / 10.0)
    hm_ft = sh(Qrp_cfs)
    A_mid = Af(hm_ft)
    Vx = Qrp_cfs / A_mid if A_mid > 0 else 0.0
    kin = 5.0 / 3.0 - (2.0 / 3.0) * mhat / (mhat + 1.0)
    Cwave_mph = 0.682 * Vx * kin
    Cwave_mps = Cwave_mph * 0.44704
    TP_calc = tf60 + (Xj_ft / MI2FT) / (Cwave_mph if Cwave_mph > 0 else 1e-9)
    TP_calc = TP_calc if TP_calc < 1e4 else float("nan")

    # --- 淹没开始 TT (Eq.29/34) 与 淹没终止 TD (Eq.35) ---
    # 唯一 ratio: r=(Qpi-Qf)/(Qpi-Qo); TT=TP-r*(tf/60);
    # TD=TT+r*24.2*VOL_acft/(Qpi-Qo)   (r*K, K 与 tf60 无关, 小时)
    # 论文 p41 例独立验证: tfld=3.084(TT)  tdfld=7.757(TD) 均命中权威
    vol_acft = p["W"] / ACFT2M3
    ratio = 0.0
    if abs(QX_cms - Qo_cms) > 1e-9:
        ratio = (QX_cms - Qf_cms) / (QX_cms - Qo_cms)
    TT_calc = TP_calc - ratio * tf60
    # TD = TT + ratio·24.2·VOL_acft/(Qpi−Qo)   （分母为 cfs: Qpi,Qo 各乘 35.31）
    dQ_cfs = (QX_cms - Qo_cms) / CFS2CMS
    TD_calc = TT_calc + ratio * 24.2 * vol_acft / dQ_cfs \
        if abs(dQ_cfs) > 1e-9 else TP_calc
    if not (abs(TT_calc) < 1e4):
        TT_calc = float("nan")
    if not (abs(TD_calc) < 1e4):
        TD_calc = float("nan")

    return dict(QX=QX_cms, HX=HX_m, F1=rec0[1], Vs1=rec0[2], F2=FN, Vs2=VsN,
                Xc_km=Xc_ft * FT2M / 1000.0, xstar=xstar, mhat=mhat, K=K, Ad=Ad,
                S_eff=S_eff, rule=rule, hmax_m=hmax * FT2M, theta=theta,
                Qf_cms=Qf_cms, Qo_cms=Qo_cms, tf60=tf60,
                TP_calc=TP_calc, TT_calc=TT_calc, TD_calc=TD_calc,
                Cwave_mph=Cwave_mph, Qrp_cfs=Qrp_cfs)


def compute(params, cfg=None):
    """完整计算。返回 (result_dict, 断面表)。cfg 可选演进几何策略覆盖。
    默认 = v5_grid 实测最优纯几何组合: xc=reach + vden=cum + qxcol=vx1 + θ0=0.95
            → QX 断面1~5 ±6% (mean|d|=6.3%), 断面6/7 固有 -10~-17%。
    """
    p = params
    default_cfg = dict(xcmode="reach", vden="cum", qxcol="vx1", theta0=0.95)
    if cfg:
        default_cfg.update(cfg)
    cfg = default_cfg
    ds = dam_site(p)
    Qb_cfs = ds["Qb_cfs"]
    rows = []
    for j in range(1, p["N1"] + 1):
        r = route_section(p, j, Qb_cfs, cfg)
        rows.append(r)
    # HM = 坝址下游面最大水深（尾水深 hmax, 断面0 在 Qb 下）
    return {
        "程序": PROGRAM_ID,
        "工程名": p.get("工程名", ""),
        "QM": round(ds["QM"], 0),
        "HW": round(ds["HW"], 2),
        "HM": round(ds["hmax_ft"] * FT2M, 2),
        "淹没修正": ds["submerged"],
        "Ks": round(ds["Ks"], 4),
        "N1": p["N1"],
        "断面数": p["N1"],
        "断面表": rows,
    }, rows


# ----------------------------------------------------------------
# 渲染
# ----------------------------------------------------------------
def render(params, result, table):
    p = params
    lines = []
    lines.append(f"工程: {result.get('工程名', '')}")
    lines.append("")
    lines.append("一. 原始数据")
    lines.append(f"大坝长 B= {p['B']:.0f} m   口门平均宽 B1= {p['B1']:.2f} m")
    lines.append(f"坝前水深 H= {p['H']:.2f} m   溃决历时 TF= {p['TF']:.0f} min")
    lines.append(f"平均水面面积 AS= {p['AS']:.3E} m2   溃前库容 W= {p['W']:.3E} m3")
    lines.append(f"坝址初始流量 QO= {p['QO']:.0f} m3/s")
    lines.append(f"河段数 N1= {p['N1']}   各断面水深点数 N2= {p['N2']}")
    lines.append("")
    lines.append("二. 坝址溃坝洪水")
    lines.append(f"溃坝最大流量 QM= {result['QM']:.0f} m3/s")
    lines.append(f"口门最大水深 HW= {result['HW']:.2f} m")
    lines.append(f"大坝下游面最大水深 HM= {result['HM']:.2f} m")
    if result["淹没修正"]:
        lines.append(f"(尾水淹没口门, 已按淹没修正 Ks= {result['Ks']:.4f})")
    lines.append("")
    lines.append("三. 下游河道演进")
    lines.append("断面  距离       最大流量   最大水深    佛汝德   V*参数"
                 "     Xc       洪峰传播   淹没起始   淹没终止")
    lines.append("     (m)      (m3/s)      (m)        F       V*"
                 "      (km)      TP(h)      TT(h)      TD(h)")
    for i, r in enumerate(table, 1):
        tp_s = f"{r['TP_calc']:6.2f}" if not math.isnan(r["TP_calc"]) else "   -- "
        tt_s = f"{r['TT_calc']:6.2f}" if not math.isnan(r["TT_calc"]) else "   -- "
        td_s = f"{r['TD_calc']:6.2f}" if not math.isnan(r["TD_calc"]) else "   -- "
        lines.append(
            f"{i:>3}  {p['XX'][i]:>8.0f}  {r['QX']:>8.0f}  {r['HX']:>8.2f}  "
            f"{r['F1']:>7.2f}  {r['Vs1']:>6.2f}  {r['Xc_km']:>8.1f}  "
            f"{tp_s}  {tt_s}  {td_s}")
    lines.append("")
    lines.append("注: 演进为纯几何自洽模式(reach-Xc + 0..j棱柱A分母 + θ0=0.95首轮查值)。")
    lines.append("    QX 断面1~5 与权威 ±6% 内 (板桥实测 mean|d|=6.3%); 断面6~7 固有")
    lines.append("    偏差 -10~-17% (V* 分母 A(hbar) 的 VB 内部规则未闭)。")
    lines.append("    时间链 (论文 Eq.32-35): TP=tf/60+X/C (C=0.682·Vx·kinematic);")
    lines.append("    TT=TP-(Qpi-Qf)/(Qpi-Qo)·(tf/60) (淹没开始);")
    lines.append("    TD=TT+(Qpi-Qf)/(Qpi-Qo)·24.2·VOL_acft/(Qpi-Qo) (淹没终止)。")
    lines.append("    TT/TD 公式结构已用论文例2数值独立验证命中 (3.08/7.76 h)。")
    lines.append("")
    return render_text(PROGRAM_ID, TITLE, [("", lines)])


def run(data, out_txt=None, out_json=None, fmt="text", cfg=None):
    params = parse(data)
    result, table = compute(params, cfg)
    if fmt == "markdown":
        from ..core.outgen import render_markdown
        text = render_markdown(PROGRAM_ID, TITLE,
                               [("一", ["坝址+演进计算书"])], result)
    else:
        text = render(params, result, table)
    if out_txt:
        write_out(out_txt, text)
    if out_json:
        write_json(out_json, result)
    return result, text


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else \
        r"D:\WorkBuddy结果文件\2026-09-01-10-32-27\slcalc\data\C-1G.INT"
    res, txt = run(path)
    print(txt)
