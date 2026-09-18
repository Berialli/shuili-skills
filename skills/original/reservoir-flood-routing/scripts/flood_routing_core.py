#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
水库双辅助曲线调洪演算核心计算引擎
=====================================
纯计算模块，无 GUI / matplotlib 依赖，可在任意环境（含无显示服务器）运行。
算法来源：原 Excel VBA 调洪演算程序（半图解法 / 双辅助曲线法）

核心能力：
  1. calc_auxiliary  : 计算调洪辅助曲线（V/dt±q/2 双辅助曲线）
  2. calc_routing    : 逐时段推求下泄流量与水库水位
  3. calc_flood      : 一站式接口（辅助曲线 + 演算 + 特征值提取）
  4. summary         : 提取调洪特征值（最大下泄流量/最高水位/最大调蓄库容）

单位约定：
  水位 m，流量 m3/s，库容 10^4 m3，时段长 h
  库容换算系数 0.36 = 3600 s/h * 10^4 m3/1e4 ... （详见 read_me）
"""

from dataclasses import dataclass, field, asdict
from typing import List, Optional, Tuple, Dict, Any
import math
import os


# ─────────────────────────────────────────────
# 数据结构
# ─────────────────────────────────────────────
@dataclass
class AuxPoint:
    """辅助曲线数据点：水位~下泄流量~总库容"""
    water_level: float = 0.0   # 水位 (m)
    discharge: float = 0.0     # 下泄流量 (m3/s)
    storage_total: float = 0.0 # 总库容 (10^4 m3)


@dataclass
class AuxResult:
    """辅助曲线计算结果"""
    points: list = field(default_factory=list)   # AuxPoint 列表
    v_above: list = field(default_factory=list)  # 溢洪水位以上库容 (10^4 m3)
    v_dt: list = field(default_factory=list)     # V/dt (m3/s)
    q_half: list = field(default_factory=list)   # q/2 (m3/s)
    sub: list = field(default_factory=list)      # V/dt - q/2 (m3/s)
    add: list = field(default_factory=list)      # V/dt + q/2 (m3/s)
    dt_hours: float = 2.0                        # 计算时段间隔 (h)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["points"] = [asdict(p) for p in self.points]
        return d


@dataclass
class InflowPoint:
    """洪水过程线数据点"""
    period: float = 0.0      # 时段序号
    dt_hours: float = 2.0    # 时段长 (h)
    inflow: float = 0.0      # 来水流量 (m3/s)


@dataclass
class RoutingResult:
    """调洪演算结果"""
    periods: list = field(default_factory=list)
    dt_list: list = field(default_factory=list)
    inflows: list = field(default_factory=list)
    avg_inflows: list = field(default_factory=list)
    sub_prev: list = field(default_factory=list)    # V/dt - q/2 (上时段)
    add_curr: list = field(default_factory=list)    # V/dt + q/2 (本时段)
    outflows: list = field(default_factory=list)
    water_levels: list = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ─────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────
def lin_interp(x1, y1, x2, y2, y_target):
    """线性内插：已知 (x1,y1),(x2,y2)，给定 y_target 求 x"""
    if abs(y2 - y1) < 1e-12:
        return x1
    return x1 + (x2 - x1) / (y2 - y1) * (y_target - y1)


def safe_float(val, default=0.0):
    """安全转换为浮点数"""
    try:
        if val is None or val == '':
            return default
        return float(val)
    except (ValueError, TypeError):
        return default


# ─────────────────────────────────────────────
# 计算引擎
# ─────────────────────────────────────────────
class FloodCalculator:
    """调洪演算核心计算（双辅助曲线法 / 半图解法）"""

    @staticmethod
    def calc_auxiliary(points: List[AuxPoint], dt_hours: float = 2.0) -> AuxResult:
        """
        计算调洪辅助曲线

        参数:
            points: 水位~流量~库容 数据点列表（水位须单调递增，至少2个点）
            dt_hours: 计算时段间隔 (h)
        返回:
            AuxResult

        计算列：
            V_溢洪以上 = V总 - V0                (10^4 m3)
            V/dt       = V_溢洪以上 / (dt*0.36)  (m3/s)
            q/2        = 下泄流量 / 2            (m3/s)
            V/dt - q/2 减辅助曲线               (m3/s)
            V/dt + q/2 加辅助曲线               (m3/s)
        """
        if not points or len(points) < 2:
            raise ValueError("辅助曲线数据至少需要2个点")
        if dt_hours <= 0:
            raise ValueError("时段间隔必须大于0")

        dt = dt_hours * 0.36  # 单位换算系数
        v0 = points[0].storage_total

        result = AuxResult(dt_hours=dt_hours)
        result.points = list(points)

        for p in points:
            v_above = round(p.storage_total - v0, 1)
            vdt = round(v_above / dt, 1) if dt > 0 else 0.0
            qh = round(p.discharge / 2, 1)
            sub_val = round(vdt - qh, 1)
            add_val = round(vdt + qh, 1)

            result.v_above.append(v_above)
            result.v_dt.append(vdt)
            result.q_half.append(qh)
            result.sub.append(sub_val)
            result.add.append(add_val)

        return result

    @staticmethod
    def calc_routing(inflows: List[InflowPoint],
                     aux: AuxResult) -> RoutingResult:
        """
        调洪演算（半图解法逐时段推求）

        参数:
            inflows: 洪水过程线数据点列表（至少2个时段）
            aux: 辅助曲线计算结果
        返回:
            RoutingResult

        递推公式：
            Qp          = (Qi + Qi-1) / 2                 平均来水流量
            (V/dt+q/2)i = Qp + (V/dt-q/2)(i-1)            加辅助曲线递推
            由 (V/dt+q/2) 查加辅助曲线 → q 和 (V/dt-q/2)
            由 q 查水位~流量关系 → 水位 G
        """
        if not inflows or len(inflows) < 2:
            raise ValueError("洪水过程线至少需要2个时段")
        if not aux or len(aux.points) < 2:
            raise ValueError("辅助曲线数据不足")

        n = len(inflows)
        res = RoutingResult()

        # 提取辅助曲线查表数组
        gx = [p.water_level for p in aux.points]
        qx = [p.discharge for p in aux.points]
        svq = list(aux.add)    # V/dt + q/2 (加辅助曲线)
        jvq = list(aux.sub)    # V/dt - q/2 (减辅助曲线)

        # 初始条件
        q_out = [0.0] * n
        qp = [0.0] * n
        sub_arr = [0.0] * n   # V/dt - q/2
        add_arr = [0.0] * n   # V/dt + q/2
        wl = [0.0] * n
        inflow_vals = [0.0] * n
        dt_vals = [0.0] * n
        period_vals = [0.0] * n

        # 初始水位 = 辅助曲线第一个水位（起调水位）
        wl[0] = gx[0]
        inflow_vals[0] = inflows[0].inflow
        dt_vals[0] = inflows[0].dt_hours
        period_vals[0] = inflows[0].period

        for i in range(1, n):
            qi = inflows[i].inflow
            inflow_vals[i] = qi
            dt_vals[i] = inflows[i].dt_hours
            period_vals[i] = inflows[i].period

            # 平均来水流量
            qp[i] = (qi + inflow_vals[i - 1]) / 2.0

            # (V/dt + q/2) = Qp + (V/dt - q/2)_prev
            qsv = qp[i] + sub_arr[i - 1]
            qsv = round(qsv, 2)
            add_arr[i] = qsv

            # 利用 (V/dt+q/2) 查加辅助曲线求 q 和 (V/dt-q/2)
            qx1 = 0.0
            jvq1 = 0.0
            gx1 = gx[0]

            found = False
            for j in range(1, len(svq)):
                if qsv <= svq[j]:
                    # 线性内插
                    qx1 = lin_interp(qx[j - 1], svq[j - 1], qx[j], svq[j], qsv)
                    jvq1 = lin_interp(jvq[j - 1], svq[j - 1], jvq[j], svq[j], qsv)

                    if abs(qx[j] - qx[j - 1]) > 1e-12:
                        gx1 = lin_interp(gx[j - 1], qx[j - 1], gx[j], qx[j], qx1)
                    else:
                        gx1 = gx[j - 1]

                    found = True
                    break

            if not found:
                # 超出辅助曲线范围，取最后一点并外推（警告由调用方提示）
                qx1 = qx[-1]
                jvq1 = jvq[-1]
                gx1 = gx[-1]

            q_out[i] = round(qx1, 2)
            sub_arr[i] = round(jvq1, 2)
            wl[i] = round(gx1, 3)

        res.periods = period_vals
        res.dt_list = dt_vals
        res.inflows = inflow_vals
        res.avg_inflows = qp
        res.sub_prev = sub_arr
        res.add_curr = add_arr
        res.outflows = q_out
        res.water_levels = wl

        return res

    @staticmethod
    def calc_flood(aux_points: List[Tuple[float, float, float]],
                   inflows: List[Tuple[float, float, float]],
                   dt_hours: float = 2.0,
                   wl_col: int = 0, q_col: int = 1, v_col: int = 2,
                   period_col: int = 0, dt_col: int = 1, qin_col: int = 2) -> Dict[str, Any]:
        """
        一站式调洪演算接口

        参数:
            aux_points: [(水位, 流量, 库容), ...] 或 [{'water_level':.., ...}]
            inflows: [(时段序号, 时段长h, 来水流量), ...] 或字典列表
            dt_hours: 时段间隔 (h)，默认 2
        返回:
            dict: {'aux': AuxResult.to_dict(), 'routing': RoutingResult.to_dict(),
                   'summary': {特征值}}
        """
        aux_pts = []
        for p in aux_points:
            if isinstance(p, (tuple, list)):
                aux_pts.append(AuxPoint(water_level=safe_float(p[wl_col]),
                                        discharge=safe_float(p[q_col]),
                                        storage_total=safe_float(p[v_col])))
            else:
                aux_pts.append(AuxPoint(water_level=safe_float(p.get('water_level')),
                                        discharge=safe_float(p.get('discharge')),
                                        storage_total=safe_float(p.get('storage_total'))))

        in_pts = []
        for p in inflows:
            if isinstance(p, (tuple, list)):
                in_pts.append(InflowPoint(period=safe_float(p[period_col]),
                                          dt_hours=safe_float(p[dt_col], 2.0),
                                          inflow=safe_float(p[qin_col])))
            else:
                in_pts.append(InflowPoint(period=safe_float(p.get('period')),
                                          dt_hours=safe_float(p.get('dt_hours'), 2.0),
                                          inflow=safe_float(p.get('inflow'))))

        aux = FloodCalculator.calc_auxiliary(aux_pts, dt_hours)
        routing = FloodCalculator.calc_routing(in_pts, aux)
        summary = FloodCalculator.summary(aux, routing)
        return {'aux': aux.to_dict(), 'routing': routing.to_dict(), 'summary': summary}

    @staticmethod
    def summary(aux: AuxResult, routing: RoutingResult) -> Dict[str, float]:
        """
        提取调洪演算特征值

        返回:
            q_max       最大下泄流量 (m3/s)
            q_max_period 最大下泄流量发生时段
            wl_max      最高水库水位 (m)
            wl_max_period 最高水位发生时段
            inflow_max  最大来水流量 (m3/s)
            v_storage_max 最大调蓄库容 (10^4 m3) ≈ Σ(Qp - q)*dt*0.36 累计峰值
        """
        outflows = routing.outflows
        wls = routing.water_levels
        inflows = routing.inflows
        qp = routing.avg_inflows
        dts = routing.dt_list

        q_max = max(outflows)
        wl_max = max(wls)
        inflow_max = max(inflows)

        # 累计蓄量过程（相对初始），取峰值即最大调蓄库容
        cum = 0.0
        cum_max = 0.0
        for i in range(1, len(qp)):
            cum += (qp[i] - outflows[i]) * (dts[i] * 0.36)
            if cum > cum_max:
                cum_max = cum
        v_storage_max = round(cum_max, 1)

        return {
            'q_max': round(q_max, 2),
            'q_max_period': outflows.index(q_max),
            'wl_max': round(wl_max, 3),
            'wl_max_period': wls.index(wl_max),
            'inflow_max': round(inflow_max, 2),
            'v_storage_max': v_storage_max,
            'dt_hours': aux.dt_hours,
            'start_wl': wls[0],
        }


# ─────────────────────────────────────────────
# Excel 读写（技能内置，无 GUI 依赖）
# ─────────────────────────────────────────────
class ExcelIO:
    """Excel 数据导入导出（.xlsx / .xls）"""

    @staticmethod
    def read_aux_data(filepath: str) -> Tuple[List[AuxPoint], float]:
        """
        从 Excel 读取辅助曲线输入数据
        返回 (AuxPoint列表, dt_hours)
        自动查找名称含"辅助"或"调洪"的工作表，否则用第一个工作表；
        数据列约定：A=水位(m)，B=下泄流量(m3/s)，C=总库容(10^4m3)
        """
        ext = os.path.splitext(filepath)[1].lower()
        if ext == '.xls':
            import xlrd
            wb = xlrd.open_workbook(filepath)
            sheet = None
            for name in wb.sheet_names():
                if '辅助' in name or '调洪' in name:
                    sheet = wb.sheet_by_name(name)
                    break
            if sheet is None:
                sheet = wb.sheet_by_index(0)
            points = []
            for row in range(5, sheet.nrows):
                wl = safe_float(sheet.cell_value(row, 0))
                q = safe_float(sheet.cell_value(row, 1))
                v = safe_float(sheet.cell_value(row, 2))
                if wl >= 0 and q >= 0 and v >= 0 and (wl > 0 or v > 0):
                    points.append(AuxPoint(water_level=wl, discharge=q, storage_total=v))
            return points, 2.0
        else:  # .xlsx
            from openpyxl import load_workbook
            wb = load_workbook(filepath, data_only=True)
            sheet = None
            for name in wb.sheetnames:
                if '辅助' in name:
                    sheet = wb[name]
                    break
            if sheet is None:
                sheet = wb.worksheets[0]
            points = []
            start_row = None
            for row in sheet.iter_rows(min_row=1, max_row=min(20, sheet.max_row), max_col=3):
                vals = [c.value for c in row]
                if all(isinstance(v, (int, float)) and v is not None for v in vals):
                    if start_row is None:
                        start_row = row[0].row
            if start_row is None:
                raise ValueError("未找到有效数据行，请检查 Excel 格式")
            for row in sheet.iter_rows(min_row=start_row, max_col=3):
                wl = safe_float(row[0].value)
                q = safe_float(row[1].value)
                v = safe_float(row[2].value)
                if wl >= 0 and q >= 0 and v >= 0 and (wl > 0 or v > 0):
                    points.append(AuxPoint(water_level=wl, discharge=q, storage_total=v))
            return points, 2.0

    @staticmethod
    def read_inflow_data(filepath: str) -> List[InflowPoint]:
        """从 Excel 读取洪水过程线数据"""
        ext = os.path.splitext(filepath)[1].lower()
        if ext == '.xls':
            import xlrd
            wb = xlrd.open_workbook(filepath)
            sheet = None
            for name in wb.sheet_names():
                if '演算' in name:
                    sheet = wb.sheet_by_name(name)
                    break
            if sheet is None and wb.nsheets > 1:
                sheet = wb.sheet_by_index(1)
            if sheet is None:
                sheet = wb.sheet_by_index(0)
            points = []
            for row in range(2, sheet.nrows):
                period = safe_float(sheet.cell_value(row, 0))
                dt_h = safe_float(sheet.cell_value(row, 1), 2.0)
                q_in = safe_float(sheet.cell_value(row, 2))
                if q_in >= 0:
                    points.append(InflowPoint(period=period, dt_hours=dt_h, inflow=q_in))
            return points
        else:  # .xlsx
            from openpyxl import load_workbook
            wb = load_workbook(filepath, data_only=True)
            sheet = None
            for name in wb.sheetnames:
                if '演算' in name:
                    sheet = wb[name]
                    break
            if sheet is None and len(wb.worksheets) > 1:
                sheet = wb.worksheets[1]
            if sheet is None:
                sheet = wb.worksheets[0]
            points = []
            start_row = None
            for row in sheet.iter_rows(min_row=1, max_row=min(20, sheet.max_row), max_col=3):
                vals = [c.value for c in row]
                if all(isinstance(v, (int, float)) and v is not None for v in vals[:2]):
                    if start_row is None:
                        start_row = row[0].row
            if start_row is None:
                raise ValueError("未找到有效数据行")
            for row in sheet.iter_rows(min_row=start_row, max_col=3):
                period = safe_float(row[0].value)
                dt_h = safe_float(row[1].value, 2.0) if row[1].value is not None else 2.0
                q_in = safe_float(row[2].value)
                if q_in >= 0:
                    points.append(InflowPoint(period=period, dt_hours=dt_h, inflow=q_in))
            return points


# ─────────────────────────────────────────────
# 示例数据（与原 Excel 调洪演算一致）
# ─────────────────────────────────────────────
def default_aux_points() -> List[AuxPoint]:
    """默认辅助曲线数据（水位~流量~库容，来自原 Excel）"""
    return [
        AuxPoint(285.0, 0.0, 3504.0),
        AuxPoint(285.5, 38.1, 3624.0),
        AuxPoint(286.0, 107.6, 3756.0),
        AuxPoint(286.5, 198.2, 3880.0),
        AuxPoint(287.0, 304.8, 4021.0),
        AuxPoint(287.5, 426.0, 4155.0),
        AuxPoint(288.0, 560.0, 4298.0),
        AuxPoint(288.5, 706.0, 4450.0),
        AuxPoint(289.0, 862.0, 4590.0),
        AuxPoint(289.5, 1032.0, 4750.0),
        AuxPoint(290.0, 1205.0, 4895.0),
        AuxPoint(290.5, 1390.0, 5070.0),
        AuxPoint(291.0, 1582.0, 5230.0),
        AuxPoint(291.5, 1795.0, 5415.0),
        AuxPoint(292.0, 1995.0, 5593.0),
    ]


def default_inflow_points() -> List[InflowPoint]:
    """默认洪水过程线数据（48时段，含两次洪水过程）"""
    raw = [
        (0, 2, 0.0), (1, 2, 116.0), (2, 2, 302.0), (3, 2, 278.0), (4, 2, 74.0),
        (5, 2, 33.0), (6, 2, 45.0), (7, 2, 700.0), (8, 2, 400.0), (9, 2, 83.0),
        (10, 2, 34.0), (11, 2, 900.0), (12, 2, 1700.0), (13, 2, 1300.0),
        (14, 2, 900.0), (15, 2, 603.0), (16, 2, 257.0), (17, 2, 50.0),
        (18, 2, 8.0), (19, 2, 0.0),
        (20, 2, 0.0), (21, 2, 0.0), (22, 2, 5.0), (23, 2, 30.0), (24, 2, 120.0),
        (25, 2, 350.0), (26, 2, 580.0), (27, 2, 820.0), (28, 2, 1100.0),
        (29, 2, 950.0), (30, 2, 680.0), (31, 2, 420.0), (32, 2, 250.0),
        (33, 2, 150.0), (34, 2, 80.0), (35, 2, 40.0), (36, 2, 20.0),
        (37, 2, 10.0), (38, 2, 5.0), (39, 2, 0.0),
        (40, 2, 0.0), (41, 2, 0.0), (42, 2, 0.0), (43, 2, 0.0), (44, 2, 0.0),
        (45, 2, 0.0), (46, 2, 0.0), (47, 2, 0.0),
    ]
    return [InflowPoint(period=p, dt_hours=dt, inflow=q) for p, dt, q in raw]


def run_demo():
    """运行示例演算并打印结果（用于验证）"""
    aux = FloodCalculator.calc_auxiliary(default_aux_points(), 2.0)
    routing = FloodCalculator.calc_routing(default_inflow_points(), aux)
    s = FloodCalculator.summary(aux, routing)
    print("=== 调洪演算结果 ===")
    print(f"最大下泄流量 q_max = {s['q_max']} m3/s (时段 {s['q_max_period']})")
    print(f"最高水位 wl_max   = {s['wl_max']} m (时段 {s['wl_max_period']})")
    print(f"最大来水流量      = {s['inflow_max']} m3/s")
    print(f"最大调蓄库容      = {s['v_storage_max']} 10^4 m3")
    print(f"起调水位          = {s['start_wl']} m")
    print(f"时段数            = {len(routing.periods)}")
    return s


if __name__ == '__main__':
    run_demo()
