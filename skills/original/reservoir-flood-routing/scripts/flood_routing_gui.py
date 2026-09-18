#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
水库双辅助曲线调洪演算计算程序 v1.0
============================================
功能：
  1. 调洪辅助曲线计算（水位~库容~流量 → V/dt±q/2 双辅助曲线）
  2. 调洪演算（利用双辅助曲线逐时段推求下泄流量和水库水位）
  3. 自动生成曲线图和过程线图
  4. 支持 Excel 导入和手动输入两种数据录入方式
  5. 计算结果导出为 Excel

算法来源：原 Excel VBA 调洪演算程序（半图解法）
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import os
import sys
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

# ─────────────────────────────────────────────
# 常量与全局设置
# ─────────────────────────────────────────────
APP_TITLE = "水库双辅助曲线调洪演算计算程序 v1.0"
APP_GEOMETRY = "1280x820"

# 尝试设置中文字体
_CN_FONTS = ['Microsoft YaHei', 'SimHei', 'FangSong', 'KaiTi',
             'WenQuanYi Micro Hei', 'Noto Sans CJK SC', 'Arial Unicode MS']

def _set_cn_font():
    """检测并设置可用的中文字体"""
    import matplotlib.font_manager as fm
    available = {f.name for f in fm.fontManager.ttflist}
    for name in _CN_FONTS:
        if name in available:
            plt.rcParams['font.sans-serif'] = [name]
            plt.rcParams['axes.unicode_minus'] = False
            return name
    # 如果没有找到，尝试用系统字体
    plt.rcParams['font.sans-serif'] = ['sans-serif']
    plt.rcParams['axes.unicode_minus'] = False
    return None

_CN_FONT_NAME = _set_cn_font()

# ttk 样式
def _apply_ttk_style():
    style = ttk.Style()
    try:
        style.theme_use('clam')
    except Exception:
        pass
    style.configure('Title.TLabel', font=('Microsoft YaHei', 14, 'bold'))
    style.configure('Section.TLabel', font=('Microsoft YaHei', 10, 'bold'))
    style.configure('Toolbar.TButton', font=('Microsoft YaHei', 9))
    style.configure('Treeview', rowheight=24, font=('Microsoft YaHei', 9))
    style.configure('Treeview.Heading', font=('Microsoft YaHei', 9, 'bold'))


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
# 数据结构
# ─────────────────────────────────────────────
@dataclass
class AuxPoint:
    """辅助曲线数据点"""
    water_level: float = 0.0   # 水位 (m)
    discharge: float = 0.0     # 下泄流量 (m3/s)
    storage_total: float = 0.0 # 总库容 (10^4 m3)

@dataclass
class AuxResult:
    """辅助曲线计算结果"""
    points: list = field(default_factory=list)   # AuxPoint 列表
    v_above: list = field(default_factory=list)  # 溢洪水位以上库容
    v_dt: list = field(default_factory=list)     # V/dt
    q_half: list = field(default_factory=list)   # q/2
    sub: list = field(default_factory=list)      # V/dt - q/2
    add: list = field(default_factory=list)      # V/dt + q/2
    dt_hours: float = 2.0

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
    sub_prev: list = field(default_factory=list)    # V/dt-q/2 (上时段)
    add_curr: list = field(default_factory=list)    # V/dt+q/2 (本时段)
    outflows: list = field(default_factory=list)
    water_levels: list = field(default_factory=list)


# ─────────────────────────────────────────────
# 计算引擎
# ─────────────────────────────────────────────
class FloodCalculator:
    """调洪演算核心计算"""

    @staticmethod
    def calc_auxiliary(points: List[AuxPoint], dt_hours: float) -> AuxResult:
        """
        计算调洪辅助曲线
        参数:
            points: 水位-流量-库容 数据点列表
            dt_hours: 计算时段间隔 (小时)
        返回:
            AuxResult 包含所有计算列
        """
        if not points or len(points) < 2:
            raise ValueError("辅助曲线数据至少需要2个点")
        if dt_hours <= 0:
            raise ValueError("时段间隔必须大于0")

        dt = dt_hours * 0.36  # 转换系数: (10^4 m3) / (h * 10000/3600) = 0.36
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
        调洪演算计算（半图解法）
        参数:
            inflows: 洪水过程线数据点列表
            aux: 辅助曲线计算结果
        返回:
            RoutingResult
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
        svq = list(aux.add)    # V/dt + q/2
        jvq = list(aux.sub)    # V/dt - q/2

        # 初始条件
        q_out = [0.0] * n
        qp = [0.0] * n
        sub_arr = [0.0] * n   # V/dt - q/2
        add_arr = [0.0] * n   # V/dt + q/2
        wl = [0.0] * n
        inflow_vals = [0.0] * n
        dt_vals = [0.0] * n
        period_vals = [0.0] * n

        # 初始水位 = 辅助曲线第一个水位
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

            # V/dt + q/2 = Qp + (V/dt - q/2)_prev
            qsv = qp[i] + sub_arr[i - 1]
            qsv = round(qsv, 2)
            add_arr[i] = qsv

            # 利用 V/dt+q/2 查辅助曲线求 q 和 V/dt-q/2
            qx1 = 0.0
            jvq1 = 0.0
            gx1 = gx[0]

            found = False
            for j in range(1, len(svq)):
                if qsv <= svq[j]:
                    # 线性内插
                    qx1 = lin_interp(qx[j-1], svq[j-1], qx[j], svq[j], qsv)
                    jvq1 = lin_interp(jvq[j-1], svq[j-1], jvq[j], svq[j], qsv)

                    if abs(qx[j] - qx[j-1]) > 1e-12:
                        gx1 = lin_interp(gx[j-1], qx[j-1], gx[j], qx[j], qx1)
                    else:
                        gx1 = gx[j-1]

                    found = True
                    break

            if not found:
                # 超出范围，取最后一个点并外推
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


# ─────────────────────────────────────────────
# Excel 读写
# ─────────────────────────────────────────────
class ExcelIO:
    """Excel 数据导入导出"""

    @staticmethod
    def read_aux_data(filepath: str) -> Tuple[List[AuxPoint], float]:
        """
        从 Excel 读取辅助曲线输入数据
        返回 (AuxPoint列表, dt_hours)
        支持 .xlsx 和 .xls
        """
        ext = os.path.splitext(filepath)[1].lower()

        if ext == '.xls':
            import xlrd
            wb = xlrd.open_workbook(filepath)
            # 查找调洪辅助曲线计算表
            sheet = None
            for name in wb.sheet_names():
                if '辅助' in name or '调洪' in name:
                    sheet = wb.sheet_by_name(name)
                    break
            if sheet is None:
                sheet = wb.sheet_by_index(0)

            points = []
            # 数据从第6行开始（索引5），表头在2-5行
            for row in range(5, sheet.nrows):
                wl = safe_float(sheet.cell_value(row, 0))
                q = safe_float(sheet.cell_value(row, 1))
                v = safe_float(sheet.cell_value(row, 2))
                if wl >= 0 and q >= 0 and v >= 0 and (wl > 0 or v > 0):
                    points.append(AuxPoint(water_level=wl, discharge=q, storage_total=v))

            # 尝试从表格中读取 dt（如果有的话）
            dt_hours = 2.0
            return points, dt_hours

        else:  # .xlsx
            from openpyxl import load_workbook
            wb = load_workbook(filepath, data_only=True)
            # 查找工作表
            sheet = None
            for name in wb.sheetnames:
                if '辅助' in name:
                    sheet = wb[name]
                    break
            if sheet is None:
                sheet = wb.worksheets[0]

            points = []
            # 自动检测数据起始行
            start_row = None
            for row in sheet.iter_rows(min_row=1, max_row=min(20, sheet.max_row), max_col=3):
                vals = [c.value for c in row]
                # 找第一个数据行（三个值都是数字）
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

            dt_hours = 2.0
            return points, dt_hours

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
            for row in range(2, sheet.nrows):  # 从第3行开始(索引2)
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

    @staticmethod
    def export_results(filepath: str, aux: AuxResult, routing: RoutingResult):
        """导出计算结果到 Excel"""
        from openpyxl import Workbook
        from openpyxl.styles import Font, Alignment, Border, Side, PatternFill

        wb = Workbook()

        # ── 辅助曲线计算表 ──
        ws1 = wb.active
        ws1.title = "调洪辅助曲线计算表"

        # 表头
        headers1 = ['水位(m)', '下泄流量(m3/s)', '总库容(10^4m3)',
                     '溢洪水位以上库容(10^4m3)', 'V/dt(m3/s)', 'q/2(m3/s)',
                     'V/dt-q/2(m3/s)', 'V/dt+q/2(m3/s)']
        header_font = Font(bold=True, size=10)
        header_fill = PatternFill(start_color='D9E1F2', end_color='D9E1F2', fill_type='solid')
        thin_border = Border(
            left=Side(style='thin'), right=Side(style='thin'),
            top=Side(style='thin'), bottom=Side(style='thin')
        )

        # 标题行
        ws1.merge_cells('A1:H1')
        ws1['A1'] = f'调洪辅助曲线计算表 (dt={aux.dt_hours}h)'
        ws1['A1'].font = Font(bold=True, size=14)
        ws1['A1'].alignment = Alignment(horizontal='center')

        for col, h in enumerate(headers1, 1):
            cell = ws1.cell(row=2, column=col, value=h)
            cell.font = header_font
            cell.fill = header_fill
            cell.border = thin_border
            cell.alignment = Alignment(horizontal='center', wrap_text=True)

        # 数据行
        for i, p in enumerate(aux.points):
            row = i + 3
            data = [p.water_level, p.discharge, p.storage_total,
                    aux.v_above[i], aux.v_dt[i], aux.q_half[i],
                    aux.sub[i], aux.add[i]]
            for col, val in enumerate(data, 1):
                cell = ws1.cell(row=row, column=col, value=val)
                cell.border = thin_border
                cell.alignment = Alignment(horizontal='center')

        # 列宽
        for col in range(1, 9):
            ws1.column_dimensions[chr(64 + col)].width = 18

        # ── 调洪演算计算表 ──
        ws2 = wb.create_sheet("调洪演算计算表")
        headers2 = ['时段序号', '时段长dt(h)', '来水流量Q(m3/s)', '平均流量Qp(m3/s)',
                     'V/dt-q/2(m3/s)', 'V/dt+q/2(m3/s)', '下泄流量q(m3/s)', '水库水位(m)']

        ws2.merge_cells('A1:H1')
        ws2['A1'] = '调洪演算计算表'
        ws2['A1'].font = Font(bold=True, size=14)
        ws2['A1'].alignment = Alignment(horizontal='center')

        for col, h in enumerate(headers2, 1):
            cell = ws2.cell(row=2, column=col, value=h)
            cell.font = header_font
            cell.fill = header_fill
            cell.border = thin_border
            cell.alignment = Alignment(horizontal='center', wrap_text=True)

        for i in range(len(routing.periods)):
            row = i + 3
            data = [routing.periods[i], routing.dt_list[i], routing.inflows[i],
                    routing.avg_inflows[i], routing.sub_prev[i], routing.add_curr[i],
                    routing.outflows[i], routing.water_levels[i]]
            for col, val in enumerate(data, 1):
                cell = ws2.cell(row=row, column=col, value=val)
                cell.border = thin_border
                cell.alignment = Alignment(horizontal='center')

        for col in range(1, 9):
            ws2.column_dimensions[chr(64 + col)].width = 18

        wb.save(filepath)


# ─────────────────────────────────────────────
# GUI 应用程序
# ─────────────────────────────────────────────
class FloodRoutingApp:
    """主应用程序 GUI"""

    def __init__(self):
        self.root = tk.Tk()
        self.root.title(APP_TITLE)
        self.root.geometry(APP_GEOMETRY)
        self.root.minsize(1024, 700)

        # 数据
        self.aux_data: List[AuxPoint] = []
        self.aux_result: Optional[AuxResult] = None
        self.inflow_data: List[InflowPoint] = []
        self.routing_result: Optional[RoutingResult] = None

        # 图表
        self.fig_aux = None
        self.fig_routing = None

        self._build_ui()
        self._load_default_data()

    # ═══════════════════════════════════════════
    # UI 构建
    # ═══════════════════════════════════════════
    def _build_ui(self):
        # 主 Notebook
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill='both', expand=True, padx=4, pady=4)

        # 4 个标签页
        self.tab_aux = ttk.Frame(self.notebook)
        self.tab_routing = ttk.Frame(self.notebook)
        self.tab_charts = ttk.Frame(self.notebook)
        self.tab_about = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_aux, text='  调洪辅助曲线计算  ')
        self.notebook.add(self.tab_routing, text='  调洪演算计算  ')
        self.notebook.add(self.tab_charts, text='  图表展示  ')
        self.notebook.add(self.tab_about, text='  程序说明  ')

        self._build_aux_tab()
        self._build_routing_tab()
        self._build_charts_tab()
        self._build_about_tab()

        # 底部状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(self.root, textvariable=self.status_var,
                               relief='sunken', anchor='w', padding=(8, 2))
        status_bar.pack(side='bottom', fill='x')

    def _build_aux_tab(self):
        """构建辅助曲线计算标签页"""
        frame = self.tab_aux

        # ── 顶部参数栏 ──
        top = ttk.Frame(frame, padding=6)
        top.pack(fill='x')

        ttk.Label(top, text="计算时段间隔 dt (h):").pack(side='left')
        self.aux_dt_var = tk.StringVar(value="2")
        dt_entry = ttk.Entry(top, textvariable=self.aux_dt_var, width=8)
        dt_entry.pack(side='left', padx=(4, 20))

        # ── 工具栏 ──
        toolbar = ttk.Frame(frame, padding=(6, 0))
        toolbar.pack(fill='x')

        btn_style = {'padding': (8, 2)}
        ttk.Button(toolbar, text="导入Excel", command=self._import_aux_excel,
                   style='Toolbar.TButton', **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="添加行", command=self._add_aux_row,
                   **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="删除末行", command=self._del_aux_row,
                   **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="清空数据", command=self._clear_aux_data,
                   **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="批量输入", command=self._batch_input_aux,
                   **btn_style).pack(side='left', padx=2)

        sep = ttk.Separator(toolbar, orient='vertical')
        sep.pack(side='left', fill='y', padx=10)

        self.btn_calc_aux = ttk.Button(toolbar, text="▶ 计算辅助曲线",
                                        command=self._calculate_aux,
                                        **btn_style)
        self.btn_calc_aux.pack(side='left', padx=2)

        # ── 数据表格 ──
        table_frame = ttk.Frame(frame, padding=4)
        table_frame.pack(fill='both', expand=True)

        columns = ('water_level', 'discharge', 'storage',
                   'v_above', 'v_dt', 'q_half', 'sub', 'add')
        col_names = ('水位\n(m)', '下泄流量\n(m3/s)', '总库容\n(10^4m3)',
                     '溢洪以上库容\n(10^4m3)', 'V/dt\n(m3/s)', 'q/2\n(m3/s)',
                     'V/dt-q/2\n(m3/s)', 'V/dt+q/2\n(m3/s)')
        col_widths = (90, 100, 100, 120, 90, 80, 110, 110)

        self.aux_tree = ttk.Treeview(table_frame, columns=columns,
                                      show='headings', selectmode='extended')
        for col, name, w in zip(columns, col_names, col_widths):
            self.aux_tree.heading(col, text=name)
            editable = col in ('water_level', 'discharge', 'storage')
            self.aux_tree.column(col, width=w, anchor='center')

        # 滚动条
        vsb = ttk.Scrollbar(table_frame, orient='vertical', command=self.aux_tree.yview)
        self.aux_tree.configure(yscrollcommand=vsb.set)

        self.aux_tree.pack(side='left', fill='both', expand=True)
        vsb.pack(side='right', fill='y')

        # 双击编辑
        self.aux_tree.bind('<Double-1>', self._on_aux_dblclick)

        # ── 提示 ──
        hint = ttk.Label(frame, text="提示：双击前3列可编辑 | 批量输入可粘贴多行 | 导入Excel支持.xlsx/.xls | 数据行数不限",
                         foreground='gray', padding=(6, 2))
        hint.pack(fill='x')

    def _build_routing_tab(self):
        """构建调洪演算计算标签页"""
        frame = self.tab_routing

        # ── 工具栏 ──
        toolbar = ttk.Frame(frame, padding=6)
        toolbar.pack(fill='x')

        btn_style = {'padding': (8, 2)}
        ttk.Button(toolbar, text="导入Excel", command=self._import_inflow_excel,
                   **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="添加行", command=self._add_inflow_row,
                   **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="删除末行", command=self._del_inflow_row,
                   **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="清空数据", command=self._clear_inflow_data,
                   **btn_style).pack(side='left', padx=2)
        ttk.Button(toolbar, text="批量输入", command=self._batch_input_inflow,
                   **btn_style).pack(side='left', padx=2)

        sep = ttk.Separator(toolbar, orient='vertical')
        sep.pack(side='left', fill='y', padx=10)

        self.btn_calc_routing = ttk.Button(toolbar, text="▶ 执行调洪演算",
                                            command=self._calculate_routing,
                                            **btn_style)
        self.btn_calc_routing.pack(side='left', padx=2)

        ttk.Button(toolbar, text="导出结果到Excel", command=self._export_results,
                   **btn_style).pack(side='left', padx=12)

        # ── 数据表格 ──
        table_frame = ttk.Frame(frame, padding=4)
        table_frame.pack(fill='both', expand=True)

        columns = ('period', 'dt', 'inflow', 'avg_inflow',
                   'sub_prev', 'add_curr', 'outflow', 'water_level')
        col_names = ('时段\n序号', '时段长\ndt(h)', '来水流量\nQ(m3/s)', '平均流量\nQp(m3/s)',
                     'V/dt-q/2\n(m3/s)', 'V/dt+q/2\n(m3/s)', '下泄流量\nq(m3/s)', '水库水位\n(m)')
        col_widths = (70, 80, 100, 100, 110, 110, 110, 90)

        self.routing_tree = ttk.Treeview(table_frame, columns=columns,
                                          show='headings', selectmode='extended')
        for col, name, w in zip(columns, col_names, col_widths):
            self.routing_tree.heading(col, text=name)
            self.routing_tree.column(col, width=w, anchor='center')

        vsb = ttk.Scrollbar(table_frame, orient='vertical', command=self.routing_tree.yview)
        self.routing_tree.configure(yscrollcommand=vsb.set)

        self.routing_tree.pack(side='left', fill='both', expand=True)
        vsb.pack(side='right', fill='y')

        self.routing_tree.bind('<Double-1>', self._on_routing_dblclick)

        # ── 提示 ──
        hint = ttk.Label(frame,
                         text="提示：双击前3列可编辑 | 批量输入可粘贴多行 | 时段数量不限 | 请先在“辅助曲线”页完成计算再执行演算",
                         foreground='gray', padding=(6, 2))
        hint.pack(fill='x')

    def _build_charts_tab(self):
        """构建图表展示标签页"""
        frame = self.tab_charts

        # 图表工具栏
        toolbar = ttk.Frame(frame, padding=6)
        toolbar.pack(fill='x')

        ttk.Button(toolbar, text="刷新图表", command=self._refresh_charts,
                   padding=(8, 2)).pack(side='left', padx=2)
        ttk.Button(toolbar, text="保存图片", command=self._save_charts,
                   padding=(8, 2)).pack(side='left', padx=2)

        # 图表区域（上下分割）
        self.chart_container = ttk.PanedWindow(frame, orient='vertical')
        self.chart_container.pack(fill='both', expand=True, padx=4, pady=4)

        # 上半部分：辅助曲线图
        self.chart_top_frame = ttk.Frame(self.chart_container)
        self.chart_container.add(self.chart_top_frame, weight=1)

        # 下半部分：演算过程线图
        self.chart_bottom_frame = ttk.Frame(self.chart_container)
        self.chart_container.add(self.chart_bottom_frame, weight=1)

        # 初始化空白图表
        self._init_blank_charts()

    def _build_about_tab(self):
        """构建程序说明标签页"""
        frame = self.tab_about

        text_widget = tk.Text(frame, wrap='word', font=('Microsoft YaHei', 11),
                              padx=20, pady=10)
        text_widget.pack(fill='both', expand=True)

        about_text = """
水库双辅助曲线调洪演算计算程序 v1.0
============================================

一、程序功能

  本程序用于水库调洪演算计算，采用"双辅助曲线法"（半图解法），
  通过 V/dt+q/2 和 V/dt-q/2 两条辅助曲线，逐时段推求水库下泄
  流量和相应水位。

二、计算原理

  辅助曲线计算：
    V = V总 - V0（溢洪水位以上库容）
    V/dt = V / (dt × 0.36)    （单位换算：10^4m3 → m3/s）
    辅助曲线1：V/dt - q/2（减辅助曲线）
    辅助曲线2：V/dt + q/2（加辅助曲线）

  调洪演算（逐时段）：
    Qp = (Qi + Qi-1) / 2      （平均来水流量）
    V/dt + q/2 = Qp + (V/dt - q/2)上时段
    由 V/dt+q/2 查辅助曲线 → 得 q 和 V/dt-q/2
    由 q 查水位~流量关系 → 得水位 G

三、使用方法

  步骤1：在"调洪辅助曲线计算"页输入（或导入）水位~流量~库容数据
  步骤2：设置计算时段间隔 dt，点击"计算辅助曲线"
  步骤3：在"调洪演算计算"页输入（或导入）洪水过程线数据
  步骤4：点击"执行调洪演算"
  步骤5：在"图表展示"页查看曲线和过程线图
  步骤6：点击"导出结果到Excel"保存计算成果

四、数据输入方式

  方式一：直接在表格中手动输入（双击单元格编辑）
  方式二：批量输入（点击"批量输入"按钮，粘贴多行文本）
  方式三：从 Excel 文件导入（支持 .xlsx 和 .xls 格式）
  时段数量不限，可输入任意多个时段。

五、注意事项

  1. 辅助曲线数据的水位必须单调递增
  2. 洪水过程线第0行为初始状态（来水流量=0，下泄流量=0）
  3. 时段长 dt 可以各时段不同（手动输入时逐行设置）
  4. 计算结果保留精度：流量2位小数，水位3位小数
  5. 时段数量无上限，支持任意长洪水过程

六、版本信息

  版本：v1.0
  基于原 Excel VBA 调洪演算程序算法
  Python + tkinter + matplotlib 实现
"""
        text_widget.insert('1.0', about_text)
        text_widget.config(state='disabled')

    def _init_blank_charts(self):
        """初始化空白图表"""
        # 上图：辅助曲线（3个子图）
        self.fig_aux = Figure(figsize=(12, 4), dpi=100)
        self.fig_aux.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.15, wspace=0.25)
        self.ax1 = self.fig_aux.add_subplot(131)
        self.ax2 = self.fig_aux.add_subplot(132)
        self.ax3 = self.fig_aux.add_subplot(133)
        self.ax1.set_title('Water Level - Storage')
        self.ax2.set_title('Water Level - Discharge')
        self.ax3.set_title('Dual Auxiliary Curves')

        self.canvas_top = FigureCanvasTkAgg(self.fig_aux, master=self.chart_top_frame)
        self.canvas_top.get_tk_widget().pack(fill='both', expand=True)

        # 下图：演算过程线（2个子图）
        self.fig_routing = Figure(figsize=(12, 4), dpi=100)
        self.fig_routing.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.15, wspace=0.25)
        self.ax4 = self.fig_routing.add_subplot(121)
        self.ax5 = self.fig_routing.add_subplot(122)
        self.ax4.set_title('Outflow Hydrograph')
        self.ax5.set_title('Water Level Process')

        self.canvas_bottom = FigureCanvasTkAgg(self.fig_routing, master=self.chart_bottom_frame)
        self.canvas_bottom.get_tk_widget().pack(fill='both', expand=True)

    # ═══════════════════════════════════════════
    # 默认数据加载
    # ═══════════════════════════════════════════
    def _load_default_data(self):
        """加载示例数据"""
        # 辅助曲线数据（来自原Excel）
        self.aux_data = [
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

        # 洪水过程线数据（来自原Excel + 扩展退水段，共48时段）
        # 时段数量不限，可任意增加
        self.inflow_data = [
            InflowPoint(0, 2, 0.0),
            InflowPoint(1, 2, 116.0),
            InflowPoint(2, 2, 302.0),
            InflowPoint(3, 2, 278.0),
            InflowPoint(4, 2, 74.0),
            InflowPoint(5, 2, 33.0),
            InflowPoint(6, 2, 45.0),
            InflowPoint(7, 2, 700.0),
            InflowPoint(8, 2, 400.0),
            InflowPoint(9, 2, 83.0),
            InflowPoint(10, 2, 34.0),
            InflowPoint(11, 2, 900.0),
            InflowPoint(12, 2, 1700.0),
            InflowPoint(13, 2, 1300.0),
            InflowPoint(14, 2, 900.0),
            InflowPoint(15, 2, 603.0),
            InflowPoint(16, 2, 257.0),
            InflowPoint(17, 2, 50.0),
            InflowPoint(18, 2, 8.0),
            InflowPoint(19, 2, 0.0),
            # ── 以下为扩展退水段（第2次洪水过程示例） ──
            InflowPoint(20, 2, 0.0),
            InflowPoint(21, 2, 0.0),
            InflowPoint(22, 2, 5.0),
            InflowPoint(23, 2, 30.0),
            InflowPoint(24, 2, 120.0),
            InflowPoint(25, 2, 350.0),
            InflowPoint(26, 2, 580.0),
            InflowPoint(27, 2, 820.0),
            InflowPoint(28, 2, 1100.0),
            InflowPoint(29, 2, 950.0),
            InflowPoint(30, 2, 680.0),
            InflowPoint(31, 2, 420.0),
            InflowPoint(32, 2, 250.0),
            InflowPoint(33, 2, 150.0),
            InflowPoint(34, 2, 80.0),
            InflowPoint(35, 2, 40.0),
            InflowPoint(36, 2, 20.0),
            InflowPoint(37, 2, 10.0),
            InflowPoint(38, 2, 5.0),
            InflowPoint(39, 2, 0.0),
            InflowPoint(40, 2, 0.0),
            InflowPoint(41, 2, 0.0),
            InflowPoint(42, 2, 0.0),
            InflowPoint(43, 2, 0.0),
            InflowPoint(44, 2, 0.0),
            InflowPoint(45, 2, 0.0),
            InflowPoint(46, 2, 0.0),
            InflowPoint(47, 2, 0.0),
        ]

        self._refresh_aux_tree()
        self._refresh_inflow_tree()
        self.status_var.set("已加载示例数据。可修改或导入新数据。")

    # ═══════════════════════════════════════════
    # 表格操作 - 辅助曲线
    # ═══════════════════════════════════════════
    def _refresh_aux_tree(self):
        """刷新辅助曲线表格"""
        for item in self.aux_tree.get_children():
            self.aux_tree.delete(item)

        if self.aux_result:
            for i, p in enumerate(self.aux_data):
                self.aux_tree.insert('', 'end', values=(
                    p.water_level, p.discharge, p.storage_total,
                    self.aux_result.v_above[i],
                    self.aux_result.v_dt[i],
                    self.aux_result.q_half[i],
                    self.aux_result.sub[i],
                    self.aux_result.add[i]
                ))
        else:
            for p in self.aux_data:
                self.aux_tree.insert('', 'end', values=(
                    p.water_level, p.discharge, p.storage_total,
                    '-', '-', '-', '-', '-'
                ))

    def _add_aux_row(self):
        """添加辅助曲线数据行"""
        # 默认值基于最后一行
        if self.aux_data:
            last = self.aux_data[-1]
            new_wl = last.water_level + 0.5
            new_q = 0.0
            new_v = last.storage_total + 100
        else:
            new_wl, new_q, new_v = 0.0, 0.0, 0.0

        self.aux_data.append(AuxPoint(new_wl, new_q, new_v))
        self.aux_result = None  # 清除计算结果
        self._refresh_aux_tree()

    def _del_aux_row(self):
        """删除辅助曲线末行"""
        if self.aux_data:
            self.aux_data.pop()
            self.aux_result = None
            self._refresh_aux_tree()

    def _clear_aux_data(self):
        """清空辅助曲线数据"""
        if messagebox.askyesno("确认", "确定清空所有辅助曲线数据？"):
            self.aux_data.clear()
            self.aux_result = None
            self._refresh_aux_tree()

    def _on_aux_dblclick(self, event):
        """双击编辑辅助曲线表格"""
        region = self.aux_tree.identify("region", event.x, event.y)
        if region != "cell":
            return

        col_id = self.aux_tree.identify_column(event.x)
        col_idx = int(col_id.replace('#', '')) - 1

        # 只允许编辑前3列
        if col_idx > 2:
            return

        item = self.aux_tree.identify_row(event.y)
        if not item:
            return

        row_idx = self.aux_tree.index(item)
        old_val = self.aux_tree.set(item, col_id)

        # 弹出输入框
        col_names = ['水位 (m)', '下泄流量 (m3/s)', '总库容 (10^4m3)']
        new_val = simpledialog.askfloat("编辑数据", f"{col_names[col_idx]}:",
                                         initialvalue=float(old_val) if old_val != '-' else 0.0,
                                         parent=self.root)
        if new_val is not None:
            if col_idx == 0:
                self.aux_data[row_idx].water_level = new_val
            elif col_idx == 1:
                self.aux_data[row_idx].discharge = new_val
            elif col_idx == 2:
                self.aux_data[row_idx].storage_total = new_val
            self.aux_result = None
            self._refresh_aux_tree()

    # ═══════════════════════════════════════════
    # 表格操作 - 洪水过程线
    # ═══════════════════════════════════════════
    def _refresh_inflow_tree(self):
        """刷新洪水过程线表格"""
        for item in self.routing_tree.get_children():
            self.routing_tree.delete(item)

        if self.routing_result:
            for i in range(len(self.routing_result.periods)):
                self.routing_tree.insert('', 'end', values=(
                    self.routing_result.periods[i],
                    self.routing_result.dt_list[i],
                    self.routing_result.inflows[i],
                    round(self.routing_result.avg_inflows[i], 2),
                    round(self.routing_result.sub_prev[i], 2),
                    round(self.routing_result.add_curr[i], 2),
                    round(self.routing_result.outflows[i], 2),
                    round(self.routing_result.water_levels[i], 3)
                ))
        else:
            for p in self.inflow_data:
                self.routing_tree.insert('', 'end', values=(
                    p.period, p.dt_hours, p.inflow,
                    '-', '-', '-', '-', '-'
                ))

    def _add_inflow_row(self):
        """添加洪水过程线行"""
        if self.inflow_data:
            last = self.inflow_data[-1]
            new_period = last.period + 1
            new_dt = last.dt_hours
        else:
            new_period, new_dt = 0, 2.0

        self.inflow_data.append(InflowPoint(new_period, new_dt, 0.0))
        self.routing_result = None
        self._refresh_inflow_tree()

    def _del_inflow_row(self):
        """删除末行"""
        if self.inflow_data:
            self.inflow_data.pop()
            self.routing_result = None
            self._refresh_inflow_tree()

    def _clear_inflow_data(self):
        """清空数据"""
        if messagebox.askyesno("确认", "确定清空所有洪水过程线数据？"):
            self.inflow_data.clear()
            self.routing_result = None
            self._refresh_inflow_tree()

    def _on_routing_dblclick(self, event):
        """双击编辑洪水过程线表格"""
        region = self.routing_tree.identify("region", event.x, event.y)
        if region != "cell":
            return

        col_id = self.routing_tree.identify_column(event.x)
        col_idx = int(col_id.replace('#', '')) - 1

        if col_idx > 2:
            return

        item = self.routing_tree.identify_row(event.y)
        if not item:
            return

        row_idx = self.routing_tree.index(item)
        old_val = self.routing_tree.set(item, col_id)

        col_names = ['时段序号', '时段长 dt (h)', '来水流量 Q (m3/s)']
        new_val = simpledialog.askfloat("编辑数据", f"{col_names[col_idx]}:",
                                         initialvalue=float(old_val) if old_val != '-' else 0.0,
                                         parent=self.root)
        if new_val is not None:
            if col_idx == 0:
                self.inflow_data[row_idx].period = new_val
            elif col_idx == 1:
                self.inflow_data[row_idx].dt_hours = new_val
            elif col_idx == 2:
                self.inflow_data[row_idx].inflow = new_val
            self.routing_result = None
            self._refresh_inflow_tree()

    # ═══════════════════════════════════════════
    # 批量输入
    # ═══════════════════════════════════════════
    def _batch_input_aux(self):
        """批量输入辅助曲线数据（粘贴多行）"""
        dialog = tk.Toplevel(self.root)
        dialog.title("批量输入 - 辅助曲线数据")
        dialog.geometry("600x500")
        dialog.transient(self.root)

        ttk.Label(dialog,
                  text="每行一个数据点，列之间用空格、Tab或逗号分隔。\n"
                       "格式：水位(m)  下泄流量(m3/s)  总库容(10^4m3)\n"
                       "示例：285.0  0.0  3504.0",
                  foreground='gray', padding=8).pack(fill='x')

        text = tk.Text(dialog, font=('Consolas', 10), wrap='none')
        text.pack(fill='both', expand=True, padx=8, pady=(0, 4))

        scrollbar = ttk.Scrollbar(text, orient='vertical', command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y')

        def do_import():
            raw = text.get("1.0", "end").strip()
            if not raw:
                messagebox.showwarning("提示", "未输入任何数据", parent=dialog)
                return
            new_points = []
            errors = []
            for line_no, line in enumerate(raw.splitlines(), 1):
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('//'):
                    continue
                parts = line.replace(',', ' ').replace('\t', ' ').split()
                if len(parts) < 3:
                    errors.append(f"第{line_no}行：需要3列数据，实际{len(parts)}列")
                    continue
                try:
                    wl, q, v = float(parts[0]), float(parts[1]), float(parts[2])
                    new_points.append(AuxPoint(wl, q, v))
                except ValueError:
                    errors.append(f"第{line_no}行：数值解析失败")

            if errors:
                messagebox.showwarning("部分行解析失败", "\n".join(errors[:10]),
                                        parent=dialog)
            if new_points:
                self.aux_data = new_points
                self.aux_result = None
                self._refresh_aux_tree()
                self.status_var.set(f"批量输入完成：{len(new_points)} 行辅助曲线数据")
                dialog.destroy()
            elif not errors:
                messagebox.showwarning("提示", "未解析到有效数据", parent=dialog)

        btn_frame = ttk.Frame(dialog, padding=8)
        btn_frame.pack(fill='x')
        ttk.Button(btn_frame, text="确认导入", command=do_import).pack(side='left', padx=4)
        ttk.Button(btn_frame, text="取消", command=dialog.destroy).pack(side='left', padx=4)

        dialog.grab_set()

    def _batch_input_inflow(self):
        """批量输入洪水过程线数据（粘贴多行）"""
        dialog = tk.Toplevel(self.root)
        dialog.title("批量输入 - 洪水过程线数据")
        dialog.geometry("600x500")
        dialog.transient(self.root)

        ttk.Label(dialog,
                  text="每行一个时段，列之间用空格、Tab或逗号分隔。\n"
                       "格式：时段序号  时段长dt(h)  来水流量Q(m3/s)\n"
                       "示例：0  2  0.0\n"
                       "       1  2  116.0\n"
                       "时段数量不限，可输入任意多行。",
                  foreground='gray', padding=8).pack(fill='x')

        text = tk.Text(dialog, font=('Consolas', 10), wrap='none')
        text.pack(fill='both', expand=True, padx=8, pady=(0, 4))

        scrollbar = ttk.Scrollbar(text, orient='vertical', command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side='right', fill='y')

        def do_import():
            raw = text.get("1.0", "end").strip()
            if not raw:
                messagebox.showwarning("提示", "未输入任何数据", parent=dialog)
                return
            new_points = []
            errors = []
            for line_no, line in enumerate(raw.splitlines(), 1):
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('//'):
                    continue
                parts = line.replace(',', ' ').replace('\t', ' ').split()
                if len(parts) < 3:
                    # 也支持只输入2列（时段长+流量），自动编号
                    if len(parts) == 2:
                        try:
                            dt_h = float(parts[0])
                            q_in = float(parts[1])
                            period = len(new_points)
                            new_points.append(InflowPoint(period, dt_h, q_in))
                            continue
                        except ValueError:
                            pass
                    errors.append(f"第{line_no}行：需要至少3列数据，实际{len(parts)}列")
                    continue
                try:
                    period = float(parts[0])
                    dt_h = float(parts[1])
                    q_in = float(parts[2])
                    new_points.append(InflowPoint(period, dt_h, q_in))
                except ValueError:
                    errors.append(f"第{line_no}行：数值解析失败")

            if errors:
                messagebox.showwarning("部分行解析失败", "\n".join(errors[:10]),
                                        parent=dialog)
            if new_points:
                self.inflow_data = new_points
                self.routing_result = None
                self._refresh_inflow_tree()
                self.status_var.set(f"批量输入完成：{len(new_points)} 个时段（时段数量不限）")
                dialog.destroy()
            elif not errors:
                messagebox.showwarning("提示", "未解析到有效数据", parent=dialog)

        btn_frame = ttk.Frame(dialog, padding=8)
        btn_frame.pack(fill='x')
        ttk.Button(btn_frame, text="确认导入", command=do_import).pack(side='left', padx=4)
        ttk.Button(btn_frame, text="取消", command=dialog.destroy).pack(side='left', padx=4)

        dialog.grab_set()

    # ═══════════════════════════════════════════
    # Excel 导入
    # ═══════════════════════════════════════════
    def _import_aux_excel(self):
        """从 Excel 导入辅助曲线数据"""
        filepath = filedialog.askopenfilename(
            title="选择辅助曲线 Excel 文件",
            filetypes=[("Excel 文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        )
        if not filepath:
            return

        try:
            points, dt = ExcelIO.read_aux_data(filepath)
            if not points:
                messagebox.showwarning("警告", "未从文件中读取到有效数据")
                return

            self.aux_data = points
            self.aux_dt_var.set(str(dt))
            self.aux_result = None
            self._refresh_aux_tree()
            self.status_var.set(f"已导入 {len(points)} 行辅助曲线数据")
            messagebox.showinfo("导入成功", f"已导入 {len(points)} 行数据")

        except Exception as e:
            messagebox.showerror("导入错误", f"读取 Excel 文件失败:\n{str(e)}")

    def _import_inflow_excel(self):
        """从 Excel 导入洪水过程线数据"""
        filepath = filedialog.askopenfilename(
            title="选择洪水过程线 Excel 文件",
            filetypes=[("Excel 文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        )
        if not filepath:
            return

        try:
            points = ExcelIO.read_inflow_data(filepath)
            if not points:
                messagebox.showwarning("警告", "未从文件中读取到有效数据")
                return

            self.inflow_data = points
            self.routing_result = None
            self._refresh_inflow_tree()
            self.status_var.set(f"已导入 {len(points)} 行洪水过程线数据")
            messagebox.showinfo("导入成功", f"已导入 {len(points)} 行数据")

        except Exception as e:
            messagebox.showerror("导入错误", f"读取 Excel 文件失败:\n{str(e)}")

    # ═══════════════════════════════════════════
    # 计算执行
    # ═══════════════════════════════════════════
    def _calculate_aux(self):
        """执行辅助曲线计算"""
        if len(self.aux_data) < 2:
            messagebox.showwarning("数据不足", "辅助曲线数据至少需要2行")
            return

        try:
            dt = float(self.aux_dt_var.get())
            if dt <= 0:
                raise ValueError("时段间隔必须大于0")
        except ValueError:
            messagebox.showerror("输入错误", "请输入有效的时段间隔（正数）")
            return

        try:
            self.aux_result = FloodCalculator.calc_auxiliary(self.aux_data, dt)
            self._refresh_aux_tree()
            self._update_aux_charts()
            self.status_var.set(
                f"辅助曲线计算完成。共 {len(self.aux_data)} 个数据点，dt={dt}h")
            messagebox.showinfo("计算完成",
                                f"辅助曲线计算完成！\n"
                                f"数据点数：{len(self.aux_data)}\n"
                                f"时段间隔：{dt} h\n"
                                f"请切换到“调洪演算计算”页执行演算。")
        except Exception as e:
            messagebox.showerror("计算错误", f"辅助曲线计算失败:\n{str(e)}")

    def _calculate_routing(self):
        """执行调洪演算"""
        if not self.aux_result:
            messagebox.showwarning("前置条件",
                                    "请先在“调洪辅助曲线计算”页完成辅助曲线计算")
            return

        if len(self.inflow_data) < 2:
            messagebox.showwarning("数据不足", "洪水过程线至少需要2个时段")
            return

        try:
            self.routing_result = FloodCalculator.calc_routing(
                self.inflow_data, self.aux_result)
            self._refresh_inflow_tree()
            self._update_routing_charts()

            # 找最大下泄流量和最高水位
            max_q = max(self.routing_result.outflows)
            max_wl = max(self.routing_result.water_levels)
            max_q_idx = self.routing_result.outflows.index(max_q)

            self.status_var.set(
                f"调洪演算完成。最大下泄流量 {max_q:.2f} m3/s (时段{max_q_idx})，"
                f"最高水位 {max_wl:.3f} m")
            messagebox.showinfo("计算完成",
                                f"调洪演算计算完成！\n\n"
                                f"最大下泄流量：{max_q:.2f} m3/s（时段 {max_q_idx}）\n"
                                f"最高水库水位：{max_wl:.3f} m\n"
                                f"计算时段数：{len(self.routing_result.periods)}")
        except Exception as e:
            messagebox.showerror("计算错误", f"调洪演算失败:\n{str(e)}")

    # ═══════════════════════════════════════════
    # 图表绘制
    # ═══════════════════════════════════════════
    def _update_aux_charts(self):
        """更新辅助曲线图表"""
        if not self.aux_result:
            return

        aux = self.aux_result
        wl = [p.water_level for p in aux.points]
        v = aux.v_above
        q = [p.discharge for p in aux.points]

        # 清图
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()

        # 图1: 水位~库容曲线
        self.ax1.plot(v, wl, 'b-o', markersize=4, linewidth=1.5)
        self.ax1.set_xlabel('Storage V (10^4 m3)')
        self.ax1.set_ylabel('Water Level (m)')
        self.ax1.set_title('Water Level - Storage Curve')
        self.ax1.grid(True, alpha=0.3)

        # 图2: 水位~流量曲线
        self.ax2.plot(q, wl, 'r-o', markersize=4, linewidth=1.5)
        self.ax2.set_xlabel('Discharge q (m3/s)')
        self.ax2.set_ylabel('Water Level (m)')
        self.ax2.set_title('Water Level - Discharge Curve')
        self.ax2.grid(True, alpha=0.3)

        # 图3: 双辅助曲线
        self.ax3.plot(aux.add, q, 'g-o', markersize=4, linewidth=1.5, label='V/dt+q/2')
        self.ax3.plot(aux.sub, q, 'b-s', markersize=4, linewidth=1.5, label='V/dt-q/2')
        self.ax3.set_xlabel('V/dt +/- q/2 (m3/s)')
        self.ax3.set_ylabel('Discharge q (m3/s)')
        self.ax3.set_title('Dual Auxiliary Curves')
        self.ax3.legend(fontsize=8)
        self.ax3.grid(True, alpha=0.3)

        # 中文标题（如果字体可用）
        if _CN_FONT_NAME:
            self.ax1.set_title('水位~库容曲线')
            self.ax1.set_xlabel('库容 (10^4 m3)')
            self.ax1.set_ylabel('水位 (m)')
            self.ax2.set_title('水位~流量曲线')
            self.ax2.set_xlabel('流量 (m3/s)')
            self.ax2.set_ylabel('水位 (m)')
            self.ax3.set_title('双辅助曲线')
            self.ax3.set_xlabel('V/dt +/- q/2 (m3/s)')
            self.ax3.set_ylabel('下泄流量 (m3/s)')

        self.fig_aux.tight_layout()
        self.canvas_top.draw()

    def _update_routing_charts(self):
        """更新演算过程线图表"""
        if not self.routing_result:
            return

        r = self.routing_result
        periods = r.periods

        self.ax4.clear()
        self.ax5.clear()

        # 图4: 下泄流量过程线 + 来水流量过程线
        self.ax4.plot(periods, r.inflows, 'b--o', markersize=3, linewidth=1,
                      alpha=0.6, label='Inflow Q' if not _CN_FONT_NAME else '来水流量Q')
        self.ax4.plot(periods, r.outflows, 'r-o', markersize=4, linewidth=1.5,
                      label='Outflow q' if not _CN_FONT_NAME else '下泄流量q')
        self.ax4.set_xlabel('Time Period' if not _CN_FONT_NAME else '时段')
        self.ax4.set_ylabel('Flow (m3/s)' if not _CN_FONT_NAME else '流量 (m3/s)')
        self.ax4.set_title('Outflow Hydrograph' if not _CN_FONT_NAME else '调洪下泄流量过程线')
        self.ax4.legend(fontsize=8)
        self.ax4.grid(True, alpha=0.3)

        # 图5: 水位过程线
        self.ax5.plot(periods, r.water_levels, 'g-o', markersize=4, linewidth=1.5)
        self.ax5.set_xlabel('Time Period' if not _CN_FONT_NAME else '时段')
        self.ax5.set_ylabel('Water Level (m)' if not _CN_FONT_NAME else '水位 (m)')
        self.ax5.set_title('Water Level Process' if not _CN_FONT_NAME else '洪水位过程线')
        self.ax5.grid(True, alpha=0.3)

        self.fig_routing.tight_layout()
        self.canvas_bottom.draw()

    def _refresh_charts(self):
        """刷新所有图表"""
        self._update_aux_charts()
        self._update_routing_charts()

    def _save_charts(self):
        """保存图表为图片"""
        if not self.aux_result and not self.routing_result:
            messagebox.showwarning("提示", "没有可保存的图表，请先执行计算")
            return

        filepath = filedialog.asksaveasfilename(
            title="保存图表",
            defaultextension=".png",
            filetypes=[("PNG 图片", "*.png"), ("所有文件", "*.*")]
        )
        if not filepath:
            return

        try:
            # 创建完整图表
            fig, axes = plt.subplots(2, 3, figsize=(18, 10))
            fig.suptitle(APP_TITLE, fontsize=16, fontweight='bold')

            if self.aux_result:
                aux = self.aux_result
                wl = [p.water_level for p in aux.points]
                v = aux.v_above
                q = [p.discharge for p in aux.points]

                axes[0, 0].plot(v, wl, 'b-o', markersize=4)
                axes[0, 0].set_title('水位~库容曲线')
                axes[0, 0].set_xlabel('库容 (10^4 m3)')
                axes[0, 0].set_ylabel('水位 (m)')
                axes[0, 0].grid(True, alpha=0.3)

                axes[0, 1].plot(q, wl, 'r-o', markersize=4)
                axes[0, 1].set_title('水位~流量曲线')
                axes[0, 1].set_xlabel('流量 (m3/s)')
                axes[0, 1].set_ylabel('水位 (m)')
                axes[0, 1].grid(True, alpha=0.3)

                axes[0, 2].plot(aux.add, q, 'g-o', markersize=4, label='V/dt+q/2')
                axes[0, 2].plot(aux.sub, q, 'b-s', markersize=4, label='V/dt-q/2')
                axes[0, 2].set_title('双辅助曲线')
                axes[0, 2].set_xlabel('V/dt +/- q/2 (m3/s)')
                axes[0, 2].set_ylabel('下泄流量 (m3/s)')
                axes[0, 2].legend(fontsize=8)
                axes[0, 2].grid(True, alpha=0.3)

            if self.routing_result:
                r = self.routing_result
                axes[1, 0].plot(r.periods, r.inflows, 'b--o', markersize=3, alpha=0.6, label='来水Q')
                axes[1, 0].plot(r.periods, r.outflows, 'r-o', markersize=4, label='下泄q')
                axes[1, 0].set_title('调洪下泄流量过程线')
                axes[1, 0].set_xlabel('时段')
                axes[1, 0].set_ylabel('流量 (m3/s)')
                axes[1, 0].legend(fontsize=8)
                axes[1, 0].grid(True, alpha=0.3)

                axes[1, 1].plot(r.periods, r.water_levels, 'g-o', markersize=4)
                axes[1, 1].set_title('洪水位过程线')
                axes[1, 1].set_xlabel('时段')
                axes[1, 1].set_ylabel('水位 (m)')
                axes[1, 1].grid(True, alpha=0.3)

                axes[1, 2].axis('off')

            fig.tight_layout(rect=[0, 0, 1, 0.96])
            fig.savefig(filepath, dpi=150, bbox_inches='tight')
            plt.close(fig)
            messagebox.showinfo("保存成功", f"图表已保存到:\n{filepath}")
        except Exception as e:
            messagebox.showerror("保存错误", f"保存图表失败:\n{str(e)}")

    # ═══════════════════════════════════════════
    # 导出
    # ═══════════════════════════════════════════
    def _export_results(self):
        """导出计算结果到 Excel"""
        if not self.aux_result:
            messagebox.showwarning("提示", "请先完成辅助曲线计算")
            return
        if not self.routing_result:
            messagebox.showwarning("提示", "请先完成调洪演算计算")
            return

        filepath = filedialog.asksaveasfilename(
            title="导出计算结果",
            defaultextension=".xlsx",
            filetypes=[("Excel 文件", "*.xlsx"), ("所有文件", "*.*")]
        )
        if not filepath:
            return

        try:
            ExcelIO.export_results(filepath, self.aux_result, self.routing_result)
            messagebox.showinfo("导出成功", f"计算结果已导出到:\n{filepath}")
            self.status_var.set(f"结果已导出: {filepath}")
        except Exception as e:
            messagebox.showerror("导出错误", f"导出失败:\n{str(e)}")

    # ═══════════════════════════════════════════
    # 运行
    # ═══════════════════════════════════════════
    def run(self):
        self.root.mainloop()


# ─────────────────────────────────────────────
# 入口
# ─────────────────────────────────────────────
if __name__ == '__main__':
    _apply_ttk_style()
    app = FloodRoutingApp()
    app.run()
