#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
水库调洪演算 CLI 命令行入口
=============================
供 AI / 脚本 / 批处理无界面调用调洪演算。

用法：
  # 使用内置示例数据
  python flood_routing_cli.py demo

  # 从 Excel 输入文件计算（参考 references/excel_format.md）
  python flood_routing_cli.py run --aux 辅助数据.xlsx --inflow 洪水过程线.xlsx --dt 2 \
      --out 结果.xlsx --json 结果.json

  # 从 CSV 输入（A列=水位, B列=流量, C列=库容；D列=时段, E列=dt, F列=来水流量）
  python flood_routing_cli.py run --aux aux.csv --inflow inflow.csv --out 结果.xlsx

退出码：0=成功, 1=参数/数据错误, 2=计算失败
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flood_routing_core import (FloodCalculator, AuxPoint, InflowPoint, ExcelIO,
                                default_aux_points, default_inflow_points)


def _read_aux_file(path, dt_hours):
    """读取辅助曲线数据文件（.xlsx / .csv）"""
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.xlsx', '.xls'):
        points, dt = ExcelIO.read_aux_data(path)
        return points, dt if dt else dt_hours
    # CSV: 每行 水位,流量,库容
    import csv
    points = []
    with open(path, 'r', encoding='utf-8-sig') as f:
        for row in csv.reader(f):
            if not row or len(row) < 3:
                continue
            try:
                points.append(AuxPoint(float(row[0]), float(row[1]), float(row[2])))
            except ValueError:
                continue
    if not points:
        raise ValueError(f"未从 {path} 读取到有效辅助曲线数据")
    return points, dt_hours


def _read_inflow_file(path):
    """读取洪水过程线数据文件（.xlsx / .csv）"""
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.xlsx', '.xls'):
        return ExcelIO.read_inflow_data(path)
    # CSV: 每行 时段序号,时段长dt,来水流量
    import csv
    points = []
    with open(path, 'r', encoding='utf-8-sig') as f:
        for row in csv.reader(f):
            if not row or len(row) < 3:
                continue
            try:
                points.append(InflowPoint(float(row[0]), float(row[1]), float(row[2])))
            except ValueError:
                continue
    if not points:
        raise ValueError(f"未从 {path} 读取到有效洪水过程线数据")
    return points


def _export_excel(filepath, aux, routing):
    """导出结果到 Excel（复用原版导出逻辑，支持中文标题）"""
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment, Border, Side, PatternFill

    wb = Workbook()
    ws1 = wb.active
    ws1.title = "调洪辅助曲线计算表"
    headers1 = ['水位(m)', '下泄流量(m3/s)', '总库容(10^4m3)',
                '溢洪水位以上库容(10^4m3)', 'V/dt(m3/s)', 'q/2(m3/s)',
                'V/dt-q/2(m3/s)', 'V/dt+q/2(m3/s)']
    header_font = Font(bold=True, size=10)
    header_fill = PatternFill(start_color='D9E1F2', end_color='D9E1F2', fill_type='solid')
    thin = Border(left=Side(style='thin'), right=Side(style='thin'),
                  top=Side(style='thin'), bottom=Side(style='thin'))
    ws1.merge_cells('A1:H1')
    ws1['A1'] = f'调洪辅助曲线计算表 (dt={aux.dt_hours}h)'
    ws1['A1'].font = Font(bold=True, size=14)
    ws1['A1'].alignment = Alignment(horizontal='center')
    for col, h in enumerate(headers1, 1):
        c = ws1.cell(row=2, column=col, value=h)
        c.font = header_font; c.fill = header_fill; c.border = thin
        c.alignment = Alignment(horizontal='center', wrap_text=True)
    for i, p in enumerate(aux.points):
        row = i + 3
        data = [p.water_level, p.discharge, p.storage_total, aux.v_above[i],
                aux.v_dt[i], aux.q_half[i], aux.sub[i], aux.add[i]]
        for col, val in enumerate(data, 1):
            c = ws1.cell(row=row, column=col, value=val)
            c.border = thin; c.alignment = Alignment(horizontal='center')
    for col in range(1, 9):
        ws1.column_dimensions[chr(64 + col)].width = 18

    ws2 = wb.create_sheet("调洪演算计算表")
    headers2 = ['时段序号', '时段长dt(h)', '来水流量Q(m3/s)', '平均流量Qp(m3/s)',
                'V/dt-q/2(m3/s)', 'V/dt+q/2(m3/s)', '下泄流量q(m3/s)', '水库水位(m)']
    ws2.merge_cells('A1:H1')
    ws2['A1'] = '调洪演算计算表'
    ws2['A1'].font = Font(bold=True, size=14)
    ws2['A1'].alignment = Alignment(horizontal='center')
    for col, h in enumerate(headers2, 1):
        c = ws2.cell(row=2, column=col, value=h)
        c.font = header_font; c.fill = header_fill; c.border = thin
        c.alignment = Alignment(horizontal='center', wrap_text=True)
    for i in range(len(routing.periods)):
        row = i + 3
        data = [routing.periods[i], routing.dt_list[i], routing.inflows[i],
                routing.avg_inflows[i], routing.sub_prev[i], routing.add_curr[i],
                routing.outflows[i], routing.water_levels[i]]
        for col, val in enumerate(data, 1):
            c = ws2.cell(row=row, column=col, value=val)
            c.border = thin; c.alignment = Alignment(horizontal='center')
    for col in range(1, 9):
        ws2.column_dimensions[chr(64 + col)].width = 18

    # 特征值表
    ws3 = wb.create_sheet("调洪特征值")
    s = FloodCalculator.summary(aux, routing)
    ws3['A1'] = '调洪演算特征值'
    ws3['A1'].font = Font(bold=True, size=14)
    rows = [
        ('起调水位 (m)', s['start_wl']),
        ('最大下泄流量 (m3/s)', s['q_max']),
        ('最大下泄流量时段', s['q_max_period']),
        ('最高水位 (m)', s['wl_max']),
        ('最高水位时段', s['wl_max_period']),
        ('最大来水流量 (m3/s)', s['inflow_max']),
        ('最大调蓄库容 (10^4 m3)', s['v_storage_max']),
        ('计算时段间隔 (h)', s['dt_hours']),
    ]
    for i, (k, v) in enumerate(rows, 2):
        ws3.cell(row=i, column=1, value=k).border = thin
        ws3.cell(row=i, column=2, value=v).border = thin
    ws3.column_dimensions['A'].width = 28
    ws3.column_dimensions['B'].width = 18

    wb.save(filepath)


def cmd_demo(args=None):
    """演示：使用内置示例数据计算"""
    aux = FloodCalculator.calc_auxiliary(default_aux_points(), 2.0)
    routing = FloodCalculator.calc_routing(default_inflow_points(), aux)
    s = FloodCalculator.summary(aux, routing)
    print("=== 水库双辅助曲线调洪演算（示例数据） ===")
    print(f"起调水位        : {s['start_wl']} m")
    print(f"最大来水流量    : {s['inflow_max']} m3/s")
    print(f"最大下泄流量    : {s['q_max']} m3/s（时段 {s['q_max_period']}）")
    print(f"最高水位        : {s['wl_max']} m（时段 {s['wl_max_period']}）")
    print(f"最大调蓄库容    : {s['v_storage_max']} 10^4 m3")
    print(f"计算时段数      : {len(routing.periods)}")
    return 0


def cmd_run(args):
    """从文件输入执行调洪演算"""
    try:
        aux_points, dt = _read_aux_file(args.aux, args.dt)
        inflows = _read_inflow_file(args.inflow)
        if args.dt is not None:
            dt = args.dt
    except Exception as e:
        print(f"[错误] 读取输入失败: {e}", file=sys.stderr)
        return 1

    try:
        aux = FloodCalculator.calc_auxiliary(aux_points, dt)
        routing = FloodCalculator.calc_routing(inflows, aux)
        s = FloodCalculator.summary(aux, routing)
    except Exception as e:
        print(f"[错误] 计算失败: {e}", file=sys.stderr)
        return 2

    # 输出
    if args.out:
        try:
            _export_excel(args.out, aux, routing)
            print(f"结果已导出: {args.out}")
        except Exception as e:
            print(f"[错误] 导出Excel失败: {e}", file=sys.stderr)
            return 2
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump({'summary': s, 'aux': aux.to_dict(),
                       'routing': routing.to_dict()}, f, ensure_ascii=False, indent=2)
        print(f"JSON已导出: {args.json}")

    print("=== 调洪演算结果 ===")
    print(f"起调水位        : {s['start_wl']} m")
    print(f"最大下泄流量    : {s['q_max']} m3/s（时段 {s['q_max_period']}）")
    print(f"最高水位        : {s['wl_max']} m（时段 {s['wl_max_period']}）")
    print(f"最大调蓄库容    : {s['v_storage_max']} 10^4 m3")
    return 0


def main():
    parser = argparse.ArgumentParser(description='水库双辅助曲线调洪演算 CLI')
    sub = parser.add_subparsers(dest='cmd')

    p_demo = sub.add_parser('demo', help='使用内置示例数据计算')
    p_demo.set_defaults(fn=cmd_demo)

    p_run = sub.add_parser('run', help='从文件输入执行调洪演算')
    p_run.add_argument('--aux', required=True, help='辅助曲线数据文件 (.xlsx/.xls/.csv)')
    p_run.add_argument('--inflow', required=True, help='洪水过程线数据文件 (.xlsx/.xls/.csv)')
    p_run.add_argument('--dt', type=float, default=None, help='计算时段间隔 (h)，默认从数据推断或2h')
    p_run.add_argument('--out', default=None, help='导出Excel结果文件路径')
    p_run.add_argument('--json', default=None, help='导出JSON结果文件路径')
    p_run.set_defaults(fn=cmd_run)

    args = parser.parse_args()
    if not hasattr(args, 'fn'):
        parser.print_help()
        return 1
    return args.fn(args)


if __name__ == '__main__':
    sys.exit(main())
