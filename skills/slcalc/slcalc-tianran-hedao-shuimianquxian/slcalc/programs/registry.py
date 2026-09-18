# -*- coding: utf-8 -*-
"""
slcalc 程序注册表
=================
程序编号 → 模块映射。新增程序内核后在此登记。
"""
from . import a01, a02, a03, a03x, c02, d06, d14a, g6

REGISTRY = {
    "A-1": {
        "name": "直线相关计算程序",
        "module": a01,
        "author": "陈丽棠（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现",
    },
    "A-2": {
        "name": "三参数幂函数曲线拟合程序",
        "module": a02,
        "author": "陈沂（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现",
    },
    "A-3": {
        "name": "水文频率计算程序",
        "module": a03,
        "author": "马明（新疆水利水电勘测设计院）",
        "priority": "★★★",
        "status": "已实现",
    },
    "A-3X": {
        "name": "水文频率计算程序（连续/不连续系列）",
        "module": a03x,
        "author": "马明（新疆水利水电勘测设计院）",
        "priority": "★★★",
        "status": "已实现",
    },
    "C-2": {
        "name": "水库调洪演算的数值解程序",
        "module": c02,
        "author": "张校正（新疆水利厅）",
        "priority": "★★★",
        "status": "已实现",
    },
    "D-6": {
        "name": "常用断面渠道水力学计算程序",
        "module": d06,
        "author": "杨志河（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现",
    },
    "D-14A": {
        "name": "推求法计算天然河道水面曲线程序",
        "module": d14a,
        "author": "张校正（新疆水利厅）",
        "priority": "★★★",
        "status": "已实现",
    },
    "G-6": {
        "name": "重力式挡土墙计算程序",
        "module": g6,
        "author": "刘以钢（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（例1/例2 有偏差，见分析报告）",
    },
}


def list_programs():
    return {k: {"name": v["name"], "status": v["status"], "priority": v["priority"]}
            for k, v in REGISTRY.items()}


def get(pid):
    pid = pid.upper()
    if pid not in REGISTRY:
        raise KeyError(f"程序 {pid} 尚未注册。可用：{', '.join(REGISTRY)}")
    return REGISTRY[pid]["module"]
