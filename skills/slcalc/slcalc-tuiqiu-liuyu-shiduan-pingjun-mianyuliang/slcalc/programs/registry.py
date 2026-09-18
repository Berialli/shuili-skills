# -*- coding: utf-8 -*-
"""
slcalc 程序注册表
=================
程序编号 → 模块映射。新增程序内核后在此登记。
"""
from . import (a01, a02, a03, a03x, a04, a05x, a06, a07, a08, a09, a10,
               a11, a12, a13, a14, a15, c01, c02, c03, d06, d14a, g6)

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
    "A-4": {
        "name": "水文系列代表性分析程序",
        "module": a04,
        "author": "孙建峰（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现",
    },
    "A-5X": {
        "name": "同频率缩放设计洪水过程线程序",
        "module": a05x,
        "author": "刘晓东（江西省水利规划设计院）",
        "priority": "★★★",
        "status": "已实现（reverse 模式：表2 57 行 ±1 命中 57/57、表3 37 点全中；表1 洪量误差<0.01%）",
    },
    "A-6": {
        "name": "推求流域时段平均面雨量程序",
        "module": a06,
        "author": "陈沂（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现",
    },
    "A-7": {
        "name": "分析经验单位线及汇流计算程序",
        "module": a07,
        "author": "谢熙曦（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（例1 C=1 单位线/SUMU/拟合 Q、例2 C=2 表流 Q 逐值验证 ±0.02 内全 PASS）",
    },
    "A-8": {
        "name": "分析马司京根法演算参数程序",
        "module": a08,
        "author": "郝福良（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（M=1/X/K/XL/KL/C0-C2/S2/O 18点逐值验证全 PASS；X/K 槽蓄 LSQ 反演与原著末位偏差 <0.1% 见模块说明）",
    },
    "A-9": {
        "name": "马司京根法分段连续演算程序",
        "module": a09,
        "author": "郝福良（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（M=2 分段折算 XL=0.5-(M/2)(1-2X)=0.0218/KL=K/M 演算：O 35点中 34点 ±0、35/35 ±1；O14 差 1 源程序舍入悬案）",
    },
    "A-10": {
        "name": "推理公式法计算洪峰流量程序",
        "module": a10,
        "author": "郝福良（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（推理公式法 Sp=HP·24^(N-1)、τ=0.278L/(mJ^(1/3)Q^(1/4))、Tc=((1-N)Sp/U)^(1/N) 判别全面/部分汇流迭代：HP=482→QM=2043、HP=339→QM=1088 两值 ±0 全 PASS）",
    },
    "A-11": {
        "name": "最大24小时洪量计算程序",
        "module": a11,
        "author": "陈沂（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（面积包围法：摘录T为时.分编码+跨日回跳→绝对时刻，梯形累计水量 C 后直线内插整点，逐时洪量=ΔC×0.006 round，24h 滑窗：逐时 35/35 ±0、W24=10142 万立米 ±0 全 PASS）",
    },
    "A-12": {
        "name": "下渗曲线产流计算程序",
        "module": a12,
        "author": "谢熙曦（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（霍顿C=1/菲利普C≠1：F~S点据51点逐点±0.01、产流时段=完整2min格×权重E、SUMR=4.98命中；28值中27/28 ±0，站1 R6 差0.01为VB单精度舍入边界）",
    },
    "A-13": {
        "name": "随机水文AR(P)模型分析程序",
        "module": a13,
        "author": "孙建峰（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（统计参数 CV=n−1无偏/CS=n−3分母、R(k)=分段去均值Pearson相关 逐位命中；正态转换 ln(Q+CC) cc=2.3116 后参数+独立性判定全命中；裁决：P 仅入标题、生成按独立抽取(AR0语义) Box-Muller，组统计公式同原始；确定性输出 68/68 PASS）",
    },
    "A-14": {
        "name": "马斯京根模型最优参数估计程序",
        "module": a14,
        "author": "叶泽纲（湖南省水文总站）",
        "priority": "★★★",
        "status": "已实现（W排队→循环差V/Y估计X=0.101→示蓄φ带截距回归K=18.21→演算 14/14 点全 PASS）",
    },
    "A-15": {
        "name": "由直方图所限定的曲线拟合程序",
        "module": a15,
        "author": "潘东海、朱凤娟（水利部天津勘测设计研究院；算法：潘东海、刘世芬《由直方图限定的曲线拟合》中国民航学院学报 1992(4)）",
        "priority": "★★★",
        "status": "已实现（段内二次+面积守恒消元+加权LSQ解节点：例1 α=0 纯C1 101/101 逐位、面积 79.976；例2 等效权重 4α:β 节点 9/9 全中、采样 100/101、面积 79986.19≈权威 79986.15）",
    },
    "C-1": {
        "name": "美国气象局溃坝洪水预报简化模型程序（SMPDBK）",
        "module": c01,
        "author": "吴媛娇（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（坝址 Eq1 复刻 + 演进纯几何自洽：板桥 QM+0.6%/HW±0.1m/QX断面1~5 ±6%/TT±0.22h/TP±0.22h 全 PASS；TD 公式破译=TT+ratio·24.2·VOL/(Qpi−Qo) 论文例2 命中 7.76h, 板桥 sec1-5 ≤1.1h。远断面 6~7 Xc 缺口 → QX -10~-17% → TD 发散已记录；Qf 精确取法未破译（非局部曼宁），21/186）",
    },
    "C-2": {
        "name": "水库调洪演算的数值解程序",
        "module": c02,
        "author": "张校正（新疆水利厅）",
        "priority": "★★★",
        "status": "已实现",
    },
    "C-3": {
        "name": "两种供水保证的水库径流调节计算程序",
        "module": c03,
        "author": "王建生（水电部天津勘测设计院）",
        "priority": "★★★",
        "status": "已实现（权威验证6组VX/VY全命中：AY/AM/BY/BM逐值核对通过）",
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
        "status": "已实现（例1 数值闭合：KC/KO 命中<1.4%，CX/CN 残差源于重心简化；例2 待后续）",
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
