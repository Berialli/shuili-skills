---
name: reservoir-flood-routing
slug: reservoir-flood-routing
displayName: 水库双辅助曲线调洪演算
description: 水库调洪演算（双辅助曲线法/半图解法）专业计算技能。当用户需要进行水库调洪演算、推求下泄流量与库水位过程、计算调洪库容、生成调洪辅助曲线（V/dt±q/2）、或处理水位-流量-库容关系与洪水过程线的组合计算时使用。提供计算引擎（Python）、CLI 批量入口、Excel 导入导出与特征值提取。算法源自 Excel VBA 半图解法，与原版程序逐点一致，适用于无闸门自由泄流水库的洪水调节计算。
version: 1.0.0
license: MIT
agent_created: true
read_when:
  - 用户要求进行水库调洪演算、洪水调节计算
  - 用户需要推求调洪下泄流量过程、库水位过程、最高洪水位
  - 用户需要计算调洪库容、校核设计洪水工况
  - 用户要求生成调洪辅助曲线（双辅助曲线法 V/dt±q/2）
  - 用户提供水位~库容~泄流能力关系与入库洪水过程线，要求联解水量平衡
  - 用户要求复现或验证既有 Excel VBA 调洪演算成果
---

# 水库双辅助曲线调洪演算

基于**双辅助曲线法（半图解法）**的水库调洪演算技能，算法与原 Excel VBA 调洪演算程序完全一致（逐点验证通过）。

## 何时使用

- 用户要求做水库调洪演算、推求下泄流量和库水位过程
- 需要最高洪水位、最大下泄流量、调洪库容等特征值
- 需要生成/校核调洪辅助曲线，或复现 Excel VBA 调洪成果
- 结合 SL 252-2017、SL 44 等规范出具调洪演算结论

## 运行环境

```bash
# 依赖（只需核心引擎时仅需 numpy 可选；Excel 读写需要 openpyxl / xlrd）
pip install openpyxl xlrd
# 有界面环境可选安装 tkinter（Windows 官方 Python 自带）
```

## 脚本清单

| 文件 | 作用 |
|---|---|
| `scripts/flood_routing_core.py` | **核心计算引擎**（无 GUI 依赖）：辅助曲线、演算、特征值、Excel 导入导出 |
| `scripts/flood_routing_cli.py` | CLI 入口：demo / run（Excel/CSV 输入 → Excel/JSON 输出） |
| `scripts/flood_routing_gui.py` | 原版 GUI 程序（tkinter + matplotlib，已修复原版嵌套引号语法错误） |
| `scripts/test_flood_routing.py` | 337 项自动化测试（含与原版逐点一致性对比） |
| `assets/示例数据_辅助曲线.xlsx` | 示例辅助曲线数据（水位~流量~库容，15 点） |
| `assets/示例数据_洪水过程线.xlsx` | 示例洪水过程线数据（48 时段） |
| `references/algorithm.md` | 算法原理、单位换算、数据格式、适用范围 |

## 使用方法

### 方式一：Python 调用（推荐，AI 批量计算）

```python
import sys
sys.path.insert(0, 'scripts')  # 按实际路径调整
from flood_routing_core import FloodCalculator

# 一站式接口：输入元组列表，返回 dict
result = FloodCalculator.calc_flood(
    aux_points=[(285.0, 0.0, 3504.0), (285.5, 38.1, 3624.0), ...],  # (水位, 流量, 库容)
    inflows=[(0, 2, 0.0), (1, 2, 116.0), ...],                       # (时段, dt, 来水)
    dt_hours=2.0,
)
summary = result['summary']
print(summary['q_max'], summary['wl_max'], summary['v_storage_max'])

# 分步调用
aux = FloodCalculator.calc_auxiliary(aux_points, dt_hours=2.0)
routing = FloodCalculator.calc_routing(inflows, aux)
s = FloodCalculator.summary(aux, routing)

# 内置示例数据
from flood_routing_core import default_aux_points, default_inflow_points
```

`calc_flood` 返回结构：
```json
{
  "aux": {"points": [...], "v_above": [...], "v_dt": [...], "q_half": [...], "sub": [...], "add": [...], "dt_hours": 2.0},
  "routing": {"periods": [...], "inflows": [...], "avg_inflows": [...], "sub_prev": [...], "add_curr": [...], "outflows": [...], "water_levels": [...]},
  "summary": {"q_max": ..., "q_max_period": ..., "wl_max": ..., "wl_max_period": ..., "inflow_max": ..., "v_storage_max": ..., "dt_hours": ..., "start_wl": ...}
}
```

### 方式二：CLI 命令行（脚本/批处理）

```bash
# 演示（内置示例数据）
python scripts/flood_routing_cli.py demo

# 从 Excel/CSV 输入计算，导出 Excel + JSON
python scripts/flood_routing_cli.py run \
    --aux 辅助曲线.xlsx --inflow 洪水过程线.xlsx --dt 2 \
    --out 结果.xlsx --json 结果.json
```

CSV 格式：辅助曲线 `水位,流量,库容`；洪水过程线 `时段序号,dt,来水流量`。
Excel 格式：见 `references/algorithm.md` §6（自动识别"辅助"/"演算"工作表）。

### 方式三：GUI 桌面程序（原版）

```bash
python scripts/flood_routing_gui.py
```

四页签：调洪辅助曲线计算 / 调洪演算计算 / 图表展示 / 程序说明。支持手动输入、批量粘贴、Excel 导入导出、图表保存。**原版存在嵌套 ASCII 引号语法错误，本技能版本已修复**。

## 计算流程规范

1. **收集数据**：水位~下泄流量~总库容关系（水位单调递增，第一点为起调水位）+ 入库洪水过程线（第 0 时段为初始状态）
2. **确定 dt**：常规 1~2 h；需与洪水过程线时段一致
3. **计算辅助曲线**：`calc_auxiliary` → V/dt±q/2 双辅助曲线
4. **执行演算**：`calc_routing` → 逐时段下泄流量与水位
5. **提取特征值**：`summary` → q_max、H_max、调洪库容、发生时段
6. **成果输出**：Excel（辅助曲线表 + 演算表 + 特征值表）或 JSON
7. **规范性引用**（如出报告）：SL 252-2017《水利水电工程等级划分及洪水标准》、SL 44《水利水电工程设计洪水计算规范》

## 质量保障

- 运行 `python scripts/test_flood_routing.py` 验证（337 项用例，含与原版逐点一致性）
- 设置环境变量 `LEGACY_FLOOD_ROUTING=原版flood_routing.py路径` 可启用与原版逐点对比
- 复核要点：起调水位、最大下泄流量发生时段、最高水位是否超过溢洪道堰顶+安全超高

## 注意事项

1. 辅助曲线数据水位必须**单调递增**；覆盖演算水位区间
2. 超出辅助曲线范围时取端点外推（结果仅供参考，建议加密数据点）
3. 本方法为半图解法，适用于无闸门自由泄流；闸门调度需改用分级泄流曲线
4. 单位：水位 m、流量 m³/s、库容 10⁴ m³、时段 h；换算系数 0.36
5. 计算精度：流量 2 位小数、水位 3 位小数（与原程序一致）


---

## 知识产权声明

本技能为**哈胜**独立原创成果。**本技能不是对任何第三方软件的改编或衍生。**

版权所有 © 哈胜，保留所有权利；未经书面授权，禁止复制、改编、再分发、转售，
或用于训练同类产品。

> 技能中引用的国家标准、行业规范、技术手册及政府公开文件，仅作为计算依据引用，
> 其著作权归各自发布机构所有。
