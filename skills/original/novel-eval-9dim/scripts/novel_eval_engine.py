# -*- coding: utf-8 -*-
"""
《天工开卷》九维小说章节评测引擎（标准版 v1.0）
=================================================
【设计目标】
  1. 9维可扩展维度表：维度定义全部在 DIMENSIONS 配置中，新增第10、11维只需改配置，引擎零改动。
  2. 动态权重：定量维度用 熵权法（数据驱动）+ 变异系数法 双法融合；
     定性维度从权重池随机抽取（模拟评委视角波动）；全部权重归一化后加权。
  3. ≥11轮评测：每轮动态权重不同，输出每轮9维得分与加权总分。
  4. 去极值统计：每轮9维得分做 去最高分+去最低分 后平均，得到该轮"稳健总分"；
     再对全部轮次的稳健总分做 二次统计（去最高最低后的平均）≥150% 才算合格。
  5. 150%口径说明（2026-08-31主理人裁定）：150% = 满分10分 × 1.5 = 15分，
     即作者写到10分(行业顶级)还不够，要在10分上再提高5分（凝练提纯区 10~15分）。
     达标线 = 15分；每轮9维去极值后平均≥15，全部轮次再去极值后平均≥15，两层都过才算合格。

【用法】
  python novel_eval_engine.py <稿件路径> [轮数N] [--json 输出.json]
  默认 N=11 轮。
"""
import argparse
import json
import math
import os
import random
import re
import sys
import statistics
from collections import OrderedDict

BASE = r"D:\WorkBuddy结果文件\2026-08-26-16-42-46"

# ======================================================================
# 一、九维维度表（可扩展：加维度=在此追加一条，引擎自动纳入）
# 每维字段：
#   id       : 唯一标识
#   name     : 维度名
#   type     : quant(定量, 需数据文件) / qual(定性, 需人工或模型打分)
#   source   : quant -> 特征数据json中的键；qual -> 'manual'(人工/模型)
#   direction: 1=越大越好 / -1=越小越好（金手指位置越小越好）
#   window   : 定量维度的检测窗口说明（仅文档用途）
#   baseline : 该维度在样本中的中位数（用于150%基准换算，自动计算，勿手填）
# ======================================================================
DIMENSIONS = [
    {"id": "hook",     "name": "钩子强度",     "type": "quant", "source": "钩子词命中(前500字)",    "direction": 1, "window": "前500字"},
    {"id": "dialog",   "name": "对话密度",     "type": "quant", "source": "对话密度(引号/千字)",    "direction": 1, "window": "全文"},
    {"id": "goal",     "name": "目标清晰度",   "type": "quant", "source": "目标词命中(前3000字)",   "direction": 1, "window": "前3000字"},
    {"id": "goldfinger","name": "金手指速度",  "type": "quant", "source": "金手指位置",            "direction": -1, "window": "第一章前3500字"},
    {"id": "humor",    "name": "幽默调节",     "type": "quant", "source": "幽默词命中(前4000字)",   "direction": 1, "window": "前4000字"},
    {"id": "writing",  "name": "文笔",         "type": "qual",  "source": "manual", "direction": 1},
    {"id": "immersion", "name": "代入感",      "type": "qual",  "source": "manual", "direction": 1},
    {"id": "emotion",  "name": "情感厚度",     "type": "qual",  "source": "manual", "direction": 1},
    {"id": "job_satisfaction", "name": "职业爽点", "type": "qual", "source": "manual", "direction": 1},
]

# 定性维度权重池（动态权重：每轮从池中随机抽取，模拟评委视角波动）
QUAL_WEIGHT_POOL = {
    "writing":  [0.10, 0.12, 0.14, 0.16],
    "immersion":[0.12, 0.14, 0.16, 0.18],
    "emotion":  [0.10, 0.12, 0.14, 0.16],
    "job_satisfaction": [0.12, 0.14, 0.16, 0.18],
}
# 定性维度每轮评分随机波动幅度（模拟多位评委的视角差异，±0.3分）
QUAL_SCORE_JITTER = 0.45
# 定量维度权重（熵权法+变异系数法融合后，再叠加一个轻微随机扰动，实现'动态'）
QUANT_JITTER = 0.02  # ±0.02 权重扰动

# 满分标准（15分制）
# 口径（主理人2026-08-31裁定）：150% = 满分10分 × 1.5 = 15分，
# 即作者写到10分(行业顶级)还不够，要在10分上再提高5分（凝练提纯区 10~15分）。
FULL_SCORE = 10.0          # 行业顶级线（百分位满分）
EXTRA_SCORE = 5.0          # 凝练提纯区上限（10分以上再提高5分）
MAX_SCORE = 15.0           # 满分 = 10 + 5 = 15（150%）
PASS_RATIO = 1.50          # 150% 达标线 = 15分


# ======================================================================
# 二、文本特征提取（与 37本样本统一算法）
# ======================================================================
HOOK_WORDS = ["突然", "竟然", "没想到", "怎么回事", "却", "猛地", "瞬间", "可怕", "诡异",
              "死了", "穿越", "重生", "系统", "金手指", "震惊", "不敢置信", "居然", "秘密",
              "杀", "血", "神秘", "未知", "意外", "不对劲", "奇怪"]
GOAL_WORDS = ["一定要", "必须", "我要", "等着", "报仇", "变强", "回家", "救", "赢", "考",
              "翻盘", "证明", "目标", "决定", "发誓", "闯", "走出", "考上", "学会", "弄明白",
              "查清楚", "治好", "活下去", "成为"]
HUMOR_WORDS = ["笑了", "好笑", "有趣", "闹", "乐了", "吐槽", "开个玩笑", "无语", "尴尬",
               "嘿嘿", "哈哈", "白痴", "傻子", "逗", "整活", "二哈", "离谱", "不是吧", "卧槽",
               "完蛋", "冤种", "沙雕"]
GOLDEN_START = re.compile(r"第[一二三四五六七八九十百0-9０-９]+[章节卷]")
GF_WORDS = ["金手指", "系统", "面板", "传承", "觉醒", "记忆", "碎片", "戒指", "老爷爷",
            "天赋", "血脉", "异能", "功法", "神器", "宝物", "炉鼎", "残印", "司", "任务",
            "修为", "灵气", "特殊能力", "复制", "扫描", "鉴定", "签到"]


def read_text(path):
    """读取文本并归一化换行：Windows CRLF 的 \r 会污染字符偏移计算
       （4000字窗口实际不足4000有效字符），统一转为 \n"""
    with open(path, "rb") as f:
        data = f.read()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        text = data.decode("gb18030", errors="replace")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def count_chars(text):
    return len(re.sub(r"\s", "", text))


def find_goldfinger_pos(text):
    first = None
    for w in GF_WORDS:
        pos = text.find(w)
        if pos != -1 and (first is None or pos < first):
            first = pos
    if first is not None and first > 3500:
        return None
    return first


def calc_features(text):
    """与 analyze_features.py 统一口径"""
    total_chars = count_chars(text)
    head500 = text[:500]
    hook = sum(1 for w in HOOK_WORDS if w in head500)
    gf = find_goldfinger_pos(text)
    quotes = (text.count('"')
              + sum(text.count(q) for q in ["\u201c", "\u201d", "\u2018", "\u2019", "\u300c", "\u300d", "\u300e", "\u300f"]))
    quote_pairs = quotes / 2.0
    dialog_density = round(quote_pairs / (total_chars / 1000.0), 1) if total_chars else 0
    head3000 = text[:3000]
    goal = sum(1 for w in GOAL_WORDS if w in head3000)
    head4000 = text[:4000]
    humor = sum(1 for w in HUMOR_WORDS if w in head4000)
    return {
        "钩子词命中(前500字)": hook,
        "问号数(前500字)": head500.count("？") + head500.count("?"),
        "金手指位置": gf,
        "对话密度(引号/千字)": dialog_density,
        "目标词命中(前3000字)": goal,
        "幽默词命中(前4000字)": humor,
    }


# ======================================================================
# 三、样本分布加载与分位评分
# ======================================================================
def load_sample_distributions():
    # 样本数据定位（发布版兼容）：
    # 1) 优先：脚本同目录的 特征数据_全量.json（SkillHub 发布包布局：SKILL.md+引擎+样本数据同目录）
    # 2) 回退：本地工作区 BASE/对标样本/ 路径（本地技能目录布局）
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "特征数据_全量.json"),
        os.path.join(BASE, "对标样本", "特征数据_全量.json"),
    ]
    sample_path = next((p for p in candidates if os.path.exists(p)), candidates[-1])
    with open(sample_path, encoding="utf-8") as f:
        data = json.load(f)
    samples = data["样本"]
    dists = {}
    for dim in DIMENSIONS:
        if dim["type"] == "quant":
            key = dim["source"]
            vals = [samples[k][key] for k in samples]
            # 金手指位置可能有None（前3500字未检测到），过滤后使用
            if dim["id"] == "goldfinger":
                vals = [v for v in vals if v is not None]
            dists[dim["id"]] = {"vals": vals,
                                "median": statistics.median(vals) if vals else 0}
    gf_vals = [samples[k]["金手指位置"] for k in samples if samples[k]["金手指位置"] is not None]
    dists["goldfinger"]["gf_vals"] = gf_vals
    return dists


def pct_score(v, vals, invert=False, full=FULL_SCORE, extra=EXTRA_SCORE):
    """百分位反映射（15分制）：
        - 0~10分：稿件在样本中的百分位（invert=True 时位置越靠前分越高）
        - 10~15分（凝练提纯区）：超过样本最大值的部分，按超出比例给0~5分
        - 达到或超过样本最大值 = 满分15分（10+5）"""
    n = len(vals)
    if n == 0:
        return full / 2
    vmax = max(vals)
    vmin = min(vals)
    pct = sum(1 for x in vals if x <= v) / n * 100.0
    if invert:
        pct = 100.0 - pct
    base = pct / 100.0 * full  # 0~10分
    # 凝练提纯区：超越样本极值
    #   direction=1（越大越好）：超过样本最大值 max 按超出比例加分
    #   direction=-1（越小越好，如金手指位置）：低于样本最小值 min 按低于比例加分
    extra_pts = 0.0
    if invert:
        if vmin > 0 and v < vmin:
            exceed = (vmin - v) / vmin
            extra_pts = max(0.0, min(extra, exceed * extra))
    else:
        if vmax > 0 and v > vmax:
            exceed = (v - vmax) / vmax
            extra_pts = max(0.0, min(extra, exceed * extra))
    return round(min(MAX_SCORE, base + extra_pts), 2)


def quant_score(feat, dim, dists):
    """定量维度打分：金手指用 invert（越靠前越好），其余越大越好"""
    key = dim["source"]
    v = feat[key]
    vals = dists[dim["id"]]["vals"]
    if v is None:
        return 2.0
    invert = (dim["direction"] == -1)
    return pct_score(v, vals, invert=invert)


# ======================================================================
# 四、动态权重计算
# ======================================================================
def entropy_weights(feat, dists):
    """熵权法（数据驱动）：信息熵越小(差异越大)权重越高"""
    e_list = []
    for dim in DIMENSIONS:
        if dim["type"] != "quant":
            continue
        vals = dists[dim["id"]]["vals"]
        # 用当前稿件值+样本值构造序列，做归一化后算熵
        series = vals + [feat[dim["source"]] if feat[dim["source"]] is not None else 0]
        series = [float(x) for x in series]
        lo, hi = min(series), max(series)
        norm = [(x - lo) / (hi - lo + 1e-9) for x in series]
        s = sum(norm) + 1e-9
        p = [x / s for x in norm]
        e = -sum(pi * math.log(pi + 1e-12) for pi in p) / math.log(len(p))
        e_list.append((dim["id"], e))
    # 熵转权重：差异大(熵小)权重大
    inv = [1 - e for _, e in e_list]
    s_inv = sum(inv) + 1e-9
    return {did: (1 - e) / s_inv for did, e in e_list}


def coefficient_variation_weights(feat, dists):
    """变异系数法：标准差/均值 越大权重越高（区分度大）"""
    cvs = []
    for dim in DIMENSIONS:
        if dim["type"] != "quant":
            continue
        vals = dists[dim["id"]]["vals"]
        m = statistics.mean(vals)
        sd = statistics.stdev(vals) if len(vals) > 1 else 0
        cv = sd / m if m else 0
        cvs.append((dim["id"], cv))
    s_cv = sum(cv for _, cv in cvs) + 1e-9
    return {did: cv / s_cv for did, cv in cvs}


def dynamic_weights(feat, dists, rng):
    """动态权重：熵权×0.5 + 变异系数×0.5 = 定量基础权重，叠加±0.02扰动；
       定性权重从池中随机抽取；全部归一化"""
    ent = entropy_weights(feat, dists)
    cvw = coefficient_variation_weights(feat, dists)
    w = {}
    quant_ids = [d["id"] for d in DIMENSIONS if d["type"] == "quant"]
    qual_ids = [d["id"] for d in DIMENSIONS if d["type"] == "qual"]
    for did in quant_ids:
        base = 0.5 * ent.get(did, 0) + 0.5 * cvw.get(did, 0)
        w[did] = base + rng.uniform(-QUANT_JITTER, QUANT_JITTER)
    for did in qual_ids:
        pool = QUAL_WEIGHT_POOL.get(did, [0.12])
        w[did] = rng.choice(pool)
    # 归一化
    s = sum(w.values()) + 1e-9
    w = {k: v / s for k, v in w.items()}
    return w


# ======================================================================
# 五、一轮评测
# ======================================================================
def one_round(feat, dists, qual_scores, rng, round_no):
    """一轮评测：9维打分 + 动态加权 + 去极值稳健总分"""
    dim_scores = {}
    for dim in DIMENSIONS:
        if dim["type"] == "quant":
            dim_scores[dim["id"]] = quant_score(feat, dim, dists)
        else:
            base = qual_scores.get(dim["id"], 5.0) * PASS_RATIO  # 10分基准×1.5 → 15分制
            jitter = rng.uniform(-QUAL_SCORE_JITTER, QUAL_SCORE_JITTER)
            dim_scores[dim["id"]] = round(max(1.0, min(MAX_SCORE, base + jitter)), 2)
    w = dynamic_weights(feat, dists, rng)
    weighted = sum(dim_scores[did] * w[did] for did in dim_scores)
    # 去极值：9维中去掉最高分和最低分，剩余7维平均
    sorted_scores = sorted(dim_scores.values())
    trimmed = sorted_scores[1:-1]
    robust_mean = statistics.mean(trimmed)
    # 150%判定（单轮）：稳健总分 >= 15分（满分10×1.5=15，凝练提纯区）
    pass_round = robust_mean >= MAX_SCORE  # 15分制：150%=15分
    return {
        "round": round_no,
        "dim_scores": dim_scores,
        "weights": {k: round(v, 4) for k, v in w.items()},
        "weighted_total": round(weighted, 3),
        "robust_mean": round(robust_mean, 3),
        "pass_150": pass_round,
    }


# ======================================================================
# 六、多轮评测 + 双层去极值统计
# ======================================================================
def run_eval(path, rounds=11, seed=42, qual_override=None):
    text = read_text(path)
    feat = calc_features(text)
    dists = load_sample_distributions()
    rng = random.Random(seed)

    # 定性评分：可用外部JSON覆盖（与旧脚本兼容），否则用默认
    qual_scores = {
        "writing": 9.5, "immersion": 9.0, "emotion": 9.0, "job_satisfaction": 9.2,
    }
    if qual_override:
        qual_scores.update(qual_override)

    rounds_out = []
    for i in range(1, rounds + 1):
        rounds_out.append(one_round(feat, dists, qual_scores, rng, i))

    # 双层去极值：
    # 第一层：每轮9维去最高最低 -> robust_mean
    # 第二层：全部轮次 robust_mean 去最高最低 -> final_mean
    robust_means = [r["robust_mean"] for r in rounds_out]
    sorted_rm = sorted(robust_means)
    trimmed_rm = sorted_rm[1:-1] if len(sorted_rm) >= 3 else sorted_rm
    final_mean = statistics.mean(trimmed_rm)
    # 二次统计：全部稳健平均分的平均（含所有轮次）
    all_mean = statistics.mean(robust_means)
    pass_final = final_mean >= MAX_SCORE  # 15分制：150%=15分

    result = {
        "path": path,
        "rounds": rounds,
        "features": feat,
        "round_details": rounds_out,
        "robust_means": robust_means,
        "trimmed_robust_means": trimmed_rm,
        "final_mean_after_trim": round(final_mean, 3),
        "overall_mean": round(all_mean, 3),
        "pass_150_final": pass_final,
        "standard": "9维可扩展评价标准 v1.0",
        "note": "150%口径（2026-08-31主理人裁定）：150% = 满分10分×1.5 = 15分，即10分(行业顶级)再提高5分(凝练提纯区)。每轮9维去极值后平均≥15分，全部轮次再去极值后平均≥15分即合格。",
    }
    return result


# ======================================================================
# 七、输出
# ======================================================================
def print_result(res):
    print("=" * 70)
    print("九维评测引擎 v1.0  轮数: %d" % res["rounds"])
    print("稿件: %s" % res["path"])
    print("=" * 70)
    print("\n[定量特征]")
    for k, v in res["features"].items():
        print("  %s: %s" % (k, v))
    print("\n[每轮评测] 稳健总分=9维去最高最低后平均 | 达标线15.0")
    for r in res["round_details"]:
        dims = " ".join("%s=%.1f" % (d, s) for d, s in r["dim_scores"].items())
        flag = "✅" if r["pass_150"] else "❌"
        print("  第%2d轮 加权%.3f 稳健%.3f %s  %s" % (r["round"], r["weighted_total"], r["robust_mean"], flag, dims))
    print("\n[双层去极值统计]")
    print("  各轮稳健总分: %s" % res["robust_means"])
    print("  去最高最低后: %s" % res["trimmed_robust_means"])
    print("  全部轮次平均(二次统计): %.3f" % res["overall_mean"])
    print("  去极值后平均: %.3f" % res["final_mean_after_trim"])
    verdict = "✅ 合格" if res["pass_150_final"] else "❌ 不合格"
    print(f"  150%判定(≥15.0): {verdict}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="稿件路径")
    ap.add_argument("rounds", nargs="?", type=int, default=11, help="评测轮数(默认11)")
    ap.add_argument("--seed", type=int, default=42, help="随机种子")
    ap.add_argument("--json", help="输出JSON路径")
    args = ap.parse_args()

    # 尝试加载同目录定性评分覆盖（兼容旧脚本约定）
    qual_override = None
    qpath = os.path.join(os.path.dirname(args.path), "定性评分_v9.json")
    if os.path.exists(qpath):
        with open(qpath, encoding="utf-8") as f:
            qual_override = json.load(f)
        print("[已加载定性评分覆盖: %s]" % qpath)

    res = run_eval(args.path, args.rounds, args.seed, qual_override)
    print_result(res)
    if args.json:
        # JSON中去除不可序列化项
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False, indent=2)
        print("\n已保存: %s" % args.json)


if __name__ == "__main__":
    main()
