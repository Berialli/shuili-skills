# -*- coding: utf-8 -*-
"""
《水利程序集》技能挂载脚本 —— scripts/calc.py
================================================
供各技能目录调用的统一计算入口。每个技能目录下的 calc.py 会
import 本文件并转发到 slcalc 内核。

用法（在技能目录内）：
  python calc.py <INT文件> [--json out.json] [--md]
  python calc.py --input '{"...": ...}' [--json out.json] [--md]

原理：
  1. 自动定位 slcalc 内核（优先找本技能目录旁的 slcalc/，再找环境变量
     SLCALC_HOME，最后找 ~/.workbuddy/skills/slcalc/）
  2. 按技能目录名（如 G-6）映射程序编号
  3. 调用 slcalc 内核计算并输出
"""
import argparse
import json
import os
import sys

# ---- 定位 slcalc 内核 ----
def find_slcalc():
    candidates = []
    here = os.path.dirname(os.path.abspath(__file__))
    # 0. 技能目录旁的 slcalc（包内自带内核，最高优先：压缩包开箱即用）
    candidates.append(os.path.join(here, "slcalc"))
    # 1. 环境变量
    if os.environ.get("SLCALC_HOME"):
        candidates.append(os.path.join(os.environ["SLCALC_HOME"], "slcalc"))
    # 2. 技能包根目录旁的 slcalc
    skill_root = os.path.dirname(here)
    candidates.append(os.path.join(skill_root, "slcalc"))
    # 3. 上一级（技能包根目录）
    candidates.append(os.path.join(os.path.dirname(skill_root), "slcalc"))
    # 4. 用户技能目录
    candidates.append(os.path.expanduser("~/.workbuddy/skills/slcalc"))
    # 5. 工作目录
    candidates.append(os.path.join(os.getcwd(), "slcalc"))
    for c in candidates:
        if os.path.isdir(os.path.join(c, "programs")):
            return c
    return None

SLCALC = find_slcalc()
if SLCALC:
    sys.path.insert(0, os.path.dirname(SLCALC))
else:
    print("错误: 未找到 slcalc 内核。请设置 SLCALC_HOME 或安装技能包。", file=sys.stderr)
    sys.exit(2)


def get_program_id(skill_dir_name=None):
    """
    由 SKILL.md 的 name 字段或技能目录名推断程序编号。
    SKILL.md frontmatter 的 name 形如 slcalc-a-1 → A-1。
    """
    if not skill_dir_name:
        skill_dir = os.path.dirname(os.path.abspath(__file__))
        # 1. 尝试读同目录 SKILL.md 的 name
        skill_md = os.path.join(skill_dir, "SKILL.md")
        if os.path.isfile(skill_md):
            try:
                with open(skill_md, encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("name:"):
                            name = line.split(":", 1)[1].strip().strip('"\'')
                            # slcalc-a-1 → 提取尾部编号 A-1
                            for seg in name.split("-"):
                                if seg.upper() in ("A", "C", "D", "E", "F", "G",
                                                   "H", "I", "J", "K", "L", "M",
                                                   "N", "P", "Q"):
                                    return f"{seg.upper()}-{name.split('-')[-1].upper()}"
                            break
            except OSError:
                pass
        name = os.path.basename(skill_dir)
    else:
        name = skill_dir_name
    name = name.upper()
    # 取最后一个下划线后的部分（如 TMP_SKILL_G-6 → G-6）
    if "_" in name:
        name = name.rsplit("_", 1)[-1]
    return name


def main():
    ap = argparse.ArgumentParser(description="《水利程序集》技能计算入口")
    ap.add_argument("data", nargs="?", help="数据文件路径")
    ap.add_argument("--input", help="JSON 格式输入")
    ap.add_argument("--json", dest="json_out", help="输出结果 JSON 文件")
    ap.add_argument("--md", action="store_true", help="输出 Markdown 计算书")
    ap.add_argument("--out", help="输出文本计算书文件")
    args = ap.parse_args()

    from slcalc.programs.registry import get

    pid = get_program_id()
    try:
        mod = get(pid)
    except KeyError as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)

    if args.input:
        data = json.loads(args.input)
    elif args.data:
        data = args.data
    else:
        print("错误: 需要提供 数据文件路径或 --input JSON", file=sys.stderr)
        sys.exit(1)

    fmt = "markdown" if args.md else "text"
    result, text = mod.run(data, out_txt=args.out, out_json=args.json_out, fmt=fmt)

    if args.out or args.json_out:
        print(f"✔ 计算完成：{pid}")
        if args.out:
            print(f"  计算书 → {args.out}")
        if args.json_out:
            print(f"  结果   → {args.json_out}")
    else:
        print(text)


if __name__ == "__main__":
    main()
