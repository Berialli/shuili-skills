# -*- coding: utf-8 -*-
"""
slcalc 命令行入口
=================
用法：
  python -m slcalc list                         # 列出已实现程序
  python -m slcalc A-1 <data.INT> [--out out.txt] [--json out.json] [--md]
  python -m slcalc A-1 --input '{"xs":[...],"ys":[...]}' [--interp "502,610"]
"""
import argparse
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from slcalc.programs.registry import list_programs, get  # noqa: E402


def main():
    ap = argparse.ArgumentParser(prog="slcalc", description="《水利程序集》Python 计算内核")
    ap.add_argument("program", nargs="?", help="程序编号，如 A-1")
    ap.add_argument("data", nargs="?", help="INT 数据文件路径")
    ap.add_argument("--input", help="JSON 格式输入（与 data 二选一）")
    ap.add_argument("--interp", help="插补 X 值列表（逗号分隔，仅部分程序）")
    ap.add_argument("--out", help="输出计算书文本文件")
    ap.add_argument("--json", dest="json_out", help="输出结果 JSON 文件")
    ap.add_argument("--md", action="store_true", help="输出 Markdown 计算书")
    ap.add_argument("--cc", type=float, help="A-13 正态转换系数 CC（权威例 2.3116）")
    args = ap.parse_args()

    if not args.program:
        ap.print_help()
        print("\n已实现程序：")
        for k, v in list_programs().items():
            print(f"  {k:<6} {v['name']}  [{v['priority']}]")
        return

    pid = args.program.upper()
    mod = get(pid)

    if args.input:
        data = json.loads(args.input)
        if args.interp and isinstance(data, dict):
            data["x_interp"] = [float(v) for v in args.interp.split(",")]
    elif args.data:
        data = args.data
    else:
        print("错误：需要提供 --input JSON 或 data INT 文件路径")
        sys.exit(1)

    fmt = "markdown" if args.md else "text"
    kw = {}
    if args.cc is not None and pid == "A-13":
        kw["cc"] = args.cc
    result, text = mod.run(data, out_txt=args.out, out_json=args.json_out,
                           fmt=fmt, **kw)

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
