# -*- coding: utf-8 -*-
"""
OUT 计算书生成器（outgen）
==========================
复刻原始程序的汉字计算书风格（A-1.OUT / C-2.OUT 等），
同时支持结构化输出（JSON）与 Markdown 计算书。

统一接口：
  gen(program_id, title, sections, results) -> str

其中 sections 为「(小标题, 行列表)」列表，行即纯文本行。
结果由各程序模块以 dict 形式给出，可一并写入 JSON。
"""
import json
import datetime

LINE = "*" * 77
THIN = "-" * 62


def _head(program_id, title):
    pad = max(2, (77 - len(title) - len(program_id) - 3) // 2)
    return (
        LINE + "\n"
        + "*" * 7 + " " * (70 - len(title) // 2) + " " + "\n"
        + f"*******{'':<3}{title} {program_id:<6}*******\n".replace("       ", "       ", 1)
    )


def render_text(program_id, title, sections, footer_note=None):
    """
    生成原始风格的文本计算书。
    sections: [(sub_title, [lines...]), ...]
    """
    out = []
    out.append(LINE)
    out.append(f"*******{title.center(60)}*******".rstrip())
    out.append(LINE)
    for sub, lines in sections:
        out.append(f"({sub})" if sub else "")
        out.extend(lines)
        out.append(THIN)
    if footer_note:
        out.append(footer_note)
    out.append("")
    return "\n".join(out)


def render_markdown(program_id, title, sections, results=None):
    """
    生成 Markdown 计算书（现代化交付格式）。
    """
    lines = [f"# {title}（{program_id}）", ""]
    for sub, rows in sections:
        if sub:
            lines.append(f"## {sub}")
            lines.append("")
        for r in rows:
            if isinstance(r, (list, tuple)):
                lines.append("| " + " | ".join(str(x) for x in r) + " |")
            else:
                lines.append(str(r))
        lines.append("")
    if results is not None:
        lines.append("## 结果汇总（JSON）")
        lines.append("")
        lines.append("```json")
        lines.append(json.dumps(results, ensure_ascii=False, indent=2))
        lines.append("```")
    return "\n".join(lines)


def write_out(path, text):
    """写出计算书文件。"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def write_json(path, data):
    """写出结构化结果 JSON。"""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def stamp():
    """计算书生成时间戳（用于页脚）。"""
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
