# -*- coding: utf-8 -*-
"""
INT 数据文件解析器（intio）
===========================
《水利程序集》的数据文件为纯文本、逗号分隔、固定参数顺序，
部分文件首行有中文标题（如「柯柯亚水库,百年一遇,5」）。

设计：
  - parse_int_file(path)  -> 原始行列表 + 数值流（兼容首行标题）
  - read_numbers(path)    -> 按逗号/换行/空白切分后的 float 列表
  - 每个程序模块再按自己的固定参数顺序消费数值流
"""
import re

_NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")


def smart_read_text(path):
    """
    智能编码读取文本文件（依次尝试 utf-8 / gbk / gb18030 / latin-1）。
    兼容原始 INT 文件的 GBK/ANSI 编码与 BOM。
    """
    with open(path, "rb") as f:
        raw = f.read()
    for enc in ("utf-8-sig", "utf-8", "gbk", "gb18030", "latin-1"):
        try:
            return raw.decode(enc)
        except (UnicodeDecodeError, UnicodeError):
            continue
    return raw.decode("latin-1", errors="replace")


def read_lines(path, encoding="gbk", errors="replace"):
    """读取 INT 文件全部行（容忍 GBK/ANSI 编码与 BOM）。"""
    text = smart_read_text(path)
    # 去除 BOM 与文件结束符
    text = text.lstrip("\ufeff").replace("\x1a", "")
    return text.splitlines()


def read_numbers(path, encoding="gbk", errors="replace"):
    """
    将 INT 文件读为数值流（float 列表）。
    自动跳过非数值文本行（如首行中文标题），
    以逗号、换行、空白为分隔符。
    """
    nums = []
    for line in read_lines(path, encoding, errors):
        # 先按逗号拆
        parts = line.split(",")
        for p in parts:
            p = p.strip()
            if not p:
                continue
            m = _NUM_RE.fullmatch(p)
            if m:
                nums.append(float(m.group(0)))
    return nums


def parse_with_header(path, encoding="gbk", errors="replace"):
    """
    解析带可选首行中文标题的 INT 文件。
    返回 (header_str_or_None, nums)
      header: 首行中第一个非数值字段（如「柯柯亚水库,百年一遇」），无则 None
      nums  : 后续数值流
    """
    lines = read_lines(path, encoding, errors)
    header = None
    nums = []
    for i, line in enumerate(lines):
        if i == 0 and not _NUM_RE.fullmatch(line.split(",")[0].strip()):
            header = line.strip()
            continue
        for p in line.split(","):
            p = p.strip()
            if not p:
                continue
            m = _NUM_RE.fullmatch(p)
            if m:
                nums.append(float(m.group(0)))
    return header, nums
