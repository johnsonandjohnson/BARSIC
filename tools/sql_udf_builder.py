#!/usr/bin/env python3
"""
Builds an SQL file from a template file that has 'include' directives for Python/Java inline UDF code.

'Include' directives are expected to be placed inside Python/Java inline code blocks (between $$ .. $$).
Example:

  -- @include python/<ModuleNameToInclude>.py

Usage:
  python tools/sql_udf_builder.py \
      --template templates/ChemCoreApi.tpl.sql \
      --output sql/ChemCoreApi.sql
"""

import argparse
import re
from pathlib import Path
from typing import Iterable

INCLUDE_RE = re.compile(r"^\s*--\s*@include\s+(?P<path>\S+)\s*$")


def process_lines(lines: Iterable[str], base_dir: Path)->Iterable[str]:
    for line in lines:
        m = INCLUDE_RE.match(line)
        if not m:
            yield line
        else:
            rel_path = m.group('path')
            py_path = (base_dir / rel_path).resolve() if not Path(rel_path).is_absolute() else Path(rel_path)
            yield from py_path.open('r', encoding='utf-8')


def build_udf(tpl_path: Path, base_dir: Path, destination: Path):
    ensure_nl = lambda line: line if line.endswith('\n') else f'{line}\n'
    with destination.open('w', encoding='utf-8') as d:
        out_lines = process_lines(tpl_path.open('r', encoding='utf-8'), base_dir)
        d.writelines(map(ensure_nl, out_lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--template', required=True, help='Path to template SQL file')
    ap.add_argument('--output', required=True, help='Path to output compiled SQL')
    args = ap.parse_args()

    tpl_path = Path(args.template)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    build_udf(tpl_path, tpl_path.parent, out_path)
    print(f"Compiled SQL written to: {out_path}")


if __name__ == '__main__':
    main()
