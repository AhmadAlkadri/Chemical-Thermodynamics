from __future__ import annotations

from typing import Any, Iterable

import numpy as np


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float, np.number)) and not isinstance(value, bool)


def _format_number(value: Any) -> str:
    return f"{float(value):.6g}"


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if _is_number(value):
        return _format_number(value)
    if isinstance(value, np.ndarray) and value.shape == ():
        return _fmt(value.item())
    if isinstance(value, (list, tuple, np.ndarray)):
        items = list(value)
        if all(_is_number(item) for item in items):
            return "[" + ", ".join(_format_number(item) for item in items) + "]"
        return "[" + ", ".join(_fmt(item) for item in items) + "]"
    return str(value)


def print_table(
    rows: Iterable[dict[str, Any]],
    columns: list[str] | None = None,
    title: str | None = None,
) -> None:
    rows = list(rows)
    if not rows:
        return
    if columns is None:
        columns = list(rows[0].keys())
    if title:
        print(title)
    formatted = [[_fmt(row.get(col)) for col in columns] for row in rows]
    widths = [max(len(col), max(len(row[i]) for row in formatted)) for i, col in enumerate(columns)]
    header = "  ".join(col.ljust(widths[i]) for i, col in enumerate(columns))
    print(header)
    print("  ".join("-" * widths[i] for i in range(len(columns))))
    for row in formatted:
        print("  ".join(row[i].ljust(widths[i]) for i in range(len(columns))))
    print()
