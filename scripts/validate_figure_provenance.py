#!/usr/bin/env python3

from __future__ import annotations

import csv
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    provenance_path = repo_root / "docs" / "FIGURE_PROVENANCE.tsv"
    if not provenance_path.exists():
        print(f"Missing: {provenance_path}", file=sys.stderr)
        return 2

    missing: list[str] = []
    checked = 0
    with provenance_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            anchors = (row.get("anchor_tables") or "").strip()
            if not anchors:
                continue
            for rel in [x.strip() for x in anchors.split(";") if x.strip()]:
                checked += 1
                path = repo_root / rel
                if not path.exists():
                    missing.append(rel)
                    continue
                if path.is_file() and path.stat().st_size == 0:
                    missing.append(rel + " (empty)")

    if missing:
        print("Missing anchor tables referenced by docs/FIGURE_PROVENANCE.tsv:", file=sys.stderr)
        for rel in sorted(set(missing)):
            print(f"- {rel}", file=sys.stderr)
        return 1

    print(f"OK: checked {checked} anchor-table references.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

