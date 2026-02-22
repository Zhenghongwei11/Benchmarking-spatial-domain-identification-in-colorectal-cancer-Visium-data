import re
import csv
from pathlib import Path

# 1. Extract citation order from manuscript
manuscript_path = Path("docs/manuscript/manuscript.md")
content = manuscript_path.read_text()
# Find all citations like [@key1; @key2]
groups = re.findall(r"\[@([^]]+)\]", content)
order = []
seen = set()
for group in groups:
    # Split by semicolon and strip whitespace and @
    keys = group.split(";")
    for k in keys:
        clean_key = k.strip().lstrip("@")
        if clean_key and clean_key not in seen:
            order.append(clean_key)
            seen.add(clean_key)

# 2. Read existing references
csv_path = Path("docs/references/doi_list.csv")
refs = {}
with open(csv_path, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        refs[row["id"]] = row

# 3. Reorder according to first appearance
ordered_refs = []
for key in order:
    if key in refs:
        ordered_refs.append(refs.pop(key))
    else:
        print(f"Warning: Citation key '{key}' in text but not in CSV.")

# Add any remaining refs that weren't cited
for key, row in refs.items():
    ordered_refs.append(row)
    print(f"Warning: Reference '{key}' in CSV but not cited in text.")

# 4. Write back to CSV
fieldnames = ["id", "doi", "notes"]
with open(csv_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(ordered_refs)

print(f"Reordered {len(ordered_refs)} references in {csv_path}.")
