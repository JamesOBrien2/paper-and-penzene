"""A release's notes, from its CHANGELOG.md section: python cmake/release-notes.py v0.9.0

Fails (exit 1) when the section is missing, so a tag can't be released without its notes.
"""
import re
import sys
from pathlib import Path

version = sys.argv[1].removeprefix("v")
changelog = Path(sys.argv[2] if len(sys.argv) > 2 else Path(__file__).parent.parent / "CHANGELOG.md").read_text(encoding="utf-8")
m = re.search(rf"^## {re.escape(version)}(?=[ \t]|$)[^\n]*\n(.*?)(?=^## |\Z)", changelog, re.M | re.S)
if not m or not m.group(1).strip():
    sys.exit(f"CHANGELOG.md has no notes for {version}: rename \"Unreleased\" to \"{version} (date)\"")
sys.stdout.reconfigure(encoding="utf-8")  # Windows would print in its code page
notes = re.sub(r"[ \t]*<!--.*?-->", "", m.group(1)).strip()
print(notes)
print(f"\nEvery release: [CHANGELOG.md](https://github.com/JamesOBrien2/penzene/blob/v{version}/CHANGELOG.md)")
