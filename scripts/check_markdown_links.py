from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import unquote, urlparse

REPO_ROOT = Path(__file__).resolve().parents[1]
IGNORED_DIRS = {
    ".git",
    ".pytest_cache",
    ".ruff_cache",
    ".venv",
    "build",
    "dist",
    "test-data",
    "tmp",
}
LINK_RE = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")


def _is_ignored(path: Path) -> bool:
    return any(part in IGNORED_DIRS for part in path.relative_to(REPO_ROOT).parts)


def _iter_markdown_lines(path: Path):
    in_code_block = False
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.lstrip().startswith("```"):
            in_code_block = not in_code_block
            continue
        if not in_code_block:
            yield line_number, line


def _target_exists(path: Path, raw_target: str) -> bool:
    target = raw_target.strip()
    if not target or target.startswith("#"):
        return True

    parsed = urlparse(target)
    if parsed.scheme or parsed.netloc:
        return True

    target_path = unquote(target.split("#", 1)[0].strip())
    if not target_path:
        return True
    if target_path.startswith("<") and target_path.endswith(">"):
        target_path = target_path[1:-1]

    return (path.parent / target_path).resolve().exists()


def main() -> None:
    failures: list[str] = []
    for path in sorted(REPO_ROOT.rglob("*.md")):
        if _is_ignored(path):
            continue
        for line_number, line in _iter_markdown_lines(path):
            for match in LINK_RE.finditer(line):
                target = match.group(1)
                if not _target_exists(path, target):
                    failures.append(f"{path.relative_to(REPO_ROOT)}:{line_number}: {target}")

    if failures:
        print("Broken Markdown links:")
        print("\n".join(failures))
        raise SystemExit(1)
    print("Markdown repo-relative links OK")


if __name__ == "__main__":
    main()
