#!/usr/bin/env python
# Create soft links in ../.venv/bin to all files in ./sys.argv[1]
# execute as:
#   > create_soft_links_in_venv_bin source
# where `source` is a subfolder of `.venv_repo`

import sys
from pathlib import Path

if __name__ == "__main__":
    assert len(sys.argv) == 2, "expected one argument (Path)"
    parent = Path(__file__).parent
    selection = parent / sys.argv[1]
    assert selection.is_dir(), f"Directory {selection} does not exist."
    destination = (parent.parent / ".venv").resolve()
    for p in selection.glob('**/*'):
        p_rel = p.relative_to(selection)
        dp = destination / p_rel
        if p.is_dir():
            print(f"{selection.name}/{p_rel} : creating folder `{dp}`")
            dp.mkdir(parents=False, exist_ok=True)
        elif p.is_file():
            (destination / p_rel).symlink_to(p.resolve())
            print(f"{selection.name}/{p_rel} : creating symlink `{dp}`")

