#!/usr/bin/env python
# Create soft links in ../.venv/bin/pythonX.Y/site-packages to all MOCCaPy/mocca
# execute as:
#   > dev_install_mocca.py

from pathlib import Path

if __name__ == "__main__":
    p_MOCCaPy = Path(__file__).parent
    while p_MOCCaPy.name != 'MOCCaPy':
        p_MOCCaPy = p_MOCCaPy.parent
    p_MOCCaPy_mocca = p_MOCCaPy / 'mocca'

    p_venv_lib = p_MOCCaPy / '.venv/lib'
    for p in p_venv_lib.glob('python*'):
        p_site_packages_mocca = (p / 'site-packages/mocca').resolve()
        p_site_packages_mocca.symlink_to(p_MOCCaPy_mocca)
        print(f"{p_site_packages_mocca} -> {p_MOCCaPy_mocca}")

