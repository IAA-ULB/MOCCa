This folder contains a few scripts for managing the Python virtual environment (beyond `uv`). 

To set up the virtual environment use `uv venv`,

# `create_soft_links_in_venv_bin.py`
```shell
> create_soft_links_in_venv_bin.py my_folder
```
This adds the directory structure of `my_folder` to the virtual environment `MOCCaPy/.venv` and creates soft links to all files under `my_folder`. This allows to isolate of additions to .venv that. 

# `dev_install_mocca.py`

This script adds a soft link to `MOCCapPy/mocca` in all `site-packages` folders under `.venv/bin` you can run `mocca` scripts and tests. 