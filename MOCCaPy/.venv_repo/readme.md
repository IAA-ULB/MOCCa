We manage the Python virtual environment with [uv](https://docs.astral.sh/uv). This folder contains a few scripts further augment the Python virtual environment:  

### `create_soft_links_in_venv_bin.py`

```shell
> create_soft_links_in_venv_bin.py my_folder
```

This adds the directory structure of `my_folder` to the virtual environment `MOCCaPy/.venv` and creates soft links to all files under `my_folder`. This allows to isolate of additions to .venv that. 

### `dev_install_mocca.py`

```shell
> dev_install_mocca.py
```
This script adds a soft link to `MOCCapPy/mocca` in all `site-packages` folders under `.venv/bin` you can run `mocca` scripts and tests. 