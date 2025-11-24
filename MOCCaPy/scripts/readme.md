```shell
(MOCCaPy) etijskens@MacOSm3@ [584] ~/workspace/tantalus_full/MOCCaPy/scripts
> ll
total 64
drwxr-xr-x@ 4 etijskens  staff    128 Nov 23 22:02 __pycache__/
-rw-r--r--@ 1 etijskens  staff  11237 Oct 19 23:13 ca40_discussion.py
drwxr-xr-x@ 7 etijskens  staff    224 Oct 16 13:23 Ca40.f90/
-rw-r--r--@ 1 etijskens  staff   1642 Nov 23 22:38 ca40.py
-rw-r--r--@ 1 etijskens  staff    601 Nov 23 22:13 insert_path.py
lrwxr-xr-x@ 1 etijskens  staff     54 Nov 23 22:41 mocca@ -> /Users/etijskens/workspace/tantalus_full/MOCCaPy/mocca
-rw-r--r--@ 1 etijskens  staff     25 Nov 23 23:14 readme.md
-rw-r--r--@ 1 etijskens  staff    180 Nov 23 22:03 sandbox.py
-rw-r--r--@ 1 etijskens  staff    193 Nov 23 21:56 tryout.py
```
the soft link `mocca` ensures that scripts in this folder can import mocca functionality. This relies on the fact that running a script in Python automatically adds the script's folder to `$PYTHONPATH`. 