# Installing ELPA on lumi

From Tor Skovgaard at support@lumi-supercomputer.eu

Hey Engelbert

The pre-installed version of ELPA does indeed look very old. If you look at
our software library you can however see that there a user-installable version
that is much newer. I would suggest that you try installing that instead. see this [link](https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Flumi-supercomputer.github.io%2FLUMI-EasyBuild-docs%2Fe%2FELPA%2F%23user-installable-modules-and-easyconfigs&data=05%7C02%7Cengelbert.tijskens%40uantwerpen.be%7C80f048dba9054d506c6d08ddc8412ab5%7C792e08fb2d544a8eaf72202548136ef6%7C0%7C0%7C638886904393514834%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=795zAsBek%2FeEKuii%2FTifBeRWENG0BvC0SG%2F9%2FywUrEY%3D&reserved=0)


Best Wishes
Tor - LUST

we will try to build:

## EasyConfig ELPA-2024.05.001-cpeGNU-24.03-CPU.eb

This will build ELPA/2024.05.001-cpeGNU-24.03-CPU. [see](https://lumi-supercomputer.github.io/LUMI-EasyBuild-docs/e/ELPA/ELPA-2024.05.001-cpeGNU-24.03-CPU/) 

the original easyconfig is here : `/appl/lumi/LUMI-EasyBuild-contrib/easybuild/easyconfigs/e/ELPA/ELPA-2024.05.001-cpeGNU-24.03-CPU.eb`

``` shell
> eb ELPA-2024.05.001-cpeGNU-24.03-CPU.eb -r
== Temporary log file in case of crash /run/user/327000497/easybuild/tmp/eb-tdl_3zp6/easybuild-5k75iz3a.log
== resolving dependencies ...
== processing EasyBuild easyconfig /appl/lumi/LUMI-EasyBuild-contrib/easybuild/easyconfigs/e/ELPA/ELPA-2024.05.001-cpeGNU-24.03-CPU.eb
== building and installing ELPA/2024.05.001-cpeGNU-24.03-CPU...
== fetching files...
== creating build dir, resetting environment...
== starting iteration #0 ...
== unpacking...
== patching...
== preparing...
== ... (took 5 secs)
== configuring...
== ... (took 30 secs)
== building...
== ... (took 3 mins 11 secs)
== testing...
== installing...
== ... (took 9 secs)
== taking care of extensions...
== creating build dir, resetting environment...
== starting iteration #1 ...
== unpacking...
== patching...
== preparing...
== ... (took 3 secs)
== configuring...
== ... (took 31 secs)
== building...
== ... (took 3 mins 18 secs)
== testing...
== installing...
== ... (took 10 secs)
== taking care of extensions...
== restore after iterating...
== postprocessing...
== sanity checking...
== ... (took 2 secs)
== cleaning up...
== creating module...
== ... (took 2 secs)
== permissions...
== packaging...
== COMPLETED: Installation ended successfully (took 8 mins 9 secs)
== Results of the build can be found in the log file(s) 
/users/entijske/EasyBuild/SW/LUMI-24.03/C/ELPA/2024.05.001-cpeGNU-24.03-CPU/easybuild/easybuild-ELPA-2024.05.001-20250721.161120.log

== Build succeeded for 1 out of 1
== [end-hook] Clearing Lmod cache directory /users/entijske/.lmod.d/.cache
== [end-hook] Clearing Lmod cache directory /users/entijske/.cache/lmod
== Temporary log file(s) /run/user/327000497/easybuild/tmp/eb-tdl_3zp6/easybuild-5k75iz3a.log* have been removed.
== Temporary directory /run/user/327000497/easybuild/tmp/eb-tdl_3zp6 has been removed.
```
That went well.

Let's try to enable the python wrapper:

``` shell
> cd workspace/tantalus_full
> cd experiments/env/lumi
> cp `/appl/lumi/LUMI-EasyBuild-contrib/easybuild/easyconfigs/e/ELPA/ELPA-2024.05.001-cpeGNU-24.03-CPU.eb` .
eb ./
```

Add `--enable-python` and `--enable-python-tests` to `local_common_configopts`

``` shell
> eb ./ELPA-2024.05.001-cpeGNU-24.03-CPU.eb -r
== Temporary log file in case of crash /run/user/327000497/easybuild/tmp/eb-_q3vuaiy/easybuild-naxue4jy.log
== ELPA/2024.05.001-cpeGNU-24.03-CPU is already installed (module found), skipping
== No easyconfigs left to be built.

== Build succeeded for 0 out of 0
== Temporary log file(s) /run/user/327000497/easybuild/tmp/eb-_q3vuaiy/easybuild-naxue4jy.log* have been removed.
== Temporary directory /run/user/327000497/easybuild/tmp/eb-_q3vuaiy has been removed.
```

OK, we probably have to rename the easyconfig file. 

Nope, same error, renaming the file didn't help. we remove instead the previous installation and run it again

still the same error. Somehow, it remembers that it is already installed.

Add `--rebuild` option

``` shell
> eb ./ELPA-2024.05.001-cpeGNU-24.03-CPU.eb -r
```

Apparently, we need not to add `--enable-python` and `--enable-python-tests` to `local_common_configopts`, but to `configopts`:

```
onfigopts = [
    local_common_configopts + ' --disable-openmp ',
    local_common_configopts + ' --enable-openmp ' ,
    local_common_configopts + ' --enable-python --enable-python-tests', # <- this line was added.
]
```
That didn't work either. pfff. After checking the build file it seems to be using the system python...

guess we'll have to resort to manual installation...

## EasyConfig ELPA-2024.05.001-cpeGNU-24.03-rocm.eb

This will build ELPA/2024.05.001-cpeGNU-24.03-rocm the gpu version