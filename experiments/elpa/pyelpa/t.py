import numpy as np
import sys
import pyscalapack
from elpa import Elpa

scalapack = pyscalapack(
    "/apps/antwerpen/skylake/rocky9/ELPA/2024.05.001-intel-2024a/lib/libelpa.so",
    "/data/antwerpen/201/vsc20170/tantalus_full/experiments/elpa/pyelpa/libelpa_object_interface.so"
)  
with Elpa(scalapack) as e:
    print(f"body")
    e.set(b'na', 8)
    