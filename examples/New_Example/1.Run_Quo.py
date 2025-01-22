import os
import pathlib
from Params1 import *
import sys
args = " ".join(sys.argv[1:])  
print(args)
command = (
    f"powershell.exe wsl env LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libmkl_def.so:"
    f"/usr/lib/x86_64-linux-gnu/libmkl_avx2.so:"
    f"/usr/lib/x86_64-linux-gnu/libmkl_core.so:"
    f"/usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so:"
    f"/usr/lib/x86_64-linux-gnu/libmkl_intel_thread.so:"
    f"/usr/lib/x86_64-linux-gnu/libiomp5.so python3.12 1.Quo_Model.py {args}"
)
os.system(command)