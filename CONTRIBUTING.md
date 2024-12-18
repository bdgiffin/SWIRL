## Changing OpenSeesPy to use a locally modified version of OpenSees

The [OpenSees documentation](https://opensees.github.io/OpenSeesDocumentation/developer/build.html) provides instructions on how to modify the local version of OpenSeesPy that gets used when running Python scripts that import the `openseespy` module. The basic steps are described below:
1. Build your local version of OpenSeesPy, following the build instructions provided in the OpenSees documentation. Within the OpenSees `build` subdirectory, this will produce the shared object library `build/lib/OpenSeesPy.dylib`. Leave this file named the way that it is (do not change the name to `opensees.so`). Navigate into this directory containing your new local version of `OpenSeesPy.dylib`, and use the `pwd` command to find out what the `FULL/SYSTEM/PATH/TO/build/lib/OpenSeesPy.dylib` is; make a note of this for use in step 4. 
2. Assuming you've already installed the release version of OpenSeesPy on your computer (using `pip install` from the PyPI), find where this version of OpenSeesPy is installed (on my Mac, it is located in `/opt/homebrew/lib/python3.13/site-packages/openseespymac/`.)
3. Navigate into the `openseespy`(`mac`) directory described in step 2. You should see a shared object library file named `opensees.so`. Do not delete/replace this file; instead, rename it to something else (e.g. `old_opensees.so`). If you want to revert your local modifications to OpenSeesPy, you only need to change the name back to the original name `opensees.so`.
4. Within this same directory, create a "symbolic link" to your local version of `OpenSeesPy.dylib` via:
   ```
   ln -s FULL/SYSTEM/PATH/TO/build/lib/OpenSeesPy.dylib opensees.so
   ```
5. Try running a Python script that imports `openseespy`. It should now be using your locally modified version of OpenSeesPy!

> NOTE: On newer Macbooks using Apple silicon with the arm64 architecture, it appears that the packaged version of OpenSeesPy is built with the x86_64 architecture (resulting in errors when attempting to run Python scripts that make use of OpenSeesPy). Replacing the version of `opensees.so` with a locally compiled version built natively on the Mac arm64 architecture (following the steps suggested above) appears to fix this error.

## Fixing errors in pyexodus version 0.1.5

The current release version of [pyexodus (0.1.5)](https://pypi.org/project/pyexodus/) installed via PyPI produces errors when run using python3.13 (the version of Python with which OpenSeesPy is built and executed). To fix these errors, you will need to directly modify the pyexodus source files installed on your computer. On my Mac, these files are located in `/opt/homebrew/lib/python3.13/site-packages/pyexodus/`. Within this folder, you will find the main pyexodus source file name `core.py` (this is the file that will need to be modified to fix the errors). Open this file, and replace all instances of `np.string_` with `np.bytes_`. For example, the following line (273) in `core.py`:
```
self._f.variables[var_name].attrs["elem_type"] = np.string_(elemType)
```
should be replaced by:
```
self._f.variables[var_name].attrs["elem_type"] = np.bytes_(elemType)
```
There should be a total of 4 such lines in `core.py` that need to be modified in this way. Once you have made these modifications, save the changes to `core.py`, and attempt to run a Python script using python3.13 and pyexodus; it should now work correctly.