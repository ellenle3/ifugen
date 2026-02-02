# ifugen
**NOTE:** This repository is a work in progress. Compiled DLLs are not provided
here due to a few known issues that may cause Zemax to crash. You are welcome to
[compile the DLLs yourself](https://optics.ansys.com/hc/en-us/articles/42661743318291-How-to-compile-a-User-Defined-DLL)
and try them at your own risk.

This repository contains Zemax DLLs for generating image slicer integral field unit
(IFUs). The mirror arrays are implemented fully sequentially, so only one configuration
is needed to specify the entire IFU. Non-sequential versions of these surfaces are
also available.

<!-- ## Installation

1. Download the latest release. See "Releases" in the sidebar to the right.
2. Navigate to the folder for Zemax data folder for DLL files on your PC. This should be something like "Documents\Zemax\DLL".
3. Copy us_slicer_std.dll and us_slicer_custom.dll into the Surfaces subdirectory.
4. Copy ImageSlicerStd.dll and ImageSlicerCustom.dll into the Objects subdirectory.

That's it! You should now be able to load the DLLs as user-defined surfaces
(or objects) in Zemax. Treat these files with the same level of caution that you
would downloading an EXE file from the internet. Alternatively, build the DLL from
source yourself. -->

## Dependencies

This project relies on the C standard library and the Windows API (windows.h).
These should already be on the machine you use for Zemax so you should not need
to install any dependencies yourself.

## Disclaimer

It is probably apparent from the source code that I am not a software engineer.
Significant work remains to be done to make the software more accessible to users.
Note that this repository uses an MIT license, so I encourage others to repurpose
the code for their own use (including commercial). If you do so, please provide
a reference to this repository.

## Citation

A paper is in preparation for this repository. Information about citing this
paper will be shared here if (when?) it becomes available.