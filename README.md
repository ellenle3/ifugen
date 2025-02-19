# ifugen

This repository contains Zemax DLLs for generating image slicer integral field unit
(IFUs). The mirror arrays are implemented fully sequentially, so only one configuration
is needed to specify the entire IFU. Non-sequential versions of these surfaces are
also available.

## Installation

1. Download the latest release. See "Releases" in the sidebar to the right.
2. Navigate to the folder for Zemax data folder for DLL files on your PC. This should be something like "Documents\Zemax\DLL".
3. Copy us_slicer_std.dll and us_slicer_custom.dll into the Surfaces subdirectory.
4. Copy ImageSlicerStd.dll and ImageSlicerCustom.dll into the Objects subdirectory.

That's it! You should now be able to load the DLLs as user-defined surfaces
(or objects) in Zemax. Treat these files with the same level of caution that you
would downloading an EXE file from the internet. Alternatively, build the DLL from
source yourself.

## Dependencies

This project relies on the C standard library and the Windows API (windows.h).
These should already be on the machine you use for Zemax so you should not need
to install any dependencies yourself.

## References

Please check the documentation for this project (not posted yet...)