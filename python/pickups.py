"""Sets pickups for image slicer IFU. The sag of the subpupil mirrors depends
on the parameters passed to the image slicer as well as the on-axis distance between
these two surfaces. The parameters can be shared through pickups. Either set them
manually or run this script to do so.

The script will iterate through each row and find pairs of "Image Slicer" and
"SubpupilMiro" surfaces. Then, it will add the appropriate pickups.

Ellen Lee
"""

import numpy as np
import matplotlib.pyplot as plt