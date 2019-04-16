#!/usr/bin/python

import subprocess
import datetime

# k = constant parameter, that helps set threshold for scale of observation
# between larger and smaller components. (large k prefers large components)

tStart = datetime.datetime.now()
p = subprocess.run("./segment 0.5 500 20 4.ppm 4out.ppm", shell=True)
tEnd = datetime.datetime.now() - tStart
secs = tEnd.seconds + 24 * 3600 * tEnd.days + 1e-6 * tEnd.microseconds
print()
print("Total test time = %.2f secs." % secs)
print()
