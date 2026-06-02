from sys import argv
import numpy as np
from genDF import genDF

currentJobs = [int(i) for i in argv[1:]]

for jobInd in currentJobs:
    genDF(jobInd)
