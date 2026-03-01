from sys import argv
import numpy as np


jobNum = int(argv[1])
rowInd = int(argv[2])
numSplits = int(argv[3])

try:
    version = int(argv[4])
except:
    version = 0

try:
    newVersion = int(argv[5])
except:
    newVersion = version + 1

if version == 0:
    suffix = ".npy"
else:
    suffix = f"_v{version}.npy"

subMats = np.load(f"IMats{jobNum}_{rowInd}" + suffix)

for i, mats in enumerate(np.array_split(subMats, numSplits)):
    np.save(f"IMats{numSplits * jobNum + i}_{rowInd}_v{newVersion}.npy", mats)

    
