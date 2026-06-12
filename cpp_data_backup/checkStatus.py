import os
from sys import argv

numJobs = int(argv[1])
rowInd = int(argv[2])

count = 0
missingJobs = []
for jobInd in range(numJobs):
    if os.path.exists(f"cpp{jobInd}_{rowInd}.tr"):
        count += 1
    else:
        missingJobs.append(jobInd)

if count == numJobs:
    print("All jobs are done!")
else:
    print(f"Completed Jobs: {count}")
    print("Missing Job Inds: ", end = '')
    for i in missingJobs:
        print(i, end = ' ')
