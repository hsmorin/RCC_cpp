import os
from sys import argv
from globalVars import *
from genSub import generateSubFile
import shutil


def genDF(jobInd, length = "L", rowMax = N):
    if not os.path.isdir(f"jobFiles"):
        os.mkdir(f"jobFiles")

    if version != 0:
        suffix = f"_v{version}.npy"
    else:
        suffix = ".npy"

    try:
        shutil.rmtree(f"./jobFiles/job{jobInd}")
    except:
        pass

    os.mkdir(f"./jobFiles/job{jobInd}")

    dagFile = open(f"./jobFiles/job{jobInd}/job{jobInd}.dag", "w")
    graphString = ""

    for rowInd in reversed(range(rowMax)):
        dagFile.write(f"JOB job{jobInd}_row{rowInd} job{jobInd}_row{rowInd}.sub\n")
        generateSubFile(jobInd, rowInd, length)
        if os.path.isfile(f"{prefix}/small_batch_data/{projectName}{jobInd}_{rowInd}{suffix}"):
            lowInd = rowInd
            dagFile.write("\n")
            break

#dagFile.write(f"SCRIPT POST ALL_NODES ../clean_jobs.sh 1003\n")
    for rowInd in range(lowInd, rowMax - 1):
        dagFile.write(f"PARENT job{jobInd}_row{rowInd} CHILD job{jobInd}_row{rowInd + 1}\n")
    dagFile.close()


if __name__ == "__main__":
    jobInd = int(argv[1])

    try:
        length = argv[2]
    except:
        length = "L"

    try:
        rowMax = int(argv[3])
    except:
        rowMax = N
    
    genDF(jobInd, length, rowMax)

