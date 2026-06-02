import os
from globalVars import *
from sys import argv

def generateSubFile(jobInd, rowInd, difficulty = "M"):
    if not os.path.isdir(f"jobFiles"):
        os.mkdir(f"jobFiles")

    if not os.path.isdir(f"jobFiles/job{jobInd}"):
        os.mkdir(f"jobFiles/job{jobInd}")

    if difficulty == "L":
        length = "Long"
        data = "60"
    elif difficulty == "M":
        length = "Medium"
        data = "15"
    else:
        length = "Medium"
        data = "15"
    
    if version != 0:
        suffix = f"_v{version}.npy"
    else:
        suffix = ".npy"
    
    subFile = open(f"jobFiles/job{jobInd}/job{jobInd}_row{rowInd}.sub", "w")
     
    subFile.write(f"""container_image = osdf:///ospool/ap40/data/henry.morin/my-container.sif

executable = ../../main.sh
transfer_input_files = {prefix}/small_batch/FindIndicatorMats.py, {prefix}/small_batch/globalVars.py, {prefix}/small_batch/data/intMat.npy, {prefix}/small_batch/data/genus.npy, {prefix}/small_batch/data/aBounds.npy, {prefix}/{projectName}_data/{projectName}{jobInd}_{rowInd}{suffix}
transfer_output_files = "{projectName}{jobInd}_{rowInd + 1}{suffix}"
transfer_output_remaps = "{projectName}{jobInd}_{rowInd + 1}{suffix} = {prefix}/{projectName}_data/{projectName}{jobInd}_{rowInd + 1}{suffix}"
arguments = {jobInd} {rowInd} {version}
 
error = ../../errors/job{jobInd}_{rowInd}.error
log = ../../logs/job{jobInd}_{rowInd}.log
output = ../../outs/job{jobInd}_{rowInd}.out


+JobDurationCategory = "{length}"
    
request_cpus    = 1 
request_memory  = {data}GB
request_disk    = {data}GB
    
queue 1""")
    
    subFile.close()

if __name__ == "__main__":
    jobInd = int(argv[1])
    rowInd = int(argv[2])
    try:
        difficulty = argv[3]
    except:
        difficulty = "m"
    generateSubFile(rowInd, jobInd, difficulty)
