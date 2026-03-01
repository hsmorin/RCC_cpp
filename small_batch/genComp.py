import os
import shutil
from sys import argv
from sys import exit
projectName = argv[1]

N = argv[2]

if not os.path.isdir(f"osdf:///ospool/ap40/data/henry.morin/{projectName}"):
    os.mkdir(f"//ospool/ap40/data/henry.morin/{projectName}")
else:
    print("Error: Name in use, try again with a new name.")
    exit()

if not os.path.isdir("logs"):
    os.mkdir("logs")
if not os.path.isdir("outs"):
    os.mkdir("outs")
if not os.path.isdir("errors"):
    os.mkdir("errors")
init = open("globalVars.py", "w")

init.write("""global projectName
global N
global version

version = 0
projectName = \"""" + projectName + """\"
N = """ + str(N))
init.close()

shutil.copyfile("globalVars.py", "osdf:///ospool/ap40/data/henry.morin/{projectName}/globalVars.py")




