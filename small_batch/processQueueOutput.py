import os
from sys import argv

filename = argv[1]

with open(filename, "r") as f:
    lines = [line for line in f.readlines() if line.startswith("henry.morin")]

with open(filename, "w") as f:
    for line in lines:
        if line[70] == "1":
            f.write(line)
