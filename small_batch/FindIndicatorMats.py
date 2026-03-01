from itertools import combinations
from sys import argv
from typing import Counter
import math
from ortools.sat.python import cp_model
import numpy as np
import os
from globalVars import *



# generates possible lists [-b_1, ..., -b_N], with (-b_i)'s in ascending order (i.e. b_i's in descending order)
# uses depth-first search to avoid creating and storing huge lists
def genBis(N, selfInt, genus, a, currentRow = [], currSelfIntSum = 0, currAdjSum = 0):

    # base cases
    if a <= 0:
        if a == 0 or (selfInt < -2 and genus == 0): # Lemma 3.4
            numOnes = 2 * a - selfInt - 1
            if 0 <= numOnes <= N - 1:
                yield [-1] * numOnes + [0] * (N - 1 - numOnes) + [-a + 1]
                return
        return
    currLen = len(currentRow)
    if currLen == N: # row is complete!
        yield [-bi for bi in currentRow]
        return
    
    # recursive case - run method for all possible next b_i's such that self-intersection and adjunction are still feasible
    numBisLeft = N - currLen
    selfIntSumLeft = a ** 2 - selfInt - currSelfIntSum
    adjSumLeft = (a - 1) * (a - 2) - 2 * genus - currAdjSum
    if selfIntSumLeft < 0 or adjSumLeft < 0:
        return
    adjLowerBound = 0 if adjSumLeft == 0 else math.sqrt(adjSumLeft / numBisLeft + 0.25) + 0.5 # fix quadratic-ness of adjunction summands (0 at both 0 and 1)
    
    nextBiMin = math.ceil(max(math.sqrt(selfIntSumLeft / numBisLeft), adjLowerBound, 0))
    nextBiMax = math.floor(min(math.sqrt(selfIntSumLeft), math.sqrt(adjSumLeft + 0.25) + 0.5, currentRow[-1] if currentRow else 1e6))

    if nextBiMin > nextBiMax:
        return
    middleBi = (nextBiMin + nextBiMax) // 2
    yield from genBis(N, selfInt, genus, a, currentRow + [middleBi], currSelfIntSum + middleBi ** 2, currAdjSum + (middleBi ** 2 - middleBi))
    biOffset = 1
    while biOffset <= (nextBiMax - nextBiMin + 1) // 2:
        nextBi = middleBi - biOffset
        if nextBi >= nextBiMin:
            yield from genBis(N, selfInt, genus, a, currentRow + [nextBi], currSelfIntSum + nextBi ** 2, currAdjSum + (nextBi ** 2 - nextBi))
        nextBi = middleBi + biOffset
        if nextBi <= nextBiMax:
            yield from genBis(N, selfInt, genus, a, currentRow + [nextBi], currSelfIntSum + nextBi ** 2, currAdjSum + (nextBi ** 2 - nextBi))
        biOffset += 1

# runs genBis for a range of values of "a" and generates solution lists in the form [a, -b_1, -b_2, ..., -b_N], with (-b_i)'s ascending
def rowCands(N, selfInt, genus, aMin, aMax):
    aLowerBound = 0 if selfInt >= -2 else math.ceil((selfInt + 1) / 2) # Lemma 3.4
    startA = max(aLowerBound, aMin)
    for a in range(startA, int(aMax + 1)):
        biVecs = genBis(N, selfInt, genus, a)
        for biVec in biVecs:
            yield [a] + biVec

# returns boolean masks to trim two lists of possible rows using the desired intersection number and the rearrangement inequality
# assumes all -b_i coefficients in the input lists are sorted in ascending order
def rearrIneqTrim(goalIntNum, coeffRowList1, coeffRowList2):
    arr1 = np.array(coeffRowList1, ndmin = 2)
    arr2 = np.array(coeffRowList2, ndmin = 2)
    pairwiseAProducts = arr1[:, 0:1] @ arr2[:, 0:1].T
    biArr1, biArr2T = arr1[:, 1:], arr2[:, 1:].T
    pairwiseMinIntNums, pairwiseMaxIntNums = pairwiseAProducts - biArr1 @ biArr2T, pairwiseAProducts - biArr1[:,::-1] @ biArr2T
    pairsValid = (pairwiseMinIntNums <= goalIntNum) & (goalIntNum <= pairwiseMaxIntNums)
    list1Valid, list2Valid = np.any(pairsValid, axis = 1), np.any(pairsValid, axis = 0)
    return list1Valid, list2Valid

# performs rearrangement inequality trimming in (chunkSize x chunkSize) chunks, saving memory when working with very long lists
def rearrIneqTrimChunked(goalIntNum, coeffRowList1, coeffRowList2):
    chunkSize = int(2000 + 5000 * math.exp(-7e-12 * (len(coeffRowList1) * len(coeffRowList2)))) # scales down slowly with list lengths
    n1, n2 = len(coeffRowList1), len(coeffRowList2)
    list1Valid = np.zeros(n1, dtype = bool)
    list2Valid = np.zeros(n2, dtype = bool)

    for i1 in range(0, n1, chunkSize):
        for i2 in range(0, n2, chunkSize):
            chunk1 = coeffRowList1[i1 : i1 + chunkSize]
            chunk2 = coeffRowList2[i2 : i2 + chunkSize]
            chunk1Valid, chunk2Valid = rearrIneqTrim(goalIntNum, chunk1, chunk2)
            list1Valid[i1 : i1 + len(chunk1)] |= chunk1Valid
            list2Valid[i2 : i2 + len(chunk2)] |= chunk2Valid
    
    return list1Valid, list2Valid

# sorts the b_i columns of a given coefficient matrix
def sortBis(mat):
    colInds = np.lexsort(mat[:, 1:])
    mat[:, 1:] = mat[:, 1:][:, colInds]
    return mat
  
# takes an array of desired dot products between the rows of prevBiMat and the elements of the counter bCtr.
# generates lists containing the elements of bCtr, representing permutations that have the desired dot products (goalDots) with prevBiMat
def validBiPerms(goalDots, prevBiMat, bCtr):
    # base cases - use rearrangement inequality to cut off recursive searching and handle 1-unique-element case
    sortedPrevBis = np.sort(prevBiMat, axis = 1)
    sortedNewBis = np.array(sorted(bCtr.elements()))
    maxDots = sortedPrevBis @ sortedNewBis
    minDots = sortedPrevBis @ np.flip(sortedNewBis)
    if np.any(goalDots < minDots) or np.any(goalDots > maxDots):
        return
    if len(bCtr) == 1: # bCtr has one unique element; only reached if (all) minDots <= goalDots <= maxDots, so dots must be correct
        yield list(bCtr.elements())
        return
    # recursive case - minDots <= goalDots <= maxDots, but this is a nontrivial range
    numCols = len(sortedNewBis)
    biToAssign = max(bCtr, key = abs) # largest absolute value -> most impact on dot product -> earlier branch pruning
    numCopies = bCtr[biToAssign]
    del bCtr[biToAssign] # remove the bi that is assigned before making recursive calls
    for colInds in combinations(range(numCols), numCopies): # indices to which copies of biToAssign are assigned

        matchedCols = prevBiMat[:,colInds]
        dotContributions = biToAssign * np.sum(matchedCols, axis = 1)
        newGoalDots = goalDots - dotContributions
        newPrevBiMat = np.delete(prevBiMat, colInds, axis = 1)

        subPerms = validBiPerms(newGoalDots, newPrevBiMat, bCtr)
        for subPerm in subPerms:
            fullPerm = subPerm.copy()
            for i in colInds: # insert left to right to restore original order
                fullPerm.insert(i, biToAssign)
            yield fullPerm
    bCtr[biToAssign] = numCopies # put back the assigned b_i's for future recursive calls

# takes in an intersection matrix, list of coefficient lists (without permutations), and the index of the row being added
# assumes lists in coeffLists are in ascending order, for proper functioning of validBiPerms
# recursively returns a set of possible submatrices (tuples of tuples) with a (lastRowInd + 1) rows, fitting the given intersection matrix
def viableSubMats(intMat, coeffLists, lastRowInd):
    if lastRowInd == 0:
        return set(tuple(tuple(row)) for row in coeffLists[0])
    prevMats = viableSubMats(intMat, coeffLists, lastRowInd - 1)
    prevMats = [np.array(mat, ndmin = 2) for mat in prevMats] # convert from set of tuple-tuples to list of np arrays
    newRowCands = np.array(coeffLists[lastRowInd], dtype = int) # dtype specification makes sure to make 2d array
    subMats = set()

    for prevMat in prevMats:
        prevBis = prevMat[:,1:]
        for newRow in newRowCands:
            newA, newBis = newRow[0], newRow[1:]
            if newA < 0 and np.any(prevMat[:,0] < 0): # Lemma 4.2
                newRowPerms = []
            else:
                goalBDots = prevMat[:,0] * newA - intMat[lastRowInd][:lastRowInd]
                newRowPerms = validBiPerms(goalBDots, prevBis, Counter(newBis))
            for perm in newRowPerms:
                newMat = sortBis(np.vstack([prevMat, np.concatenate([[newA], perm])]))
                subMats.add(tuple(tuple(row) for row in newMat))

    # if lastRowInd == 2:
    #     for j, mat in enumerate(subMats):
    #         print(f"Matrix {j}: ")
    #         print(*mat, sep = "\n")
    return subMats

# takes in a list of submatrices as arrays, the self-intersection, genus, and a-bounds of the new surface, and the mutual intersections between
# the new row and the existing ones, and returns a list of possible matrix extensions, i.e. matrices with one more row than subMats
def matExtensions(subMats, interNums, selfInt, genus, aMin, aMax):
    newMatSet = set()
    N = len(subMats[0][0]) - 1

    aLowerBound = 0 if selfInt >= -2 else math.ceil((selfInt + 1) / 2) # Lemma 3.4
    startA = max(aLowerBound, aMin)
    for a in range(startA, int(aMax + 1)):
        biVecs = genBis(N, selfInt, genus, a)
        for biVec in biVecs:
            for subMat in subMats:
                if a < 0 and np.any(subMat[:,0] < 0): # Lemma 4.2
                    newRowPerms = []
                else:
                    goalBDots = subMat[:,0] * a - interNums
                    newRowPerms = validBiPerms(goalBDots, subMat[:,1:], Counter(biVec))
                for perm in newRowPerms:
                    newMat = np.vstack([subMat, np.zeros((1, N + 1))])
                    newMat[-1, 0], newMat[-1, 1:] = a, perm
                    newMatSet.add(tuple(tuple(row) for row in sortBis(newMat)))

    return [np.array(mat, ndmin = 2) for mat in newMatSet] # convert to list of numpy arrays

# uses CP-SAT solver to compile list of rows [a, -b_1, ..., -b_N] (with a from 1 to aMax inclusive) that have intersection numbers interNums with the rows of subMat
# solver finds b_i's, then solutions are flipped to -b_i's as they are compiled, and solutions with a <= 0 are manually added
def solveNextRow(subMat: np.ndarray, interNums: np.ndarray, selfInt: int, genus: int, aMax: int):

    numPrevRows, N = len(subMat), len(subMat[0]) - 1

    model = cp_model.CpModel()

    a = model.NewIntVar(1, aMax, "a")
    biUB = min(3 * aMax - selfInt + 2 * genus - 2, int(math.sqrt(aMax ** 2 - selfInt)))
    b = [model.NewIntVar(0, biUB, f"b_{i}") for i in range(N)]
    
    # intersection number constraints: for all rowInd, subMat[rowInd] . [a, b_1, ..., b_N] = interNums[rowInd]
    for rowInd in range(numPrevRows):
        expr = int(subMat[rowInd, 0]) * a
        for i in range(N):
            expr += int(subMat[rowInd, i + 1]) * b[i]
        model.Add(expr == int(interNums[rowInd]))

    # adjunction formula
    model.Add(sum(b) == 3 * a - selfInt + 2 * genus - 2)

    # self-intersection
    bisSquared = []
    for i in range(N):
        bSquared = model.NewIntVar(0, biUB ** 2, f"b_{i}_squared")
        model.AddMultiplicationEquality(bSquared, [b[i], b[i]])
        bisSquared.append(bSquared)

    aSquared = model.NewIntVar(1, aMax ** 2, "a_squared")
    model.AddMultiplicationEquality(aSquared, [a, a])

    model.Add(sum(bisSquared) == aSquared - selfInt)

    # symmetry breaking for repeated columns in subMat; assumes that subMat has columns in sortBis() order
    sameAsRightNeighbor = [] # saving this for later use in a <= 0 case
    for i in range(N - 1):
        if np.array_equal(subMat[:, i+1], subMat[:, i+2]):  # columns for b[i] vs b[i+1]
            model.Add(b[i] >= b[i+1])
            sameAsRightNeighbor.append(i)

    # solution collection

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True

    sols = []

    class SolutionCollector(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            sols.append([self.Value(a)] + [-self.Value(b[i]) for i in range(N)])
    
    collector = SolutionCollector()
    solver.Solve(model, collector)

    # handle a <= 0 case, with only one solution (plus permutations) per a
    aMin = 0 if selfInt >= -2 else math.ceil((selfInt + 1) / 2) # Lemma 3.4
    for currA in range(aMin, 1):
        numOnes = 2 * currA - selfInt - 1
        if 0 <= numOnes <= N - 1:
            goalBDots = subMat[:,0] * currA - interNums
            negBCounter = Counter({-1: numOnes, 0: N - 1 - numOnes, (- currA + 1): 1})
            for perm in validBiPerms(goalBDots, subMat[:,1:], negBCounter):
                discard = False
                for i in sameAsRightNeighbor:
                    if perm[i] < perm[i + 1]:
                        discard = True
                        break
                if not discard:
                    sols.append([currA] + perm)

    return sols

# 4 possible key functions (1-4 in aggOrder) for sorting surfaces using coeffLists, besides length of list

# returns average number of unique b_i values in a list in coeffList
def avgUniqueVals(coeffList):
    return sum(len(set(row[1:])) for row in coeffList) / len(coeffList)

# returns average absolute value of b_i's in coeffList
def avgAbsVal(coeffList):
    return sum(sum(abs(bi) for bi in row[1:]) for row in coeffList) / (len(coeffList) * (len(coeffList[0]) - 1))

# returns average absolute value of UNIQUE b_i's in coeffList (i.e. unweighted within each list)
def avgUniqueAbs(coeffList):
    biSets = [set(row[1:]) for row in coeffList]
    return sum(sum(abs(bi) for bi in biSet) / len(biSet) for biSet in biSets) / len(coeffList)

# returns NEGATIVE of average variance of b_i's in a list in coeffList (negative because higher variance should be earlier)
def avgBiVar(coeffList):
    return -1 * sum(np.var(row[1:]) for row in coeffList) / len(coeffList)

# aggregation of key functions (see above) L, L1, L2, L3, L4, 1, 2, 3, 4, and default by taking an average of the rankings they induce;
# L1, L2, 1, 2 are usually good
def aggOrder(coeffLists):
    invIndsL = np.argsort(sorted(range(len(coeffLists)), key = lambda i: len(coeffLists[i])))
    invInds1 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: avgUniqueVals(coeffLists[i])))
    invInds2 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: avgAbsVal(coeffLists[i])))
    invInds3 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: avgUniqueAbs(coeffLists[i])))
    invInds4 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: avgBiVar(coeffLists[i])))
    invIndsL1 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: (len(coeffLists[i]),avgUniqueVals(coeffLists[i]))))
    invIndsL2 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: (len(coeffLists[i]),avgAbsVal(coeffLists[i]))))
    invIndsL3 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: (len(coeffLists[i]),avgUniqueAbs(coeffLists[i]))))
    invIndsL4 = np.argsort(sorted(range(len(coeffLists)), key = lambda i: (len(coeffLists[i]),avgBiVar(coeffLists[i]))))
    return sorted(range(len(coeffLists)),
                  key = lambda i: 0 * invInds1[i]  + 0 * invInds2[i]  + 0 * invInds3[i]  + 0 * invInds4[i] 
                                + 1 * invIndsL1[i] + 1 * invIndsL2[i] + 0 * invIndsL3[i] + 0 * invIndsL4[i]
                                + 0 * invIndsL[i]  + 0 * i)

# returns list of possible coefficient matrices that agree with given intersection matrix, Lemma 3.3, Lemma 3.4, and adjunction formula
# excludes permutations of b_i (non-first) columns
def FindIMats(N, intMat, genera, aMaxes):
    numSurfs = len(genera)
    rowOrdInds = sorted(range(numSurfs), key = lambda i: aMaxes[i] ** 2 - intMat[i][i], reverse = False)
    intMatReOrd = intMat[rowOrdInds][:,rowOrdInds]
    coeffListsReOrd = list()
    generaReOrd = genera[rowOrdInds]
    aMaxesReOrd = aMaxes[rowOrdInds]
    
    for i in range(numSurfs):
        sols = list(rowCands(N, intMatReOrd[i,i], generaReOrd[i], -1e6, aMaxesReOrd[i]))
        coeffListsReOrd.append(sols)
        for j in range(i):
            listIValid, listJValid = rearrIneqTrimChunked(intMatReOrd[i][j], coeffListsReOrd[i], coeffListsReOrd[j])
            coeffListsReOrd[i] = [coeffListsReOrd[i][ind] for ind in range(len(coeffListsReOrd[i])) if listIValid[ind]]
            coeffListsReOrd[j] = [coeffListsReOrd[j][ind] for ind in range(len(coeffListsReOrd[j])) if listJValid[ind]]

    restoreOrder = np.argsort(rowOrdInds)
    coeffLists = [coeffListsReOrd[i] for i in restoreOrder]

    # key function(s) for reordering rows is some combination of avgUniqueVals, avgAbsVal, avgUniqueAbs, and avgBiVar
    rowOrdInds = aggOrder(coeffLists)
    intMatReOrd = intMat[rowOrdInds][:,rowOrdInds]
    coeffListsReOrd = [coeffLists[i] for i in rowOrdInds]
    solsReOrd = viableSubMats(intMatReOrd, coeffListsReOrd, numSurfs - 1)
    restoreOrder = np.argsort(rowOrdInds)
    finalSols = [np.array(solMat)[restoreOrder] for solMat in solsReOrd]
    return finalSols

def main():
   
    try:
        rowInd = int(argv[2])
    except:
        rowInd = -1

    try:
        version = int(argv[3])
    except:
        pass

    if version != 0:
        suffix = f"_v{version}.npy"
    else:
        suffix = ".npy"

    print(f"Job for row index {rowInd} and version {version} with suffix {suffix}.")
    Q = np.load("./intMat.npy")
    N = len(Q)
    aMaxes = np.load("./aBounds.npy")
    
    rowOrdInds = sorted(range(N), key = lambda i: aMaxes[i] ** 2 - Q[i][i])
    intMatReOrd = Q[rowOrdInds][:,rowOrdInds]
    aMaxesReOrd = aMaxes[rowOrdInds]
    generaReOrd = np.load("./genus.npy")[rowOrdInds]
    # [np.array(firstRow).reshape(1, -1) for firstRow in rowCands(N, intMatReOrd[0, 0], generaReOrd[0], aMin = -1e6, aMax = aMaxesReOrd[0])]
    if rowInd == -1:
        np.save(f"{projectName}0.npy", np.array([np.array(firstRow).reshape(1, -1) for firstRow in rowCands(N - 1, intMatReOrd[0, 0], generaReOrd[0], aMin = -1e6, aMax = aMaxesReOrd[0])]))
    else:
        jobInd = int(argv[1])
        newSubMatSet = set()
        subMats = np.load(f"./{projectName}{jobInd}_{rowInd}{suffix}")
        for subMat in subMats:
            nextRows = solveNextRow(subMat, interNums = intMatReOrd[rowInd + 1, :rowInd + 1], selfInt = intMatReOrd[rowInd + 1, rowInd + 1], genus = 0, aMax = aMaxesReOrd[rowInd + 1])
            newMats = [np.vstack([subMat, newRow]) for newRow in nextRows]
            for newMat in newMats:
                newSubMatSet.add(tuple(tuple(row) for row in sortBis(newMat)))
        subMats = [np.array(mat, ndmin = 2) for mat in newSubMatSet]
        np.save(f"{projectName}{jobInd}_{rowInd + 1}{suffix}", subMats)
        print(f"Job was sucessful with {len(subMats)} new matrices saved to {projectName}{jobInd}_{rowInd + 1}{suffix}")

if __name__ == "__main__":
    main()
