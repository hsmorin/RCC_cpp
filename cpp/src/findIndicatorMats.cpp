#include "absl/base/log_severity.h"
#include "absl/log/globals.h"
#include "ortools/base/init_google.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/util/sorted_interval_list.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <cmath>
#include <generator>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>
using namespace std;

void displayMatrix(const vector<vector<int>> &matrix) {
  for (const auto &row : matrix) {
    for (int val : row) {
      cout << val << " ";
    }
    cout << "\n";
  }
}

void displayMatrix(const vector<array<int, 11>> &matrix) {
  for (const auto &row : matrix) {
    for (int val : row) {
      cout << val << " ";
    }
    cout << "\n";
  }
}

void displayArray(const array<int, 11> &array) {
  for (const auto &val : array) {
    cout << val << " ";
  }
  cout << "\n";
}

void displayVector(const vector<int> &array) {
  for (const auto &val : array) {
    cout << val << " ";
  }
  cout << "\n";
}

vector<int> append(const vector<int> &v, int x) {
  vector<int> result = v;
  result.push_back(x);
  return result;
}

vector<vector<int>> sliceMat(const vector<vector<int>> &matrix,
                             int rowStart = 0, int rowEnd = -1,
                             int colnStart = 0, int colnEnd = -1) {
  vector<vector<int>> result;
  if (rowEnd == -1) {
    rowEnd = matrix.size();
  }
  if (colnEnd == -1) {
    colnEnd = matrix[0].size();
  }
  for (int rowInd = rowStart; rowInd < rowEnd; rowInd++) {
    vector<int> currentRow = matrix[rowInd];
    vector<int> arraySlice = {};
    for (int colnInd = colnStart; colnInd < colnEnd; colnInd++) {
      arraySlice.push_back(currentRow[colnInd]);
    }
    result.push_back(arraySlice);
  }
  return result;
}
vector<vector<int>> sliceMat(const vector<array<int, 11>> &matrix,
                             int rowStart = 0, int rowEnd = -1,
                             int colnStart = 0, int colnEnd = -1) {
  vector<vector<int>> result;
  if (rowEnd == -1) {
    rowEnd = matrix.size();
  }
  if (colnEnd == -1) {
    colnEnd = 11;
  }
  for (int rowInd = rowStart; rowInd < rowEnd; rowInd++) {
    array<int, 11> currentRow = matrix[rowInd];
    vector<int> arraySlice = {};
    for (int colnInd = colnStart; colnInd < colnEnd; colnInd++) {
      arraySlice.push_back(currentRow[colnInd]);
    }
    result.push_back(arraySlice);
  }
  return result;
}

generator<array<int, 10>> genBis(int N, int selfInt, int genus, int a,
                                 vector<int> currentRow = {},
                                 int currSelfIntSum = 0, int currAdjSum = 0) {

  // base cases
  if (a <= 0) {
    if ((a == 0) || ((selfInt < -2) && (genus == 0))) {
      int numOnes = 2 * a - selfInt - 1;
      if ((0 <= numOnes) && (numOnes <= N - 1)) {
        vector<int> bis(numOnes, -1);
        bis.insert(bis.end(), N - 1 - numOnes, 0);
        bis.push_back(1 - a);
        array<int, 10> bisArray;
        copy(bis.begin(), bis.end(), bisArray.begin());
        co_yield bisArray;
        co_return;
      }
    }
    co_return;
  }

  int currLen = currentRow.size();
  if (currLen == N) {
    array<int, 10> bisArray;
    copy(currentRow.begin(), currentRow.end(), bisArray.begin());
    transform(bisArray.begin(), bisArray.end(), bisArray.begin(),
              negate<int>());
    co_yield bisArray;
    co_return;
  }

  int numBisLeft = N - currLen;
  int selfIntSumLeft = a * a - selfInt - currSelfIntSum;
  int adjSumLeft = (a - 1) * (a - 2) - 2 * genus - currAdjSum;
  if ((selfIntSumLeft < 0) || (adjSumLeft < 0))
    co_return;
  int adjLowerBound = 0;
  if (adjSumLeft != 0)
    adjLowerBound = sqrt(((double)adjSumLeft / numBisLeft) + 0.25) + 0.5;

  int nextBiMin = ceil(max(
      {sqrt((double)selfIntSumLeft / numBisLeft), (double)adjLowerBound, 0.0}));
  double placeholderUpperBound;
  if (currentRow.empty())
    placeholderUpperBound = 1000000.0;
  else {
    placeholderUpperBound = static_cast<double>(currentRow.back());
  }
  int nextBiMax = floor(
      min({sqrt((double)selfIntSumLeft), sqrt((double)adjSumLeft + 0.25) + 0.5,
           placeholderUpperBound}));

  if (nextBiMin > nextBiMax)
    co_return;

  int middleBi = (nextBiMin + nextBiMax) / 2;
  for (const auto &row :
       genBis(N, selfInt, genus, a, append(currentRow, middleBi),
              currSelfIntSum + middleBi * middleBi,
              currAdjSum + middleBi * middleBi - middleBi)) {
    array<int, 10> bisArray;
    copy(row.begin(), row.end(), bisArray.begin());
    co_yield bisArray;
  }

  int biOffset = 1;
  while (biOffset <= (nextBiMax - nextBiMin + 1) / 2) {
    int nextBi = middleBi - biOffset;
    if (nextBi >= nextBiMin) {
      for (const auto &row :
           genBis(N, selfInt, genus, a, append(currentRow, nextBi),
                  currSelfIntSum + nextBi * nextBi,
                  currAdjSum + (nextBi * nextBi - nextBi))) {
        array<int, 10> bisArray;
        copy(row.begin(), row.end(), bisArray.begin());
        co_yield bisArray;
      }
    }
    nextBi = middleBi + biOffset;
    if (nextBi <= nextBiMax) {
      for (const auto &row :
           genBis(N, selfInt, genus, a, append(currentRow, nextBi),
                  currSelfIntSum + nextBi * nextBi,
                  currAdjSum + (nextBi * nextBi - nextBi))) {
        array<int, 10> bisArray;
        copy(row.begin(), row.end(), bisArray.begin());
        co_yield bisArray;
      }
    }

    biOffset++;
  }
}

array<int, 11> arrayConcat(const array<int, 1> &a, const array<int, 10> &bis) {
  array<int, 11> Fk;
  copy(a.begin(), a.end(), Fk.begin());
  copy(bis.begin(), bis.end(), Fk.begin() + 1);
  return Fk;
}

generator<array<int, 11>> rowCands(int N, int selfInt, int genus, int aMin,
                                   int aMax) {

  int aLowerBound = 0;
  if (selfInt < -2)
    aLowerBound = ceil((float)(selfInt + 1) / 2);

  int startA = max(aLowerBound, aMin);
  for (int a = startA; a <= aMax; a++) {
    generator<array<int, 10>> biVecs = genBis(N, selfInt, genus, a);
    array<int, 1> a_array = {a};
    for (const auto &bis : biVecs) {
      // TO DO, Array concat
      co_yield arrayConcat(a_array, bis);
    }
  }
}

// Helper: get sorted elements from map counter
vector<int> counterElements(const map<int, int> &ctr) {
  vector<int> elems;
  for (const auto &[val, count] : ctr) {
    for (int i = 0; i < count; i++)
      elems.push_back(val);
  }
  sort(elems.begin(), elems.end());
  return elems;
}

// Helper: dot product of a row of prevBiMat with a vector
int dotProduct(const vector<vector<int>> &mat, int row,
               const vector<int> &vec) {
  int sum = 0;
  for (int j = 0; j < vec.size(); j++)
    sum += mat[row][j] * vec[j];
  return sum;
}

// Helper: sum of a column subset of a matrix row
int sumCols(const vector<vector<int>> &mat, int row,
            const vector<int> &colInds) {
  int sum = 0;
  for (int col : colInds)
    sum += mat[row][col];
  return sum;
}

// Helper: delete columns from matrix
vector<vector<int>> deleteCols(const vector<vector<int>> &mat,
                               const vector<int> &colInds) {
  vector<vector<int>> result;
  for (const auto &row : mat) {
    vector<int> newRow;
    for (int j = 0; j < row.size(); j++) {
      if (find(colInds.begin(), colInds.end(), j) == colInds.end())
        newRow.push_back(row[j]);
    }
    result.push_back(newRow);
  }
  return result;
}

// Helper: generate combinations of r indices from 0..n-1
generator<vector<int>> combinations(int n, int r) {
  vector<int> indices(r);
  iota(indices.begin(), indices.end(), 0);
  co_yield indices;
  while (true) {
    int i = r - 1;
    while (i >= 0 && indices[i] == i + n - r)
      i--;
    if (i < 0)
      co_return;
    indices[i]++;
    for (int j = i + 1; j < r; j++)
      indices[j] = indices[j - 1] + 1;
    co_yield indices;
  }
}

generator<vector<int>> validBiPerms(vector<int> goalDots,
                                    vector<vector<int>> prevBiMat,
                                    map<int, int> bCtr) {
  // Sort prevBiMat rows
  vector<vector<int>> sortedPrevBis = prevBiMat;
  for (auto &row : sortedPrevBis)
    sort(row.begin(), row.end());

  vector<int> sortedNewBis = counterElements(bCtr); // already sorted
  vector<int> flippedNewBis(sortedNewBis.rbegin(), sortedNewBis.rend());

  int numRows = prevBiMat.size();

  // Compute maxDots and minDots
  vector<int> maxDots(numRows), minDots(numRows);
  for (int i = 0; i < numRows; i++) {
    maxDots[i] = dotProduct(sortedPrevBis, i, sortedNewBis);
    minDots[i] = dotProduct(sortedPrevBis, i, flippedNewBis);
  }

  // Check bounds
  for (int i = 0; i < numRows; i++) {
    if (goalDots[i] < minDots[i] || goalDots[i] > maxDots[i])
      co_return;
  }

  // Base case: one unique element
  if (bCtr.size() == 1) {
    co_yield counterElements(bCtr);
    co_return;
  }

  int numCols = sortedNewBis.size();

  // Find element with largest absolute value
  int biToAssign =
      max_element(bCtr.begin(), bCtr.end(), [](const auto &a, const auto &b) {
        return abs(a.first) < abs(b.first);
      })->first;
  int numCopies = bCtr[biToAssign];
  bCtr.erase(biToAssign);

  for (const auto &colInds : combinations(numCols, numCopies)) {
    // Compute dot contributions and new goal dots
    vector<int> newGoalDots(numRows);
    for (int i = 0; i < numRows; i++)
      newGoalDots[i] =
          goalDots[i] - biToAssign * sumCols(prevBiMat, i, colInds);

    vector<vector<int>> newPrevBiMat = deleteCols(prevBiMat, colInds);

    for (auto subPerm : validBiPerms(newGoalDots, newPrevBiMat, bCtr)) {
      vector<int> fullPerm = subPerm;
      // Insert left to right to restore original order
      for (int i = 0; i < colInds.size(); i++)
        fullPerm.insert(fullPerm.begin() + colInds[i], biToAssign);
      co_yield fullPerm;
    }
  }

  bCtr[biToAssign] = numCopies; // restore counter
}

namespace operations_research {
namespace sat {
vector<array<int, 11>> solveNextRow(const vector<array<int, 11>> &subMat,
                                    vector<vector<int>> &interNums, int selfInt,
                                    int genus, int aMax) {

  int numPrevRows = subMat.size();
  int N = 10;

  CpModelBuilder model;

  const IntVar a = model.NewIntVar(Domain(1, aMax)).WithName("a");
  int biUB = (int)min((double)(3 * aMax - selfInt + 2 * genus - 2),
                      trunc(sqrt((double)(aMax * aMax - selfInt))));
  vector<IntVar> b;
  for (int i = 0; i < N; i++) {
    b.push_back(model.NewIntVar(Domain(0, biUB)).WithName("b_" + to_string(i)));
  }

  for (int rowInd = 0; rowInd < numPrevRows; rowInd++) {
    LinearExpr expr = LinearExpr::Term(a, subMat[rowInd][0]);
    for (int i = 0; i < N; i++) {
      expr += LinearExpr::Term(b[i], subMat[rowInd][i + 1]);
    }
    model.AddEquality(expr, interNums[0][rowInd]);
  }

  model.AddEquality(LinearExpr::Sum(b),
                    LinearExpr::Term(a, 3) - selfInt + 2 * genus - 2);

  vector<IntVar> bisSquared = {};
  for (int i = 0; i < N; i++) {
    IntVar bSquared = model.NewIntVar(Domain(0, biUB * biUB))
                          .WithName("b_" + to_string(i) + "_squared");
    model.AddMultiplicationEquality(bSquared, {b[i], b[i]});
    bisSquared.push_back(bSquared);
  }

  IntVar aSquared =
      model.NewIntVar(Domain(1, aMax * aMax)).WithName("a_squared");
  model.AddMultiplicationEquality(aSquared, {a, a});
  model.AddEquality(LinearExpr::Sum(bisSquared),
                    LinearExpr::Term(aSquared, 1) - selfInt);

  vector<int> sameAsRightNeighbor = {};
  for (int i = 0; i < N - 1; i++) {
    // TO DO: Conditional involving slicing subMats
    if (sliceMat(subMat, 0, -1, i + 1, i + 2) ==
        sliceMat(subMat, 0, -1, i + 2, i + 3)) {
      model.AddGreaterOrEqual(b[i], b[i + 1]);
      sameAsRightNeighbor.push_back(i);
    }
  }

  SatParameters parameters;
  parameters.set_enumerate_all_solutions(true);

  vector<array<int, 11>> sols = {};

  Model orModel;
  orModel.Add(NewSatParameters(parameters));
  orModel.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse &r) {
    array<int, 11> sol;
    sol[0] = SolutionIntegerValue(r, a);
    for (int i = 0; i < N; i++) {
      sol[i + 1] = -SolutionIntegerValue(r, b[i]);
    }
    sols.push_back(sol);
  }));

  SolveCpModel(model.Build(), &orModel);
  /* Debug Tools
std::cout << model.Build().DebugString() << std::endl << flush;
std::cout << model.Build().variables_size() << std::endl << flush;
std::cout << model.Build().constraints_size() << std::endl << flush;*/
  // a <= 0 case

  /*
   int aMin = 0;
   if (selfInt < -2)
     aMin = ceil((double)(selfInt + 1) / 2);
   for (int currA = aMin; currA < 1; currA++) {
     int numOnes = 2 * currA - selfInt - 1;
     if ((0 <= numOnes) && (numOnes <= N - 1)) {
       vector<int> goalBDots;
       vector<vector<int>> firstCol = sliceMat(subMat, 0, -1, 0, 1);
       for (int i = 0; i < interNums[0].size(); i++) {
         if (i < firstCol[0].size()) {
           goalBDots.push_back(firstCol[0][i] * currA - interNums[0][i]);
         } else {
           goalBDots.push_back(-interNums[0][i]);
         }
       }
       map<int, int> negBCounter;
       negBCounter[-1] = numOnes;
       negBCounter[0] = N - 1 - numOnes;
       negBCounter[-currA + 1] = 1;

       for (const auto &perm : validBiPerms(
                goalBDots, sliceMat(subMat, 0, -1, 1, -1), negBCounter)) {
         bool discard = false;
         for (const auto &i : sameAsRightNeighbor) {
           if (perm[i] < perm[i + 1]) {
             discard = true;
             break;
           }
         }
         if (!discard) {
           array<int, 1> currA_array = {currA};
           array<int, 11> sol;
           sol[0] = currA;
           copy(perm.begin(), perm.end(), sol.begin() + 1);
           // sols.push_back(sol);
         }
       }
     }
}*/

  return sols;
}
} // namespace sat
} // namespace operations_research

int main(int argc, char *argv[]) {

  int jobInd = 0;
  int rowInd = 0;
  if (argc != 0) {
    stringstream convert{argv[2]};
    if (!(convert >> rowInd)) {
      rowInd = 0;
    }

    if (rowInd != 0) {
      stringstream convert{argv[3]};

      if (!(convert >> jobInd)) {
        jobInd = 0;
      }
    }
  }

  const int N = 10;
  vector<array<int, 11>> intMat = {
      {-3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, -3, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0}, {1, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1},
      {1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1}};

  array<int, 11> aMaxes = {11, 36, 18, 12, 9, 6, 12, 6, 35, 4, 4};
  array<int, 11> genera = {0};

  array<int, 11> keyArray;

  for (int i = 0; i < 11; i++) {
    keyArray[i] = aMaxes[i] * aMaxes[i] - intMat[i][i];
  }

  array<int, 11> indices;
  iota(indices.begin(), indices.end(), 0);
  sort(indices.begin(), indices.end(),
       [&keyArray](int i, int j) { return keyArray[i] < keyArray[j]; });

  vector<array<int, 11>> intMatReOrd(11);
  array<int, 11> aMaxesReOrd;
  array<int, 11> generaReOrd;

  for (int i = 0; i < 11; i++) {
    aMaxesReOrd[i] = aMaxes[indices[i]];
    generaReOrd[i] = genera[indices[i]];
    for (int j = 0; j < 11; j++) {
      intMatReOrd[i][j] = intMat[indices[i]][indices[j]];
    }
  }

  // if (rowInd == -1) {
  generator<array<int, 11>> firstRows =
      rowCands(N, intMatReOrd[0][0], generaReOrd[0], -1000000, aMaxesReOrd[0]);

  rowInd = 0;
  vector<vector<array<int, 11>>> newSols;
  vector<vector<int>> intMatSlice =
      sliceMat(intMatReOrd, rowInd + 1, rowInd + 2, 0, rowInd + 1);
  int count = 0;
  for (const auto &firstRow : firstRows) {
    count++;
    vector<array<int, 11>> firstRowVect = {firstRow};
    vector<array<int, 11>> nextRows = operations_research::sat::solveNextRow(
        firstRowVect, intMatSlice, intMatReOrd[rowInd + 1][rowInd + 1],
        generaReOrd[rowInd + 1], aMaxesReOrd[rowInd + 1]);

    sort(nextRows.begin(), nextRows.end());
    nextRows.erase(unique(nextRows.begin(), nextRows.end()), nextRows.end());
    for (const auto &newRow : nextRows) {
      vector<array<int, 11>> newSol;
      newSol.insert(newSol.end(), firstRowVect.begin(), firstRowVect.end());
      newSol.push_back(newRow);
      newSols.push_back(newSol);
    }
  }
  cout << count << "\n";

  rowInd++;
  for (const auto &sol : newSols) {
    displayMatrix(sol);
    cout << "\n";
  }

  while (rowInd < 11) {
    cout << newSols.size() << "\n" << flush;
    vector<vector<array<int, 11>>> oldSols = newSols;
    newSols.clear();
    vector<vector<int>> intMatSlice =
        sliceMat(intMatReOrd, rowInd + 1, rowInd + 2, 0, rowInd + 2);
    for (const auto &subMat : oldSols) {
      vector<array<int, 11>> nextRows = operations_research::sat::solveNextRow(
          subMat, intMatSlice, intMatReOrd[rowInd + 1][rowInd + 1],
          generaReOrd[rowInd + 1], aMaxesReOrd[rowInd + 1]);
      sort(nextRows.begin(), nextRows.end());
      nextRows.erase(unique(nextRows.begin(), nextRows.end()), nextRows.end());
      for (const auto &newRow : nextRows) {
        vector<array<int, 11>> newSol;
        newSol.insert(newSol.end(), subMat.begin(), subMat.end());
        newSol.push_back(newRow);
        newSols.push_back(newSol);
      }
    }
    rowInd++;
  }

  cout << newSols.size();
  // TO DO : Figure out how to slice matrices (<vector<array<int, 11>> to
  // vector<array<int, 11>> nextRows = solveNextRow(subMat, intMat,
  // intMat[1, 1], 0, aMaxes[1]);
}
