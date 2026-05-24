#include "absl/base/log_severity.h"
#include "absl/log/globals.h"
#include "ortools/base/init_google.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
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

void displayMatrix(vector<array<int, 11>> &matrix) {
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
}

vector<int> append(const vector<int> &v, int x) {
  vector<int> result = v;
  result.push_back(x);
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
        // TO DO: Write Concat Statment
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
  for (int a = startA; a < trunc(aMax + 1); a++) {
    generator<array<int, 10>> biVecs = genBis(N, selfInt, genus, a);
    array<int, 1> a_array = {a};
    for (const auto &bis : biVecs) {
      // TO DO, Array concat
      co_yield arrayConcat(a_array, bis);
    }
  }
}

/* validBiPerms(goalDots, prevBiMat, bCtr) {

     }

namespace operations_research {
namespace sat {
vector<vector<array<int, 11>>>
solveNextRow(const vector<array<int, 11>> &subMat,
             vector<vector<int>> &interNums, int selfInt, int genus) {
  int numPrevRows = subMat.size();
  int N = 11;

  CpModelBuilder model;

  const intVar a = model.NewIntVar(1, aMax, "a");
  int biUB = min((double)(3 * aMax - selfInt + 2 * genus - 2),
                 trunc(sqrt((double)(aMax * aMax - selfInt))));
  vector<IntVar> b;
  for (int i = 0; i < N; i++) {
    b.push_back(model.NewIntVar(0, biUB, "b_" + to_string(i)));
  }

  for (int rowInd = 0; rowInd < numPrevRows; rowInd++) {
    int expr = trunc(subMat[rowInt][0]) * a;
    for (int i = 0; i < N; i++) {
      expr += trunc(subMat[rowInd][i + 1]);
    }
    model.Add(expr == trunc(interNums[rowInd]));
  }

  model.Add(accumulate(b.begin(), b.end(), 0) ==
            3 * a - selfInt + 2 * genus - 2);

  vector<intVar> bisSquared = {};
  for (i = 0; i < N; i++) {
    intVar bSquared = model.NewIntVar(0, biUB * *2, format("b_{}_squared", i));
    model.AddMultiplicationEquality({bSquared[i], bSquared[i]});
    bisSquared.push_back(bSquared);
  }

  intVar aSquared = model.NewIntVar(1, aMax * *2, "a_squared");
  model.AddMultiplicationEquality(aSquared, {a, a});
  model.Add(accumulate(bisSquared.begin(), bisSquared.end(), 0) ==
            aSquared - selfInt);

  vector<int> sameAsRightNeighbor = {};
  for (i = 0; i < N - 1; i++) {
    // TO DO: Conditional involving slicing subMats
    if (subMat[:][i + 1] == subMat[:][i + 2]) {
      model.Add(b[i] >= b[i + 1]);
      sameAsRightNeighbor.push_back(i);
    }
  }

  Model model;
  // TO DO: Fix Solutions Initialization
  solver = cp_model.CpSolver();
  SatParameters parameters;
  solver.parameters.enumerate_all_solutions(true);
  model.Add(NewSatParameters(parameters));

  vector<array<int, 11>> sols = {};

  class SolutionCollector : public cp_model.CpSolverSolutionCallback {
  public:
    void on_solution_callback()
        :

          array<int, 1> a_array = Value(a);
    array<int, 10> bis_array;
    for (int i = 0; i < 10; i++) {
      bis_array[i] = -1 * Value(b[i])
    }
    sols.push_back(append(a_array, bis_array))
  }

  int aMin = 0;
  if selfInt
    < -2 aMin = ceil((double)(selfInt + 1) / 2);
  for (currA = aMin; currA < 1, currA++) {
    numOnes = 2 * currA - selfInt - 1;
    if ((0 <= numOnes) && (numOnes <= N - 1)) {
      goalBDots = subMat[:, 0] * currA - interNums;
      // TO DO: Counter Stuff
      for (const auto &perm :
           validBiPerms(goalBDots, subMat[:] [1:], negBCounter)) {
        bool discard = false;
        for (const auto &i : sameAsRightNeighbor) {
          if (perm[i] < perm[i + 1]) {
            discard = true;
            break;
          }
          if !(discard) {
            array<int, 1> currA_array = {currA};
            sols.push_back(append(currA_array, perm));
          }
        }
      }
    }
  }

  return sols;
}
} // namespace sat
} // namespace operations_research */

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

  int aMaxes[11] = {11, 36, 18, 12, 9, 6, 12, 6, 35, 4, 4};
  int genera[11] = {0};

  // if (rowInd == -1) {
  generator<array<int, 11>> firstRows =
      rowCands(N, intMat[0][0], genera[0], -1000000, aMaxes[0]);
  for (const auto &firstRow : firstRows) {
    displayArray(firstRow);
    cout << "\n";
    //}
  }

  // TO DO : Figure out how to slice matrices (<vector<array<int, 11>> to
  // <vector<vector<int>>)
  // vector<array<int, 11>> nextRows = solveNextRow(subMat, intMat, intMat[1,
  // 1], 0, aMaxes[1]);
}
