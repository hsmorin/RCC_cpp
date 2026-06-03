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
#include "tree.hh"
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

/* IMPORTANT WARNING! This program computes and saves the rows of the indicator
matrices in a different order to the input to improve the efficiency of the
computation. Remember to change the order of the rows of the matrices back
before using them!!!*/
// Before use update the problem data choose which main function to use

// Problem Data (update this)
constexpr int N = 22; // N as in n - 1
constexpr int M = 23; // n (number of rows of the intersection matrix)
const vector<array<int, M>> intMat = {
    {-12, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
    {0, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 1, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -6, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -1}}; // Intersection Matrix
const array<int, M> aMaxes = {
    155, 3,   3, 3, 3,  3,  3,   3,  3,  3,  3, 95,
    480, 239, 3, 3, 96, 46, 160, 80, 53, 40, 32}; // a_maxes of F_k
const array<int, M> genera = {0};                 // Genera of F_k
// End of Problem Data (Only Change which Main function to use after this)

void displayMatrix(const vector<vector<int>> &matrix) {
  for (const auto &row : matrix) {
    for (int val : row) {
      cout << val << " ";
    }
    cout << "\n";
  }
}

void displayMatrix(const vector<array<int, M>> &matrix) {
  for (const auto &row : matrix) {
    for (int val : row) {
      cout << val << " ";
    }
    cout << "\n";
  }
}

void displayArray(const array<int, M> &array) {
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

// Save Matrix Tree (To Store Data)
void saveTree(const tree<array<int, M>> &matTree, const string &filename) {
  ofstream out(filename, ios::binary);
  for (auto it = matTree.begin(); it != matTree.end(); ++it) {
    if (!matTree.is_valid(it))
      continue;
    int d = matTree.depth(it);
    if (d == 0)
      continue; // skip sentinel root
    //
    out.write(reinterpret_cast<const char *>(&d), sizeof(int));
    out.write(reinterpret_cast<const char *>((*it).data()), M * sizeof(int));
  }
}

// Load Matrix Tree (To Read Data)
tree<array<int, M>> loadTree(const string &filename) {
  tree<array<int, M>> matTree;
  array<int, M> sentinel = {};
  auto root = matTree.insert(matTree.begin(), sentinel);

  vector<tree<array<int, M>>::iterator> lastAtDepth;
  lastAtDepth.push_back(root); // depth 0 = sentinel

  ifstream in(filename, ios::binary);
  int depth;
  while (in.read(reinterpret_cast<char *>(&depth), sizeof(int))) {
    array<int, M> row;
    in.read(reinterpret_cast<char *>(row.data()), M * sizeof(int));

    while ((int)lastAtDepth.size() <= depth)
      lastAtDepth.push_back(root);

    auto parent = lastAtDepth[depth - 1];
    bool found = false;
    for (auto sib = matTree.begin(parent); sib != matTree.end(parent); ++sib) {
      if (*sib == row) {
        lastAtDepth[depth] = sib;
        found = true;
        break;
      }
    }

    if (!found) {
      auto node = matTree.append_child(parent, row);
      lastAtDepth[depth] = node;
    }
  }
  return matTree;
}

void pruneShallowLeaves(tree<array<int, M>> &matTree, int targetDepth) {
  bool changed = true;
  // Repeat until no more shallow leaves exist
  // (pruning a leaf may expose its parent as a new shallow leaf)
  while (changed) {
    changed = false;
    vector<tree<array<int, M>>::iterator> toErase;
    for (auto it = matTree.begin_leaf(); it != matTree.end_leaf(); ++it) {
      if (matTree.depth(it) < targetDepth) {
        toErase.push_back(it);
      }
    }
    for (auto it = toErase.rbegin(); it != toErase.rend(); ++it) {
      matTree.erase(*it);
      changed = true;
    }
  }
}

// Reconstructs Matrix from Leaf Itterator
vector<array<int, M>> reconstructMatrix(const tree<array<int, M>> &matTree,
                                        tree<array<int, M>>::iterator leaf) {

  vector<array<int, M>> matrix;
  auto it = leaf;

  while (matTree.is_valid(it) && it != matTree.begin()) {
    matrix.push_back(*it);
    it = matTree.parent(it);
  }

  reverse(matrix.begin(), matrix.end());
  return matrix;
}

void saveSubtree(const tree<array<int, M>> &matTree,
                 const vector<tree<array<int, M>>::iterator> &leaves, int start,
                 int end, const string &filename) {
  ofstream out(filename, ios::binary);

  set<void *> written;
  for (int i = start; i < end; i++) {
    // Reconstruct path for each leaf
    vector<tree<array<int, M>>::iterator> path;
    auto it = leaves[i];
    while (matTree.depth(it) > 0) {
      path.push_back(it);
      it = matTree.parent(it);
    }
    reverse(path.begin(), path.end());

    for (auto &nodeIt : path) {
      void *ptr = nodeIt.node;
      if (written.find(ptr) == written.end()) {
        int depth = matTree.depth(nodeIt);
        out.write(reinterpret_cast<const char *>(&depth), sizeof(int));
        out.write(reinterpret_cast<const char *>((*nodeIt).data()),
                  M * sizeof(int));
        written.insert(ptr);
      }
    }
  }
}

void saveTreeChunks(const tree<array<int, M>> &matTree, string filePrefix,
                    int rowInd, int numSplits) {
  // Collect all leaves
  vector<tree<array<int, M>>::iterator> leaves;
  for (auto it = matTree.begin_leaf(); it != matTree.end_leaf(); ++it) {
    if (matTree.depth(it) > 0)
      leaves.push_back(it);
  }

  int total = leaves.size();
  int chunkSize = (total + numSplits - 1) / numSplits; // ceiling division

  for (int i = 0; i < numSplits; i++) {
    int start = i * chunkSize;
    int end = min(start + chunkSize, total);
    if (start >= total)
      break; // fewer leaves than chunks

    string filename =
        filePrefix + to_string(i) + "_" + to_string(rowInd) + ".tr";
    saveSubtree(matTree, leaves, start, end, filename);
    cout << "Saved chunk " << i << " (" << (end - start) << " leaves) to "
         << filename << "\n";
  }
}

// Slices a Generic Matrix Type
template <typename MatrixType>
vector<vector<int>> sliceMat(const MatrixType &matrix, int rowStart = 0,
                             int rowEnd = -1, int colnStart = 0,
                             int colnEnd = -1) {
  vector<vector<int>> result;
  if (matrix.empty())
    return result;
  if (rowEnd == -1) {
    rowEnd = matrix.size();
  }
  if (colnEnd == -1) {
    colnEnd = matrix[0].size();
  }

  result.reserve(rowEnd - rowStart);

  for (int rowInd = rowStart; rowInd < rowEnd; rowInd++) {
    const auto &currentRow = matrix[rowInd];
    vector<int> arraySlice(currentRow.begin() + colnStart,
                           currentRow.begin() + colnEnd);
    result.push_back(std::move(arraySlice));
  }
  return result;
}

// Overloading of sliceMat
vector<vector<int>> sliceMat(const vector<array<int, M>> &matrix,
                             int rowStart = 0, int rowEnd = -1,
                             int colnStart = 0, int colnEnd = -1) {
  vector<vector<int>> result;
  if (rowEnd == -1) {
    rowEnd = matrix.size();
  }
  if (colnEnd == -1) {
    colnEnd = M;
  }
  for (int rowInd = rowStart; rowInd < rowEnd; rowInd++) {
    array<int, M> currentRow = matrix[rowInd];
    vector<int> arraySlice = {};
    for (int colnInd = colnStart; colnInd < colnEnd; colnInd++) {
      arraySlice.push_back(currentRow[colnInd]);
    }
    result.push_back(arraySlice);
  }
  return result;
}

// Helper Function for Row Cands
generator<array<int, N>> genBis(int selfInt, int genus, int a,
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
        array<int, N> bisArray;
        copy(bis.begin(), bis.end(), bisArray.begin());
        co_yield bisArray;
        co_return;
      }
    }
    co_return;
  }

  int currLen = currentRow.size();
  if (currLen == N) {
    array<int, N> bisArray;
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
  for (const auto &row : genBis(selfInt, genus, a, append(currentRow, middleBi),
                                currSelfIntSum + middleBi * middleBi,
                                currAdjSum + middleBi * middleBi - middleBi)) {
    array<int, N> bisArray;
    copy(row.begin(), row.end(), bisArray.begin());
    co_yield bisArray;
  }

  int biOffset = 1;
  while (biOffset <= (nextBiMax - nextBiMin + 1) / 2) {
    int nextBi = middleBi - biOffset;
    if (nextBi >= nextBiMin) {
      for (const auto &row :
           genBis(selfInt, genus, a, append(currentRow, nextBi),
                  currSelfIntSum + nextBi * nextBi,
                  currAdjSum + (nextBi * nextBi - nextBi))) {
        array<int, N> bisArray;
        copy(row.begin(), row.end(), bisArray.begin());
        co_yield bisArray;
      }
    }
    nextBi = middleBi + biOffset;
    if (nextBi <= nextBiMax) {
      for (const auto &row :
           genBis(selfInt, genus, a, append(currentRow, nextBi),
                  currSelfIntSum + nextBi * nextBi,
                  currAdjSum + (nextBi * nextBi - nextBi))) {
        array<int, N> bisArray;
        copy(row.begin(), row.end(), bisArray.begin());
        co_yield bisArray;
      }
    }

    biOffset++;
  }
}

// Helper Function for rowCands
array<int, M> arrayConcat(const array<int, 1> &a, const array<int, N> &bis) {
  array<int, M> Fk;
  copy(a.begin(), a.end(), Fk.begin());
  copy(bis.begin(), bis.end(), Fk.begin() + 1);
  return Fk;
}

// Computes the possible first rows of the matrices
generator<array<int, M>> rowCands(int selfInt, int genus, int aMin, int aMax) {

  int aLowerBound = 0;
  if (selfInt < -2)
    aLowerBound = ceil((float)(selfInt + 1) / 2);

  int startA = max(aLowerBound, aMin);
  for (int a = startA; a <= aMax; a++) {
    generator<array<int, N>> biVecs = genBis(selfInt, genus, a);
    array<int, 1> a_array = {a};
    for (const auto &bis : biVecs) {
      co_yield arrayConcat(a_array, bis);
    }
  }
}

// Helper for solveNextRow: get sorted elements from map counter
vector<int> counterElements(const map<int, int> &ctr) {
  vector<int> elems;
  for (const auto &[val, count] : ctr) {
    for (int i = 0; i < count; i++)
      elems.push_back(val);
  }
  sort(elems.begin(), elems.end());
  return elems;
}

// Helper for solveNextRow: dot product of a row of prevBiMat with a vector
int dotProduct(const vector<vector<int>> &mat, int row,
               const vector<int> &vec) {
  int sum = 0;
  for (int j = 0; j < vec.size(); j++)
    sum += mat[row][j] * vec[j];
  return sum;
}

// Helper for solveNextRow: sum of a column subset of a matrix row
int sumCols(const vector<vector<int>> &mat, int row,
            const vector<int> &colInds) {
  int sum = 0;
  for (int col : colInds)
    sum += mat[row][col];
  return sum;
}

// Helper for SolveNextRow: delete columns from matrix
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

// Helper for SolveNextRow: generate combinations of r indices from 0..n-1
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

// Helper for SolveNextRow
generator<vector<int>> validBiPerms(const vector<int> &goalDots,
                                    const vector<vector<int>> &prevBiMat,
                                    map<int, int> &bCtr) {

  int numRows = prevBiMat.size();

  // Sort prevBiMat rows
  vector<vector<int>> sortedPrevBis = prevBiMat;
  for (auto &row : sortedPrevBis)
    sort(row.begin(), row.end());

  vector<int> sortedNewBis = counterElements(bCtr); // already sorted
  vector<int> flippedNewBis(sortedNewBis.rbegin(), sortedNewBis.rend());

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
  auto maxIt =
      max_element(bCtr.begin(), bCtr.end(), [](const auto &a, const auto &b) {
        return abs(a.first) < abs(b.first);
      });
  int biToAssign = maxIt->first;
  int numCopies = maxIt->second;

  // Backtracking: remove for recursion, restore afterwords
  bCtr.erase(maxIt);

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

// Finds the next possible rows given a subMatrix
// Call via operations_research::sat::solveNextRow()
vector<array<int, M>> solveNextRow(const vector<array<int, M>> &subMat,
                                   vector<vector<int>> &interNums, int selfInt,
                                   int genus, int aMax) {

  int numPrevRows = subMat.size();

  // a > 0 case
  //  Solves the problem using CP_SAT with constraints
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

  vector<array<int, M>> sols = {};

  Model orModel;
  orModel.Add(NewSatParameters(parameters));
  orModel.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse &r) {
    array<int, M> sol;
    sol[0] = SolutionIntegerValue(r, a);
    for (int i = 0; i < N; i++) {
      sol[i + 1] = -SolutionIntegerValue(r, b[i]);
    }
    sols.push_back(sol);
  }));

  SolveCpModel(model.Build(), &orModel);

  // a <= 0 case
  int aMin = 0;
  if (selfInt < -2)
    aMin = ceil((double)(selfInt + 1) / 2);
  for (int currA = aMin; currA < 1; currA++) {
    int numOnes = 2 * currA - selfInt - 1;
    if (0 <= numOnes && numOnes <= (N - 1)) {
      vector<int> goalBDots(numPrevRows);
      for (int i = 0; i < numPrevRows; i++) {
        goalBDots[i] = subMat[i][0] * currA - interNums[0][i];
      }
      map<int, int> negBCounter;
      if (numOnes > 0)
        negBCounter[-1] = numOnes;
      if (N - 1 - numOnes > 0)
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
          array<int, M> sol;
          sol[0] = currA;
          copy(perm.begin(), perm.end(), sol.begin() + 1);
          sols.push_back(sol);
        }
      }
    }
  }

  return sols;
}
} // namespace sat
} // namespace operations_research

// Main function using trees
// Computes everthing at once
int mainTr() {

  int rowInd = 0;

  // Sort matrix rows by key
  array<int, M> keyArray;

  for (int i = 0; i < M; i++) {
    keyArray[i] = aMaxes[i] * aMaxes[i] - intMat[i][i];
  }

  array<int, M> indices;
  iota(indices.begin(), indices.end(), 0);
  sort(indices.begin(), indices.end(),
       [&keyArray](int i, int j) { return keyArray[i] < keyArray[j]; });

  vector<array<int, M>> intMatReOrd(M);
  array<int, M> aMaxesReOrd;
  array<int, M> generaReOrd;

  for (int i = 0; i < M; i++) {
    aMaxesReOrd[i] = aMaxes[indices[i]];
    generaReOrd[i] = genera[indices[i]];
    for (int j = 0; j < M; j++) {
      intMatReOrd[i][j] = intMat[indices[i]][indices[j]];
    }
  }

  tree<array<int, M>> matTree;
  // Empty nodes are the zero vector by default
  array<int, M> firstNode;
  auto root = matTree.insert(matTree.begin(), firstNode);

  // if (rowInd == -1) {
  generator<array<int, M>> firstRows =
      rowCands(intMatReOrd[0][0], generaReOrd[0], -1000000, aMaxesReOrd[0]);

  rowInd = 0;
  vector<vector<int>> intMatSlice =
      sliceMat(intMatReOrd, rowInd + 1, rowInd + 2, 0, rowInd + 1);
  int count = 0;
  int countFirstRows = 0;

  for (const auto &firstRow : firstRows) {
    vector<array<int, M>> firstRowVect = {firstRow};
    vector<array<int, M>> nextRows = operations_research::sat::solveNextRow(
        firstRowVect, intMatSlice, intMatReOrd[rowInd + 1][rowInd + 1],
        generaReOrd[rowInd + 1], aMaxesReOrd[rowInd + 1]);

    auto firstRowNode = matTree.append_child(root, firstRow);
    countFirstRows++;
    for (const auto &newRow : nextRows) {
      count++;
      matTree.append_child(firstRowNode, newRow);
    }
  }

  rowInd++;
  cout << countFirstRows << "\n" << flush;
  while (rowInd < N) {
    cout << count << "\n" << flush;
    count = 0;
    vector<tree<array<int, M>>::iterator> leaves;
    for (auto it = matTree.begin_leaf(); it != matTree.end_leaf(); it++) {
      leaves.push_back(it);
    }

    vector<vector<int>> intMatSlice =
        sliceMat(intMatReOrd, rowInd + 1, rowInd + 2, 0, rowInd + 2);

    for (const auto &leafIt : leaves) {
      vector<array<int, M>> subMat = reconstructMatrix(matTree, leafIt);

      if (subMat.size() != rowInd + 1) {
        matTree.erase(leafIt);
        continue;
      }

      vector<array<int, M>> nextRows = operations_research::sat::solveNextRow(
          subMat, intMatSlice, intMatReOrd[rowInd + 1][rowInd + 1],
          generaReOrd[rowInd + 1], aMaxesReOrd[rowInd + 1]);

      if (nextRows.empty()) {
        matTree.erase(leafIt);
        continue;
      }

      for (const auto &newRow : nextRows) {
        matTree.append_child(leafIt, newRow);
        count++;
      }
    }
    rowInd++;

    pruneShallowLeaves(matTree, rowInd);
    saveTree(matTree, "matTree_" + to_string(rowInd) + ".tr");
    matTree = loadTree("matTree_" + to_string(rowInd) + ".tr");
  }

  cout << count << "\n" << flush;
  // TO DO : Figure out how to slice matrices (<vector<array<int, 11>> to
  // vector<array<int, 11>> nextRows = solveNextRow(subMat, intMat,
  // intMat[1, 1], 0, aMaxes[1]);
  //
  return 0;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " filePrefix rowIndex NumberOfSplits\n";
    return 1;
  }

  string filePrefix = argv[1];
  int rowInd;
  int numSplits;

  stringstream convert2{argv[2]};
  stringstream convert3{argv[3]};

  if (!(convert2 >> rowInd)) {
    cerr << "Usage: " << argv[0] << " filePrefix rowIndex NumberOfSplits\n";
    return 1;
  }

  if (!(convert3 >> numSplits)) {
    cerr << "Usage: " << argv[0] << " filePrefix rowIndex NumberOfSplits\n";
    return 1;
  }

  string inputFilename = filePrefix + to_string(rowInd) + ".tr";
  auto matTree = loadTree(inputFilename);
  saveTreeChunks(matTree, filePrefix, rowInd, numSplits);
  return 0;
}

int mainTR() {

  int count = 0;
  auto matTree = loadTree("matTree_4.tr");
  vector<tree<array<int, M>>::iterator> leaves;
  for (auto it = matTree.begin_leaf(); it != matTree.end_leaf(); it++) {
    leaves.push_back(it);
  }

  for (int jobInd = 0; jobInd < 10; jobInd++) {
    auto matSubTree = loadTree("matTree_" + to_string(jobInd) + "_4.tr");

    vector<tree<array<int, M>>::iterator> subLeaves;
    for (auto it = matSubTree.begin_leaf(); it != matSubTree.end_leaf(); it++) {
      subLeaves.push_back(it);
    }

    for (auto leaf : subLeaves) {
      auto mat1 = reconstructMatrix(matSubTree, leaf);

      auto mat2 = reconstructMatrix(matTree, leaves[count]);

      if (mat1 == mat2) {
        cout << ":)";
      } else {
        displayMatrix(mat1);
        displayMatrix(mat2);
        cout << ":(";
        abort();
      }
      count += 1;
      cout << "\n";
    }
  }
  cout << "Count: " << count;
  return 0;
}
