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
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace std;

/*
IMPORTANT WARNING! This program computes and saves the rows of the indicator
matrices in a different order to the input to improve the efficiency of the
computation. Remember to change the order of the rows of the matrices back
before using them.
*/

// Problem Data (update this)
// Degree 10
constexpr int N = 22;  // N as in n - 1
constexpr int M = 23;  // n, the number of rows of the intersection matrix
static_assert(M == N + 1, "This code assumes one a-coordinate plus N b-coordinates.");
static_assert(sizeof(int) == 4, "The binary .tr format assumes 32-bit int values.");

using Row = array<int, M>;
using BisRow = array<int, N>;
using DenseMatrix = vector<vector<int>>;
using RowMatrix = vector<Row>;
using MatrixTree = tree<Row>;

constexpr int kDefaultOuterWorkers = 8;
constexpr int kMaxOuterWorkers = 8;
constexpr int kCpSatWorkersPerSolve = 1;
constexpr uint32_t kTreeMagic = 0x31434352u;  // "RCC1" on little-endian machines.
constexpr uint32_t kTreeVersion = 1u;

const RowMatrix intMat = {
    Row{-12, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
    Row{0, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{1, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 1, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -3, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0},
    Row{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -6, 0, 0, 0, 0, 1},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1, 0},
    Row{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2, 1},
    Row{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -1}};

const Row aMaxes = Row{155, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 95, 480, 239, 3, 3, 96, 46, 160, 80, 53, 40, 32};
const Row genera = Row{};  // all zero for the current degree-10 data
// End of Problem Data

void displayMatrix(const DenseMatrix &matrix) {
    for (const auto &row : matrix) {
        for (int val : row) cout << val << ' ';
        cout << '\n';
    }
}

void displayMatrix(const RowMatrix &matrix) {
    for (const auto &row : matrix) {
        for (int val : row) cout << val << ' ';
        cout << '\n';
    }
}

void displayArray(const Row &arrayValues) {
    for (int val : arrayValues) cout << val << ' ';
    cout << '\n';
}

void displayVector(const vector<int> &arrayValues) {
    for (int val : arrayValues) cout << val << ' ';
    cout << '\n';
}

int parsePositiveEnv(const char *name) {
    const char *raw = getenv(name);
    if (raw == nullptr || *raw == '\0') return 0;

    char *end = nullptr;
    long value = strtol(raw, &end, 10);
    if (end == raw || *end != '\0' || value <= 0) return 0;
    return static_cast<int>(min<long>(value, kMaxOuterWorkers));
}

int requestedOuterWorkers() {
    if (int n = parsePositiveEnv("RCC_NUM_THREADS"); n > 0) return n;
    if (int n = parsePositiveEnv("OMP_NUM_THREADS"); n > 0) return n;
    return kDefaultOuterWorkers;
}

void writeOrThrow(ostream &out, const char *data, streamsize size, const string &what) {
    out.write(data, size);
    if (!out) throw runtime_error("failed while writing " + what);
}

bool readExact(istream &in, char *data, streamsize size) {
    in.read(data, size);
    return static_cast<bool>(in);
}

void insertTreeRecord(MatrixTree &matTree,
                      vector<MatrixTree::iterator> &lastAtDepth,
                      MatrixTree::iterator root,
                      int depth,
                      const Row &row) {
    if (depth <= 0 || depth > M) {
        throw runtime_error("invalid depth in .tr file: " + to_string(depth));
    }
    while (static_cast<int>(lastAtDepth.size()) <= depth) lastAtDepth.push_back(root);

    auto parent = lastAtDepth[depth - 1];
    if (!matTree.is_valid(parent)) {
        throw runtime_error("invalid parent while reconstructing .tr tree at depth " + to_string(depth));
    }

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

    // Do not let stale deeper-path iterators from a previous branch become parents.
    if (static_cast<int>(lastAtDepth.size()) > depth + 1) {
        lastAtDepth.resize(depth + 1);
    }
}

void readDepthRowRecords(istream &in,
                         MatrixTree &matTree,
                         vector<MatrixTree::iterator> &lastAtDepth,
                         MatrixTree::iterator root) {
    while (true) {
        int depth = 0;
        Row row{};

        in.read(reinterpret_cast<char *>(&depth), sizeof(depth));
        if (in.eof()) break;
        if (!in) throw runtime_error("truncated .tr file while reading a node depth");
        if (!readExact(in, reinterpret_cast<char *>(row.data()), M * sizeof(int))) {
            throw runtime_error("truncated .tr file while reading a matrix row");
        }
        insertTreeRecord(matTree, lastAtDepth, root, depth, row);
    }
}

// Save Matrix Tree (To Store Data). The format is:
//   uint32 magic "RCC1", uint32 version, then repeated (int depth, Row row) records.
void saveTree(const MatrixTree &matTree, const string &filename) {
    ofstream out(filename, ios::binary | ios::trunc);
    if (!out) throw runtime_error("could not open output file: " + filename);

    writeOrThrow(out, reinterpret_cast<const char *>(&kTreeMagic), sizeof(kTreeMagic), "tree magic");
    writeOrThrow(out, reinterpret_cast<const char *>(&kTreeVersion), sizeof(kTreeVersion), "tree version");

    for (auto it = matTree.begin(); it != matTree.end(); ++it) {
        if (!matTree.is_valid(it)) continue;
        const int depth = matTree.depth(it);
        if (depth == 0) continue;  // skip sentinel root

        writeOrThrow(out, reinterpret_cast<const char *>(&depth), sizeof(depth), "node depth");
        writeOrThrow(out, reinterpret_cast<const char *>((*it).data()), M * sizeof(int), "matrix row");
    }
}

// Load Matrix Tree (To Read Data). Accepts the new magic-header format and the
// older depth+row format. Row-only legacy files cannot be reconstructed because
// they do not contain parent/depth information.
MatrixTree loadTree(const string &filename) {
    ifstream in(filename, ios::binary);
    if (!in) throw runtime_error("could not open input file: " + filename);

    in.seekg(0, ios::end);
    const streamoff fileSize = in.tellg();
    in.seekg(0, ios::beg);

    MatrixTree matTree;
    Row sentinel{};
    auto root = matTree.insert(matTree.begin(), sentinel);
    vector<MatrixTree::iterator> lastAtDepth;
    lastAtDepth.push_back(root);  // depth 0 = sentinel root

    if (fileSize == 0) return matTree;

    uint32_t magic = 0;
    if (fileSize >= static_cast<streamoff>(2 * sizeof(uint32_t))) {
        in.read(reinterpret_cast<char *>(&magic), sizeof(magic));
        if (!in) throw runtime_error("failed while reading .tr magic");
    }

    if (magic == kTreeMagic) {
        uint32_t version = 0;
        if (!readExact(in, reinterpret_cast<char *>(&version), sizeof(version))) {
            throw runtime_error("truncated .tr file while reading version");
        }
        if (version != kTreeVersion) {
            throw runtime_error("unsupported .tr version: " + to_string(version));
        }
        readDepthRowRecords(in, matTree, lastAtDepth, root);
        return matTree;
    }

    // Legacy file without a magic header. It must be exactly a sequence of
    // (int depth, Row row) records.
    in.clear();
    in.seekg(0, ios::beg);
    const streamoff depthRowRecordSize = static_cast<streamoff>(sizeof(int) + M * sizeof(int));
    const streamoff rowOnlyRecordSize = static_cast<streamoff>(M * sizeof(int));
    if (fileSize % depthRowRecordSize != 0) {
        if (fileSize % rowOnlyRecordSize == 0) {
            throw runtime_error(
                "input appears to be a legacy row-only .tr file. It cannot be "
                "reconstructed because node depths were not saved; regenerate it "
                "with the fixed saveTree() format.");
        }
        throw runtime_error("input .tr file has an invalid byte length");
    }

    readDepthRowRecords(in, matTree, lastAtDepth, root);
    return matTree;
}

void pruneShallowLeaves(MatrixTree &matTree, int targetDepth) {
    bool changed = true;
    while (changed) {
        changed = false;
        vector<MatrixTree::iterator> toErase;
        for (auto it = matTree.begin_leaf(); it != matTree.end_leaf(); ++it) {
            if (matTree.depth(it) < targetDepth) toErase.push_back(it);
        }
        for (auto it = toErase.rbegin(); it != toErase.rend(); ++it) {
            matTree.erase(*it);
            changed = true;
        }
    }
}

// Reconstructs Matrix from Leaf Iterator.
RowMatrix reconstructMatrix(const MatrixTree &matTree, MatrixTree::iterator leaf) {
    RowMatrix matrix;
    auto it = leaf;
    while (matTree.is_valid(it) && it != matTree.begin()) {
        matrix.push_back(*it);
        it = matTree.parent(it);
    }
    reverse(matrix.begin(), matrix.end());
    return matrix;
}

// Slices a matrix-like object whose rows support begin()/end(). Returns vectors
// because OR-Tools expressions and permutation routines use dynamically-sized rows.
template <typename MatrixType>
DenseMatrix sliceMat(const MatrixType &matrix,
                     int rowStart = 0,
                     int rowEnd = -1,
                     int colStart = 0,
                     int colEnd = -1) {
    DenseMatrix result;
    if (matrix.empty()) return result;

    const int matrixRows = static_cast<int>(matrix.size());
    if (rowEnd == -1) rowEnd = matrixRows;
    rowStart = max(rowStart, 0);
    rowEnd = min(rowEnd, matrixRows);
    if (rowStart >= rowEnd) return result;

    const int matrixCols = static_cast<int>(matrix[rowStart].size());
    if (colEnd == -1) colEnd = matrixCols;
    colStart = max(colStart, 0);
    colEnd = min(colEnd, matrixCols);
    if (colStart >= colEnd) return result;

    result.reserve(rowEnd - rowStart);
    for (int rowInd = rowStart; rowInd < rowEnd; ++rowInd) {
        const auto &currentRow = matrix[rowInd];
        result.emplace_back(currentRow.begin() + colStart, currentRow.begin() + colEnd);
    }
    return result;
}

// Helper Function for Row Candidates.
void genBisForEach(int selfInt,
                   int genus,
                   int a,
                   vector<int> &currentRow,
                   int currSelfIntSum,
                   int currAdjSum,
                   const function<void(const BisRow &)> &emit) {
    // Base cases for a <= 0.
    if (a <= 0) {
        if ((a == 0) || ((selfInt < -2) && (genus == 0))) {
            const int numOnes = 2 * a - selfInt - 1;
            if ((0 <= numOnes) && (numOnes <= N - 1)) {
                vector<int> bis(numOnes, -1);
                bis.insert(bis.end(), N - 1 - numOnes, 0);
                bis.push_back(1 - a);

                BisRow bisArray{};
                copy(bis.begin(), bis.end(), bisArray.begin());
                emit(bisArray);
            }
        }
        return;
    }

    const int currLen = static_cast<int>(currentRow.size());
    if (currLen == N) {
        BisRow bisArray{};
        for (int i = 0; i < N; ++i) bisArray[i] = -currentRow[i];
        emit(bisArray);
        return;
    }

    const int numBisLeft = N - currLen;
    const int selfIntSumLeft = a * a - selfInt - currSelfIntSum;
    const int adjSumLeft = (a - 1) * (a - 2) - 2 * genus - currAdjSum;
    if ((selfIntSumLeft < 0) || (adjSumLeft < 0)) return;

    int adjLowerBound = 0;
    if (adjSumLeft != 0) {
        adjLowerBound = static_cast<int>(sqrt((static_cast<double>(adjSumLeft) / numBisLeft) + 0.25) + 0.5);
    }

    const int nextBiMin = static_cast<int>(ceil(max({sqrt(static_cast<double>(selfIntSumLeft) / numBisLeft),
                                                       static_cast<double>(adjLowerBound),
                                                       0.0})));
    const double placeholderUpperBound = currentRow.empty()
                                             ? 1000000.0
                                             : static_cast<double>(currentRow.back());
    const int nextBiMax = static_cast<int>(floor(min({sqrt(static_cast<double>(selfIntSumLeft)),
                                                       sqrt(static_cast<double>(adjSumLeft) + 0.25) + 0.5,
                                                       placeholderUpperBound})));
    if (nextBiMin > nextBiMax) return;

    auto recurseWith = [&](int nextBi) {
        currentRow.push_back(nextBi);
        genBisForEach(selfInt,
                      genus,
                      a,
                      currentRow,
                      currSelfIntSum + nextBi * nextBi,
                      currAdjSum + nextBi * nextBi - nextBi,
                      emit);
        currentRow.pop_back();
    };

    const int middleBi = (nextBiMin + nextBiMax) / 2;
    recurseWith(middleBi);

    int biOffset = 1;
    while (biOffset <= (nextBiMax - nextBiMin + 1) / 2) {
        int nextBi = middleBi - biOffset;
        if (nextBi >= nextBiMin) recurseWith(nextBi);

        nextBi = middleBi + biOffset;
        if (nextBi <= nextBiMax) recurseWith(nextBi);

        ++biOffset;
    }
}

// Computes the possible first rows of the matrices.
void rowCandsForEach(int selfInt,
                     int genus,
                     int aMin,
                     int aMax,
                     const function<void(const Row &)> &emit) {
    int aLowerBound = 0;
    if (selfInt < -2) aLowerBound = static_cast<int>(ceil(static_cast<double>(selfInt + 1) / 2.0));
    const int startA = max(aLowerBound, aMin);

    vector<int> currentRow;
    currentRow.reserve(N);
    for (int a = startA; a <= aMax; ++a) {
        genBisForEach(selfInt, genus, a, currentRow, 0, 0, [&](const BisRow &bis) {
            Row row{};
            row[0] = a;
            copy(bis.begin(), bis.end(), row.begin() + 1);
            emit(row);
        });
    }
}

// Helper for solveNextRow: get sorted elements from map counter.
vector<int> counterElements(const map<int, int> &ctr) {
    vector<int> elems;
    for (const auto &[val, count] : ctr) {
        for (int i = 0; i < count; ++i) elems.push_back(val);
    }
    sort(elems.begin(), elems.end());
    return elems;
}

// Helper for solveNextRow: dot product of a row of prevBiMat with a vector.
int dotProduct(const DenseMatrix &mat, int row, const vector<int> &vec) {
    int sum = 0;
    for (int j = 0; j < static_cast<int>(vec.size()); ++j) sum += mat[row][j] * vec[j];
    return sum;
}

// Helper for solveNextRow: sum of a column subset of a matrix row.
int sumCols(const DenseMatrix &mat, int row, const vector<int> &colInds) {
    int sum = 0;
    for (int col : colInds) sum += mat[row][col];
    return sum;
}

bool columnsEqual(const DenseMatrix &mat, int col1, int col2) {
    for (const auto &row : mat) {
        if (row[col1] != row[col2]) return false;
    }
    return true;
}

// Helper for solveNextRow: delete columns from matrix.
DenseMatrix deleteCols(const DenseMatrix &mat, const vector<int> &colInds) {
    DenseMatrix result;
    if (mat.empty()) return result;

    vector<char> deleteCol(mat[0].size(), 0);
    for (int col : colInds) {
        if (0 <= col && col < static_cast<int>(deleteCol.size())) deleteCol[col] = 1;
    }

    result.reserve(mat.size());
    for (const auto &row : mat) {
        vector<int> newRow;
        newRow.reserve(row.size() - colInds.size());
        for (int j = 0; j < static_cast<int>(row.size()); ++j) {
            if (!deleteCol[j]) newRow.push_back(row[j]);
        }
        result.push_back(std::move(newRow));
    }
    return result;
}

// Helper for solveNextRow: generate combinations of r indices from 0..n-1.
void combinationsForEach(int n, int r, const function<void(const vector<int> &)> &emit) {
    if (r < 0 || r > n) return;
    vector<int> indices(r);
    iota(indices.begin(), indices.end(), 0);
    emit(indices);

    while (true) {
        int i = r - 1;
        while (i >= 0 && indices[i] == i + n - r) --i;
        if (i < 0) return;
        ++indices[i];
        for (int j = i + 1; j < r; ++j) indices[j] = indices[j - 1] + 1;
        emit(indices);
    }
}

// Helper for solveNextRow.
void validBiPermsForEach(const vector<int> &goalDots,
                         const DenseMatrix &prevBiMat,
                         map<int, int> &bCtr,
                         const function<void(const vector<int> &)> &emit) {
    const int numRows = static_cast<int>(prevBiMat.size());

    // Sort prevBiMat rows.
    DenseMatrix sortedPrevBis = prevBiMat;
    for (auto &row : sortedPrevBis) sort(row.begin(), row.end());

    const vector<int> sortedNewBis = counterElements(bCtr);
    const vector<int> flippedNewBis(sortedNewBis.rbegin(), sortedNewBis.rend());

    // Compute maxDots and minDots.
    vector<int> maxDots(numRows), minDots(numRows);
    for (int i = 0; i < numRows; ++i) {
        maxDots[i] = dotProduct(sortedPrevBis, i, sortedNewBis);
        minDots[i] = dotProduct(sortedPrevBis, i, flippedNewBis);
    }

    // Check bounds.
    for (int i = 0; i < numRows; ++i) {
        if (goalDots[i] < minDots[i] || goalDots[i] > maxDots[i]) return;
    }

    // Base case: one unique element.
    if (bCtr.size() == 1) {
        emit(counterElements(bCtr));
        return;
    }

    const int numCols = static_cast<int>(sortedNewBis.size());

    // Find element with largest absolute value.
    auto maxIt = max_element(bCtr.begin(), bCtr.end(), [](const auto &a, const auto &b) {
        return abs(a.first) < abs(b.first);
    });
    const int biToAssign = maxIt->first;
    const int numCopies = maxIt->second;

    // Backtracking: remove for recursion, restore afterwards.
    bCtr.erase(maxIt);
    combinationsForEach(numCols, numCopies, [&](const vector<int> &colInds) {
        vector<int> newGoalDots(numRows);
        for (int i = 0; i < numRows; ++i) {
            newGoalDots[i] = goalDots[i] - biToAssign * sumCols(prevBiMat, i, colInds);
        }

        DenseMatrix newPrevBiMat = deleteCols(prevBiMat, colInds);
        validBiPermsForEach(newGoalDots, newPrevBiMat, bCtr, [&](const vector<int> &subPerm) {
            vector<int> fullPerm = subPerm;
            // Insert left-to-right to restore original order.
            for (int i = 0; i < static_cast<int>(colInds.size()); ++i) {
                fullPerm.insert(fullPerm.begin() + colInds[i], biToAssign);
            }
            emit(fullPerm);
        });
    });
    bCtr[biToAssign] = numCopies;
}

namespace operations_research {
namespace sat {

// Finds the next possible rows given a subMatrix.
// Call via operations_research::sat::solveNextRow().
RowMatrix solveNextRow(const RowMatrix &subMatRows,
                       const DenseMatrix &interNums,
                       int selfInt,
                       int genus,
                       int aMax,
                       int cpSatWorkers = kCpSatWorkersPerSolve) {
    const DenseMatrix subMat = sliceMat(subMatRows);
    const int numPrevRows = static_cast<int>(subMat.size());
    RowMatrix sols;

    vector<int> sameAsRightNeighbor;
    for (int i = 0; i < N - 1; ++i) {
        if (columnsEqual(subMat, i + 1, i + 2)) sameAsRightNeighbor.push_back(i);
    }

    // a > 0 case: solve using CP-SAT.
    if (aMax >= 1) {
        CpModelBuilder model;
        const IntVar a = model.NewIntVar(Domain(1, aMax)).WithName("a");

        int biUB = static_cast<int>(min(static_cast<double>(3 * aMax - selfInt + 2 * genus - 2),
                                        trunc(sqrt(static_cast<double>(aMax * aMax - selfInt)))));
        biUB = max(0, biUB);

        vector<IntVar> b;
        b.reserve(N);
        for (int i = 0; i < N; ++i) {
            b.push_back(model.NewIntVar(Domain(0, biUB)).WithName("b_" + to_string(i)));
        }

        for (int rowInd = 0; rowInd < numPrevRows; ++rowInd) {
            LinearExpr expr = LinearExpr::Term(a, subMat[rowInd][0]);
            for (int i = 0; i < N; ++i) {
                expr += LinearExpr::Term(b[i], subMat[rowInd][i + 1]);
            }
            model.AddEquality(expr, interNums[0][rowInd]);
        }

        model.AddEquality(LinearExpr::Sum(b), LinearExpr::Term(a, 3) - selfInt + 2 * genus - 2);

        vector<IntVar> bisSquared;
        bisSquared.reserve(N);
        for (int i = 0; i < N; ++i) {
            IntVar bSquared = model.NewIntVar(Domain(0, biUB * biUB)).WithName("b_" + to_string(i) + "_squared");
            model.AddMultiplicationEquality(bSquared, {b[i], b[i]});
            bisSquared.push_back(bSquared);
        }

        IntVar aSquared = model.NewIntVar(Domain(1, aMax * aMax)).WithName("a_squared");
        model.AddMultiplicationEquality(aSquared, {a, a});
        model.AddEquality(LinearExpr::Sum(bisSquared), LinearExpr::Term(aSquared, 1) - selfInt);

        for (int i : sameAsRightNeighbor) {
            model.AddGreaterOrEqual(b[i], b[i + 1]);
        }

        SatParameters parameters;
        parameters.set_num_search_workers(max(1, cpSatWorkers));
        parameters.set_enumerate_all_solutions(true);

        Model orModel;
        orModel.Add(NewSatParameters(parameters));
        orModel.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse &r) {
            Row sol{};
            sol[0] = static_cast<int>(SolutionIntegerValue(r, a));
            for (int i = 0; i < N; ++i) {
                sol[i + 1] = -static_cast<int>(SolutionIntegerValue(r, b[i]));
            }
            sols.push_back(sol);
        }));
        SolveCpModel(model.Build(), &orModel);
    }

    // a <= 0 case.
    int aMin = 0;
    if (selfInt < -2) aMin = static_cast<int>(ceil(static_cast<double>(selfInt + 1) / 2.0));

    for (int currA = aMin; currA < 1; ++currA) {
        const int numOnes = 2 * currA - selfInt - 1;
        if (0 <= numOnes && numOnes <= (N - 1)) {
            vector<int> goalBDots(numPrevRows);
            for (int i = 0; i < numPrevRows; ++i) {
                goalBDots[i] = subMat[i][0] * currA - interNums[0][i];
            }

            map<int, int> negBCounter;
            if (numOnes > 0) negBCounter[-1] = numOnes;
            if (N - 1 - numOnes > 0) negBCounter[0] = N - 1 - numOnes;
            negBCounter[-currA + 1] = 1;

            validBiPermsForEach(goalBDots, sliceMat(subMat, 0, -1, 1, -1), negBCounter, [&](const vector<int> &perm) {
                bool discard = false;
                for (int i : sameAsRightNeighbor) {
                    if (perm[i] < perm[i + 1]) {
                        discard = true;
                        break;
                    }
                }
                if (!discard) {
                    Row sol{};
                    sol[0] = currA;
                    copy(perm.begin(), perm.end(), sol.begin() + 1);
                    sols.push_back(sol);
                }
            });
        }
    }

    return sols;
}

}  // namespace sat
}  // namespace operations_research

struct LeafResult {
    bool invalidDepth = false;
    RowMatrix nextRows;
};

// Main function using trees. Computes everything at once; kept for local checkpoint generation.
int mainAll() {
    int rowInd = 0;

    // Sort matrix rows by key.
    Row keyArray{};
    for (int i = 0; i < M; ++i) keyArray[i] = aMaxes[i] * aMaxes[i] - intMat[i][i];

    Row indices{};
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&keyArray](int i, int j) { return keyArray[i] < keyArray[j]; });

    RowMatrix intMatReOrd(M);
    Row aMaxesReOrd{};
    Row generaReOrd{};
    for (int i = 0; i < M; ++i) {
        aMaxesReOrd[i] = aMaxes[indices[i]];
        generaReOrd[i] = genera[indices[i]];
        for (int j = 0; j < M; ++j) intMatReOrd[i][j] = intMat[indices[i]][indices[j]];
    }

    MatrixTree matTree;
    Row firstNode{};
    auto root = matTree.insert(matTree.begin(), firstNode);

    RowMatrix intMatRowSlice = {intMatReOrd[rowInd + 1]};
    DenseMatrix intMatSlice = sliceMat(intMatRowSlice, 0, 1, 0, rowInd + 1);
    int count = 0;
    int countFirstRows = 0;

    rowCandsForEach(intMatReOrd[0][0], generaReOrd[0], -1000000, aMaxesReOrd[0], [&](const Row &firstRow) {
        RowMatrix firstRowVect = {firstRow};
        RowMatrix nextRows = operations_research::sat::solveNextRow(firstRowVect,
                                                                     intMatSlice,
                                                                     intMatReOrd[rowInd + 1][rowInd + 1],
                                                                     generaReOrd[rowInd + 1],
                                                                     aMaxesReOrd[rowInd + 1]);
        auto firstRowNode = matTree.append_child(root, firstRow);
        ++countFirstRows;
        for (const auto &newRow : nextRows) {
            ++count;
            matTree.append_child(firstRowNode, newRow);
        }
    });

    ++rowInd;
    cout << countFirstRows << '\n' << flush;

    while (rowInd < N) {
        cout << count << '\n' << flush;
        count = 0;

        vector<MatrixTree::iterator> leaves;
        for (auto it = matTree.begin_leaf(); it != matTree.end_leaf(); ++it) leaves.push_back(it);

        RowMatrix nextIntMatRowSlice = {intMatReOrd[rowInd + 1]};
        intMatSlice = sliceMat(nextIntMatRowSlice, 0, 1, 0, rowInd + 2);

        for (const auto &leafIt : leaves) {
            RowMatrix subMat = reconstructMatrix(matTree, leafIt);
            if (static_cast<int>(subMat.size()) != rowInd + 1) {
                matTree.erase(leafIt);
                continue;
            }

            RowMatrix nextRows = operations_research::sat::solveNextRow(subMat,
                                                                         intMatSlice,
                                                                         intMatReOrd[rowInd + 1][rowInd + 1],
                                                                         generaReOrd[rowInd + 1],
                                                                         aMaxesReOrd[rowInd + 1]);
            if (nextRows.empty()) {
                matTree.erase(leafIt);
                continue;
            }

            for (const auto &newRow : nextRows) {
                matTree.append_child(leafIt, newRow);
                ++count;
            }
        }

        ++rowInd;
        pruneShallowLeaves(matTree, rowInd + 1);
        saveTree(matTree, "matTree_" + to_string(rowInd) + ".tr");
        matTree = loadTree("matTree_" + to_string(rowInd) + ".tr");
    }

    cout << count << '\n' << flush;
    return 0;
}

// Computes one step at a time and saves the output.
// Usage: ./main inputFileName.tr outputFileName.tr
int main(int argc, char *argv[]) {
    try {
        if (argc != 3) {
            cerr << "Usage: " << argv[0] << " inputFileName.tr outputFileName.tr\n";
            return 1;
        }

        const string inputFileName = argv[1];
        const string outputFileName = argv[2];

        // Sort problem data by keyArray.
        Row keyArray{};
        for (int i = 0; i < M; ++i) keyArray[i] = aMaxes[i] * aMaxes[i] - intMat[i][i];

        Row indices{};
        iota(indices.begin(), indices.end(), 0);
        sort(indices.begin(), indices.end(), [&keyArray](int i, int j) { return keyArray[i] < keyArray[j]; });

        RowMatrix intMatReOrd(M);
        Row aMaxesReOrd{};
        Row generaReOrd{};
        for (int i = 0; i < M; ++i) {
            aMaxesReOrd[i] = aMaxes[indices[i]];
            generaReOrd[i] = genera[indices[i]];
            for (int j = 0; j < M; ++j) intMatReOrd[i][j] = intMat[indices[i]][indices[j]];
        }

        MatrixTree matTree = loadTree(inputFileName);

        vector<MatrixTree::iterator> leaves;
        for (auto it = matTree.begin_leaf(); it != matTree.end_leaf(); ++it) leaves.push_back(it);

        // A file containing only the sentinel root has one leaf, namely the root.
        leaves.erase(remove_if(leaves.begin(), leaves.end(), [&](const auto &it) { return it == matTree.begin(); }), leaves.end());

        if (leaves.empty()) {
            cerr << argv[0] << ' ' << inputFileName << ' ' << outputFileName
                 << " ran into an error: MatTree is empty; please select a different file.\n";
            return 1;
        }

        const RowMatrix firstMat = reconstructMatrix(matTree, leaves.front());
        const int currentDepth = static_cast<int>(firstMat.size());
        const int nextRowIndex = currentDepth;
        if (nextRowIndex >= M) {
            cerr << "Input tree already appears complete at depth " << currentDepth << ". No next row exists.\n";
            return 1;
        }

        RowMatrix intMatRowSlice = {intMatReOrd[nextRowIndex]};
        const DenseMatrix intMatSlice = sliceMat(intMatRowSlice, 0, 1, 0, nextRowIndex + 1);
        const int selfInt = intMatReOrd[nextRowIndex][nextRowIndex];
        const int genus = generaReOrd[nextRowIndex];
        const int aMax = aMaxesReOrd[nextRowIndex];

        const int numThreads = max(1, min({requestedOuterWorkers(), static_cast<int>(leaves.size()), kMaxOuterWorkers}));
        cout << "Processing " << leaves.size() << " leaves at depth " << currentDepth
             << " with " << numThreads << " worker threads; CP-SAT workers per solve = "
             << kCpSatWorkersPerSolve << ".\n" << flush;

        vector<LeafResult> results(leaves.size());
        atomic<size_t> nextLeaf{0};
        mutex treeReadMutex;
        mutex exceptionMutex;
        exception_ptr firstException = nullptr;

        auto worker = [&]() {
            while (true) {
                const size_t idx = nextLeaf.fetch_add(1, memory_order_relaxed);
                if (idx >= leaves.size()) break;

                try {
                    RowMatrix subMat;
                    {
                        lock_guard<mutex> lock(treeReadMutex);
                        subMat = reconstructMatrix(matTree, leaves[idx]);
                    }

                    if (static_cast<int>(subMat.size()) != currentDepth) {
                        results[idx].invalidDepth = true;
                        continue;
                    }

                    results[idx].nextRows = operations_research::sat::solveNextRow(subMat,
                                                                                   intMatSlice,
                                                                                   selfInt,
                                                                                   genus,
                                                                                   aMax,
                                                                                   kCpSatWorkersPerSolve);
                } catch (...) {
                    lock_guard<mutex> lock(exceptionMutex);
                    if (firstException == nullptr) firstException = current_exception();
                }
            }
        };

        vector<thread> threads;
        threads.reserve(numThreads);
        for (int i = 0; i < numThreads; ++i) threads.emplace_back(worker);
        for (auto &thread : threads) thread.join();

        if (firstException != nullptr) rethrow_exception(firstException);

        int count = 0;
        vector<MatrixTree::iterator> toErase;
        for (size_t idx = 0; idx < leaves.size(); ++idx) {
            auto leafIt = leaves[idx];
            if (results[idx].invalidDepth || results[idx].nextRows.empty()) {
                toErase.push_back(leafIt);
                continue;
            }

            for (const auto &newRow : results[idx].nextRows) {
                matTree.append_child(leafIt, newRow);
                ++count;
            }
        }

        for (auto it = toErase.rbegin(); it != toErase.rend(); ++it) matTree.erase(*it);

        const int newDepth = currentDepth + 1;
        cout << "Number of leaves at depth " << newDepth << ": " << count << '\n' << flush;

        pruneShallowLeaves(matTree, newDepth);
        saveTree(matTree, outputFileName);
        return 0;
    } catch (const exception &e) {
        cerr << "Fatal error: " << e.what() << '\n';
        return 2;
    }
}
