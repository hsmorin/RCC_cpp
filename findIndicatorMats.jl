import NPZ as np
using ResumableFunctions
import ORTools_jll
using JuMP, SCIP
import SCIP as scp
import MiniZinc as MiniZinc

#key = "7c786db7-e071-4b37-a314-0a8a8bb7a604"
#run(`$(Gurobi_jll.grbgetkey()) $key`)

include("globalVars.jl")

@resumable function genBis(N_local::Int64, selfInt::Int64, genus::Int64, a, currentRow = [], currSelfIntSum = 0, currAdjSum = 0)

# base cases
    if a <= 0
        if a == 0 || (selfInt < -2 && genus == 0) # Lemma 3.4
            numOnes = 2 * a - selfInt - 1
            if 0 <= numOnes <= N_local - 1
                @yield vcat([-1 for i in 1:numOnes], zeros(Int64, N_local - 1 - numOnes), [1 - a])
                return
            end
        return
        end
    end

    currLen = length(currentRow)
    if currLen == N_local # row is complete!
        @yield [-bi for bi in currentRow]
        return
    end 

    # recursive case - run method for all possible next b_i's such that self-intersection and adjunction are still feasible
    numBisLeft = N_local - currLen
    selfIntSumLeft = a^2 - selfInt - currSelfIntSum
    adjSumLeft = (a - 1) * (a - 2) - 2 * genus - currAdjSum
    if selfIntSumLeft < 0 || adjSumLeft < 0
        return
    end

    if adjSumLeft == 0 
        adjLowerBound = 0  
    else 
        adjLowerBound = sqrt((adjSumLeft / numBisLeft) + 0.25) + 0.5 # fix quadratic-ness of adjunction summands (0 at both 0 and 1)
    end
    
    nextBiMin = ceil(Int64, maximum([sqrt(selfIntSumLeft / numBisLeft), adjLowerBound, 0]))
    
    if !isempty(currentRow)
        upperBound = last(currentRow)
    else
        upperBound = 1e6
    end

    nextBiMax = minimum([sqrt(selfIntSumLeft), sqrt(adjSumLeft + 0.25) + 0.5, upperBound])
    nextBiMax = floor(Int64, nextBiMax)
    
    #TO DO: Yield from?
    if nextBiMin > nextBiMax
        return
    end

    middleBi = (nextBiMin + nextBiMax) ÷ 2
     
    for i in genBis(N_local, selfInt, genus, a, vcat(currentRow, [middleBi]), currSelfIntSum + middleBi^2, currAdjSum + (middleBi^2 - middleBi))
        @yield i 
    end
    biOffset = 1
    while biOffset <= (nextBiMax - nextBiMin + 1) ÷ 2
        nextBi = middleBi - biOffset
        if nextBi >= nextBiMin
            for i in genBis(N_local, selfInt, genus, a, vcat(currentRow, [nextBi]), currSelfIntSum + nextBi^2, currAdjSum + (nextBi^2 - nextBi))
                @yield i 
            end
        end
        nextBi = middleBi + biOffset
        if nextBi <= nextBiMax
            for i in genBis(N_local, selfInt, genus, a, vcat(currentRow, [nextBi]), currSelfIntSum + nextBi^2, currAdjSum + (nextBi^2 - nextBi))
                @yield i
            end
        end
        biOffset += 1
    end
end


# runs genBis for a range of values of "a" and generates solution lists in the form [a, -b_1, -b_2, ..., -b_N_local], with (-b_i)'s ascending
@resumable function rowCands(N_local, selfInt, genus, aMin, aMax)

    if selfInt >= -2
        aLowerBound = 0
    else
        aLowerBound = ceil((selfInt + 1) / 2) 
    end

    startA = maximum(Int64, [aLowerBound, aMin])
    for a in startA:trunc(Int64, aMax)
        biVecs = genBis(N_local, selfInt, genus, a)
        for biVec in biVecs
            @yield vcat([a], biVec)
        end
    end
end    

function solveNextRow(subMat, interNums, selfInt::Int64, genus::Int64, aMax::Int64)
    
    if length(subMat[1]) == 1
        numPrevRows, N_local = 1, length(subMat) - 1
        subMat = reshape(subMat, (1,:))
    else
        numPrevRows, N_local = length(subMat)[1], length(subMat[1]) - 1
    end

    path = joinpath(ORTools_jll.artifact_dir, "share", "minizinc", "solvers", "cp-sat.msc") 
    model = Model(() -> MiniZinc.Optimizer{Float64}("highs"))
    #model = Model(SCIP.Optimizer)
    #set_silent(model)

    @variable(model, 1 <= a <= aMax, Int)

    biUB = minimum(Int64, [3 * aMax - selfInt + 2 * genus - 2, trunc(Int64, sqrt(aMax^2 - selfInt))])
    @variable(model, 0 <= b[i=1:N_local] <= biUB, Int)

    # intersection number constraints: for all rowInd, subMat[rowInd] . [a, b_1, ..., b_N] = interNums[rowInd]
    for rowInd in 1:numPrevRows
        @constraint(model, 
                    trunc(Int64, subMat[rowInd, 1]) * a + sum([trunc(Int64, subMat[rowInd, i + 1]) * b[i] for i in 1:N_local]) == trunc(Int64, interNums[rowInd]))
    end

    # adjunction formula
    @constraint(model, 
                sum([b[i] for i in 1:N_local]) == 3 * a - trunc(Int64, selfInt + 2 * genus - 2))

    # self-intersection
    @constraint(model,
                sum([b[i]^2 for i in 1:N_local]) == a^2 - selfInt)

    # symmetry breaking for repeated columns in subMat; assumes that subMat has columns in sortBis() order
    sameAsRightNeighbor = [] # saving this for later use in a <= 0 case
    for i in 1:(N_local - 1)
        if subMat[:,i+1] == subMat[:,i+2]  # columns for b[i] vs b[i+1]
            @constraint(model, b[i] >= b[i + 1])
            sameAsRightNeighbor = vcat(sameAsRightNeighbor, i)
        end
    end
    
    # collect all solutions via no-good cuts
    all_solutions = []
    
    function add_no_good_cut(model, all_solutions)

        # get current solution
        a_val = round(Int, value(a))
        b_vals = [round(Int, value(b[i])) for i in 1:N_local]
        sol = vcat(a_val, b_vals)
        push!(all_solutions, sol)
        println("Found solution: $sol")
        
        # add constraint ruling out this exact solution
        # at least one variable must differ from its current value
        @constraint(model,
            (a - a_val)^2 + sum((b[i] - b_vals[i])^2 for i in 1:N_local) >= 1)
    end

    # iteratively solve and exclude
    while true
        optimize!(model)
        status = termination_status(model)
        println(status) 
        if status ∈ (MOI.OPTIMAL, MOI.FEASIBLE_POINT, MOI.LOCALLY_SOLVED)
            add_no_good_cut(model, all_solutions)
        else
            break
        end
    end
end

function main_cluster()
    
    jobNum = parse(Int64, ARGS[1])
    rowInd = parse(Int64, ARGS[2])
    version = parse(Int64, ARGS[3])

    if version != 0
        suffix = "_v$(version).npy"
    else 
        suffix = ".npy"
    end

    println("Job for row index $(rowInd) and version $(version) with suffix $(suffix).")

    Q = np.npzread("./data/intMat.npy") 
    aMaxes = np.npzread("./data/aBounds.npy")
    genera = np.npzread("./data/genus.npy")
    println(N)
    println(size(Q)[1])

    if N != size(Q)[1]
        throw("Intersection matrix not compatable with N.")
    end

    if N != length(aMaxes)
        throw("Intersection matrix not compatable with N.")
    end

    if N != length(genera)
        throw("Intersection matrix not compatable with N.")
    end

    rowOrdInds = sortperm([aMaxes[i]^2 - Q[i, i] for i in range(1, N)])
    intMatReOrd = Q[rowOrdInds,rowOrdInds]
    aMaxesReOrd = aMaxes[rowOrdInds]
    generaReOrd = genera[rowOrdInds]
    
    #if rowInd == -1
        #TO DO: Parse
        #np.npzwrite("./data/$(projectName)0.npy", [np.array(firstRow).reshape(1, -1) for firstRow in rowCands(N - 1, intMatReOrd[0, 0], generaReOrd[0], aMin = -1e6, aMax = aMaxesReOrd[0])])
    #else
        newSubMatSet = Set()
        #subMats = npz.read(f"./{projectName}{jobInd}_{rowInd}{suffix}")
    #end
end

function main_local()
    Q = np.npzread("./data/intMat.npy") 
    aMaxes = np.npzread("./data/aBounds.npy")
    genera = np.npzread("./data/genus.npy")

    if N != size(Q)[1] - 1
        throw("Intersection matrix not compatable with N.")
    end

    if N != length(aMaxes) - 1
        throw("aMaxes not compatable with N.")
    end

    if N != length(genera) - 1
        throw("genera not compatable with N.")
    end

    #rowOrdInds = sortperm([aMaxes[i]^2 - Q[i, i] for i in 1:(N + 1)])
    rowOrdInds = 1:(N + 1)
    intMatReOrd = Q[rowOrdInds,rowOrdInds]
    aMaxesReOrd = aMaxes[rowOrdInds]
    generaReOrd = genera[rowOrdInds]

    
    firstRows = Set(rowCands(N, intMatReOrd[1,1], generaReOrd[1], -1e6, aMaxesReOrd[1]))
    
    rowInd = 1
    newSubMatSet = Set() 
    subMats = firstRows

    for subMat in subMats
        println(subMat)
        nextRows = solveNextRow(subMat, intMatReOrd[rowInd + 1, 1:rowInd + 1], intMatReOrd[rowInd + 1, rowInd + 1], generaReOrd[rowInd + 1], aMaxesReOrd[rowInd + 1])
        #newMats = [vcat(subMat, newRow) for newRow in nextRows]
        #for newMat in newMats
            #TO DO SortBis
            #newSubMatSet.add(newMat)
        #subMats = [mat for mat in newSubMatSet]
        #end
    end
end

main_local()
