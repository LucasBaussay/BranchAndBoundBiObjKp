
include("dataStruct.jl")

@enum Neighborhood KP_1_1 KP_1_2 KP_2_1 KP_0_1
@enum Dominance DOMINATES DOMINATED NOTHING

function d(solA::Solution, solB::Solution)
    return sqrt((solB.y[1] - solA.y[1])^2 + (solB.y[2] - solA.y[2])^2)
end

function d(yA::Vector{Float64}, yB::Vector{Float64})
    #@assert length(yA) == 2 && length(yB) == 2 "We only work with 2D"
    return sqrt((yB[1] - yA[1])^2 + (yB[2] - yA[2])^2)
end

function d(yA::Vector{Float64}, yB::Vector{Float64})
    #@assert length(yA) == 2 && length(yB) == 2 "We only work with 2D"
    return sqrt((yB[1] - yA[1])^2 + (yB[2] - yA[2])^2)
end

function middle(yA::Vector{Float64}, yB::Vector{Float64})
    #@assert length(yA) == 2 && length(yB) == 2 "We only work with 2D"
    return [(yA[1]+yB[1])/2, (yA[2]+yB[2])/2]
end

function evaluate(y::Vector{Float64}, consecutiveSet::Vector{PairOfSolution}; fastEvaluate = false)
    #@assert length(y) == 2 "We only work with 2D"
    minEval = Inf
    for pairOfSol in consecutiveSet
        solL,solR = pairOfSol.solL, pairOfSol.solR
        middlePoint = middle(solL.y,solR.y)
        if !fastEvaluate
            minLengthTriangle = d(solL,solR)
            nadir = [solL.y[1],solR.y[2]]
            if d(solL.y,nadir) < minLengthTriangle
                minLengthTriangle = d(solL.y,nadir)
            end
            if d(solR.y,nadir) < minLengthTriangle
                minLengthTriangle = d(solR.y,nadir)
            end
            w = 1/minLengthTriangle
            evalSol = d(y,middlePoint)*w
        else
            evalSol = d(y,middlePoint)
        end
        if evalSol < minEval
            minEval = evalSol
        end
    end
    return minEval
end

function bestKP01(sol::Solution, consecutiveSet::Vector{PairOfSolution}, prob::Problem; fastEvaluate = false)
    bestKP = [0]
    bestEval = evaluate(sol.y,consecutiveSet,fastEvaluate=fastEvaluate)
    evalBuffer = 0.0
    for up in 1:prob.nbVar
        if sol.x[up] == 0
            # test is the kp is feasible
            if sol.w + prob.weights[up] <= prob.maxWeight # feasible
                #println("feasible : $(sol.w) + $(prob.weights[up]) <= $(prob.maxWeight)")
                # evaluate new solution
                y = [sol.y[1]+prob.profits[1,up],sol.y[2]+prob.profits[2,up]]
                evalBuffer = evaluate(y,consecutiveSet,fastEvaluate=fastEvaluate)
                if evalBuffer < bestEval
                    bestEval = evalBuffer
                    bestKP[1] = up
                end
            end
        end
    end
    return bestKP
end

function firstKP01(sol::Solution, consecutiveSet::Vector{PairOfSolution}, prob::Problem; fastEvaluate = false)
    bestKP = [0]
    bestEval = evaluate(sol.y,consecutiveSet,fastEvaluate=fastEvaluate)
    evalBuffer = 0.0
    improved = false
    up = 0
    while !improved && up < prob.nbVar
        up += 1
        if sol.x[up] == 0
            # test is the kp is feasible
            if sol.w + prob.weights[up] <= prob.maxWeight # feasible
                #println("feasible : $(sol.w) + $(prob.weights[up]) <= $(prob.maxWeight)")
                # evaluate new solution
                y = [sol.y[1]+prob.profits[1,up],sol.y[2]+prob.profits[2,up]]
                evalBuffer = evaluate(y,consecutiveSet,fastEvaluate=fastEvaluate)
                if evalBuffer < bestEval
                    improved = true
                    bestEval = evalBuffer
                    bestKP[1] = up
                end
            end
        end
    end
    return bestKP
end

function bestKP12(sol::Solution, consecutiveSet::Vector{PairOfSolution}, prob::Problem; fastEvaluate = false)
    bestKP = [0,0,0]
    bestEval = evaluate(sol.y,consecutiveSet,fastEvaluate=fastEvaluate)
    evalBuffer = 0.0
    for down in 1:prob.nbVar
        if sol.x[down] == 1
            for up in 1:prob.nbVar
                if sol.x[up] == 0
                    # test if the kp could be feasible
                    if sol.w - prob.weights[down] + prob.weights[up] < prob.maxWeight
                        for up2 in 1:prob.nbVar
                            if sol.x[up2] == 0
                                # test is the kp is feasible
                                if sol.w + prob.weights[up] + prob.weights[up2] - prob.weights[down] <= prob.maxWeight # feasible
                                    # evaluate new solution
                                    y = [sol.y[1]+prob.profits[1,up]+prob.profits[1,up2]-prob.profits[1,down],sol.y[2]+prob.profits[2,up]+prob.profits[2,up2]-prob.profits[2,down]]
                                    evalBuffer = evaluate(y,consecutiveSet,fastEvaluate=fastEvaluate)
                                    if evalBuffer < bestEval
                                        bestEval = evalBuffer
                                        bestKP[1] = down
                                        bestKP[2] = up
                                        bestKP[3] = up2
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return bestKP
end

function firstKP12(sol::Solution, consecutiveSet::Vector{PairOfSolution}, prob::Problem; fastEvaluate = false)
    bestKP = [0,0,0]
    bestEval = evaluate(sol.y,consecutiveSet,fastEvaluate=fastEvaluate)
    evalBuffer = 0.0
    improved = false
    down = 0
    while !improved && down < prob.nbVar
        down += 1
        if sol.x[down] == 1
            up = 0
            while !improved && up < prob.nbVar
                up += 1
                if sol.x[up] == 0
                    # test if the kp could be feasible
                    if sol.w - prob.weights[down] + prob.weights[up] < prob.maxWeight
                        up2 = 0
                        while !improved && up2 < prob.nbVar
                            up2 += 1
                            if sol.x[up2] == 0
                                # test is the kp is feasible
                                if sol.w + prob.weights[up] + prob.weights[up2] - prob.weights[down] <= prob.maxWeight # feasible
                                    # evaluate new solution
                                    y = [sol.y[1]+prob.profits[1,up]+prob.profits[1,up2]-prob.profits[1,down],sol.y[2]+prob.profits[2,up]+prob.profits[2,up2]-prob.profits[2,down]]
                                    evalBuffer = evaluate(y,consecutiveSet,fastEvaluate=fastEvaluate)
                                    if evalBuffer < bestEval
                                        improved = true
                                        bestEval = evalBuffer
                                        bestKP[1] = down
                                        bestKP[2] = up
                                        bestKP[3] = up2
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return bestKP
end

function bestKP11(sol::Solution, consecutiveSet::Vector{PairOfSolution}, prob::Problem, tEval::Float64, tE2::Float64; fastEvaluate = false)
    startFun = time()
    bestKP = [0,0]
    bestEval = evaluate(sol.y,consecutiveSet,fastEvaluate=fastEvaluate)
    evalBuffer = 0.0
    for down in 1:prob.nbVar
        if sol.x[down] == 1
            for up in 1:prob.nbVar
                if sol.x[up] == 0
                    # test is the kp is feasible
                    if sol.w + prob.weights[up] - prob.weights[down] <= prob.maxWeight # feasible
                        # evaluate new solution
                        y = [sol.y[1]+prob.profits[1,up]-prob.profits[1,down],sol.y[2]+prob.profits[2,up]-prob.profits[2,down]]
                        start = time()
                        evalBuffer = evaluate(y,consecutiveSet,fastEvaluate=fastEvaluate)
                        tEval += time() - start
                        if evalBuffer < bestEval
                            bestEval = evalBuffer
                            bestKP[1] = down
                            bestKP[2] = up
                        end
                    end
                end
            end
        end
    end
    return bestKP, tEval, tE2 + time() - startFun
end

function firstKP11(sol::Solution, consecutiveSet::Vector{PairOfSolution}, prob::Problem; fastEvaluate = false)
    bestKP = [0,0]
    bestEval = evaluate(sol.y,consecutiveSet,fastEvaluate=fastEvaluate)
    evalBuffer = 0.0
    improved = false
    down = 0
    while !improved && down < prob.nbVar
        down += 1
        if sol.x[down] == 1
            up = 0
            while !improved && up < prob.nbVar
                up += 1
                if sol.x[up] == 0
                    if sol.w + prob.weights[up] - prob.weights[down] <= prob.maxWeight # feasible
                        # evaluate new solution
                        y = [sol.y[1]+prob.profits[1,up]-prob.profits[1,down],sol.y[2]+prob.profits[2,up]-prob.profits[2,down]]
                        evalBuffer = evaluate(y,consecutiveSet,fastEvaluate=fastEvaluate)
                        if evalBuffer < bestEval
                            improved = true
                            bestEval = evalBuffer
                            bestKP[1] = down
                            bestKP[2] = up
                        end
                    end
                end
            end
        end
    end
    return bestKP
end

function updateSolByKP(sol::Solution,kp::Vector{Int}, neighborhood::Neighborhood, prob::Problem)
    newX = copy(sol.x)
    if neighborhood == KP_1_1
        newX[kp[1]] = 0
        newX[kp[2]] = 1
        newY = [sol.y[1]+prob.profits[1,kp[2]]-prob.profits[1,kp[1]],sol.y[2]+prob.profits[2,kp[2]]-prob.profits[2,kp[1]]]
        newW = sol.w + prob.weights[kp[2]] - prob.weights[kp[1]]
    elseif neighborhood == KP_0_1
        newX[kp[1]] = 1
        newY = [sol.y[1]+prob.profits[1,kp[1]],sol.y[2]+prob.profits[2,kp[1]]]
        newW = sol.w + prob.weights[kp[1]]
    elseif neighborhood == KP_1_2
        newX[kp[1]] = 0
        newX[kp[2]] = 1
        newX[kp[3]] = 1
        newY = [sol.y[1]+prob.profits[1,kp[2]]+prob.profits[1,kp[3]]-prob.profits[1,kp[1]],sol.y[2]+prob.profits[2,kp[2]]+prob.profits[2,kp[3]]-prob.profits[2,kp[1]]]
        newW = sol.w + prob.weights[kp[2]] + prob.weights[kp[3]] - prob.weights[kp[1]]
    end
    
    return Solution(newX,newY,newW)
end

function isNewSol(y::Vector{Float64}, kp::Vector{Int}, neighborhood::Neighborhood, allSolsI::Vector{Solution})
    for oldSol in allSolsI
        if oldSol.y == y
            return false
        end
    end
end

"""
    returns the index in solsInTriangles where sol should be added (triangle idex)
"""
function getTriangleIndex(solsInTriangles::Vector{Vector{Solution}}, sol::Solution)
    for i in 1:length(solsInTriangles)
        minZ1 = solsInTriangles[i][1].y[1]
        maxZ1 = solsInTriangles[i][2].y[1]
        if sol.y[1] >= minZ1 && sol.y[1] <= maxZ1
            return i
        end
    end
    @assert false "No triangle is fitting the solution"
end

function secureList(solsInTriangles::Vector{Vector{Solution}})
    for i in 1:length(solsInTriangles)
        for j in 3:length(solsInTriangles[i])
            for i1 in 1:length(solsInTriangles)
                for j1 in 1:length(solsInTriangles[i1])
                    if i1 != i && j != j1
                        @assert solsInTriangles[i][j].y != solsInTriangles[i1][j1].y "pas d'unicité dans la liste"
                    end
                end
            end
        end
    end
end

function whichDominates(sol1, sol2)
    sol1Dominates = true
    sol2Dominates = true

    for i in 1:2
        if sol1.y[i] < sol2.y[i]
            sol1Dominates = false
        end
        if sol1.y[i] > sol2.y[i]
            sol2Dominates = false
        end
    end

    if sol1Dominates && sol2Dominates # the Δs are the same
        @assert false "Those sols shouldn't be the same"
    elseif !sol1Dominates && !sol2Dominates # we can't compare the two kp-exchanges
        return NOTHING
    elseif sol1Dominates
        return DOMINATES
    else
        return DOMINATED
    end
end

function size2DList(list::Vector{Vector{Solution}})
    size = 0
    for i in 1:length(list)
        size += length(list[i])
    end
    return size
end

function cleanList(solsInTriangles::Vector{Vector{Solution}})
    newList = Vector{Vector{Solution}}(undef,length(solsInTriangles))
    for i in 1:length(newList)
        newList[i] = [solsInTriangles[i][1],solsInTriangles[i][2]]
    end

    fnished = false
    for i in 1:length(solsInTriangles)
        for j in 3:length(solsInTriangles[i])
            finished = false
            for i1 in 1:length(solsInTriangles)
                for j1 in 1:length(solsInTriangles[i1])
                    if i1 != i || j != j1
                        dominance_indicator = whichDominates(solsInTriangles[i][j], solsInTriangles[i1][j1])
                        if dominance_indicator == DOMINATED
                            finished = true
                            break
                        end
                        if finished
                            break
                        end
                    end
                    if finished
                        break
                    end
                end
                if finished
                    break
                end
            end

            if !finished
                push!(newList[i],solsInTriangles[i][j])
            end
        end
    end
    return newList
end

function descentNeighborhood(sol::Solution, consecutiveSet::Vector{PairOfSolution}, neighborhood::Neighborhood, prob::Problem, allSolsI::Vector{Solution}, solsInTriangles::Vector{Vector{Solution}};steepest=true,fastEvaluate=fastEvaluate)
    #println("descentNeighborhood")
    improved = false
    finished = false
    cpt = 0
    tEval = 0.
    tFun = 0.
    if neighborhood == KP_1_1
        if steepest
            bestKP,tEval,tFun = bestKP11(sol,consecutiveSet,prob,tEval,tFun,fastEvaluate=fastEvaluate)
        else
            bestKP = firstKP11(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
        end
    elseif neighborhood == KP_0_1
        if steepest
            bestKP = bestKP01(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
        else
            bestKP = firstKP01(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
        end
    elseif neighborhood == KP_1_2
        if steepest
            bestKP = bestKP12(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
        else
            bestKP = firstKP12(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
        end
    else
        @assert false "unknown neighborhood"
    end
    push!(allSolsI, Solution(copy(sol.x),copy(sol.y),sol.w))
    while bestKP[1] != 0
        cpt += 1
        sol = updateSolByKP(sol,bestKP,neighborhood,prob)
        push!(allSolsI, Solution(copy(sol.x),copy(sol.y),sol.w))
        # test if the solution should be added to solsInTriangles and if some old ones should be deleted
        indexTriangle = getTriangleIndex(solsInTriangles,sol)
        for i in 1:length(solsInTriangles[indexTriangle])
            if solsInTriangles[indexTriangle][i].y == sol.y # the solution is not new
                finished = true
                break
            end
        end

        if !finished
            push!(solsInTriangles[indexTriangle],sol)
            improved = true
            secureList(solsInTriangles)
        end

        if finished
            bestKP[1] = 0
        elseif neighborhood == KP_1_1
            if steepest
                bestKP,tEval, tFun = bestKP11(sol,consecutiveSet,prob,tEval,tFun,fastEvaluate=fastEvaluate)
            else
                bestKP = firstKP11(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
            end
        elseif neighborhood == KP_0_1
            if steepest
                bestKP = bestKP01(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
            else
                bestKP = firstKP01(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
            end
        elseif neighborhood == KP_1_2
            if steepest
                bestKP = bestKP12(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
            else
                bestKP = firstKP12(sol,consecutiveSet,prob,fastEvaluate=fastEvaluate)
            end
        else
            @assert false "Unknwon neighborhood"
        end
    end

    #println("tEval = $tEval")
    #println("tFun = $tFun")

    return sol, cpt, improved
end

function descent(sol::Solution, consecutiveSet::Vector{PairOfSolution}, neighborhoods::Vector{Neighborhood}, prob::Problem, allSolsI::Vector{Solution}, solsInTriangles::Vector{Vector{Solution}}; steepest = steepest, fastEvaluate = fastEvaluate)
    #println("descent")
    cptGlobal = zeros(Int,length(neighborhoods))
    improving = false
    for i in 1:length(neighborhoods)
        sol, cpt, improved = descentNeighborhood(sol,consecutiveSet,neighborhoods[i],prob,allSolsI,solsInTriangles,steepest=steepest,fastEvaluate=fastEvaluate)
        cptGlobal[i] += cpt
        improving = improving || improved
    end

    #println("   <cptGlobal = $cptGlobal>")

    return sol
end

function improveLowerBoundSet(lowerBoundSet::Vector{Solution},consecutiveSet::Vector{PairOfSolution},prob::Problem;steepest = true, fastEvaluate = true, neighborhoods = [KP_1_1])
    #println("steepest = $steepest, fastEval = $fastEvaluate, neighborhoods = $neighborhoods")
    if length(lowerBoundSet) == 0
        return lowerBoundSet, consecutiveSet
    else
        @assert length(lowerBoundSet) == length(consecutiveSet)+1 "length(lowerBoundSet) != length(consecutiveSet)+1"

        solsPrime = Vector{Solution}(undef,length(lowerBoundSet))
        allSols = Vector{Vector{Solution}}(undef,length(lowerBoundSet))
        solsInTriangles = Vector{Vector{Solution}}(undef,length(consecutiveSet))

        for i in 1:length(consecutiveSet)
            solsInTriangles[i] = [consecutiveSet[i].solL, consecutiveSet[i].solR]
        end

        for i in 1:length(lowerBoundSet)
            allSols[i] = Vector{Solution}()
        end

        improving = true
        while improving
            improving = false
            for i in 1:length(lowerBoundSet)
                solsPrime[i], improved = descent(lowerBoundSet[i],consecutiveSet, neighborhoods, prob, allSols[i], solsInTriangles, steepest = steepest, fastEvaluate = fastEvaluate)
                improving = improving || improved
            end

            # setting the structures for a new iteration
            solsClean = cleanList(solsInTriangles)
            newLowerBoundSet = Vector{Solution}()

            for i in 1:length(solsClean)
                push!(newLowerBoundSet,solsClean[i][1])
                for j in 3:length(solsClean[i])
                    push!(newLowerBoundSet,solsClean[i][j])
                end
                if i == length(solsClean)
                    push!(newLowerBoundSet,solsClean[i][2])
                end
            end

            sort!(newLowerBoundSet, by = sol -> sol.y[1])

            newConsecutiveSet = Vector{PairOfSolution}(undef,length(newLowerBoundSet)-1)
            for i in 1:length(consecutiveSet)
                consecutiveSet[i] = PairOfSolution(newLowerBoundSet[i],newLowerBoundSet[i+1])
            end
            # at this point we have constructed our new structures


        end

        #println("allSols = $allSols")

        #println("solsPrime = $solsPrime")

        #println("solsInTriangles = $solsInTriangles")
        #println("size2DList(solsInTriangles) = $(size2DList(solsInTriangles))")

        solsClean = cleanList(solsInTriangles)
        lb = nil(Solution)
        lb = cons(solsClean[length(solsClean)][2],lb)

        for i in 1:length(solsClean)
            print(solsClean[i])
        end

        lbVect = Vector{Solution}()

        for i in 1:length(solsClean)
            push!(lbVect,solsClean[i][1])
            for j in 3:length(solsClean[i])
                push!(lbVect,solsClean[i][j])
            end
            if i == length(solsClean)
                push!(lbVect,solsClean[i][2])
            end
        end

        println("lbVect = $lbVect")

        sort!(lbVect, by = sol -> sol.y[1])

        println("lbVect = $lbVect")

        lb = nil(Solution)
        lb = cons(lbVect[length(lbVect)],lb)

        for i in length(lbVect)-1:-1:1
            lb = cons(lbVect[i],lb)
        end

        nadirPoints = getNadirPoints(lb) #Jules

        println("nadirPoints = $nadirPoints")

        #println("solsClean = $solsClean")
        #println("length(lowerBoundSet) = $(length(lowerBoundSet))")
        #println("size2DList(solsClean) = $(size2DList(solsClean) - length(lowerBoundSet) + 2)")

        #plotRunBis("Inst200",lowerBoundSet,solsPrime, allSols, consecutiveSet)

        #plotRun3("Inst200",lowerBoundSet,solsPrime, allSols, consecutiveSet, solsClean)

        return lb, nadirPoints, size2DList(solsClean) - 2*length(lowerBoundSet) + 2, solsInTriangles
    end
end

function plotRun3(iname, sols, newSols, allSols, consecutiveSet, solsInTriangle)
    figure("Run Steepest",figsize=(6,6)) # Create a new figure
    title("Steepest-Descent-2OKP | " * iname)
    xlabel("z1")
    ylabel("z2")

    #println("consecutiveSet = $consecutiveSet")
    for i in 1:length(solsInTriangle)
        for j in 3:length(solsInTriangle[i])
            sol = solsInTriangle[i][j]

            plot(sol.y[1], sol.y[2], ms = "4", marker = "o", color = "black")
        end
    end
end

function plotRunBis(iname, sols, newSols, allSols, consecutiveSet)
    figure("Run Steepest",figsize=(6,6)) # Create a new figure
    #title("Steepest-Descent-2OKP | " * iname)
    xlabel("z1")
    ylabel("z2")

    x = Vector{Float64}(undef,length(sols))
    y = Vector{Float64}(undef,length(sols))
    newX = Vector{Float64}(undef,length(newSols))
    newY = Vector{Float64}(undef,length(newSols))
    for i in 1:length(sols)
        x[i] = sols[i].y[1]
        y[i] = sols[i].y[2]
    end

    for i in 1:length(newSols)
        newX[i] = newSols[i].y[1]
        newY[i] = newSols[i].y[2]
    end

    X = Vector{Vector{Float64}}(undef,length(allSols))
    Y = Vector{Vector{Float64}}(undef,length(allSols))
    for i in 1:length(allSols)
        X[i] = Vector{Float64}()
        Y[i] = Vector{Float64}()
        for j in 1:length(allSols[i])
            push!(X[i],allSols[i][j].y[1])
            push!(Y[i],allSols[i][j].y[2])
        end
    end

    XT = Vector{Float64}(undef,3*length(consecutiveSet))
    YT = Vector{Float64}(undef,3*length(consecutiveSet))
    #println("consecutiveSet = $consecutiveSet")
    for i in 1:length(consecutiveSet)
        pair = consecutiveSet[i]
        solL = pair.solL
        solR = pair.solR

        plot([solL.y[1], solL.y[1], solR.y[1], solL.y[1]], [solL.y[2], solR.y[2], solR.y[2], solL.y[2]], color = "orange", ls = "-", linewidth=1)

        XT[3*(i-1)+1] = solL.y[1]
        XT[3*(i-1)+2] = solL.y[1]
        XT[3*(i-1)+3] = solR.y[1]

        YT[3*(i-1)+1] = solL.y[2]
        YT[3*(i-1)+2] = solR.y[2]
        YT[3*(i-1)+3] = solR.y[2]
    end

    for i in 1:length(allSols)
        plot(X[i],Y[i],ls="-", marker = "o", ms = 3)
    end

    #plot(XT,YT,ls="-",marker="o",ms=3,color="blue",label="Triangles")
    plot(x,y,ls="",marker="o",ms=4,color="orange",label="LB")
    #plot(newX,newY,ls="",marker="x",ms=10,color="red",label="nouvelles")
    #legend(loc=4, fontsize ="small")
end

function plotRun(iname, sols, newSols)
    figure("Run",figsize=(6,6)) # Create a new figure
    #title("Steepest-Descent-2OKP | " * iname)
    xlabel("z1")
    ylabel("z2")

    x = Vector{Float64}(undef,length(sols))
    y = Vector{Float64}(undef,length(sols))
    newX = Vector{Float64}(undef,length(newSols))
    newY = Vector{Float64}(undef,length(newSols))
    for i in 1:length(sols)
        x[i] = sols[i].y[1]
        y[i] = sols[i].y[2]
    end

    for i in 1:length(newSols)
        newX[i] = newSols[i].y[1]
        newY[i] = newSols[i].y[2]
    end

    #plot(x,y,ls="",marker="+",ms=5,color="green",label="LB")
    #plot(newX,newY,ls="",marker="x",ms=5,color="red",label="nouvelles")
    #legend(loc=4, fontsize ="small")
end