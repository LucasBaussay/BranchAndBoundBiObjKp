
"""
    jules's functions
"""

"""
    redifinition of data structures and empty functions
"""
struct Problem
	nbVar::Int
	nbObj::Int
	profits::Array{Float64, 2}
	weights::Array{Float64, 1}
	maxWeight::Float64
end

struct Solution
	x::Vector{Bool}
	y::Vector{Float64}
	
	w::Float64
end

struct PairOfSolution
	solL::Solution
	solR::Solution
end

function isLDominatingR(yL::Vector{Float64}, yR::Vector{Float64})
    @assert length(yL) == 2 && length(yR) == 2 "We only support two objectives"
    solLDominates = true
    solRDominates = true

    for i in 1:2
        if yL[i] < yR[i]
            return false
        end
        if yL[i] > yR[i]
            solRDominates = false
        end
    end

    if solLDominates && !solRDominates
        return true
    else
        return false
    end
end

@enum Dominance DOMINATES DOMINATED NOTHING

function whichDominates(yL::Vector{Float64}, yR::Vector{Float64})
    @assert length(yL) == 2 && length(yR) == 2 "We only support two objectives"
    solLDominates = true
    solRDominates = true

    for i in 1:2
        if yL[i] < yR[i]
            solLDominates = false
        end
        if yL[i] > yR[i]
            solRDominates = false
        end
    end

    if solLDominates && solRDominates # the Δs are the same
        @assert false "Those sols shouldn't be the same"
    elseif !solLDominates && !solRDominates # we can't compare the two kp-exchanges
        return NOTHING
    elseif solLDominates
        return DOMINATES
    else
        return DOMINATED
    end
end

function replaceSolByOtherSol(oldSol::Solution, newSol::Solution)
    oldSol.x = newSol.x
    oldSol.y = newSol.y
    oldSol.w = oldSol.w
end

"""
    mainJules(;verbose=false)

    @verbose : print comments if verbose equals true
"""
function mainJules(;verbose=false)
    verbose && println("Welcome to the Jules's section :\n\nFunctions available :\n   getNadirPoints(LB::Vector{Solution})\n   updateLowerBound(LB::Vector{Solution}, nadirPoints::Vector{PairOfSolution})\n")
end

"""
    getNadirPoints(LB::Vector{Solution})

    returns a the nadir points (ordered by first objective value) in a vector of pair of solutions (Vector{PairOfSolution})

    IMPORTANT : the solutions inside each a PairOfSolution are the solutions of the LB parameter. We didn't copy them but directly paste their reference.
"""
function getNadirPoints(LB::Vector{Solution})
    nbNadirPoints = length(LB) - 1
    nadirPoints = Vector{PairOfSolution}(undef,nbNadirPoints)
    for i in 1:nbNadirPoints
        nadirPoints[i] = PairOfSolution(LB[i],LB[i+1]) # we put the solutions by reference not by copy
    end
    return nadirPoints
end

"""
    updateLowerBound(LB::Vector{Solution}, nadirPoints::Vector{PairOfSolution}, subLB::Vector{Solution})

    update the LB parameter and returns a new nadirPoints vector.

    @LB : current best LB for the main problem (solutions oredered by first objective value)
    @nadirPoints : nadirPoints of the LB
    @subLB : 
"""
function updateLowerBound(LB::Vector{Solution}, nadirPoints::Vector{PairOfSolution}, subLB::Vector{Solution})
    # check if solutions of subLB dominate solutions of LB
    indFirstDominated = zeros(Int,length(subLB))
    indFirstNotDominated = zeros(Int,length(subLB))
    for i in 1:length(subLB)
        deb = 1
        if i > 1 # we skip the first part of the LB to optimize the process
            if indFirstNotDominated[i-1] > 0
                deb = indFirstNotDominated[i-1]
            end
        end
        for j in deb:length(LB)
            if isLDominatingR(subLB[i].y,LB[j].y)
                if indFirstDominated[i] == 0
                    indFirstDominated[i] = j
                end
            elseif indFirstDominated[i] > 0 && j == length(LB)
                indFirstNotDominated[i] = j+1
                break
            elseif indFirstDominated[i] > 0 && indFirstNotDominated[i] == 0
                indFirstNotDominated[i] = j
                break
            end
        end
    end
    # modify LB and nadirPoints in consequence
    decay = 0
    itStart = 0
    itEnd = 0
    for i in 1:length(subLB)
        if indFirstDominated[i] > 0
            replaceSolByOtherSol(LB[indFirstDominated[i]-decay],subLB[i])
        elseif indFirstNotDominated[i] - indFirstDominated[i] > 1
            itStart = indFirstDominated[i]+1 # 3
            itEnd = indFirstNotDominated[i] # 5
            for k in 1:(indFirstNotDominated[i]-indFirstDominated[i]-1) # 1:2 (deux it)
                if itEnd <= length(LB)
                    replaceSolByOtherSol(LB[itStart-decay], LB[itEnd]) # LB[3] <- LB[5], LB[4] <- LB[6]
                    itStart += 1
                    itEnd += 1
                end
                decay += 1
            end
        end
    end

    for i in 1:decay
        pop!(LB)
    end
end
