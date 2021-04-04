
using JuMP, GLPK, MathOptInterface, PyPlot

const MOI = MathOptInterface

"""
    jules's functions
"""

"""
    redifinition of data structures and empty functions
"""
@enum Pruned optimality infeasibility dominance none

struct Problem
	nbVar::Int
	nbObj::Int
	profits::Array{Float64, 2}
	weights::Array{Float64,1}
	weightMax::Int
end

function Problem()
	return Problem(
					6,
					2,
					[11 2 8 10 9 1; 2 7 8 4 1 3],
					[4, 4, 6, 4, 3, 2],
					11
				)
end

mutable struct Solution
	x::Vector{Float64}
	y::Vector{Float64}
	w::Float64
end

struct Assignment
	assign::Vector{Float64}
	profit::Vector{Float64}
	weight::Float64
	assignEndIndex::Int
end

function Assignment()
	return Assignment(
			Vector{Float64}(),
			Vector{Float64}(),
			0,
			0)
end

function Assignment(prob::Problem)
	return(
		ones(Float64, prob.nbVar)*-1,
		zeros(Float64, prob.nbObj),
		0,
		0)
end

struct PairOfSolutions
	solL::Solution
	solR::Solution
end

struct DualSet
	A::Array{Float64, 2}
	b::Array{Float64, 1}
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

function solve1OKP(prob::Problem, assignment::Assignment)
    @assert prob.nbObj == 1 "solve1OKP only supports one objective function"

	model = Model(GLPK.Optimizer)
	x = @variable(model, x[1:(prob.nbVar-assignment.assignEndIndex)], Bin)
	@constraint(model, Weights, sum(x .* prob.weights[assignment.assignEndIndex+1:end]) + assignment.weight <= prob.maxWeight)
	@objective(model, Max, sum(x .* prob.profits[1, assignment.assignEndIndex+1:end]) + assignment.profit)
    
	optimize!(model)

    X = append!(assignment.assign[1:assignment.assignEndIndex],(Float64).(value.(x)))

    if MOI.get(GLPK.Optimizer, TerminationStatus()) == success
        return Solution(X, [objective_value(model)], sum(value.(x) .* prob.weights)), true
    else
        return Solution(X, [objective_value(model)], sum(value.(x) .* prob.weights)), false
    end
end

#=
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
        pop!(nadirPoints)
    end
end
=#