
"""
    jules's functions
"""

"""
    redifinition of data structures
"""
struct Problem
	nbVar::Int
	nbObj::Int
	profits::Array{Int, 2}
	weights::Array{Int, 1}
	maxWeight::Int
end

struct Solution
	x::Vector{Bool}
	y::Vector{Int}
	
	w::Int
end

struct PairOfSolution
	solL::Solution
	solR::Solution
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
"""
function getNadirPoints(LB::Vector{Solution})
    nbNadirPoints = length(LB) - 1
    nadirPoints = Vector{PairOfSolution}(undef,nbNadirPoints)
    for i in 1:nbNadirPoints
        nadirPoints[i] = PairOfSolution(LB[i],LB[i+1]) # we put the solutions by reference not by copy
    end
    return nadirPoints
end

function updateLowerBound(LB::Vector{Solution}, nadirPoints::Vector{PairOfSolution})

end