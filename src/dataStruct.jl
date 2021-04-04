
@enum PRUNED optimality infeasibility dominance none

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

struct DualSetLinear
	A::Array{Float64, 2}
	b::Array{Float64, 1}
end