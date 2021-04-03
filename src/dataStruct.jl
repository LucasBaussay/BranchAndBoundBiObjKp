
@enum PRUNED optimality infeasibility dominance none

struct Problem
	nbVar::Int
	nbObj::Int
	profits::Array{Float64, 2}
	weights::Array{Float64,1}
	weightMax::Int
end

mutable struct Solution
	x::Array{Float64, 1}
	y::Array{Float64, 1}
	w::Float64
end

struct Assignment
	assign::Array{Float64, 1}
	profit::Float64
	weight::Float64
	assignEndIndex:: Int
end

struct PairOfSolution
	sol1::Solution
	sol2::Solution
end
