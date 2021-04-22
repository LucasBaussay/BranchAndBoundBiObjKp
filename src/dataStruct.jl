
@enum Pruned optimality infeasibility dominance none

struct Problem
	nbVar::Int
	nbObj::Int
	profits::Array{Float64, 2}
	weights::Array{Float64,1}
	maxWeight::Float64
	
	isInteger::Bool
end

function Problem()
	return Problem(
					6,
					2,
					[11 2 8 10 9 1; 2 7 8 4 1 3],
					[4, 4, 6, 4, 3, 2],
					11, 
					true
				)
end

function Problem(nbVar::Int, nbObj::Int, profits::Array{Float64, 2}, weights::Array{Float64, 1}, maxWeight::Float64)
	return Problem(
					nbVar,
					nbObj,
					profits,
					weights,
					maxWeight,
					false
					)
end

function Problem(fname::String)
	f = open(fname)

	nbVar = parse(Int, readline(f))
	nbObj = parse(Int, readline(f))
	nbConst = parse(Int, readline(f))

	@assert nbConst == 1 "We only care about mono-dimensional Knapsack"

	profits = zeros(Float64, nbObj, nbVar)
	weights = zeros(Float64, nbVar)

	for iter = 1:nbObj
		profits[iter, 1:end] = parse.(Float64, split(readline(f)))
	end

	weights = parse.(Float64, split(readline(f)))

	maxWeight = parse(Float64, readline(f))

	close(f)

	return Problem(
					nbVar,
					nbObj,
					profits,
					weights,
					maxWeight
					)
end

mutable struct Solution
	x::Vector{Float64}
	y::Vector{Float64}
	w::Float64
	
	isBinary::Bool
end

mutable struct Assignment
	assign::Vector{Int}
	profit::Vector{Float64}
	weight::Float64
	assignEndIndex::Int
end

function Assignment()
	return Assignment(
			Vector{Int}(),
			Vector{Float64}(),
			0,
			0)
end

function Assignment(prob::Problem)
	return Assignment(
		ones(Int, prob.nbVar)*-1,
		zeros(Float64, prob.nbObj),
		0,
		0)
end

struct PairOfSolution
	solL::Solution
	solR::Solution
end

struct DualSet
	A::Array{Float64, 2}
	b::Array{Float64, 1}
end

mutable struct Compteur
	value::Int
end

function Compteur()
	return Compteur(0)
end

function Solution()
	return Solution(
			Vector{Float64}(), 
			Vector{Float64}(),
			0.,
			false
			)
end

function Solution(x::Vector{Float64}, y::Vector{Float64}, w::Float64)
	return Solution(x,
					y,
					w,
					true)
end

function PairOfSolution()
	return PairOfSolution(
				Solution(),
				Solution()
			)
end

import Base.:(==), Base.:(!=)

Base.:(==)(sol1::Solution, sol2::Solution) = sol1.y == sol2.y
Base.:(==)(pair1::PairOfSolution, pair2::PairOfSolution) = pair1.solL == pair2.solL && pair1.solR == pair2.solR

Base.:(!=)(pair1::PairOfSolution, pair2::PairOfSolution) = !(pair1 == pair2)
