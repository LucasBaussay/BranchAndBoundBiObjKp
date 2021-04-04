
using JuMP, GLPK, MathOptInterface, PyPlot

const MOI = MathOptInterface

include("dataStruct.jl")

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
    solve1OKP(prob::Problem, assignment::Assignment)

    returns a tuple (Solution,Bool). The boolean is false if the solver couldn't solve the problem, and true otherwise.

    usage :
    > sol, optimal_status = solve1OKP(myProb, myAssignment)
    > if !optimal_status
            return empty_sol() // for exemple
      else
            return sol
      end
"""
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
