using JuMP, GLPK

include("dataStruct.jl")

function weightedRelax(prob::Problem, λ::Vector{Float64})

	obj = zeros(Float64, 1, prob.nbVar)
	for iter = 1:prob.nbVar
		obj[iter] = sum( λ .* prob.profits[1:end, iter])
	end

	return Problem(
			prob.nbVar,
			1,
			obj,
			prob.weights,
			prob.weightMax)
end

#temporaire c'est pour tester frère
function solve1OKP(prob::Problem)

	model = Model(GLPK.Optimizer)
	x = @variable(model, x[1:prob.nbVar], Bin)
	@constraint(model, Weights, sum(x .* prob.weights) <= prob.weightMax)
	@objective(model, Max, sum(x .* prob.profits[1, 1:end]))

	optimize!(model)

	return Solution(Float64.(value.(x)), [objective_value(model)], sum(value.(x) .* prob.weights))
end

function evaluate(prob::Problem, x::Vector{Float64})
    y = zeros(Float64, prob.nbObj)
    weight = 0

    for iterObj = 1:prob.nbObj
        for iter = 1:prob.nbVar
            if x[iter] > 0
                y[iterObj] += prob.profits[iterObj, iter] * x[iter]
				if iterObj == 1
                	weight += prob.weights[iter] * x[iter]
				end
            end
        end
    end

    return Solution(x, y, weight)
end

function dichoSearch(prob::Problem, assignment::Assignment = Assignment(), M = 1000., verbose = true)

    # @assert

    assignment.assign == [] && (assignment = Assignment(prob))
    lowerBound = Vector{Solution}()
	toStudy = Vector{PairOfSolution}()

    λ = [1, M]
    leftSol = solve1OKP(weightedRelax(prob, λ))
	leftSol = evaluate(prob, leftSol.x)
    λ = [M, 1]
    rightSol = solve1OKP(weightedRelax(prob, λ))
	rightSol = evaluate(prob, rightSol.x)

	verbose && println("Two found solutions : $leftSol && $rightSol")

    # In the event that the two solutions found are identical,
    # this solution is optimal and, thus, added to the lower bound
    if leftSol.y == rightSol.y
        lowerBound = [leftSol]
		verbose && println("The two solutions are identical")
    else
		# lowerBound doit être trié
		# A REVOIR
		push!(lowerBound, leftSol)
		push!(lowerBound, rightSol)

		verbose && println("LB = $lowerBound")

        push!(toStudy, PairOfSolution(leftSol, rightSol))

		verbose && println("toStudy origin: $toStudy")

		#while some pairs of solutions are found by the dichotomic search, we keep going
        while toStudy != []

			currPair = pop!(toStudy)
			verbose && println("toStudy : $toStudy")
			leftSol = currPair.solL
			rightSol = currPair.solR

			verbose && println("Found solutions : $leftSol && $rightSol")

	        λ = [leftSol.y[2] - rightSol.y[2], rightSol.y[1] - leftSol.y[1]]
	        midSol = solve1OKP(weightedRelax(prob, λ))
			midSol = evaluate(prob, midSol.x)

			# if the solution dominates one of the other, it's added to the LB
			if sum(λ .* midSol.y) > sum(λ .* rightSol.y)
				push!(lowerBound, midSol)
				push!(toStudy, PairOfSolution(leftSol, midSol))
				push!(toStudy, PairOfSolution(midSol, rightSol))
			end
        end

		sort!(lowerBound, by=x->x.y[1])
    end

	return lowerBound

end

function pruningTest(lengthSubLB::Int, listPointsNadir::Vector{PairOfSolution}, subUpperBound::DualSet, verbose = true)

	lengthSubLB == 1 && return optimality, Vector{PairOfSolution}()
	lengthSubLB == 0 && return infeasibility, Vector{PairOfSolution}()

	domi = true
	noneDomiNadirs = Vector{PairOfSolution}()

	for pairNadir in listPointsNadir

		verbose && println("pairNadir : $pairNadir")
		nadir = [pairNadir.solL.y[1], pairNadir.solR.y[2]]
		nadirA = subUpperBound.A * nadir

		iter = 1
		first = true

		while iter <= length(subUpperBound.b) && (domi == true || first == true)
			if nadirA[iter] <= subUpperBound.b[iter]
				if first
					push!(noneDomiNadirs, pairNadir)
					first = false
				end
				domi = false
			end
			iter += 1
		end
	end

	if domi == true
		return dominance, Vector{PairOfSolution}()
	else
		return none, noneDomiNadirs
	end

end
