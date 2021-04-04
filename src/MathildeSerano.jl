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
	toStudy = Vector{PairOfSolutions}()

    λ = [1, M]
    leftSol = solveKP(weightedRelax(prob, λ))
	leftSol = evaluate(prob, leftSol.x)
    λ = [M, 1]
    rightSol = solveKP(weightedRelax(prob, λ))
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

        push!(toStudy, PairOfSolutions(leftSol, rightSol))

		verbose && println("toStudy origin: $toStudy")

		#while some pairs of solutions are found by the dichotomic search, we keep going
        while toStudy != []

			currPair = pop!(toStudy)
			verbose && println("toStudy : $toStudy")
			leftSol = currPair.solL
			rightSol = currPair.solR

			verbose && println("Found solutions : $leftSol && $rightSol")

	        λ = [leftSol.y[2] - rightSol.y[2], rightSol.y[1] - leftSol.y[1]]
	        midSol = solveKP(weightedRelax(prob, λ))
			midSol = evaluate(prob, midSol.x)

			# if the solution dominates one of the other, it's added to the LB
			if sum(λ .* midSol.y) > sum(λ .* rightSol.y)
				push!(lowerBound, midSol)
				push!(toStudy, PairOfSolutions(leftSol, midSol))
				push!(toStudy, PairOfSolutions(midSol, rightSol))
			end
        end
    end

	return lowerBound

end

function pruningTest(lengthSubLB::Int, listPointsNadir::Vector{PairOfSolution}, subUpperBound::DualSet)

	lengthSubLB == 1 && return optimality, Vector{PairOfSolutions}()
	lengthSubLB == 0 && return infeasibility, Vector{PairOfSolutions}()

	for pairNadir in listPointsNadir

		nadir = pairNadir.solL.y[1], pairNadir.solR.y[2]
		nadirA = nadir * subUpperBound.A

		# domi = true
		noneDomiNadirs = Vector{PairOfSolutions}()

		while iter <= length(subUpperBound) # && domi = true
			if nadirA[iter] <= subUpperBound.b[iter]
				# domi = false
				push!(noneDomiNadirs, pairNadir)
			end
			iter += 1

		end
	end

	if domi == true
		return dominance, Vector{PairOfSolutions}()
	else
		return none, noneDomiNadirs
	end

end
