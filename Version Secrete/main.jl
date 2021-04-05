using JuMP, GLPK
using DataStructures

include("dataStruct.jl")

function addLowerBound!(lowerBound::T, sol::Solution) where T<:LinkedList{Solution}
	test = false
	while lowerBound != nil(Solution)
		if lowerBound.tail != nil(Solution) && !test && lowerBound.tail.head.y[1] > sol.y[1]
			lowerBound.tail = cons(sol, lowerBound.tail)
			test = true
		elseif lowerBound.tail == nil(Solution) && !test
			lowerBound.tail = cons(sol, nil(Solution))
			test = true
			@assert true "Normalement je ne passe pas là"
		end
		lowerBound = lowerBound.tail
	end
end

function getNadirPoints(LB::T) where T<:LinkedList{Solution}
	nbNadirPoints = length(LB) - 1
    nadirPoints = Vector{PairOfSolution}(undef,nbNadirPoints)
    
    iter = 1
    while LB.tail != nil(Solution)
    	nadirPoints[iter] = PairOfSolution(LB.head, LB.tail.head)
    	LB = LB.tail
    	iter += 1
    end
    
    return nadirPoints
end

function solve1OKP(prob::Problem, assignment::Assignment)
	#Passer Lambda en paramètre
	@assert prob.nbObj == 1 "solve1OKP only supports one objective function"

	model = Model(GLPK.Optimizer)
	x = @variable(model, x[1:(prob.nbVar-assignment.assignEndIndex)], Bin)
	@constraint(model, Weights, sum(x .* prob.weights[(assignment.assignEndIndex+1):end]) + assignment.weight <= prob.maxWeight)
	@objective(model, Max, sum(x .* prob.profits[1, (assignment.assignEndIndex+1):end]))

	optimize!(model)

	X = append!(assignment.assign[1:assignment.assignEndIndex],(Float64).(value.(x)))
	
	termStatus = termination_status(model)

	if termStatus == MOI.OPTIMAL
		return Solution(X, [sum(X .* prob.profits[1,1:end])], sum(X .* prob.weights)), true
	elseif termStatus == MOI.INFEASIBLE
		return Solution(), false
	else
		error("Cette Fois c'est MOI qui renvoie de la merde")
	end
end

function weightedRelax(prob::Problem, λ::Vector{Float64})::Problem

	obj = zeros(Float64, 1, prob.nbVar)
	for iter = 1:prob.nbVar
		obj[iter] = sum( λ .* prob.profits[1:end, iter])
	end

	return Problem(
			prob.nbVar,
			1,
			obj,
			prob.weights,
			prob.maxWeight)
end

function evaluate(prob::Problem, x::Vector{Float64})::Solution
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

function dichoSearch(prob::Problem, assignment::Assignment; M::Float64 = 1000., verbose::Bool = false)

    lowerBound = nil(Solution)
	toStudy = Vector{PairOfSolution}()

    λ = [1, M]
    leftSol, testLeft = solve1OKP(weightedRelax(prob, λ), assignment)
	leftSol = evaluate(prob, leftSol.x)
    λ = [M, 1]
    rightSol, testRight = solve1OKP(weightedRelax(prob, λ), assignment)
	rightSol = evaluate(prob, rightSol.x)
	
	@assert testLeft && testRight "Bah écoute, on choppe des solutions infaisables, de mieux en mieux"

	verbose && println("Two found solutions : $leftSol && $rightSol")

    # In the event that the two solutions found are identical,
    # this solution is optimal and, thus, added to the lower bound
    if leftSol.y == rightSol.y
        lowerBound = cons(leftSol, lowerBound)
		verbose && println("The two solutions are identical")
    else
		# lowerBound doit être trié
		# A REVOIR - C'est fait :P
		lowerBound = cons(leftSol, cons(rightSol, nil(Solution)))

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
	        midSol, testMid = solve1OKP(weightedRelax(prob, λ), assignment)
			midSol = evaluate(prob, midSol.x)
			
			@assert testMid "Bah écoute, on choppe des solutions infaisables, de mieux en mieux"

			# if the solution dominates one of the other, it's added to the LB
			if sum(λ .* midSol.y) > sum(λ .* rightSol.y)
				addLowerBound!(lowerBound, midSol)
				push!(toStudy, PairOfSolution(leftSol, midSol))
				push!(toStudy, PairOfSolution(midSol, rightSol))
			end
        end
    end

	return lowerBound

end

function computeUpperBound(lowerBound::T) where T<:LinkedList{Solution}

	nbConstraints = length(lowerBound) - 1 + 2

	A = zeros(Float64, nbConstraints, 2)
	b = zeros(Float64, nbConstraints)

	max1 = 0
	max2 = 0
	
	iter = 1
	while lowerBound.tail != nil(Solution)

		pair = PairOfSolution(lowerBound.head, lowerBound.tail.head)

		A[iter, 1] = pair.solL.y[2] - pair.solR.y[2]
		A[iter, 2] = pair.solR.y[1] - pair.solL.y[1]

		b[iter] = (pair.solL.y[2] - pair.solR.y[2]) * pair.solR.y[1] + (pair.solR.y[1] - pair.solL.y[1]) * pair.solR.y[2]

		if pair.solL.y[1] > max1
			max1 = pair.solL.y[1]
		end
		if pair.solR.y[1] > max1
			max1 = pair.solR.y[1]
		end
		if pair.solL.y[2] > max2
			max2 = pair.solL.y[2]
		end
		if pair.solR.y[2] > max2
			max2 = pair.solR.y[2]
		end
		
		lowerBound = lowerBound.tail
		iter += 1

	end

	A[nbConstraints-1, 1] = 1
	A[nbConstraints-1, 2] = 0

	b[nbConstraints-1] = max1

	A[nbConstraints, 1] = 0
	A[nbConstraints, 2] = 1

	b[nbConstraints] = max2

	return DualSet(A, b)

end

function dominate(sol1::Solution, sol2::Solution)
	return (sol1.y[1] >= sol2.y[1] && sol1.y[2] > sol2.y[2]) || (sol1.y[1] > sol2.y[1] && sol1.y[2] >= sol2.y[2])
end

function dominate(sol::Solution, pair::PairOfSolution)
	return sol.y[1] > pair.solL.y[1] && sol.y[2] > pair.solR.y[2]
end

function updateLowerBound!(lowerBound::T, nadirPoints::Vector{PairOfSolution}, listOfPoint::LinkedList{Solution}) where T<:LinkedList{Solution}

	for sol in listOfPoint
		
		studiedLowerBound = lowerBound
		
		testFirstIn = false
		testFirstOut = false
		
		eltBeforeIn = nil(Solution)
		eltAfterOut = nil(Solution)
		
		
		while studiedLowerBound != nil(Solution) && !testFirstOut && !(testFirstIn && studiedLowerBound.head.y[1] > sol.y[1])
			
			if !testFirstIn && studiedLowerBound.tail != nil(Solution) && dominate(sol, studiedLowerBound.tail.head)
				testFirstIn = true
				eltBeforeIn = studiedLowerBound
			elseif testFirstIn && !dominate(sol, studiedLowerBound.head)
				testFirstOut = true
				eltAfterOut = studiedLowerBound
			end
			
			studiedLowerBound = studiedLowerBound.tail
		end
		
		if testFirstIn && testFirstOut
			eltBeforeIn.tail = cons(sol, eltAfterOut)
			
			indBeforeIn = 0
			indAfterOut = 0
			iterNadir = 1
			while iterNadir <= length(nadirPoints) && indAfterOut == 0
				if indBeforeIn == 0
					if nadirPoints[iterNadir].solL.y == eltBeforeIn.head.y
						indBeforeIn = iterNadir - 1
					end
				elseif nadirPoints[iterNadir].solR.y == eltAfterOut.head.y
					indAfterOut = iterNadir + 1
				end
			end
			nadirPoints = append!(push!(nadirPoints[1:indBeforeIn], PairOfSolution(nadirPoints[indBeforeIn+1].solL, sol), PairOfSolution(sol, nadirPoints[indAfterOut-1].solR)), nadirPoints[indAfterOut:end])
			
		elseif testFirstIn || testFirstOut
			error("Putain de Merde")
		else
			
			studiedLowerBound = lowerBound
			testDominNadir = false
			pairDomin = PairOfSolution()
			
			while studiedLowerBound.tail != nil(Solution) && !testDominNadir
				pair = PairOfSolution(studiedLowerBound.head, studiedLowerBound.tail.head)
				if dominate(sol, pair)
					studiedLowerBound.tail = cons(sol, studiedLowerBound.tail)
					testDominNadir = true
					pairDomin = pair
				end
				
				studiedLowerBound = studiedLowerBound.tail
			end
			
			if testDominNadir
				iterPair = 1
				while iterPair <= length(nadirPoints) && nadirPoints[iterPair] != pairDomin
					iterPair += 1
				end
				nadirPoints = append!(push!(nadirPoints[1:(iterPair-1)], PairOfSolution(pairDomin.solL, sol), PairOfSolution(sol, pairDomin.solR)), nadirPoints[(iterPair+1):end])
			end
			
		end
		
	end
	
	return nadirPoints

end



function pruningTest(lowerBound::T, listPointsNadir::Vector{PairOfSolution}, subUpperBound::DualSet, verbose = false) where T<:LinkedList{Solution}

	lowerBound == nil(Solution) && return infeasibility, Vector{PairOfSolution}()
	lowerBound.tail == nil(Solution) && return optimality, Vector{PairOfSolution}()
	
	
	domi = true
	MyDude = Vector{PairOfSolution}()

	for pairNadir in listPointsNadir

		verbose && println("pairNadir : $pairNadir")
		nadir = [pairNadir.solL.y[1], pairNadir.solR.y[2]]
		nadirA = subUpperBound.A * nadir
		
		println("$nadirA - $(subUpperBound.b)")

		iter = 1
		
		nonDomiNadir = false
	
		while iter <= length(subUpperBound.b) && !nonDomiNadir
			
			nonDomiNadir = nonDomiNadir || nadirA[iter] > subUpperBound.b[iter]
			
			iter += 1
		end
		if !nonDomiNadir
			push!(MyDude, pairNadir)
		end
		
		domi = domi && nonDomiNadir
	end

	if domi
		return dominance, Vector{PairOfSolution}()
	else
		return none, MyDude
	end

end

function addVarAssignment!(assignment::Assignment, prob::Problem)
	assignment.assignEndIndex += 1
	assignment.assign[assignment.assignEndIndex] = 1

	assignment.weight += prob.weights[assignment.assignEndIndex]
	
	assignment.profit += prob.profits[1:end, assignment.assignEndIndex]
	assignment
end

function removeVarAssignment!(assignment::Assignment, prob::Problem, testAddVar::Bool)
	if !testAddVar
		assignment.assignEndIndex += 1
	else
		assignment.weight -= prob.weights[assignment.assignEndIndex]
		assignment.profit -= prob.profits[1:end, assignment.assignEndIndex]
	end

	assignment.assign[assignment.assignEndIndex] = 0
end

function returnParentAssignment!(assignment::Assignment, prob::Problem)
	assignment.assign[assignment.assignEndIndex] = -1
	assignment.assignEndIndex -= 1
end

function branchAndBound!(lowerBound::T, prob::Problem, assignment::Assignment, nadirPointsToStudy::Vector{PairOfSolution}; M::Float64 = 1000., num::Int = 1) where T<:LinkedList{Solution}
	
	println(num)
	
	subLowerBound = dichoSearch(prob, assignment, M = M) #Nathalie
	subUpperBound = computeUpperBound(lowerBound) #Lucas

	prunedType, newNadirPoints = pruningTest(subLowerBound, nadirPointsToStudy, subUpperBound) #Nathalie

	println(prunedType)

	if prunedType == optimality || prunedType == none
		newNadirPoints = updateLowerBound!(lowerBound, newNadirPoints, subLowerBound) #Lucas/Jules/Nathalie
	end

	if prunedType == none && assignment.assignEndIndex < prob.nbVar
		testAddVar = (assignment.weight + prob.weights[assignment.assignEndIndex + 1] <= prob.maxWeight)

		if testAddVar
			addVarAssignment!(assignment, prob) #Lucas
			branchAndBound!(lowerBound, prob, assignment, newNadirPoints, M = M, num = num + 1) #Lucas
		end
		removeVarAssignment!(assignment, prob, testAddVar) #Lucas
		branchAndBound!(lowerBound, prob, assignment, newNadirPoints, M = M, num = num + 1) #Lucas

		returnParentAssignment!(assignment, prob) #Lucas
	end

end

function main(prob::Problem; withLinear::Bool = false, M::Float64 = 1000.)

	@assert prob.nbObj == 2 "This Branch and Bound supports only a Bio-objective problem"

	assignment = Assignment(prob) #Nathalie



	lowerBound = dichoSearch(prob, assignment, M = M) #Nathalie
	nadirPoints = getNadirPoints(lowerBound) #Jules

	branchAndBound!(lowerBound, prob, assignment, nadirPoints, M = M) #Lucass

	return lowerBound

end

function main(fname::String; withLinear::Bool = false, M::Float64 = 1000.)
	prob = Problem(fname)
	return main(prob, withLinear = withLinear, M = M)
end

function main(;withLinear::Bool = false, M::Float64 = 1000.)
	prob = Problem()
	return main(prob, withLinear = withLinear, M = M)
end
