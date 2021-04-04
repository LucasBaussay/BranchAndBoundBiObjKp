function computeUpperBound(lowerBound::Vector{Solution})

	nbConstraints = length(lowerBound) - 1 + 2

	A = zeros(Float64, nbConstraints, 2)
	b = zeros(Float64, nbConstraints)

	max1 = 0
	max2 = 0

	for iter = 1:(nbConstraints-2)

		pair = PairOfSolution(lowerBound[iter], lowerBound[iter+1])
		
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

	end

	A[nbConstraints-1, 1] = 1
	A[nbConstraints-1, 2] = 0

	b[nbConstraints-1] = max1

	A[nbConstraints, 1] = 0
	A[nbConstraints, 2] = 1

	b[nbConstraints] = max2

	return DualSet(A, b)

end

function addVarAssignment!(assignment::Assignment, prob::Problem)
	assignment.assignEndIndex += 1
	assignment.assign[assignment.assignEndIndex] = 1
	
	assignment.weight += prob.weights[assignment.assignEndIndex]
	assignment.profit += prob.profits[1:end, assignment.assignEndIndex]
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

function branchAndBound!(lowerBound::Vector{Solution}, assignment::Assignment, nadirPointsToStudy::Vector{PairOfSolution}; M::Int = 1000)
	
	subLowerBound = dichoSearch(prob, assignment, M = M) #Nathalie
	subUpperBound = computeUpperBound(lowerBound) #Lucas
	
	prunedType, newNadirPoints = pruningTest(subUpperBound, nadirPointsTuStudy) #Nathalie
	
	if prunedType == oprimality || prunedType == none
		updateLowerBound!(lowerBound, newNadirPoints, subLowerBound) #Jules
	end
	
	if prunedType == none && assignment.assignEndIndex < prob.nbVar
		testAddVar = (assignment.weight + prob.weights[assignment.assignEndIndex] <= prob.maxWeight)
		
		if testAddVar
			addVarAssignment!(assignment, prob) #Lucas
			branchAndABound!(lowerBound, assignment, newNadirPoints, M = M) #Lucas
		end
		removeVarAssignment!(assignment, prob, testAddVar) #Lucas
		branchAndBound!(lowerBound, newNadirPoints, M = M) #Lucas
		
		returnParentAssignment!(assignment, prob) #Lucas
	end
	
	

end

function main(prob::Problem; withLinear::Bool = false, M::Int = 1000)

	@assert prob.nbObj == 2 "This Branch and Bound supports only a Bio-objective problem"

	assignment = Assignment(prob) #Nathalie



	lowerBound = dichoSearch(prob, assignment, M = M) #Nathalie
	nadirPoints = getNadirPoint(lowerBound) #Jules

	branchAndBound!(lowerBound, assignment, nadirPoints, M = M) #Lucass

	return lowerBound

end
