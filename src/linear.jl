using DataStructures, Random
include("dataStruct.jl")

function addLowerBound!(lowerBound::T, sol::Solution) where T<:LinkedList{Solution}
	test = false
	
	while lowerBound.tail != nil(Solution) && !test
	
		if lowerBound.tail.head.y[1] > sol.y[1]
			lowerBound.tail = cons(sol, lowerBound.tail)
			test = true
		end
		
		lowerBound = lowerBound.tail
		
	end
	
	if !test
		lowerBound.tail = cons(sol, nil(Solution))
		
		@assert false "Normalement je ne passe jamais par là\n $(sol.y)"
	end
end

function merge(listNadirA::Vector{PairOfSolution}, listNadirB::Vector{PairOfSolution})
	lengthA = length(listNadirA)
	lengthB = length(listNadirB)
	
	iterA = 1
	iterB = 1
	iter = 1
	
	finalList = Vector{PairOfSolution}(undef, lengthA + lengthB)
	
	while iterA <= lengthA && iterB <= lengthB
		if listNadirA[iterA].solL.y[1] < listNadirB[iterB].solL.y[1]
			finalList[iter] = listNadirA[iterA]
			iterA += 1
		else
			finalList[iter] = listNadirB[iterB]
			iterB += 1
		end
		iter += 1
	end
	
	if iterA > lengthA
		finalList[iter:end] = listNadirB[iterB:end]
	else
		finalList[iter:end] = listNadirA[iterA:end]
	end
	
	return finalList
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

function solve1OKPLinear(prob::Problem, λ::Vector{Float64}, assignment::Assignment)
	
	utilities = [sum(λ .* prob.profits[1:end, iter]) for iter = (assignment.assignEndIndex+1):prob.nbVar] ./ prob.weights[(assignment.assignEndIndex+1):end]
	permList = sortperm(utilities, rev = true)
	
	subLength = length(utilities)
	
	weight = assignment.weight
	profit = assignment.profit[1:end]
	x = Float64.(assignment.assign[1:end])
	for iter = 1:subLength
		permList[iter] += assignment.assignEndIndex
		x[assignment.assignEndIndex + iter] = 0.
	end
	
	iter = 1
	isBinary = false
	while iter <= subLength && weight + prob.weights[permList[iter]] <= prob.maxWeight 
		weight += prob.weights[permList[iter]]
		x[permList[iter]] = 1.
		profit += prob.profits[1:end, permList[iter]]
		
		iter += 1
	end
	
	if iter <= subLength && (prob.maxWeight - weight) > 0.
		x[permList[iter]] = (prob.maxWeight - weight) / prob.weights[permList[iter]]
		
		profit += x[permList[iter]] * prob.profits[1:end, permList[iter]]
		weight += x[permList[iter]] * prob.weights[permList[iter]]
		
	else
		isBinary = true
	end
	
	return Solution(x, profit, weight, isBinary)
		
	
end

function dichoSearch(prob::Problem, assignment::Assignment; pointsAlreadyCompute = nil(Solution), M::Float64 = 1000., compteur = nothing) where T <: LinkedList{Solution}

    lowerBound = copy(pointsAlreadyCompute)
	toStudy = Vector{PairOfSolution}()

    λ = [1, M]
	leftSol = solve1OKPLinear(prob, λ, assignment)
	
    λ = [M, 1]
	rightSol = solve1OKPLinear(prob, λ, assignment)

    if round.(leftSol.y, digits = 7) == round.(rightSol.y, digits = 7)
        lowerBound = cons(leftSol, lowerBound)
    else
		# lowerBound doit être trié
		# A REVOIR - C'est fait :P
		
		if pointsAlreadyCompute == nil(Solution) || pointsAlreadyCompute.head != leftSol
			lowerBound = cons(leftSol, lowerBound)
		end
		tmpLowerBound = lowerBound
		while tmpLowerBound.tail != nil(Solution)
			push!(toStudy, PairOfSolution(tmpLowerBound.head, tmpLowerBound.tail.head))
			tmpLowerBound = tmpLowerBound.tail
		end
		if tmpLowerBound.head != rightSol
			tmpLowerBound.tail = cons(rightSol, nil(Solution))
			push!(toStudy, PairOfSolution(tmpLowerBound.head, rightSol))
		end
		#verbose && println("toStudy origin: $toStudy")

		#while some pairs of solutions are found by the dichotomic search, we keep going
        while toStudy != []
			currPair = pop!(toStudy)
			#verbose && println("toStudy : $toStudy")
			leftSol = currPair.solL
			rightSol = currPair.solR

	        λ = [leftSol.y[2] - rightSol.y[2], rightSol.y[1] - leftSol.y[1]]
			midSol = solve1OKPLinear(prob, λ, assignment)
			# if the solution dominates one of the other, it's added to the LB
			if round(sum(λ .* midSol.y), digits = 7) - round(sum(λ .* rightSol.y), digits = 7) > 0.00001
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
	return sol.y[1] >= (pair.solL.y[1] + 1) && sol.y[2] >= (pair.solR.y[2]+1)
end

function updateLowerBound!(lowerBound::T, nadirPoints::Vector{PairOfSolution}, listOfPoint::LinkedList{Solution}; compteur = nothing) where T<:LinkedList{Solution}

	
	for sol in listOfPoint
		
		
		if sol.isBinary
		
			studiedLowerBound = lowerBound
			
			testFirstIn = false
			testFirstOut = false
			
			eltBeforeIn = nil(Solution)
			eltAfterOut = nil(Solution)
			
			
			while studiedLowerBound != nil(Solution) && !testFirstOut && !(!testFirstIn && studiedLowerBound.head.y[1] > sol.y[1])
				
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
					end
					if nadirPoints[iterNadir].solR.y == eltAfterOut.head.y
						indAfterOut = iterNadir + 1
					end
					iterNadir += 1
				end
				nadirPoints = append!(push!(nadirPoints[1:indBeforeIn], PairOfSolution(nadirPoints[indBeforeIn+1].solL, sol), PairOfSolution(sol, nadirPoints[indAfterOut-1].solR)), nadirPoints[indAfterOut:end])
				
			elseif testFirstIn || testFirstOut
				error("testFirstIn : $testFirstIn - testFirstOut : $testFirstOut\n lowerBound : $(map(sol->sol.y, lowerBound))")
			else
				iterPair = 1
				testDomiNadir = false
				pairDomin = PairOfSolution()
				
				while iterPair <= length(nadirPoints) && !testDomiNadir
					pairNadir = nadirPoints[iterPair]
					
					if dominate(sol, pairNadir)
						testDomiNadir = true
						nadirPoints = append!(push!(nadirPoints[1:(iterPair-1)], PairOfSolution(pairNadir.solL, sol), PairOfSolution(sol, pairNadir.solR)), nadirPoints[(iterPair+1):end])
						pairDomin = pairNadir
						
					end
					iterPair += 1
				end
				
				if testDomiNadir
					studiedLowerBound = lowerBound
					testSolFound = false
					
					while studiedLowerBound != nil(Solution) && !testSolFound
						
						if studiedLowerBound.head == pairDomin.solL
							studiedLowerBound.tail = cons(sol, studiedLowerBound.tail)
							
							
							if studiedLowerBound.tail == nil(Solution) || studiedLowerBound.tail.tail.head != pairDomin.solR
								
								@assert false "J'en peux plus des bugs ! $(studiedLowerBound.tail.tail.head.y) - $(pairDomin.solR.y)"
							end
							
							testSolFound = true
						end
						
						studiedLowerBound = studiedLowerBound.tail
					end
				end
				
			end
		end
	end
	
	return nadirPoints

end

function pruningTest(lowerBound::T, listPointsNadir::Vector{PairOfSolution}, subUpperBound::DualSet, verbose = false; compteur = nothing) where T<:LinkedList{Solution}

	lowerBound == nil(Solution) && return infeasibility, Vector{PairOfSolution}(), listPointsNadir
	lowerBound.tail == nil(Solution) && lowerBound.head.isBinary && return optimality, listPointsNadir, Vector{PairOfSolution}()
	
	lowerBound.tail == nil(Solution) && return none, listPointsNadir, Vector{PairOfSolution}()
	
	subIsDominated = true
	MyDude = Vector{PairOfSolution}()
	dominatedPoints = Vector{PairOfSolution}()

	for pairNadir in listPointsNadir

		verbose && println("pairNadir : $pairNadir")
		nadir = [pairNadir.solL.y[1], pairNadir.solR.y[2]]
		nadirA = subUpperBound.A * nadir

		iter = 1
		
		nadirInUpper = true
	
		while iter <= length(subUpperBound.b) && nadirInUpper
			
			nadirInUpper = nadirInUpper && nadirA[iter] <= subUpperBound.b[iter]
			iter += 1
		end
		if nadirInUpper
			push!(MyDude, pairNadir)
		else
			push!(dominatedPoints, pairNadir)
		end
		
		subIsDominated = subIsDominated && !nadirInUpper
	end

	if subIsDominated
		return dominance, Vector{PairOfSolution}(), listPointsNadir
	else
		return none, MyDude, dominatedPoints
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

function computeAlreadyPoints(subLowerBound::T, index::Int, value::Int) where T <: LinkedList{Solution}
	
	newList = nil(Solution)
	solInNewList = nil(Solution)
	testFirstVarIn = false
	
	while subLowerBound != nil(Solution)
		if !testFirstVarIn && subLowerBound.head.x[index] == value
			testFirstVarIn = true
			newList = cons(subLowerBound.head, nil(Solution))
			solInNewList = newList
		elseif subLowerBound.head.x[index] == value
			solInNewList.tail = cons(subLowerBound.head, nil(Solution))
			solInNewList = solInNewList.tail
		end
		
		subLowerBound = subLowerBound.tail
	end
	
	return newList

end

function branchAndBound!(lowerBound::T, prob::Problem, assignment::Assignment, nadirPointsToStudy::Vector{PairOfSolution}; pointsAlreadyCompute = nil(Solution), M::Float64 = 1000., compteur = nothing, timeMax = nothing, start = nothing) where T<:LinkedList{Solution}
	
	if timeMax != nothing && time() - start > timeMax
		return nadirPointsToStudy, false
	end
	testTime = false
	
	compteur != nothing && (compteur.value += 1)
	compteur.value%10000 == 0 && println(compteur.value)	
	
	subLowerBound = dichoSearch(prob, assignment, pointsAlreadyCompute = pointsAlreadyCompute, M = M, compteur = compteur) #Nathalie	
	subUpperBound = computeUpperBound(subLowerBound) #Lucas
	

	prunedType, newNadirPoints, dominatedNadir = pruningTest(subLowerBound, nadirPointsToStudy, subUpperBound, compteur = compteur)


	if prunedType == none
		newNadirPoints = updateLowerBound!(lowerBound, newNadirPoints, subLowerBound, compteur = compteur) #Lucas/Jules/Nathalie
	elseif prunedType == optimality
		newNadirPoints = updateLowerBound!(lowerBound, newNadirPoints, subLowerBound, compteur = compteur)
	end
	

	if prunedType == none && assignment.assignEndIndex < prob.nbVar
		testAddVar = (assignment.weight + prob.weights[assignment.assignEndIndex + 1] <= prob.maxWeight)


		if testAddVar
			addVarAssignment!(assignment, prob) #Lucas
			pointsAlreadyCompute = computeAlreadyPoints(subLowerBound, assignment.assignEndIndex, 1)
			newNadirPoints, testTime = branchAndBound!(lowerBound, 
												prob, 
												assignment, 
												newNadirPoints, 
												pointsAlreadyCompute = pointsAlreadyCompute, 
												M = M,
												compteur = compteur, 
												timeMax = timeMax,
												start = start) #Lucas
		end
		removeVarAssignment!(assignment, prob, testAddVar) #Lucas
		pointsAlreadyCompute = computeAlreadyPoints(subLowerBound, assignment.assignEndIndex, 0)
		newNadirPoints, testTime = branchAndBound!(lowerBound, 
											prob, 
											assignment, 
											newNadirPoints,  
											pointsAlreadyCompute = pointsAlreadyCompute, 
											M = M,
											compteur = compteur, 
											timeMax = timeMax,
											start = start) #Lucas

		returnParentAssignment!(assignment, prob) #Lucas
		return merge(newNadirPoints, dominatedNadir), (timeMax == nothing || time() - start < timeMax)
		
	elseif prunedType == optimality
		return merge(newNadirPoints, dominatedNadir), (timeMax == nothing || time() - start < timeMax)
	else
		return nadirPointsToStudy, (timeMax == nothing || time() - start < timeMax)
	end

end

function main(prob::Problem; M::Float64 = 1000., timeMax = nothing, start = nothing, withHeuristic::Bool = false)

	@assert prob.nbObj == 2 "This Branch and Bound supports only a Bio-objective problem"

	assignment = Assignment(prob) #Nathalie

	compt = Compteur()
	
	lowerBound = dichoSearch(prob, assignment, M = M, compteur = compt) #Nathalie
	
	if withHeuristic
		
		n = length(lowerBound)
		lowerBoundSet = Vector{Solution}(undef, n)
		consecutiveSet = Vector{PairOfSolution}(undef, n-1)
		
		lowerBoundBis = lowerBound
		iter = 1
		while lowerBoundBis.tail != nil(Solution)
			lowerBoundSet[iter] = lowerBoundBis.head
			consecutiveSet[iter] = PairOfSolution(lowerBoundBis.head, lowerBoundBis.tail.head)
			
			lowerBoundBis = lowerBoundBis.tail
			iter += 1
		end
		lowerBoundSet[iter] = lowerBoundBis.head
	
		lowerBound, tmp2 = improveLowerBoundSet(lowerBoundSet, consecutiveSet, prob)
		
	end

	nadirPoints = getNadirPoints(lowerBound) #Jules
	
	list, testTime = branchAndBound!(lowerBound, prob, assignment, nadirPoints, pointsAlreadyCompute = nil(Solution), M = M, compteur = compt, timeMax = timeMax, start = start) #Lucass

	println("N° Assignement : $(compt.value)")

	return lowerBound, testTime, compt

end

function experimentation(timeMax::Float64 = 600.)
	
	Random.seed!(123456)
	
	problems = Vector{Vector{Problem}}(undef, 11)
	for nbVar = 50:10:150
		instances = Vector{Problem}(undef, 5)
		for iter = 1:5
			instances[iter] = Problem(nbVar)
		end
		problems[Int(nbVar/10 - 4)] = copy(instances)
	end
	
	println(problems[4][2])
	return 0
	
	println("Début ConvexeBase")
	for instance in problems[4:end]
		
		main(instance[1], timeMax = 5., start = time())

		println("$(instance[1].nbVar) : ")		

		start = time()
		iterInstance = 1
		file = open("../Experimentation/LinearBaseBis", "a")
		while time() - start <= timeMax && iterInstance <= 5
			
			println("Instance n° $iterInstance - $(time()-start)")

			prob = instance[iterInstance]
		
			startInstance = time()
			lb, testTime, compt = main(prob, timeMax = timeMax, start = start)
			endInstance = time()
			if time() - start < timeMax
				write(file, "$(prob.nbVar);$iterInstance;$(endInstance-startInstance);$(compt.value)\n")
			end
			iterInstance += 1
		end
		close(file)
	end
	
end
