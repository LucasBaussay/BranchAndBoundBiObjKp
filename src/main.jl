
using DataStructures, Random, JuMP, GLPK

const comboPath = joinpath(@__DIR__,"..","deps","libcombo.so")

include("dataStruct.jl")
include("parserCombo.jl")
include("steepestDescent.jl")

function invertPerm(perm::Vector{Int})
	n = 0
	for item in perm
		if item > n
			n = item
		end
	end
	
	rev = Vector{Int}(undef, n)
	for iter = 1:length(perm)
		rev[perm[iter]] = iter
	end
	
	return rev
end

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
	
	return Solution(x, profit, weight, isBinary), true
		
	
end

function solve1OKP(prob::Problem, λ::Vector{Float64}, assignment::Assignment)
	#Passer Lambda en paramètre

	model = Model(GLPK.Optimizer)
	x = @variable(model, x[1:(prob.nbVar-assignment.assignEndIndex)], Bin)
	@constraint(model, Weights, sum(x .* prob.weights[(assignment.assignEndIndex+1):end]) + assignment.weight <= prob.maxWeight)
	@objective(model, Max, λ[1] * sum(x .* prob.profits[1, (assignment.assignEndIndex+1):end]) + λ[2] * sum(x .* prob.profits[2, (assignment.assignEndIndex+1):end]))

	optimize!(model)

	X = append!(assignment.assign[1:assignment.assignEndIndex],(Float64).(value.(x)))
	
	termStatus = termination_status(model)

	if termStatus == MOI.OPTIMAL
		return Solution(X, [sum(X .* prob.profits[1,1:end]), sum(X .* prob.profits[2,1:end])], sum(X .* prob.weights), true), true
	elseif termStatus == MOI.INFEASIBLE
		return Solution(), false
	else
		error("Cette Fois c'est MOI qui renvoie de la merde")
	end
end

function dichoSearch(prob::Problem, assignment::Assignment; withFirstConvex::Bool = false, pointsAlreadyCompute = nil(Solution), M::Float64 = 1000., withLinear::Bool = false, verbose::Bool = false, compteur = nothing, linearAmelioration::Bool = false, solutionList::Union{Nothing, Vector{AmeliorateLinear}} = nothing) where T <: LinkedList{Solution}

    lowerBound = copy(pointsAlreadyCompute)
	toStudy = Vector{PairOfSolution}()

    λ = [1, M]
    
    verbose && println("\n $(broadcast(sol->sol.solution.y, solutionList))")
    verbose && println(" $(broadcast(sol->sol.lambdaValue, solutionList))\n")
    
    if withFirstConvex
    	
    	leftSol, testLeft = solve1OKP(prob, λ, assignment)
    	
    elseif linearAmelioration
    	
    	lambdaValue = λ[1] / sum(λ)
    	
    	iter = 1
    	while lambdaValue < solutionList[iter].lambdaValue && iter < length(solutionList)
    		iter += 1
    	end
    	
    	leftSol, testLeft = copy(solutionList[iter].solution), true
	else
		
		if withLinear
			leftSol, testLeft = solve1OKPLinear(prob, λ, assignment)
		else
			leftSol, testLeft = solve_monoBinary(prob, λ, assignment)
		end
	
	end
	
    λ = [M, 1]
    
     if withFirstConvex
    	
    	rightSol, testRight = solve1OKP(prob, λ, assignment)
    	
    elseif linearAmelioration
    	
    	lambdaValue = λ[1] / sum(λ)
    	
    	iter = 1
    	while iter <= length(solutionList) && lambdaValue < solutionList[iter].lambdaValue
    		iter += 1
    	end
    	
    	rightSol, testRight = copy(solutionList[iter].solution), true
	else
		
		if withLinear
			rightSol, testRight = solve1OKPLinear(prob, λ, assignment)
		else
			rightSol, testRight = solve_monoBinary(prob, λ, assignment)
		end
	
	end
	
	@assert testLeft && testRight "Bah écoute, on choppe des solutions infaisables, de mieux en mieux"

	verbose && println("Two found solutions : $(leftSol.y) && $(rightSol.y)")

    # In the event that the two solutions found are identical,
    # this solution is optimal and, thus, added to the lower bound
    if round.(leftSol.y, digits = 7) - round.(rightSol.y, digits = 7) < 0.00001
        lowerBound = cons(leftSol, lowerBound)
		verbose && println("The two solutions are identical")
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

		verbose && println("pointsAlreadyCompute = $(map(sol->sol.y, pointsAlreadyCompute))") #Aled

		#verbose && println("toStudy origin: $toStudy")

		#while some pairs of solutions are found by the dichotomic search, we keep going
        while toStudy != []
			verbose && println("LB = $(map(sol->sol.y, lowerBound))")
			currPair = pop!(toStudy)
			#verbose && println("toStudy : $toStudy")
			leftSol = currPair.solL
			rightSol = currPair.solR

			verbose && println("Found solutions : $(leftSol.y) && $(rightSol.y)")

	        λ = [leftSol.y[2] - rightSol.y[2], rightSol.y[1] - leftSol.y[1]]
	        if withFirstConvex
				
				midSol, testMid = solve1OKP(prob, λ, assignment)
				
			elseif linearAmelioration
				
				lambdaValue = λ[1] / sum(λ)
				
				iter = 1
				while iter <= length(solutionList) && lambdaValue < solutionList[iter].lambdaValue
					iter += 1
				end
				
				midSol, testMid = copy(solutionList[iter].solution), true
				testMid = true
				
				
			else
				
				if withLinear
					midSol, testMid = solve1OKPLinear(prob, λ, assignment)
				else
					midSol, testMid = solve_monoBinary(prob, λ, assignment)
				end
			
			end
			
			verbose && println("Found solution : $(midSol.y)")
			
			@assert testMid "Bah écoute, on choppe des solutions infaisables, de mieux en mieux"

			# if the solution dominates one of the other, it's added to the LB
			if round(sum(λ .* midSol.y), digits = 7) - round(sum(λ .* rightSol.y), digits = 7) > 0.00001
				verbose && println("It ameliorates")
				verbose && println("LB = $(map(sol->sol.y, lowerBound))")
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

function updateLowerBound!(lowerBound::T, nadirPoints::Vector{PairOfSolution}, listOfPoint::LinkedList{Solution}; compteur = nothing, withLinear = false) where T<:LinkedList{Solution}

	
	for sol in listOfPoint
		
		
		if !withLinear || sol.isBinary
		
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



function pruningTest(lowerBound::T, listPointsNadir::Vector{PairOfSolution}, subUpperBound::DualSet, verbose = false; compteur = nothing, avecLesFigures = false) where T<:LinkedList{Solution}

	lowerBound == nil(Solution) && return infeasibility, Vector{PairOfSolution}(), listPointsNadir
	lowerBound.tail == nil(Solution) && lowerBound.head.isBinary && return optimality, listPointsNadir, Vector{PairOfSolution}()
	
	lowerBound.tail == nil(Solution) && return none, listPointsNadir, Vector{PairOfSolution}()
	
	if compteur != nothing && avecLesFigures
		
		fig = figure()
		ax = fig.add_subplot(111)
		
		lengthLower = length(lowerBound)
		solX = Vector{Float64}(undef, lengthLower)
		solY = Vector{Float64}(undef, lengthLower)
		
		iter = 1
		for sol in lowerBound
			solX[iter] = sol.y[1]
			solY[iter] = sol.y[2]
			
			iter += 1
		end
		
		ax.plot(solX, solY, "b-.")
		ax.plot(broadcast(pair->pair.solL.y[1], listPointsNadir), broadcast(pair->pair.solR.y[2], listPointsNadir), "r.")
		ax.set_xlabel("z_1(x)")
		ax.set_ylabel("z_2(x)")
		ax.grid(true)
		fig.savefig("SaveFig/lowerBound_$(compteur.value).png")
		close(fig)
	end
	
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


function linearPretreatment(permObjBis::Vector{Int}, permLambda::Vector{Int}, lambdaList::Vector{LambdaChange}, indEndLambdaList::Int, prob::Problem, assignment::Assignment)

	permObj = permObjBis[1:end]
	revPermObj = invertPerm(permObj)
	
	solutionList = Vector{AmeliorateLinear}(undef, indEndLambdaList+1)
	indEndSolutionList = 1
	
	nbVar = prob.nbVar - assignment.assignEndIndex
	
	
	indSol = 0
	weight = assignment.weight
	profit = assignment.profit[1:end]
	#=x = Vector{Float64}(undef, prob.nbVar)
	for iter = 1:prob.nbVar
		if assignment.assign[iter] == -1
			x[iter] = 0.
		else
			x[iter] = assignment.assign[iter]
		end
	end=#
	x = Float64.(assignment.assign[1:end])
	isBinary = true
	while indSol < nbVar && weight < prob.maxWeight
		indSol += 1
		if weight + prob.weights[permObj[indSol]] <= prob.maxWeight
			x[permObj[indSol]] = 1.
			profit += prob.profits[1:end, permObj[indSol]]
			weight += prob.weights[permObj[indSol]]
		else
			x[permObj[indSol]] = (prob.maxWeight - weight)/prob.weights[permObj[indSol]]
			weight = prob.maxWeight
			profit += x[permObj[indSol]] * prob.profits[1:end, permObj[indSol]]
			isBinary = false
		end
		
	end
	
	if (x[permObj[indSol]] == 1.)
		indSol += 1
	end
	
	indSol == 0 ? indSol = 1 : nothing
	indSol == nbVar+1 ? indSol = nbVar : nothing
	
	sol = Solution(copy(x), copy(profit), weight, isBinary)
	
	solutionList[indEndSolutionList] = AmeliorateLinear(sol, indEndLambdaList >= 1 ? lambdaList[permLambda[1]].value : 0.)
	
	
	for indLambda = 1:indEndLambdaList
		nextVal = indEndLambdaList >= (indLambda+1) ? lambdaList[permLambda[indLambda+1]].value : 0.
		
		index1 = revPermObj[lambdaList[permLambda[indLambda]].index[1]] #Indice dans la liste triée
		index2 = revPermObj[lambdaList[permLambda[indLambda]].index[2]]
		
		if index2 < index1
			index1, index2 = index2, index1
		end
		
		permObj[index1], permObj[index2] = permObj[index2], permObj[index1]
		revPermObj[permObj[index1]], revPermObj[permObj[index2]] = revPermObj[permObj[index2]], revPermObj[permObj[index1]]
		
		#Liste triée  -> Liste réelle : perm
		#Liste réelle -> Liste triée  : revPerm
		
		if index1 < indSol
		
			if index2 < indSol
				nothing
			elseif index2 == indSol #IL EST VALIDE ++
				floatPart = x[permObj[index1]]
				
				x[permObj[index1]] = 1.
				profit += (1. - floatPart) * prob.profits[1:end, permObj[index1]]
				weight += (1. - floatPart) * prob.weights[permObj[index1]]
				
				while weight > prob.maxWeight && indSol > 0
					x[permObj[indSol]] = max(0., (prob.maxWeight - (weight - prob.weights[permObj[indSol]]))/prob.weights[permObj[indSol]])
					
					profit -= (1 - x[permObj[indSol]]) * prob.profits[1:end, permObj[indSol]]
					weight -= (1 - x[permObj[indSol]]) * prob.weights[permObj[indSol]]
					
					if (x[permObj[indSol]] == 0.)
						indSol -= 1
					end
					
					indSol == 0 ? indSol = 1 : nothing
				end
				
				
				
				if lambdaList[indLambda].value != nextVal
					indEndSolutionList += 1
					solutionList[indEndSolutionList] = AmeliorateLinear(Solution(copy(x), copy(profit), weight, isBinary), nextVal)
				end
			else #Valide ++
				
				x[permObj[index1]], x[permObj[index2]] = 1., 0.
				
				profit += prob.profits[1:end, permObj[index1]] - prob.profits[1:end, permObj[index2]]
				weight += prob.weights[permObj[index1]] - prob.weights[permObj[index2]]
				
				if weight > prob.maxWeight
				
					profit -= x[permObj[indSol]] * prob.profits[1:end, permObj[indSol]]
					weight -= x[permObj[indSol]] * prob.weights[permObj[indSol]]
					x[permObj[indSol]] = 1.
					profit += prob.profits[1:end, permObj[indSol]]
					weight += prob.weights[permObj[indSol]]
					
					while indSol <= nbVar && weight > prob.maxWeight && indSol > 0
						
						x[permObj[indSol]] = max(0., (prob.maxWeight - (weight - prob.weights[permObj[indSol]]))/prob.weights[permObj[indSol]])
						profit -= (1. - x[permObj[indSol]]) * prob.profits[1:end, permObj[indSol]]
						weight -= (1. - x[permObj[indSol]]) * prob.weights[permObj[indSol]]
						
						if (x[permObj[indSol]] == 0.)
							indSol -= 1
						end
					end
					
					indSol == 0 ? indSol = 1 : nothing
					
					isBinary = (x[permObj[indSol]] == 1. || x[permObj[indSol]] == 0.)
					
					if lambdaList[indLambda].value != nextVal
						indEndSolutionList += 1
						solutionList[indEndSolutionList] = AmeliorateLinear(Solution(copy(x), copy(profit), weight, isBinary), nextVal)
					end
					
				elseif weight < prob.maxWeight
					
					profit -= x[permObj[indSol]] * prob.profits[1:end, permObj[indSol]]
					weight -= x[permObj[indSol]] * prob.weights[permObj[indSol]]
					x[permObj[indSol]] = 0.
					
					while indSol <= nbVar && weight < prob.maxWeight && indSol > 0
						x[permObj[indSol]] = min(1., (prob.maxWeight - (weight - prob.weights[permObj[indSol]]))/prob.weights[permObj[indSol]])
						profit += x[permObj[indSol]] * prob.profits[1:end, permObj[indSol]]
						weight += x[permObj[indSol]] * prob.weights[permObj[indSol]]
						
						if (x[permObj[indSol]] == 1.)
							indSol += 1
						end
					end
					
					indSol == nbVar+1 ? indSol = nbVar : nothing
					
					isBinary = (x[permObj[indSol]] == 1. || x[permObj[indSol]] == 0.)
					
					if lambdaList[indLambda].value != nextVal
						indEndSolutionList += 1
						solutionList[indEndSolutionList] = AmeliorateLinear(Solution(copy(x), copy(profit), weight, isBinary), nextVal)
					end
					
				else
					nothing
				end
				
			end
		
		elseif index1 == indSol
		
			floatPart = x[permObj[index2]]
		
			x[permObj[index1]], x[permObj[index2]] = 1., 0.
			profit += prob.profits[1:end, permObj[index1]] - floatPart * prob.profits[1:end, permObj[index2]]
			weight += prob.weights[permObj[index1]] - floatPart * prob.weights[permObj[index2]]
			####################""
			
			if weight > prob.maxWeight
				
				while weight > prob.maxWeight && indSol > 0
					x[permObj[indSol]] = max(0., (prob.maxWeight - (weight - prob.weights[permObj[indSol]]))/prob.weights[permObj[indSol]])
					profit -= (1. - x[permObj[indSol]]) * prob.profits[1:end, permObj[indSol]]
					weight -= (1. - x[permObj[indSol]]) * prob.weights[permObj[indSol]]
					
					if (x[permObj[indSol]] == 0.)
						indSol -= 1
					end
				end
				
				indSol == 0 ? indSol = 1 : nothing
				
				isBinary = (x[permObj[indSol]] == 1. || x[permObj[indSol]] == 0.)
				
				if lambdaList[indLambda].value != nextVal
					indEndSolutionList += 1
					solutionList[indEndSolutionList] = AmeliorateLinear(Solution(copy(x), copy(profit), weight, isBinary), nextVal)
				end
				
			elseif weight < prob.maxWeight
				
				indSol += 1
				
				while indSol <= nbVar && weight < prob.maxWeight
					x[permObj[indSol]] = min(1., (prob.maxWeight - (weight - prob.weights[permObj[indSol]]))/prob.weights[permObj[indSol]])
					profit += x[permObj[indSol]] * prob.profits[1:end, permObj[indSol]]
					weight += x[permObj[indSol]] * prob.weights[permObj[indSol]]
					
					if (x[permObj[indSol]] == 1.)
						indSol += 1
					end
				end
				
				indSol == nbVar+1 ? indSol = nbVar : nothing
				
				isBinary = (x[permObj[indSol]] == 1. || x[permObj[indSol]] == 0.)
				
				if lambdaList[indLambda].value != nextVal
					indEndSolutionList += 1
					solutionList[indEndSolutionList] = AmeliorateLinear(Solution(copy(x), copy(profit), weight, isBinary), nextVal)
				end
				
			else
				nothing
			end
		
		else
			nothing
		end
		
	end
		
	return solutionList[1:indEndSolutionList]
	
end


function branchAndBound!(lowerBound::T, prob::Problem, assignment::Assignment, nadirPointsToStudy::Vector{PairOfSolution}, permObj, permLambda, lambdaList, indEndLambdaList::Int; withFirstConvex::Bool = false, pointsAlreadyCompute = nil(Solution), M::Float64 = 1000., num::Int = 1, compteur = nothing, avecLesFigures = false, withLinear::Bool = false, linearAmelioration::Bool = false, timeMax = nothing, start = nothing) where T<:LinkedList{Solution}
	
	if timeMax != nothing && time() - start > timeMax
		return nadirPointsToStudy, false
	end
	testTime = false
	
	compteur != nothing && (compteur.value += 1)
	compteur.value%10000 == 0 && println(compteur.value)
	
	solutionList = nothing
	if linearAmelioration 
	
		if assignment.assignEndIndex != prob.nbVar
			
			solutionList = linearPretreatment(permObj, permLambda, lambdaList, indEndLambdaList, prob, assignment)
		else
			solutionList = [AmeliorateLinear(Solution(assignment.assign[1:end], assignment.profit[1:end], assignment.weight, true), 0.)]
		end
		
		for iter = 1:length(solutionList)
			println("$(solutionList[iter].solution.y) - $(solutionList[iter].solution.x) - $(solutionList[iter].lambdaValue)")
		end
		println()
		
	end
	
	
	
	subLowerBound = dichoSearch(prob, assignment, withFirstConvex = withFirstConvex, pointsAlreadyCompute = pointsAlreadyCompute, withLinear = withLinear, M = M, compteur = compteur, linearAmelioration = linearAmelioration, solutionList = solutionList) #Nathalie	
	subUpperBound = computeUpperBound(subLowerBound) #Lucas

	prunedType, newNadirPoints, dominatedNadir = pruningTest(subLowerBound, nadirPointsToStudy, subUpperBound, compteur = compteur, avecLesFigures = avecLesFigures) #Nathalie

	#println("lowerBound : $(map(sol->sol.y, lowerBound))")
	#println("Pruned : $prunedType")
	#println()
	#println("subLowerBound : $(map(sol->sol.y, subLowerBound))")
	#println("nadirPoints : $(broadcast(pair->(pair.solL.y[1], pair.solR.y[2]), nadirPointsToStudy)), newNadirPoints : $(broadcast(pair->(pair.solL.y[1], pair.solR.y[2]), newNadirPoints))")
	#println("nadirPoints : $(broadcast(pair->(pair.solL.y, pair.solR.y), nadirPointsToStudy)), newNadirPoints : $(broadcast(pair->(pair.solL.y, pair.solR.y), newNadirPoints))")


	if prunedType == none
		newNadirPoints = updateLowerBound!(lowerBound, newNadirPoints, subLowerBound, withLinear = withLinear, compteur = compteur) #Lucas/Jules/Nathalie
	elseif prunedType == optimality
		newNadirPoints = updateLowerBound!(lowerBound, newNadirPoints, subLowerBound, withLinear = withLinear, compteur = compteur)
	end

	if prunedType == none && assignment.assignEndIndex < prob.nbVar
		testAddVar = (assignment.weight + prob.weights[assignment.assignEndIndex + 1] <= prob.maxWeight)

		newPermObj = nothing
		newPermLambda = nothing
		newIndEndLambdaList = 0

		if linearAmelioration
			newPermObj = filter(x->x!=assignment.assignEndIndex+1, permObj)
			
			newPermLambda = filter( x-> lambdaList[x].index[1] != assignment.assignEndIndex+1 && lambdaList[x].index[2] != assignment.assignEndIndex+1, permLambda)
			newIndEndLambdaList = length(newPermLambda)
		end

		if testAddVar
			addVarAssignment!(assignment, prob) #Lucas
			pointsAlreadyCompute = computeAlreadyPoints(subLowerBound, assignment.assignEndIndex, 1)
			newNadirPoints, testTime = branchAndBound!(lowerBound, 
												prob, 
												assignment, 
												newNadirPoints, 
												newPermObj, 
												newPermLambda,
												lambdaList,
												newIndEndLambdaList, 
												withFirstConvex = withFirstConvex,
												pointsAlreadyCompute = pointsAlreadyCompute, 
												M = M, 
												num = num + 1, 
												withLinear = withLinear, 
												compteur = compteur, 
												avecLesFigures = avecLesFigures, 
												linearAmelioration = linearAmelioration,
												timeMax = timeMax,
												start = start) #Lucas
		end
		removeVarAssignment!(assignment, prob, testAddVar) #Lucas
		pointsAlreadyCompute = computeAlreadyPoints(subLowerBound, assignment.assignEndIndex, 0)
		newNadirPoints, testTime = branchAndBound!(lowerBound, 
											prob, 
											assignment, 
											newNadirPoints,  
											newPermObj, 
											newPermLambda,
											lambdaList,
											newIndEndLambdaList,
											withFirstConvex = withFirstConvex,
											pointsAlreadyCompute = pointsAlreadyCompute, 
											M = M,
											num = num + 1, 
											compteur = compteur, 
											withLinear = withLinear, 
											avecLesFigures = avecLesFigures, 
											linearAmelioration = linearAmelioration,
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


function initializeLinearPretreatment(prob::Problem)


	ratio = zeros(Float64, 2, prob.nbVar)
	lambdaList = Vector{LambdaChange}(undef, (prob.nbVar^2 - prob.nbVar))
	indEndLambdaList = 0
	
	for iter = 1:prob.nbVar
		ratio[1, iter] = prob.profits[1, iter] / prob.weights[iter]
		ratio[2, iter] = prob.profits[2, iter] / prob.weights[iter]
	end
	
	for iter = 1:prob.nbVar
		for iter2 = 1:prob.nbVar
			if iter2 > iter
				value = (ratio[1, iter] - ratio[2, iter] - ratio[1, iter2] + ratio[2, iter2])
				if value != 0.
					val = (ratio[2, iter2] - ratio[2, iter]) / value
					if val < 1 && val > 0
						indEndLambdaList += 1
						lambdaList[indEndLambdaList] = LambdaChange((iter, iter2), val) #Indice réel
					end
				end
			end
		end
	end

	permLambda = sortperm(lambdaList[1:indEndLambdaList], by = x->x.value, rev = true)
	
	permObj = sortperm(1:prob.nbVar, lt = (x, y) -> if ratio[1, x] == ratio[1, y]
														isless(ratio[2, x], ratio[2, y])
													else
														isless(ratio[1, x], ratio[1, y])
													end,
													rev = true)
													
	return permObj, permLambda, lambdaList, indEndLambdaList
													
end

function main(prob::Problem; withFirstConvex::Bool = false, withLinear::Bool = false, M::Float64 = 1000., avecLesFigures::Bool = false, linearAmelioration::Bool = false, timeMax = nothing, start = nothing, withHeuristic::Bool = true)

	@assert prob.nbObj == 2 "This Branch and Bound supports only a Bio-objective problem"

	assignment = Assignment(prob) #Nathalie

	compt = Compteur()
	
	solutionList = nothing
	permObj = nothing
	permLambda = nothing
	lambdaList = nothing
	indEndLambdaList = 0
	if linearAmelioration
		permObj, permLambda, lambdaList, indEndLambdaList = initializeLinearPretreatment(prob)
		solutionList = linearPretreatment(permObj, permLambda, lambdaList, indEndLambdaList, prob, assignment)
	end
	
	lowerBound = dichoSearch(prob, assignment, M = M, compteur = compt, withFirstConvex = withFirstConvex) #Nathalie
	
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
		
		println(lowerBoundSet)
		println(consecutiveSet)
	
		lowerBound, tmp2 = improveLowerBoundSet(lowerBoundSet, consecutiveSet, prob)
		
	end

	nadirPoints = getNadirPoints(lowerBound) #Jules
	
	list, testTime = branchAndBound!(lowerBound, prob, assignment, nadirPoints, permObj, permLambda, lambdaList, indEndLambdaList, withFirstConvex = withFirstConvex, pointsAlreadyCompute = nil(Solution), M = M, compteur = compt, avecLesFigures = avecLesFigures, withLinear = withLinear, linearAmelioration = linearAmelioration, timeMax = timeMax, start = start) #Lucass

	println("N° Assignement : $(compt.value)")

	return lowerBound, testTime, compt

end

function main(fname::String; withLinear::Bool = false, withFirstConvex::Bool = false, withHeuristic::Bool = true, M::Float64 = 1000., avecLesFigures::Bool = false, linearAmelioration::Bool = false)
	prob = Problem(fname)
	return main(prob, withLinear = withLinear, withFirstConvex = withFirstConvex, withHeuristic = withHeuristic, M = M, avecLesFigures = avecLesFigures, linearAmelioration = linearAmelioration)
end

function main(;withLinear::Bool = false, withFirstConvex::Bool = false, withHeuristic::Bool = true, M::Float64 = 1000., avecLesFigures::Bool = false, linearAmelioration::Bool = false)
	prob = Problem()
	return main(prob, withLinear = withLinear, withFirstConvex = withFirstConvex, withHeuristic = withHeuristic, M = M, avecLesFigures = avecLesFigures, linearAmelioration = linearAmelioration)
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
	println(problems[3][2])
	
	println("Début ConvexeBase")
	for instance in problems
		
		main(instance[1], withLinear = true, timeMax = 5., start = time())

		println("$(instance[1].nbVar) : ")		

		start = time()
		iterInstance = 1
		file = open("../Experimentation/LinearBase", "a")
		while time() - start <= timeMax && iterInstance <= 5
			
			println("Instance n° $iterInstance - $(time()-start)")

			prob = instance[iterInstance]
		
			startInstance = time()
			lb, testTime, compt = main(prob, withLinear = true, timeMax = timeMax, start = start)
			endInstance = time()
			if time() - start < timeMax
				write(file, "$(prob.nbVar);$iterInstance;$(endInstance-startInstance);$(compt.value)\n")
			end
			iterInstance += 1
		end
		close(file)
	end
	
#=
	println("Début : LinearBase")
			
	file = open("../Experimentation/LinearBase", "a")
	for instance in problems
		#instance = problems[3]
		println("$(instance[1].nbVar) : ")
		main(instance[1], withLinear = true, timeMax = 5., start = time())
		
		start = time()
		iterInstance = 1
	
		while time() - start <= timeMax && iterInstance <= 5
			println("Instance n°$iterInstance - $(time()-start)")
			prob = instance[iterInstance]
		
			startInstance = time()
			lb, testTime, compt = main(prob, withLinear = true, timeMax = timeMax, start = start)
			endInstance = time()
			if time() - start < timeMax
				write(file, "$(prob.nbVar);$iterInstance;$(endInstance-startInstance);$(compt.value)\n")
			end
			iterInstance += 1
		end
	end
	close(file)
	
	file = open("../Experimentation/ConvexBetter", "a")
	for instance in problems
	
		main(instance[1], timeMax = 5., start = time())
	
		start = time()
		iterInstance = 1
	
		while time() - start <= timeMax && iterInstance <= 5
		
			prob = instance[iterInstance]
		
			startInstance = time()
			lb, testTime, compt = main(prob, timeMax = timeMax, start = start)
			endInstance = time()
			if time() - start < timeMax
				write(file, "$(prob.nbVar);$iterInstance;$(endInstance-startInstance);$(compt.value)\n")
			end
			iterInstance += 1
		end
	end
	close(file)
	
	file = open("../Experimentation/LinearBaseHeuristic", "a")
	for instance in problems
	
		main(instance[1], withLinear = true, withHeuristic = true, timeMax = 5., start = time())
	
		start = time()
		iterInstance = 1
	
		while time() - start <= timeMax && iterInstance <= 5
		
			prob = instance[iterInstance]
		
			startInstance = time()
			lb, testTime, compt = main(prob, withHeuristic = true, withLinear = true, timeMax = timeMax, start = start)
			endInstance = time()
			if time() - start < timeMax
				write(file, "$(prob.nbVar);$iterInstance;$(endInstance-startInstance);$(compt.value)\n")
			end
			iterInstance += 1
		end
	end
	close(file)
	
	file = open("../Experimentation/ConvexBetterHeuristic", "a")
	for instance in problems
	
		main(instance[1], withHeuristic = true, timeMax = 5., start = time())
	
		start = time()
		iterInstance = 1
	
		while time() - start <= timeMax && iterInstance <= 5
		
			prob = instance[iterInstance]
		
			startInstance = time()
			lb, testTime, compt = main(prob, withHeuristic = true, timeMax = timeMax, start = start)
			endInstance = time()
			if time() - start < timeMax
				write(file, "$(prob.nbVar);$iterInstance;$(endInstance-startInstance);$(compt.value)\n")
			end
			iterInstance += 1
		end
	end
	close(file)=#
end
