struct combo_item
    p::Clonglong
    w::Clonglong
    x::Cint
    i::Cint
end

function solve_monoBinary(prob::Problem, λ::Vector{Float64}, assignment::Assignment)
	
	weightedProfits = [sum(λ .* prob.profits[1:end, iter]) for iter=1:prob.nbVar]

	if prob.nbVar - assignment.assignEndIndex == 0
		return Solution(Float64.(assignment.assign), assignment.profit, assignment.weight, true), true
	else
	
		test = true
		for iter = (assignment.assignEndIndex+1):prob.nbVar
			test = test && prob.weights[iter] > (prob.maxWeight - assignment.weight)
		end
		
		if test
			return Solution(append!(Float64.(assignment.assign[1:assignment.assignEndIndex]), zeros(Float64, prob.nbVar-assignment.assignEndIndex)), 
							assignment.profit,
							assignment.weight,
							true), true
		elseif sum(prob.weights[(assignment.assignEndIndex+1):end]) <= (prob.maxWeight - assignment.weight)
			
			return Solution(append!(Float64.(assignment.assign[1:assignment.assignEndIndex]), ones(Float64, prob.nbVar-assignment.assignEndIndex)), 
							assignment.profit + [sum(prob.profits[iter, (assignment.assignEndIndex+1):end]) for iter=1:prob.nbObj],
							assignment.weight + sum(prob.weights[(assignment.assignEndIndex+1):end]),
							true), true
		end
	end
	
	items = [combo_item(weightedProfits[iter], prob.weights[iter], 0, iter) for iter = (assignment.assignEndIndex+1):prob.nbVar]
	
	z = ccall((:solve, comboPath),
    Clonglong, 
    (Ref{combo_item}, Cint, Clonglong, Clonglong, Clonglong),
    items, (prob.nbVar - assignment.assignEndIndex), (prob.maxWeight - assignment.weight), 0, 0)
	

    if z == 0
        z = dot(weightedProfits[(assignment.assignEndIndex+1):end], broadcast(it->it.x, items))
    end
    
    x = Float64.(append!(assignment.assign[1:assignment.assignEndIndex], falses(prob.nbVar-assignment.assignEndIndex)))
    for it in items
    	x[it.i] = Float64(it.x)
    end
    
    return Solution(x,
    				assignment.profit + [sum(x[(assignment.assignEndIndex+1):end] .* prob.profits[iter, (assignment.assignEndIndex+1):end]) for iter = 1:prob.nbObj],
    				assignment.weight + sum(prob.weights[(assignment.assignEndIndex+1):end] .* x[(assignment.assignEndIndex+1):end]),
    				true), true
end
