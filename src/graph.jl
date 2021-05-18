using PyPlot

function draw3Graph()

	file = open("../Experimentation/"*"ConvexeBase", "r")
	
	moyTimeConvexeBase = zeros(Float64, 11)
	moyNbIterConvexeBase = zeros(Float64, 11)
	nbInstConvexeBase = zeros(Int, 11)

	x = readline(file)
	while x != ""
		informations = split(x, ";")
		nbVar = parse(Int, informations[1])
		numInst = parse(Int, informations[2])
		timeTaken = parse(Float64, informations[3])
		nbIter = parse(Int, informations[4])
		
		moyTimeConvexeBase[Int(nbVar/10 - 4)] += timeTaken
		moyNbIterConvexeBase[Int(nbVar/10 - 4)] += nbIter
		nbInstConvexeBase[Int(nbVar/10 - 4)] += 1
	
		x = readline(file)
	end
	moyTimeConvexeBase ./= nbInstConvexeBase
	moyNbIterConvexeBase ./= nbInstConvexeBase
	
	close(file)
	
	file = open("../Experimentation/LinearBase", "r")
	
	moyTimeLinearBase = zeros(Float64, 11)
	moyNbIterLinearBase = zeros(Float64, 11)
	nbInstLinearBase = zeros(Int, 11)

	x = readline(file)
	while x != ""
		informations = split(x, ";")
		nbVar = parse(Int, informations[1])
		numInst = parse(Int, informations[2])
		timeTaken = parse(Float64, informations[3])
		nbIter = parse(Int, informations[4])
		
		moyTimeLinearBase[Int(nbVar/10 - 4)] += timeTaken
		moyNbIterLinearBase[Int(nbVar/10 - 4)] += nbIter
		nbInstLinearBase[Int(nbVar/10 - 4)] += 1
	
		x = readline(file)
	end
	moyTimeLinearBase ./= nbInstLinearBase
	moyNbIterLinearBase ./= nbInstLinearBase
	
	close(file)
	
	file = open("../Experimentation/Combo", "r")
	
	moyTimeCombo = zeros(Float64, 11)
	moyNbIterCombo = zeros(Float64, 11)
	nbInstCombo = zeros(Int, 11)

	x = readline(file)
	while x != ""
		informations = split(x, ";")
		nbVar = parse(Int, informations[1])
		numInst = parse(Int, informations[2])
		timeTaken = parse(Float64, informations[3])
		nbIter = parse(Int, informations[4])
		
		moyTimeCombo[Int(nbVar/10 - 4)] += timeTaken
		moyNbIterCombo[Int(nbVar/10 - 4)] += nbIter
		nbInstCombo[Int(nbVar/10 - 4)] += 1
	
		x = readline(file)
	end
	moyTimeCombo ./= nbInstCombo
	moyNbIterCombo ./= nbInstCombo
	
	close(file)
	
	file = open("../Experimentation/ComboHeuristic", "r")
	
	moyTimeComboHeuristic = zeros(Float64, 11)
	moyNbIterComboHeuristic = zeros(Float64, 11)
	nbInstComboHeuristic = zeros(Int, 11)

	x = readline(file)
	while x != ""
		informations = split(x, ";")
		nbVar = parse(Int, informations[1])
		numInst = parse(Int, informations[2])
		timeTaken = parse(Float64, informations[3])
		nbIter = parse(Int, informations[4])
		
		moyTimeComboHeuristic[Int(nbVar/10 - 4)] += timeTaken
		moyNbIterComboHeuristic[Int(nbVar/10 - 4)] += nbIter
		nbInstComboHeuristic[Int(nbVar/10 - 4)] += 1
	
		x = readline(file)
	end
	moyTimeComboHeuristic ./= nbInstComboHeuristic
	moyNbIterComboHeuristic ./= nbInstComboHeuristic
	
	close(file)
	
	file = open("../Experimentation/LinearHeuristic", "r")
	
	moyTimeLinearHeuristic = zeros(Float64, 11)
	moyNbIterLinearHeuristic = zeros(Float64, 11)
	nbInstLinearHeuristic = zeros(Int, 11)

	x = readline(file)
	while x != ""
		informations = split(x, ";")
		nbVar = parse(Int, informations[1])
		numInst = parse(Int, informations[2])
		timeTaken = parse(Float64, informations[3])
		nbIter = parse(Int, informations[4])
		
		moyTimeLinearHeuristic[Int(nbVar/10 - 4)] += timeTaken
		moyNbIterLinearHeuristic[Int(nbVar/10 - 4)] += nbIter
		nbInstLinearHeuristic[Int(nbVar/10 - 4)] += 1
	
		x = readline(file)
	end
	moyTimeLinearHeuristic ./= nbInstLinearHeuristic
	moyNbIterLinearHeuristic ./= nbInstLinearHeuristic
	
	close(file)
	
	
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeConvexeBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeConvexeBase[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeConvexeBase.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeConvexeBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterConvexeBase[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Nombre de sous-problème moyen étudié")
	ax.grid(true)
	fig.savefig("SaveFig/moyNbIterConvexeBase.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeConvexeBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeConvexeBase[Int.((solX ./ 10) .- 4)] ./ moyNbIterConvexeBase[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen d'un sous-problème (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeByIterConvexeBase.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeLinearBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearBase[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeLinearBase.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeLinearBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterLinearBase[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Nombre de sous-problème moyen étudié")
	ax.grid(true)
	fig.savefig("SaveFig/moyNbIterLinearBase.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeLinearBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearBase[Int.((solX ./ 10) .- 4)] ./ moyNbIterLinearBase[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen d'un sous-problème (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeByIterLinearBase.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeCombo[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeCombo[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeCombo.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeCombo[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterCombo[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Nombre de sous-problème moyen étudié")
	ax.grid(true)
	fig.savefig("SaveFig/moyNbIterCombo.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeCombo[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeCombo[Int.((solX ./ 10) .- 4)] ./ moyNbIterCombo[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen d'un sous-problème (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeByIterCombo.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeComboHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeComboHeuristic[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeComboHeuristic.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeComboHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterComboHeuristic[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Nombre de sous-problème moyen étudié")
	ax.grid(true)
	fig.savefig("SaveFig/moyNbIterComboHeuristic.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeComboHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeComboHeuristic[Int.((solX ./ 10) .- 4)] ./ moyNbIterComboHeuristic[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen d'un sous-problème (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeByIterComboHeuristic.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeLinearHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearHeuristic[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeLinearHeuristic.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeLinearHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterLinearHeuristic[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Nombre de sous-problème moyen étudié")
	ax.grid(true)
	fig.savefig("SaveFig/moyNbIterLinearHeuristic.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeLinearHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearHeuristic[Int.((solX ./ 10) .- 4)] ./ moyNbIterLinearHeuristic[Int.((solX ./ 10) .- 4)]
	
	ax.plot(solX, solY, "b.-")
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen d'un sous-problème (en secondes)")
	ax.grid(true)
	fig.savefig("SaveFig/moyTimeByIterLinearHeuristic.png")
	close(fig)
	
	###############################################################################""
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeConvexeBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeConvexeBase[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Convexe JuMP")
	
	solX = filter(x->moyTimeLinearBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearBase[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Relaxation Linéaire")
	
	solX = filter(x->moyTimeCombo[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeCombo[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Combo")
	
	solX = filter(x->moyTimeComboHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeComboHeuristic[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Combo + Heuristique Primale")
	
	solX = filter(x->moyTimeLinearHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearHeuristic[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Relaxation linéaire + Heuristique Primale")
	
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen (en secondes)")
	ax.grid(true)
	ax.legend()
	fig.savefig("SaveFig/moyTimeGeneral.png")
	close(fig)
	
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyNbIterConvexeBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterConvexeBase[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Convexe JuMP")
	
	solX = filter(x->moyNbIterLinearBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterLinearBase[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Relaxation Linéaire")
	
	solX = filter(x->moyNbIterCombo[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterCombo[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Combo")
	
	solX = filter(x->moyNbIterComboHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterComboHeuristic[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Combo + Heuristique Primale")
	
	solX = filter(x->moyNbIterLinearHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyNbIterLinearHeuristic[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Relaxation linéaire + Heuristique Primale")
	
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Nombre de sous-problème étudié")
	ax.grid(true)
	ax.legend()
	fig.savefig("SaveFig/moyNbIterGeneral.png")
	close(fig)
	
	fig = figure()
	ax = fig.add_subplot(111)
	
	solX = filter(x->moyTimeConvexeBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeConvexeBase[Int.((solX ./ 10) .- 4)] ./ moyNbIterConvexeBase[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Convexe JuMP")
	
	solX = filter(x->moyTimeLinearBase[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearBase[Int.((solX ./ 10) .- 4)] ./ moyNbIterLinearBase[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Relaxation Linéaire")
	
	solX = filter(x->moyTimeCombo[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeCombo[Int.((solX ./ 10) .- 4)] ./ moyNbIterCombo[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Combo")
	
	solX = filter(x->moyTimeComboHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeComboHeuristic[Int.((solX ./ 10) .- 4)] ./ moyNbIterComboHeuristic[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Combo + Heuristique Primale")
	
	solX = filter(x->moyTimeLinearHeuristic[Int(x/10 - 4)]!=0, 50:10:150)
	solY = moyTimeLinearHeuristic[Int.((solX ./ 10) .- 4)] ./ moyNbIterLinearHeuristic[Int.((solX ./ 10) .- 4)]
	ax.plot(solX, solY, ".-", label = "Relaxation linéaire + Heuristique Primale")
	
	ax.set_xlabel("Nombre d'objets")
	ax.set_ylabel("Temps de résolution moyen d'un sous-problème (en secondes)")
	ax.grid(true)
	ax.legend()
	fig.savefig("SaveFig/moyTimeByIterGeneral.png")
	close(fig)
	
end
