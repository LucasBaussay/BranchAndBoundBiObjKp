using Profile, PProf
include("main.jl")
include("dataStruct.jl")

main("Inst50.dat")
@profile main("Inst50.dat")
pprof()
