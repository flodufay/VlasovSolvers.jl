using VlasovSolvers

α = 0.01

print("\n one waves with α =", α, " : \n")

include("find_parameters_landau_stable.jl")

print("\n two waves with α =", α, " : \n")

include("find_parameters_landau_instable.jl")

print("\n three waves with α =", α, " : \n")

include("find_parameters_three_stream_instability.jl")


α = 0.001

print("\n one waves with α =", α, " : \n")

include("find_parameters_landau_stable.jl")

print("\n two waves with α =", α, " : \n")

include("find_parameters_landau_instable.jl")

print("\n three waves with α =", α, " : \n")

include("find_parameters_three_stream_instability.jl")

α = 0.0001

print("\n one waves with α =", α, " : \n")

include("find_parameters_landau_stable.jl")

print("\n two waves with α =", α, " : \n")

include("find_parameters_landau_instable.jl")

print("\n three waves with α =", α, " : \n")

include("find_parameters_three_stream_instability.jl")

α = 0.00001

print("\n one waves with α =", α, " : \n")

include("find_parameters_landau_stable.jl")

print("\n two waves with α =", α, " : \n")

include("find_parameters_landau_instable.jl")

print("\n three waves with α =", α, " : \n")

include("find_parameters_three_stream_instability.jl")
