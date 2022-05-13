using VlasovSolvers, SpecialFunctions


dev = CPU()                  # device
nx, nv = 128, 128              # grid resolution
stepper = StrangSplitting()  # timestepper
dt = 0.1                    # timestep
nsteps = 500                  # total number of time-steps

vmin, vmax = -10, 10           # V Domain length (m)
kx = 0.4               # Wave number of perturbation
xmin, xmax = 0, 2π/kx           # X Domain length (m)
α = 0.001         
t_start = 5          # Optimization Domain
ind_start = Int(t_start/dt)
window = [ind_start, Int(50/dt)]   # Display domain 
#a_book, b_book, c_book, phi_book = -0.1533, 1.4156, 0.004*0.3677 * sqrt(pi/kx), -0.536245
a_book, b_book, c_book, phi_book = -0.0661, 1.285, 0.002*0.4246, -0.3357
l = [[a_book, b_book, c_book, phi_book]]
a, b ,c, d = -0.5, 1.7, 0.01, -1


T = LinRange(0, nsteps * dt, nsteps)                

 #   sol = log(c) .+ a * T .+ log.(abs.(cos.(b * T .+ d)))
sol=c.*exp.(a*T).*(cos.(b * T .+ d))

#pour initialiser automatiquement, il faut encore une initialisation pour trouver les 0.
function D(x)
    w = x[2] + im * x[1] 
    abs(1 + 1/(kx*kx) * (1 + w/kx * sqrt(pi/2) * exp(-w*w/(2kx*kx)) * ( im - erfi(w/(sqrt(2) *kx)))))
end


find_parameters(l, sol, D, T, window, ind_start)