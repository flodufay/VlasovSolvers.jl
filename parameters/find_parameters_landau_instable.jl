   using VlasovSolvers, SpecialFunctions


dev = CPU()                  # device
nx, nv = 12, 12              # grid resolution
stepper = StrangSplitting()  # timestepper
dt = 0.1                    # timestep
nsteps = 600                  # total number of time-steps

vmin, vmax = -10, 10           # V Domain length (m)
kx = 0.2               # Wave number of perturbation
xmin, xmax = 0, 2π/kx           # X Domain length (m)
v0 = 2.3
α = 0.00001         
t_start = dt       # Optimization Domain
ind_start = Int(t_start/dt)
window = [ind_start, Int(60/dt)]   # Display domain 
#a_book, b_book, c_book, phi_book, linear_book = -    0.00104, 1.1648, 0.002*0.42466, -0.33577, 0.02115
a_book, b_book, c_book, phi_book, linear_book = - 0.00242, 1.339, 0.002*0.42466, -0.33577, 0.2258
l = [[a_book, b_book, c_book, phi_book], [linear_book, 0, c_book, 0]]

T = LinRange(0, nsteps * dt, nsteps)                

    xgrid = OneDGrid(dev, nx, xmin, xmax)
    vgrid = OneDGrid(dev, nv, vmin, vmax)

    f = DistributionFunction( xgrid, vgrid )
    
    landau_instable!(f, α / 2, kx, v0)

    prob = VlasovProblem(f, Fourier(xgrid, vgrid), dev)

    sol = map(exp , solve!(prob, stepper, dt, nsteps ))


#pour initialiser automatiquement, il faut encore une initialisation pour trouver les 0.
function D(x)
    w = x[2] + im * x[1] 
    abs(1 + 1/(2*kx*kx) * (2 + (w/kx - v0) * (sqrt(pi/2) * exp(-(((w/kx)-v0)^2)/2) * ( im - erfi(((w/kx) - v0)/(sqrt(2))))) + (w/kx + v0) * sqrt(pi/2) * exp(-(((w/kx)+v0)^2)/2) * ( im - erfi(((w/kx)+ v0)/(sqrt(2)))) ))
end

find_parameters(l, sol, D,T, window, ind_start)