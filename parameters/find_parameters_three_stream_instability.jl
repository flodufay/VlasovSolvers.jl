using VlasovSolvers, SpecialFunctions


dev = CPU()                  # device
nx, nv = 96, 2 * 96              # grid resolution
stepper = StrangSplitting()  # timestepper
dt = 0.1                    # timestep
nsteps = 500                  # total number of time-steps

vmin, vmax = -10, 10           # V Domain length (m)
kx = 0.5               # Wave number of perturbation
xmin, xmax = 0, 2π/kx           # X Domain length (m)
v0 = 3.4
b = 0.2          
Tc = 0.1
t_start = 5       # Optimization Domain
ind_start = Int(t_start/dt)
window = [ind_start, Int(nsteps)]         # Display domain 
#a_book, b_book, c_book, phi_book, linear_book = -    0.00104, 1.1648, 0.002*0.42466, -0.33577, 0.02115
a_book, b_book, c_book, phi_book =  0.05, 1.06, 0.002*0.42466, 1.0894
l = [[a_book, b_book, c_book, phi_book]]

T = LinRange(0, nsteps * dt, nsteps)                

    xgrid = OneDGrid(dev, nx, xmin, xmax)
    vgrid = OneDGrid(dev, nv, vmin, vmax)

    f = DistributionFunction( xgrid, vgrid )
    
    three_stream_instability!(f, α, kx, v0, b, Tc)

    prob = VlasovProblem(f, Fourier(xgrid, vgrid), dev)

    sol = map(exp , solve!(prob, stepper, dt, nsteps ))


#pour initialiser automatiquement, il faut encore une initialisation pour trouver les 0.
function D(x)
    w = x[2] + im * x[1] 
    abs(1 + 1/(kx*kx) * (((1 - b)/Tc) * (1 + w/(kx * sqrt(Tc)) * sqrt(pi/2) * exp(-(((w/kx))^2)/(2 * Tc)) * ( im - erfi(((w/kx))/(sqrt(2 * Tc))))) + b/2 * (2 + (w/kx - v0) * (sqrt(pi/2) * exp(-(((w/kx)-v0)^2)/2) * ( im - erfi(((w/kx) - v0)/(sqrt(2))))) + (w/kx + v0) * sqrt(pi/2) * exp(-(((w/kx)+v0)^2)/2) * ( im - erfi(((w/kx)+ v0)/(sqrt(2)))))))
end

find_parameters(l, sol, D, T, window, ind_start)