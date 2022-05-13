using FFTW, LinearAlgebra, Statistics, Plots, Optim, SpecialFunctions


export Parameters

function construct_optimum(l, T1)
    res = zeros(length(T1))
    i = 0
    while i <= length(l)- 4
        res += ((l[i + 3] .* exp.(l[i + 1] .* T1) .*abs.((cos.(l[i + 2] .* T1 .+ l[i + 4])))))
        i+= 4
    end
 #   res = log.(abs(l[3])) .+ l[1].*T1 .+ log.(abs.(cos.(l[2] .* T1 .+ l[4])))
    res
end

function concat(l :: Array{Array{Float64,1},1})
    res = l[1]
    if length(l) > 1
        for w in l[2:length(l)]
            res = [res;w]
        end
    end
    return res
end

function unconcat(l)
    res = Array{Array{Float64,1},1}(undef, 0)
    i = 0
    while i <= length(l)- 4
        res = append!(res, [[l[i+1], l[i+2], l[i+3], l[i+4]]])
        i +=4
    end
    return res
end



function optimum(sol, init :: Array{Array{Float64,1},1}, T, ind_start)
    function f(x) 
        return (construct_optimum(x, T[ind_start : length(T)]) .- sol[ind_start : length(T)])
    end

    function gradient!(G,x) 
        tab = f(x)
        i = 0
        d = 1
        T2 = T[ind_start : length(T)]
        while i <= length(x)- 4
            G[i + 1] = 0
            G[i + 2] = 0
            G[i + 3] = 0
            G[i + 4] = 0
            for t in 1 : length(T2)
                if (cos(x[i + 2] * T[2] + x[i + 4])) >= 0
                    d = 1
                else
                    d = -1
                end
                G[i + 1] += tab[t]/norm(tab) * T2[t] * x[i + 3] * exp(T2[t] * x[i + 1]) * abs(cos(x[i + 2] * T[2] + x[i + 4]))
                G[i + 2] += tab[t]/norm(tab) * (- 1 * d) * T2[t] * x[i + 3] * exp(T2[t] * x[i + 1]) * sin(x[i + 2] * T[2] + x[i + 4])
                G[i + 3] += tab[t]/norm(tab) * exp(T2[t] * x[i + 1]) * abs(cos(x[i + 2] * T[2] + x[i + 4]))
                G[i + 4] += tab[t]/norm(tab) * (- 1 * d) * x[i + 3] * exp(T2[t] * x[i + 1]) * sin(x[i + 2] * T[2] + x[i + 4])
                end
            i+= 4
        end
    end

    function dist(x) 
        return norm(f(x))
    end

#    res = optimize(dist, concat(init), ConjugateGradient())
    res = optimize(dist, concat(init))
 #res = optimize(dist, gradient!,concat(init),BFGS())
    print("actualisé \n")
    
    l, min = Optim.minimizer(res), Optim.minimum(res)
    return unconcat(l), min
end


function initialization(l :: Array{Array{Float64,1},1}, D, sol, ind_start)
    phi0 = 0
    c0 = sol[ind_start]/length(l)
    lower = [-10, -0.1]
    upper = [10, 10]
    init = copy(l)
    for i in 1: length(l)
        inner_optimizer = GradientDescent()
        opt = optimize(D, lower, upper, [l[i][1], l[i][2]], Fminbox(inner_optimizer))
        print("D(w) = ")
        print(Optim.minimum(opt))
        print("\n")
        init[i] = [Optim.minimizer(opt)[1], Optim.minimizer(opt)[2], c0, phi0]
    end
    return init
end

export find_parameters


function find_parameters(l :: Array{Array{Float64,1},1}, sol, D, T, window, ind_start, p = false)
    init = initialization(l, D, sol, ind_start)
    #init = l
    parameters, dist = optimum(sol, init, T, ind_start)

    if p 

        Y = construct_optimum(concat(parameters), T)
        T_tronc = T[window[1] : window[2]]
        sol_tronc = sol[window[1] : window[2]]
        Y_tronc = construct_optimum(concat(parameters), T_tronc)
        test_tampon = construct_optimum(concat(l), T_tronc)
        test_tampon = test_tampon * sol_tronc[1]/test_tampon[1]

        display(plot(T_tronc, [log.(abs.(sol_tronc)), log.(abs.(Y_tronc))] , label = ["simulation" "approximation of the simulation" ], legend = :topleft ))
    end
    print("\n our initialization : ")
    print(init)
    print("\n book initialization (adapted) : ")
    print(l)
    print("\n parameters result and distance to the simulation : ")
    print(parameters,", ", dist, "\n")
    return function f(t)
            res = 0
            for w in parameters
                res += ((w[3] * exp(w[1] * t) * abs(cos(w[2] * t + w[4]))))
            end
            return res
        end
    
end

