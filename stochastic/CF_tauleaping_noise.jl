using Distributions
using Random
using KernelDensity
using LinearAlgebra
using Dates
using CSV
using Tables
using ColorSchemes
using Plots
using DataFrames

cd("/Users/katherinedixon/Documents/StuffINeed/_Research/UTokyo/Code/")

function run_tau_noise(τ, sdC::Float64 = 0.0, sdF::Float64 = 0.0, cor::Float64 = 0.0)

    years = 1000

    tmax = years*KC
    toss = years*0.25*KC
    N = zeros(Int64(years*KC/τ) + 1,2)
    t = zeros(Int64(years*KC/τ) + 1,1)
    #rCsave = zeros(Int64(years*KC/τ) + 1,1)
    #rFsave = zeros(Int64(years*KC/τ) + 1,1)

    μ = [0.0, 0.0]::Vector{Float64}
    σ = [sdC, sdF]::Vector{Float64} 
    ρ = [1.0 cor;
         cor 1.0]::Matrix{Float64}

    cov_matrix = Symmetric(Diagonal(σ) * ρ * Diagonal(σ))

    if !isposdef(cov_matrix)
        cov_matrix = cov_matrix + I * 1e-8
    end

    cor_norm = MvNormal(μ,cov_matrix)

    stoch = rand(cor_norm, Int64(years))::Matrix{Float64}

    rCt = rC*exp.(stoch[1,:])::Vector{Float64}
    rFt = rF*exp.(stoch[2,:])::Vector{Float64}

    C0 = round(C_eq, digits = (ndigits(KC)-1)) #round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    F0 = round(F_eq, digits = (ndigits(KC)-1)) #round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    #C0 = round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    #F0 = round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    
    N[1,:] = [C0, F0]

    row = 1::Int64

    event_C = [1/KC, -1/KC, 0.0, 0.0, 0.0]::Vector{Float64}
    event_F = [0.0, 0.0, 1/KF, -1/KF, -1/KF]::Vector{Float64}

    event = hcat(event_C, event_F)::Matrix{Float64}

    extinct_counter = 0::Int64

    while t[row] < tmax #&& extinct_counter <= 500 #round((now()-runtime).value/(60000),digits = 2) < timeout && birthC >= 0 && birthF >=0 #&& minimum(N[row,:]) > 0 
        
        yr = Int(floor(t[row]/KC)) + 1

        rCyr = rCt[yr]
        rFyr = rFt[yr]

        birthC = rCyr*N[row,1]*(1 - N[row,1])::Float64
        deathC = dC*N[row,1]/(1 + θ*N[row,2])::Float64

        birthF = rFyr*N[row,2]*(1 - N[row,2])::Float64
        deathFc = rFyr*N[row,2]* α*N[row,1]::Float64
        deathFe = dF*N[row,2]/(1 + θ*N[row,2])::Float64

        rates = [birthC, deathC, birthF, deathFc, deathFe]::Vector{Float64}

        counts = [rand(Poisson(λ * τ)) for λ in rates]

        updates = vec(sum(counts .* event,dims = 1))

        N[row + 1, :] = N[row, :] .+ updates

        if N[row + 1,1] < 1/KC
            N[row + 1,1] = 0.0
            extinct_counter = extinct_counter + 1     
        end

        if N[row + 1,2] < 1/KF
            N[row + 1,2] = 0.0
            extinct_counter = extinct_counter + 1
        end

        t[row + 1] = t[row] + τ
        #rCsave[row] = rCyr
        #rFsave[row] = rFyr

        row = row + 1

    end

    sim_end = row - 1

    sim_start = 1#argmin(abs.(t .- t[row]/2))[1]

    # if extinct_counter >= 1
    #     N2 = N[sim_start:sim_end,:]
    #     t = t[sim_start:sim_end]
    # else
    #     N2 = N[sim_start:sim_end,:]
    #     t = t[sim_start:sim_end]
    # end

    N2 = N[sim_start:sim_end,:]
    t = t[sim_start:sim_end]

    end_t = size(N2)[1]
    buffer = 5000

    pop_avg = vec(mean(N2[(end_t - buffer):end_t,:], dims = 1))

    pop_end = N2[end_t,:]

    first_ext_C = findfirst(x -> x == 0, N2[:,1]) 
    first_ext_F = findfirst(x -> x == 0, N2[:,2]) 

    if isnothing(first_ext_C)
        if isnothing(first_ext_F)
            first_ext = maximum(t)/KC
        else
            first_ext = t[first_ext_F]/KC
        end
    else
        if isnothing(first_ext_F)
            first_ext = t[first_ext_C]/KC
        else
            first_ext = minimum([t[first_ext_C], t[first_ext_F]])/KC
        end
    end

    return(pop_avg, pop_end, first_ext)
end

#d = Uniform(0.2,0.6)

γ = 2.5 # stress
α = 2.0 # interspecific competition aKC/KF
θ = 20 # facilitation

if θ == 20
    case = "case1"
elseif θ == 40
    case = "case2"
else
    case = "unknown"
end

rC = 1.0 #1.1
rF = 2.2 #2.5

rF/rC

dC = γ
dF = γ

F_eq = (θ*(1- α) - 1 + sqrt((1 + θ*(α - 1))^2 - 4*θ*(α - 1 + dF/rF - α*dC/rC)))/(2*θ)
C_eq = 1 - dC/(rC*(1 + θ*F_eq))

KC = 10000
KF = 10000

stdev = 0.15
correlation = 0.0

atime = now()

extinction_df = DataFrame()
for sig in [0.0, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.40, 0.45, 0.5]
    for i in 1:1000
        pop_avg, pop_end, first_ext = run_tau_noise(100, sig, sig, correlation)

        temp = DataFrame(C_avg = pop_avg[1], F_avg = pop_avg[2], C_end = pop_end[1], F_end = pop_end[2], t_ext = first_ext, sigma = sig, cor = correlation, rep = i, d = γ, theta = θ)

        append!(extinction_df, temp)

        # if i % 10 == 0
        #     println("Finishing run "*string(i))
        # end
    end
    println("Finishing sigma = "*string(sig))
end

println("Finished in "*string(round((now()-atime).value/(60000),digits = 2))*" minutes")

CSV.write("data/extinctions/tau_"*case*"_"*string(correlation)*"_"*string(γ)*"_extinctions_1k.csv", extinction_df)


pop_avg, first_ext = run_tau_noise(100, 0.4,0.4,correlation)

p1 = plot(t,N2[:,1], label = "Competitor", color = :blue);
plot!(p1, t,N2[:,2], label = "Facilitator", color = :orange)

pop_avg
first_ext/KC

end_t = size(N2)[1]
buffer = 5000 
vec(mean(N2[(end_t - buffer):end_t,:], dims = 1))

p1 = plot(t[(end_t - buffer):end_t],N2[(end_t - buffer):end_t,1], label = "Competitor", color = :blue);
plot!(p1, t[(end_t - buffer):end_t],N2[(end_t - buffer):end_t,2], label = "Facilitator", color = :orange)
