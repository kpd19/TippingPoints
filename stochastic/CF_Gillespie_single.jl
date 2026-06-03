using Distributions
using Random
using KernelDensity
using LinearAlgebra
using Dates
using CSV
using Tables
using ColorSchemes
using Plots

cd("/Users/katherinedixon/Documents/StuffINeed/_Research/UTokyo/Code/")

d = Uniform(0.1,0.9)

γ = 2.5 # stress
α = 2.0 # interspecific competition aKC/KF
θ = 20 # facilitation

rC = 1.0 #1.1
rF = 2.2 #2.5

rF/rC # rho

dC = γ
dF = γ

KC = 10000
KF = 10000

F_eq = (θ*(1- α) - 1 + sqrt((1 + θ*(α - 1))^2 - 4*θ*(α - 1 + dF/rF - α*dC/rC)))/(2*θ)
C_eq = 1 - dC/(rC*(1 + θ*F_eq))

function run_gillespie(sdC::Float64 = 0.0, sdF::Float64 = 0.0, cor::Float64 = 0.0)

    tmax = 1000*KC
    toss = 250*KC
    N = zeros(5*tmax,2)
    t = zeros(5*tmax,1)
    rCsave = zeros(5*tmax,1)
    rFsave = zeros(5*tmax,1)

    μ = [0.0, 0.0]::Vector{Float64}
    σ = [sdC, sdF]::Vector{Float64} 
    ρ = [1.0 cor;
         cor 1.0]::Matrix{Float64}

    cov_matrix = Symmetric(Diagonal(σ) * ρ * Diagonal(σ))

    if !isposdef(cov_matrix)
        cov_matrix = cov_matrix + I * 1e-8
    end

    cor_norm = MvNormal(μ,cov_matrix)

    C0 = C_eq #round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    F0 = F_eq #round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    N[1,:] = [C0, F0]

    row = 1::Int64

    event_C = [1/KC, -1/KC, 0.0, 0.0, 0.0]::Vector{Float64}
    event_F = [0.0, 0.0, 1/KF, -1/KF, -1/KF]::Vector{Float64}

    event = hcat(event_C, event_F)::Matrix{Float64}

    extinct_counter = 0::Int64

    while t[row] < tmax && extinct_counter <= 1 #round((now()-runtime).value/(60000),digits = 2) < timeout && birthC >= 0 && birthF >=0 #&& minimum(N[row,:]) > 0 
        
        stoch = rand(cor_norm, 1)::Matrix{Float64}

        rCt = rC*exp(stoch[1,1])::Float64
        rFt = rF*exp(stoch[2,1])::Float64

        birthC = rCt*N[row,1]*(1 - N[row,1])::Float64
        deathC = dC*N[row,1]/(1 + θ*N[row,2])::Float64

        birthF = rFt*N[row,2]*(1 - N[row,2])::Float64
        deathFc = rFt*N[row,2]* α*N[row,1]::Float64
        deathFe = dF*N[row,2]/(1 + θ*N[row,2])::Float64

        ak = [birthC, deathC, birthF, deathFc, deathFe]::Vector{Float64}
        atot = sum(ak) # total propensity

        τ = -1/atot*log(rand(Uniform(0,1)))::Float64

        r2 = rand(Uniform(0,1))::Float64

        idx = findfirst(x -> x == 1, 1/atot*cumsum(ak) .> r2)::Int64

        N[row + 1,:] = N[row,:] + event[idx,:]

        if N[row + 1,1] <= 0
            N[row + 1,1] = 0.0
            extinct_counter = extinct_counter + 1     
        end

        if N[row + 1,2] <= 0
            N[row + 1,2] = 0.0
            extinct_counter = extinct_counter + 1
        end

        t[row + 1] = t[row] + τ
        rCsave[row + 1] = rCt
        rFsave[row + 1] = rFt

        row = row + 1

    end

    sim_end = row-1

    sim_start = argmin(abs.(t .- t[row]/2))[1]

    N2 = N[sim_start:sim_end,:]
    t = t[sim_start:sim_end]
    rCsave = rCsave[sim_start:sim_end,:]
    rFsave = rFsave[sim_start:sim_end,:]

    return(N2, t, C0, F0, extinct_counter, rCsave, rFsave)
end

N2, t, _, _, _, rCt, rFt = run_gillespie(0.5,0.5, 0.0);

count(iszero,N2)

p1 = plot(t,N2[:,1], label = "Competitor", color = :blue, ylimits = [0,0.7]);
plot!(p1, t,N2[:,2], label = "Facilitator", color = :orange)

CSV.write("data/stoch_0.5_N2.csv", Tables.table(N2))
CSV.write("data/stoch_0.5_t.csv", Tables.table(t))
CSV.write("data/stoch_0.5_rCt.csv", Tables.table(rCt))
CSV.write("data/stoch_0.5_rFt.csv", Tables.table(rFt))


nplt = plot(N2[:,1], N2[:,2], xlabel = "Competitor", ylabel = "Facilitator", color = :black, legend = false);

dataC = N2[:,1]
dataF = N2[:,2]

xvals = range(0, stop=1, length=100) 
yvals = range(0, stop=1, length=100) 

k2 = kde((dataF, dataC), (xvals, yvals))

density_2d = k2.density # Matrix of values

heat = heatmap(xvals, yvals, density_2d, 
        color = cgrad(:roma, rev = true),
        xlabel = "Competitor", ylabel = "Facilitator")


sum(density_2d)*diff(xvals)[1]*diff(yvals)[1]

# histC = histogram(N2[:,1], color = :blue)
# histF = histogram(N2[:,2], color = :orange)
# histt = histogram(diff(t), color = :pink)

p0 = []
push!(p0, p1)
push!(p0, nplt)
push!(p0, heat)

plot(p0..., layout = (1,3), size = (1200, 400))
savefig("gillespie_example5.pdf")


