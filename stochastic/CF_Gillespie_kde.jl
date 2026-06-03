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

function run_gillespie(sdC::Float64 = 0.0, sdF::Float64 = 0.0, cor::Float64 = 0.0)

    tmax = 1000*KC
    toss = 250*KC
    N = zeros(5*tmax,2)
    t = zeros(5*tmax,1)
    
    μ = [0.0, 0.0]::Vector{Float64}
    σ = [sdC, sdF]::Vector{Float64} 
    ρ = [1.0 cor;  
         cor 1.0]::Matrix{Float64}

    cov_matrix = Symmetric(Diagonal(σ) * ρ * Diagonal(σ))

    if !isposdef(cov_matrix)
        cov_matrix = cov_matrix + I * 1e-8
    end

    cor_norm = MvNormal(μ,cov_matrix)

    C0 = round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    F0 = round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    N[1,:] = [C0, F0]

    row = 1::Int64
    #timeout = 10
    #runtime = now()

    #event_C = [1/KC, -1/KC, 0.0, 0.0, 0.0, 1/KC, 0.0]::Vector{Float64}
    #event_F = [0.0, 0.0, 1/KF, -1/KF, -1/KF, 0.0, 1/KF]::Vector{Float64}

    event_C = [1/KC, -1/KC, 0.0, 0.0, 0.0]::Vector{Float64}
    event_F = [0.0, 0.0, 1/KF, -1/KF, -1/KF]::Vector{Float64}

    event = hcat(event_C, event_F)::Matrix{Float64}

    #birthC = 1.0
    #birthF = 1.0

    extinct_counter = 0::Int64

    while t[row] < tmax && extinct_counter <= 1 #round((now()-runtime).value/(60000),digits = 2) < timeout && birthC >= 0 && birthF >=0 #&& minimum(N[row,:]) > 0 
        
        stoch = rand(cor_norm, 1)::Matrix{Float64}

        rCt = rC*exp(stoch[1,1])::Float64
        rFt = rF*exp(stoch[2,1])::Float64

        birthC = rCt*N[row,1]*(1 - N[row,1])::Float64
        deathC = dC*N[row,1]/(1 + θ*N[row,2])::Float64
        #iC = δ::Float64

        birthF = rFt*N[row,2]*(1 - N[row,2])::Float64
        deathFc = rFt*N[row,2]* α*N[row,1]::Float64
        deathFe = dF*N[row,2]/(1 + θ*N[row,2])::Float64
        #iF = δ::Float64

        ak = [birthC, deathC, birthF, deathFc, deathFe]::Vector{Float64}
        # ak = [birthC, deathC, birthF, deathFc, deathFe, iC, iF]::Vector{Float64}
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
        row = row + 1

    end

    sim_end = row-1 #findfirst(x -> x == 0, t[2:size(t)[1]]) 

    sim_start = argmin(abs.(t .- t[row]/2))[1]# Int64(round(row*0.5))

    if extinct_counter >= 1
        N2 = nothing
        t = nothing
    else
        N2 = N[sim_start:sim_end,:]
        t = t[sim_start:sim_end]
    end

    if isnothing(N2)
        println("Population went extinct- rerun")
    end

    return(N2, t, C0, F0, extinct_counter)
end

d = Uniform(0.2,0.6)

γ = 2.5 # stress
α = 2.0 # interspecific competition aKC/KF
θ = 20 # facilitation

rC = 1.0 #1.1
rF = 2.2 #2.5

rF/rC

dC = γ
dF = γ

KC = 10000
KF = 10000

stdev = 0.25
correlation = -0.5

atime = now()

NCs = Vector{Float64}()
NFs = Vector{Float64}()
ts = Vector{Float64}()
rlzn = Vector{Float64}()

for i in 1:10
    N2, t, _, _, _, = run_gillespie(stdev,stdev,correlation)

    if isnothing(N2)
        ntries = 0
        while isnothing(N2)
            N2, t, _, _, _, = run_gillespie(stdev,stdev,correlation)
            ntries = ntries + 1
        end
        println("Successful rerun after "*string(ntries)*" tries")
    end

    NCs = vcat(NCs, N2[:,1])
    NFs = vcat(NFs, N2[:,2])

    println("Finishing run "*string(i))

end

println("Finished in "*string(round((now()-atime).value/(60000),digits = 2))*" minutes")

xvals = range(0, stop=1, length=100); 
yvals = range(0, stop=1, length=100); 

dataC = NCs;
dataF = NFs;

k2 = kde((dataF, dataC), (xvals, yvals));

density_2d = k2.density; # Matrix of values

heat = heatmap(yvals, xvals, density_2d, 
        xlabel = "Competitor", ylabel = "Facilitator")

CSV.write("data/kde_case1_"*string(correlation)*"_g2.5_10k_r2.csv", Tables.table(density_2d))