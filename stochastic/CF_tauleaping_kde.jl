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

function run_tau_leaping(τ, sdC::Float64 = 0.0, sdF::Float64 = 0.0, cor::Float64 = 0.0)

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

    while t[row] < tmax && extinct_counter <= 1 #round((now()-runtime).value/(60000),digits = 2) < timeout && birthC >= 0 && birthF >=0 #&& minimum(N[row,:]) > 0 
        
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

    sim_start = argmin(abs.(t .- t[row]/2))[1]

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

    return(N2, t)
end

d = Uniform(0.2,0.6)

γ = 2.5 # stress
α = 2.0 # interspecific competition aKC/KF
θ = 20 # facilitation

if θ == 20
    case = "case1"
elseif θ == 40
    case == "case2"
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

stdev = 0.05
correlation = 0.9

atime = now()

NCs = Vector{Float64}()
NFs = Vector{Float64}()
ts = Vector{Float64}()
rlzn = Vector{Float64}()

for i in 1:50
    N2, t = run_tau_leaping(100, stdev,stdev,correlation)

    if isnothing(N2)
        ntries = 0
        while isnothing(N2)
            N2, t = run_tau_leaping(100, stdev,stdev,correlation)
            ntries = ntries + 1
        end
        println("Successful rerun after "*string(ntries)*" tries")
    end

    NCs = vcat(NCs, N2[:,1])
    NFs = vcat(NFs, N2[:,2])


    if i % 10 == 0
        println("Finishing run "*string(i))
    end

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

CSV.write("data/tau_"*case*"_"*string(stdev)*"_"*string(correlation)*"_g2.5_10k_r1.csv", Tables.table(density_2d))
