using Distributions
using Random
using KernelDensity
using LinearAlgebra
using Dates
using CSV
using Tables
using ColorSchemes
using Plots


function run_tau_leaping(τ, sdC::Float64 = 0.0, sdF::Float64 = 0.0, cor::Float64 = 0.0)

    years = 1000

    tmax = years*KC
    toss = years*0.25*KC
    N = zeros(Int64(years*KC/τ) + 1,2)
    t = zeros(Int64(years*KC/τ) + 1,1)
    rCsave = zeros(Int64(years*KC/τ) + 1,1)
    rFsave = zeros(Int64(years*KC/τ) + 1,1)

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

    C0 = round(C_eq, digits = (ndigits(KC)-1)) 
    F0 = round(F_eq, digits = (ndigits(KC)-1)) 
   
    #C0 = round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    #F0 = round(rand(d)[1], digits = (ndigits(KC)-1))::Float64
    
    N[1,:] = [C0, F0]

    row = 1::Int64

    event_C = [1/KC, -1/KC, 0.0, 0.0, 0.0]::Vector{Float64}
    event_F = [0.0, 0.0, 1/KF, -1/KF, -1/KF]::Vector{Float64}

    #event_C = [1/KC, -1/KC, 0.0, 0.0]::Vector{Float64}
    #event_F = [0.0, 0.0, 1/KF, -1/KF]::Vector{Float64}

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
        #rates = [birthC, deathC, birthF, deathFe]::Vector{Float64}

        #rates = max.(rates, 0)

        counts = [rand(Poisson(λ * τ)) for λ in rates]

        updates = vec(sum(counts .* event,dims = 1))

        N[row + 1, :] = N[row, :] .+ updates

        if N[row + 1,1] <= 0
            N[row + 1,1] = 0.0
        end

        if N[row + 1,2] <= 0
            N[row + 1,2] = 0.0
        end

        if N[row + 1,1] == 0.0 && N[row + 1, 2] == 0.0
            extinct_counter = extinct_counter + 1
        end


        t[row + 1] = t[row] + τ
        rCsave[row] = rCyr
        rFsave[row] = rFyr

        row = row + 1

    end

    sim_end = row - 1

    sim_start = 1 #argmin(abs.(t .- t[row]/2))[1]

    N2 = N[sim_start:sim_end,:]
    t = t[sim_start:sim_end]
    rCsave = rCsave[sim_start:sim_end,:]
    rFsave = rFsave[sim_start:sim_end,:]

    return(N2, t, C0, F0, extinct_counter, rCsave, rFsave)
end

d = Uniform(0.2,0.6)

γ = 8.5 # stress
α = 2.0 # interspecific competition aKC/KF
θ = 40 # facilitation

rC = 1.0 #1.1
rF = 2.2 #2.5

rF/rC

dC = γ
dF = γ

F_eq = (θ*(1- α) - 1 + sqrt((1 + θ*(α - 1))^2 - 4*θ*(α - 1 + dF/rF - α*dC/rC)))/(2*θ)
C_eq = 1 - dC/(rC*(1 + θ*F_eq))

KC = 10000
KF = 10000

stdev = 0.1
correlation = -0.5

N2, t, _, _, _, rCt, rFt = run_tau_leaping(50, stdev, stdev, correlation);

count(iszero,N2)

p1 = plot(t,N2[:,1], label = "Competitor", color = :blue);
plot!(p1, t,N2[:,2], label = "Facilitator", color = :orange)

maximum(abs.(diff(N2[:,2])))/mean(N2[:,2])
maximum(abs.(diff(N2[:,1])))/mean(N2[:,1])

# tau needs to be small enough so that the propensities don't change over the time period
# only supposed to change 0.01-0.03

CSV.write("data/tau100_ann/stoch_"*string(stdev)*"_N2.csv", Tables.table(N2))
CSV.write("data/tau100_ann/stoch_"*string(stdev)*"_t.csv", Tables.table(t))
#CSV.write("data/tau100_ann/stoch_"*string(stdev)*"_rCt.csv", Tables.table(rCt))
#CSV.write("data/tau100_ann/stoch_"*string(stdev)*"_rFt.csv", Tables.table(rFt))

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




clip = 1:1000
p1 = plot(t[clip],N2[clip,1], ylabel = "Scaled density", label = "Competitor", color = :blue);
plot!(p1, t[clip],N2[clip,2], label = "Facilitator", color = :orange)
vline!(p1, [10000:10000:100000], label = "", color = :gray)

p2 = plot(t[clip], rCt[clip], ylabel = "rC", label = "rC", )
p3 = plot(t[clip], rFt[clip], ylabel = "rF", label = "rF")

p0 = []
push!(p0, p1)
push!(p0, p2)
push!(p0, p3)
plot(p0..., layout = (3,1), size = (800, 600))
savefig("annual_stochasticity_tau100.pdf")

