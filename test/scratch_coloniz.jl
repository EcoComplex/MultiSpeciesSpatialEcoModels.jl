## Initialization
using Revise
using MultiSpeciesSpatialEcoModels
using Plots

## Initialize model (alt-enter runs all the code between ##)

fil = 100
col = 100
numsp = 4
densIni = fill(0.0, numsp)
# densIni = [ 0.1 * i for i in 1:numsp]

λ = fill(2,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)
α = fill(0.05, numsp)

m = init_multispecies(fil,col,densIni)
calc_species_density(m,numsp)
plot_multispecies(m,4)
step_birthdeathcol!(m,λ,δ,α)
##
# run 100 steps and return the densities
#
vs = run_birthdeathcol!(m,λ,δ,α,200)
plot_multispecies(m,length(densIni))

## Plot time series
#
plot(vs)

## 
# Initialize a 100 species model 
#
using Distributions
numsp = 10
fil = 100
col = 100

#  
#
λ = rand(Normal(1.7,.1),numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = rand(Normal(1,.02),numsp)
α = rand(Normal(0.01,0.001), numsp)
α = fill(0.0, numsp)
all(α .> 0)
#  
#
densIni = fill(0.1,numsp)
# densIni = [ 0.1 * i for i in 1:numsp]

m = init_multispecies(fil,col,densIni)

##  run 100 steps and generate a gif 
fint = 400
dens = zeros(Float64,fint,numsp)

anim = @gif for i in 1:fint
    plt1 = plot_multispecies(m,numsp)
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,.16),label=nothing) #label=permutedims(["sp $i" for i=1:numsp])
    dens[i,:] = calc_species_density(m,numsp)
    step_birthdeathcol!(m,λ,δ,α)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2))
end

## gif(anim, "macro10sp.gif")

# Number of species at the end of simulation
#
count(calc_species_density(m,numsp).>0)
bar(calc_species_density(m,numsp))

## Compare with birthdeath model setting colonization to 0
#
numsp=10
λ = rand(Normal(1.7,.1),numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = rand(Normal(1,.02),numsp)
α = fill(0.0, numsp)
densIni = fill(0.1,numsp)
m = init_multispecies(fil,col,densIni)
vs = run_birthdeathcol!(m,λ,δ,α,fint)
count(calc_species_density(m,numsp).>0)
bar(calc_species_density(m,numsp))

## Find the surviving species

surv = unique(filter(x->x>0,m))

# Check the parameters of surviving species
#
[ round(λ[i]/δ[i],digits=4) for i in surv]

# Check the parameters of extinct species
#
[ round(λ[i]/δ[i],digits=4) for i in 1:length(λ) if i ∉ surv ]


