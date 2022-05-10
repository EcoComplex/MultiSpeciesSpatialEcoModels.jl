## Initialization
using MultiSpeciesSpatialEcoModels
using Plots

## Initialize model (alt-enter runs all the code between ##)

fil = 100
col = 100
numsp = 4
densIni = fill(0.1, numsp)
densIni = [ 0.1 * i for i in 1:numsp]

λ = fill(1.7,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)
λ = fill(1,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
λ[1]=2
m = init_multispecies(fil,col,densIni)
calc_species_density(m,numsp)
plot_multispecies(m,4)

##
# run 100 steps and return the densities
#
vs = run_birthdeath!(m,λ,δ,100)

## Plot time series
#
plot(vs)

# Se puede plotear luego de 100 pasos
#
plot_multispecies(m,length(densIni))

## 
# Initialize a 100 species model 
#
using Distributions
numsp = 100
fil = 100
col = 100

#  
#
λ = rand(Normal(1.7,.04),numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = rand(Normal(1,.02),numsp)

#  
#
λ = fill(1.7,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1,numsp)


densIni = fill(0.1,numsp)
densIni = [ 0.1 * i for i in 1:numsp]

m = init_multispecies(fil,col,densIni)

##  run 100 steps and generate a gif 
fint = 200
dens = zeros(Float64,fint,numsp)

anim = @gif for i in 1:fint
    plt1 = plot_multispecies(m,numsp)
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,.1),label=nothing) #label=permutedims(["sp $i" for i=1:numsp])
    dens[i,:] = calc_species_density(m,numsp)
    step_birthdeath!(m,λ,δ)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2))
end

## gif(anim, "macro10sp.gif")

# Number of species at the end of simulation
#
count(calc_species_density(m,numsp).>0)

##
# Plot each 10 steps --> time 1000
#
anim = @gif for i in 1:fint
    plt1 = plot_multispecies(m,numsp)
    plt2=  plot(range(start=0,step=10,length=fint),dens,xlim=(0,fint*10),ylim=(0,.1),label=nothing) #label=permutedims(["sp $i" for i=1:numsp])
    dens[i,:] = calc_species_density(m,numsp)
    run_birthdeath!(m,λ,δ,10)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2))
end

##
# Find the surviving species
#
surv = unique(filter(x->x>0,m))

# Check the parameters of surviving species
#
[ round(λ[i]/δ[i],digits=4) for i in surv]

# Check the parameters of extinct species
#
[ round(λ[i]/δ[i],digits=4) for i in 1:length(λ) if i ∉ surv ]