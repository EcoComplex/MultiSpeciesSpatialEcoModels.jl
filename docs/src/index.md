# MultiSpeciesSpatialEcoModel.jl 

## Overview

This is an implementation of several models that represent individuals as a discrete units spatially located in a two dimensional grid. 
The most simple spatially explicit model in the context of ecology, is the contact process which is a straithforward implementation of the logistic model in 2D @Gastner2009, here we developed the multispecies version of this model with the additon of: 1) colonization 2) positive or negative interactions between species 3) Spatial heterogeneity in terms of one dimensional niche. 

## The multiple species contact model

An natural extension of the one species contact model is the multiple species contact model, which automatically adds competition for space to the system. Considering $n$ species denoted $S_i$ where $i=(1... n)$, the processes of birth and death could be represented by the following transitions:

$S_i + Ø \xrightarrow{\lambda_i} S_i + S_i$

$S_i \xrightarrow{\delta_i} Ø$

The first equation indicates that an individual of species $S_i$ gives birth to another individual on an empty site $Ø$ in a local neighborhod at a rate $\lambda_i$. The second equation shows that an individual $S_i$ dies at a rate $\delta_i$ to give an empty site $Ø$.

### Simulations

To simulate the birth-death multiple species contact model we could start initializing the model:

```julia

using MultiSpeciesSpatialEcoModels
using Plots

## Initialize model (alt-enter runs all the code between ##)

fil = 100               # number of rows
col = 100               # number of cols
numsp = 4               # number of species 

densIni = fill(0.1, numsp) # Initial density of species 

λ = fill(1.7,numsp)    # Groth rate
δ = fill(1, numsp)     # Death rate 

m = init_multispecies(fil,col,densIni)  # Generates a matrix (landscape) with species distributed at random 
                                        # and densities chosen from a Multinomial distribution

calc_species_density(m,numsp)           # calculates species densities 

plot_multispecies(m,4)                  # Plot the landscape using GR.heatmap

```

Then we could run the model and plot the time-series and the landscape

```julia

vs = run_birthdeath!(m,λ,δ,200)
plot(vs)

plot_multispecies(m,4)                  

```

Or we could generate a gif


```julia
fint = 400
m = init_multispecies(fil,col,densIni)

dens = zeros(Float64,fint,numsp)

anim = @gif for i in 1:fint
    plt1 = plot_multispecies(m,numsp)
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,.16),label=permutedims(["sp $i" for i=1:numsp])
    dens[i,:] = calc_species_density(m,numsp)
    step_birthdeath!(m,λ,δ)
    plot(plt1, plt2, layout = (1, 2))
end

```

The model is using absorving boundary conditions. 


# References

@Martin2020
@Black2012 
@Neuhauser1998 

@Levine2017

@Keymer2000
