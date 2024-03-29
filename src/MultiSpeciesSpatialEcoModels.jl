module MultiSpeciesSpatialEcoModels

export set_neighborhood!
export get_neighborhood
export step_birthdeath!
export step_birthdeathcol!
export calc_species_density
export init_multispecies
export run_birthdeath!
export run_birthdeathcol!
export plot_multispecies
export max_patch_size
export max_patch_highlight

#
# Modelo estocastico espacial de una población
#
using Plots
using Statistics
using Distributions
using StatsBase
using Random
using Images

include("birthdeathcol.jl")  # step function for the birth-death-colonization model
include("maxpatchstats.jl")

#
# Vecindad por defecto de Von Neuman = 4 vecinos mas cercanos
#
neighborhood = ((0,1),(0,-1),(1,0),(-1,0) )

"""
    randvecino_2D(vecindad::Tuple)

Una vecindad aleatoria en coordenadas 2D basado en una Tuple con coordenadas
para otros tipos de vecindades debería tomar una función
"""
function randvecino_2D(vecindad::Tuple)
    rand(vecindad)    
end



function set_neighborhood!(nbh)
    neighborhood = nbh
end

function get_neighborhood()
    return neighborhood 
end


"""
    randvecino_2D(alpha)

A random neighborhood in 2D coordinates based on an inverse potential distribution 
with x_min=1 and exponent alpha, equivalent to a Pareto distribution with exponent alpha-1

"""
function randvecino_2D(alpha::AbstractFloat)
    r = rand(Pareto(alpha-1,1))
    theta = rand()*360
    x , y = Int(round(r*cos(theta))), Int(round(r*sin(theta)))
end

    


"""
    step_birthdeath!(landscape,λ,δ;vecindad=neighborhood)

Performs a step of the model a birth-death multispecies model, where `landscape` is the model grid, and 
λ and δ are vectors of the birth and death rate parameters. These vectors must have the same length that also defines the number of competing 
species. The length of λ and δ is the number of species. The parameter `vecindad` determines the dispersal, if it is a tuple uses a random 
position withing the tuple as dispersal neighborhood, if it is a real it uses a power-law dispersal with exponent `vecindad` 

The parameters λ and δ could change in any step of the model.


"""
function step_birthdeath!(landscape,λ,δ;vecindad=neighborhood)
    if length(λ) != length(δ)
        error("Length of λ must be = to length of δ")
    end

    R = sum(λ) + sum(δ) 

    λ = λ ./ R 
    δ = δ ./ R 

    # Vector with acummulated probabilities of each event 
    #
    acumProb = zeros(length(λ)+ length(δ))

    j = 2 
    acumProb[1] =  λ[1]
    acumProb[2] = acumProb[1] + δ[1]

    for i in 3:2:length(acumProb)
        acumProb[i] +=  acumProb[i-1] + λ[j]
        acumProb[i+1] = acumProb[i] + δ[j]
        j+=1
    end

    # Check that for 3 species the computations was done right
    #
    # (λ[1],λ[1]+δ[1],λ[1]+δ[1]+λ[2], λ[1]+δ[1]+λ[2]+δ[2], λ[1]+δ[1]+λ[2]+δ[2]+λ[3], λ[1]+δ[1]+λ[2]+δ[2]+λ[3]+δ[3]) .≈ acumProb
    fil = size(landscape,1)
    col = size(landscape,2)
    N = fil*col*R
    pushfirst!(acumProb,0.0)
    for z in 1:N
        i= rand(1:fil) 
        j= rand(1:col)
        #
        # For each site any event could happen but if no species are present no event will happen
        #
        if landscape[i,j]!=0

            sp = landscape[i,j]
            #@info "Landscape $i $j not $sp"
            z = sp * 2
            rnd = rand()
            if acumProb[z-1] < rnd 
                if rnd ≤ acumProb[z]
                    #@info "Birth event $acumProb z=$z $rnd"
                    x=(i,j) .+  randvecino_2D(vecindad)

                    if !(x[1] > fil || x[1] <1 || x[2] > col || x[2]<1)
                        if(landscape[x[1],x[2]] == 0 )                       # if is empty sp reproduces!
                            landscape[x[1],x[2]] = sp
                        end
                    end
                elseif rnd ≤ acumProb[z+1]
                    #@info "Mortality event  $acumProb z=$z $rnd"
                    landscape[i,j]=0
                end
            end
        end
    end
end






"""
    calc_species_density(landscape, numSpecies) 

Calculate the density of the species in the `landscape` matrix
"""
function calc_species_density(landscape, numSpecies) 
    return proportions(landscape, 1:numSpecies )
end

"""
    plot_multispecies(landscape,numsp)

plot the `landscape` matrix of a multispecies model as a heatmap with `numsp` species
"""
function plot_multispecies(landscape,numsp) 
    heatmap(landscape,aspect_ratio=:equal,legend=:none,xticks=:none,yticks=:none,framestyle=:none,clims=(0,numsp))
end

"""
    plot_multispecies(landscape,numsp,colors)

plot the `landscape` matrix of a multispecies model as a heatmap with `numsp` species, with a color palette defined by `colors`
"""
function plot_multispecies(landscape,numsp,colors) 
    c = palette(colors , numsp+1)

    heatmap(landscape,aspect_ratio=:equal,legend=:none,xticks=:none,yticks=:none,framestyle=:none,clims=(0,numsp),color=c)
end



"""
    init_multispecies(fil, col, densidadIni)

Initialize a multispecies model landscape of size `fil x col` with densities and number of species given by the vector `densidadIni`
and returns a matrix of type `Int16`.
"""
function init_multispecies(fil,col,densidadIni)

    landscape = zeros(Int16,fil,col)
    numsp = length(densidadIni)
    numIni = zeros(Int32,numsp)
    #@info "Sum densidad $(sum(densidadIni))"
    if sum(densidadIni) > 1
        densidadIni = densidadIni ./ sum(densidadIni)
        numIni = rand(Multinomial(fil*col,densidadIni))
    else
        # Si hay espacio vacio hay que agregarlo al vector densidadIni para usar la distribucion Multinomial
        #
        vacio = 1 - sum(densidadIni)
        push!(densidadIni,vacio)
        numIni = rand(Multinomial(fil*col,densidadIni))
        #
        # Luego borramos el ultimo elemento de ambos vectores
        #
        pop!(numIni)
        pop!(densidadIni)
    end
    #@info "Sum densidad $(sum(densidadIni)) length $(length(densidadIni))"
    
    #@info "numIni $numIni"

    mins = 0
    maxs = 0
    for s in 1:numsp
        mins = maxs + 1
        maxs += numIni[s]
        
        for m in mins:maxs
            landscape[m] = s                   # Especie s
        end
    end
    shuffle!(landscape)
    return landscape
end


"""
    run_birthdeath!(landscape,λ,δ, steps; vecindad=neighborhood )

Runs the model with `landscape` grid and parameters vectors λ and δ (birth and death rate) during `steps` time
The model must be initialized previously, it modifies the `landscape` matrix and returns a matrix of populations over time.
"""
function run_birthdeath!(landscape,λ,δ, steps; vecindad=neighborhood)
    if length(λ) != length(δ)
        error("Length of λ must be = to length of δ")
    end

    di = zeros(steps,length(δ))

    for j in 1:steps
        step_birthdeath!(landscape,λ,δ;vecindad)

        di[j,:] = calc_species_density(landscape,length(δ))
    end

    return di
end


    
end