module MultiSpeciesSpatialEcoModels

export set_neighborhood
export step_birthdeath!
export calc_species_density
export init_multispecies
export run_birthdeath!
export plot_multispecies
#
# Modelo estocastico espacial de una población
#
using Plots
using Statistics
using Distributions
using StatsBase
using Random

# Vecindad por defecto 
#
neighborhood = ((0,1),(0,-1),(1,0),(-1,0) )

function set_neighborhood(nbh)
    neighborhood = nbh
end

"""
    step_birthdeath!(landscape, λ,δ)

Performs a step of the model a birth-death multispecies model, where `landscape` is the model grid, and 
λ and δ are vectors of the birth and death rate parameters. These vectors must have the same length that also defines the number of competing 
species. The length of λ and δ is the number of species. 

The parameters λ and δ could change in any step of the model.


"""
function step_birthdeath!(landscape,λ,δ)
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

                    x=[i,j] .+  rand(neighborhood)
                    if !(x[1] > fil || x[1] <1 || x[2] > col || x[2]<1)
                        landscape[x[1],x[2]] = sp
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

# Una funcion para Inicializar
# Pensar la escala del sitio 30x30 cm ponele!
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
    run_birthdeath!(landscape,λ,δ, steps )

Runs the model with `landscape` grid and parameters vectors λ and δ (birth and death rate) during `steps` time
The model must be initialized previously, it modifies the `landscape` matrix and returns a matrix of populations over time.
"""
function run_birthdeath!(landscape,λ,δ, steps )
    if length(λ) != length(δ)
        error("Length of λ must be = to length of δ")
    end

    di = zeros(steps,length(δ))

    for j in 1:steps
        step_birthdeath!(landscape,λ,δ)
        di[j,:] = calc_species_density(landscape,length(δ))
    end

    return di
end



    
end