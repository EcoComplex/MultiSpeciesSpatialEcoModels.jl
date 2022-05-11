"""
    step_birthdeathcol!(landscape, λ,δ,α)

Performs a step of the model a birth-death-colonization multispecies model, where `landscape` is the model grid, and 
λ and δ are vectors of the birth and death rate parameters. These vectors must have the same length that also defines the number of competing 
species. The length of λ and δ is the number of species. 

The parameters λ and δ could change in any step of the model.


"""
function step_birthdeathcol!(landscape,λ,δ,α)
    if !(length(α)==length(λ)==length(δ))
        error("Lengths of λ,δ,α must be equal")
    end

    R = sum(λ) + sum(δ) + sum(α)

    λ = λ ./ R 
    δ = δ ./ R 
    α = α ./ R 

    # Vector with acummulated probabilities of each event 
    #
    acumProb = zeros(length(λ)+ length(δ)+length(α))

    j = 2 
    acumProb[1] =  λ[1]
    acumProb[2] = acumProb[1] + δ[1]
    acumProb[3] = acumProb[2] + α[1]


    for i in 4:3:length(acumProb)
        acumProb[i] +=  acumProb[i-1] + λ[j]
        acumProb[i+1] = acumProb[i] + δ[j]
        acumProb[i+2] = acumProb[i+1] + α[j]
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
            z = sp * 3
            rnd = rand()
            if acumProb[z-2] < rnd 
                if rnd ≤ acumProb[z-1]
                    #@info "Birth event $acumProb z=$z $rnd"

                    x=[i,j] .+  rand(neighborhood)
                    if !(x[1] > fil || x[1] <1 || x[2] > col || x[2]<1)
                        if(landscape[x[1],x[2]] == 0 )                         # if is empty sp reproduces!
                            landscape[x[1],x[2]] = sp
                        end
                    end
                elseif rnd ≤ acumProb[z]
                    #@info "Mortality event  $acumProb z=$z $rnd"
                    landscape[i,j]=0
                end
            end
        else 
            rnd = rand()
            for sp in 1:length(δ)
                z = sp * 3
                if acumProb[z] < rnd ≤ acumProb[z+1]
                    #@info "Colonization event  $acumProb z=$z $rnd"
                    if (i == fil || i==1 || j==col || j==1)           # Colonization only occurs at the borders
                       landscape[i,j]=sp
                    end
                end
            end
        end
    end
end

"""
    run_birthdeathcol!(landscape,λ,δ, steps )

Runs the model with `landscape` grid and parameters vectors λ δ α (birth death  and colonization rate) during `steps` time
The model must be initialized previously, it modifies the `landscape` matrix and returns a matrix of populations over time.
"""
function run_birthdeathcol!(landscape,λ,δ,α, steps )

    di = zeros(steps,length(δ))

    for j in 1:steps
        step_birthdeathcol!(landscape,λ,δ,α)
        di[j,:] = calc_species_density(landscape,length(δ))
    end

    return di
end
