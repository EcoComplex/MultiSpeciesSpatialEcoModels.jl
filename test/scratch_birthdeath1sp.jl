## Initialization
using Revise
using MultiSpeciesSpatialEcoModels
using Plots
using Statistics

## Initialize model (alt-enter runs all the code between ##)

fil = 100
col = 100
numsp = 1
densIni = fill(0.001, numsp)
λ = fill(2.5,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)
# λ = fill(1,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
# λ[1]=2
m = init_multispecies(fil,col,densIni)
calc_species_density(m,numsp)
plot_multispecies(m,4)

##
# run 100 steps and return the densities
#
vs = run_birthdeath!(m,λ,δ,200, vecindad=2.1)
vs = run_birthdeath!(m,λ,δ,200)

# Plot time series
#
plot(vs)
mean(vs)

# Se puede plotear luego de 100 pasos
#
plot_multispecies(m,numsp,[:black,:green])

## test power law dispersal
#
m = zeros(Int16,fil,col)
m[50:51,50:51] .= 1
plot_multispecies(m,4)
step_birthdeath!(m,λ,δ,vecindad=2.1)
plot_multispecies(m,4)
step_birthdeath!(m,λ,δ)
plot_multispecies(m,4)

## Testing patch size 
#
using Images
m = zeros(Int16,10,10)
m[1,1]=m[10,1]=m[1,10]=m[10,10]=1
plot_multispecies(m,4)
step_birthdeath!(m,λ,δ)
plot_multispecies(m,4)
# Detect patches
labels = label_components(m)            
plot_multispecies(labels,maximum(labels))
# Get Patch size
patch_size = component_lengths(labels)
# Skip first element, background = 0
maxS = maximum(patch_size[2:end])
findmax(patch_size[2:end])[2]+1


plot_multispecies(max_patch_highlight(m),2,[:black,:green])
#
## Create a gif with max patch dynamics and density 
#
densIni = fill(0.6,numsp)
m = init_multispecies(fil,col,densIni)
fint = 400
dens = zeros(Float64,fint,numsp)
λ = fill(3.0,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)

anim = @animate for i in 1:fint
    mmax = max_patch_highlight(m)
    plt1 = plot_multispecies(mmax,2,[:black,:green])
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,1),label=nothing, annotate=(100,1,string("λ: ", round(λ[1],digits=2)))) #label=permutedims(["sp $i" for i=1:numsp])
    dens[i,:] = calc_species_density(m,numsp)
    if i % 100 == 0
        λ[1] -= .5
    end
    # step_birthdeath!(m,λ,δ, vecindad=2.1)
    step_birthdeath!(m,λ,δ)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2))
end every 10
gif(anim,fps=2)

#
## Create a gif with max patch dynamics and max patch size
#
densIni = fill(0.6,numsp)
m = init_multispecies(fil,col,densIni)
fint = 400
dens = zeros(Float64,fint)
λ = fill(3.0,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)

anim = @animate for i in 1:fint
    mmax = max_patch_highlight(m)
    msize = max_patch_size(m)
    dens[i]=msize/(fil*col)
    plt1 = plot_multispecies(mmax,2,[:black,:green])
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,1),label=nothing, annotate=(100,1,string("λ: ", round(λ[1],digits=2)))) #label=permutedims(["sp $i" for i=1:numsp])
    dens[i,:] = calc_species_density(m,numsp)
    if i % 100 == 0
        λ[1] -= .5
    end
    # step_birthdeath!(m,λ,δ, vecindad=2.1)
    step_birthdeath!(m,λ,δ)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2))
end every 10
gif(anim,fps=2)



#
##  run 100 steps and generate a gif 
#
densIni = fill(0.6,numsp)
m = init_multispecies(fil,col,densIni)
fint = 400
dens = zeros(Float64,fint,numsp)
λ = fill(2.6,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)

anim = @gif for i in 1:fint
    plt1 = plot_multispecies(m,numsp,[:black,:green])
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,1),label=nothing, annotate=(100,1,string("λ: ", round(λ[1],digits=2)))) #label=permutedims(["sp $i" for i=1:numsp])
    dens[i,:] = calc_species_density(m,numsp)
    if i % 100 == 0
        λ[1] -= .5
    end
    step_birthdeath!(m,λ,δ, vecindad=2.1)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2))
end 
## gif(anim, "macro10sp.gif")

m = init_multispecies(fil,col,densIni)
fint = 400
dens = Float64[]
λ = fill(3.0,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)
densIni = fill(0.7,numsp)

for i in 1:fint
    plt1 = plot_multispecies(m,numsp,[:black,:green])
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,1),label=nothing, annotate=(100,1,string("λ: ", round(λ[1],digits=2)))) #label=permutedims(["sp $i" for i=1:numsp])
    if i % 100 == 0
        λ[1] -= 0.50
    end
    step_birthdeath!(m,λ,δ)
    push!(dens, calc_species_density(m,numsp)[1])
    # push!(pob,sum(mm)/N)         #
    if i % 50 == 0
        display(plot(plt1, plt2, layout = (1, 2)))
    end 
end 
## gif(anim, "macro10sp.gif")

## Calculate the mean density for different lambda
numsp=1
fint = 200
λ = fill(3.0,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)
densIni = fill(0.5,numsp)
dens = Float64[]

for lamb in 1.2:0.1:3
    @show lamb
    m = init_multispecies(fil,col,densIni)
    vs = run_birthdeath!(m,lamb,δ,fint)
    
    push!(dens,mean(vs[100:fint,1]))

end
plot(1.2:0.1:3,dens, label=nothing)

