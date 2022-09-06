## Initialization
using Revise
using MultiSpeciesSpatialEcoModels
using Plots
using Statistics

## Initialize model (alt-enter runs all the code between ##)

fil = 100
col = 100
numsp = 1

#
##  Animation Plot of density with \lambda 2.6 - 1.2  power-law dispersion 
#
densIni = fill(0.6,numsp)
m = init_multispecies(fil,col,densIni)
fint = 400
dens = zeros(Float64,fint,numsp)
λ = fill(2.6,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)

# Plot 2 steps 
#
plot_multispecies(m,numsp,[:black,:green])
step_birthdeath!(m,λ,δ)
plot_multispecies(m,numsp,[:black,:green])

# Generate animations
#
anim = @gif for i in 1:fint
    plt1 = plot_multispecies(m,numsp,[:black,:green])
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,1),label=nothing, annotate=(100,1,string("λ: ", round(λ[1],digits=2))))
    dens[i,:] = calc_species_density(m,numsp)
    if i % 100 == 0
        λ[1] -= .5
    end
    step_birthdeath!(m,λ,δ, vecindad=2.1)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2))
end 
#
# 
#
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

#
## Calculate the mean density for different lambda
#
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
    plot(plt1, plt2, layout = (1, 2),plot_title)
end every 10
gif(anim,fps=2)

#
## Create a gif with max patch dynamics and max patch size
#
densIni = fill(0.6,numsp)
m = init_multispecies(fil,col,densIni)
fint = 500
dens = zeros(Float64,fint)
λ = fill(3.5,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)

anim = @animate for i in 1:fint
    mmax = max_patch_highlight(m)
    msize = max_patch_size(m)
    plt1 = plot_multispecies(mmax,2,[:black,:green])
    plt2=  plot(dens,xlim=(0,fint),ylim=(0,1),label=nothing) #label=permutedims(["sp $i" for i=1:numsp])
    dens[i]=msize/(fil*col)
    if i % 100 == 0
        λ[1] -= .5
    end
    # step_birthdeath!(m,λ,δ, vecindad=2.1)
    step_birthdeath!(m,λ,δ)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2),plot_title=string("λ: ", round(λ[1],digits=2), " - Smax"))
end every 10
gif(anim, fps=2)
## gif(anim, "macro10sp.gif")
# "MaxPatchSizeLambda3-1.5.gif"

#
## Plot of delta Smax vs lambda 
#
densIni = fill(0.6,numsp)
m = init_multispecies(fil,col,densIni)
fint = 500
dens = zeros(Float64,fint)
smax = zeros(Float64,fint)
Δsmax= zeros(Float64,fint)
λ = fill(3.5,numsp)   # Todas las especies por encima del punto critico λ/δ > 1.5
δ = fill(1, numsp)
smax[1] = max_patch_size(m)/(fil*col)

anim = @animate for i in 2:fint
    mmax = max_patch_highlight(m)
    msize = max_patch_size(m)
    smax[i]=msize/(fil*col)
    dens[i]=calc_species_density(m,numsp)[1]
    Δsmax[i] = smax[i]-smax[i-1]
    plt1 = plot_multispecies(mmax,2,[:black,:green])
    #plt1=  plot(dens,xlim=(0,fint),ylim=(0,1),label=nothing, annotate=(100,1,string("λ: ", round(λ[1],digits=2)))) 
    plt2=  plot(Δsmax, xlim=(0,fint),ylim=(-0.4,0.4),label=nothing ) 
    if i % 100 == 0
        λ[1] -= .5
    end
    # step_birthdeath!(m,λ,δ, vecindad=2.1)
    step_birthdeath!(m,λ,δ)
    # push!(pob,sum(mm)/N)         #
    plot(plt1, plt2, layout = (1, 2),plot_title=string("λ: ", round(λ[1],digits=2), " - Δsmax"))
end every 10
gif(anim,fps=2)

# z = repeat([1, 2, 3, 4],inner=100)
# plot(Δsmax,group=z, cgrad=:viridis)