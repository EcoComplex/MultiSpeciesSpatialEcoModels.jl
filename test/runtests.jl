using MultiSpeciesSpatialEcoModels
using Test
using Statistics
using Random

Random.seed!(1707)
fil=100
col=100

@testset "Initialization 4 sp 0.0 density" begin
    densIni=fill(0.0, 4)
    m = init_multispecies(fil,col,densIni)
    spd = calc_species_density(m,length(densIni))
    @test spd == zeros(4)
end

@testset "Initialization 4 sp 0.7 density" begin
    densIni=fill(0.7, 4)
    m = init_multispecies(fil,col,densIni)
    spd = calc_species_density(m,length(densIni))
    @test sum(spd) == 1.0
end

@testset "Initialization 4 sp 0.1 density" begin
    densIni=fill(0.1, 4)
    m = init_multispecies(fil,col,densIni)
    spd = calc_species_density(m,length(densIni))
    
    @test round.(spd, digits=1) == densIni   # initialization
    
end

@testset "birthdeath 1 species survival, 3 extinct" begin
    numsp = 4
    densIni = fill(0.1, numsp)
    δ = fill(1, numsp)
    λ = fill(1,numsp)   # Only species 1 above the critical point λ/δ > 1.5
    λ[1]=2              # Species 1 should survive

    m = init_multispecies(fil,col,densIni)
    vs = run_birthdeath!(m,λ,δ,100)             # run 100 time steps 
    pop =mean(vs[50:100,:], dims=1)             # mean of the last 50 steps

    @test pop[1] > 0.3                          # species 1 density > 3
    @test all( pop[2:numsp] .== 0)              # other species densities == 0 
    
end

