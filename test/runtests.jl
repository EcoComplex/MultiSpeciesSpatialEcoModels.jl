using MultiSpeciesSpatialEcoModels
using Test

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
