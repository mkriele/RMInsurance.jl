using RMInsurance
using DataFrames
using Test
using CSV

### TESTS
curr_dir = @__DIR__


println("start test ChainLadder_Test.jl...")
# We use the example from Mack's 1993-paper
# Cumulated claim amounts

## First code cell
cum_triangle = CSV.File(curr_dir * "/ChainLadderMack_Input.csv", 
                        header=false, ignoreemptylines=true) |> DataFrame
cum_triangle = Matrix(cum_triangle)
cum_triangle = convert.(Float64, cum_triangle)

## Next code cell
triangle =cum2claims(cum_triangle)      
@test cum_triangle ≈ claims2cum(triangle)

mack = Mack(cum_triangle, cum=true)

# Test: development factor
@test  [3.491, 1.747, 1.457, 1.174, 1.104, 1.086, 1.054,
        1.077,  1.018] ≈ round.(mack.f, digits=3)
# Test: future claims
@test sum(mack.futureclaims) ≈ mack.tot_res
# Test: standard error
@test [0.80, 0.26, 0.19, 0.27, 0.29, 0.26,
       0.22, 0.23, 0.29] ≈ round.(sqrt.(mack.mse) ./ mack.res, digits=2)[2:end]
# Test: total standard error
@test 0.13 ≈ round(√(mack.tot_mse)/mack.tot_res, digits = 2)

println("...end test ChainLadder_Test.jl")

