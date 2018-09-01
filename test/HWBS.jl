using Distributions
using DataFrames
using Test
using RMInsurance
using Dierckx
using CSV


import RMInsurance.zb

### TESTS
curr_dir = dirname(@__FILE__())

println("start HWBS.jl...")

## All Data sourced from Swiss National Bank (www.snb.ch)
# Read the raw data for market prices.  The raw data only
# constains spot rates for different measurement dates
# (columns :xyyyy_mm, where yyyy is the year and mm the month)
# and the corresponding Durations (column :Duration).
df_chf_spots = CSV.read( curr_dir * "/HWBS_Input_CHF_Spot.csv",
                          header = true)
# contains various short rates ("SARON", "Overnight",
# "MM_3months") for different measurement dates
df_chf_shorts = CSV.read(curr_dir * "/HWBS_Input_CHF_Short.csv",
                          header = true)
# parameters for the Hull-White Black-Scholes model
hwpar = CSV.read(curr_dir * "/HWBS_Input_Parameters.csv",
                  header = true)
# Asset portfolio
assets = CSV.read( curr_dir * "/HWBS_Input_Portfolio.csv",
                    header = true)

# t₀ refers to the time at which we will base the spot curve
# it is given as a column identifyer for df_chf_spots_sparse
t₀ =  Symbol("2016-09")  #:x2016_09
# We add the 3month money market debt rate (Switzerland - CHF -
# Money market debt register claims of the Swiss Confederation,
# 3-month) as a proxi for the spot rate of duration 0.
# This is used for extrapolation to spot rates with
# duration < 1 year
spot₀ =
  df_chf_shorts[findall(x->x=="SARON",df_chf_shorts[:Type])[1], t₀]

df_spot_coarse = DataFrame()
df_spot_coarse[:t] = vcat(0.0, df_chf_spots[:Duration])
df_spot_coarse[:spot] = vcat(spot₀, df_chf_spots[t₀])



dict_hwpar =
  Dict{Symbol, Real}( Symbol( hwpar[i, :name]) =>
                      hwpar[i, :value] for i = 1:nrow(hwpar))

function spot_coarse2fine(nʸ, df_spot_coarse)
  δt = 1 / nʸ
  nᵀ = maximum(df_spot_coarse[:t]) * nʸ
  func_spot =
    Spline1D( convert(Array, df_spot_coarse[:t]),
              convert(Array, df_spot_coarse[:spot]))
  # Durations and corresponding spot rates
  DataFrame(t = collect(1:nᵀ) * δt,
            spot = func_spot(collect(1:nᵀ) * δt))
end



hw = HWBS(dict_hwpar, spot_coarse2fine(dict_hwpar[:nʸ],
                                       df_spot_coarse))

asset_id = Array{Symbol}(undef, nrow(assets))
for 𝑟𝑜𝑤 ∈ 1:nrow(assets)
  if assets[𝑟𝑜𝑤,:Type] in ["zb"]
    asset_id[𝑟𝑜𝑤] =
      Symbol(assets[𝑟𝑜𝑤,:Type] *"_" * string(assets[𝑟𝑜𝑤,:Maturity]))
  else
    asset_id[𝑟𝑜𝑤] = Symbol(assets[𝑟𝑜𝑤,:Type])
  end
end

n_scen = 100_000
t₁ = 1
V₀ = zeros(nrow(assets))            ## Value of assets at time 0
V₁ = zeros(nrow(assets), n_scen)    ## Value of assets at time 1
## loss in value from 0 to 1, discounted to time 0
loss = DataFrame()
## random variables for capital market model at time 1
cm₁ = rand_rw(hw, t₁, n_scen)
r₁ = cm₁[1,:]  ## random interest rates at time 1
S₁ = cm₁[2,:]  ## random index values at time 1

discount = zb(hw, 1) ## discount rate from 0 to 1

for 𝑟𝑜𝑤 ∈ 1:nrow(assets)
  if assets[𝑟𝑜𝑤, :Type] == "index"
    V₀[𝑟𝑜𝑤] = assets[𝑟𝑜𝑤, :Nominal]
    V₁[𝑟𝑜𝑤, :] = V₀[𝑟𝑜𝑤] * S₁
  else
    V₀[𝑟𝑜𝑤] =
      zb(hw, assets[𝑟𝑜𝑤, :Maturity]) * assets[𝑟𝑜𝑤, :Nominal]
    zb_tmp(r) = zb(hw, r, t₁, assets[𝑟𝑜𝑤, :Maturity])
    V₁[𝑟𝑜𝑤, :] = zb_tmp.(r₁) * assets[𝑟𝑜𝑤, :Nominal]
  end
  loss[asset_id[𝑟𝑜𝑤]] = V₀[𝑟𝑜𝑤] .- discount * V₁[𝑟𝑜𝑤,:]
end
loss[:total] = sum(loss[asset_id[𝑟𝑜𝑤]] for 𝑟𝑜𝑤 ∈ 1:nrow(assets))

es_assets = es(loss[:total], 0.99)


println("...end HWBS.jl")
