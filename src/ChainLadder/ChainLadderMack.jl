using DataFrames
using LinearAlgebra

export Mack, claims2cum, cum2claims, cum2futureclaims, replacena!
"""
A reserve triangle with chain ladder information
"""
mutable struct Mack
  "Length of history"
  I::Int
  "Cummulative claims (upper triangle)"
  cum::Matrix{Float64}
  "Time vector of development factors"
  f::Vector{Float64}
  "Time vector of future claims payments"
  futureclaims::Vector{Float64}
  "Vector of reseves"
  res::Vector{Float64}
  """
  Total reserves: ``\\mathit{tot\\_res} =\\sum_{i=1}^I \\mathit{res}_i``
  """
  tot_res::Float64
  "Vector of mean square errors"
  mse::Vector{Float64}
  "Total mean square error"
  tot_mse::Float64
end


# function Mack(df_triang::DataFrame; cum::Bool = false)
# Mack(convert(DataArray{Real, 2}, df_triang); cum = cum)
# end
#
# function Mack(da_triang::DataArray{Real, 2}; cum::Bool = false)
#   Mack(convert(Array{Real,2}, da_triang, 0); cum = cum)
# end

function Mack(triang::Matrix{Float64}; cum::Bool = false) 
  if cum
    C = deepcopy(triang)
  else
    C = claims2cum(triang)
  end
  # Number of accident / development years
  I = size(C)[1]
  # Development factors
  f = [sum(C[1:(I-𝑘),𝑘+1])/sum(C[1:(I-𝑘),𝑘]) for 𝑘 ∈ 1:(I-1)]
  # Ultimate claim amounts
  for 𝑘  ∈ 2:I
    for 𝑙 ∈ (I+1-𝑘+1):I
      C[𝑘,𝑙] = C[𝑘,I+1-𝑘] * prod(f[(I+1-𝑘):( 𝑙 -1)])
    end
  end
  futclaims = cum2futureclaims(C)
  # Reserves
  R = [C[𝑘,I]-C[𝑘,I+1-𝑘] for 𝑘 ∈ 1:I]
  σ² =
    [1/(I-𝑘-1) *
      C[1:(I-𝑘),𝑘] ⋅ (C[1:(I-𝑘),𝑘+1] ./ C[1:(I-𝑘),𝑘] .- f[𝑘]).^2
    for 𝑘 ∈ 1:(I-1)]
  σ²[I-1] = min(σ²[I-2]^2 / σ²[I-3], min(σ²[I-3],σ²[I-2]))
  # mean square error
  mse = zeros(Float64, I)
  for 𝑖 ∈ 2:I
    mse[𝑖] = 0
    for 𝑘 ∈ (I+1-𝑖):(I-1)
      mse[𝑖] += σ²[𝑘]/f[𝑘]^2 * (1/C[𝑖,𝑘] + 1/ sum(C[1:I-𝑘,𝑘]))
    end
    mse[𝑖] *= C[𝑖,I]^2
  end
  # total mean square error
  mse_total = 0.0
  for 𝑖 ∈ 2:I
    tmp = 0.0
    for 𝑘 ∈ collect(I+1-𝑖 : I-1)
      tmp += (2σ²[𝑘] / f[𝑘]^2) / sum(C[1 : I-𝑘, 𝑘])
    end
    tmp *= C[𝑖,I] * sum(C[𝑖+1 : I, I])
    tmp += mse[𝑖]
    mse_total += tmp
  end
  return Mack(I, C, f, futclaims, R, sum(R), mse, mse_total)
end

"""
replacena!(df::DataFrame, replacement::Any)

Replaces all 'NA'-values with 'replacement'
"""
function replacena!(df::DataFrame, replacement::Any)
  nrows, ncols = size(df)
  for 𝑗 ∈ 1:ncols; for 𝑖 ∈ 1:nrows
    if ismissing(df[𝑖,𝑗]); df[𝑖,𝑗] = replacement; end
    end
  end
end

"""
upperleft(quadrat_mat::Array{Real,2})

Sets all values in the strict lower right triangle to 0.
"""
function upperleft(quadrat_mat::Matrix{T}) where T <: Real 
  c = deepcopy(quadrat_mat)
  n = size(c,1)
  if n ≠ size(c,2)
    error("known: Matrix is not quadratic")
  end
  for 𝑗 ∈ 1:n; for 𝑖 ∈ 1:n
    if 𝑖 > n-𝑗+1; c[𝑖,𝑗] = 0 end
  end; end
  c
end

"""
lowerright(quadrat_mat::Array{Real,2})

Sets all values in the upperleft triangle (incl. diagonal) to 0.
"""
function lowerright(quadrat_mat::Matrix{T}) where T <: Real 
  quadrat_mat - upperleft(quadrat_mat)
end

"""
    claims2cum(c::Matrix)

    Convert a quadratic upper left triangular reserve matrix to
    the correponding cumulative triangular reserve matrix.
    The (strictly) lower triangular part is set to 0
"""
function claims2cum(c::Matrix{T}) where T <: Real 
  if size(c,1) ≠ size(c,2)
    error("claims2cum: Matrix is not quadratic")
  end
 # cum = zeros(Real,size(c))
  #for i in 1:size(c, 1)
    #cum[i,1:size(c, 1)-i+1] = cumsum(c[i,1:size(c, 1)-i+1], 2)
  #end
  upperleft(cumsum(c, dims=2))
end

"""
    cum2claims(c::Matrix)

    Convert a quadratic cumulative claims matrix into a
    claims matrix.
"""
function cum2claims(c::Matrix{T}) where T <: Real 
  I = size(c,2)
  if size(c,1) ≠ I
    error("cum2claims: Matrix is not quadratic")
  end
  claims = zero(c)
  for 𝑖 ∈ 1:I
    claims[𝑖,1] = c[𝑖,1]
    for 𝑘 ∈ 2:I
     claims[𝑖,𝑘] = c[𝑖,𝑘]-c[𝑖,𝑘-1]
    end
  end
  claims
end

"""
    cum2futclaims(c::Array{Real,2})

    Calculate the future claims vector from the quadratic,
    (strictly) lower triangle of a projected cumulative claims
    matrix
"""
function cum2futureclaims(c::Matrix{T}) where T <: Real 
  I = size(c,2)
  if size(c,1) ≠ I
    error("cumclaims: Matrix is not quadratic")
  end
  claims = cum2claims(c)
  futureclaims = zero(claims[1:(I-1), 1])
   for 𝑘 ∈ 2:I
    futureclaims[𝑘-1] = sum(claims[I-𝑘+2:I,𝑘])
  end
  futureclaims
end
