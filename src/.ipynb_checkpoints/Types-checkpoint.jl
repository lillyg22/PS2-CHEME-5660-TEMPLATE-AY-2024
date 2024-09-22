abstract type AbstractPriceTreeModel end
abstract type AbstractProbabilityMeasure end

mutable struct MyTrinomialLatticeNodeModel

    # data -
    price::Float64
    probability::Float64
    path::Array{Char,1}

    # constructor -
    MyTrinomialLatticeNodeModel() = new();
end

mutable struct MyAdjacencyBasedTrinomialSharePriceTreeModel <: AbstractPriceTreeModel

    # Parameters set by the user
    p::Float64;            # what is the probability of an up move
    p̄::Float64;            # what is the probability of a down move
    q::Float64;            # what is the probability of an unch. move
    u::Float64;            # what is the up move factor
    d::Float64;            # what is the down move factor
    ϵ::Float64;            # what is the ϵ-margin around zero
    Δt::Float64;           # what is the time step for the tree

    # Properties of the tree computed by the populate method
    data::Dict{Int64, MyTrinomialLatticeNodeModel}
    connectivity::Dict{Int64,Array{Int64,1}}
    levels::Union{Nothing, Dict{Int64,Array{Int64,1}}}

    # constructor 
    MyAdjacencyBasedTrinomialSharePriceTreeModel() = new()
end

"""
   struct RealWorldTrinomialProbabilityMeasure <: AbstractProbabilityMeasure

Immutable type that represents the real-world probability measure. 
This type is passed as an argument to various functions to indicate that the real-world probability measure should be used in calculations.   
"""
struct RealWorldTrinomialProbabilityMeasure <: AbstractProbabilityMeasure
    RealWorldTrinomialProbabilityMeasure() = new()
end

mutable struct MyRealWorldTrinomialSharePriceTreeParameters

    # data -
    p::Float64;            # what is the probability of an up move
    p̄::Float64;            # what is the probability of a down move
    q::Float64;            # what is the probability of an unch. move
    u::Float64;            # what is the up move factor
    d::Float64;            # what is the down move factor
    ϵ::Float64;            # what is the ϵ-margin around zero
    Δt::Float64;           # what is the time step for the tree

    # constructor -
    MyRealWorldTrinomialSharePriceTreeParameters() = new(); # empty constructor
end