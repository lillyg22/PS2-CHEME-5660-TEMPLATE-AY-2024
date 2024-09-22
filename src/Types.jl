abstract type AbstractPriceTreeModel end
abstract type AbstractProbabilityMeasure end

"""
    mutable struct MyTrinomialLatticeNodeModel

The `MyTrinomialLatticeNodeModel` type is a mutable struct that holds the trinomial lattice node.

### Fields
- `price::Float64`: The price at the node.
- `probability::Float64`: The probability of the node.
- `path::Array{Char,1}`: The path to the node from the root.

We use the code:
- `a::Char` to represent an up move.
- `b::Char` to represent an unchanged move.
- `c::Char` to represent a down move.
"""
mutable struct MyTrinomialLatticeNodeModel

    # data -
    price::Float64
    probability::Float64
    path::Array{Char,1}

    # constructor -
    MyTrinomialLatticeNodeModel() = new();
end

"""
    mutable struct MyAdjacencyBasedTrinomialSharePriceTreeModel <: AbstractPriceTreeModel

The `MyAdjacencyBasedTrinomialSharePriceTreeModel` type is a mutable struct that holds the trinomial share price tree.

### Fields
- `p::Float64`: The probability of an up move.
- `p̄::Float64`: The probability of a down move.
- `q::Float64`: The probability of an unchanged move.
- `u::Float64`: The up move factor.
- `d::Float64`: The down move factor.
- `ϵ::Float64`: The ϵ-margin around zero.
- `Δt::Float64`: The time step for the tree.

### Computed Fields
- `data::Dict{Int64, MyTrinomialLatticeNodeModel}`: The data for the tree (price, probability, path).
- `connectivity::Dict{Int64,Array{Int64,1}}`: The connectivity of the tree.
- `levels::Union{Nothing, Dict{Int64,Array{Int64,1}}}`: The levels of the tree.
"""
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

"""
    mutable struct MyRealWorldTrinomialSharePriceTreeParameters

The `MyRealWorldTrinomialSharePriceTreeParameters` type is a mutable struct that holds the parameters required to build a trinomial share price tree.

### Fields
- `p::Float64`: The probability of an up move.
- `p̄::Float64`: The probability of a down move.
- `q::Float64`: The probability of an unchanged move.
- `u::Float64`: The up move factor.
- `d::Float64`: The down move factor.
- `ϵ::Float64`: The ϵ-margin around zero.
- `Δt::Float64`: The time step for the tree.

"""
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