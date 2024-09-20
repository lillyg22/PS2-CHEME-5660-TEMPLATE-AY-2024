abstract type AbstractPriceTreeModel end
abstract type AbstractProbabilityMeasure end



mutable struct MyAdjacencyBasedTernarySharePriceTree <: AbstractPriceTreeModel

    # data -


    # computed properties -
    data::Dict{Int64,Float64}
    connectivity::Dict{Int64,Array{Int64,1}}
    levels::Union{Nothing, Dict{Int64,Array{Int64,1}}}

    # constructor 
    MyAdjacencyBasedTernarySharePriceTree() = new()
end

"""
   struct RealWorldTrinomialProbabilityMeasure <: AbstractProbabilityMeasure

Immutable type that represents the real-world probability measure. 
This type is passed as an argument to various functions to indicate that the real-world probability measure should be used in calculations.   
"""
struct RealWorldTernaryProbabilityMeasure <: AbstractProbabilityMeasure
    RealWorldTernaryProbabilityMeasure() = new()
end