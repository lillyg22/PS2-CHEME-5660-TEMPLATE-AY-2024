# --- PRIVATE METHODS BELOW HERE -------------------------------------------------------------------------------------- #
function _build_nodes_level_dictionary(levels::Int64)::Dict{Int64,Array{Int64,1}}

    # initialize -
    index_dict = Dict{Int64, Array{Int64,1}}()

    counter = 0
    for l = 0:levels
        
        # create index set for this level -
        index_array = Array{Int64,1}()
        for _ = 1:(3^l)
            counter = counter + 1
            push!(index_array, counter)
        end

        index_dict[l] = (index_array .- 1) # zero based
    end

    # return -
    return index_dict
end
# --- PRIVATE METHODS ABOVE HERE -------------------------------------------------------------------------------------- #

# --- PUBLIC METHODS BELOW HERE --------------------------------------------------------------------------------------- #
"""
    function build(modeltype::Type{MyAdjacencyBasedTernarySharePriceTree}, 
        parameters::MyRealWorldTernarySharePriceTreeParameters) -> MyAdjacencyBasedTernarySharePriceTree

This function builds and initializes a MyAdjacencyBasedTernarySharePriceTree object from a MyRealWorldTernarySharePriceTreeParameters instance.

### Arguments
- `modeltype::Type{MyAdjacencyBasedTernarySharePriceTree}`: The type of the model to build, i.e. a MyAdjacencyBasedTernarySharePriceTree instance
- `parameters::MyRealWorldTernarySharePriceTreeParameters`: The parameters to use to build the model. 

### Returns
- `MyAdjacencyBasedTernarySharePriceTree`: The initialized MyAdjacencyBasedTernarySharePriceTree instance. 
"""
function build(modeltype::Type{MyAdjacencyBasedTrinomialSharePriceTreeModel}, 
    parameters::MyRealWorldTrinomialSharePriceTreeParameters)::MyAdjacencyBasedTrinomialSharePriceTreeModel

    # initialize -
    tree = modeltype();

    # set the constant parameters on the tree -
    tree.p = parameters.p;
    tree.p̄ = parameters.p̄;
    tree.q = parameters.q;
    tree.u = parameters.u;
    tree.d = parameters.d;
    tree.ϵ = parameters.ϵ;
    tree.Δt = parameters.Δt;

    # return -
    return tree;
end
# --- PUBLIC METHODS ABOVE HERE --------------------------------------------------------------------------------------- #