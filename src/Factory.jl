"""
    build(type::Type{AdjacencyBasedTernaryCommodityPriceTree};
        h::Int64 = 1, price::Float64 = 1.0, u::Float64 = 0.02, d::Float64 = 0.01) -> AdjacencyBasedTernaryCommodityPriceTree
"""
function build(modeltype::Type{MyAdjacencyBasedTernarySharePriceTree}, 
    data::NamedTuple)::MyAdjacencyBasedTernarySharePriceTree

    # get data required to build the tree from the NamedTuple -
    h = data.h;             # how many levels in the tree
    price = data.price;     # what is the price at the root of the tree
    p = data.p;             # what is the probability of an up move
    p̄ = data.p̄;             # what is the probability of a down move
    u = data.u;             # what is the up move factor
    d = data.d;             # what is the down move factor
    ϵ = data.epsilon        # what is the ϵ-margin around zero
    Δt = data.Δt            # what is the time step for the tree

    # compute the probability of unch. -> save this as q
    q = 1 - p - p̄;

    # initialize -
    model = modeltype(); # build an empty tree model
    Nₕ = sum([3^i for i ∈ 0:h]) # compute how many nodes do we have in the tree
    P = Dict{Int64, Float64}() # use Dict for zero-based array hack. Hold price information
    connectivity = Dict{Int64, Array{Int64,1}}() # holds tree connectivity information

    # setup Δ - the amount the price moves up, unch, or down
    Δ = [u, 0, -d];

    
    P[0] = price; # set the price at the root of the tree

    # build connectivity -
    for i ∈ 0:(Nₕ - 3^h - 1)
        
        # what is the *parent* price
        Pᵢ = P[i]

        # Compute the children for this node -
        Cᵢ = [j for j ∈ (3*i+1):(3*i+3)]; 
        connectivity[i] = Cᵢ # stores the children indices of node i

        # cmpute the prices at the child nodes
        for c ∈ 1:3 # for each node (no matter what i) we have three children

            # what is the child index?
            child_index = Cᵢ[c]

            # compute the new price for the child node
            P[child_index] = Pᵢ*exp(Δ[c]*Δt);
        end
    end

    # compute the levels -


    # set the data, and connectivity for the model -
    model.data = P;
    model.connectivity = connectivity;


    # return -
    return model;
end