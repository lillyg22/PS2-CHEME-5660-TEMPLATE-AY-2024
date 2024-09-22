# --- PRIVATE METHODS BELOW HERE -------------------------------------------------------------------------------------- #
_𝔼(X::Array{Float64,1}, p::Array{Float64,1}) = sum(X.*p)

"""
    𝕍(model::MyAdjacencyBasedTrinomialSharePriceTreeModel; level::Int = 0) -> Float64
"""
function _𝕍(model::MyAdjacencyBasedTrinomialSharePriceTreeModel; level::Int = 0)::Float64

    # initialize -
    variance_value = 0.0;
    X = Array{Float64,1}();
    p = Array{Float64,1}();

    # get the levels dictionary -
    levels = model.levels;
    nodes_on_this_level = levels[level]
    for i ∈ nodes_on_this_level
 
        # grab the node -
        node = model.data[i];
         
        # get the data -
        x_value = node.price;
        p_value = node.probability;
 
        # store unique data -
        if ((x_value ∈ X) == false && (p_value ∈ p) == false)
            push!(X,x_value);
            push!(p,p_value);
        end
    end

    # update -
    variance_value = (_𝔼(X.^2,p) - (_𝔼(X,p))^2)

    # return -
    return variance_value;
end

"""
    𝕍(model::MyAdjacencyBasedTrinomialSharePriceTreeModel, levels::Array{Int64,1}; startindex::Int64 = 0) -> Array{Float64,2}

Computes the variance of the model simulation. Takes a model::MyAdjacencyBasedTrinomialSharePriceTreeModel instance and a vector of
tree levels, i.e., time steps and returns a variance array where the first column is the time and the second column is the variance.
Each row is a time step.
"""
function _𝕍(model::MyAdjacencyBasedTrinomialSharePriceTreeModel, levels::Array{Int64,1}; startindex::Int64 = 0)::Array{Float64,2}

    # initialize -
    number_of_levels = length(levels);
    variance_value_array = Array{Float64,2}(undef, number_of_levels, 2);

    # loop -
    for i ∈ 0:(number_of_levels - 1)
        level = levels[i+1];
        variance_value = _𝕍(model, level=level);
        variance_value_array[i+1,1] = level + startindex
        variance_value_array[i+1,2] = variance_value;
    end

    # return -
    return variance_value_array;
end

function _𝔼(model::MyAdjacencyBasedTrinomialSharePriceTreeModel, levels::Array{Int64,1}; 
    startindex::Int64 = 0)::Array{Float64,2}

    # initialize -
    number_of_levels = length(levels);
    expected_value_array = Array{Float64,2}(undef, number_of_levels, 2);

    # loop -
    for i ∈ 0:(number_of_levels-1)

        # get the level -
        level = levels[i+1];

        # get the expected value -
        expected_value = _𝔼(model, level = level);
        expected_value_array[i+1,1] = level + startindex;
        expected_value_array[i+1,2] = expected_value;
    end

    # return -
    return expected_value_array;
end

function _𝔼(model::MyAdjacencyBasedTrinomialSharePriceTreeModel; level::Int = 0)::Float64

    # initialize -
    expected_value = 0.0;
    X = Array{Float64,1}();
    p = Array{Float64,1}();

    # get the levels dictionary -
    levels = model.levels;
    nodes_on_this_level = levels[level]
    for i ∈ nodes_on_this_level

        # grab the node -
        node = model.data[i];
        
        # get the data -
        x_value = node.price;
        p_value = node.probability;

        # store only the unique values -
        if ((x_value ∈ X) == false && (p_value ∈ p) == false)
            push!(X,x_value);
            push!(p,p_value);
        end
    end

    # compute -
    expected_value = _𝔼(X,p) # inner product

    # return -
    return expected_value
end

# --- PRIVATE METHODS ABOVE HERE -------------------------------------------------------------------------------------- #

# --- PUBLIC METHODS BELOW HERE --------------------------------------------------------------------------------------- #

# Short-hand notation for the variance computation -
Var(model::MyAdjacencyBasedTrinomialSharePriceTreeModel, levels::Array{Int64,1}; 
    startindex::Int64 = 0) = _𝕍(model, levels, startindex = startindex)

𝔼(model::MyAdjacencyBasedTrinomialSharePriceTreeModel, levels::Array{Int64,1}; 
    startindex::Int64 = 0) = _𝔼(model, levels, startindex = startindex)

# Short-hand notation for lattice parameter computation -
# TODO: you need to implement the _analyze_real_world_ternary_single_asset function in Studdent.jl
(m::RealWorldTrinomialProbabilityMeasure)(R::Array{Float64,1};  
    Δt::Float64 = (1.0/252.0), ϵ::Float64 = 0.0)::MyRealWorldTrinomialSharePriceTreeParameters = _analyze_real_world_ternary_single_asset(R, Δt, ϵ)
# --- PUBLIC METHODS ABOVE HERE --------------------------------------------------------------------------------------- #