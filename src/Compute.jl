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


function _analyze_real_world_ternary_single_asset(R::Array{Float64,1}, Δt::Float64, ϵ::Float64)::MyRealWorldTrinomialSharePriceTreeParameters

    # initialize -
    N = length(R);
    p = 0.0; # probability of up move
    p̄ = 0.0; # probability of down move
    q = 0.0; # probability of unch. move
    u = Array{Float64,1}(); # up move factor array
    d = Array{Float64,1}(); # down move factor array

    # get indicies for the up and down moves -
    i₊ = findall(x -> x > ϵ, R);
    i₋ = findall(x -> x < -ϵ, R);
    iₒ = findall(x -> abs(x) <= ϵ, R);

    # compute the probabilities -
    p = length(i₊)/N; # up
    p̄ = length(i₋)/N; # down
    q = length(iₒ)/N; # unch.

    # compute the up and down move factors -
    for i ∈ i₊
        push!(u, exp(R[i]*Δt)); # up move factors
    end

    for i ∈ i₋
        push!(d, exp(R[i]*Δt)); # down move factors
    end

    # build parameter wrapper -
    parameters = MyRealWorldTrinomialSharePriceTreeParameters();
    parameters.p = p;
    parameters.p̄ = p̄;
    parameters.q = q;
    parameters.u = mean(u);
    parameters.d = mean(d);
    parameters.ϵ = ϵ;
    parameters.Δt = Δt;

    # return -
    return parameters;
end

# --- PRIVATE METHODS ABOVE HERE -------------------------------------------------------------------------------------- #

# --- PUBLIC METHODS BELOW HERE --------------------------------------------------------------------------------------- #

# Short-hand notation for the variance computation -
Var(model::MyAdjacencyBasedTrinomialSharePriceTreeModel, levels::Array{Int64,1}; 
    startindex::Int64 = 0) = _𝕍(model, levels, startindex = startindex)

𝔼(model::MyAdjacencyBasedTrinomialSharePriceTreeModel, levels::Array{Int64,1}; 
    startindex::Int64 = 0) = _𝔼(model, levels, startindex = startindex)

# Short-hand notation for lattice parameter computation -
(m::RealWorldTrinomialProbabilityMeasure)(R::Array{Float64,1};  
    Δt::Float64 = (1.0/252.0), ϵ::Float64 = 0.0)::MyRealWorldTrinomialSharePriceTreeParameters = _analyze_real_world_ternary_single_asset(R, Δt, ϵ)


function populate(model::MyAdjacencyBasedTrinomialSharePriceTreeModel; 
    Sₒ::Float64 = 100.0, h::Int = 1)::MyAdjacencyBasedTrinomialSharePriceTreeModel

    # get data required to build the tree from the NamedTuple -
    p = model.p;             # what is the probability of an up move
    p̄ = model.p̄;             # what is the probability of a down move
    q = model.q;             # what is the probability of an unch. move
    u = model.u;             # what is the up move factor
    d = model.d;             # what is the down move factor
    ϵ = model.ϵ              # what is the ϵ-margin around zero
    Δt = model.Δt            # what is the time step for the tree

    # initialize -
    Nₕ = sum([3^i for i ∈ 0:h]) # compute how many nodes do we have in the tree
    P = Dict{Int64, MyTrinomialLatticeNodeModel}() # use Dict for zero-based array hack. Hold price and probability information
    connectivity = Dict{Int64, Array{Int64,1}}() # holds tree connectivity information

    # setup Δ - the amount the price moves up, unch, or down
    Δ = [u, 0, -d];

    rootnode = MyTrinomialLatticeNodeModel();
    rootnode.price = Sₒ;
    rootnode.probability = 1.0;
    rootnode.path = [];
    P[0] = rootnode; # set the price at the root of the tree

    # build connectivity -
    for i ∈ 0:(Nₕ - 3^h - 1)
        
        # what is the *parent* price
        Pᵢ = P[i].price;
        path_parent = P[i].path;

        # Compute the children for this node -
        Cᵢ = [j for j ∈ (3*i+1):(3*i+3)]; 
        connectivity[i] = Cᵢ # stores the children indices of node i

        # cmpute the prices at the child nodes
        for c ∈ 1:3 # for each node (no matter what i) we have three children

            my_path = copy(path_parent);
            if (c == 1)
                push!(my_path, 'a');
            elseif (c == 2)
                push!(my_path, 'b');
            else (c == 3)
                push!(my_path, 'c');
            end
        
            # build empty node model -
            node = MyTrinomialLatticeNodeModel();

            # what is the child index?
            child_index = Cᵢ[c]

            # to build a tree node, we need two things: the price and the probability
            node.price = Pᵢ*exp(Δ[c]*Δt);
            node.probability = 0.1; # placeholder
            node.path = my_path;

            # compute the new price for the child node
            P[child_index] = node;
        end
    end

    # for each node, compute the probability of seeing that node -
    for (_ , v) ∈ P

        patharray = v.path;
        ntotal = length(patharray);

        # compute the n up moves, n down moves, and n unchanged moves
        n_up = count(x -> x == 'a', patharray);
        n_unch = count(x -> x == 'b', patharray);
        n_down = count(x -> x == 'c', patharray);

        # compute the probability of seeing this path -
        factor = (factorial(ntotal))/(factorial(n_up)*factorial(n_unch)*factorial(n_down));
        v.probability = factor*(p^n_up)*(q^n_unch)*(p̄^n_down);
    end

    # set the data, and connectivity for the model -
    model.data = P;
    model.connectivity = connectivity;
    model.levels = _build_nodes_level_dictionary(h);

    # return -
    return model;
end


"""
    function log_growth_matrix(dataset::Dict{String, DataFrame}, 
                firms::Array{String,1}; Δt::Float64 = (1.0/252.0), risk_free_rate::Float64 = 0.0) -> Array{Float64,2}

The `log_growth_matrix` function computes the excess log growth matrix for a given set of firms where we define the log growth as:

```math
    \\mu_{t,t-1}(r_{f}) = \\frac{1}{\\Delta t} \\log\\left(\\frac{S_{t}}{S_{t-1}}\\right) - r_f
```

where ``S_t`` is the volume weighted average price (units: USD/share) at time `t`, ``\\Delta t`` is the time increment (in years), and ``r_f`` is the annual risk-free rate (units: 1/years) assuming
continuous compounding.

### Arguments
- `dataset::Dict{String, DataFrame}`: A dictionary of data frames where the keys are the firm ticker symbols and the values are the data frames holding price data. We use the `volume_weighted_average_price` column to compute the log growth by default.
- `firms::Array{String,1}`: An array of firm ticker symbols for which we want to compute the log growth matrix.
- `Δt::Float64`: The time increment used to compute the log growth. The default value is `1/252`, i.e., one trading day in units of years.
- `risk_free_rate::Float64`: The risk-free rate used to compute the log growth. The default value is `0.0`.
- `keycol::Symbol`: The column in the data frame to use to compute the log growth. The default value is `:volume_weighted_average_price`.
- `testfirm::String`: The firm ticker symbol to use to determine the number of trading days. By default, we use "AAPL".

### Returns
- `Array{Float64,2}`: An array of the excess log growth values for the given set of firms. The time series is the rows and the firms are the columns. The columns are ordered according to the order of the `firms` array.

### See:
* The `DataFrame` type (and methods for working with data frames) is exported from the [DataFrames.jl package](https://dataframes.juliadata.org/stable/)
"""
function log_growth_matrix(dataset::Dict{String, DataFrame}, 
    firms::Array{String,1}; Δt::Float64 = (1.0/252.0), risk_free_rate::Float64 = 0.0, 
    testfirm="AAPL", keycol::Symbol = :volume_weighted_average_price)::Array{Float64,2}

    # initialize -
    number_of_firms = length(firms);
    number_of_trading_days = nrow(dataset[testfirm]);
    return_matrix = Array{Float64,2}(undef, number_of_trading_days-1, number_of_firms);

    # main loop -
    for i ∈ eachindex(firms) 

        # get the firm data -
        firm_index = firms[i];
        firm_data = dataset[firm_index];

        # compute the log returns -
        for j ∈ 2:number_of_trading_days
            S₁ = firm_data[j-1, keycol];
            S₂ = firm_data[j, keycol];
            return_matrix[j-1, i] = (1/Δt)*(log(S₂/S₁)) - risk_free_rate;
        end
    end

    # return -
    return return_matrix;
end
# --- PUBLIC METHODS ABOVE HERE --------------------------------------------------------------------------------------- #