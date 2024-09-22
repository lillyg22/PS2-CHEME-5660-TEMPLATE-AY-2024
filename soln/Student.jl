"""
    function _analyze_real_world_ternary_single_asset(R::Array{Float64,1}, 
                Δt::Float64, ϵ::Float64) -> MyRealWorldTrinomialSharePriceTreeParameters

The `_analyze_real_world_ternary_single_asset` function computes the real-world trinomial share price tree parameters for a single asset. 
The function takes as input the log growth rate for the asset, the time increment, and the margin around zero. 
The function returns a `MyRealWorldTrinomialSharePriceTreeParameters` object that contains the computed parameters.


### Arguments
- `R::Array{Float64,1}`: An array of log growth rates for the asset.
- `Δt::Float64`: The time increment used to compute the log growth. The default value is `1/252`, i.e., one trading day in units of years.
- `ϵ::Float64`: The margin around zero used to determine the up, down, and unchanged moves.

### Returns
- `MyRealWorldTrinomialSharePriceTreeParameters`: A `MyRealWorldTrinomialSharePriceTreeParameters` object that contains the computed parameters.
"""
function _analyze_real_world_ternary_single_asset(R::Array{Float64,1}, 
    Δt::Float64, ϵ::Float64)::MyRealWorldTrinomialSharePriceTreeParameters

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

"""
    function populate(model::MyAdjacencyBasedTrinomialSharePriceTreeModel; 
        Sₒ::Float64 = 100.0, h::Int = 1) -> MyAdjacencyBasedTrinomialSharePriceTreeModel

The `populate` function builds the trinomial share price tree for a given model. 
The function takes as input the model, the initial price of the asset and the number of time steps (levels) in the tree.
Each node in the tree is represented by a `MyTrinomialLatticeNodeModel` object that holds the price, probability, and path information.
This reference implementation builds a full trinomial tree with the given number of levels. Thus, there are duplicates of the same price and probability information in the tree.


### Arguments
- `model::MyAdjacencyBasedTrinomialSharePriceTreeModel`: The model for which we want to build the trinomial share price tree.
- `Sₒ::Float64`: The initial price of the asset. The default value is `100.0`.
- `h::Int`: The number of time steps (levels) in the tree. The default value is `1`.

### Returns
- `MyAdjacencyBasedTrinomialSharePriceTreeModel`: The model with the trinomial share price tree populated.
"""
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