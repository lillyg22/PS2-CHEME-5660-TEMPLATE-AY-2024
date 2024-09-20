# --- PRIVATE METHODS BELOW HERE -------------------------------------------------------------------------------------- #
function _jld2(path::String)::Dict{String,Any}
    return load(path);
end
# --- PRIVATE METHODS ABOVE HERE -------------------------------------------------------------------------------------- #


# --- PUBLIC METHODS BELOW HERE --------------------------------------------------------------------------------------- #
MyTrainingMarketPriceDataSet() = _jld2(joinpath(_PATH_TO_DATA, "SP500-Daily-OHLC-1-3-2018-to-12-29-2023.jld2"));
# --- PUBLIC METHODS ABOVE HERE --------------------------------------------------------------------------------------- #