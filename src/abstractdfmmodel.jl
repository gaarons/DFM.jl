abstract AbstractModel{T}

# I dont want to define my own Parameter types at all.  My parameters are not as complex as DSGE, simple steady-state or time-invariant only
description(m::AbstractModel) = "DFM Model $(spec(m)) (subspec $(subspec(m)))"

function Base.show{T<:AbstractModel}(io::IO, m::T)
    @printf io "Dynamic Factor Model (Style of BGR 2008)\n"
    @printf io "%s\n" T
    @printf io "description:\n %s\n"          description(m)
    @printf io "model frequency:           %s\n" m_frequency(m)
    @printf io "data vintage:           %s\n" data_vintage(m)
    @printf io "number of baskets:           %s\n" n_baskets(m)
    @printf io "number of series:           %s\n" n_series(m)


# Loop through the basket names and print them out with quantities of fac extracted and number series in basket
    # @printf io "factors per basket (vector):           %s\n" n_baskets_fac(m)
    # @printf io "series per basket (vector):           %s\n" n_baskets_series(m)
# Write these looping get_setting to construct a string below....?
    # How many factors drawn from each basket, will be a vector of length equal to :num_baskets
    # n_baskets_fac(m::AbstractModel)                = get_setting(m, :num_baskets_fac)

end

@inline function Base.getindex(m::AbstractModel, i::Integer)
    if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]
    end
end

# need to define like this so we can disable bounds checking
@inline function Base.getindex(m::AbstractModel, k::Symbol)
    i = m.keys[k]
    @inbounds if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]
    end
end


@inline function Base.setindex!{T<:Number}(m::AbstractModel, value::T, i::Integer)
    if i <= (j = length(m.parameters))
        param = m.parameters[i]
        param.value = value
        if isa(param, ScaledParameter)
            param.scaledvalue = param.scaling(value)
        end
        return param
    else
        steady_state_param = m.steady_state[i-j]
        steady_state_param.value = value
        return steady_state_param
    end
end

Base.setindex!(m::AbstractModel, value, k::Symbol) = Base.setindex!(m, value, m.keys[k])

# Index into data matrix of first period forecasted without full data
function fcst_start_index(m::AbstractModel)
    fcst_start_quarter = get_setting(m, :date_fcstregime_start)
    presample_start_quarter = get_setting(m, :date_presample_start)
    quarters = subtract_quarters(fcst_start_quarter, presample_start_quarter)
    index = quarters+1
    return index
end

# Number of presample periods
function n_presample_periods(m::AbstractModel)
    mainsample_start_quarter = get_setting(m, :date_mainsample_start)
    presample_start_quarter = get_setting(m, :date_presample_start)
    quarters = subtract_quarters(mainsample_start_quarter, presample_start_quarter)
    return quarters
end

# Number of a few things that are useful
n_observables(m::AbstractModel)            = length(m.observables)
# n_parameters(m::AbstractModel)             = length(m.parameters)
n_baskets(m::AbstractModel)                = get_setting(m, :num_baskets)
m_frequency(m::AbstractModel)                = get_setting(m, :frequency)

# Interface for I/O settings
spec(m::AbstractModel)         = m.spec
subspec(m::AbstractModel)      = m.subspec
saveroot(m::AbstractModel)     = get_setting(m, :saveroot)
dataroot(m::AbstractModel)     = get_setting(m, :dataroot)

# Interface for data
data_vintage(m::AbstractModel) = get_setting(m, :data_vintage)
# cond_vintage(m::AbstractModel) = get_setting(m, :cond_vintage)

# Interface for general computation settings
use_parallel_workers(m::AbstractModel)    = get_setting(m, :use_parallel_workers)

#=
Build paths to where input/output/results data are stored.

Description:
Creates the proper directory structure for input and output files, treating the DSGE/save
    directory as the root of a savepath directory subtree. Specifically, the following
    structure is implemented:

    dataroot/

    savepathroot/
                 output_data/<spec>/<subspec>/log/
                 output_data/<spec>/<subspec>/<out_type>/raw/
                 output_data/<spec>/<subspec>/<out_type>/work/
                 output_data/<spec>/<subspec>/<out_type>/tables/
                 output_data/<spec>/<subspec>/<out_type>/figures/

Note: we refer to the savepathroot/output_data/<spec>/<subspec>/ directory as saveroot.
=#
"""
```
logpath(model)
```
Returns path to log file. Path built as
```
<output root>/output_data/<spec>/<subspec>/log/log_<filestring>.log
```
"""
function logpath(m::AbstractModel)
    return savepath(m, "log", "log.log")
end

strs = [:work, :raw, :tables, :figures]
fns = [symbol(x, "path") for x in strs]
for (str, fn) in zip(strs, fns)
    @eval begin
        # First eval function
        function $fn{T<:AbstractString}(m::AbstractModel,
                                             out_type::T,
                                             file_name::T="",
                                             filestring_addl::Vector{T}=Vector{T}())
            return savepath(m, out_type, $(string(str)), file_name, filestring_addl)
        end

        # Then, add docstring to it
        @doc $(
        """
        ```
        $fn{T<:AbstractString}(m::AbstractModel, out_type::T, file_name::T="")
        ```

        Returns path to specific $str output file, creating containing directory as needed. If
        `file_name` not specified, creates and returns path to containing directory only. Path built
        as
        ```
        <output root>/output_data/<spec>/<subspec>/<out_type>/$str/<file_name>_<filestring>.<ext>
        ```
        """
        ) $fn
    end
end

# Not exposed to user. Actually create path and insert model string to file name.
function savepath{T<:AbstractString}(m::AbstractModel,
                                     out_type::T,
                                     sub_type::T,
                                     file_name::T="",
                                     filestring_addl::Vector{T}=Vector{T}())
    # Containing dir
    path = joinpath(saveroot(m), "output_data", spec(m), subspec(m), out_type, sub_type)
    if !isdir(path)
        mkpath(path)
    end

    # File with model string inserted
    if !isempty(file_name)
        if isempty(filestring_addl)
            myfilestring = filestring(m)
        else
            myfilestring = filestring(m, filestring_addl)
        end
        (base, ext) = splitext(file_name)
        file_name_detail = base * myfilestring * ext
        path = joinpath(path, file_name_detail)
    end

    return path
end

# Input data handled slightly differently, because it is not model-specific.
"""
```
inpath{T<:AbstractString}(m::AbstractModel, in_type::T, file_name::T="")
```

Returns path to specific input data file, creating containing directory as needed. If
`file_name` not specified, creates and returns path to containing directory only. Valid
`in_type` includes:

* `"data"`: recorded data
* `"cond"`: conditional data - nowcasts for the current forecast quarter, or related
* `"user"`: user-supplied data for starting parameter vector, hessian, or related

Path built as
```
<data root>/<in_type>/<file_name>
```
"""
function inpath{T<:AbstractString}(m::AbstractModel, in_type::T, file_name::T="")
    path = dataroot(m)
    # Normal cases.
    if in_type == "data" || in_type == "cond"
        path = joinpath(path, in_type)
    # User-provided inputs. May treat this differently in the future.
    elseif in_type == "user"
        path = joinpath(path, "user")
    else
        error("Invalid in_type: ", in_type)
    end

    # Containing dir
    if !isdir(path)
        mkpath(path)
    end

    # If file_name provided, return full path
    if !isempty(file_name)
        path = joinpath(path, file_name)
    end

    return path
end

filestring(m::AbstractModel) = filestring(m, Vector{AbstractString}())
filestring(m::AbstractModel, d::AbstractString) = filestring(m, [d])
function filestring{T<:AbstractString}(m::AbstractModel,
                                        d::Vector{T})
    if !m.testing
        filestrings = Vector{T}()
        for (skey, sval) in m.settings
            if sval.print
                push!(filestrings, to_filestring(sval))
            end
        end
        append!(filestrings, d)
        sort!(filestrings)
        return "_"*join(filestrings, "_")
    else
        return "_test"
    end
end
