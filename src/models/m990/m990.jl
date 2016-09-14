"""
```
Model990{T} <: AbstractModel{T}
```

The `Model990` type defines the structure of the FRBNY DSGE model.

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed as a function of elements of
  `parameters`.

* `keys::Dict{Symbol,Int}`: Maps human-readable names for all model parameters and
  steady-states to their indices in `parameters` and `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readible names to row and column
indices in the matrix representations of of the measurement equation and equilibrium
conditions.

* `endogenous_states::Dict{Symbol,Int}`: Maps each state to a column in the measurement and
  equilibrium condition matrices.

* `exogenous_shocks::Dict{Symbol,Int}`: Maps each shock to a column in the measurement and
  equilibrium condition matrices.

* `expected_shocks::Dict{Symbol,Int}`: Maps each expected shock to a column in the
  measurement and equilibrium condition matrices.

* `equilibrium_conditions::Dict{Symbol,Int}`: Maps each equlibrium condition to a row in the
  model's equilibrium condition matrices.

* `endogenous_states_augmented::Dict{Symbol,Int}`: Maps lagged states to their columns in
  the measurement and equilibrium condition equations. These are added after Gensys solves the
  model.

* `observables::Dict{Symbol,Int}`: Maps each observable to a row in the model's measurement
  equation matrices.

#### Model Specifications and Settings

* `spec::AbstractString`: The model specification identifier, "m990", cached here for
  filepath computation.

* `subspec::AbstractString`: The model subspecification number, indicating that some
  parameters from the original model spec ("ss0") are initialized differently. Cached here for
  filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation without changing
  the economic or mathematical setup of the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. Can be is seeded to ensure
  reproducibility in algorithms that involve randomness (such as Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If `true`, settings from
  `m.test_settings` are used in place of those in `m.settings`.

* `data_series::Dict{Symbol,Vector{Symbol}}`: A dictionary that
  stores data sources (keys) and lists of series mnemonics
  (values). DSGE.jl will fetch data from the Federal Reserve Bank of
  St. Louis's FRED database; all other data must be downloaded by the
  user. See `load_data` for further details.
"""

# Overwrite anywhere we see Model990 with mname = split(basename(@__FILE__),'.')[1]
# This is of course, if we are allowed to dynamically name a concrete type <:AbstractModel .....????
type Model990{T} <: AbstractModel{T}
    # parameters::ParameterVector{T}                  # vector of all time-invariant model parameters
    # steady_state::ParameterVector{T}                # model steady-state values
    # keys::Dict{Symbol,Int}                          # human-readable names for all the model
                                                    # parameters and steady-n_states
    observables::Dict{Symbol,Int}                   #
    transf::Dict{Symbol,AbstractString}                   #
    nfrequency::Dict{Symbol,AbstractString}                   #
    startDate::Dict{Symbol,AbstractString}                   #
    endDate::Dict{Symbol,AbstractString}                   #
    blockNames::Array{AbstractString,2}                   #
    # block_1::Dict{Symbol,Int}                  #
    # block_2::Dict{Symbol,Int}                  #
    # block_3::Dict{Symbol,Int}                  #
    # block_4::Dict{Symbol,Int}                  #
    # block_5::Dict{Symbol,Int}                  #
    dataBasket::Dict{Symbol, Array{Int,2}}

    p         # number of lags in the SKF and FIS kalman estimation
    r         # Number of factors in each basket


    spec::ASCIIString                               # Model specification number (eg "m990")
    subspec::ASCIIString                            # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                  # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}             # Settings/flags for testing mode
    rng::MersenneTwister                            # Random number generator
    testing::Bool                                   # Whether we are in testing mode or not

    data_series::Dict{Symbol,Vector{Symbol}}       # Keys = data sources, values = vector of series mnemonics
    data_transforms::OrderedDict{Symbol,Function}  # functions to transform raw data into input matrix
end


# initialize empty Res if one not on disk yet


description(m::Model990) = "DFM Model m990, $(m.subspec)"
# description(m::Model990) = "DFM Model $(split(basename(@__FILE__),'.')[1]), $(m.subspec)"

# default value of "ss2" there is a way to call the model with a different subspec from the command line
# m = Model990("ss1")

# Will loop through vintages in this way
# m <= Setting(:data_vintage, "yymmdd")
function Model990(subspec::AbstractString="ss2")

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister()
    testing            = false

    # Raw Pull everything and anything that we might be interested in
    # Only the variables that appear in observables are passed to later parts of program
    # Set up data sources and series
    fred_series        = [:GDP, :GDPCTPI, :PCE, :FPI, :CNP16OV, :CE16OV, :PRS85006013,
                          :UNRATE, :AWHNONAG, :DFF, :BAA, :GS10, :PRS85006063, :CES0500000030, :CLF16OV,
                          :PCEPILFE, :COMPNFB, :THREEFYTP10]
    spf_series         = [:ASACX10]
    fernald_series     = [:TFPJQ, :TFPKQ]
    longrate_series    = [:FYCZZA]
    # ois data taken care of in load_data
    # Eventually want to load in BlsData and Google trends information and geographic information (regional)
    # Boils down to writing future packages for those datasources

    data_series = Dict{Symbol,Vector{Symbol}}(:fred => fred_series, :spf => spf_series,
                                              :fernald => fernald_series, :longrate => longrate_series)

    # set up data transformations
    # data_transforms = Dict{Symbol,Function}()
    data_transforms = Dict{Symbol,Function}()

    # initialize empty model
    m = Model990{Float64}(
      # # Anything that has an <m.WHATEVER> declaration should be initialized to types here!
      # # Come back to this section to ammend the initialized types
      #       # model dictionaries for 1. transformations, 2. native frequency, 3. StartDate, 4. EndDate, 5.  and basket string definition
      #       # the following 5 Dict{Symbol,Int}() allow the programmer to define up to 5 different types of factor baskets and ensure that the variables basketed are never lost or confused!
      #       # Dictionary of factor strings to block numbers
      #

        # Ordering here does not matter, so long as you have initialized proper quantities
            # observables, transf, nfrequency, startdates, enddates,
            Dict{Symbol,Int}(), Dict{Symbol,AbstractString}(), Dict{Symbol,AbstractString}(), Dict{Symbol,AbstractString}(), Dict{Symbol,AbstractString}(),

            # # Blocknames,
            Array{AbstractString,2}(),
            # up to 5 blocks - for now going to run with the dataBasket and work on basketing from there
            #Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(),

            # # full databasket binary
            Dict{Symbol, Array{Int,2}}(),

            # parameters for extraction (lags in kalman, and number of factors per basket)
            Int64[], Int64[],

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            data_series,
            data_transforms)

            # settings_m990!(m)


    # Set settings
    default_settings!(m)  # removed variables in here
    # default_test_settings!(m) # havent touched this --> will need to later

    # Set data to collect and transformations to apply
    init_subspec!(m)  # In subspec need to store the number of factor types, the count used in each basket, the allocation decisions.

    init_transformation!(m)

# Do I even want factor variables broken down like this?  I thought I just needed a matrix of 0 & 1's
    # factor_baskets = Dict{Symbol,Vector{Symbol}}(:f_global => observables, :f_macro => obs_macro,
    #                                              :f_finance => obs_financial)

    return m
end # function

"""
`init_model_indices!(m::Model990)`

Arguments:
`m:: Model990`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_transformation!(m::Model990)
  # Standard data transformations - Leave out 0 [user-defined] and 1 [levels]
  default_transforms = Dict{AbstractString, Function}("2"=>diff,"4"=>log,"5"=>difflog)

  # In case any of the raw data types pulled need to be combined in some derivative manner
  # Loop through the symbols
  for obs in keys(m.observables)
      if get(default_transforms, m.transf[obs], 0) != 0
        m.data_transforms[obs] = default_transforms[m.transf[obs]]
      elseif m.transf[obs] == "1"
        m.data_transforms[obs] = function (levels)
          return levels # do nothing, leave in levels the way it was downloaded raw
        end
      # not a default_transforms - anything that combines multiple series in unique way
      elseif obs==:obs_wages
        # 1. Real wage growth
        m.data_transforms[:obs_wages] = function (levels)
        # FROM: Nominal compensation per hour (:COMPNFB from FRED)
        # TO: quarter to quarter percent change of real compensation (using GDP deflator)
        oneqtrpctchange(nominal_to_real(:COMPNFB, levels))
        end # transformation function
      elseif obs==:obs_ADD
        m.data_transforms[:obs_ADD] = function (levels)
        end
      else
        error("$(obs)--> no user-defined function for transformation (need String!=0 in subspec.jl or a user-defined function in the m###.jl file)")
      end
  end # obs


  # Look to truncate the set at a later point
  # Prior to actually conducting the transformations, this is just writing them out
    # for obs in keys(m.observables)
    # if m.startDate[obs] != ""
      # transform all the values prior to a specific date to NaNs

  # m.data_transforms[obs] = function (levels)
  #   return levels # do nothing, leave in levels the way it was downloaded raw
  # end

end # function
