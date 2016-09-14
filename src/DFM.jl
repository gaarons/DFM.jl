isdefined(Base, :__precompile__) && __precompile__()

module DFM
    # using FredData, DataFrames, Base.Dates
    # using Interpolations
    # using Distributions, Roots.fzero, HDF5
    # using DataStructures: SortedDict, insert!, ForwardOrdering, OrderedDict
    # using QuantEcon: solve_discrete_lyapunov

    using Distributions, Roots.fzero, HDF5
    using DataStructures: SortedDict, insert!, ForwardOrdering, OrderedDict
    using FredData, DataFrames, Base.Dates
    using QuantEcon: solve_discrete_lyapunov
    import Calculus
    import Optim
    using Optim: OptimizationTrace, OptimizationState, MultivariateOptimizationResults


    import Calculus
    import Optim
    using Optim: OptimizationTrace, OptimizationState, MultivariateOptimizationResults


# Keep getting warnings to explicitly import the julia Base package
    # import Base

# Things you might want to use outside the module, when calling module DFM
    export
    #     # abstractdfmmodel.jl
        AbstractModel, description, spec, subspec,
        dataroot, saveroot, inpath, workpath, rawpath, tablespath, figurespath, logpath,
        data_vintage,
    # models/
        Model990,
    # settings.jl
        Setting, get_setting, default_settings!, default_test_settings!,

    # estimate/
        estimateDFM, em_converged!, InitCond, EMstep,
        remNaNs_spline!,
        runKF!,SKF!,FIS!,MissData!,

    # types, not sure if i need these or if can be imbedded in the function definition
        KalmanType, ResultsType,


    #     kalman_filter, likelihood, posterior, posterior!,
    # optimize!, csminwel, hessian!, estimate, proposal_distribution,
    # metropolis_hastings, compute_parameter_covariance, compute_moments,
    # find_density_bands, prior,

    # data/
        load_data, load_data_levels, load_fred_data, transform_data, save_data,
        df_to_matrix, hpfilter, difflog, quartertodate, percapita, nominal_to_real,
        hpadjust, oneqtrpctchange, annualtoquarter


    # transform_to_model_space!, transform_to_real_line!,
    #     n_states, n_states_augmented, n_shocks_exogenous, n_shocks_expectational,
    #     n_equilibrium_conditions, n_observables, n_parameters, n_parameters_steady_state,
    #     n_parameters_free, n_anticipated_shocks,
    #
    #     reoptimize, calculate_hessian,
    #     n_mh_blocks, n_mh_simulations, n_mh_burn, mh_thin, specify_mh_start,
    #
    #     specify_mode!, specify_hessian, load_parameters_from_file,
    #     use_population_forecast,

    #
    #     # distributions_ext.jl
    #     BetaAlt, GammaAlt,
    #
    #
    #
    #     # parameters.jl
    #     parameter, Transform, NullablePrior, AbstractParameter,
    #     Parameter, ParameterVector, ScaledParameter,
    #     UnscaledParameter, SteadyStateParameter, transform_to_real_line, transform_to_model_space,
    #     update, update!, transform_to_model_space, transform_to_real_line, Interval, ParamBoundsError,
    #


    #     # solve/
    #     gensys, solve,
    #

    const VERBOSITY = Dict(:none => 0, :low => 1, :high => 2)
    const DFM_DATE_FORMAT = "yyyymmdd"


# Establish what the model types and basic functions are
    include("abstractdfmmodel.jl")
    include("settings.jl")
    include("defaults.jl")
    include("util.jl")

# Get the data into dataframes and prepared for use
    include("data/load_data.jl")
    include("data/fred_data.jl")
    include("data/transformations.jl")
    include("data/transform_data.jl")
    include("data/util.jl")

# Preferred initialization for the DFM algorithm to converge to Quasi-Maximum likelihood results
    include("estimate/remNaNs_spline.jl")


    include("models/m990/m990.jl")
    include("models/m990/subspecs.jl")

    include("estimate/runKF.jl")
    include("estimate/estimateDFM.jl")





    # include("estimate/posterior.jl")
    # include("estimate/optimize.jl")
    # include("estimate/csminwel.jl")
    # include("estimate/hessian.jl")
    # include("estimate/hessizero.jl")
    # include("estimate/estimate.jl")
    # include("estimate/moments.jl")
    #
    #  include("models/m990/eqcond.jl")
    # include("models/m990/measurement.jl")
    # include("models/m990/augment_states.jl")
    #
    # include("models/smets_wouters/smets_wouters.jl")
    # include("models/smets_wouters/subspecs.jl")
    # include("models/smets_wouters/eqcond.jl")
    # include("models/smets_wouters/measurement.jl")
    # include("models/smets_wouters/augment_states.jl")

end
