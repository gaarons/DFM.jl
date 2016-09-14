"""
```
default_settings!(m::AbstractModel)
```

Default Settings are constructed, initialized and added to `m.settings`.
"""
function default_settings!(m::AbstractModel)

    # I/O File locations
    saveroot = normpath(joinpath(dirname(@__FILE__), "..","save"))
    datapath = normpath(joinpath(dirname(@__FILE__), "..","save","input_data"))

    m <= Setting(:saveroot, saveroot, "Root of data directory structure")
    m <= Setting(:dataroot, datapath, "Input data directory path")

    # Data settings for released and conditional data. Default behavior is to set vintage
    # of data to today's date.
    vint = Dates.format(now(), DFM_DATE_FORMAT)
    m <= Setting(:data_vintage, vint, true, "vint", "Vintage of data")



    # Default to run a Monthly model with 1 basket == global
    m <= Setting(:frequency, "m", "Frequency of the model")
    # Complicated setting, as this is a product of the model code we are running
    # m <= Setting(:num_baskets, 1, "Number of factor baskets")
    # fully determined in the subspec

    # Dates
    m <= Setting(:date_sample_start, quartertodate("1960-Q1"), "Start date of main sample")
    m <= Setting(:date_sample_end, Dates.lastdayofquarter(Dates.today()-Dates.Month(3)), "End date of main sample")
    m <= Setting(:date_fcst_start, Dates.lastdayofquarter(Dates.today()), "Start date of forecast period")
    m <= Setting(:date_fcst_end, Dates.lastdayofquarter(Dates.today()+Dates.Month(60*3)), "End date of forecast period")

    # General computation
    m <= Setting(:use_parallel_workers, true, "Use available parallel workers in computations")
end
