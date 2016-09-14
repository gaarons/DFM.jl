"""
`prev_quarter(q::TimeType = now())`

Returns Date identifying last day of the previous quarter
"""
function prev_quarter(q::TimeType = now())
    Date(lastdayofquarter(firstdayofquarter(q)-Dates.Day(1)))
end

"""
`next_quarter(q::TimeType = now())`

Returns Date identifying last day of the next quarter
"""
function next_quarter(q::TimeType = now())
    Date(lastdayofquarter(lastdayofquarter(q)+Dates.Day(1)))
end

"""
`get_quarter_ends(start_date::Date,end_date::Date)`

Returns a DataArray of quarter end dates between `start_date` and `end_date`.
"""
function get_quarter_ends(start_date::Date,end_date::Date)
    map(lastdayofquarter, collect(start_date:Dates.Month(3):end_date))
end

"""
`stringstodates(stringarray)`

Converts a collection of strings in "y-m-d" format to Dates.
"""
function stringstodates(stringarray)
    n = length(stringarray)
    dates = Vector{Date}(n)
    for i = 1:n
        dates[i] = Date(stringarray[i], "y-m-d")
    end

    dates
end

"""
`quartertodate(string::AbstractString)`

Convert `string` in the form "YYqX", "YYYYqX", or "YYYY-qX" to a Date of the end of
the indicated quarter. "X" is in `{1,2,3,4}` and the case of "q" is ignored.
"""
function quartertodate(string::AbstractString)
    if ismatch(r"^[0-9]{2}[qQ][1-4]$", string)
        year = "20"*string[1:2]
        quarter = string[end]
    elseif ismatch(r"^[0-9]{4}[qQ][1-4]$", string)
        year = string[1:4]
        quarter = string[end]
    elseif ismatch(r"^[0-9]{4}-[qQ][1-4]$", string)
        year = string[1:4]
        quarter = string[end]
    else
        throw(ParseError("Invalid format: $string"))
    end

    year = parse(Int, year)
    quarter = parse(Int, quarter)
    month = 3*quarter
    day = 1

    return lastdayofquarter(Date(year, month, day))
end

"""
`subtract_quarters(t1::Date, t0::Date)`

Compute the number of quarters between t1 and t0, including t0 and excluding t1.
"""
function subtract_quarters(t1::Date, t0::Date)
    days = t1 - t0
    quarters = round(days.value / 365.25 * 4.0)
    return convert(Int, quarters)
end


"""
`format_dates!(col, df)`

Change column `col` of dates in `df` from String to Date, and map any dates given in the
interior of a quarter to the last day of the quarter.
"""
function format_dates!(col::Symbol, df::DataFrame)
    df[col] = stringstodates(df[col])
    map!(lastdayofquarter, df[col], df[col])
end

"""
```
na2nan!(df::DataFrame)
```

Convert all NAs in a DataFrame to NaNs.
"""
function na2nan!(df::DataFrame)
    for col in names(df)
        df[isna(df[col]), col] = NaN
    end
end

"""
```
na2nan!(df::DataFrame)
```

Convert all NAs in a DataFrame to NaNs.
"""
function na2nan!(v::DataArray)
    for i = 1:length(v)
        v[i] = isna(v[i]) ?  NaN : v[i]
    end
end
