using Dates, Loess, Interpolations, Missings, TimeSeries

# # These are used only during testing
using Plots, Distributions

"""
    calcHoL(timestamps)

Calculates the hours on-line for a series of timestamps, assuming the first timestamp = 0 hours on-line.
"""
function calcHoL(timestamps)
    starttime = timestamps[begin]
    endtime = timestamps[end]

    return (timestamps .- starttime) ./ Hour(1)
end


# Generate dummy data for testing
function gendata(timestamps, period, fracfilled, fracdouble)
    # We use HoL since we need the x-axis to run from 0.0
    # The timestamps are also converted into a Float64 with hours since start
    HoL = calcHoL(timestamps)
    data = zeros(Union{Float64, Missing}, length(times))

    for i in eachindex(HoL)
        data[i] = sin(π*period*(HoL[i]/(Hour(endtime - starttime)/Hour(1))))
        norm = Normal(0, 0.2*abs(data[i]))
        data[i] += rand(norm)
    end

    len = length(data)
    totalmissing = round(Int, (1 - fracfilled)*len)

    # Add missing data
    nummissing = 0
    numdouble = 0

    while nummissing < totalmissing
        idx = rand(1:len)
        if ismissing(data[idx])
            continue
        end
        data[idx] = missing
        nummissing += 1
        if idx < len && rand() < fracdouble
            data[idx+1] = missing            
            numdouble += 1
        end
    end
    
    return data
end


"""
    filldata(raw; allsmoothed=false, denoise=false, sensitivity = 0.1, loessspan=0.3, 
        suggest_start=false, startvals=Float64[], suggest_end=false, endvals=Float64[])

Fill a time series using LOESS with suggested start and end values or linear extrapolations.
If `suggest_start = true`, the values in `startvals` will be used as the start values, if these are missing.
If `suggest_end = true`, the values in `endvals` will be used as the end values, if these are missing.
If start or end values are missing and suggested values not supplied, linear extrapolation is used to fill them.

If `denoise` is true, datapoints will be replaced with the smoothed value when `abs(smoothed - original) > sensitivity * abs(smoothed)`.
If `allsmoothed` is true, all values are smoothed using LOESS, otherwise only missings are filled.

"""
function filldata(raw; allsmoothed=false, denoise=false, sensitivity = 0.1, loessspan=0.3, 
    suggest_start=false, startvals=Float64[], suggest_end=false, endvals=Float64[])

    HoL = calcHoL(timestamp(raw))
    fulldata = similar(values(raw))
    data = similar(values(raw))

    for (i, col) in enumerate(colnames(raw))
        rawvals = values(raw[col])
        # Drop missing data entries in this column
        _HoL, _data = collect.(skipmissings(HoL, rawvals))
        
        # Prepare a linear interpolator that will extrapolate to start and end values
        extrap = linear_interpolation(_HoL, _data, extrapolation_bc = Line())

        # Add start and end values, if missing
        if !any(x -> x == HoL[begin], _HoL)
            if suggest_start
                pushfirst!(_HoL, HoL[begin])
                pushfirst!(_data, startvals[i])
            else   
                pushfirst!(_HoL, HoL[begin])
                pushfirst!(_data, extrap(HoL[begin]))
            end
        end
        if !any(x -> x == HoL[end], _HoL)
            if suggest_end
                push!(_HoL, HoL[end])
                push!(_data, endvals[i])
            else
                push!(_HoL, HoL[end])
                push!(_data, extrap(HoL[end]))
            end
        end

        # Fit a loess model to the data
        model = loess(_HoL, _data; span=loessspan, degree=2)
        
        fulldata[:, i] .= predict(model, HoL)

        if allsmoothed
            data[:, i] .= fulldata[:, i]
        else
            # Fill in missing data only
            data[:, i] .= coalesce.(rawvals, fulldata[:, i])
        end

        if denoise
            # Replace points that are too different from the smoothed value
                data[:, i] .= ifelse.(abs.(fulldata[:, i] .- data[:, i]) .< sensitivity .* abs.(fulldata[:, i]), data[:, i], fulldata[:, i])
        end
    end
    
    return TimeArray(timestamp(raw), data, colnames(raw))
end


# Generate dummy data with missing values
starttime = DateTime(2023, 1, 1, 0, 0)
endtime = DateTime(2023, 1, 12, 24, 0)
times = starttime:Hour(6):endtime

# raw = TimeArray(times, [gendata(times, 2, 0.75) gendata(times, 2, 0.75)], [:raw1, :raw2])
raw = TimeArray(times, gendata(times, 2, 0.75, 0.5), [:raw1])

sf1 = filldata(raw)
scatter(sf1, leg=:bottomleft)
scatter!(raw)

sf2 = filldata(raw, denoise = true, suggest_start=true, suggest_end=true, startvals=[0.0, 0.0], endvals=[0.0, 0.0])
scatter(sf2)
scatter!(raw)

sf3 = filldata(raw, allsmoothed = true, suggest_start=true, suggest_end=true, startvals=[0.0, 0.0], endvals=[0.0, 0.0])
scatter(sf3)
scatter!(raw)
