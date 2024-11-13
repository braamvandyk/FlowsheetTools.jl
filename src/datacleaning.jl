using Dates, Loess, Interpolations, Missings, TimeSeries, Statistics

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
    basedata = zeros(Float64, length(times))
    data = zeros(Union{Float64, Missing}, length(times))

    for i in eachindex(HoL)
        basedata[i] = sin(π*period*(HoL[i]/(Hour(endtime - starttime)/Hour(1))))
        norm = Normal(0, 0.2*abs(basedata[i]))
        data[i] = basedata[i] + rand(norm)
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
    
    return data, basedata
end



"""
    filldata(raw; allsmoothed=false, denoise=false, threshold = 2, α=0.3, 
        suggest_start=false, startvals=Float64[], suggest_end=false, endvals=Float64[])

Fill a time series using LOESS with suggested start and end values or linear extrapolations.
If `suggest_start = true`, the values in `startvals` will be used as the start values, if these are missing.
If `suggest_end = true`, the values in `endvals` will be used as the end values, if these are missing.
If start or end values are missing and suggested values not supplied, linear extrapolation is used to fill them.

If `denoise` is true, datapoints will be replaced with the smoothed value when `abs(smoothed - original) > sensitivity * abs(smoothed)`.
If `allsmoothed` is true, all values are smoothed using LOESS, otherwise only missings are filled.

"""
function filldata(raw; fullsmooth=false, denoise=false, threshold = 2, α=0.3, 
    suggest_start=false, startvals=Float64[], suggest_end=false, endvals=Float64[])

    HoL = calcHoL(timestamp(raw))
    fulldata = similar(values(raw))
    data = zeros(nonmissingtype(eltype(values(raw))), length(values(raw)))

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
        model = loess(_HoL, _data; span=α, degree=2)
        
        fulldata[:, i] .= predict(model, HoL)

        if fullsmooth
            data[:, i] .= fulldata[:, i]
        else
            # Fill in missing data only
            data[:, i] .= coalesce.(rawvals, fulldata[:, i])
        end

        if denoise
            # Replace points that are too different from the smoothed value
            # Find the standard deviation for the difference between smoothed and original
            σ = std(fulldata[:, i] .- data[:, i])
            data[:, i] .= ifelse.(abs.(fulldata[:, i] .- data[:, i]) .< (threshold * σ), data[:, i], fulldata[:, i])
        end
    end
    
    return TimeArray(timestamp(raw), data, colnames(raw))
end


# Generate dummy data with missing values
starttime = DateTime(2023, 1, 1, 0, 0)
endtime = DateTime(2023, 1, 12, 24, 0)
times = starttime:Hour(6):endtime

rawvals, purevals = gendata(times, 2, 0.75, 0.5)
raw = TimeArray(times, rawvals, [:raw])
pure = TimeArray(times, purevals, [:pure])
# raw = TimeArray(times, genstepdata(times, 20, 0.75, 0.5, 10) .+ 50.0, [:raw1])

sf1 = filldata(raw)
rename!(sf1, [:default02])

sf2 = filldata(raw, denoise = true)
rename!(sf2, [:denoise02])

sf3 = filldata(raw, fullsmooth = true, suggest_start=true, suggest_end=true, startvals=[0.0, 0.0], endvals=[0.0, 0.0])
rename!(sf3, [:fullsmooth02])

sf4 = filldata(raw, α=0.5)
rename!(sf4, [:default05])

sf5 = filldata(raw, denoise = true, α=0.5)
rename!(sf5, [:denoise05])

sf6 = filldata(raw, fullsmooth = true, α=0.5, suggest_start=true, suggest_end=true, startvals=[0.0, 0.0], endvals=[0.0, 0.0])
rename!(sf6, [:fullsmooth05])

# https://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf

# begin
#     scatter(raw, ms=6, label="raw")
#     scatter!(sf1, label="default", marker=:square)
#     scatter!(sf2, label="denoise", marker=:diamond)
#     scatter!(sf3, label="allsmoothed", ms=2)
# end

begin
    alldata = TimeSeries.merge(sf1, sf2, sf3, sf4, sf5, sf6)
    l = @layout [a b c d;
                 e f g h]

    pltraw = let
        scatter(alldata[:raw], leg=:bottomleft, size=(640, 480));
        plot!(pure[:pure])
    end;

    pltdef02 = let 
        scatter(alldata[:default02], leg=:bottomleft, size=(640, 480))
        plot!(pure[:pure])
    end;

    pltnoise02 = let 
        scatter(alldata[:denoise02], leg=:bottomleft, size=(640, 480))
        plot!(pure[:pure])
    end;

    pltsmooth02 = let 
        scatter(alldata[:fullsmooth02], leg=:bottomleft, size=(640, 480))
        plot!(pure[:pure])
    end;

    pltdef05 = let 
        scatter(alldata[:default05], leg=:bottomleft, size=(640, 480))
        plot!(pure[:pure])
    end;

    pltnoise05 = let 
        scatter(alldata[:denoise05], leg=:bottomleft, size=(640, 480))
        plot!(pure[:pure])
    end;

    pltsmooth05 = let 
        scatter(alldata[:fullsmooth05], leg=:bottomleft, size=(640, 480))
        plot!(pure[:pure])
    end;

    plot(pltraw, pltdef02, pltnoise02, pltsmooth02, pltraw, pltdef05, pltnoise05, pltsmooth05, layout = l, size=(2560, 960))
end
savefig("cleandemo.png")