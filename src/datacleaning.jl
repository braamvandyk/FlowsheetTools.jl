# # These are used only during testing
# using Dates, Loess, Interpolations, Missings, TimeSeries, Statistics
# using Plots, Distributions

"""
    calcHoL(timestamps)

Calculates the hours on-line for a series of timestamps, assuming the first timestamp = 0 hours on-line.
"""
function calcHoL(timestamps)
    starttime = timestamps[begin]
    endtime = timestamps[end]

    return (timestamps .- starttime) ./ Hour(1)
end


@enum FillMethod begin
    Default
    Denoise
    FullSmooth
end

"""
    filldata(raw; method=Default, threshold = 2, α=0.3, startvals=Float64[], endvals=Float64[])

Fill a time series using LOESS with suggested start and end values or linear extrapolations.
Values in `startvals`, if provided, will be used as the start values, if these are missing.
Calues in `endvals`, if provided, will be used as the end values, if these are missing.
If start or end values are missing and suggested values not supplied, linear extrapolation is used to fill them.

If method == Default, only missing values are filled.
If method == Denoise, only values that are significantly different from the smoothed value are replaced with the smoothed value, where significant is defined as `abs(smoothed - original) > threshold * std(smoothed - original)`.
If method == FullSmooth, all values are smoothed using LOESS
"""
function filldata(raw::TimeArray; fullsmooth=false, denoise=false, threshold = 2, α=0.3, 
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
            if length(startvals) >= 1
                pushfirst!(_HoL, HoL[begin])
                pushfirst!(_data, startvals[i])
            else   
                pushfirst!(_HoL, HoL[begin])
                pushfirst!(_data, extrap(HoL[begin]))
            end
        end

        if !any(x -> x == HoL[end], _HoL)
            if length(endvals) >= 1
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

        if method == Default
            fullsmooth = false
            denoise = false
        elseif method == Denoise
            fullsmooth = false
            denoise = true
        elseif method == FullSmooth
            fullsmooth = true
            denoise = false
        end


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

# function filldata(rawstream::Stream, basis=:mass; fullsmooth=false, denoise=false, threshold = 2, α=0.3, 
#     suggest_start=false, startvals=Float64[], suggest_end=false, endvals=Float64[])

#     if basis == :mass
#         rawvals = values(rawstream.massflows)
#         colnames = colnames(rawstream.massflows)
#     else # basis == :mole
#         rawvals = values(rawstream.moleflows)
#         colnames = colnames(rawstream.moleflows)
#     end

#     rawdata = TimeArray(timestamp(rawstream), rawvals, colnames)
#     filled = filldata(rawdata; fullsmooth=fullsmooth, denoise=denoise, threshold=threshold, α=α, suggest_start=suggest_start, startvals=startvals,
#         suggest_end=suggest_end, endvals=endvals)

# end


# # Generate dummy data with missing values
# starttime = DateTime(2023, 1, 1, 0, 0)
# endtime = DateTime(2023, 1, 12, 24, 0)
# times = starttime:Hour(6):endtime

# rawvals, purevals = gendata(times, 2, 0.75, 0.5)
# raw = TimeArray(times, rawvals, [:raw])
# pure = TimeArray(times, purevals, [:pure])
# # raw = TimeArray(times, genstepdata(times, 20, 0.75, 0.5, 10) .+ 50.0, [:raw1])

# sf1 = filldata(raw)
# rename!(sf1, [:default03])

# sf2 = filldata(raw, method = Denoise)
# rename!(sf2, [:denoise03])

# sf3 = filldata(raw, method = FullSmooth, startvals=[0.0, 0.0], endvals=[0.0, 0.0])
# rename!(sf3, [:fullsmooth03])

# sf4 = filldata(raw, α=0.5)
# rename!(sf4, [:default05])

# sf5 = filldata(raw, method = Denoise, α=0.5)
# rename!(sf5, [:denoise05])

# sf6 = filldata(raw, method = FullSmooth, α=0.5, startvals=[0.0, 0.0], endvals=[0.0, 0.0])
# rename!(sf6, [:fullsmooth05])

# alldata = TimeSeries.merge(sf1, sf2, sf3, sf4, sf5, sf6, method=:left)

# # https://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf


# begin
#     alldata = TimeSeries.merge(sf1, sf2, sf3, sf4, sf5, sf6)
#     l = @layout [a b c d;
#                  e f g h]

#     pltraw = let
#         scatter(raw, leg=:bottomleft, size=(640, 480));
#         plot!(pure[:pure])
#     end;

#     pltdef03 = let 
#         scatter(alldata[:default03], leg=:bottomleft, size=(640, 480))
#         plot!(pure[:pure])
#     end;

#     pltnoise03 = let 
#         scatter(alldata[:denoise03], leg=:bottomleft, size=(640, 480))
#         plot!(pure[:pure])
#     end;

#     pltsmooth03 = let 
#         scatter(alldata[:fullsmooth03], leg=:bottomleft, size=(640, 480))
#         plot!(pure[:pure])
#     end;

#     pltdef05 = let 
#         scatter(alldata[:default05], leg=:bottomleft, size=(640, 480))
#         plot!(pure[:pure])
#     end;

#     pltnoise05 = let 
#         scatter(alldata[:denoise05], leg=:bottomleft, size=(640, 480))
#         plot!(pure[:pure])
#     end;

#     pltsmooth05 = let 
#         scatter(alldata[:fullsmooth05], leg=:bottomleft, size=(640, 480))
#         plot!(pure[:pure])
#     end;

#     plot(pltraw, pltdef03, pltnoise03, pltsmooth03, pltraw, pltdef05, pltnoise05, pltsmooth05, layout = l, size=(2560, 960))
# end
# savefig("cleandemo.png")

# end # module
