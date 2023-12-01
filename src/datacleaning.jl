# using Dates, Loess, Interpolations, Missings, TimeSeries

# # These are used only during testing
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


# Generate dummy data for testing
# function gendata(timestamps, period, fracfilled)
#     HoL = calcHoL(timestamps)
#     data = zeros(Union{Float64, Missing}, length(times))

#     for i in eachindex(HoL)
#         if rand() > fracfilled
#             data[i] = missing
#         else
#             data[i] = sin(π*period*(HoL[i]/(Hour(endtime - starttime)/Hour(1))))
#             norm = Normal(0, 0.1*abs(data[i]))
#             data[i] += rand(norm)
#         end
#     end

#     return data
# end

"""
    smooth_and_fill(raw; loessspan=0.3, suggest_start=false, startval=0.0, suggest_end=false, endval=0.0)

Smooth and fill a time series using loess, with suggested start and end values or linear extrapolations.
"""
function smooth_and_fill(raw; loessspan=0.3, suggest_start=false, startvals=Float64[], suggest_end=false, endvals=Float64[])
    HoL = calcHoL(timestamp(raw))
    data = similar(values(raw))
    
    for (i, col) in enumerate(colnames(raw))
        # Drop missing data entries
        _HoL, _data = collect.(skipmissings(HoL, values(raw[col])))
        
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
        data[:, i] .= predict(model, HoL)
    end

    return TimeArray(timestamp(raw), data, colnames(raw))
end


# # Generate dummy data with missing values
# starttime = DateTime(2023, 1, 1, 0, 0)
# endtime = DateTime(2023, 1, 12, 24, 0)
# times = starttime:Hour(6):endtime

# raw = TimeArray(times, [gendata(times, 2, 0.95) gendata(times, 2, 0.95)], [:raw1, :raw2])
# plot(raw, leg=:bottomleft)

# sf1 = smooth_and_fill(raw, suggest_start=true, suggest_end=true, startvals=[0.0, 0.0], endvals=[0.0, 0.0])
# sf2 = smooth_and_fill(raw)
# plot!(sf1)
# plot!(sf2)
