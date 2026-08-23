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
    filldata(raw::TimeArray; method=Default, threshold = 2, α=0.3, startvals=Float64[], endvals=Float64[])

Fill a time series using LOESS with suggested start and end values or linear extrapolations.
Values in `startvals`, if provided, will be used as the start values, if these are missing.
Calues in `endvals`, if provided, will be used as the end values, if these are missing.
If start or end values are missing and suggested values not supplied, linear extrapolation is used to fill them.

If method == Default, only missing values are filled.
If method == Denoise, only values that are significantly different from the smoothed value are replaced with the smoothed value, where significant is defined as `abs(smoothed - original) > threshold * std(smoothed - original)`.
If method == FullSmooth, all values are smoothed using LOESS
"""
function filldata(raw::TimeArray; method::FillMethod=Default, threshold = 2, α=0.3, nonzeroonly=true, startvals=Float64[], endvals=Float64[])

    HoL = calcHoL(timestamp(raw))
    fulldata = similar(values(raw))
    data = zeros(nonmissingtype(eltype(values(raw))), size(values(raw)))

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
            data[:, i] = coalesce.(rawvals, fulldata[:, i])
        end

        if denoise
            # Replace points that are too different from the smoothed value
            # Find the standard deviation for the difference between smoothed and original
            σ = std(fulldata[:, i] .- data[:, i])
            data[:, i] .= ifelse.(abs.(fulldata[:, i] .- data[:, i]) .< (threshold * σ), data[:, i], fulldata[:, i])
        end
    end
    
    nonzeroonly && clamp!(data, 0, Inf) # Clamp negative values to zero
    return TimeArray(timestamp(raw), data, colnames(raw))
end


function filldata!(input_filename, output_filename; method::FillMethod=Default, threshold = 2, α=0.3, nonzeroonly=true, startvals=Float64[], endvals=Float64[])

    data, header = readdlm(lowercase(input_filename), ',', '\n', header=true)
    comps = string.(header[2:end]) # readdlm returns an array of AbstractStrings for some reason
    
    temp = Array{Union{Float64, Missing}}(undef, size(data, 1))

    for i in eachindex(data)
        if data[i] == ""
            data[i] = missing
        end
    end

    data = sortslices(data, dims=1, lt=(x,y)->isless(x[1],y[1])) #sort by timestamps
    timestamps = DateTime.(data[:, 1], "yyyy/mm/dd HH:MM")
    flows = data[:, 2:end]
    flows = convert(Matrix{Union{Float64, Missing}}, flows)
    raw = TimeArray(timestamps, flows, comps)

    filled = filldata(raw; method=method, threshold=threshold, α=α, nonzeroonly=nonzeroonly, startvals=startvals, endvals=endvals)
    writetimearray(filled, lowercase(output_filename), format="yyyy/mm/dd HH:MM", delim=',')

    return nothing
end
