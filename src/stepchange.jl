# TODO Kill this code and replace with ChangePointDetection.changepoints(ts; threshold = 0.5, window = 150)...

"""
    function findsteps(data)

Finds step changes in the data. When reconciling data, you should first do step change detection,
then reconcile between steps to make sure you can calculate meaningful uncertainties in the data
"""
function findsteps(data; window=3, noiseperc=0.9)
    fwdmean = similar(data)
    revmean = similar(data)

    numdata = length(data)
    @assert numdata > 2*window "Too few data points ($numdata) for the selectived window size ($window)."
    
    # Calculate the forward and reverse moving averages `window` points from i
    # i.e. the actual window is `window + 1` points wide
    for i ∈ eachindex(data)
        if i <= window
            revmean[i] = mean(data[1:i])
            fwdmean[i] = mean(data[i:(i+window)])
        elseif window < i <= (numdata - window)
            revmean[i] = mean(data[(i-window):i])
            fwdmean[i] = mean(data[i:(i+window)])
        else
            revmean[i] = mean(data[(i-window):i])
            fwdmean[i] = mean(data[i:numdata])
        end
    end

    # Take the difference between forwards and reverse moving averages
    # then kill everthing smaller than the specified quantile to remove noise
    deltas = abs.(fwdmean .- revmean)
    cutoff = quantile(deltas, noiseperc)

    for i ∈ eachindex(deltas)
        if deltas[i] < cutoff
            deltas[i] = 0.0
        end
    end

    # Now look for the largest deltas in forward and reverse moving windows again
    # This allows an arbitrary number of local maxima to be found
    fwdmax = similar(data)
    revmax = similar(data)

    for i ∈ eachindex(deltas)
        if i <= window
            revmax[i] = maximum(deltas[1:i])
            fwdmax[i] = maximum(deltas[i:(i+window)])
        elseif window < i <= (numdata - window)
            revmax[i] = maximum(deltas[(i-window):i])
            fwdmax[i] = maximum(deltas[i:(i+window)])
        else
            revmax[i] = maximum(deltas[(i-window):i])
            fwdmax[i] = maximum(deltas[i:numdata])
        end
    end

    steps = falses(numdata)

    # And the step changes are at the points where the delta = the local maxima and
    # at least one of the local maxima is > 0
    for i ∈ 1:numdata
        if (revmax[i] + fwdmax[i]) > 0.0
            if deltas[i] == max(revmax[i], fwdmax[i])
                steps[i] = true
            end
        end
    end

    # Find all the steps - reported as the index of the first point AFTER the step,
    # i.e. the first point of the new group
    results = findall(steps)

    return results
end