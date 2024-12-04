using Dates, TimeSeries, Statistics, HypothesisTests
using Distributions, Missings




function genstepdata(timestamps, period, fracfilled, fracdouble, stepsize=0.0)
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
    
    stepidx = rand(1:len)
    data[stepidx:end] .= data[stepidx:end] .+= stepsize # Add a step change to the end of the data
    return data
end


# function findchangepoints(data, window=72)

    function getwindow(point)
        allleft = findall(x -> x < HoL[point] - window, HoL)
        left =  length(allleft) + 1
        checkright = findfirst(x -> x > HoL[point] + window, HoL)
        if isnothing(checkright)
            right = length(HoL)
        else
            right = checkright - 1
        end

        return left, right
    end

    
    starttime = DateTime(2023, 1, 1, 0, 0)
    endtime = DateTime(2023, 1, 12, 24, 0)
    times = starttime:Hour(6):endtime
    raw = TimeArray(times, genstepdata(times, 20, 0.75, 0.5, 10) .+ 50.0, [:raw1])
    clean = filldata(raw)
    scatter(raw)
    scatter!(clean, ms=2)
    
    HoL = calcHoL(timestamp(clean))
    vals = values(clean)

    pvals = Float64[]
    for point in eachindex(HoL)
    # point = 25
    # window = 72
        if point == 1 || point == length(HoL)
            continue
        end

        left, right = getwindow(point)
        leftdata = vals[left:point]
        rightdata = vals[point:right]
        leftmean = mean(leftdata)
        leftstd = std(leftdata)
        rightmean = mean(rightdata)
        rightstd = std(rightdata)

        pval = pvalue(UnequalVarianceTTest(leftdata, rightdata))
        push!(pvals, pval)

        println("$point \t $left \t $right \t $(leftmean) \t $(rightmean) \t $pval")
    end


# end



raw = TimeArray(times, genstepdata(times, 20, 0.75, 0.5, 2) .+ 50.0, [:raw1])
data = filldata(raw)
# scatter(data)
# scatter!(raw)