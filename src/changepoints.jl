"""
    function mad(data)

Return the median absolute deviation of `data`.
"""
function mad(data)
    return median(abs.(data .- median(data)))
end



"""


Wild Bindary Segmentation as per Fryzlewicz, P. (2014) Wild binary segmentation for multiple change-point detection, Annals of Statistics 42(6), 2243-2281

M: Number of samples, default at 5000
C: Threshold parameter
"""
function WBS(data, M = 5000, C=1.0, σ = mad(data))
    T = length(data)
    basesum = cumsum(data)

    function CUSUM(b, s, e)
        n = e - s + 1
        sum1 = basesum[b]
        sum2 = basesum[e] - sum1
    
        return abs(sqrt((e-b)/(n*(b-s+1)))*sum1 - sqrt((b-s+1)/(n*(e-b)))*sum2)/σ
    end

    function sSIC()
        n = length(data)
         return (T/2) * log(median(abs.(data[s:t] .- mean(data[s:t])))^2)
        return cost
    end
    

