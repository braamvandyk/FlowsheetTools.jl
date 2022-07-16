"""
    function prettyround(x, n=3)

Rounds to n significant digits for x < 1, or n digits for x >= 1.
Used to pretty up the display of flows.
"""
function prettyround(x, n=3)
    if x < 1
        return round(x, sigdigits = n)
    else
        return round(x, digits = n)
    end
end