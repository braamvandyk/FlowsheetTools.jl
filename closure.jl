#----Closure-----------------

"""
    function closemb(boundary::BalanceBoundary; totalweight=1.0, elementweight=1.0)
    function closemb(boundary::BalanceBoundaryHistory; totalweight=1.0, elementweight=1.0)

Basic mass balance reconciliation. Error in total mass closure and average element closures are
weighted by *totalweight* and *elementweight* respectively and the squared weighted error is minimized.

Results are returned as a dict of streams and corrections to their flows.
"""
function closemb(boundary::BalanceBoundary; totalweight=1.0, elementweight=1.0)         
    corrections = Dict{String, Float64}()
    inlets = boundary.inlets
    outlets = boundary.outlets

    ins = length(inlets)
    outs = length(outlets)
    factors = ones(ins + outs)

   
    function f(factors)
        total_in = sum(factors[1:ins] .* inlets)
        total_out = sum(factors[ins+1:end] .* outlets)
        masserror = abs2(total_out.totalmassflow/total_in.totalmassflow - 1.0)

        atomclosures = Dict{String, Float64}()
        for atom in keys(total_in.atomflows)
            in = total_in.atomflows[atom]
            out = total_out.atomflows[atom]
            atomclosures[atom] = out/in
        end
        numatoms = length(atomclosures)
        atomerror = abs2(sum(values(atomclosures)) - numatoms)
        return totalweight*masserror + elementweight*atomerror
    end

    res = optimize(f, factors)
    optfactors = res.minimizer

    i = 1
    for stream in inlets
        corrections[stream.name] = optfactors[i]
        i += 1
    end
    for stream in outlets
        corrections[stream.name] = optfactors[i]
        i += 1
    end

    return corrections
end


function closemb(boundary::BalanceBoundaryHistory; totalweight=1.0, elementweight=1.0)         
    corrections = Dict{String, Float64}()
    inlets = boundary.inlets
    outlets = boundary.outlets

    ins = length(inlets)
    outs = length(outlets)
    factors = ones(ins + outs)

    numdata = inlets[1].numdata
   
    function f(factors)
        # Since inlets and outlets are arrays of StreamHistory, summing them produces StreamHistory objects
        total_in = sum(factors[1:ins] .* inlets)
        total_out = sum(factors[ins+1:end] .* outlets)
        
        masserror = 0.0
        atomerror = 0.0
        for datum = 1:numdata
            masserror += abs2(total_out.totalmassflow[datum]/total_in.totalmassflow[datum] - 1.0)
            
            atomclosures = Dict{String, Float64}()
            for atom in keys(total_in.atomflows[datum])
                in = total_in.atomflows[datum][atom]
                out = total_out.atomflows[datum][atom]
                atomclosures[atom] = out/in
            end
            numatoms = length(atomclosures)
            atomerror += abs2(sum(values(atomclosures)) - numatoms)
        end
        return totalweight*masserror + elementweight*atomerror
    end

    res = optimize(f, factors)
    optfactors = res.minimizer

    i = 1
    for stream in inlets
        corrections[stream.name] = optfactors[i]
        i += 1
    end
    for stream in outlets
        corrections[stream.name] = optfactors[i]
        i += 1
    end

    return corrections
end
