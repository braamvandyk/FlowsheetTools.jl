#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------

struct Flowsheet
    unitops::Vector{UnitOp}
    order::Vector{Int}
end

function (fs::Flowsheet)(neworder = nothing)
    if isnothing(neworder)
        for i in fs.order
            fs.unitops[i]()
        end
    else
        for i in neworder
            u = fs.unitops[i]
            u()
        end
    end

    return nothing
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


function addunitop!(fs::Flowsheet, u::UnitOp)
    push!(fs.unitops, u)
    push!(fs.order, length(fs.order) + 1)

    return nothing
end


function addunitop!(fs::Flowsheet, us::Vector{UnitOp})
    for u in us
        push!(fs.unitops, u)
        push!(fs.order, length(fs.order) + 1)
    end

    return nothing
end


function setorder!(fs::Flowsheet, neworder)
    numold = length(fs.order)
    numnew = length(neworder)

    if numold == numnew
        for i in 1:numold
            fs.order[i] = neworder[i]
        end
    elseif numold < numnew
        for i in 1:numold
            fs.order[i] = neworder[i]
        end            
        for j = 1:(numnew - numold)
            push!(fs.order, neworder[numold+j])
        end
    else
        for i in 1:numold
            fs.order[i] = neworder[i]
        end            
        for j = 1:(numold - numnew)
            pop!(fs.order)
        end
    end

    return nothing
end


# function generateBFD(fs:Flowsheet)
#     # 1. Find all the feed and product streams. We need to assign nodes to them for the diagram, olr they won't get rendered
#     # 2. Iterate through the blocks and add them to the 
    
#     for i in fs.order
#         u = fs.unitops[i]
