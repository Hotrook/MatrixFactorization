# Sebastian Mroz 
# 221433 
function iterDFS(j::Int64, postCounter::Int64, head::Array{Int64, 1}, next::Array{Int64, 1}, post::Array{Int64, 1}, n::Int64, stack::Array{Int64, 1})
    stack[ 1 ] = j
    top = 1

    while top >= 1
        parent = stack[ top ]
        child = head[ parent ]
        if child == 0
            top = top - 1
            post[ postCounter ] = parent
            postCounter = postCounter + 1
        else
            head[ parent ] = next[ child ]
            top = top + 1
            stack[ top ] = child
        end
    end

    return postCounter
end
