# Sebastian Mroz 
# 221433 

function postordering(parent::Array{Int, 1})
    n = size(parent)[1]
    head = zeros(Int, n)
    next = zeros(Int, n)
    post = zeros(Int, n)
    stack = zeros( Int, n )
    
    postCounter = 1

    for j = n:-1:1
        if parent[ j ] == 0
            continue
        end
        next[ j ] = head[ parent[ j ] ]
        head[ parent[ j ] ] = j
    end

    for j = 1:n
        if parent[ j ] != 0
            continue
        end
        postCounter = iterDFS(j, postCounter, head, next, post, n, stack)
    end

    return post
end
