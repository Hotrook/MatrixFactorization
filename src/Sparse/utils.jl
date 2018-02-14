# Sebastian Mroz 
# 221433 
function FLIP( i )
    -(i)-2
end

function UNFLIP( i )
    if i < 0
        return FLIP(i)
    else
        return i
    end
end

function MARKED(w, j)
    w[ j ] < 0
end

function MARK(w, j)
    w[ j ] = FLIP( w[ j ] )
end
