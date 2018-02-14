# Sebastian Mroz 
# 221433 
function isOfDiag( i, j, aij )
    return i != j
end

function fkeep0indexed( A::CCSparseMatrix )
    nz = 0

    for j = 0 : (A.n - 1)
        p = A.columnPointer[ j + 1 ]
        A.columnPointer[ j + 1 ] = nz
        while p < A.columnPointer[ j + 2 ]
            # if filter( A.rowIndex[ p + 1 ], j, A.value[ p + 1 ] )
            if A.rowIndex[ p + 1 ] != j
                A.rowIndex[ nz + 1 ] = A.rowIndex[ p + 1 ]
                nz += 1
            end
            p += 1
        end
    end
    A.columnPointer[ A.n + 1 ] = nz
    A.nzmax = nz
    return nz
end

function iterDFS0Indexed(j::Int64, postCounter::Int64, head::Array{Int64, 1}, next::Array{Int64, 1}, post::Array{Int64, 1}, n::Int64, stack::Array{Int64, 1})
    stack[ 1 ] = j
    top = 1

    while top >= 1
        parent = stack[ top ]
        child = head[ parent ]
        if child == -1
            top -= 1
            post[ postCounter ] = parent
            postCounter += 1
        else
            head[ parent ] = next[ child + 1 ]
            top += 1
            stack[ top ] = child + 1
        end
    end

    return postCounter
end

function downgrade!( A::CCSparseMatrix )
    for i = 1:(A.columnPointer[A.n+1]-1)
        A.rowIndex[ i ] -=1
    end
    for i = 1:(A.n+1)
        A.columnPointer[ i ] -= 1
    end
end



function amd( order::Integer, A::CCSparseMatrix)

    counter = 1
    dense = d = dk = dext = e = elenk = eln = i = j = k = k1 = k2 = k3 = jlast = ln = dens = nzmax = nvi = nvj = nvk = mark = wnvi = ok = cnz = p = p1 = p2 = p3 = p4 = pj = pk = pk1 = pk2 = pn = q = t = 0
    nel = 0
    mindeg = 0 # originally starts frm 0
    lemax = 0
    C = CCSparseMatrix(1,1,1)

    AT = transposeWithoutValues(A)
    m = A.m
    n = A.n

    dense = max(16, 10 * Int64(floor(sqrt( n ))))
    dense = min( n-2 , dense )

    if order == 1 && n == m
        C = addWithoutValues(A, AT, 0.0, 0.0)
    elseif order == 2
        p2 = 0
        for j = 0 : (m-1)
            p = AT.columnPointer[ j + 1 ]
            AT.columnPointer[ j + 1 ] = p2 + 1
            if AT.columnPointer[ j + 2 ] - p > dense
                continue
            end
            while p < AT.columnPointer[ j + 2 ]
                AT.rowIndex[ p2 + 1 ] = AT.rowIndex[ p ]
                p2 += 1
                p += 1
            end
        end
        AT.columnPointer[ m + 1 ] = p2 + 1
        A2 = transposeWithoutValues( AT )
        C = multiplyWithoutValues( AT, A2 )
    else
        C = multiplyWithoutValues( AT, A ) # QR
    end
    downgrade!( C )

    fkeep0indexed( C )
    cnz = C.columnPointer[ n + 1 ]
    t = (cnz) + Int64(floor((cnz)/5)) + 2*n
    resize!(C.rowIndex, t)
    C.nzmax = t

    P = fill( -1, n + 1 )
    len = zeros( Int64, n+1 )
    nv = ones( Int64, n+1 )
    next = fill( -1, n+1 )
    head = fill( -1, n+1 )
    elen = zeros( Int64, n+1 )

    degree = zeros( Int64, n+1 )
    w = ones( Int64, n+1 )
    hhead = fill( -1, n+1 )
    last = P
    for k = 1:n
        degree[ k ] = len[ k ] = C.columnPointer[ k + 1 ] - C.columnPointer[ k ]
    end
    degree[ n + 1 ] = len[ n + 1 ] = 0
    nzmax = C.nzmax

    mark = wclear(0, 0, w, n)

    elen[ n + 1 ] = -2
    C.columnPointer[ n + 1 ] = -1
    w[ n + 1 ] = 0

    for i = 0 : (n-1)
        d = degree[ i + 1 ]
        if d == 0
            elen[ i + 1  ] = -2
            nel += 1
            C.columnPointer[ i + 1 ] = -1
            w[ i + 1 ] = 0
        elseif d > dense
            nv[ i + 1 ] = 0
            elen[ i + 1 ] = -1
            nel += 1
            C.columnPointer[ i + 1 ] = FLIP( n ) # TODO: analyze argument of FLIP - n
            nv[ n+1 ] += 1
        else
            if head[ d+1 ] != -1
                 last[ head[ d + 1 ] + 1 ] = i
            end
            next[ i + 1 ] = head[ d + 1 ]
            head[ d + 1 ] = i
        end
    end




    while nel < n
        k = -1
        while mindeg < n && (k = head[ mindeg + 1 ]  ) == -1
            mindeg += 1
        end
        if next[ k + 1 ] != -1
            last[ next[ k + 1 ] + 1 ] = -1
        end
        head[ mindeg + 1 ] = next[ k + 1 ]
        elenk = elen[ k + 1 ]
        nvk = nv[ k + 1 ]
        nel += nvk


        if elenk > 0  && cnz + mindeg  >= nzmax
            for j =  0 : (n-1)
                if ( p = C.columnPointer[ j + 1 ] ) >= 0
                    C.columnPointer[ j + 1 ] = C.rowIndex[ p + 1 ]
                    C.rowIndex[ p + 1 ] = FLIP( j )
                end
            end
            q = 0 # change 0 to 1
            p = 0 # change 0 to 1
            while p < cnz
                tempp = p
                p += 1
                if (j = FLIP( C.rowIndex[ tempp + 1 ] ) ) >= 0
                    C.rowIndex[ q + 1 ] = C.columnPointer[ j + 1 ]
                    C.columnPointer[ j + 1 ] = q
                    q += 1
                    k3 = 0
                    for k3 = 0 : (len[ j + 1 ]- 2)
                        C.rowIndex[ q + 1 ] = C.rowIndex[ p + 1 ]
                        q += 1
                        p += 1
                    end
                end
            end
            cnz = q
        end

        dk = 0
        nv[ k + 1 ] = -nvk
        p = C.columnPointer[ k + 1 ]
        pk1 = cnz
        if elenk == 0
            pk1 = p
        end
        pk2 = pk1

        for k1 = 1 : (elenk + 1)
            if k1 > elenk
                e = k
                pj = p
                ln = len[ k + 1 ] - elenk
            else
                e = C.rowIndex[ p + 1 ]
                p += 1
                pj = C.columnPointer[ e + 1 ]
                ln = len[ e + 1 ]
            end

            for k2 = 1 : ln
                i = C.rowIndex[ pj + 1 ]
                pj += 1
                if ( nvi = nv[ i + 1 ] ) <= 0
                    continue
                end
                dk += nvi
                nv[ i + 1 ] = -nvi
                C.rowIndex[ pk2 + 1 ] = i
                pk2 += 1
                if next[ i + 1 ] != -1
                    last[ next[ i + 1 ] + 1 ] = last[ i + 1 ]
                end
                if last[ i + 1 ] != -1
                    next[ last[ i + 1 ] + 1 ] = next[ i + 1 ]
                else
                    head[ degree[ i + 1 ] + 1 ] = next[ i + 1 ]
                end

            end
            if e != k
                C.columnPointer[ e + 1 ] = FLIP( k )
                w[ e + 1] = 0
            end

        end

        if elenk != 0
            cnz = pk2
        end

        degree[ k + 1 ] = dk
        C.columnPointer[ k + 1 ] = pk1
        len[ k + 1 ] = pk2 - pk1
        elen[ k + 1 ] = -2
        mark = wclear( mark, lemax, w, n )

        for pk = pk1 : ( pk2 - 1 )
            i = C.rowIndex[ pk + 1 ]
            if (eln = elen[ i + 1  ]) <= 0
                continue
            end
            nvi = -nv[ i + 1 ]
            wnvi = mark - nvi

            for p = C.columnPointer[ i + 1  ] : ( C.columnPointer[ i + 1 ] + eln - 1 )
                e = C.rowIndex[ p + 1 ]
                if w[ e + 1 ] >= mark
                    w[ e + 1 ] -= nvi
                elseif w[ e + 1 ] != 0
                    w[ e + 1  ] = degree[ e + 1  ] + wnvi
                end
            end
        end


        for pk = pk1 : ( pk2 - 1 )
            i = C.rowIndex[ pk + 1 ]
            p1 = C.columnPointer[ i + 1 ]
            p2 = p1 + elen[ i + 1 ] - 1
            pn = p1
            h = 0
            d = 0

            for p = p1 : p2
                e = C.rowIndex[ p + 1  ]
                if w[ e + 1  ] != 0
                    dext = w[ e + 1 ] - mark
                    if dext > 0
                        d += dext
                        C.rowIndex[ pn + 1 ] = e
                        pn += 1
                        h += e
                    else
                        C.columnPointer[ e + 1 ] = FLIP( k )
                        w[ e + 1 ] = 0
                    end
                end
            end
            elen[ i + 1  ] = pn - p1  + 1
            p3 = pn
            p4 = p1 + len[ i + 1  ]

            for p = (p2 + 1) : (p4 -1 )
                j = C.rowIndex[ p + 1 ]
                if (nvj = nv[ j + 1 ] ) <= 0
                    continue
                end
                d += nvj
                C.rowIndex[ pn + 1 ] = j
                pn += 1
                h += j
            end

            if d == 0
                C.columnPointer[ i + 1 ] = FLIP( k ) # check
                nvi = -nv[ i + 1 ]
                dk -= nvi
                nvk += nvi
                nel += nvi
                nv[ i + 1  ] = 0
                elen[ i + 1 ] = -1
            else
                degree[ i + 1 ] = min( degree[ i + 1 ], d )
                C.rowIndex[ pn + 1 ] = C.rowIndex[ p3 + 1 ]
                C.rowIndex[ p3 + 1 ] = C.rowIndex[ p1 + 1 ]
                C.rowIndex[ p1 + 1 ] = k
                len[ i + 1 ] = pn - p1 + 1
                h = abs(h) % n
                next[ i + 1 ] = hhead[ h + 1 ]
                hhead[ h + 1 ] = i
                last[ i + 1 ] = h
            end

        end

        degree[ k + 1 ] = dk
        lemax = max( lemax, dk ) # TODO: check max
        mark = wclear( mark + lemax, lemax, w, n )
        # -- Supernode detection --- #
        for pk = pk1 : (pk2-1)
            i = C.rowIndex[ pk + 1 ]
            if nv[ i + 1 ] >= 0 # check
                continue
            end
            h = last[ i + 1 ]
            i = hhead[ h + 1  ]
            hhead[ h + 1  ] = -1

            while i != -1 && next[ i + 1 ] != -1
                ln = len[ i + 1]
                eln = elen[ i + 1]
                for p = (C.columnPointer[ i + 1 ] + 1):(C.columnPointer[ i + 1 ] + ln -1)
                    w[ C.rowIndex[ p + 1 ] + 1 ] = mark
                end

                jlast = i
                j = next[ i + 1 ]
                while j != -1
                    ok = ( len[ j + 1 ] == ln ) && ( elen[ j + 1 ] == eln )
                    p = C.columnPointer[ j + 1 ] + 1
                    while ok && p <= (C.columnPointer[ j + 1 ] + ln - 1)
                        if w[ C.rowIndex[ p + 1 ] + 1 ] != mark
                            ok = false
                        end
                        p +=1
                    end
                    if ok
                        C.columnPointer[ j + 1 ] = FLIP( i )
                        nv[ i + 1 ] += nv[ j + 1 ]
                        nv[ j + 1 ] = 0
                        elen[ j + 1 ] = -1
                        j = next[ j + 1]
                        next[ jlast + 1 ] = j
                    else
                        jlast = j
                        j = next[ j + 1 ]
                    end
                end
                i = next[ i + 1 ]
                mark += 1
            end
        end
        # finaize new element
        p = pk1

        for pk = pk1:(pk2-1)
            i = C.rowIndex[ pk + 1 ]
            if ( nvi =  -nv[ i + 1] ) <= 0
                continue
            end
            nv[ i + 1  ] = nvi
            d = degree[ i + 1 ] + dk - nvi
            d = min( d, n - nel - nvi ) # TODO: check min

            if head[ d + 1  ] != -1
                last[ head[ d + 1  ] + 1 ] = i
            end

            next[ i + 1 ] = head[ d + 1 ]
            last[ i + 1 ] = -1
            head[ d + 1 ] = i
            mindeg = min( mindeg, d )
            degree[ i + 1 ] = d
            C.rowIndex[ p + 1 ] = i
            p += 1
            pk += 1
        end
        nv[ k + 1 ] = nvk
        if (len[ k + 1 ] = (p - pk1) ) == 0
            C.columnPointer[ k + 1 ] = -1
            w[ k + 1 ] = 0
        end
        if elenk != 0
            cnz = p
        end
    end

    for i = 1 : n
        C.columnPointer[ i ] = FLIP( C.columnPointer[ i ] )
    end

    for j = 0 : n
        head[ j + 1 ] = -1
    end
    for j = n : -1 : 0
        if nv[ j + 1 ] > 0
            continue
        end
        next[ j + 1 ] = head[ C.columnPointer[ j + 1 ] + 1 ]
        head[ C.columnPointer[ j + 1 ] + 1 ] = j
    end
    for e = n : -1 : 0
        if nv[ e + 1 ] <= 0
            continue
        end
        if C.columnPointer[ e + 1 ] != -1
            next[ e + 1 ] = head[ C.columnPointer[ e + 1 ] + 1 ]
            head[ C.columnPointer[ e + 1 ] + 1 ] = e
        end
    end

    k = 1
    stack = zeros( Int, n )
    for i = 0 : n
        if C.columnPointer[ i + 1 ] == -1
            k = iterDFS0Indexed( i+1, k, head, next, P, n, stack)
        end
    end

    return P
end
