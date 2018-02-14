directory = "testresults/netlib/pwr/lu"
tuples = []

function processDirectory( dir )
    files = readdir( dir )
    for file_name in files
        if file_name != ".DS_Store"
            if contains( file_name, ".result")
                file = open(dir * "/" * file_name)
                lines = readlines( file )
                n = parse(Int64, lines[ 1 ])
                density = parse(Float64, lines[ 4 ])
                ccsparse_time = parse(Float64, lines[ 5 ])
                julia_time = parse(Float64, lines[ 6 ])
                ccsparse_mem = lines[ 7 ]
                julia_mem = lines[ 8 ]

                tuple = (n, density, ccsparse_time, julia_time, ccsparse_mem, julia_mem )
                push!( tuples, tuple )
                close(file)
            else
                new_dir = dir * "/" * file_name
                processDirectory( new_dir )
            end
        end
    end
end


processDirectory(directory)

sort!(tuples, by = x -> x[1])

for tuple in tuples
    @printf("%5d & %.5f & %.5f & %s & %.5f & %s\\\\ \\hline\n", tuple[ 1 ], tuple[ 2 ], tuple[ 3 ], tuple[ 5 ], tuple[ 4 ], tuple[ 6 ])
end
