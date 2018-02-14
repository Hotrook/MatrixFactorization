using Gadfly

file = "test0.0001-0.01"
algorithm = "cholesky"
param = "time"
path = "testresults/" * algorithm * "/" * file * "/"

first_description = ""
second_description = "Algorytm z języka Julia\n dla macierzy rzadkich"
third_description = "Algorytm z języka Julia\n dla macierzy gęstych"

if algorithm == "cholesky"
    first_description = "Algortym up-looking przy\n użyciu struktury CCSparseMatrix"
else
    first_description = "Algortym left-looking przy\n użyciu struktury CCSparseMatrix"
end

enumeration_file_name = path * "enumeration.txt"
symbolic_results_file_name = path * "symbolic_results.txt"
ccsparse_results_file_name = path * "ccsparse_results.txt"
julia_sparse_results_file_name = path * "julia_sparse_results.txt"
julia_results_file_name = path * "julia_results.txt"

type PlotData
    time::Array{Float64}
    mem::Array{Int64}

    function PlotData(n)
        new(zeros(Float64, n), zeros(Int64, n))
    end
end

function readData(filename, n::Int64)
    result = PlotData(n)
    file = open(filename)
    lines = readlines(file)
    result.time = map( x -> parse(Float64, split(x)[1]), lines)
    result.mem = map( x -> parse(Int64, split(x)[2]), lines)
    close(file)
    result
end

function readEnumeration()
    file = open( enumeration_file_name )
    lines = readlines( file )
    close(file)
    map( x -> parse( Float64, x), lines )
end


println("test")
enumeration = readEnumeration()
n = length(enumeration)
scope = n

symbolic_data = readData(symbolic_results_file_name, n)
julia_sparse_data = readData(julia_sparse_results_file_name, n)
julia_data = readData(julia_results_file_name, n)
ccsparse_data = readData(ccsparse_results_file_name, n)

ccsparse_layer = julia_sparse_layer = julia_layer = layer()
ylabel = ""
if param == "time"
    ccsparse_layer = layer(x = enumeration[1:scope], y =  symbolic_data.time[1:scope] + ccsparse_data.time[1:scope], Geom.line, Theme(default_color=color("green")))
    julia_sparse_layer = layer(x = enumeration[1:scope], y = julia_sparse_data.time[1:scope], Geom.line, Theme(default_color=color("red")))
    julia_layer = layer(x = enumeration[1:scope], y = julia_data.time[1:scope], Geom.line, Theme(default_color=color("purple")))
    ylabel = "Czas trwania [s]"

else
    ccsparse_layer = layer(x = enumeration[1:scope], y =  symbolic_data.mem[1:scope] + ccsparse_data.mem[1:scope], Geom.line, Theme(default_color=color("green")))
    julia_sparse_layer = layer(x = enumeration[1:scope], y = julia_sparse_data.mem[1:scope], Geom.line, Theme(default_color=color("red")))
    julia_layer = layer(x = enumeration[1:scope], y = julia_data.mem[1:scope], Geom.line, Theme(default_color=color("purple")))
    ylabel = "Zużycie pamięci [b]"
end

 p = plot(
    ccsparse_layer,
    julia_sparse_layer,
    julia_layer,
    Guide.manual_color_key(
        "Algorytm",
        [first_description, second_description, third_description],
        ["green", "red", "purple"]),
    Theme(
        background_color="white",
        key_label_font_size=14pt,
        minor_label_font_size=14pt,
        major_label_font_size=14pt,
        ),
    Guide.xlabel("Gęstość macierzy"),
    Guide.ylabel(ylabel)
    )

    img = PNG(file * param * ".png", 12inch, 6inch)
    draw(img, p)
