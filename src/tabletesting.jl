using JuliaDB, BenchmarkTools, CategoricalArrays, DataFrames, WeakRefStrings, Query

function tabletest()
    names = StringArray{WeakRefString}(fill.(["red","blue","green"], 150000) |> Base.Iterators.flatten |> collect)
    loci = "loc".*StringArray{WeakRefString}(string.(collect(1:1500)))
    loci = fill.(loci, 300) |> Base.Iterators.flatten |> collect
    genotypes = fill((1,2), 450000)

    return table((names = names, loci = loci, genotypes = genotypes), pkey = :names)
end

function df_test_ws()
    names = StringArray{WeakRefString}(fill.(["red","blue","green"], 150000) |> Base.Iterators.flatten |> collect)
    loci = "loc".*StringArray{WeakRefString}(string.(collect(1:1500)))
    loci = fill.(loci, 300) |> Base.Iterators.flatten |> collect
    genotypes = fill((1,2), 450000)

    DataFrame([names, loci, genotypes])
end

@benchmark tabletest()

function table_table()
    names = StringArray{WeakRefString}(fill.(["red","blue","green"], 150000) |> Base.Iterators.flatten |> collect)
    loci = "loc".*StringArray{WeakRefString}(string.(collect(1:1500)))
    loci = fill.(loci, 300) |> Base.Iterators.flatten |> collect
    genotypes = fill((1,2), 450000)
    
    JuliaDB.table(Tables.columntable((names, loci, genotypes)))
end


====================
function dftest()
    names = fill.(["red","blue","green"], 150000) |> Base.Iterators.flatten |> collect
    loci = "loc".*string.(collect(1:1500))
    loci = fill.(loci, 300) |> Base.Iterators.flatten |> collect
    genotypes = fill((1,2), 450000)

    return DataFrame([names, loci, genotypes])
end

function dftest_cat()
    names = categorical(fill.(["red","blue","green"], 150000) |> Base.Iterators.flatten |> collect, compress = true)
    loci = "loc".*string.(collect(1:1500))
    loci = categorical(fill.(loci, 300) |> Base.Iterators.flatten |> collect, compress = true)
    genotypes = fill((1,2), 450000) |> categorical

    return DataFrame([names, loci, genotypes])
end

function dftest_weak()
    names = StringArray{WeakRefString}(fill.(["red","blue","green"], 150000) |> Base.Iterators.flatten |> collect)
    loci = "loc".*string.(collect(1:1500))
    loci = StringArray{WeakRefString}(fill.(loci, 300) |> Base.Iterators.flatten |> collect)
    genotypes = fill((1,2), 450000) |> categorical

    return DataFrame([names, loci, genotypes])
end

@benchmark  dftest()

function df_tbl_test()
    names = fill.(["red","blue","green"], 150000) |> Base.Iterators.flatten |> collect
    loci = "loc".*string.(collect(1:150000))
    loci = fill.(loci, 3) |> Base.Iterators.flatten |> collect
    genotypes = fill((1,2), 450000)

    df = DataFrame([names, loci, genotypes])
    return table((names = df.x1, loci = df.x2, genotypes = df.x3), pkey = :names)
end

@btime df_tbl_test()