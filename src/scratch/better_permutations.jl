function perm_samp(data::PopData)
    popcounts = countmap(data.metadata.sampleinfo.population)
    pops, counts = keys(popcounts), values(popcounts)
    gdf = groupby(data.genodata, :name)
    perm_idx = shuffle(Xoroshiro128Star(), keys(gdf))
    
    #return gdf

end