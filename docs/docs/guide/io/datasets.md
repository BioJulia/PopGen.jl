# Provided datasets
PopGen.jl provides two datasets as examples, `nancycats` and `gulfsharks`. The datasets can be retrieved using the `dataset` function, or their names as commands without arguments (e.g. `gulfsharks()`). 

::: tip identitcal methods
The methods are identical (one is a wrapper for the other), but the benefit of calling the datasets directly by name is that you get the luxury of tab auto-completion :grin:
:::

## datasets
```julia
dataset(::String)
```
Returns a `PopData` object of the dataset you would like to retrieve by calling the dataset as a string by name.

**Example:**
```julia
sharks = dataset("gulfsharks")
cats = dataset("nancycats")
```
### nancycats

We include the familiar nancycats microsatellite data, as featured in [adegenet](http://adegenet.r-forge.r-project.org), for easy importing into PopGen.jl as `PopData`. As an alternative to `datasets`, you can invoke the `nancycats()`  command without any arguments.

```
julia> ncats = nancycats() ; summary(ncats)
PopData Object
  Marker type: Microsatellite
  Ploidy: 2
  Number of individuals: 237
  Number of loci: 9
  Populations: 17
  Longitude: absent
  Latitude: absent
```

The spatial coordinates provided for the dataset in `adegenet` are completely unfamiliar to us (and some geospatial folks we spoke to), so they have been omitted.  If you recognize what coordinate system has 485.111 appear in Nancy, France, please let us know!

### gulfsharks

We also include the SNP dataset used in Dimens *et al.* 2019 "[A genomic assessment of movement and gene flow around the South Florida vicariance zone in the migratory coastal blacknose shark, *Carcharhinus acronotus*](https://link.springer.com/article/10.1007/s00227-019-3533-1)" since it was already on hand. Like `nancycats`, we provide a convenient function to load these data into PopGen.jl as `PopData`. As an alternative to `datasets`, you can invoke the `gulfsharks()` command without any arguments. 

```jullia
julia> sharks = gulfsharks() ; summary(sharks)
PopData Object
  Marker type: SNP
  Ploidy: 2
  Number of individuals: 212
  Number of loci: 2213
  Populations: 7
  Longitude: present with 0 missing
  Latitude: present with 0 missing
```

