# Modifying a genepop file for `Genepop` Isolation by Distance calculations

**Before we get into the nitty-gritty, allow me a moment to gripe:**

So, the (beloved) genepop file format specification exists to be used by the software (or R package) `Genepop`. That's all well and good (though honestly XML or JSON is probably a better format for SNP data with metadata), but if you want to perform Isolation by Distance calculations using Genepop, they ask that you reformat your genepop file for a completely wacky format, where you preface each sample name with x and y coordinates.

**Example:**

Normal sample row for genepop file:

```
thing_001, 001001 002001 001002 002002
```

With location data:

```
-43.111 12.221 thing_001, 001001 002001 001002 002002
```

The name `thing_001` is optional under this format (to my understanding), but it doesn't hurt to have it.

**But**

If you want to perform an individual-based model, you must also separate each sample as its own "population" using the word `POP`.

```
POP
-43.111 12.221 thing_001, 001001 002001 001002 002002
POP
-43.121 21.261 thing_002, 002001 002001 001002 002001
POP
-43.411 11.211 thing_002, 001002 001001 001002 001002
POP
-43.177 13.271 thing_002, 001001 002001 002001 002002
```



As you can imagine, the amount of work to do this correctly _quickly_ adds up if you are editing the file by hand. Modern SNP datasets have hundreds of individuals, and that reformatting is an absurd amount of work to do, and every time you manually edit something, it introduces the possibility of human error :frowning:.

Hopefully there will be a time when PopGen.jl will incorporate Isolation-by-Distance calculations (looking for volunteers :wink:) and you won't need to format your `PopObj` whatsoever, but until then we encourage you to continue using `Genepop` (available [here](https://cran.r-project.org/package=genepop)) or whatever your preferred method is. 

In the meantime, using the most basic parts of PopGen.jl and Julia, we can reformat a genepop file really quickly and relatively painlessly because most of the work is already done for you in the code snippets below. If you set up those first three variables `gfile`, `loc_data`, and `outfile` correctly and your coordinate data is in a two-column X/Y format, the rest of this code _should_ work to generate your modified Isolation-by-Distance ready genepop file. Make sure your `outfile` doesn't exist yet (or is empty) or else the output will append to whatever is in that file.

### Conversion for by-individual

Along with prefixing the sample names with the location data, this snippet will force a `POP` to separate individuals are shown above.

```julia
using PopGen, CSV

gfile = "PATH/TO/YOUR/GENEPOP/FILE"
loc_data = "PATH/TO/YOUR/LAT-LONG/FILE"
outfile = "PATH/AND/NAME/OF/INTENDED/OUTPUT/FILE"

#load in your genepop as a PopObj
gpop_file = genepop(gfile)
# load in your location data as a table
loc_df = DataFrame(CSV.File(loc_data))

oldnames = sample_names(gpop_file)  # this is the only reason we need PopGen.jl here
newnames = Vector{String}()

# merge the x,y coordinates with the original sample name
# make sure your longitude (x) and latitude (y) are in the correct order
# this loop generates the format "x y sample_name" as a string
for (i,j,k) in zip(loc_df[!,1], loc_df[!,2],oldnames)
  push!(newnames, string("$i\t$j\t$k"))
end

# open the original genepop file and parse each line as a string 
popgenfile = open(readlines, gfile)

# Create the output file (if it doesn't already exist) and parse the original genepop file to replace the sample names with the new names we created above and print to outfile
open(outfile, "w") do f
  name_num = 1
  len = length(popgenfile)
  for line in 1:len
    if occursin("POP", popgenfile[line]) == true
      # A little condition to add the first POP to the file
      if name_num == 1
        println(f, "POP")
      end
      # skip POP lines
      continue
    elseif occursin(",",  popgenfile[line]) != true
      # output locus names
      println(f,  popgenfile[line])
      continue
    else
      eachline_mod = replace(popgenfile[line], oldnames[name_num] => newnames[name_num])
        if line != len
          println(f, eachline_mod, "\nPOP")
          name_num += 1
        else
          print(f, "$eachline_mod")
        end
    end
  end
end
```

### Conversion for by-population

This snippet will prefix the location information to each sample name, but preserve the original POP tags.

```julia
using PopGen, CSV

gfile = "PATH/TO/YOUR/GENEPOP/FILE"
loc_data = "PATH/TO/YOUR/LAT-LONG/FILE"
outfile = "PATH/AND/NAME/OF/INTENDED/OUTPUT/FILE"

#load in your genepop as a PopObj
gpop_file = genepop(gfile)
# load in your location data as a table
loc_df = DataFrame(CSV.File(loc_data))

oldnames = sample_names(gpop_file) # this is the only reason we need PopGen.jl here
newnames = Vector{String}()

# merge the x,y coordinates with the original sample name
# make sure your longitude (x) and latitude (y) are in the correct order
# this loop generates the format "x y sample_name" as a string
for (i,j,k) in zip(loc_df[!,1], loc_df[!,2],oldnames)
  push!(newnames, string("$i\t$j\t$k"))
end

# open the original genepop file and parse each line as a string 
popgenfile = open(readlines, gfile)

# Create the output file (if it doesn't already exist) and parse the original genepop file to replace the sample names with the new names we created above  and print to outfile
open(outfile, "w") do f
  name_num = 1
  len = length(popgenfile)
  for line in 1:len
    if occursin("POP", popgenfile[line]) == true
      println(f, "POP")
      continue
    elseif occursin(",",  popgenfile[line]) != true
      # output locus names
      println(f,  popgenfile[line])
      continue
    else
      eachline_mod = replace(popgenfile[line], oldnames[name_num] => newnames[name_num])
      println(f, eachline_mod)
      name_num += 1
    end
  end
end
```

