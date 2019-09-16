<<<<<<< HEAD
## Format

Variant Call Format (or *VCF*) files already follow a format standard, and while there is some wiggle-room for optional values, PopGen.jl only requires the core/mandatory components of a VCF, meaning problems should hopefully not arise regardless of which variant caller you are using (although we use `Freebayes` ourselves). Please open an issue if they do, or reach out to us on the community Slack.

!!! failure "filter VCF files beforehand"
    Keep in mind, VCF files need to be filtered **before** importing them into PopGen.jl. There is no and will be no VCF-filtering functionality to this package, as it is outside of the purpose of PopGen.jl. Refer to `vcftools` and `vcflib` to filter your sequence data. 



## Import a VCF file as a `PopObj`

PopGen.jl provides a simple command `vcf` to import a VCF file as a `PopObj`, which requires only the name of the file and nothing else.



```julia
potato = vcf("/home/data/russet_potatoes.vcf")
```



## What VCF files lack

Due to the nature of the file format, importing VCF files will provide:

- sample names
- ploidy of each sample
- locus names
- genotypes

but they will **not** provide:

- population information
- latitude or longitude



This means you will need to add that information separately afterwards. Location data (which is optional!) can be added to the `PopObj` directly with `.samples.latitude` or `.samples.longitude` or with the `locations!` command. Population names (mandatory!) can be added by overwriting `.samples.population` with an array of population names.  
=======
Next on th To-Do list!
>>>>>>> master
