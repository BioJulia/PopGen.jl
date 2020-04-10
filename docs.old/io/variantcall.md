## Import a BCF/VCF file as a `PopObj`

PopGen.jl provides a simple command `vcf` to import a VCF file as a `PopObj`, which requires only the name of the file and nothing else. While not strictly necessary, if you have polyploid or mixed-ploidy samples, this import method may be most efficient.

```julia
potato = vcf("/home/data/russet_potatoes.vcf")
```

and similar for `bcf` formatted files:

```julia
cabbage = bcf("/home/data/nappa_cabbage.bcf")
```

??? warning "Windows users"
    make sure to change your backslashes "\" to forward slashes "/" 

## What BCF/VCF files contain and lack

Due to the nature of the file format, importing variant call files **will** provide:

- sample names
- ploidy of each sample
- locus names
- genotypes

but they **will not** provide:

- population information
- latitude or longitude



This means you will need to add that information separately afterwards. Location data (which is optional!) can be added to the `PopObj` directly to `.samples.latitude` or `.samples.longitude` or with the `locations!` command. Population names (mandatory!) can be added using `populations!()`

## Format

Variant Call Format (or *VCF*) files already follow a format standard, and while there is some wiggle-room for optional values, PopGen.jl only requires the core/mandatory components of a VCF, meaning problems should hopefully not arise regardless of which variant caller you are using (although we use `Freebayes` ourselves). Please open an issue if they do, or reach out to us on the community Slack.

!!! failure ""
    Keep in mind, BCF/VCF files need to be filtered **before** importing them into PopGen.jl. There is no and will be no VCF-filtering functionality to this package, as it is outside of the purpose of PopGen.jl. Refer to `vcftools` and `vcflib` to filter your sequence data. 



## Acknowledgements

The majority of the BCF/VCF reader is thanks to the tremendous efforts of Ben J. Ward and the BioJulia contributors involved in [GeneticVariation.jl](https://github.com/BioJulia/GeneticVariation.jl), which we use to parse your files into `PopObj` format. If you love the file importer, then give those folks your thanks. If something is wrong and/or you hate the importer, blame us (and please [open up an issue](https://github.com/pdimens/PopGen.jl/issues) :sweat_smile:).