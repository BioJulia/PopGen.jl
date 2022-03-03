module.exports = {
  "docs": {
    "Getting Started": [
      "gettingstarted/install",
      "gettingstarted/juliaprimer",
      "gettingstarted/tips",
      "gettingstarted/comparison",
      "gettingstarted/datasets"
    ],
    "File I/O": [
      "io/readingdata",
      "io/writingdata",
      {
        "File Formats": [
          "io/delimited",
          "io/genepop",
          "io/structure",
          "io/vcf"
        ]
      },
    ],
    "Working with PopData": [
      "popdata/popdata",
      "popdata/workingwithpopdata",
      "popdata/viewdata",
      "popdata/addingdata",
      "popdata/exclusion",
      "popdata/conditionals",
      "popdata/populationdata",
      "popdata/locationdata",
      "popdata/dataexploration",
      "popdata/advancedindexing"
    ],
    "Analyses":[
      "analyses/hardyweinberg",
      "analyses/relatedness",
      "analyses/fstatistics"
    ],
    "Simulations":[
      "simulations/simulate_samples",
      "simulations/breedingcrosses",
      "simulations/sibship_simulations",
    ],    
    "API": [
      "api/api",
      {
        "PopGen":[
          "api/PopGen/dataexploration",
          "api/PopGen/fstbylocus",
          "api/PopGen/fstglobal",
          "api/PopGen/fstpermutations",
          "api/PopGen/hardyweinberg",
          "api/PopGen/heterozygosity",
          "api/PopGen/pairwisefst",
          "api/PopGen/pairwisekinship",
          "api/PopGen/kinshipmoments",
          "api/PopGen/kinshipposthocs",
          "api/PopGen/summaryinfo",
          "api/PopGen/utils",
        ]
      },
      {
        "PopGenCore":[
          "api/PopGenCore/allelefreq",
          "api/PopGenCore/conditionals",
          "api/PopGenCore/datasets",
          "api/PopGenCore/delimited",
          "api/PopGenCore/genepop",
          "api/PopGenCore/generalutils",
          "api/PopGenCore/genofreq",
          "api/PopGenCore/genotypeutils",
          "api/PopGenCore/ioutils",
          "api/PopGenCore/iterators",
          "api/PopGenCore/manipulate",
          "api/PopGenCore/mathutils",
          "api/PopGenCore/missingutils",
          "api/PopGenCore/permutations",
          "api/PopGenCore/popdatawrappers",
          "api/PopGenCore/read",
          "api/PopGenCore/structure",
          "api/PopGenCore/types",
          "api/PopGenCore/variantcall",
        ]
      },
      {
      "PopGenSims":[
        "api/PopGenSims/popgensims_cross",
        "api/PopGenSims/popgensims_samples",
        "api/PopGenSims/popgensims_sibship",
        "api/PopGenSims/popgensims_utils",
      ]
      }
    ]
  }
}
