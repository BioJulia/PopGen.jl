"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[9369],{3905:function(e,n,t){t.d(n,{Zo:function(){return d},kt:function(){return k}});var l=t(7294);function a(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function i(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);n&&(l=l.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,l)}return t}function r(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?i(Object(t),!0).forEach((function(n){a(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):i(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function o(e,n){if(null==e)return{};var t,l,a=function(e,n){if(null==e)return{};var t,l,a={},i=Object.keys(e);for(l=0;l<i.length;l++)t=i[l],n.indexOf(t)>=0||(a[t]=e[t]);return a}(e,n);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(l=0;l<i.length;l++)t=i[l],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(a[t]=e[t])}return a}var p=l.createContext({}),m=function(e){var n=l.useContext(p),t=n;return e&&(t="function"==typeof e?e(n):r(r({},n),e)),t},d=function(e){var n=m(e.components);return l.createElement(p.Provider,{value:n},e.children)},u="mdxType",s={inlineCode:"code",wrapper:function(e){var n=e.children;return l.createElement(l.Fragment,{},n)}},c=l.forwardRef((function(e,n){var t=e.components,a=e.mdxType,i=e.originalType,p=e.parentName,d=o(e,["components","mdxType","originalType","parentName"]),u=m(t),c=a,k=u["".concat(p,".").concat(c)]||u[c]||s[c]||i;return t?l.createElement(k,r(r({ref:n},d),{},{components:t})):l.createElement(k,r({ref:n},d))}));function k(e,n){var t=arguments,a=n&&n.mdxType;if("string"==typeof e||a){var i=t.length,r=new Array(i);r[0]=c;var o={};for(var p in n)hasOwnProperty.call(n,p)&&(o[p]=n[p]);o.originalType=e,o[u]="string"==typeof e?e:a,r[1]=o;for(var m=2;m<i;m++)r[m]=t[m];return l.createElement.apply(null,r)}return l.createElement.apply(null,t)}c.displayName="MDXCreateElement"},7450:function(e,n,t){t.r(n),t.d(n,{assets:function(){return d},contentTitle:function(){return p},default:function(){return k},frontMatter:function(){return o},metadata:function(){return m},toc:function(){return u}});var l=t(7462),a=t(3366),i=(t(7294),t(3905)),r=["components"],o={id:"plink",title:"Plink.jl",sidebar_label:"Plink.jl"},p=void 0,m={unversionedId:"api/PopGenCore/plink",id:"api/PopGenCore/plink",title:"Plink.jl",description:"PopGenCore.jl/src/io/Plink.jl",source:"@site/docs/api/PopGenCore/Plink.md",sourceDirName:"api/PopGenCore",slug:"/api/PopGenCore/plink",permalink:"/PopGen.jl/docs/api/PopGenCore/plink",draft:!1,editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGenCore/Plink.md",tags:[],version:"current",lastUpdatedAt:1658771861,formattedLastUpdatedAt:"Jul 25, 2022",frontMatter:{id:"plink",title:"Plink.jl",sidebar_label:"Plink.jl"},sidebar:"docs",previous:{title:"Permutations.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/permutations"},next:{title:"PopDataWrappers.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/popdatawrappers"}},d={},u=[{value:"PopGenCore.jl/src/io/Plink.jl",id:"popgencorejlsrcioplinkjl",level:2},{value:"\ud83d\udce6 _plinkindex",id:"-_plinkindex",level:3},{value:"\ud83d\udce6 _SNP",id:"-_snp",level:3},{value:"\ud83d\udce6 _plinkped",id:"-_plinkped",level:3},{value:"\ud83d\udce6 _plinkbed",id:"-_plinkbed",level:3},{value:"\ud83d\udce6 _genoconversion",id:"-_genoconversion",level:3},{value:"\ud83d\udfea\ud83d\udd35 plink",id:"-plink",level:3},{value:"Keyword Arguments",id:"keyword-arguments",level:3},{value:"Example",id:"example",level:2},{value:"Example",id:"example-1",level:2}],s={toc:u},c="wrapper";function k(e){var n=e.components,t=(0,a.Z)(e,r);return(0,i.kt)(c,(0,l.Z)({},s,t,{components:n,mdxType:"MDXLayout"}),(0,i.kt)("h2",{id:"popgencorejlsrcioplinkjl"},"PopGenCore.jl/src/io/Plink.jl"),(0,i.kt)("table",null,(0,i.kt)("thead",{parentName:"table"},(0,i.kt)("tr",{parentName:"thead"},(0,i.kt)("th",{parentName:"tr",align:"center"},"\ud83d\udce6  not exported"),(0,i.kt)("th",{parentName:"tr",align:"center"},"\ud83d\udfea  exported by PopGenCore.jl"),(0,i.kt)("th",{parentName:"tr",align:"center"},"\ud83d\udd35  exported by PopGen.jl")))),(0,i.kt)("h3",{id:"-_plinkindex"},"\ud83d\udce6 _plinkindex"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"_plinkindex(s::Matrix{UInt8}, i::Integer, j::Integer)\n_plinkindex(s::Matrix{UInt8})\n")),(0,i.kt)("p",null,"Tthis function copies the getindex tool OpenMendel/SnpArrays.jl uses\nto pull out the byte values from the compressed hex genotypes\nRepo: ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/OpenMendel/SnpArrays.jl"},"https://github.com/OpenMendel/SnpArrays.jl")),(0,i.kt)("hr",null),(0,i.kt)("h3",{id:"-_snp"},"\ud83d\udce6 _SNP"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"_SNP(genotype::UInt8)\n_SNP(genomatrix::AbstractArray{UInt8})\n")),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"00\tHomozygous for first allele (0x00)"),(0,i.kt)("li",{parentName:"ul"},"01\tMissing genotype (0x01)"),(0,i.kt)("li",{parentName:"ul"},"10\tHeterozygous  (0x02)"),(0,i.kt)("li",{parentName:"ul"},"11\tHomozygous for second allele in .bim file (0x03)")),(0,i.kt)("hr",null),(0,i.kt)("h3",{id:"-_plinkped"},"\ud83d\udce6 _plinkped"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"_plinkped(infile::String, keepfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)\n")),(0,i.kt)("hr",null),(0,i.kt)("h3",{id:"-_plinkbed"},"\ud83d\udce6 _plinkbed"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"_plinkbed(infile::String, famfields::Union{Symbol,Vector{Symbol}} = :all, bimfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)\n")),(0,i.kt)("hr",null),(0,i.kt)("h3",{id:"-_genoconversion"},"\ud83d\udce6 _genoconversion"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'_genoconversion(genotype::T) where T<:Genotype = join(genotype, " ")\n_genoconversion(genotype::Missing) = "0 0"\n')),(0,i.kt)("hr",null),(0,i.kt)("h3",{id:"-plink"},"\ud83d\udfea\ud83d\udd35 plink"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"plink(infile::String; famfields::Union{Symbol,Vector{Symbol}} = :all, bimfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)\n")),(0,i.kt)("p",null,"Read a PLINK ",(0,i.kt)("inlineCode",{parentName:"p"},".ped")," or binary ",(0,i.kt)("inlineCode",{parentName:"p"},".bed")," file into memory as a ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," object.\nRequires an accompanying ",(0,i.kt)("inlineCode",{parentName:"p"},".fam")," file in the same directory, but an accompanying ",(0,i.kt)("inlineCode",{parentName:"p"},".bim")," file is optional."),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"infile::String")," : path to ",(0,i.kt)("inlineCode",{parentName:"li"},".ped")," or ",(0,i.kt)("inlineCode",{parentName:"li"},".bed")," file")),(0,i.kt)("h3",{id:"keyword-arguments"},"Keyword Arguments"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"famfields::Symbol|Vector{Symbol}")," : which additional fields to import from the ",(0,i.kt)("inlineCode",{parentName:"li"},".fam")," file",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},":all")," ","[default]"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},":none")),(0,i.kt)("li",{parentName:"ul"},"any one or combination of ",(0,i.kt)("inlineCode",{parentName:"li"},"[:sire, :dam, :sex, :phenotype]")))),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"bimfields::Symbol|Vector{Symbol}")," : which additional fields to import from the optional ",(0,i.kt)("inlineCode",{parentName:"li"},".bim")," file",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},":all")," ","[default]"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},":none")),(0,i.kt)("li",{parentName:"ul"},"any one or combination of ",(0,i.kt)("inlineCode",{parentName:"li"},"[:chromosome, :cm, :bp]")))),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"silent::Bool"),"   : whether to print file information during import (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"false"),")")),(0,i.kt)("h2",{id:"example"},"Example"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'parakeet = plink("datadir/parakeet.ped", famfields = :sex)\nparrot = plink("datadir/parrot.bed", famfields = [:sire, :dam], bimfields = :chromosome)\n')),(0,i.kt)("hr",null),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"plink(data::PopData; filename::String)\n")),(0,i.kt)("p",null,"Write a biallelic ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," object to PLINK ",(0,i.kt)("inlineCode",{parentName:"p"},".ped")," format with an accompanying\n",(0,i.kt)("inlineCode",{parentName:"p"},".fam")," file. Genotypes are coded by the PLINK standard:"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Integers are the alleles"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"0")," encodes missing"),(0,i.kt)("li",{parentName:"ul"},"After column 6, every two numbers indicate a diploid genotype.")),(0,i.kt)("h2",{id:"example-1"},"Example"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'sharks = dropmultiallelic(@gulfsharks) ;\nplink(sharks, filename = "biallelic_sharks.ped")\n')))}k.isMDXComponent=!0}}]);