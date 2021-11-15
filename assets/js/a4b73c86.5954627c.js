"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[3962],{3905:function(e,t,a){a.d(t,{Zo:function(){return m},kt:function(){return u}});var n=a(7294);function i(e,t,a){return t in e?Object.defineProperty(e,t,{value:a,enumerable:!0,configurable:!0,writable:!0}):e[t]=a,e}function o(e,t){var a=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),a.push.apply(a,n)}return a}function r(e){for(var t=1;t<arguments.length;t++){var a=null!=arguments[t]?arguments[t]:{};t%2?o(Object(a),!0).forEach((function(t){i(e,t,a[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(a)):o(Object(a)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(a,t))}))}return e}function l(e,t){if(null==e)return{};var a,n,i=function(e,t){if(null==e)return{};var a,n,i={},o=Object.keys(e);for(n=0;n<o.length;n++)a=o[n],t.indexOf(a)>=0||(i[a]=e[a]);return i}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(n=0;n<o.length;n++)a=o[n],t.indexOf(a)>=0||Object.prototype.propertyIsEnumerable.call(e,a)&&(i[a]=e[a])}return i}var p=n.createContext({}),s=function(e){var t=n.useContext(p),a=t;return e&&(a="function"==typeof e?e(t):r(r({},t),e)),a},m=function(e){var t=s(e.components);return n.createElement(p.Provider,{value:t},e.children)},d={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},c=n.forwardRef((function(e,t){var a=e.components,i=e.mdxType,o=e.originalType,p=e.parentName,m=l(e,["components","mdxType","originalType","parentName"]),c=s(a),u=i,f=c["".concat(p,".").concat(u)]||c[u]||d[u]||o;return a?n.createElement(f,r(r({ref:t},m),{},{components:a})):n.createElement(f,r({ref:t},m))}));function u(e,t){var a=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var o=a.length,r=new Array(o);r[0]=c;var l={};for(var p in t)hasOwnProperty.call(t,p)&&(l[p]=t[p]);l.originalType=e,l.mdxType="string"==typeof e?e:i,r[1]=l;for(var s=2;s<o;s++)r[s]=a[s];return n.createElement.apply(null,r)}return n.createElement.apply(null,a)}c.displayName="MDXCreateElement"},8211:function(e,t,a){a.r(t),a.d(t,{frontMatter:function(){return l},contentTitle:function(){return p},metadata:function(){return s},toc:function(){return m},default:function(){return c}});var n=a(7462),i=a(3366),o=(a(7294),a(3905)),r=["components"],l={id:"vcf",title:"Variant Call Format",sidebar_label:"Variant Call Format"},p=void 0,s={unversionedId:"io/vcf",id:"io/vcf",isDocsHomePage:!1,title:"Variant Call Format",description:"Import a BCF/VCF file as PopData",source:"@site/docs/io/variantcall.md",sourceDirName:"io",slug:"/io/vcf",permalink:"/PopGen.jl/docs/io/vcf",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/io/variantcall.md",tags:[],version:"current",lastUpdatedAt:1635451805,formattedLastUpdatedAt:"10/28/2021",frontMatter:{id:"vcf",title:"Variant Call Format",sidebar_label:"Variant Call Format"},sidebar:"docs",previous:{title:"Structure",permalink:"/PopGen.jl/docs/io/structure"},next:{title:"The PopData type",permalink:"/PopGen.jl/docs/workingwithpopdata/popdata"}},m=[{value:"Import a BCF/VCF file as <code>PopData</code>",id:"import-a-bcfvcf-file-as-popdata",children:[{value:"Arguments",id:"arguments",children:[],level:3},{value:"Keyword Arguments",id:"keyword-arguments",children:[],level:3},{value:"Example",id:"example",children:[],level:3},{value:"Mixed-Ploidy data",id:"mixed-ploidy-data",children:[],level:3},{value:"Format",id:"format",children:[],level:3},{value:"Allele encodings",id:"allele-encodings",children:[],level:3},{value:"What BCF/VCF files contain",id:"what-bcfvcf-files-contain",children:[],level:3},{value:"What BCF/VCF files lack",id:"what-bcfvcf-files-lack",children:[],level:3}],level:2},{value:"Acknowledgements",id:"acknowledgements",children:[],level:2}],d={toc:m};function c(e){var t=e.components,a=(0,i.Z)(e,r);return(0,o.kt)("wrapper",(0,n.Z)({},d,a,{components:t,mdxType:"MDXLayout"}),(0,o.kt)("h2",{id:"import-a-bcfvcf-file-as-popdata"},"Import a BCF/VCF file as ",(0,o.kt)("inlineCode",{parentName:"h2"},"PopData")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"vcf(infile::String; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)\nbcf(infile::String; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)\n")),(0,o.kt)("p",null,"PopGen.jl provides the commands ",(0,o.kt)("inlineCode",{parentName:"p"},"vcf")," and ",(0,o.kt)("inlineCode",{parentName:"p"},"bcf")," to import a variant call format files into ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData"),". The reader also accepts files that are gzipped. "),(0,o.kt)("h3",{id:"arguments"},"Arguments"),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"infile::String")," : path to file, in quotes. ",(0,o.kt)("strong",{parentName:"li"},"must end in ",(0,o.kt)("inlineCode",{parentName:"strong"},".gz")," if gzipped"))),(0,o.kt)("h3",{id:"keyword-arguments"},"Keyword Arguments"),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"rename_loci::Bool"),": whether to simplify loci names to ",(0,o.kt)("inlineCode",{parentName:"li"},"snp_#")," (default: ",(0,o.kt)("inlineCode",{parentName:"li"},"false"),")"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"allow_monomorphic::Bool")," : whether to keep monomorphic loci (default: ",(0,o.kt)("inlineCode",{parentName:"li"},"false"),")"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"silent::Bool")," : whether to print file information during import (default: ",(0,o.kt)("inlineCode",{parentName:"li"},"false"),")")),(0,o.kt)("h3",{id:"example"},"Example"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'cabbage = bcf("/home/data/nappa_cabbage.bcf", rename_loci = true, silent = true)\npotato = vcf("/home/data/russet_potatoes.vcf.gz", allow_monomorphic = true)\n')),(0,o.kt)("h3",{id:"mixed-ploidy-data"},"Mixed-Ploidy data"),(0,o.kt)("p",null,"In the event your variant call file is for mixed-ploidy data (where ploidy is not the same across all samples, e.g. PoolSeq), you will need to perform an additional step after reading in your data as ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," to convert the ",(0,o.kt)("inlineCode",{parentName:"p"},".genodata.genotype")," column into a ",(0,o.kt)("inlineCode",{parentName:"p"},"GenoArray"),":"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'julia> mydata = bcf("path/to/file.bcf", silent = true, rename_loci = true) ;\n\njulia> mydata.genodata.genotype =  mydata.genodata.genotype |> Array{Union{Missing, NTuple}}\n')),(0,o.kt)("div",{className:"admonition admonition-caution alert alert--warning"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"16",height:"16",viewBox:"0 0 16 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M8.893 1.5c-.183-.31-.52-.5-.887-.5s-.703.19-.886.5L.138 13.499a.98.98 0 0 0 0 1.001c.193.31.53.501.886.501h13.964c.367 0 .704-.19.877-.5a1.03 1.03 0 0 0 .01-1.002L8.893 1.5zm.133 11.497H6.987v-2.003h2.039v2.003zm0-3.004H6.987V5.987h2.039v4.006z"}))),"WIP")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},"The extra step required by mixed-ploidy data is a work in progress. Feel free to submit a PR if you have ideas!"))),(0,o.kt)("h3",{id:"format"},"Format"),(0,o.kt)("p",null,"Variant Call Format files follow a format standard, and while there is some wiggle-room for optional values, PopGen.jl only requires the core/mandatory components of a BCF/VCF, meaning problems should hopefully not arise regardless of which variant caller you are using (although we use ",(0,o.kt)("inlineCode",{parentName:"p"},"Freebayes")," ourselves). Please open an issue if they do, or reach out to us on the community Slack."),(0,o.kt)("h3",{id:"allele-encodings"},"Allele encodings"),(0,o.kt)("p",null,"When converting to ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData"),", the nucleotides will be recoded according to the table below. Note that this system differs slightly from\nhow PGDSpider2 recodes alleles (the 3 and 4 are switched)."),(0,o.kt)("table",null,(0,o.kt)("thead",{parentName:"table"},(0,o.kt)("tr",{parentName:"thead"},(0,o.kt)("th",{parentName:"tr",align:"left"},"Base"),(0,o.kt)("th",{parentName:"tr",align:"center"},"A"),(0,o.kt)("th",{parentName:"tr",align:"center"},"T"),(0,o.kt)("th",{parentName:"tr",align:"center"},"C"),(0,o.kt)("th",{parentName:"tr",align:"center"},"G"))),(0,o.kt)("tbody",{parentName:"table"},(0,o.kt)("tr",{parentName:"tbody"},(0,o.kt)("td",{parentName:"tr",align:"left"},(0,o.kt)("strong",{parentName:"td"},"Allele")),(0,o.kt)("td",{parentName:"tr",align:"center"},"1"),(0,o.kt)("td",{parentName:"tr",align:"center"},"2"),(0,o.kt)("td",{parentName:"tr",align:"center"},"3"),(0,o.kt)("td",{parentName:"tr",align:"center"},"4")))),(0,o.kt)("div",{className:"admonition admonition-caution alert alert--warning"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"16",height:"16",viewBox:"0 0 16 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M8.893 1.5c-.183-.31-.52-.5-.887-.5s-.703.19-.886.5L.138 13.499a.98.98 0 0 0 0 1.001c.193.31.53.501.886.501h13.964c.367 0 .704-.19.877-.5a1.03 1.03 0 0 0 .01-1.002L8.893 1.5zm.133 11.497H6.987v-2.003h2.039v2.003zm0-3.004H6.987V5.987h2.039v4.006z"}))),"Filter files beforehand")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},"Keep in mind, BCF/VCF files need to be filtered ",(0,o.kt)("strong",{parentName:"p"},"before")," importing them into PopGen.jl. There is no and will be no VCF-filtering functionality to this package, as it is outside of the purpose of PopGen.jl. Refer to ",(0,o.kt)("inlineCode",{parentName:"p"},"vcftools"),", ",(0,o.kt)("inlineCode",{parentName:"p"},"bcftools"),", and ",(0,o.kt)("inlineCode",{parentName:"p"},"vcflib")," to filter your sequence data. "))),(0,o.kt)("h3",{id:"what-bcfvcf-files-contain"},"What BCF/VCF files contain"),(0,o.kt)("p",null,"Due to the nature of the file format, importing variant call files ",(0,o.kt)("strong",{parentName:"p"},"will")," provide:"),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},"sample names"),(0,o.kt)("li",{parentName:"ul"},"ploidy of each sample"),(0,o.kt)("li",{parentName:"ul"},"locus names"),(0,o.kt)("li",{parentName:"ul"},"genotypes")),(0,o.kt)("h3",{id:"what-bcfvcf-files-lack"},"What BCF/VCF files lack"),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},"population information"),(0,o.kt)("li",{parentName:"ul"},"geographical coordinate information")),(0,o.kt)("p",null,"This means you will need to add that information separately afterwards. Location data (which is ",(0,o.kt)("em",{parentName:"p"},"optional"),") can be added to the ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," with the ",(0,o.kt)("inlineCode",{parentName:"p"},"locations!")," command. Population names (",(0,o.kt)("em",{parentName:"p"},"mandatory"),") can be added using ",(0,o.kt)("inlineCode",{parentName:"p"},"populations!()")),(0,o.kt)("h2",{id:"acknowledgements"},"Acknowledgements"),(0,o.kt)("p",null,"The heavy lifting of the BCF/VCF reader is thanks to the tremendous efforts of the contributors involved with\n",(0,o.kt)("a",{parentName:"p",href:"https://github.com/BioJulia/GeneticVariation.jl"},"GeneticVariation.jl"),", and its successor ",(0,o.kt)("a",{parentName:"p",href:"https://github.com/rasmushenningsson/VariantCallFormat.jl"},"VariantCallFormat.jl"),"\nwhich we use to parse files into ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," format. More specifically, the two packages use a file parser created from ",(0,o.kt)("a",{parentName:"p",href:"https://github.com/BioJulia/Automa.jl"},"Automa.jl"),". If you love the file importer, then give those folks your thanks. If something is wrong and/or you hate the importer, blame us\nfirst (and please ",(0,o.kt)("a",{parentName:"p",href:"https://github.com/biojulia/PopGenCore.jl/issues"},"open up an issue")," \ud83d\ude05)."))}c.isMDXComponent=!0}}]);