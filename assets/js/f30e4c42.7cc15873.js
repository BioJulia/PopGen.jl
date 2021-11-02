"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[9215],{3905:function(e,t,r){r.d(t,{Zo:function(){return s},kt:function(){return d}});var n=r(7294);function a(e,t,r){return t in e?Object.defineProperty(e,t,{value:r,enumerable:!0,configurable:!0,writable:!0}):e[t]=r,e}function l(e,t){var r=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),r.push.apply(r,n)}return r}function i(e){for(var t=1;t<arguments.length;t++){var r=null!=arguments[t]?arguments[t]:{};t%2?l(Object(r),!0).forEach((function(t){a(e,t,r[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(r)):l(Object(r)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(r,t))}))}return e}function o(e,t){if(null==e)return{};var r,n,a=function(e,t){if(null==e)return{};var r,n,a={},l=Object.keys(e);for(n=0;n<l.length;n++)r=l[n],t.indexOf(r)>=0||(a[r]=e[r]);return a}(e,t);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);for(n=0;n<l.length;n++)r=l[n],t.indexOf(r)>=0||Object.prototype.propertyIsEnumerable.call(e,r)&&(a[r]=e[r])}return a}var u=n.createContext({}),p=function(e){var t=n.useContext(u),r=t;return e&&(r="function"==typeof e?e(t):i(i({},t),e)),r},s=function(e){var t=p(e.components);return n.createElement(u.Provider,{value:t},e.children)},c={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},m=n.forwardRef((function(e,t){var r=e.components,a=e.mdxType,l=e.originalType,u=e.parentName,s=o(e,["components","mdxType","originalType","parentName"]),m=p(r),d=a,f=m["".concat(u,".").concat(d)]||m[d]||c[d]||l;return r?n.createElement(f,i(i({ref:t},s),{},{components:r})):n.createElement(f,i({ref:t},s))}));function d(e,t){var r=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var l=r.length,i=new Array(l);i[0]=m;var o={};for(var u in t)hasOwnProperty.call(t,u)&&(o[u]=t[u]);o.originalType=e,o.mdxType="string"==typeof e?e:a,i[1]=o;for(var p=2;p<l;p++)i[p]=r[p];return n.createElement.apply(null,i)}return n.createElement.apply(null,r)}m.displayName="MDXCreateElement"},3946:function(e,t,r){r.r(t),r.d(t,{frontMatter:function(){return o},contentTitle:function(){return u},metadata:function(){return p},toc:function(){return s},default:function(){return m}});var n=r(7462),a=r(3366),l=(r(7294),r(3905)),i=["components"],o={id:"structure",title:"Structure.jl",sidebar_label:"Structure.jl"},u=void 0,p={unversionedId:"api/PopGenCore/structure",id:"api/PopGenCore/structure",isDocsHomePage:!1,title:"Structure.jl",description:"PopGenCore.jl/src/io/Structure.jl",source:"@site/docs/api/PopGenCore/Structure.md",sourceDirName:"api/PopGenCore",slug:"/api/PopGenCore/structure",permalink:"/PopGen.jl/docs/api/PopGenCore/structure",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGenCore/Structure.md",tags:[],version:"current",lastUpdatedAt:1635534271,formattedLastUpdatedAt:"10/29/2021",frontMatter:{id:"structure",title:"Structure.jl",sidebar_label:"Structure.jl"},sidebar:"docs",previous:{title:"ReadWrite.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/read"},next:{title:"PopData.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/types"}},s=[{value:"PopGenCore.jl/src/io/Structure.jl",id:"popgencorejlsrciostructurejl",children:[{value:"\u2757phase_structure",id:"phase_structure",children:[{value:"Example",id:"example",children:[],level:4}],level:3},{value:"\ud83d\udfea structure",id:"-structure",children:[{value:"Keyword Arguments",id:"keyword-arguments",children:[],level:4},{value:"File must follow this Structure format:",id:"file-must-follow-this-structure-format",children:[],level:4},{value:"Structure file example:",id:"structure-file-example",children:[],level:4},{value:"fastStructure file format:",id:"faststructure-file-format",children:[],level:4},{value:"fastStructure file example:",id:"faststructure-file-example",children:[],level:4},{value:"keyword arguments",id:"keyword-arguments-1",children:[],level:4}],level:3}],level:2}],c={toc:s};function m(e){var t=e.components,r=(0,a.Z)(e,i);return(0,l.kt)("wrapper",(0,n.Z)({},c,r,{components:t,mdxType:"MDXLayout"}),(0,l.kt)("h2",{id:"popgencorejlsrciostructurejl"},"PopGenCore.jl/src/io/Structure.jl"),(0,l.kt)("p",null,"\u2757 => not exported |\n\ud83d\udfea => exported by PopGenCore.jl |\n\ud83d\udd35 => exported by PopGen.jl"),(0,l.kt)("h3",{id:"phase_structure"},"\u2757phase_structure"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"phase_structure(datatype::DataType, args...)\n")),(0,l.kt)("p",null,"Takes a DataType (such as ",(0,l.kt)("inlineCode",{parentName:"p"},"Int8"),") and a series of integers to return\na sorted Tuple of those integers converted to that DataType. i.e. takes\na series of alleles and returns a genotype. Returns ",(0,l.kt)("inlineCode",{parentName:"p"},"missing")," if args are\n",(0,l.kt)("inlineCode",{parentName:"p"},"missing"),". Used internally in PopGen.structure file reader."),(0,l.kt)("h4",{id:"example"},"Example"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"phase_structure(Int8, 1,2,3,4,3,4,6,1)\n(1, 1, 2, 3, 3, 4, 4, 6)\n\nphase_structure(Int16, missing, missing)\nmissing\n")),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-structure"},"\ud83d\udfea structure"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"structure(infile::String; kwargs...)\n")),(0,l.kt)("p",null,"Load a Structure format file into memory as a PopData object."),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"infile::String")," : path to Structure file")),(0,l.kt)("h4",{id:"keyword-arguments"},"Keyword Arguments"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"extracols::Integer"),": how many additional optional columns there are beyond Stucture's POPDATA the reader needs to ignore (default: ",(0,l.kt)("inlineCode",{parentName:"li"},"0"),")",(0,l.kt)("ul",{parentName:"li"},(0,l.kt)("li",{parentName:"ul"},"these include POPFLAG, LOCDATA, or anything else you might have added"))),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"extrarows::Integer")," : how many additional optional rows there are beyond the first row of locus names (default: ",(0,l.kt)("inlineCode",{parentName:"li"},"0"),")"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"missingval::String"),"  : the value used to identify missing values in the data (default: ",(0,l.kt)("inlineCode",{parentName:"li"},'"-9"'),")"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"silent::Bool"),"   : whether to print file information during import (default: ",(0,l.kt)("inlineCode",{parentName:"li"},"false"),")"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"allow_monomorphic::Bool")," : whether to keep monomorphic loci in the dataset (default: ",(0,l.kt)("inlineCode",{parentName:"li"},"false"),")"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"faststructure::Bool"),": whether the file is fastStructure format (default: ",(0,l.kt)("inlineCode",{parentName:"li"},"false"),")")),(0,l.kt)("h4",{id:"file-must-follow-this-structure-format"},"File must follow this Structure format:"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"the file is ",(0,l.kt)("inlineCode",{parentName:"li"},"tab")," or ",(0,l.kt)("inlineCode",{parentName:"li"},"space")," delimited ",(0,l.kt)("strong",{parentName:"li"},"but not both")),(0,l.kt)("li",{parentName:"ul"},"first row is locus names separated by the delimiter",(0,l.kt)("ul",{parentName:"li"},(0,l.kt)("li",{parentName:"ul"},"leading/trailing whitespaces are tolerated"),(0,l.kt)("li",{parentName:"ul"},"optional rows allowed ",(0,l.kt)("strong",{parentName:"li"},"after")," the locus names"))),(0,l.kt)("li",{parentName:"ul"},"number of rows per sample = ploidy",(0,l.kt)("ul",{parentName:"li"},(0,l.kt)("li",{parentName:"ul"},"e.g. if diploid, that sample would have 2 rows"),(0,l.kt)("li",{parentName:"ul"},"multi-column variant not supported"))),(0,l.kt)("li",{parentName:"ul"},"first data column is sample name"),(0,l.kt)("li",{parentName:"ul"},"second data column is population ID",(0,l.kt)("ul",{parentName:"li"},(0,l.kt)("li",{parentName:"ul"},"optional columns allowed ",(0,l.kt)("strong",{parentName:"li"},"after")," the population ID (2nd) column"))),(0,l.kt)("li",{parentName:"ul"},"remaining columns are the genotype for that individual for that locus")),(0,l.kt)("h4",{id:"structure-file-example"},"Structure file example:"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"locus_1 locus_2 locus_3 locus_4 locus_5\nwalnut_01   1   -9  145 66  0   92\nwalnut_01   1   -9  -9  64  0   94\nwalnut_02   1   106 142 68  1   92\nwalnut_02   1   106 148 64  0   94\nwalnut_03   2   110 145 -9  0   92\nwalnut_03   2   110 148 66  1   -9\n")),(0,l.kt)("h4",{id:"faststructure-file-format"},"fastStructure file format:"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"the file is ",(0,l.kt)("inlineCode",{parentName:"li"},"tab")," or ",(0,l.kt)("inlineCode",{parentName:"li"},"space")," delimited ",(0,l.kt)("strong",{parentName:"li"},"but not both")),(0,l.kt)("li",{parentName:"ul"},"no first row of loci names"),(0,l.kt)("li",{parentName:"ul"},"number of rows per sample = ploidy",(0,l.kt)("ul",{parentName:"li"},(0,l.kt)("li",{parentName:"ul"},"e.g. if diploid, that sample would have 2 rows"))),(0,l.kt)("li",{parentName:"ul"},"first data column is sample name"),(0,l.kt)("li",{parentName:"ul"},"second data column is population ID"),(0,l.kt)("li",{parentName:"ul"},"remaining columns are the genotype for that individual for that locus"),(0,l.kt)("li",{parentName:"ul"},"usually, first 6 colums are empty (but not necessary)"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("strong",{parentName:"li"},"no")," extra rows or columns.")),(0,l.kt)("h4",{id:"faststructure-file-example"},"fastStructure file example:"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"chestnut_01 1   -9  145 66  0   92\nchestnut_01 1   -9  -9  64  0   94\nchestnut_02 1   106 142 68  1   92\nchestnut_02 1   106 148 64  0   94\nchestnut_03 2   110 145 -9  0   92\nchestnut_03 2   110 148 66  1   -9\n")),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'walnuts = structure("juglans_nigra.str", extracols = 0, extrarows = 0)\n')),(0,l.kt)("hr",null),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"structure(data::PopData; filename::String, faststructure::Bool, delim::String)\n")),(0,l.kt)("p",null,"Write a ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," object to a Stucture format file"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"data"),": the ",(0,l.kt)("inlineCode",{parentName:"li"},"PopData")," object you wish to convert to a Structure file")),(0,l.kt)("h4",{id:"keyword-arguments-1"},"keyword arguments"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"filename"),": a ",(0,l.kt)("inlineCode",{parentName:"li"},"String")," of the output filename"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"delim")," : a ",(0,l.kt)("inlineCode",{parentName:"li"},"String")," of either ",(0,l.kt)("inlineCode",{parentName:"li"},'"tab"')," or ",(0,l.kt)("inlineCode",{parentName:"li"},'"space"')," indicating the delimiter (default: ",(0,l.kt)("inlineCode",{parentName:"li"},'"tab"'),")"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"faststructure"),": true/false of whether the output should be formatted for fastStructure (default: ",(0,l.kt)("inlineCode",{parentName:"li"},"false"),")")),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'cats = @nancycats;\nfewer_cats = omit(cats, name = samples(cats)[1:10]);\nstructure(fewer_cats, filename = "filtered_nancycats.str", faststructure = true)\n')))}m.isMDXComponent=!0}}]);