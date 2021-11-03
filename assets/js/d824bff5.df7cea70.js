"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[3271],{3905:function(e,n,t){t.d(n,{Zo:function(){return d},kt:function(){return c}});var a=t(7294);function o(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function l(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);n&&(a=a.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,a)}return t}function i(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?l(Object(t),!0).forEach((function(n){o(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):l(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function r(e,n){if(null==e)return{};var t,a,o=function(e,n){if(null==e)return{};var t,a,o={},l=Object.keys(e);for(a=0;a<l.length;a++)t=l[a],n.indexOf(t)>=0||(o[t]=e[t]);return o}(e,n);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);for(a=0;a<l.length;a++)t=l[a],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(o[t]=e[t])}return o}var p=a.createContext({}),s=function(e){var n=a.useContext(p),t=n;return e&&(t="function"==typeof e?e(n):i(i({},n),e)),t},d=function(e){var n=s(e.components);return a.createElement(p.Provider,{value:n},e.children)},u={inlineCode:"code",wrapper:function(e){var n=e.children;return a.createElement(a.Fragment,{},n)}},m=a.forwardRef((function(e,n){var t=e.components,o=e.mdxType,l=e.originalType,p=e.parentName,d=r(e,["components","mdxType","originalType","parentName"]),m=s(t),c=o,k=m["".concat(p,".").concat(c)]||m[c]||u[c]||l;return t?a.createElement(k,i(i({ref:n},d),{},{components:t})):a.createElement(k,i({ref:n},d))}));function c(e,n){var t=arguments,o=n&&n.mdxType;if("string"==typeof e||o){var l=t.length,i=new Array(l);i[0]=m;var r={};for(var p in n)hasOwnProperty.call(n,p)&&(r[p]=n[p]);r.originalType=e,r.mdxType="string"==typeof e?e:o,i[1]=r;for(var s=2;s<l;s++)i[s]=t[s];return a.createElement.apply(null,i)}return a.createElement.apply(null,t)}m.displayName="MDXCreateElement"},6066:function(e,n,t){t.r(n),t.d(n,{frontMatter:function(){return r},contentTitle:function(){return p},metadata:function(){return s},toc:function(){return d},default:function(){return m}});var a=t(7462),o=t(3366),l=(t(7294),t(3905)),i=["components"],r={id:"manipulate",title:"Manipulate.jl",sidebar_label:"Manipulate.jl"},p=void 0,s={unversionedId:"api/PopGenCore/manipulate",id:"api/PopGenCore/manipulate",isDocsHomePage:!1,title:"Manipulate.jl",description:"PopGenCore.jl/src/Manipulate.jl",source:"@site/docs/api/PopGenCore/Manipulate.md",sourceDirName:"api/PopGenCore",slug:"/api/PopGenCore/manipulate",permalink:"/PopGen.jl/docs/api/PopGenCore/manipulate",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGenCore/Manipulate.md",tags:[],version:"current",lastUpdatedAt:1635540406,formattedLastUpdatedAt:"10/29/2021",frontMatter:{id:"manipulate",title:"Manipulate.jl",sidebar_label:"Manipulate.jl"},sidebar:"docs",previous:{title:"Iterators.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/iterators"},next:{title:"MathUtils.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/mathutils"}},d=[{value:"PopGenCore.jl/src/Manipulate.jl",id:"popgencorejlsrcmanipulatejl",children:[{value:"\ud83d\udfea\ud83d\udd35 sampleinfo!",id:"-sampleinfo",children:[{value:"Arguments",id:"arguments",children:[],level:4},{value:"Keyword Arguments",id:"keyword-arguments",children:[],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 locusinfo!",id:"-locusinfo",children:[{value:"Arguments",id:"arguments-1",children:[],level:4},{value:"Keyword Arguments",id:"keyword-arguments-1",children:[],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 locationdata!",id:"-locationdata",children:[{value:"Formatting requirements",id:"formatting-requirements",children:[],level:4},{value:"Formatting requirements",id:"formatting-requirements-1",children:[{value:"NOTE",id:"note",children:[],level:5}],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 populations!",id:"-populations",children:[{value:"Rename using a Dictionary",id:"rename-using-a-dictionary",children:[],level:4},{value:"Rename using a Vector of Strings",id:"rename-using-a-vector-of-strings",children:[],level:4},{value:"Reassign using samples and new population assignments",id:"reassign-using-samples-and-new-population-assignments",children:[],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 exclude!",id:"-exclude",children:[{value:"Keyword Arguments",id:"keyword-arguments-2",children:[],level:4},{value:"<code>locus</code>",id:"locus",children:[],level:4},{value:"<code>population</code>",id:"population",children:[],level:4},{value:"<code>name</code>",id:"name",children:[],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 exclude",id:"-exclude-1",children:[{value:"Keyword Arguments",id:"keyword-arguments-3",children:[],level:4},{value:"<code>locus</code>",id:"locus-1",children:[],level:4},{value:"<code>population</code>",id:"population-1",children:[],level:4},{value:"<code>name</code>",id:"name-1",children:[],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 keep!",id:"-keep",children:[{value:"Keyword Arguments",id:"keyword-arguments-4",children:[],level:4},{value:"<code>locus</code>",id:"locus-2",children:[],level:4},{value:"<code>population</code>",id:"population-2",children:[],level:4},{value:"<code>name</code>",id:"name-2",children:[],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 keep",id:"-keep-1",children:[{value:"Keyword Arguments",id:"keyword-arguments-5",children:[],level:4},{value:"<code>locus</code>",id:"locus-3",children:[],level:4},{value:"<code>population</code>",id:"population-3",children:[],level:4},{value:"<code>name</code>",id:"name-3",children:[],level:4}],level:3},{value:"\ud83d\udfea\ud83d\udd35 filter",id:"-filter",children:[],level:3},{value:"\ud83d\udfea\ud83d\udd35 filter!",id:"-filter-1",children:[],level:3}],level:2}],u={toc:d};function m(e){var n=e.components,t=(0,o.Z)(e,i);return(0,l.kt)("wrapper",(0,a.Z)({},u,t,{components:n,mdxType:"MDXLayout"}),(0,l.kt)("h2",{id:"popgencorejlsrcmanipulatejl"},"PopGenCore.jl/src/Manipulate.jl"),(0,l.kt)("p",null,"\u2757 => not exported |\n\ud83d\udfea => exported by PopGenCore.jl |\n\ud83d\udd35 => exported by PopGen.jl"),(0,l.kt)("h3",{id:"-sampleinfo"},"\ud83d\udfea\ud83d\udd35 sampleinfo!"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"sampleinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)\nsampleinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)\n")),(0,l.kt)("p",null,"Add an additional sample information to ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," metadata. Mutates ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," in place. Metadata\nmust be in the same order as the samples in ",(0,l.kt)("inlineCode",{parentName:"p"},"sampleinfo(popdata)"),"."),(0,l.kt)("h4",{id:"arguments"},"Arguments"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"metadata")," : A Pair of :ColumnName => ","[Values]")),(0,l.kt)("h4",{id:"keyword-arguments"},"Keyword Arguments"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"categorical"),' : Boolean of whether the metadata being added is categorical aka "factors" (default: ',(0,l.kt)("inlineCode",{parentName:"li"},"false"),")\n",(0,l.kt)("strong",{parentName:"li"},"Example"))),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'cats = @nancycats\nsampleinfo!(cats, :whiskerlength => rand(cats.metadata.samples))\nsampleinfo!(cats, "tailcolor" => rand(["orange", "brown"], metadata(cats).samples), categorical = true)\ncats\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 237\n  Populations: 17\n  Other Info: ["whiskerlength", "tailcolor"]\n')),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-locusinfo"},"\ud83d\udfea\ud83d\udd35 locusinfo!"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"locusinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)\nlocusinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)\n")),(0,l.kt)("p",null,"Add an additional locus information to ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," metadata. Mutates ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," in place. Metadata\nmust be in the same order as the samples in ",(0,l.kt)("inlineCode",{parentName:"p"},"locusinfo(PopData)"),"."),(0,l.kt)("h4",{id:"arguments-1"},"Arguments"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"metadata")," : A Pair of ",(0,l.kt)("inlineCode",{parentName:"li"},":ColumnName => [Values]"))),(0,l.kt)("h4",{id:"keyword-arguments-1"},"Keyword Arguments"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},"categorical"),' : Boolean of whether the metadata being added is categorical aka "factors" (default: ',(0,l.kt)("inlineCode",{parentName:"li"},"false"),")\n",(0,l.kt)("strong",{parentName:"li"},"Example"))),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'cats = @nancycats\nlocusinfo!(cats, :quality => rand(metadata(cats).loci))\ncats\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 237\n  Populations: 17\n  Other Info: ["quality"]\n')),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-locationdata"},"\ud83d\udfea\ud83d\udd35 locationdata!"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"locationdata!(data::PopData; longitude::Vector{Float64}, latitude::Vector{Float64})\n")),(0,l.kt)("p",null,"Replaces existing ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," geographic coordinate data.\nTakes ",(0,l.kt)("strong",{parentName:"p"},"decimal degrees")," as a ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector")," of any ",(0,l.kt)("inlineCode",{parentName:"p"},"AbstractFloat"),"."),(0,l.kt)("h4",{id:"formatting-requirements"},"Formatting requirements"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"Decimal Degrees format: ",(0,l.kt)("inlineCode",{parentName:"li"},"-11.431")),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("strong",{parentName:"li"},"Must")," use negative sign ",(0,l.kt)("inlineCode",{parentName:"li"},"-")," instead of cardinal directions"),(0,l.kt)("li",{parentName:"ul"},"Location data must be in the order that samples appear in your ",(0,l.kt)("inlineCode",{parentName:"li"},"PopData")),(0,l.kt)("li",{parentName:"ul"},"Missing data should be coded as ",(0,l.kt)("inlineCode",{parentName:"li"},"missing")," values of type ",(0,l.kt)("inlineCode",{parentName:"li"},"Missing")," (can be accomplished with ",(0,l.kt)("inlineCode",{parentName:"li"},"replace!()"),")")),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"ncats = @nancycats ;\nx = rand(237) ; y = rand(237)\nlocationdata!(ncats, longitude = x, latitude = y)\n")),(0,l.kt)("hr",null),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"locationdata!(data::PopData; longitude::Vector{String}, latitude::Vector{String})\n")),(0,l.kt)("p",null,"Replaces existing ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," geographic coordinate data. Takes\n",(0,l.kt)("strong",{parentName:"p"},"decimal minutes")," or ",(0,l.kt)("strong",{parentName:"p"},"degrees minutes seconds")," format as a ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector")," of ",(0,l.kt)("inlineCode",{parentName:"p"},"String"),". Recommended to use ",(0,l.kt)("inlineCode",{parentName:"p"},"CSV.read"),"\nfrom ",(0,l.kt)("inlineCode",{parentName:"p"},"CSV.jl")," to import your spatial coordinates from a text file."),(0,l.kt)("h4",{id:"formatting-requirements-1"},"Formatting requirements"),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},"Coordinates as a ",(0,l.kt)("inlineCode",{parentName:"li"},"String")," separated by spaces (",(0,l.kt)("inlineCode",{parentName:"li"},'"11 43 41"'),") or colons (",(0,l.kt)("inlineCode",{parentName:"li"},'"11:43:41"'),")"),(0,l.kt)("li",{parentName:"ul"},"Must use negative sign (",(0,l.kt)("inlineCode",{parentName:"li"},'"-11 43.52"'),") or single-letter cardinal direction (",(0,l.kt)("inlineCode",{parentName:"li"},'"11 43.52W"'),")"),(0,l.kt)("li",{parentName:"ul"},"Missing data should be coded as the string ",(0,l.kt)("inlineCode",{parentName:"li"},'"missing"')," (can be accomplished with ",(0,l.kt)("inlineCode",{parentName:"li"},"replace!()"),")"),(0,l.kt)("li",{parentName:"ul"},"Can mix colons and spaces (although it's bad practice)")),(0,l.kt)("h5",{id:"note"},"NOTE"),(0,l.kt)("p",null,"If you read in the coordinate data as 4 vectors (longitude degrees, longitude minutes, latitude degrees, latitude minutes),\nthen the easiest course of action would be to merge them into two vectors of strings\n(one for longitude, one for latitude):"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'long_string = string.(lat_deg, " ", lat_min)\nlat_string = string.(long_deg, " ", long_min)\n')),(0,l.kt)("p",null,"and use these as inputs into ",(0,l.kt)("inlineCode",{parentName:"p"},"locations!")),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'ncats = @nancycats;\nx = fill("11 22.33W", 237) ; y = fill("-41 31.52", 237)\nlocationdata!(ncats, longitude = x, latitude = y)\n')),(0,l.kt)("hr",null),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"locationdata!(data::PopData, longitude::Vector{String}, latitude::Vector{String})\nlocationdata!(data::PopData; kwargs...)\n")),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-populations"},"\ud83d\udfea\ud83d\udd35 populations!"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"# Replace by matching\npopulations!(data::PopData, rename::Dict)\npopulations!(data::PopData, rename::Vector{String})\npopulations!(data::PopData, samples::Vector{String}, populations::Vector{String})\n")),(0,l.kt)("p",null,"Multiple methods to rename or reassign population names to a ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("h4",{id:"rename-using-a-dictionary"},"Rename using a Dictionary"),(0,l.kt)("p",null,"Rename existing population ID's of ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," using a ",(0,l.kt)("inlineCode",{parentName:"p"},"Dict")," of\n",(0,l.kt)("inlineCode",{parentName:"p"},"population_name => replacement")),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'potatopops = Dict("1" => "Idaho", "2" => "Russet")\npopulations!(potatoes, potatopops)\n')),(0,l.kt)("h4",{id:"rename-using-a-vector-of-strings"},"Rename using a Vector of Strings"),(0,l.kt)("p",null,"If the number of new names is equal to the number of current unique population names,\nthe method will rename the existing populations in the order with which they appear\nvia ",(0,l.kt)("inlineCode",{parentName:"p"},"unique()"),". If the number of new population names is equal to the number of samples,\nthe method will instead assign new population names to every sample in the order with which they appear in ",(0,l.kt)("inlineCode",{parentName:"p"},"sampleinfo(popdata)"),"."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'# rename (2) existing populations\npotatopops = ["Idaho", "Russet"]\npopulations!(potatoes, potatopops)\n# assign new names to all [44] samples\npotatopops = repeat(["Idaho", "Russet"], inner = 22) ;\npopulations!(potatoes, potatopops)\n')),(0,l.kt)("h4",{id:"reassign-using-samples-and-new-population-assignments"},"Reassign using samples and new population assignments"),(0,l.kt)("p",null,"Completely reassign populations for each individual. Takes two vectors of strings\nas input: one of the sample names, and the other with their new corresponding\npopulation name. This can be useful to change population names for only some individuals."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'populations!(potatoes, ["potato_1", "potato_2"], ["north", "south"])\n')),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-exclude"},"\ud83d\udfea\ud83d\udd35 exclude!"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"exclude!(data::PopData, kwargs...)\n")),(0,l.kt)("p",null,"Edit a ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," object in-place by excluding all occurences of the specified information.\nThe keywords can be used in any combination. Synonymous with ",(0,l.kt)("inlineCode",{parentName:"p"},"omit!")," and ",(0,l.kt)("inlineCode",{parentName:"p"},"remove!"),". All\nvalues are converted to ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," for filtering, so ",(0,l.kt)("inlineCode",{parentName:"p"},"Symbol")," and numbers will also work.\nThis can be considered a simpler and more rudimentary syntax for subsetting\nor filtering PopData."),(0,l.kt)("h4",{id:"keyword-arguments-2"},"Keyword Arguments"),(0,l.kt)("h4",{id:"locus"},(0,l.kt)("inlineCode",{parentName:"h4"},"locus")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of loci you want to remove from the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),".\nThe keyword ",(0,l.kt)("inlineCode",{parentName:"p"},"loci")," also works."),(0,l.kt)("h4",{id:"population"},(0,l.kt)("inlineCode",{parentName:"h4"},"population")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of populations you want to remove from the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"\nThe keyword ",(0,l.kt)("inlineCode",{parentName:"p"},"populations")," also works."),(0,l.kt)("h4",{id:"name"},(0,l.kt)("inlineCode",{parentName:"h4"},"name")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of samples you want to remove from the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"\nThe keywords ",(0,l.kt)("inlineCode",{parentName:"p"},"names"),", ",(0,l.kt)("inlineCode",{parentName:"p"},"sample"),", and ",(0,l.kt)("inlineCode",{parentName:"p"},"samples")," also work."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Examples")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'cats = @nancycats;\nexclude!(cats, name = "N100", population = 1:5)\nexclude!(cats, name = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])\nexclude!(cats, name = "N102", locus = :fca8, population = "3")\n')),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"const omit! = exclude!\nconst remove! = exclude!\n")),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-exclude-1"},"\ud83d\udfea\ud83d\udd35 exclude"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"exclude(data::PopData, kwargs...)\n")),(0,l.kt)("p",null,"Returns a new ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," object excluding all occurrences of the specified keywords.\nThe keywords can be used in any combination. Synonymous with ",(0,l.kt)("inlineCode",{parentName:"p"},"omit")," and ",(0,l.kt)("inlineCode",{parentName:"p"},"remove"),". All\nvalues are converted to ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," for filtering, so ",(0,l.kt)("inlineCode",{parentName:"p"},"Symbol")," and numbers will also work.\nThis can be considered a simpler and more rudimentary syntax for subsetting\nor filtering PopData."),(0,l.kt)("h4",{id:"keyword-arguments-3"},"Keyword Arguments"),(0,l.kt)("h4",{id:"locus-1"},(0,l.kt)("inlineCode",{parentName:"h4"},"locus")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of loci you want to remove from the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("h4",{id:"population-1"},(0,l.kt)("inlineCode",{parentName:"h4"},"population")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of populations you want to remove from the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("h4",{id:"name-1"},(0,l.kt)("inlineCode",{parentName:"h4"},"name")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of samples you want to remove from the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Examples")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'cats = @nancycats;\nexclude(cats, name = "N100", population = 1:5)\nexclude(cats, name = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])\nexclude(cats, name = "N102", locus = :fca8, population = "3")\n')),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"const omit = exclude\nconst remove = exclude\n")),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-keep"},"\ud83d\udfea\ud83d\udd35 keep!"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"keep!(data::PopData, kwargs...)\n")),(0,l.kt)("p",null,"Edit a ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),' object in-place by keeping only the occurrences of the specified keywords.\nIf using multiple fields, they will be chained together as "',(0,l.kt)("inlineCode",{parentName:"p"},"or"),'" statements.\nAll values are converted to ',(0,l.kt)("inlineCode",{parentName:"p"},"String")," for filtering, so ",(0,l.kt)("inlineCode",{parentName:"p"},"Symbol")," and numbers will also work.\nThis can be considered a simpler and more rudimentary syntax for subsetting\nor filtering PopData."),(0,l.kt)("h4",{id:"keyword-arguments-4"},"Keyword Arguments"),(0,l.kt)("h4",{id:"locus-2"},(0,l.kt)("inlineCode",{parentName:"h4"},"locus")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of loci you want to keep in the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("h4",{id:"population-2"},(0,l.kt)("inlineCode",{parentName:"h4"},"population")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of populations you want to keep in the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("h4",{id:"name-2"},(0,l.kt)("inlineCode",{parentName:"h4"},"name")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of samples you want to keep in the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Examples")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'cats = @nancycats;\nkeep!(cats, population = 1:5)\n# keep 4 populations and 3 specific samples\nkeep!(cats, name = ["N100", "N102", "N211"])\n# keep 2 loci, 2 populations, and 10 specific individuals\nkeep!(cats, locus = [:fca8, "fca37"], population = [7,8], name = samplenames(cats)[1:10])\n')),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-keep-1"},"\ud83d\udfea\ud83d\udd35 keep"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"keep(data::PopData, kwargs...)\n")),(0,l.kt)("p",null,"Returns a new ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," object keeping only the occurrences of the specified keyword.\nUnlike ",(0,l.kt)("inlineCode",{parentName:"p"},"exclude()"),". only one keyword can be used at a time. All values are\nconverted to ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," for filtering, so ",(0,l.kt)("inlineCode",{parentName:"p"},"Symbol")," and numbers will also work.\nThis can be considered a simpler and more rudimentary syntax for subsetting\nor filtering PopData."),(0,l.kt)("h4",{id:"keyword-arguments-5"},"Keyword Arguments"),(0,l.kt)("h4",{id:"locus-3"},(0,l.kt)("inlineCode",{parentName:"h4"},"locus")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of loci you want to keep in the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("h4",{id:"population-3"},(0,l.kt)("inlineCode",{parentName:"h4"},"population")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of populations you want to keep in the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("h4",{id:"name-3"},(0,l.kt)("inlineCode",{parentName:"h4"},"name")),(0,l.kt)("p",null,"A ",(0,l.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,l.kt)("inlineCode",{parentName:"p"},"Vector{String}")," of samples you want to keep in the ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),"."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Examples")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'cats = @nancycats;\nkeep(cats, population = 1:5)\n# equivalent to cats[cats.genodata.population .\u2208 Ref(string.(1:5)), :]\nkeep(cats, name = ["N100", "N102", "N211"])\n# equivalent to cats[cats.genodata.name .\u2208 Ref(["N100", "N102", "N211"]), :]\nkeep(cats, locus = [:fca8, "fca37"])\n# equivalent to cats[cats.genodata.locus .\u2208 Ref(["fca8", "fca37"]), :]\n')),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-filter"},"\ud83d\udfea\ud83d\udd35 filter"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"filter(data::PopData, args...)\n")),(0,l.kt)("p",null,"A drop-in replacement for DataFrames.filter where ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," is the first\nargument and the filtering conditions are the second argument. Returns a\nnew ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),". ",(0,l.kt)("strong",{parentName:"p"},"Note")," the argument order is opposite of that from DataFrames.jl.\n",(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"x = @nancycats ;\ny = filter(x, :name => i -> i \u2208 samplenames(x)[1:10]) ;\nshow(x)\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 237\n  Populations: 17\nshow(y)\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 10\n  Populations: 1\n")),(0,l.kt)("hr",null),(0,l.kt)("h3",{id:"-filter-1"},"\ud83d\udfea\ud83d\udd35 filter!"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"filter!(data::PopData, args...)\n")),(0,l.kt)("p",null,"A drop-in replacement for the DataFrames.filter! where ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," is the first\nargument and the filtering conditions are the second argument. Mutates the\n",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," in place and returns it. ",(0,l.kt)("strong",{parentName:"p"},"Note")," the argument order is opposite of that from DataFrames.jl."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"x = @nancycats ;\nfilter!(x, :name => i -> i \u2208 samplenames(x)[1:10]) ;\nshow(x)\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 10\n  Populations: 1\n")))}m.isMDXComponent=!0}}]);