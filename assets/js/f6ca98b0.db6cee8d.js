"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[9535],{3905:function(e,l,n){n.d(l,{Zo:function(){return u},kt:function(){return d}});var t=n(7294);function r(e,l,n){return l in e?Object.defineProperty(e,l,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[l]=n,e}function a(e,l){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var t=Object.getOwnPropertySymbols(e);l&&(t=t.filter((function(l){return Object.getOwnPropertyDescriptor(e,l).enumerable}))),n.push.apply(n,t)}return n}function o(e){for(var l=1;l<arguments.length;l++){var n=null!=arguments[l]?arguments[l]:{};l%2?a(Object(n),!0).forEach((function(l){r(e,l,n[l])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):a(Object(n)).forEach((function(l){Object.defineProperty(e,l,Object.getOwnPropertyDescriptor(n,l))}))}return e}function p(e,l){if(null==e)return{};var n,t,r=function(e,l){if(null==e)return{};var n,t,r={},a=Object.keys(e);for(t=0;t<a.length;t++)n=a[t],l.indexOf(n)>=0||(r[n]=e[n]);return r}(e,l);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(t=0;t<a.length;t++)n=a[t],l.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var i=t.createContext({}),c=function(e){var l=t.useContext(i),n=l;return e&&(n="function"==typeof e?e(l):o(o({},l),e)),n},u=function(e){var l=c(e.components);return t.createElement(i.Provider,{value:l},e.children)},s={inlineCode:"code",wrapper:function(e){var l=e.children;return t.createElement(t.Fragment,{},l)}},f=t.forwardRef((function(e,l){var n=e.components,r=e.mdxType,a=e.originalType,i=e.parentName,u=p(e,["components","mdxType","originalType","parentName"]),f=c(n),d=r,m=f["".concat(i,".").concat(d)]||f[d]||s[d]||a;return n?t.createElement(m,o(o({ref:l},u),{},{components:n})):t.createElement(m,o({ref:l},u))}));function d(e,l){var n=arguments,r=l&&l.mdxType;if("string"==typeof e||r){var a=n.length,o=new Array(a);o[0]=f;var p={};for(var i in l)hasOwnProperty.call(l,i)&&(p[i]=l[i]);p.originalType=e,p.mdxType="string"==typeof e?e:r,o[1]=p;for(var c=2;c<a;c++)o[c]=n[c];return t.createElement.apply(null,o)}return t.createElement.apply(null,n)}f.displayName="MDXCreateElement"},2844:function(e,l,n){n.r(l),n.d(l,{frontMatter:function(){return p},contentTitle:function(){return i},metadata:function(){return c},toc:function(){return u},default:function(){return f}});var t=n(7462),r=n(3366),a=(n(7294),n(3905)),o=["components"],p={id:"allelefreq",title:"AlleleFreq.jl",sidebar_label:"AlleleFreq.jl"},i=void 0,c={unversionedId:"api/PopGenCore/allelefreq",id:"api/PopGenCore/allelefreq",isDocsHomePage:!1,title:"AlleleFreq.jl",description:"PopGenCore.jl/src/AlleleFreq.jl",source:"@site/docs/api/PopGenCore/AlleleFreq.md",sourceDirName:"api/PopGenCore",slug:"/api/PopGenCore/allelefreq",permalink:"/PopGen.jl/docs/api/PopGenCore/allelefreq",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGenCore/AlleleFreq.md",tags:[],version:"current",lastUpdatedAt:1635534271,formattedLastUpdatedAt:"10/29/2021",frontMatter:{id:"allelefreq",title:"AlleleFreq.jl",sidebar_label:"AlleleFreq.jl"},sidebar:"docs",previous:{title:"Utils.jl",permalink:"/PopGen.jl/docs/api/PopGen/utils"},next:{title:"Conditionals.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/conditionals"}},u=[{value:"PopGenCore.jl/src/AlleleFreq.jl",id:"popgencorejlsrcallelefreqjl",children:[{value:"\ud83d\udfea allelefreq",id:"-allelefreq",children:[],level:3},{value:"\ud83d\udfea avg_allelefreq",id:"-avg_allelefreq",children:[],level:3},{value:"\ud83d\udfea allelefreq_vec",id:"-allelefreq_vec",children:[],level:3}],level:2}],s={toc:u};function f(e){var l=e.components,n=(0,r.Z)(e,o);return(0,a.kt)("wrapper",(0,t.Z)({},s,n,{components:l,mdxType:"MDXLayout"}),(0,a.kt)("h2",{id:"popgencorejlsrcallelefreqjl"},"PopGenCore.jl/src/AlleleFreq.jl"),(0,a.kt)("p",null,"\u2757 => not exported |\n\ud83d\udfea => exported by PopGenCore.jl |\n\ud83d\udd35 => exported by PopGen.jl"),(0,a.kt)("h3",{id:"-allelefreq"},"\ud83d\udfea allelefreq"),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},"allelefreq(allele::Int, genos::GenoArray)\n")),(0,a.kt)("p",null,"Return the frequency of a specific ",(0,a.kt)("inlineCode",{parentName:"p"},"allele")," from a vector of Genotypes ",(0,a.kt)("inlineCode",{parentName:"p"},"genos"),"."),(0,a.kt)("p",null,(0,a.kt)("strong",{parentName:"p"},"Example")),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre"},'using DataFramesMeta\nncats = @nancycats;\nncats_sub @where(ncats.genodata, :locus .== "fca8", :genotype .!== missing)\npop_grp = groupby(ncats_sub, :population)\nDataFrames.combine(pop_grp, :genotype => (geno -> allelefreq(137, geno)) => :freq_137)\n')),(0,a.kt)("hr",null),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},"allelefreq(geno::Genotype)\n")),(0,a.kt)("p",null,"Return a ",(0,a.kt)("inlineCode",{parentName:"p"},"Dict")," of allele frequencies of the alleles within a single Genotype in a ",(0,a.kt)("inlineCode",{parentName:"p"},"PopData"),"\nobject."),(0,a.kt)("hr",null),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},"allelefreq(locus::T) where T<:GenotypeArray\n")),(0,a.kt)("p",null,"Return a ",(0,a.kt)("inlineCode",{parentName:"p"},"Dict")," of allele frequencies of a single locus in a ",(0,a.kt)("inlineCode",{parentName:"p"},"PopData"),"\nobject."),(0,a.kt)("hr",null),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},"allelefreq(data::PopData, locus::String; population::Bool = false)\n")),(0,a.kt)("p",null,"Return a ",(0,a.kt)("inlineCode",{parentName:"p"},"Dict")," of allele frequencies of a single locus in a ",(0,a.kt)("inlineCode",{parentName:"p"},"PopData"),"\nobject. Use ",(0,a.kt)("inlineCode",{parentName:"p"},"population = true")," to return a table of allele frequencies\nof that locus per population."),(0,a.kt)("p",null,(0,a.kt)("strong",{parentName:"p"},"Example")),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},'cats = @nancycats\nallelefreq(cats, "fca8")\nallelefreq(cats, "fca8", population = true)\n')),(0,a.kt)("hr",null),(0,a.kt)("h3",{id:"-avg_allelefreq"},"\ud83d\udfea avg_allelefreq"),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},"avg_allelefreq(allele_dicts::AbstractVector{Dict{T, Float64}}, power::Int = 1) where T<:Signed  \n")),(0,a.kt)("p",null,"Takes a Vector of Dicts generated by ",(0,a.kt)("inlineCode",{parentName:"p"},"allelefreq")," and returns a Dict of the average\nallele frequencies raised to the ",(0,a.kt)("inlineCode",{parentName:"p"},"power")," (exponent) specified (default: ",(0,a.kt)("inlineCode",{parentName:"p"},"1"),").\nThis is typically done to calculate average allele frequencies across populations."),(0,a.kt)("p",null,(0,a.kt)("strong",{parentName:"p"},"Example")),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre"},"cats = @nancycats;\nalleles_df = DataFrames.combine(\n    groupby(cats.genodata, [:locus, :population]),\n    :genotype => allelefreq => :alleles\n);\nDataFrames.combine(\n    groupby(alleles_df, :locus),\n    :alleles => (i -> sum(avg_allelefreq(i, 2))) => :avg_freq\n)\n")),(0,a.kt)("hr",null),(0,a.kt)("h3",{id:"-allelefreq_vec"},"\ud83d\udfea allelefreq_vec"),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},"allelefreq_vec(locus::T) where T<:GenotypeArray\n")),(0,a.kt)("p",null,"Return a Vector of allele frequencies of a single locus in a ",(0,a.kt)("inlineCode",{parentName:"p"},"PopData")," object. Similar to ",(0,a.kt)("inlineCode",{parentName:"p"},"allelefreq()"),", except it returns only the frequencies, without the allele names, meaning they can be in any order. This can be useful for getting the expected genotype frequencies."))}f.isMDXComponent=!0}}]);