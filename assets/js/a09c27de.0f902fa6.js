"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[1506],{3905:function(e,n,t){t.d(n,{Zo:function(){return u},kt:function(){return d}});var l=t(7294);function i(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function r(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);n&&(l=l.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,l)}return t}function a(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?r(Object(t),!0).forEach((function(n){i(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):r(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function s(e,n){if(null==e)return{};var t,l,i=function(e,n){if(null==e)return{};var t,l,i={},r=Object.keys(e);for(l=0;l<r.length;l++)t=r[l],n.indexOf(t)>=0||(i[t]=e[t]);return i}(e,n);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);for(l=0;l<r.length;l++)t=r[l],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(i[t]=e[t])}return i}var o=l.createContext({}),p=function(e){var n=l.useContext(o),t=n;return e&&(t="function"==typeof e?e(n):a(a({},n),e)),t},u=function(e){var n=p(e.components);return l.createElement(o.Provider,{value:n},e.children)},c={inlineCode:"code",wrapper:function(e){var n=e.children;return l.createElement(l.Fragment,{},n)}},m=l.forwardRef((function(e,n){var t=e.components,i=e.mdxType,r=e.originalType,o=e.parentName,u=s(e,["components","mdxType","originalType","parentName"]),m=p(t),d=i,g=m["".concat(o,".").concat(d)]||m[d]||c[d]||r;return t?l.createElement(g,a(a({ref:n},u),{},{components:t})):l.createElement(g,a({ref:n},u))}));function d(e,n){var t=arguments,i=n&&n.mdxType;if("string"==typeof e||i){var r=t.length,a=new Array(r);a[0]=m;var s={};for(var o in n)hasOwnProperty.call(n,o)&&(s[o]=n[o]);s.originalType=e,s.mdxType="string"==typeof e?e:i,a[1]=s;for(var p=2;p<r;p++)a[p]=t[p];return l.createElement.apply(null,a)}return l.createElement.apply(null,t)}m.displayName="MDXCreateElement"},381:function(e,n,t){t.r(n),t.d(n,{frontMatter:function(){return s},contentTitle:function(){return o},metadata:function(){return p},toc:function(){return u},default:function(){return m}});var l=t(7462),i=t(3366),r=(t(7294),t(3905)),a=["components"],s={id:"genotypeutils",title:"GenotypeUtils.jl",sidebar_label:"GenotypeUtils.jl"},o=void 0,p={unversionedId:"api/PopGenCore/genotypeutils",id:"api/PopGenCore/genotypeutils",isDocsHomePage:!1,title:"GenotypeUtils.jl",description:"PopGenCore.jl/src/Utils/GenotypeUtils.jl",source:"@site/docs/api/PopGenCore/GenotypeUtils.md",sourceDirName:"api/PopGenCore",slug:"/api/PopGenCore/genotypeutils",permalink:"/PopGen.jl/docs/api/PopGenCore/genotypeutils",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGenCore/GenotypeUtils.md",tags:[],version:"current",lastUpdatedAt:1636029729,formattedLastUpdatedAt:"11/4/2021",frontMatter:{id:"genotypeutils",title:"GenotypeUtils.jl",sidebar_label:"GenotypeUtils.jl"},sidebar:"docs",previous:{title:"GenoFreq.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/genofreq"},next:{title:"ioUtils.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/ioutils"}},u=[{value:"PopGenCore.jl/src/Utils/GenotypeUtils.jl",id:"popgencorejlsrcutilsgenotypeutilsjl",children:[{value:"\ud83d\udfea alleles",id:"-alleles",children:[],level:3},{value:"\ud83d\udfea alleles",id:"-alleles-1",children:[],level:3},{value:"\ud83d\udfea alleles",id:"-alleles-2",children:[],level:3},{value:"\ud83d\udfea uniquealleles",id:"-uniquealleles",children:[],level:3},{value:"\ud83d\udfea locidataframe",id:"-locidataframe",children:[],level:3},{value:"\ud83d\udfea locimatrix",id:"-locimatrix",children:[],level:3},{value:"\ud83d\udfea phasedmatrix",id:"-phasedmatrix",children:[],level:3}],level:2}],c={toc:u};function m(e){var n=e.components,t=(0,i.Z)(e,a);return(0,r.kt)("wrapper",(0,l.Z)({},c,t,{components:n,mdxType:"MDXLayout"}),(0,r.kt)("h2",{id:"popgencorejlsrcutilsgenotypeutilsjl"},"PopGenCore.jl/src/Utils/GenotypeUtils.jl"),(0,r.kt)("p",null,"\ud83d\udce6  => not exported |\n\ud83d\udfea => exported by PopGenCore.jl |\n\ud83d\udd35 => exported by PopGen.jl"),(0,r.kt)("h3",{id:"-alleles"},"\ud83d\udfea alleles"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"allelecount(locus::T) where T<:GenoArray\n")),(0,r.kt)("p",null,"Return the number of unique alleles present at a locus."),(0,r.kt)("hr",null),(0,r.kt)("h3",{id:"-alleles-1"},"\ud83d\udfea alleles"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"alleles(locus::T) where T<:GenoArray\n")),(0,r.kt)("p",null,"Return an array of all the non-missing alleles of a locus."),(0,r.kt)("hr",null),(0,r.kt)("h3",{id:"-alleles-2"},"\ud83d\udfea alleles"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"alleles(locus::T, miss::Bool = false) where T<:GenoArray\n")),(0,r.kt)("p",null,"Return an array of all the non-missing alleles of a locus. Use the second positional\nargument as ",(0,r.kt)("inlineCode",{parentName:"p"},"true")," to include missing values."),(0,r.kt)("hr",null),(0,r.kt)("h3",{id:"-uniquealleles"},"\ud83d\udfea uniquealleles"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"uniquealleles(locus::T) where T<:GenoArray\n")),(0,r.kt)("p",null,"Return an array of all the unique non-missing alleles of a locus."),(0,r.kt)("hr",null),(0,r.kt)("h3",{id:"-locidataframe"},"\ud83d\udfea locidataframe"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"locidataframe(data::PopData)\n")),(0,r.kt)("p",null,"Return a wide ",(0,r.kt)("inlineCode",{parentName:"p"},"DataFrame")," of samples as columns, ommitting population information."),(0,r.kt)("p",null,(0,r.kt)("strong",{parentName:"p"},"Example")),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"julia> locidataframe(@nancycats)\n9\xd7237 DataFrame. Omitted printing of 232 columns\n\u2502 Row \u2502 N215       \u2502 N216       \u2502 N217       \u2502 N218       \u2502 N219       \u2502\n\u2502     \u2502 Tuple\u2026?    \u2502 Tuple\u2026?    \u2502 Tuple\u2026?    \u2502 Tuple\u2026?    \u2502 Tuple\u2026?    \u2502\n\u251c\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524\n\u2502 1   \u2502 missing    \u2502 missing    \u2502 (135, 143) \u2502 (133, 135) \u2502 (133, 135) \u2502\n\u2502 2   \u2502 (136, 146) \u2502 (146, 146) \u2502 (136, 146) \u2502 (138, 138) \u2502 (140, 146) \u2502\n\u2502 3   \u2502 (139, 139) \u2502 (139, 145) \u2502 (141, 141) \u2502 (139, 141) \u2502 (141, 145) \u2502\n\u2502 4   \u2502 (116, 120) \u2502 (120, 126) \u2502 (116, 116) \u2502 (116, 126) \u2502 (126, 126) \u2502\n\u2502 5   \u2502 (156, 156) \u2502 (156, 156) \u2502 (152, 156) \u2502 (150, 150) \u2502 (152, 152) \u2502\n\u2502 6   \u2502 (142, 148) \u2502 (142, 148) \u2502 (142, 142) \u2502 (142, 148) \u2502 (142, 148) \u2502\n\u2502 7   \u2502 (199, 199) \u2502 (185, 199) \u2502 (197, 197) \u2502 (199, 199) \u2502 (193, 199) \u2502\n\u2502 8   \u2502 (113, 113) \u2502 (113, 113) \u2502 (113, 113) \u2502 (91, 105)  \u2502 (113, 113) \u2502\n\u2502 9   \u2502 (208, 208) \u2502 (208, 208) \u2502 (210, 210) \u2502 (208, 208) \u2502 (208, 208) \u2502\n")),(0,r.kt)("hr",null),(0,r.kt)("h3",{id:"-locimatrix"},"\ud83d\udfea locimatrix"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"locimatrix(data::PopData)\n")),(0,r.kt)("p",null,"Return a matrix of genotypes with dimensions ",(0,r.kt)("inlineCode",{parentName:"p"},"samples \xd7 loci"),".\nRows are samples and columns are loci. Will return an error if ploidy varies between samples."),(0,r.kt)("p",null,(0,r.kt)("strong",{parentName:"p"},"Example")),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"julia> locimatrix(@nancycats)\n237\xd79 Array{Union{Missing, Tuple{Int16,Int16}},2}:\n missing     (136, 146)  (139, 139)  \u2026  (199, 199)  (113, 113)  (208, 208)\n missing     (146, 146)  (139, 145)     (185, 199)  (113, 113)  (208, 208)\n (135, 143)  (136, 146)  (141, 141)     (197, 197)  (113, 113)  (210, 210)\n (133, 135)  (138, 138)  (139, 141)     (199, 199)  (91, 105)   (208, 208)\n (133, 135)  (140, 146)  (141, 145)     (193, 199)  (113, 113)  (208, 208)\n (135, 143)  (136, 146)  (145, 149)  \u2026  (193, 195)  (91, 113)   (208, 208)\n (135, 135)  (136, 146)  (139, 145)     (199, 199)  (105, 113)  (208, 208)\n (135, 143)  (136, 146)  (135, 149)     (193, 197)  (91, 91)    (208, 212)\n (137, 143)  (136, 146)  (139, 139)     (197, 197)  (105, 113)  (208, 212)\n (135, 135)  (132, 132)  (141, 145)     (197, 197)  (91, 105)   (208, 208)\n (137, 141)  (130, 136)  (137, 145)  \u2026  (193, 199)  (91, 91)    (182, 182)\n (129, 133)  (130, 136)  (135, 145)     (193, 199)  (91, 113)   (182, 208)\n \u22ee                                   \u22f1                          \n (133, 135)  (136, 136)  (135, 139)  \u2026  (199, 199)  (113, 113)  (182, 182)\n (133, 141)  (136, 136)  (135, 139)     (197, 197)  (113, 113)  (182, 208)\n (133, 141)  (130, 146)  (141, 141)     (191, 199)  missing     (208, 208)\n (123, 133)  (138, 138)  (141, 145)     (191, 197)  missing     (208, 208)\n (123, 133)  (138, 138)  (139, 139)     (197, 199)  missing     (208, 208)\n (133, 141)  (136, 146)  (139, 139)  \u2026  (197, 197)  missing     (208, 208)\n (133, 141)  (130, 136)  (139, 145)     (191, 199)  missing     (208, 208)\n (133, 141)  (136, 146)  (139, 145)     (199, 199)  missing     (208, 220)\n (133, 143)  (130, 130)  (135, 145)     (197, 197)  missing     (208, 208)\n (135, 141)  (136, 144)  (143, 143)     (191, 197)  (113, 117)  (208, 208)\n (137, 143)  (130, 136)  (135, 145)  \u2026  (193, 199)  (113, 117)  (208, 208)\n (135, 141)  (130, 146)  (135, 139)     (197, 197)  missing     (208, 208)\n")),(0,r.kt)("hr",null),(0,r.kt)("h3",{id:"-phasedmatrix"},"\ud83d\udfea phasedmatrix"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"phasedmatrix(data::PopData)\n")),(0,r.kt)("p",null,"Return a ",(0,r.kt)("inlineCode",{parentName:"p"},"Vector")," of length ",(0,r.kt)("inlineCode",{parentName:"p"},"ploidy")," composed of allele matrices with dimensions ",(0,r.kt)("inlineCode",{parentName:"p"},"samples \xd7 loci"),".\nRows are samples and columns are loci. Will return an error if ploidy varies between samples. "),(0,r.kt)("p",null,(0,r.kt)("strong",{parentName:"p"},"Example")),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"julia> mtx = phasedmatrix(@nancycats)\n2-element Array{Array{Union{Missing, Int16},2},1}:\n [missing 136 \u2026 113 208; missing 146 \u2026 113 208; \u2026 ; 137 130 \u2026 113 208; 135 130 \u2026 missing 208]\n [missing 146 \u2026 113 208; missing 146 \u2026 113 208; \u2026 ; 143 136 \u2026 117 208; 141 146 \u2026 missing 208]\njulia> mtx[1]\n237\xd79 Array{Union{Missing, Int16},2}:\n    missing  136  139  116         156  142  199  113         208\n    missing  146  139  120         156  142  185  113         208\n 135         136  141  116         152  142  197  113         210\n 133         138  139  116         150  142  199   91         208\n 133         140  141  126         152  142  193  113         208\n 135         136  145  120         150  148  193   91         208\n 135         136  139  116         152  142  199  105         208\n 135         136  135  120         154  142  193   91         208\n 137         136  139  116         150  142  197  105         208\n 135         132  141  120         150  148  197   91         208\n 137         130  137  128         152  142  193   91         182\n 129         130  135  126         144  140  193   91         182\n   \u22ee                                      \u22ee                   \n 133         136  135     missing  146  142  199  113         182\n 133         136  135     missing  150  142  197  113         182\n 133         130  141     missing  148  142  191     missing  208\n 123         138  141     missing  148  142  191     missing  208\n 123         138  139     missing  150  142  197     missing  208\n 133         136  139     missing  150  142  197     missing  208\n 133         130  139     missing  152  142  191     missing  208\n 133         136  139     missing  150  142  199     missing  208\n 133         130  135     missing  148  142  197     missing  208\n 135         136  143     missing  144  142  191  113         208\n 137         130  135     missing  150  142  193  113         208\n 135         130  135     missing  150  142  197     missing  208\n")))}m.isMDXComponent=!0}}]);