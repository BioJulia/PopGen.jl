"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[7578],{3905:function(e,t,n){n.d(t,{Zo:function(){return s},kt:function(){return m}});var r=n(7294);function a(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function o(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);t&&(r=r.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,r)}return n}function i(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?o(Object(n),!0).forEach((function(t){a(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):o(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,r,a=function(e,t){if(null==e)return{};var n,r,a={},o=Object.keys(e);for(r=0;r<o.length;r++)n=o[r],t.indexOf(n)>=0||(a[n]=e[n]);return a}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(r=0;r<o.length;r++)n=o[r],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(a[n]=e[n])}return a}var p=r.createContext({}),u=function(e){var t=r.useContext(p),n=t;return e&&(n="function"==typeof e?e(t):i(i({},t),e)),n},s=function(e){var t=u(e.components);return r.createElement(p.Provider,{value:t},e.children)},c={inlineCode:"code",wrapper:function(e){var t=e.children;return r.createElement(r.Fragment,{},t)}},d=r.forwardRef((function(e,t){var n=e.components,a=e.mdxType,o=e.originalType,p=e.parentName,s=l(e,["components","mdxType","originalType","parentName"]),d=u(n),m=a,f=d["".concat(p,".").concat(m)]||d[m]||c[m]||o;return n?r.createElement(f,i(i({ref:t},s),{},{components:n})):r.createElement(f,i({ref:t},s))}));function m(e,t){var n=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var o=n.length,i=new Array(o);i[0]=d;var l={};for(var p in t)hasOwnProperty.call(t,p)&&(l[p]=t[p]);l.originalType=e,l.mdxType="string"==typeof e?e:a,i[1]=l;for(var u=2;u<o;u++)i[u]=n[u];return r.createElement.apply(null,i)}return r.createElement.apply(null,n)}d.displayName="MDXCreateElement"},4707:function(e,t,n){n.r(t),n.d(t,{frontMatter:function(){return l},contentTitle:function(){return p},metadata:function(){return u},toc:function(){return s},default:function(){return d}});var r=n(7462),a=n(3366),o=(n(7294),n(3905)),i=["components"],l={id:"hardyweinberg",title:"HardyWeinberg.jl",sidebar_label:"HardyWeinberg.jl"},p=void 0,u={unversionedId:"api/PopGen/hardyweinberg",id:"api/PopGen/hardyweinberg",isDocsHomePage:!1,title:"HardyWeinberg.jl",description:"PopGen.jl/src/HardyWeinberg.jl",source:"@site/docs/api/PopGen/HardyWeinberg.md",sourceDirName:"api/PopGen",slug:"/api/PopGen/hardyweinberg",permalink:"/PopGen.jl/docs/api/PopGen/hardyweinberg",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGen/HardyWeinberg.md",tags:[],version:"current",lastUpdatedAt:1636029729,formattedLastUpdatedAt:"11/4/2021",frontMatter:{id:"hardyweinberg",title:"HardyWeinberg.jl",sidebar_label:"HardyWeinberg.jl"},sidebar:"docs",previous:{title:"FstPermutations.jl",permalink:"/PopGen.jl/docs/api/PopGen/fstpermutations"},next:{title:"Heterozygosity.jl",permalink:"/PopGen.jl/docs/api/PopGen/heterozygosity"}},s=[{value:"PopGen.jl/src/HardyWeinberg.jl",id:"popgenjlsrchardyweinbergjl",children:[{value:"\ud83d\udce6 _chisqlocus",id:"-_chisqlocus",children:[],level:3},{value:"\ud83d\udd35 hwetest",id:"-hwetest",children:[],level:3}],level:2}],c={toc:s};function d(e){var t=e.components,n=(0,a.Z)(e,i);return(0,o.kt)("wrapper",(0,r.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,o.kt)("h2",{id:"popgenjlsrchardyweinbergjl"},"PopGen.jl/src/HardyWeinberg.jl"),(0,o.kt)("p",null,"\ud83d\udce6  => not exported |\n\ud83d\udd35 => exported by PopGen.jl"),(0,o.kt)("h3",{id:"-_chisqlocus"},"\ud83d\udce6 _chisqlocus"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"_chisqlocus(locus::T) where T <: GenotypeArray\n")),(0,o.kt)("p",null,"Calculate the chi square statistic and p-value for a locus\nReturns a tuple with chi-square statistic, degrees of freedom, and p-value."),(0,o.kt)("hr",null),(0,o.kt)("h3",{id:"-hwetest"},"\ud83d\udd35 hwetest"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'hwetest(data::PopData; by_pop::Bool = false; correction = "none")\n')),(0,o.kt)("p",null,"Calculate chi-squared test of HWE for each locus and returns observed and\nexpected heterozygosity with chi-squared, degrees of freedom and p-values for each locus. Use ",(0,o.kt)("inlineCode",{parentName:"p"},"by_pop = true")," to perform this separately for each population (default: by_pop = false) and return a NamedTuple of DataFrames with the names corresponding to the population names. Use ",(0,o.kt)("inlineCode",{parentName:"p"},"correction =")," to specify a P-value adjustment method for multiple testing."),(0,o.kt)("p",null,(0,o.kt)("strong",{parentName:"p"},"correction methods (case insensitive)")),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"bonferroni"')," : Bonferroni adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"holm"')," : Holm adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"hochberg"')," : Hochberg adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"bh"')," : Benjamini-Hochberg adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"by"')," : Benjamini-Yekutieli adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"bl"'),"  : Benjamini-Liu adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"hommel"')," : Hommel adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"sidak"')," : \u0160id\xe1k adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"forwardstop"')," or ",(0,o.kt)("inlineCode",{parentName:"li"},'"fs"')," : Forward-Stop adjustment"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},'"bc"')," : Barber-Cande\u0300s adjustment")),(0,o.kt)("p",null,(0,o.kt)("strong",{parentName:"p"},"Example")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'hwetest(@gulfsharks, correction = "bh")\nhwetest(@gulfsharks, by_pop = true, correction = "bh")\n')))}d.isMDXComponent=!0}}]);