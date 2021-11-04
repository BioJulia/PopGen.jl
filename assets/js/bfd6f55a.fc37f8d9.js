"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[4952],{3905:function(e,t,n){n.d(t,{Zo:function(){return s},kt:function(){return m}});var r=n(7294);function a(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function l(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);t&&(r=r.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,r)}return n}function i(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?l(Object(n),!0).forEach((function(t){a(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):l(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function o(e,t){if(null==e)return{};var n,r,a=function(e,t){if(null==e)return{};var n,r,a={},l=Object.keys(e);for(r=0;r<l.length;r++)n=l[r],t.indexOf(n)>=0||(a[n]=e[n]);return a}(e,t);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);for(r=0;r<l.length;r++)n=l[r],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(a[n]=e[n])}return a}var p=r.createContext({}),u=function(e){var t=r.useContext(p),n=t;return e&&(n="function"==typeof e?e(t):i(i({},t),e)),n},s=function(e){var t=u(e.components);return r.createElement(p.Provider,{value:t},e.children)},c={inlineCode:"code",wrapper:function(e){var t=e.children;return r.createElement(r.Fragment,{},t)}},d=r.forwardRef((function(e,t){var n=e.components,a=e.mdxType,l=e.originalType,p=e.parentName,s=o(e,["components","mdxType","originalType","parentName"]),d=u(n),m=a,f=d["".concat(p,".").concat(m)]||d[m]||c[m]||l;return n?r.createElement(f,i(i({ref:t},s),{},{components:n})):r.createElement(f,i({ref:t},s))}));function m(e,t){var n=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var l=n.length,i=new Array(l);i[0]=d;var o={};for(var p in t)hasOwnProperty.call(t,p)&&(o[p]=t[p]);o.originalType=e,o.mdxType="string"==typeof e?e:a,i[1]=o;for(var u=2;u<l;u++)i[u]=n[u];return r.createElement.apply(null,i)}return r.createElement.apply(null,n)}d.displayName="MDXCreateElement"},9833:function(e,t,n){n.r(t),n.d(t,{frontMatter:function(){return o},contentTitle:function(){return p},metadata:function(){return u},toc:function(){return s},default:function(){return d}});var r=n(7462),a=n(3366),l=(n(7294),n(3905)),i=["components"],o={id:"utils",title:"Utils.jl",sidebar_label:"Utils.jl"},p=void 0,u={unversionedId:"api/PopGen/utils",id:"api/PopGen/utils",isDocsHomePage:!1,title:"Utils.jl",description:"PopGen.jl/src/Utils.jl",source:"@site/docs/api/PopGen/Utils.md",sourceDirName:"api/PopGen",slug:"/api/PopGen/utils",permalink:"/PopGen.jl/docs/api/PopGen/utils",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGen/Utils.md",tags:[],version:"current",lastUpdatedAt:1636029729,formattedLastUpdatedAt:"11/4/2021",frontMatter:{id:"utils",title:"Utils.jl",sidebar_label:"Utils.jl"},sidebar:"docs",previous:{title:"SummaryInfo.jl",permalink:"/PopGen.jl/docs/api/PopGen/summaryinfo"},next:{title:"AlleleFreq.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/allelefreq"}},s=[{value:"PopGen.jl/src/Utils.jl",id:"popgenjlsrcutilsjl",children:[{value:"\ud83d\udce6 _adjacency_matrix",id:"-_adjacency_matrix",children:[],level:3},{value:"\ud83d\udce6 _p_adjust",id:"-_p_adjust",children:[],level:3},{value:"\ud83d\udce6 feature_req",id:"-feature_req",children:[],level:3}],level:2}],c={toc:s};function d(e){var t=e.components,n=(0,a.Z)(e,i);return(0,l.kt)("wrapper",(0,r.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,l.kt)("h2",{id:"popgenjlsrcutilsjl"},"PopGen.jl/src/Utils.jl"),(0,l.kt)("p",null,"\ud83d\udce6  => not exported |\n\ud83d\udd35 => exported by PopGen.jl"),(0,l.kt)("h3",{id:"-_adjacency_matrix"},"\ud83d\udce6 _adjacency_matrix"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"_adjacency_matrix(data::PopData)\n")),(0,l.kt)("h3",{id:"-_p_adjust"},"\ud83d\udce6 _p_adjust"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"_p_adjust(pvals::Vector{T}, method::String) where T <: Union{Missing, <:AbstractFloat}\n")),(0,l.kt)("p",null,"Modification to ",(0,l.kt)("inlineCode",{parentName:"p"},"MultipleTesting.adjust")," to include ",(0,l.kt)("inlineCode",{parentName:"p"},"missing")," values in the\nreturned array. See MultipleTesting.jl docs for full more detailed information."),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},"Example")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},'julia> _p_adjust([0.1, 0.01, 0.005, 0.3], "bh")\n')),(0,l.kt)("p",null,(0,l.kt)("strong",{parentName:"p"},(0,l.kt)("inlineCode",{parentName:"strong"},"correction")," methods (case insensitive)")),(0,l.kt)("ul",null,(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"bonferroni"')," : Bonferroni adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"holm"')," : Holm adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"hochberg"')," : Hochberg adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"bh"')," : Benjamini-Hochberg adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"by"')," : Benjamini-Yekutieli adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"bl"')," : Benjamini-Liu adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"hommel"')," : Hommel adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"sidak"')," : \u0160id\xe1k adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"forwardstop"')," or ",(0,l.kt)("inlineCode",{parentName:"li"},'"fs"')," : Forward-Stop adjustment"),(0,l.kt)("li",{parentName:"ul"},(0,l.kt)("inlineCode",{parentName:"li"},'"bc"'),' : Barber-Cande\u0300s adjustment\n"""')),(0,l.kt)("h3",{id:"-feature_req"},"\ud83d\udce6 feature_req"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"feature_req()\n")),(0,l.kt)("p",null,"Returns the text: ",(0,l.kt)("inlineCode",{parentName:"p"},"Please open an Issue or Pull Request on https://www.github.com/biojulia/PopGen.jl if you would like this feature implemented")))}d.isMDXComponent=!0}}]);