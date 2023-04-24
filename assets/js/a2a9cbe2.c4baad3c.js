"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[6633],{3905:function(e,t,n){n.d(t,{Zo:function(){return m},kt:function(){return f}});var a=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function i(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function l(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?i(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):i(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function o(e,t){if(null==e)return{};var n,a,r=function(e,t){if(null==e)return{};var n,a,r={},i=Object.keys(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var p=a.createContext({}),s=function(e){var t=a.useContext(p),n=t;return e&&(n="function"==typeof e?e(t):l(l({},t),e)),n},m=function(e){var t=s(e.components);return a.createElement(p.Provider,{value:t},e.children)},u="mdxType",d={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},c=a.forwardRef((function(e,t){var n=e.components,r=e.mdxType,i=e.originalType,p=e.parentName,m=o(e,["components","mdxType","originalType","parentName"]),u=s(n),c=r,f=u["".concat(p,".").concat(c)]||u[c]||d[c]||i;return n?a.createElement(f,l(l({ref:t},m),{},{components:n})):a.createElement(f,l({ref:t},m))}));function f(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var i=n.length,l=new Array(i);l[0]=c;var o={};for(var p in t)hasOwnProperty.call(t,p)&&(o[p]=t[p]);o.originalType=e,o[u]="string"==typeof e?e:r,l[1]=o;for(var s=2;s<i;s++)l[s]=n[s];return a.createElement.apply(null,l)}return a.createElement.apply(null,n)}c.displayName="MDXCreateElement"},7296:function(e,t,n){n.r(t),n.d(t,{assets:function(){return m},contentTitle:function(){return p},default:function(){return f},frontMatter:function(){return o},metadata:function(){return s},toc:function(){return u}});var a=n(7462),r=n(3366),i=(n(7294),n(3905)),l=["components"],o={id:"tsne",title:"TSNE.jl",sidebar_label:"TSNE.jl"},p=void 0,s={unversionedId:"api/PopGen/tsne",id:"api/PopGen/tsne",title:"TSNE.jl",description:"PopGen.jl/src/SummaryInfo.jl",source:"@site/docs/api/PopGen/TSNE.md",sourceDirName:"api/PopGen",slug:"/api/PopGen/tsne",permalink:"/PopGen.jl/docs/api/PopGen/tsne",draft:!1,editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGen/TSNE.md",tags:[],version:"current",lastUpdatedAt:1658766707,formattedLastUpdatedAt:"Jul 25, 2022",frontMatter:{id:"tsne",title:"TSNE.jl",sidebar_label:"TSNE.jl"}},m={},u=[{value:"PopGen.jl/src/SummaryInfo.jl",id:"popgenjlsrcsummaryinfojl",level:2},{value:"\ud83d\udd35 tsne",id:"-tsne",level:3}],d={toc:u},c="wrapper";function f(e){var t=e.components,n=(0,r.Z)(e,l);return(0,i.kt)(c,(0,a.Z)({},d,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("h2",{id:"popgenjlsrcsummaryinfojl"},"PopGen.jl/src/SummaryInfo.jl"),(0,i.kt)("table",null,(0,i.kt)("thead",{parentName:"table"},(0,i.kt)("tr",{parentName:"thead"},(0,i.kt)("th",{parentName:"tr",align:"center"},"\ud83d\udce6  not exported"),(0,i.kt)("th",{parentName:"tr",align:"center"},"\ud83d\udd35  exported by PopGen.jl")))),(0,i.kt)("h3",{id:"-tsne"},"\ud83d\udd35 tsne"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"tsne(data::PopData, args...; kwargs...)\n")),(0,i.kt)("p",null,"Perform t-SNE (t-Stochastic Neighbor Embedding) on a PopData object, returning a DataFrame. Converts the\nPopData object into a matrix of allele frequencies with missing values replaced with\nthe global mean frequency of that allele. First performs PCA on that matrix, retaining\n",(0,i.kt)("inlineCode",{parentName:"p"},"reduce_dims")," dimensions of the PCA prior to t-SNE analysis. The positional and keyword arguments\nare the same as ",(0,i.kt)("inlineCode",{parentName:"p"},"tsne")," from ",(0,i.kt)("inlineCode",{parentName:"p"},"TSne.jl"),"."),(0,i.kt)("p",null,(0,i.kt)("strong",{parentName:"p"},"Arguments")),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"data"),": a ",(0,i.kt)("inlineCode",{parentName:"li"},"PopData")," object"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"ndims"),": Dimension of the embedded space (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"2"),")"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"reduce_dims")," the number of the first dimensions of X PCA to use for t-SNE, if 0, all available dimension are used (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"0"),")"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"max_iter"),": Maximum number of iterations for the optimization (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"1000"),")"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"perplexity"),": The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. Different values can result in significantly different results (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"30"),")")),(0,i.kt)("p",null,(0,i.kt)("strong",{parentName:"p"},"Keyword Arguments (optional)")),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"distance"),": type ",(0,i.kt)("inlineCode",{parentName:"li"},"Function")," or ",(0,i.kt)("inlineCode",{parentName:"li"},"Distances.SemiMetric"),", specifies the function to\nuse for calculating the distances between the rows"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"pca_init"),": whether to use the first ",(0,i.kt)("inlineCode",{parentName:"li"},"ndims")," of the PCA as the initial t-SNE layout,\nif ",(0,i.kt)("inlineCode",{parentName:"li"},"false")," (the default), the method is initialized with the random layout"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"max_iter"),": how many iterations of t-SNE to do"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"verbose"),": output informational and diagnostic messages"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"progress"),": display progress meter during t-SNE optimization (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"true"),")"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"min_gain"),": ",(0,i.kt)("inlineCode",{parentName:"li"},"eta"),": ",(0,i.kt)("inlineCode",{parentName:"li"},"initial_momentum"),", ",(0,i.kt)("inlineCode",{parentName:"li"},"final_momentum"),", ",(0,i.kt)("inlineCode",{parentName:"li"},"momentum_switch_iter"),",\n",(0,i.kt)("inlineCode",{parentName:"li"},"stop_cheat_iter"),": ",(0,i.kt)("inlineCode",{parentName:"li"},"cheat_scale")," low-level parameters of t-SNE optimization"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"extended_output"),": if ",(0,i.kt)("inlineCode",{parentName:"li"},"true"),", returns a tuple of embedded coordinates matrix,\npoint perplexities and final Kullback-Leibler divergence")))}f.isMDXComponent=!0}}]);