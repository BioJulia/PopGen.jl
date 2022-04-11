"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[7356],{3905:function(e,n,t){t.d(n,{Zo:function(){return p},kt:function(){return d}});var a=t(7294);function r(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function i(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);n&&(a=a.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,a)}return t}function s(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?i(Object(t),!0).forEach((function(n){r(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):i(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function o(e,n){if(null==e)return{};var t,a,r=function(e,n){if(null==e)return{};var t,a,r={},i=Object.keys(e);for(a=0;a<i.length;a++)t=i[a],n.indexOf(t)>=0||(r[t]=e[t]);return r}(e,n);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(a=0;a<i.length;a++)t=i[a],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(r[t]=e[t])}return r}var l=a.createContext({}),u=function(e){var n=a.useContext(l),t=n;return e&&(t="function"==typeof e?e(n):s(s({},n),e)),t},p=function(e){var n=u(e.components);return a.createElement(l.Provider,{value:n},e.children)},c={inlineCode:"code",wrapper:function(e){var n=e.children;return a.createElement(a.Fragment,{},n)}},m=a.forwardRef((function(e,n){var t=e.components,r=e.mdxType,i=e.originalType,l=e.parentName,p=o(e,["components","mdxType","originalType","parentName"]),m=u(t),d=r,f=m["".concat(l,".").concat(d)]||m[d]||c[d]||i;return t?a.createElement(f,s(s({ref:n},p),{},{components:t})):a.createElement(f,s({ref:n},p))}));function d(e,n){var t=arguments,r=n&&n.mdxType;if("string"==typeof e||r){var i=t.length,s=new Array(i);s[0]=m;var o={};for(var l in n)hasOwnProperty.call(n,l)&&(o[l]=n[l]);o.originalType=e,o.mdxType="string"==typeof e?e:r,s[1]=o;for(var u=2;u<i;u++)s[u]=t[u];return a.createElement.apply(null,s)}return a.createElement.apply(null,t)}m.displayName="MDXCreateElement"},7101:function(e,n,t){t.r(n),t.d(n,{frontMatter:function(){return o},contentTitle:function(){return l},metadata:function(){return u},assets:function(){return p},toc:function(){return c},default:function(){return d}});var a=t(7462),r=t(3366),i=(t(7294),t(3905)),s=["components"],o={id:"kmeans",title:"K-Means Clustering",sidebar_label:"K-Means Clustering"},l=void 0,u={unversionedId:"analyses/kmeans",id:"analyses/kmeans",title:"K-Means Clustering",description:"Usually the beginning of a study without prior population information requires guesstimating the number of clusters present in the data.",source:"@site/docs/analyses/kmeans.md",sourceDirName:"analyses",slug:"/analyses/kmeans",permalink:"/PopGen.jl/docs/analyses/kmeans",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/analyses/kmeans.md",tags:[],version:"current",lastUpdatedAt:1649705603,formattedLastUpdatedAt:"4/11/2022",frontMatter:{id:"kmeans",title:"K-Means Clustering",sidebar_label:"K-Means Clustering"},sidebar:"docs",previous:{title:"Hardy-Weinberg Equilibrium",permalink:"/PopGen.jl/docs/analyses/hardyweinberg"},next:{title:"Pairwise F-Statistics",permalink:"/PopGen.jl/docs/analyses/fstatistics"}},p={},c=[{value:"<code>kmeans</code>",id:"kmeans",level:3}],m={toc:c};function d(e){var n=e.components,t=(0,r.Z)(e,s);return(0,i.kt)("wrapper",(0,a.Z)({},m,t,{components:n,mdxType:"MDXLayout"}),(0,i.kt)("p",null,"Usually the beginning of a study without prior population information requires guesstimating the number of clusters present in the data.\nOne way to accomplish this is with K-means clustering, an unsupervised clustering algorithm that clusters data based on similarity. The\nPopGen.jl implementation of K-means clustering uses the K-means ++ algorithm (",(0,i.kt)("a",{parentName:"p",href:"http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf"},"Arthur & Vassilvitskii 2007"),") under the hood, as implemented in ",(0,i.kt)("inlineCode",{parentName:"p"},"Clustering.jl")," (",(0,i.kt)("a",{parentName:"p",href:"https://github.com/JuliaStats/Clustering.jl"},"link"),")."),(0,i.kt)("h3",{id:"kmeans"},(0,i.kt)("inlineCode",{parentName:"h3"},"kmeans")),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"# krange as a vector : [1,3,6,7]\nkmeans(data::PopData; krange::Vector{Int64}, iterations::Int64 = 100)\n# krange as a range: 2:9\nkmeans(data::PopData; krange::UnitRange{Int64}, iterations::Int64 = 100)\n# without keywords\nkmeans(data::PopData, krange::Union{UnitRange{Int64},Vector{Int64}}, iterations::Int64 = 100)\n")),(0,i.kt)("p",null,"Perform PCA-based Kmeans clustering (using Kmeans++) on a ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," object. Returns a ",(0,i.kt)("inlineCode",{parentName:"p"},"KMeansResults"),"\nobject storing the results of each ",(0,i.kt)("inlineCode",{parentName:"p"},"kmeans()")," run across a range of K values (",(0,i.kt)("inlineCode",{parentName:"p"},"krange"),", default: ",(0,i.kt)("inlineCode",{parentName:"p"},"2:nsamples/10"),").\nUse the keyword argument ",(0,i.kt)("inlineCode",{parentName:"p"},"iterations")," (default: 100) to set the maximum number of iterations allowed to\nachieve convergence. Interally, kmeans clustering is performed on the principal components of the scaled allele frequencies\nmatrix, where ",(0,i.kt)("inlineCode",{parentName:"p"},"missing")," values are replaced by the global mean allele frequency."),(0,i.kt)("p",null,(0,i.kt)("strong",{parentName:"p"},"Example")),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"julia> cats = @nancycats ;\njulia> km = kmeans(cats, krange = 2:7)\nK-Means(++) Clustering Results\n  K values:          2, 3, 4, 5, 6, 7\n  Iterations per K:  2, 4, 6, 5, 5, 7\n  Convergence per K: T, T, T, T, T, T\n  Total Cost per K:  91.41, 90.679, 89.54, 89.287, 88.284, 87.957\n  Available fields to inspect: assignments, costs, centers, other\n  \njulia> km.assignments\n237\xd77 DataFrame\n Row \u2502 name     2      3      4      5      6      7     \n     \u2502 String   Int64  Int64  Int64  Int64  Int64  Int64\n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 N215         1      1      1      2      3      3\n   2 \u2502 N216         1      2      1      2      3      6\n   3 \u2502 N217         1      3      1      1      3      3\n   4 \u2502 N218         1      3      1      1      3      6\n   5 \u2502 N219         1      2      1      2      3      3\n   6 \u2502 N220         1      2      1      2      3      1\n  \u22ee  \u2502    \u22ee       \u22ee      \u22ee      \u22ee      \u22ee      \u22ee      \u22ee\n 233 \u2502 N296         1      1      1      2      3      6\n 234 \u2502 N297         1      1      3      2      4      3\n 235 \u2502 N281         1      2      1      2      3      3\n 236 \u2502 N289         1      3      3      2      3      4\n 237 \u2502 N290         1      1      1      2      3      6\n")))}d.isMDXComponent=!0}}]);