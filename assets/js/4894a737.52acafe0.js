"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[6116],{3905:function(e,t,n){n.d(t,{Zo:function(){return c},kt:function(){return h}});var a=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function i(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?i(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):i(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function s(e,t){if(null==e)return{};var n,a,r=function(e,t){if(null==e)return{};var n,a,r={},i=Object.keys(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var l=a.createContext({}),p=function(e){var t=a.useContext(l),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},c=function(e){var t=p(e.components);return a.createElement(l.Provider,{value:t},e.children)},d={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},u=a.forwardRef((function(e,t){var n=e.components,r=e.mdxType,i=e.originalType,l=e.parentName,c=s(e,["components","mdxType","originalType","parentName"]),u=p(n),h=r,y=u["".concat(l,".").concat(h)]||u[h]||d[h]||i;return n?a.createElement(y,o(o({ref:t},c),{},{components:n})):a.createElement(y,o({ref:t},c))}));function h(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var i=n.length,o=new Array(i);o[0]=u;var s={};for(var l in t)hasOwnProperty.call(t,l)&&(s[l]=t[l]);s.originalType=e,s.mdxType="string"==typeof e?e:r,o[1]=s;for(var p=2;p<i;p++)o[p]=n[p];return a.createElement.apply(null,o)}return a.createElement.apply(null,n)}u.displayName="MDXCreateElement"},6362:function(e,t,n){n.r(t),n.d(t,{frontMatter:function(){return s},contentTitle:function(){return l},metadata:function(){return p},toc:function(){return c},default:function(){return u}});var a=n(7462),r=n(3366),i=(n(7294),n(3905)),o=["components"],s={id:"othertypes",title:"Other data types",sidebar_label:"Other data types"},l=void 0,p={unversionedId:"gettingstarted/othertypes",id:"gettingstarted/othertypes",isDocsHomePage:!1,title:"Other data types",description:"While not strictly their own composite types, we also define aliases for genotypes and vectors of genotypes, as their explicit types can get a little unwieldy to use. The types shown below in the code blocks include their name and type (all types are of type DataType) on the first line, and what the alias actually defines on the second line.",source:"@site/docs/gettingstarted/othertypes.md",sourceDirName:"gettingstarted",slug:"/gettingstarted/othertypes",permalink:"/PopGen.jl/docs/gettingstarted/othertypes",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/gettingstarted/othertypes.md",tags:[],version:"current",lastUpdatedAt:1635451805,formattedLastUpdatedAt:"10/28/2021",frontMatter:{id:"othertypes",title:"Other data types",sidebar_label:"Other data types"},sidebar:"docs",previous:{title:"The PopData type",permalink:"/PopGen.jl/docs/gettingstarted/popdata"},next:{title:"Provided datasets",permalink:"/PopGen.jl/docs/gettingstarted/datasets"}},c=[{value:"Genotype",id:"genotype",children:[{value:"SNP and Msat",id:"snp-and-msat",children:[],level:4}],level:3},{value:"GenoArray",id:"genoarray",children:[],level:3}],d={toc:c};function u(e){var t=e.components,n=(0,r.Z)(e,o);return(0,i.kt)("wrapper",(0,a.Z)({},d,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("p",null,"While not strictly their own composite types, we also define aliases for genotypes and vectors of genotypes, as their explicit types can get a little unwieldy to use. The types shown below in the code blocks include their name and type (all types are of type ",(0,i.kt)("inlineCode",{parentName:"p"},"DataType"),") on the first line, and what the alias actually defines on the second line."),(0,i.kt)("h3",{id:"genotype"},"Genotype"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"Genotype::DataType\nNTuple{N,<:Signed} where N\n")),(0,i.kt)("p",null,"An ",(0,i.kt)("inlineCode",{parentName:"p"},"NTuple")," is itself an alias for a ",(0,i.kt)("inlineCode",{parentName:"p"},"Tuple{Vararg{}}")," , but you can think of it as Tuple of ",(0,i.kt)("inlineCode",{parentName:"p"},"N")," length composed of items of a particular type, in this case it's items that are subtypes of ",(0,i.kt)("inlineCode",{parentName:"p"},"Signed")," (the integer types). The length of the tuple (",(0,i.kt)("inlineCode",{parentName:"p"},"N"),") will vary based on the ploidy of the sample, and the element ",(0,i.kt)("inlineCode",{parentName:"p"},"Type")," will vary whether the markers are snps (",(0,i.kt)("inlineCode",{parentName:"p"},"Int8"),") or microsatellites (",(0,i.kt)("inlineCode",{parentName:"p"},"Int16"),"), making this a pretty flexible (but immutable) structure."),(0,i.kt)("h4",{id:"snp-and-msat"},"SNP and Msat"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"Snp::NTuple{N,Int8} where N\nMsat::NTuple{N,Int16} where N\n")),(0,i.kt)("p",null,"These are convenience aliases for the two main kinds of NTuples of genotypes you will see.\nThese are typically used internally."),(0,i.kt)("h3",{id:"genoarray"},"GenoArray"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"GenoArray::DataType\nAbstractVector{S} where S<:Union{Missing,Genotype}\n")),(0,i.kt)("p",null,"As you can guess from the name, this Type wraps ",(0,i.kt)("inlineCode",{parentName:"p"},"Genotype")," into a Vector, while permitting ",(0,i.kt)("inlineCode",{parentName:"p"},"missing")," values (what's genetics without missing data!?). By using ",(0,i.kt)("inlineCode",{parentName:"p"},"AbstractVector")," (rather than ",(0,i.kt)("inlineCode",{parentName:"p"},"Vector"),"), we also have the flexibility of functions working on things like ",(0,i.kt)("inlineCode",{parentName:"p"},"SubArrays")," out of the box. "),(0,i.kt)("div",{className:"admonition admonition-note alert alert--secondary"},(0,i.kt)("div",{parentName:"div",className:"admonition-heading"},(0,i.kt)("h5",{parentName:"div"},(0,i.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,i.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,i.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.3 5.69a.942.942 0 0 1-.28-.7c0-.28.09-.52.28-.7.19-.18.42-.28.7-.28.28 0 .52.09.7.28.18.19.28.42.28.7 0 .28-.09.52-.28.7a1 1 0 0 1-.7.3c-.28 0-.52-.11-.7-.3zM8 7.99c-.02-.25-.11-.48-.31-.69-.2-.19-.42-.3-.69-.31H6c-.27.02-.48.13-.69.31-.2.2-.3.44-.31.69h1v3c.02.27.11.5.31.69.2.2.42.31.69.31h1c.27 0 .48-.11.69-.31.2-.19.3-.42.31-.69H8V7.98v.01zM7 2.3c-3.14 0-5.7 2.54-5.7 5.68 0 3.14 2.56 5.7 5.7 5.7s5.7-2.55 5.7-5.7c0-3.15-2.56-5.69-5.7-5.69v.01zM7 .98c3.86 0 7 3.14 7 7s-3.14 7-7 7-7-3.12-7-7 3.14-7 7-7z"}))),"why bother defining these aliases?")),(0,i.kt)("div",{parentName:"div",className:"admonition-content"},(0,i.kt)("p",{parentName:"div"},"Getting the most out of Julia and demonstrating good practices means making sure functions work on the things they're supposed to, and give informative error messages when the input isn't suitable for the function (a rare case of ",(0,i.kt)("em",{parentName:"p"},"wanting")," MethodErrors). Without these aliases, functions would either have vague definitions like ",(0,i.kt)("inlineCode",{parentName:"p"},"f(x,y,z) where x <: AbstractArray")," and potentially cause errors, or overly complicated definitions like ",(0,i.kt)("inlineCode",{parentName:"p"},"f(x::AbstractVector{S},y,z) where {N, T<:Signed,S<:NTuple{N,T}}")," and not be very legible. Instead, functions are written as ",(0,i.kt)("inlineCode",{parentName:"p"},"f(x,y,z) where x<:GenotypeArray"),", and that seems like a good compromise of getting the latter while looking like the former."))))}u.isMDXComponent=!0}}]);