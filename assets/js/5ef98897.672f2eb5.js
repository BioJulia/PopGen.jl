"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[724],{3905:function(e,t,a){a.d(t,{Zo:function(){return p},kt:function(){return m}});var n=a(7294);function r(e,t,a){return t in e?Object.defineProperty(e,t,{value:a,enumerable:!0,configurable:!0,writable:!0}):e[t]=a,e}function o(e,t){var a=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),a.push.apply(a,n)}return a}function i(e){for(var t=1;t<arguments.length;t++){var a=null!=arguments[t]?arguments[t]:{};t%2?o(Object(a),!0).forEach((function(t){r(e,t,a[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(a)):o(Object(a)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(a,t))}))}return e}function s(e,t){if(null==e)return{};var a,n,r=function(e,t){if(null==e)return{};var a,n,r={},o=Object.keys(e);for(n=0;n<o.length;n++)a=o[n],t.indexOf(a)>=0||(r[a]=e[a]);return r}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(n=0;n<o.length;n++)a=o[n],t.indexOf(a)>=0||Object.prototype.propertyIsEnumerable.call(e,a)&&(r[a]=e[a])}return r}var l=n.createContext({}),d=function(e){var t=n.useContext(l),a=t;return e&&(a="function"==typeof e?e(t):i(i({},t),e)),a},p=function(e){var t=d(e.components);return n.createElement(l.Provider,{value:t},e.children)},c={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},u=n.forwardRef((function(e,t){var a=e.components,r=e.mdxType,o=e.originalType,l=e.parentName,p=s(e,["components","mdxType","originalType","parentName"]),u=d(a),m=r,f=u["".concat(l,".").concat(m)]||u[m]||c[m]||o;return a?n.createElement(f,i(i({ref:t},p),{},{components:a})):n.createElement(f,i({ref:t},p))}));function m(e,t){var a=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var o=a.length,i=new Array(o);i[0]=u;var s={};for(var l in t)hasOwnProperty.call(t,l)&&(s[l]=t[l]);s.originalType=e,s.mdxType="string"==typeof e?e:r,i[1]=s;for(var d=2;d<o;d++)i[d]=a[d];return n.createElement.apply(null,i)}return n.createElement.apply(null,a)}u.displayName="MDXCreateElement"},3006:function(e,t,a){a.r(t),a.d(t,{frontMatter:function(){return s},contentTitle:function(){return l},metadata:function(){return d},toc:function(){return p},default:function(){return u}});var n=a(7462),r=a(3366),o=(a(7294),a(3905)),i=["components"],s={id:"datasets",title:"Provided datasets",sidebar_label:"Provided datasets"},l=void 0,d={unversionedId:"gettingstarted/datasets",id:"gettingstarted/datasets",isDocsHomePage:!1,title:"Provided datasets",description:"PopGen.jl (via PopGenCore.jl) provides two datasets as examples, nancycats and gulfsharks. The datasets can be retrieved using the dataset function, or their names as macros  (e.g. @gulfsharks).",source:"@site/docs/gettingstarted/datasets.md",sourceDirName:"gettingstarted",slug:"/gettingstarted/datasets",permalink:"/PopGen.jl/docs/gettingstarted/datasets",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/gettingstarted/datasets.md",tags:[],version:"current",lastUpdatedAt:1635454733,formattedLastUpdatedAt:"10/28/2021",frontMatter:{id:"datasets",title:"Provided datasets",sidebar_label:"Provided datasets"},sidebar:"docs",previous:{title:"Other data types",permalink:"/PopGen.jl/docs/gettingstarted/othertypes"},next:{title:"Reading data",permalink:"/PopGen.jl/docs/io/readingdata"}},p=[{value:"datasets",id:"datasets",children:[{value:"nancycats",id:"nancycats",children:[],level:3},{value:"gulfsharks",id:"gulfsharks",children:[],level:3}],level:2}],c={toc:p};function u(e){var t=e.components,a=(0,r.Z)(e,i);return(0,o.kt)("wrapper",(0,n.Z)({},c,a,{components:t,mdxType:"MDXLayout"}),(0,o.kt)("p",null,"PopGen.jl (via PopGenCore.jl) provides two datasets as examples, ",(0,o.kt)("inlineCode",{parentName:"p"},"nancycats")," and ",(0,o.kt)("inlineCode",{parentName:"p"},"gulfsharks"),". The datasets can be retrieved using the ",(0,o.kt)("inlineCode",{parentName:"p"},"dataset")," function, or their names as macros  (e.g. ",(0,o.kt)("inlineCode",{parentName:"p"},"@gulfsharks"),")."),(0,o.kt)("div",{className:"admonition admonition-info alert alert--info"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"}))),"identitcal methods")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},"The methods are identical (one is a wrapper for the other), but the benefit of calling the datasets directly by name is that you get the luxury of tab auto-completion \ud83d\ude01"))),(0,o.kt)("h2",{id:"datasets"},"datasets"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"PopGen.dataset(::String)\n")),(0,o.kt)("p",null,"Returns a ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," object of the dataset you would like to retrieve by calling the dataset as a string by name."),(0,o.kt)("p",null,(0,o.kt)("strong",{parentName:"p"},"Example:")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'sharks = PopGen.dataset("gulfsharks")\ncats = PopGen.dataset("nancycats")\n')),(0,o.kt)("h3",{id:"nancycats"},"nancycats"),(0,o.kt)("p",null,"We include the familiar nancycats microsatellite data, as featured in ",(0,o.kt)("a",{parentName:"p",href:"http://adegenet.r-forge.r-project.org"},"adegenet"),", for easy importing into PopGen.jl as ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData"),". As an alternative to ",(0,o.kt)("inlineCode",{parentName:"p"},"datasets"),", you can invoke the ",(0,o.kt)("inlineCode",{parentName:"p"},"@nancycats")," macro."),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"julia> cats = @nancycats\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 237\n  Populations: 17\n")),(0,o.kt)("p",null,"The spatial coordinates provided for the dataset in ",(0,o.kt)("inlineCode",{parentName:"p"},"adegenet")," are completely unfamiliar to us (and some geospatial folks we spoke to), so they have been omitted. If you recognize what coordinate system has 485.111 appear in Nancy, France, please let us know!"),(0,o.kt)("h3",{id:"gulfsharks"},"gulfsharks"),(0,o.kt)("p",null,"We also include the SNP dataset used in ",(0,o.kt)("a",{parentName:"p",href:"https://link.springer.com/article/10.1007/s00227-019-3533-1"},"Dimens ",(0,o.kt)("em",{parentName:"a"},"et al.")," 2019")," since it was already on hand. Like ",(0,o.kt)("inlineCode",{parentName:"p"},"nancycats"),", we provide a convenient function to load these data into PopGen.jl as ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData"),". As an alternative to ",(0,o.kt)("inlineCode",{parentName:"p"},"dataset"),", you can invoke the ",(0,o.kt)("inlineCode",{parentName:"p"},"@gulfsharks")," macro. "),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'julia> sharks = @gulfsharks\nPopData{Diploid, 2209 SNP loci}\n  Samples: 212\n  Populations: 7\n  Other Info: ["longitude", "latitude"]\n')))}u.isMDXComponent=!0}}]);