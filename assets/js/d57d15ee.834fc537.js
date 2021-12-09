"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[3117],{3905:function(e,t,n){n.d(t,{Zo:function(){return m},kt:function(){return c}});var a=n(7294);function i(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function r(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function l(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?r(Object(n),!0).forEach((function(t){i(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):r(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function o(e,t){if(null==e)return{};var n,a,i=function(e,t){if(null==e)return{};var n,a,i={},r=Object.keys(e);for(a=0;a<r.length;a++)n=r[a],t.indexOf(n)>=0||(i[n]=e[n]);return i}(e,t);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);for(a=0;a<r.length;a++)n=r[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(i[n]=e[n])}return i}var s=a.createContext({}),d=function(e){var t=a.useContext(s),n=t;return e&&(n="function"==typeof e?e(t):l(l({},t),e)),n},m=function(e){var t=d(e.components);return a.createElement(s.Provider,{value:t},e.children)},p={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},u=a.forwardRef((function(e,t){var n=e.components,i=e.mdxType,r=e.originalType,s=e.parentName,m=o(e,["components","mdxType","originalType","parentName"]),u=d(n),c=i,h=u["".concat(s,".").concat(c)]||u[c]||p[c]||r;return n?a.createElement(h,l(l({ref:t},m),{},{components:n})):a.createElement(h,l({ref:t},m))}));function c(e,t){var n=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var r=n.length,l=new Array(r);l[0]=u;var o={};for(var s in t)hasOwnProperty.call(t,s)&&(o[s]=t[s]);o.originalType=e,o.mdxType="string"==typeof e?e:i,l[1]=o;for(var d=2;d<r;d++)l[d]=n[d];return a.createElement.apply(null,l)}return a.createElement.apply(null,n)}u.displayName="MDXCreateElement"},3919:function(e,t,n){function a(e){return!0===/^(\w*:|\/\/)/.test(e)}function i(e){return void 0!==e&&!a(e)}n.d(t,{b:function(){return a},Z:function(){return i}})},4996:function(e,t,n){n.d(t,{C:function(){return r},Z:function(){return l}});var a=n(2263),i=n(3919);function r(){var e=(0,a.Z)().siteConfig,t=(e=void 0===e?{}:e).baseUrl,n=void 0===t?"/":t,r=e.url;return{withBaseUrl:function(e,t){return function(e,t,n,a){var r=void 0===a?{}:a,l=r.forcePrependBaseUrl,o=void 0!==l&&l,s=r.absolute,d=void 0!==s&&s;if(!n)return n;if(n.startsWith("#"))return n;if((0,i.b)(n))return n;if(o)return t+n;var m=n.startsWith(t)?n:t+n.replace(/^\//,"");return d?e+m:m}(r,n,e,t)}}}function l(e,t){return void 0===t&&(t={}),(0,r().withBaseUrl)(e,t)}},8215:function(e,t,n){var a=n(7294);t.Z=function(e){var t=e.children,n=e.hidden,i=e.className;return a.createElement("div",{role:"tabpanel",hidden:n,className:i},t)}},6396:function(e,t,n){n.d(t,{Z:function(){return u}});var a=n(7462),i=n(7294),r=n(2389),l=n(9443);var o=function(){var e=(0,i.useContext)(l.Z);if(null==e)throw new Error('"useUserPreferencesContext" is used outside of "Layout" component.');return e},s=n(9521),d=n(6010),m="tabItem_1uMI";function p(e){var t,n,a,r=e.lazy,l=e.block,p=e.defaultValue,u=e.values,c=e.groupId,h=e.className,k=i.Children.map(e.children,(function(e){if((0,i.isValidElement)(e)&&void 0!==e.props.value)return e;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})),f=null!=u?u:k.map((function(e){var t=e.props;return{value:t.value,label:t.label}})),g=(0,s.lx)(f,(function(e,t){return e.value===t.value}));if(g.length>0)throw new Error('Docusaurus error: Duplicate values "'+g.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.');var v=null===p?p:null!=(t=null!=p?p:null==(n=k.find((function(e){return e.props.default})))?void 0:n.props.value)?t:null==(a=k[0])?void 0:a.props.value;if(null!==v&&!f.some((function(e){return e.value===v})))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+v+'" but none of its children has the corresponding value. Available values are: '+f.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");var N=o(),b=N.tabGroupChoices,y=N.setTabGroupChoices,w=(0,i.useState)(v),C=w[0],j=w[1],T=[],x=(0,s.o5)().blockElementScrollPositionUntilNextRender;if(null!=c){var M=b[c];null!=M&&M!==C&&f.some((function(e){return e.value===M}))&&j(M)}var _=function(e){var t=e.currentTarget,n=T.indexOf(t),a=f[n].value;a!==C&&(x(t),j(a),null!=c&&y(c,a))},P=function(e){var t,n=null;switch(e.key){case"ArrowRight":var a=T.indexOf(e.currentTarget)+1;n=T[a]||T[0];break;case"ArrowLeft":var i=T.indexOf(e.currentTarget)-1;n=T[i]||T[T.length-1]}null==(t=n)||t.focus()};return i.createElement("div",{className:"tabs-container"},i.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,d.Z)("tabs",{"tabs--block":l},h)},f.map((function(e){var t=e.value,n=e.label;return i.createElement("li",{role:"tab",tabIndex:C===t?0:-1,"aria-selected":C===t,className:(0,d.Z)("tabs__item",m,{"tabs__item--active":C===t}),key:t,ref:function(e){return T.push(e)},onKeyDown:P,onFocus:_,onClick:_},null!=n?n:t)}))),r?(0,i.cloneElement)(k.filter((function(e){return e.props.value===C}))[0],{className:"margin-vert--md"}):i.createElement("div",{className:"margin-vert--md"},k.map((function(e,t){return(0,i.cloneElement)(e,{key:t,hidden:e.props.value!==C})}))))}function u(e){var t=(0,r.Z)();return i.createElement(p,(0,a.Z)({key:String(t)},e))}},9443:function(e,t,n){var a=(0,n(7294).createContext)(void 0);t.Z=a},8619:function(e,t,n){n.r(t),n.d(t,{frontMatter:function(){return m},contentTitle:function(){return p},metadata:function(){return u},toc:function(){return c},default:function(){return k}});var a=n(7462),i=n(3366),r=(n(7294),n(3905)),l=n(6396),o=n(8215),s=n(4996),d=["components"],m={id:"relatedness",title:"Relatedness (Kinship)",sidebar_label:"Relatedness (Kinship)"},p=void 0,u={unversionedId:"analyses/relatedness",id:"analyses/relatedness",isDocsHomePage:!1,title:"Relatedness (Kinship)",description:"Background",source:"@site/docs/analyses/relatedness.md",sourceDirName:"analyses",slug:"/analyses/relatedness",permalink:"/PopGen.jl/docs/analyses/relatedness",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/analyses/relatedness.md",tags:[],version:"current",lastUpdatedAt:1627398780,formattedLastUpdatedAt:"7/27/2021",frontMatter:{id:"relatedness",title:"Relatedness (Kinship)",sidebar_label:"Relatedness (Kinship)"},sidebar:"docs",previous:{title:"Hardy-Weinberg Equilibrium",permalink:"/PopGen.jl/docs/analyses/hardyweinberg"},next:{title:"Pairwise F-Statistics",permalink:"/PopGen.jl/docs/analyses/fstatistics"}},c=[{value:"Background",id:"background",children:[],level:2},{value:"Calculate Relatedness",id:"calculate-relatedness",children:[{value:"Arguments",id:"arguments",children:[],level:4},{value:"Keyword Arguments",id:"keyword-arguments",children:[],level:4},{value:"Arguments",id:"arguments-1",children:[],level:4},{value:"Keyword Arguments",id:"keyword-arguments-1",children:[],level:4},{value:"Methods",id:"methods",children:[{value:"Examples",id:"examples",children:[],level:4}],level:3}],level:2},{value:"Posthoc analyses",id:"posthoc-analyses",children:[{value:"Arguments",id:"arguments-2",children:[],level:4},{value:"Keyword Arguments",id:"keyword-arguments-2",children:[],level:4}],level:2},{value:"Acknowledgements",id:"acknowledgements",children:[],level:2}],h={toc:c};function k(e){var t=e.components,m=(0,i.Z)(e,d);return(0,r.kt)("wrapper",(0,a.Z)({},h,m,{components:t,mdxType:"MDXLayout"}),(0,r.kt)("link",{rel:"stylesheet",href:(0,s.Z)("katex/katex.min.css")}),(0,r.kt)("h2",{id:"background"},"Background"),(0,r.kt)("p",null,"Sometimes you want or need to know the relatedness of individuals in a sample. Relatedness is exactly what its name implies: how related individuals are given some provided genetic information (e.g. full siblings, half siblings, etc.). Relatedness can be used in quantitative genetics to estimate heritability, additive genetic variances, and covariances. It can also be used in population genetics to study isolation-by-distance or population structure."),(0,r.kt)("p",null,'The goal of calculating relatedness from molecular markers is to accurately estimate the proportion of the genome which is identical by descent between two individuals. With a pedigree this is "relatively" straightforward. However, for large, natural, populations pedigrees tend not to exist and some brilliant minds have developed various equations to estimate the relatedness from a set of molecular markers. Given two diploid individuals, there are 9 "identity by descent" models available between them (',(0,r.kt)("a",{parentName:"p",href:"https://www.springer.com/gp/book/9783642884177"},"Jacquard 1975"),", paywall), as shown below (from ",(0,r.kt)("a",{parentName:"p",href:"https://www.genetics.org/content/163/3/1153.full"},"Milligan 2003"),"):"),(0,r.kt)("p",null,(0,r.kt)("img",{alt:"Jacquard IBD",src:n(4422).Z})),(0,r.kt)("p",null,"Broadly speaking there are two different ways of estimating genetic relatedness using molecular markers: methods of moments, and likelihood estimators. Generally, moments estimators will be faster but aren't constrained to being between the theoretical minimum and maximum values of 0 and 1. The likelihood estimators use likelihood functions derived from the set of Jacquard Identity States (above) to determine the most likely inheritance pattern. One difference between the two classes is (generally) moments estimators require an assumption of no inbreeding, while that assumption isn't necessarily required for likelihood estimators (though it does simplify the math). It is increasingly common to use multiple estimators on pairs, simulated from your molecular markers, with known relationships to determine the most appropriate estimator to use with your data."),(0,r.kt)("p",null,"PopGen.jl implements a wide variety of moments-based estimators: Blouin, Li & Horvitz, Loiselle, Lynch, Lynch/Li, Lynch & Ritland, Moran, Queller & Goodnight, Ritland, and Wang. Along with these, we provide an option to estimate mean, median, standard error, and confidence intervals using bootstrapping. Have a look at ",(0,r.kt)("a",{parentName:"p",href:"/blog/relatedness"},"the guide")," we provide detailing how to perform a basic relatedness analysis."),(0,r.kt)("div",{className:"admonition admonition-note alert alert--secondary"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.3 5.69a.942.942 0 0 1-.28-.7c0-.28.09-.52.28-.7.19-.18.42-.28.7-.28.28 0 .52.09.7.28.18.19.28.42.28.7 0 .28-.09.52-.28.7a1 1 0 0 1-.7.3c-.28 0-.52-.11-.7-.3zM8 7.99c-.02-.25-.11-.48-.31-.69-.2-.19-.42-.3-.69-.31H6c-.27.02-.48.13-.69.31-.2.2-.3.44-.31.69h1v3c.02.27.11.5.31.69.2.2.42.31.69.31h1c.27 0 .48-.11.69-.31.2-.19.3-.42.31-.69H8V7.98v.01zM7 2.3c-3.14 0-5.7 2.54-5.7 5.68 0 3.14 2.56 5.7 5.7 5.7s5.7-2.55 5.7-5.7c0-3.15-2.56-5.69-5.7-5.69v.01zM7 .98c3.86 0 7 3.14 7 7s-3.14 7-7 7-7-3.12-7-7 3.14-7 7-7z"}))),"A note about removing kin")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"There are reasons for removing kin in population genetics datasets. For one, there are no siblings/kin or mixed-generations in a Hardy-Weinberg Equilibrium population, and the inclusion of siblings/kin in analyses that rely on HWE assumptions ","[technically]"," violate those assumptions. However, there are also arguments to keep kin/siblings in your data, those data are important for effective population size, linkage disequilibrium, etc. "),(0,r.kt)("p",{parentName:"div"},"see  ",(0,r.kt)("a",{parentName:"p",href:"https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14022"},"Waples and Anderson (2017)")))),(0,r.kt)("h2",{id:"calculate-relatedness"},"Calculate Relatedness"),(0,r.kt)(l.Z,{block:!0,defaultValue:"a",values:[{label:"All vs. All",value:"a"},{label:"Specific Samples",value:"s"}],mdxType:"Tabs"},(0,r.kt)(o.Z,{value:"a",mdxType:"TabItem"},(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"relatedness(data::PopData; method::F, iterations::Int64, interval::Tuple(Float64, Float64))\n")),(0,r.kt)("p",null,"Return a dataframe of pairwise relatedness estimates for all individuals in a ",(0,r.kt)("inlineCode",{parentName:"p"},"PopData")," object using\nmethod ",(0,r.kt)("inlineCode",{parentName:"p"},"F")," (see below). To calculate means, median, standard error, and confidence intervals using bootstrapping,\nset ",(0,r.kt)("inlineCode",{parentName:"p"},"iterations = n")," where ",(0,r.kt)("inlineCode",{parentName:"p"},"n")," is an integer greater than ",(0,r.kt)("inlineCode",{parentName:"p"},"0")," (the default) corresponding to the number\nof bootstrap iterations you wish to perform for each pair. The default confidence interval is ",(0,r.kt)("inlineCode",{parentName:"p"},"(0.275, 0.975)")," (i.e. 95%), however that can be changed by supplying a ",(0,r.kt)("inlineCode",{parentName:"p"},"Tuple{Float64, Float64}")," of ",(0,r.kt)("inlineCode",{parentName:"p"},"(low, high)")," to the keyword ",(0,r.kt)("inlineCode",{parentName:"p"},"interval"),". ",(0,r.kt)("strong",{parentName:"p"},"Note:")," samples must be diploid."),(0,r.kt)("h4",{id:"arguments"},"Arguments"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"data")," : A PopData object")),(0,r.kt)("h4",{id:"keyword-arguments"},"Keyword Arguments"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"method")," : A method function or vector of method functions (see below)"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"iterations")," : The number of iterations to perform bootstrapping (default: ",(0,r.kt)("inlineCode",{parentName:"li"},"0"),", will not perform bootstrapping)"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"interval")," : A Tuple of (low, high) indicating the confidence intervals you would like for bootstrapping (default: ",(0,r.kt)("inlineCode",{parentName:"li"},"(0.275, 0.975)"),", i.e. 95%)"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"inbreeding")," : true/false of whether to consider inbreeding in the calculations (default: ",(0,r.kt)("inlineCode",{parentName:"li"},"false"),"). Only used in ",(0,r.kt)("inlineCode",{parentName:"li"},"dyadML"))),(0,r.kt)("p",null,(0,r.kt)("strong",{parentName:"p"},"Examples")),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"julia> cats = @nancycats;\n\njulia> relatedness(cats, method = Ritland)\n27966\xd74 DataFrame\n    Row \u2502 sample_1  sample_2  n_loci  Ritland\n        \u2502 String    String    Int64   Float64?\n \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n      1 \u2502 N215      N216           8   0.258824\n      2 \u2502 N215      N217           8   0.193238\n      3 \u2502 N215      N218           8   0.127497\n      4 \u2502 N215      N219           8   0.0453471\n   \u22ee   \u2502    \u22ee         \u22ee        \u22ee          \u22ee\n  27963 \u2502 N297      N290           7   0.189647\n  27964 \u2502 N281      N289           8   0.0892068\n  27965 \u2502 N281      N290           7   0.104614\n  27966 \u2502 N289      N290           7   0.0511663\n                                27958 rows omitted\n"))),(0,r.kt)(o.Z,{value:"s",mdxType:"TabItem"},(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"relatedness(data::PopData, sample_names::Vector{String}; method::F, iterations::Int64, interval::Tuple(Float64, Float64))\n")),(0,r.kt)("p",null,"Return a dataframe of pairwise relatedness estimates for all pairs of the supplied sample names in a ",(0,r.kt)("inlineCode",{parentName:"p"},"PopData")," object using\nmethod ",(0,r.kt)("inlineCode",{parentName:"p"},"F")," (see below). To calculate means, median, standard error, and confidence intervals using bootstrapping,\nset ",(0,r.kt)("inlineCode",{parentName:"p"},"iterations = n")," where ",(0,r.kt)("inlineCode",{parentName:"p"},"n")," is an integer greater than ",(0,r.kt)("inlineCode",{parentName:"p"},"0")," (the default) corresponding to the number\nof bootstrap iterations you wish to perform for each pair. The default confidence interval is ",(0,r.kt)("inlineCode",{parentName:"p"},"(0.275, 0.975)")," (i.e. 95%),\nhowever that can be changed by supplying a ",(0,r.kt)("inlineCode",{parentName:"p"},"Tuple{Float64, Float64}")," of ",(0,r.kt)("inlineCode",{parentName:"p"},"(low, high)")," to the keyword ",(0,r.kt)("inlineCode",{parentName:"p"},"interval"),".\n",(0,r.kt)("strong",{parentName:"p"},"Note:")," samples must be diploid."),(0,r.kt)("h4",{id:"arguments-1"},"Arguments"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"data")," : A PopData object"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"sample_names")," : A list of samples names to calculate relatedness for")),(0,r.kt)("h4",{id:"keyword-arguments-1"},"Keyword Arguments"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"method")," : A method function or vector of method functions (see below)"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"iterations")," : The number of iterations to perform bootstrapping (default: ",(0,r.kt)("inlineCode",{parentName:"li"},"0"),", will not perform bootstrapping)"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"interval")," : A Tuple of (low, high) indicating the confidence intervals you would like for bootstrapping (default: ",(0,r.kt)("inlineCode",{parentName:"li"},"(0.275, 0.975)"),", i.e. 95%)"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"inbreeding")," : true/false of whether to consider inbreeding in the calculations (default: ",(0,r.kt)("inlineCode",{parentName:"li"},"false"),"). Only used in ",(0,r.kt)("inlineCode",{parentName:"li"},"dyadML"))),(0,r.kt)("p",null,(0,r.kt)("strong",{parentName:"p"},"Examples")),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},'julia> cats = @nancycats;\n\njulia> relatedness(cats, ["N7", "N111", "N115"], method = [Ritland, Wang])\n3\xd75 DataFrame\n\u2502 Row \u2502 sample_1 \u2502 sample_2 \u2502 n_loci \u2502 Ritland    \u2502 Wang      \u2502\n\u2502     \u2502 String   \u2502 String   \u2502 Int64  \u2502 Float64?   \u2502 Float64?  \u2502\n\u251c\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524\n\u2502 1   \u2502 N7       \u2502 N111     \u2502 9      \u2502 -0.129432  \u2502 -0.395806 \u2502\n\u2502 2   \u2502 N7       \u2502 N115     \u2502 9      \u2502 -0.0183925 \u2502 0.0024775 \u2502\n\u2502 3   \u2502 N111     \u2502 N115     \u2502 9      \u2502 0.0240152  \u2502 0.183966  \u2502\n\n\njulia> relatedness(cats, ["N7", "N111", "N115"], method = [Loiselle, Moran], iterations = 100, interval = (0.025, 0.975))\n3\xd713 DataFrame. Omitted printing of 7 columns\n\u2502 Row \u2502 sample_1 \u2502 sample_2 \u2502 n_loci \u2502 Loiselle   \u2502 Loiselle_mean \u2502 Loiselle_median \u2502\n\u2502     \u2502 String   \u2502 String   \u2502 Int64  \u2502 Float64?   \u2502 Float64?      \u2502 Float64?        \u2502\n\u251c\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524\n\u2502 1   \u2502 N7       \u2502 N111     \u2502 9      \u2502 -0.101618  \u2502 0.0141364     \u2502 0.00703954      \u2502\n\u2502 2   \u2502 N7       \u2502 N115     \u2502 9      \u2502 -0.0428898 \u2502 0.0743497     \u2502 0.0798708       \u2502\n\u2502 3   \u2502 N111     \u2502 N115     \u2502 9      \u2502 0.13681    \u2502 0.266043      \u2502 0.257748        \u2502\n')))),(0,r.kt)("h3",{id:"methods"},"Methods"),(0,r.kt)("p",null,"There are several estimators available and are listed below. ",(0,r.kt)("inlineCode",{parentName:"p"},"relatedness")," takes the\nfunction names as arguments (",(0,r.kt)("strong",{parentName:"p"},"case sensitive"),"), therefore do not use quotes or colons\nin specifying the methods. Methods can be supplied as a vector. "),(0,r.kt)("table",null,(0,r.kt)("thead",{parentName:"table"},(0,r.kt)("tr",{parentName:"thead"},(0,r.kt)("th",{parentName:"tr",align:"left"},"Method"),(0,r.kt)("th",{parentName:"tr",align:"left"},"Type"),(0,r.kt)("th",{parentName:"tr",align:"left"},"Method Call"))),(0,r.kt)("tbody",{parentName:"table"},(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.1996.00094.x"},"Blouin 1996")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"Blouin"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1716461/"},"Li & Horvitz 1953")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"LiHorvitz"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://bsapubs.onlinelibrary.wiley.com/doi/abs/10.1002/j.1537-2197.1995.tb12679.x"},"Loiselle 1995")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"Loiselle"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://pubmed.ncbi.nlm.nih.gov/3193879/"},"Lynch 1988")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"Lynch"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://pubmed.ncbi.nlm.nih.gov/8514326/"},"Lynch/Li 1993")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"LynchLi"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://www.genetics.org/content/152/4/1753.short"},"Lynch & Ritland 1999")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"LynchRitland"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://www.jstor.org/stable/2332142?origin=crossref&seq=1#metadata_info_tab_contents"},"Moran 1950")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"Moran"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1989.tb04226.x"},"Queller & Goodnight 1989")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"QuellerGoodnight"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1996.tb02347.x"},"Ritland 1996")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("inlineCode",{parentName:"td"},"Ritland"))),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://pubmed.ncbi.nlm.nih.gov/12663552/"},"Milligan 2003"),' "DyadML"'),(0,r.kt)("td",{parentName:"tr",align:"left"},"maximum-likelihood"),(0,r.kt)("td",{parentName:"tr",align:"left"},"incomplete*")),(0,r.kt)("tr",{parentName:"tbody"},(0,r.kt)("td",{parentName:"tr",align:"left"},(0,r.kt)("a",{parentName:"td",href:"https://www.genetics.org/content/160/3/1203.short"},"Wang 2002")),(0,r.kt)("td",{parentName:"tr",align:"left"},"moments-based"),(0,r.kt)("td",{parentName:"tr",align:"left"},"incomplete*")))),(0,r.kt)("div",{className:"admonition admonition-note alert alert--secondary"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.3 5.69a.942.942 0 0 1-.28-.7c0-.28.09-.52.28-.7.19-.18.42-.28.7-.28.28 0 .52.09.7.28.18.19.28.42.28.7 0 .28-.09.52-.28.7a1 1 0 0 1-.7.3c-.28 0-.52-.11-.7-.3zM8 7.99c-.02-.25-.11-.48-.31-.69-.2-.19-.42-.3-.69-.31H6c-.27.02-.48.13-.69.31-.2.2-.3.44-.31.69h1v3c.02.27.11.5.31.69.2.2.42.31.69.31h1c.27 0 .48-.11.69-.31.2-.19.3-.42.31-.69H8V7.98v.01zM7 2.3c-3.14 0-5.7 2.54-5.7 5.68 0 3.14 2.56 5.7 5.7 5.7s5.7-2.55 5.7-5.7c0-3.15-2.56-5.69-5.7-5.69v.01zM7 .98c3.86 0 7 3.14 7 7s-3.14 7-7 7-7-3.12-7-7 3.14-7 7-7z"}))),"*more relatedness")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"Contact us or submit a pull request if you're interested in contributing to the relatedness methods. DyadML and Wang (2002) estimators are currently incomplete and the others\ncould use some optimizations. More help is always welcomed! Our wishlist also includes the KING method \ud83d\ude04"))),(0,r.kt)("h4",{id:"examples"},"Examples"),(0,r.kt)(l.Z,{block:!0,defaultValue:"s",values:[{label:"Single Method",value:"s"},{label:"Multiple Methods",value:"m"},{label:"With Bootstrapping",value:"b"}],mdxType:"Tabs"},(0,r.kt)(o.Z,{value:"s",mdxType:"TabItem"},(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"julia> cats = @nancycats;\n\njulia> cat_kin = relatendess(cats, samples(cats)[1:10], method = Ritland)\n"))),(0,r.kt)(o.Z,{value:"m",mdxType:"TabItem"},(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"julia> cats = @nancycats;\n\njulia> cat_kin = relatendess(cats, samples(cats)[1:10], method = [Moran, QuellerGoodnight])\n"))),(0,r.kt)(o.Z,{value:"b",mdxType:"TabItem"},(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"julia> cats = @nancycats;\n\njulia> cat_kin = relatendess(cats, method = [Ritland, Wang], iterations = 100)\n")))),(0,r.kt)("div",{className:"admonition admonition-tip alert alert--success"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"12",height:"16",viewBox:"0 0 12 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.5 0C3.48 0 1 2.19 1 5c0 .92.55 2.25 1 3 1.34 2.25 1.78 2.78 2 4v1h5v-1c.22-1.22.66-1.75 2-4 .45-.75 1-2.08 1-3 0-2.81-2.48-5-5.5-5zm3.64 7.48c-.25.44-.47.8-.67 1.11-.86 1.41-1.25 2.06-1.45 3.23-.02.05-.02.11-.02.17H5c0-.06 0-.13-.02-.17-.2-1.17-.59-1.83-1.45-3.23-.2-.31-.42-.67-.67-1.11C2.44 6.78 2 5.65 2 5c0-2.2 2.02-4 4.5-4 1.22 0 2.36.42 3.22 1.19C10.55 2.94 11 3.94 11 5c0 .66-.44 1.78-.86 2.48zM4 14h5c-.23 1.14-1.3 2-2.5 2s-2.27-.86-2.5-2z"}))),"autocompletion")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"Since the methods correspond to function names, they will tab-autocomplete when\ninputting them. For more information on a specific method, please see the respective docstring (e.g. ",(0,r.kt)("inlineCode",{parentName:"p"},"?Loiselle"),")."))),(0,r.kt)("h2",{id:"posthoc-analyses"},"Posthoc analyses"),(0,r.kt)("p",null,"There are several different kinds of things you can do with kinship information (e.g. network analysis), and one that's provided is lovingly\ncalled ",(0,r.kt)("inlineCode",{parentName:"p"},"relatedness_posthoc()"),", which performs a permutation analysis to\ntest if within-population relatedness is significantly greater than between-population relatedness for each population. We recommend that you\ncorrect for multiple testing using ",(0,r.kt)("inlineCode",{parentName:"p"},"MultipleTesting.jl"),"."),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"relatedness_posthoc(data::PopData, results::Union{DataFrame, NamedTuple}; iterations::Int)\n")),(0,r.kt)("h4",{id:"arguments-2"},"Arguments"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"data")," : A PopData object"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"results")," : the DataFrame or NamedTuple results from ",(0,r.kt)("inlineCode",{parentName:"li"},"relatedness()"))),(0,r.kt)("h4",{id:"keyword-arguments-2"},"Keyword Arguments"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"iterations")," : number of iterations for the permutation tests (default: ",(0,r.kt)("inlineCode",{parentName:"li"},"20000"),")")),(0,r.kt)("div",{className:"admonition admonition-tip alert alert--success"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"12",height:"16",viewBox:"0 0 12 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.5 0C3.48 0 1 2.19 1 5c0 .92.55 2.25 1 3 1.34 2.25 1.78 2.78 2 4v1h5v-1c.22-1.22.66-1.75 2-4 .45-.75 1-2.08 1-3 0-2.81-2.48-5-5.5-5zm3.64 7.48c-.25.44-.47.8-.67 1.11-.86 1.41-1.25 2.06-1.45 3.23-.02.05-.02.11-.02.17H5c0-.06 0-.13-.02-.17-.2-1.17-.59-1.83-1.45-3.23-.2-.31-.42-.67-.67-1.11C2.44 6.78 2 5.65 2 5c0-2.2 2.02-4 4.5-4 1.22 0 2.36.42 3.22 1.19C10.55 2.94 11 3.94 11 5c0 .66-.44 1.78-.86 2.48zM4 14h5c-.23 1.14-1.3 2-2.5 2s-2.27-.86-2.5-2z"}))),"not a great name")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},'We admit "relatedness_posthoc" is not a great name for this function. Please\ncontact us with better ideas! \ud83d\ude01'))),(0,r.kt)("p",null,(0,r.kt)("strong",{parentName:"p"},"Example")),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"julia> cats = @nancycats ;\n\njulia> rel_out = relatedness(cats, method = [Ritland, Moran], iterations = 100);\n\njulia> relatedness_posthoc(cats, rel_out)\n17x3 DataFrame\n Row \u2502 population  Ritland_P  Moran_P\n     \u2502 String      Float64    Float64\n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 1              5.0e-5   5.0e-5\n   2 \u2502 2              5.0e-5   5.0e-5\n   3 \u2502 3              5.0e-5   5.0e-5\n   4 \u2502 4              5.0e-5   5.0e-5\n   5 \u2502 5              5.0e-5   5.0e-5\n   6 \u2502 6              5.0e-5   5.0e-5\n   7 \u2502 7              5.0e-5   5.0e-5\n   8 \u2502 8              5.0e-5   5.0e-5\n   9 \u2502 9              5.0e-5   5.0e-5\n  10 \u2502 10             5.0e-5   5.0e-5\n  11 \u2502 11             5.0e-5   5.0e-5\n  12 \u2502 12             5.0e-5   5.0e-5\n  13 \u2502 13             5.0e-5   5.0e-5\n  14 \u2502 14             5.0e-5   5.0e-5\n  15 \u2502 15             5.0e-5   5.0e-5\n  16 \u2502 16             5.0e-5   5.0e-5\n  17 \u2502 17             5.0e-5   5.0e-5\n")),(0,r.kt)("hr",null),(0,r.kt)("h2",{id:"acknowledgements"},"Acknowledgements"),(0,r.kt)("p",null,"The relatedness methods were dutifully written and verified against R analogues by Jason Selwyn. These anaylses can take a while, especially if performing bootstrapping; we provide a progress bar via ",(0,r.kt)("inlineCode",{parentName:"p"},"ProgressMeter.jl")," so you can move on and focus on other things in the meantime."))}k.isMDXComponent=!0},4422:function(e,t,n){t.Z=n.p+"assets/images/jacquard_identitiies-dd9e2b1bc1f371819ea078c6ad3fc86c.jpg"}}]);