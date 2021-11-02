"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[3944],{3905:function(e,n,t){t.d(n,{Zo:function(){return p},kt:function(){return d}});var a=t(7294);function i(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function r(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);n&&(a=a.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,a)}return t}function l(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?r(Object(t),!0).forEach((function(n){i(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):r(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function o(e,n){if(null==e)return{};var t,a,i=function(e,n){if(null==e)return{};var t,a,i={},r=Object.keys(e);for(a=0;a<r.length;a++)t=r[a],n.indexOf(t)>=0||(i[t]=e[t]);return i}(e,n);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);for(a=0;a<r.length;a++)t=r[a],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(i[t]=e[t])}return i}var s=a.createContext({}),u=function(e){var n=a.useContext(s),t=n;return e&&(t="function"==typeof e?e(n):l(l({},n),e)),t},p=function(e){var n=u(e.components);return a.createElement(s.Provider,{value:n},e.children)},m={inlineCode:"code",wrapper:function(e){var n=e.children;return a.createElement(a.Fragment,{},n)}},c=a.forwardRef((function(e,n){var t=e.components,i=e.mdxType,r=e.originalType,s=e.parentName,p=o(e,["components","mdxType","originalType","parentName"]),c=u(t),d=i,f=c["".concat(s,".").concat(d)]||c[d]||m[d]||r;return t?a.createElement(f,l(l({ref:n},p),{},{components:t})):a.createElement(f,l({ref:n},p))}));function d(e,n){var t=arguments,i=n&&n.mdxType;if("string"==typeof e||i){var r=t.length,l=new Array(r);l[0]=c;var o={};for(var s in n)hasOwnProperty.call(n,s)&&(o[s]=n[s]);o.originalType=e,o.mdxType="string"==typeof e?e:i,l[1]=o;for(var u=2;u<r;u++)l[u]=t[u];return a.createElement.apply(null,l)}return a.createElement.apply(null,t)}c.displayName="MDXCreateElement"},8215:function(e,n,t){var a=t(7294);n.Z=function(e){var n=e.children,t=e.hidden,i=e.className;return a.createElement("div",{role:"tabpanel",hidden:t,className:i},n)}},5064:function(e,n,t){t.d(n,{Z:function(){return c}});var a=t(7462),i=t(7294),r=t(2389),l=t(9443);var o=function(){var e=(0,i.useContext)(l.Z);if(null==e)throw new Error('"useUserPreferencesContext" is used outside of "Layout" component.');return e},s=t(3039),u=t(6010),p="tabItem_1uMI";function m(e){var n,t,a,r=e.lazy,l=e.block,m=e.defaultValue,c=e.values,d=e.groupId,f=e.className,v=i.Children.map(e.children,(function(e){if((0,i.isValidElement)(e)&&void 0!==e.props.value)return e;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})),h=null!=c?c:v.map((function(e){var n=e.props;return{value:n.value,label:n.label}})),b=(0,s.lx)(h,(function(e,n){return e.value===n.value}));if(b.length>0)throw new Error('Docusaurus error: Duplicate values "'+b.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.');var g=null===m?m:null!=(n=null!=m?m:null==(t=v.find((function(e){return e.props.default})))?void 0:t.props.value)?n:null==(a=v[0])?void 0:a.props.value;if(null!==g&&!h.some((function(e){return e.value===g})))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+g+'" but none of its children has the corresponding value. Available values are: '+h.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");var y=o(),k=y.tabGroupChoices,w=y.setTabGroupChoices,N=(0,i.useState)(g),_=N[0],j=N[1],P=[],T=(0,s.o5)().blockElementScrollPositionUntilNextRender;if(null!=d){var O=k[d];null!=O&&O!==_&&h.some((function(e){return e.value===O}))&&j(O)}var S=function(e){var n=e.currentTarget,t=P.indexOf(n),a=h[t].value;a!==_&&(T(n),j(a),null!=d&&w(d,a))},x=function(e){var n,t=null;switch(e.key){case"ArrowRight":var a=P.indexOf(e.currentTarget)+1;t=P[a]||P[0];break;case"ArrowLeft":var i=P.indexOf(e.currentTarget)-1;t=P[i]||P[P.length-1]}null==(n=t)||n.focus()};return i.createElement("div",{className:"tabs-container"},i.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,u.Z)("tabs",{"tabs--block":l},f)},h.map((function(e){var n=e.value,t=e.label;return i.createElement("li",{role:"tab",tabIndex:_===n?0:-1,"aria-selected":_===n,className:(0,u.Z)("tabs__item",p,{"tabs__item--active":_===n}),key:n,ref:function(e){return P.push(e)},onKeyDown:x,onFocus:S,onClick:S},null!=t?t:n)}))),r?(0,i.cloneElement)(v.filter((function(e){return e.props.value===_}))[0],{className:"margin-vert--md"}):i.createElement("div",{className:"margin-vert--md"},v.map((function(e,n){return(0,i.cloneElement)(e,{key:n,hidden:e.props.value!==_})}))))}function c(e){var n=(0,r.Z)();return i.createElement(m,(0,a.Z)({key:String(n)},e))}},9443:function(e,n,t){var a=(0,t(7294).createContext)(void 0);n.Z=a},7678:function(e,n,t){t.r(n),t.d(n,{frontMatter:function(){return u},contentTitle:function(){return p},metadata:function(){return m},toc:function(){return c},default:function(){return f}});var a=t(7462),i=t(3366),r=(t(7294),t(3905)),l=t(5064),o=t(8215),s=["components"],u={id:"simulate_samples",title:"Simulating Samples",sidebar_label:"Simulating Samples"},p=void 0,m={unversionedId:"simulations/simulate_samples",id:"simulations/simulate_samples",isDocsHomePage:!1,title:"Simulating Samples",description:"To perfom simulations, you will need add and import the package PopGenSims.jl (available here).",source:"@site/docs/simulations/simulations.md",sourceDirName:"simulations",slug:"/simulations/simulate_samples",permalink:"/PopGen.jl/docs/simulations/simulate_samples",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/simulations/simulations.md",tags:[],version:"current",lastUpdatedAt:1635529158,formattedLastUpdatedAt:"10/29/2021",frontMatter:{id:"simulate_samples",title:"Simulating Samples",sidebar_label:"Simulating Samples"},sidebar:"docs",previous:{title:"Pairwise F-Statistics",permalink:"/PopGen.jl/docs/analyses/fstatistics"},next:{title:"Breeding Crosses",permalink:"/PopGen.jl/docs/simulations/breedingcrosses"}},c=[{value:"Simulate samples within populations",id:"simulate-samples-within-populations",children:[],level:2}],d={toc:c};function f(e){var n=e.components,t=(0,i.Z)(e,s);return(0,r.kt)("wrapper",(0,a.Z)({},d,t,{components:n,mdxType:"MDXLayout"}),(0,r.kt)("div",{className:"admonition admonition-note alert alert--secondary"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.3 5.69a.942.942 0 0 1-.28-.7c0-.28.09-.52.28-.7.19-.18.42-.28.7-.28.28 0 .52.09.7.28.18.19.28.42.28.7 0 .28-.09.52-.28.7a1 1 0 0 1-.7.3c-.28 0-.52-.11-.7-.3zM8 7.99c-.02-.25-.11-.48-.31-.69-.2-.19-.42-.3-.69-.31H6c-.27.02-.48.13-.69.31-.2.2-.3.44-.31.69h1v3c.02.27.11.5.31.69.2.2.42.31.69.31h1c.27 0 .48-.11.69-.31.2-.19.3-.42.31-.69H8V7.98v.01zM7 2.3c-3.14 0-5.7 2.54-5.7 5.68 0 3.14 2.56 5.7 5.7 5.7s5.7-2.55 5.7-5.7c0-3.15-2.56-5.69-5.7-5.69v.01zM7 .98c3.86 0 7 3.14 7 7s-3.14 7-7 7-7-3.12-7-7 3.14-7 7-7z"}))),"Requires PopGenSims.jl")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"To perfom simulations, you will need add and import the package ",(0,r.kt)("inlineCode",{parentName:"p"},"PopGenSims.jl")," (available ",(0,r.kt)("a",{parentName:"p",href:"https://github.com/pdimens/PopGenSims.jl"},"here"),")."))),(0,r.kt)("h2",{id:"simulate-samples-within-populations"},"Simulate samples within populations"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"simulate(data::PopData; n::Int = 100)\n")),(0,r.kt)("p",null,"If you want to generate simulated data of a certain number of individuals per population, you can do so with the ",(0,r.kt)("inlineCode",{parentName:"p"},"simulate()")," function, which takes a ",(0,r.kt)("inlineCode",{parentName:"p"},"PopData")," object and simulates ",(0,r.kt)("inlineCode",{parentName:"p"},"n")," number of individuals per population using the allele frequencies of each population. This returns\nnew ",(0,r.kt)("inlineCode",{parentName:"p"},"PopData"),". The simulated samples will have the naming convention ",(0,r.kt)("inlineCode",{parentName:"p"},"sim_#")," where ",(0,r.kt)("inlineCode",{parentName:"p"},"#")," is a number from 1:",(0,r.kt)("inlineCode",{parentName:"p"},"n"),". These simulations return samples with the same ploidy as the source ",(0,r.kt)("inlineCode",{parentName:"p"},"PopData"),", but will ",(0,r.kt)("strong",{parentName:"p"},"not")," work on mixed-ploidy data. "),(0,r.kt)("p",null,"In the example below, we simulate 100 individuals per\npopulation using the nancycats data, which has 17 populations, therefore the resulting ",(0,r.kt)("inlineCode",{parentName:"p"},"PopData")," will have 1700 samples (100 samples x 17 populations)"),(0,r.kt)("p",null,(0,r.kt)("strong",{parentName:"p"},"Example")),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-julia"},"cats = @nancycats;\n\njulia> sims = simulate(cats , n = 100)\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 1700\n  Populations: 17\n")),(0,r.kt)("p",null,"Here is a look inside the ",(0,r.kt)("inlineCode",{parentName:"p"},"PopData")," to verify everything looks as expected."),(0,r.kt)(l.Z,{block:!0,defaultValue:"s",values:[{label:"sampleinfo",value:"s"},{label:"genodata",value:"g"}],mdxType:"Tabs"},(0,r.kt)(o.Z,{value:"s",mdxType:"TabItem"},(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"julia> sampleinfo(sims)\n\n  1700\xd75 DataFrame\n  Row \u2502 name      population  ploidy   \n      \u2502 String    String      Int8      \n\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n    1 \u2502 sim_1     1                2    \n    2 \u2502 sim_2     1                2    \n    3 \u2502 sim_3     1                2    \n    4 \u2502 sim_4     1                2    \n    5 \u2502 sim_5     1                2    \n  \u22ee   \u2502    \u22ee          \u22ee         \u22ee \n 1697 \u2502 sim_1697  17               2  \n 1698 \u2502 sim_1698  17               2  \n 1699 \u2502 sim_1699  17               2  \n 1700 \u2502 sim_1700  17               2  \n                                         1691 rows omitted \n"))),(0,r.kt)(o.Z,{value:"g",mdxType:"TabItem"},(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"julia> genodata(sims)\n15300\xd74 DataFrame\n   Row \u2502 name      population  locus   genotype   \n       \u2502 String    String      String  Tuple\u2026?    \n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n     1 \u2502 sim_1     1           fca8    (135, 143)\n     2 \u2502 sim_1     1           fca23   (136, 146)\n     3 \u2502 sim_1     1           fca43   (141, 145)\n     4 \u2502 sim_1     1           fca45   (120, 126)\n     5 \u2502 sim_1     1           fca77   (156, 156)\n   \u22ee   \u2502    \u22ee          \u22ee         \u22ee         \u22ee\n 15297 \u2502 sim_1700  17          fca78   (150, 150)\n 15298 \u2502 sim_1700  17          fca90   (197, 197)\n 15299 \u2502 sim_1700  17          fca96   (113, 113)\n 15300 \u2502 sim_1700  17          fca37   (208, 208)\n                                15291 rows omitted\n")))))}f.isMDXComponent=!0}}]);