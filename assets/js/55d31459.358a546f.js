"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[829],{4137:function(e,t,n){n.d(t,{Zo:function(){return s},kt:function(){return f}});var a=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function i(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function l(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?i(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):i(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function o(e,t){if(null==e)return{};var n,a,r=function(e,t){if(null==e)return{};var n,a,r={},i=Object.keys(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var u=a.createContext({}),p=function(e){var t=a.useContext(u),n=t;return e&&(n="function"==typeof e?e(t):l(l({},t),e)),n},s=function(e){var t=p(e.components);return a.createElement(u.Provider,{value:t},e.children)},c="mdxType",m={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},d=a.forwardRef((function(e,t){var n=e.components,r=e.mdxType,i=e.originalType,u=e.parentName,s=o(e,["components","mdxType","originalType","parentName"]),c=p(n),d=r,f=c["".concat(u,".").concat(d)]||c[d]||m[d]||i;return n?a.createElement(f,l(l({ref:t},s),{},{components:n})):a.createElement(f,l({ref:t},s))}));function f(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var i=n.length,l=new Array(i);l[0]=d;var o={};for(var u in t)hasOwnProperty.call(t,u)&&(o[u]=t[u]);o.originalType=e,o[c]="string"==typeof e?e:r,l[1]=o;for(var p=2;p<i;p++)l[p]=n[p];return a.createElement.apply(null,l)}return a.createElement.apply(null,n)}d.displayName="MDXCreateElement"},425:function(e,t,n){n.d(t,{Z:function(){return l}});var a=n(7294),r=n(6010),i={tabItem:"tabItem_Ymn6"};function l(e){var t=e.children,n=e.hidden,l=e.className;return a.createElement("div",{role:"tabpanel",className:(0,r.Z)(i.tabItem,l),hidden:n},t)}},3992:function(e,t,n){n.d(t,{Z:function(){return N}});var a=n(7462),r=n(7294),i=n(6010),l=n(2957),o=n(6550),u=n(5238),p=n(3609),s=n(2560);function c(e){return function(e){var t,n;return null!=(t=null==(n=r.Children.map(e,(function(e){if(!e||(0,r.isValidElement)(e)&&(t=e.props)&&"object"==typeof t&&"value"in t)return e;var t;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})))?void 0:n.filter(Boolean))?t:[]}(e).map((function(e){var t=e.props;return{value:t.value,label:t.label,attributes:t.attributes,default:t.default}}))}function m(e){var t=e.values,n=e.children;return(0,r.useMemo)((function(){var e=null!=t?t:c(n);return function(e){var t=(0,p.l)(e,(function(e,t){return e.value===t.value}));if(t.length>0)throw new Error('Docusaurus error: Duplicate values "'+t.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.')}(e),e}),[t,n])}function d(e){var t=e.value;return e.tabValues.some((function(e){return e.value===t}))}function f(e){var t=e.queryString,n=void 0!==t&&t,a=e.groupId,i=(0,o.k6)(),l=function(e){var t=e.queryString,n=void 0!==t&&t,a=e.groupId;if("string"==typeof n)return n;if(!1===n)return null;if(!0===n&&!a)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return null!=a?a:null}({queryString:n,groupId:a});return[(0,u._X)(l),(0,r.useCallback)((function(e){if(l){var t=new URLSearchParams(i.location.search);t.set(l,e),i.replace(Object.assign({},i.location,{search:t.toString()}))}}),[l,i])]}function g(e){var t,n,a,i,l=e.defaultValue,o=e.queryString,u=void 0!==o&&o,p=e.groupId,c=m(e),g=(0,r.useState)((function(){return function(e){var t,n=e.defaultValue,a=e.tabValues;if(0===a.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(n){if(!d({value:n,tabValues:a}))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+n+'" but none of its children has the corresponding value. Available values are: '+a.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");return n}var r=null!=(t=a.find((function(e){return e.default})))?t:a[0];if(!r)throw new Error("Unexpected error: 0 tabValues");return r.value}({defaultValue:l,tabValues:c})})),k=g[0],v=g[1],h=f({queryString:u,groupId:p}),b=h[0],w=h[1],N=(t=function(e){return e?"docusaurus.tab."+e:null}({groupId:p}.groupId),n=(0,s.Nk)(t),a=n[0],i=n[1],[a,(0,r.useCallback)((function(e){t&&i.set(e)}),[t,i])]),y=N[0],C=N[1],O=function(){var e=null!=b?b:y;return d({value:e,tabValues:c})?e:null}();return(0,r.useLayoutEffect)((function(){O&&v(O)}),[O]),{selectedValue:k,selectValue:(0,r.useCallback)((function(e){if(!d({value:e,tabValues:c}))throw new Error("Can't select invalid tab value="+e);v(e),w(e),C(e)}),[w,C,c]),tabValues:c}}var k=n(1048),v={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};function h(e){var t=e.className,n=e.block,o=e.selectedValue,u=e.selectValue,p=e.tabValues,s=[],c=(0,l.o5)().blockElementScrollPositionUntilNextRender,m=function(e){var t=e.currentTarget,n=s.indexOf(t),a=p[n].value;a!==o&&(c(t),u(a))},d=function(e){var t,n=null;switch(e.key){case"Enter":m(e);break;case"ArrowRight":var a,r=s.indexOf(e.currentTarget)+1;n=null!=(a=s[r])?a:s[0];break;case"ArrowLeft":var i,l=s.indexOf(e.currentTarget)-1;n=null!=(i=s[l])?i:s[s.length-1]}null==(t=n)||t.focus()};return r.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,i.Z)("tabs",{"tabs--block":n},t)},p.map((function(e){var t=e.value,n=e.label,l=e.attributes;return r.createElement("li",(0,a.Z)({role:"tab",tabIndex:o===t?0:-1,"aria-selected":o===t,key:t,ref:function(e){return s.push(e)},onKeyDown:d,onClick:m},l,{className:(0,i.Z)("tabs__item",v.tabItem,null==l?void 0:l.className,{"tabs__item--active":o===t})}),null!=n?n:t)})))}function b(e){var t=e.lazy,n=e.children,a=e.selectedValue,i=(Array.isArray(n)?n:[n]).filter(Boolean);if(t){var l=i.find((function(e){return e.props.value===a}));return l?(0,r.cloneElement)(l,{className:"margin-top--md"}):null}return r.createElement("div",{className:"margin-top--md"},i.map((function(e,t){return(0,r.cloneElement)(e,{key:t,hidden:e.props.value!==a})})))}function w(e){var t=g(e);return r.createElement("div",{className:(0,i.Z)("tabs-container",v.tabList)},r.createElement(h,(0,a.Z)({},e,t)),r.createElement(b,(0,a.Z)({},e,t)))}function N(e){var t=(0,k.Z)();return r.createElement(w,(0,a.Z)({key:String(t)},e))}},1243:function(e,t,n){n.r(t),n.d(t,{assets:function(){return m},contentTitle:function(){return s},default:function(){return k},frontMatter:function(){return p},metadata:function(){return c},toc:function(){return d}});var a=n(7462),r=n(3366),i=(n(7294),n(4137)),l=n(3992),o=n(425),u=["components"],p={id:"genepop",title:"Genepop",sidebar_label:"Genepop"},s=void 0,c={unversionedId:"io/genepop",id:"io/genepop",title:"Genepop",description:"Import a genepop file as PopData",source:"@site/docs/io/genepop.md",sourceDirName:"io",slug:"/io/genepop",permalink:"/PopGen.jl/docs/io/genepop",draft:!1,editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/io/genepop.md",tags:[],version:"current",lastUpdatedAt:1652451030,formattedLastUpdatedAt:"May 13, 2022",frontMatter:{id:"genepop",title:"Genepop",sidebar_label:"Genepop"},sidebar:"docs",previous:{title:"Delimited",permalink:"/PopGen.jl/docs/io/delimited"},next:{title:"Plink",permalink:"/PopGen.jl/docs/io/plink"}},m={},d=[{value:"Import a genepop file as <code>PopData</code>",id:"import-a-genepop-file-as-popdata",level:2},{value:"Arguments",id:"arguments",level:3},{value:"Keyword Arguments",id:"keyword-arguments",level:3},{value:"Example",id:"example",level:3},{value:"Format",id:"format",level:3},{value:"Writing to a Genepop file",id:"writing-to-a-genepop-file",level:2},{value:"Arguments",id:"arguments-1",level:3},{value:"Keyword arguments",id:"keyword-arguments-1",level:3},{value:"Example",id:"example-1",level:3},{value:"Acknowledgements",id:"acknowledgements",level:2}],f={toc:d},g="wrapper";function k(e){var t=e.components,n=(0,r.Z)(e,u);return(0,i.kt)(g,(0,a.Z)({},f,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("h2",{id:"import-a-genepop-file-as-popdata"},"Import a genepop file as ",(0,i.kt)("inlineCode",{parentName:"h2"},"PopData")),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"genepop(infile; kwargs...)\n")),(0,i.kt)("h3",{id:"arguments"},"Arguments"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"infile::String")," : path to genepop file, in quotes")),(0,i.kt)("h3",{id:"keyword-arguments"},"Keyword Arguments"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"digits::Integer"),": number of digits denoting each allele (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"3"),")"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"popsep::String")," : word that separates populations in ",(0,i.kt)("inlineCode",{parentName:"li"},"infile"),' (default: "POP")'),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"diploid::Bool"),"  : whether samples are diploid for parsing optimizations (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"true"),")"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"silent::Bool"),"   : whether to print file information during import (default: ",(0,i.kt)("inlineCode",{parentName:"li"},"false"),")")),(0,i.kt)("admonition",{title:"population names",type:"info"},(0,i.kt)("p",{parentName:"admonition"},"By default, the file reader will assign numbers as population ID's (as Strings) in order of appearance in the genepop file. Use the ",(0,i.kt)("inlineCode",{parentName:"p"},"populations!")," function to rename these with your own population ID's.")),(0,i.kt)("h3",{id:"example"},"Example"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'julia> wasp_data = genepop("/data/wasp_hive.gen", digits = 3, popsep = "POP")\n')),(0,i.kt)("h3",{id:"format"},"Format"),(0,i.kt)("p",null,"Files must follow standard Genepop formatting:"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"First line is a comment (and skipped)"),(0,i.kt)("li",{parentName:"ul"},"Loci are listed after first line as one-per-line without commas or in single comma-separated row"),(0,i.kt)("li",{parentName:"ul"},"A line with a particular and consistent keyword must delimit populations",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("strong",{parentName:"li"},"must")," be the same word each time and not a unique population name"))),(0,i.kt)("li",{parentName:"ul"},"File is ",(0,i.kt)("strong",{parentName:"li"},"tab")," delimited or ",(0,i.kt)("strong",{parentName:"li"},"space")," delimited, but not both")),(0,i.kt)(l.Z,{block:!0,defaultValue:"v",values:[{label:"genepop w/loci stacked vertically",value:"v"},{label:"genepop w/loci stacked horizontally",value:"h"}],mdxType:"Tabs"},(0,i.kt)(o.Z,{value:"v",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},"Wasp populations in New York\nLocus1\nLocus2\nLocus3\nPOP\nOneida_01,  250230  564568  110100\nOneida_02,  252238  568558  100120\nOneida_03,  254230  564558  090100\nPOP\nNewcomb_01, 254230  564558  080100\nNewcomb_02, 000230  564558  090080\nNewcomb_03, 254230  000000  090100\nNewcomb_04, 254230  564000  090120\n"))),(0,i.kt)(o.Z,{value:"h",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},"Wasp populations in New York\nLocus1,Locus2,Locus3\nPOP\nOneida_01,  250230  564568  110100\nOneida_02,  252238  568558  100120\nOneida_03,  254230  564558  090100\nPOP\nNewcomb_01, 254230  564558  080100\nNewcomb_02, 000230  564558  090080\nNewcomb_03, 254230  000000  090100\nNewcomb_04, 254230  564000  090120\n")))),(0,i.kt)("h2",{id:"writing-to-a-genepop-file"},"Writing to a Genepop file"),(0,i.kt)("p",null,"All file writing options can be performed using ",(0,i.kt)("inlineCode",{parentName:"p"},"PopGen.write()"),", which calls ",(0,i.kt)("inlineCode",{parentName:"p"},"genpop")," when writing to a Genepop file."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'genepop(data::PopData; filename::String = "output.gen", digits::Int = 3, format::String = "vertical", miss::Int = 0)\n')),(0,i.kt)("p",null,"Writes a ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," object to a Genepop-formatted file."),(0,i.kt)("h3",{id:"arguments-1"},"Arguments"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"data"),": the ",(0,i.kt)("inlineCode",{parentName:"li"},"PopData")," object you wish to convert to a Genepop file")),(0,i.kt)("h3",{id:"keyword-arguments-1"},"Keyword arguments"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"filename::String"),": the output filename"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"digits::Integer"),": how many digits to format each allele",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},"e.g. ",(0,i.kt)("inlineCode",{parentName:"li"},"digits = 3")," will turn ",(0,i.kt)("inlineCode",{parentName:"li"},"(1, 2)")," into ",(0,i.kt)("inlineCode",{parentName:"li"},"001002")," "))),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"format::String")," : the way loci should be formatted ",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},"vertically (",(0,i.kt)("inlineCode",{parentName:"li"},'"v"')," or ",(0,i.kt)("inlineCode",{parentName:"li"},'"vertical"'),")"),(0,i.kt)("li",{parentName:"ul"},"hortizontally (",(0,i.kt)("inlineCode",{parentName:"li"},'"h"'),", or ",(0,i.kt)("inlineCode",{parentName:"li"},'"horizontal"'),")"),(0,i.kt)("li",{parentName:"ul"},"isolation-by-distance (",(0,i.kt)("inlineCode",{parentName:"li"},'"ibd"'),") where each sample is a population with coordinate data prepended"))),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"miss::Integer")," : how you would like missing values written ",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"0")," : as a genotype represented as a number of zeroes equal to ",(0,i.kt)("inlineCode",{parentName:"li"},"digits \xd7 ploidy")," like ",(0,i.kt)("inlineCode",{parentName:"li"},"000000")," (default) "),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"-9")," : as a single value ",(0,i.kt)("inlineCode",{parentName:"li"},"-9"))))),(0,i.kt)("h3",{id:"example-1"},"Example"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'cats = @nancycats;\nfewer_cats = omit(cats, name = samplenames(cats)[1:10]);\njulia> genepop(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")\n')),(0,i.kt)("hr",null),(0,i.kt)("h2",{id:"acknowledgements"},"Acknowledgements"),(0,i.kt)("p",null,"The original implementations of the importing parser were written using only Base Julia, and while the speed was fantastic, the memory footprint involved seemed unusually high (~650mb RAM to parse ",(0,i.kt)("inlineCode",{parentName:"p"},"gulfsharks"),", which is only 3.2mb in size). However, thanks to the efforts of ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/JuliaData/CSV.jl"},"CSV.jl"),", we leverage that package to preserve the speed and reduce the memory footprint quite a bit!"))}k.isMDXComponent=!0}}]);