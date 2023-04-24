"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[3628],{3905:function(e,n,t){t.d(n,{Zo:function(){return p},kt:function(){return d}});var a=t(7294);function r(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function o(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);n&&(a=a.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,a)}return t}function i(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?o(Object(t),!0).forEach((function(n){r(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):o(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,a,r=function(e,n){if(null==e)return{};var t,a,r={},o=Object.keys(e);for(a=0;a<o.length;a++)t=o[a],n.indexOf(t)>=0||(r[t]=e[t]);return r}(e,n);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(a=0;a<o.length;a++)t=o[a],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(r[t]=e[t])}return r}var s=a.createContext({}),u=function(e){var n=a.useContext(s),t=n;return e&&(t="function"==typeof e?e(n):i(i({},n),e)),t},p=function(e){var n=u(e.components);return a.createElement(s.Provider,{value:n},e.children)},c="mdxType",f={inlineCode:"code",wrapper:function(e){var n=e.children;return a.createElement(a.Fragment,{},n)}},m=a.forwardRef((function(e,n){var t=e.components,r=e.mdxType,o=e.originalType,s=e.parentName,p=l(e,["components","mdxType","originalType","parentName"]),c=u(t),m=r,d=c["".concat(s,".").concat(m)]||c[m]||f[m]||o;return t?a.createElement(d,i(i({ref:n},p),{},{components:t})):a.createElement(d,i({ref:n},p))}));function d(e,n){var t=arguments,r=n&&n.mdxType;if("string"==typeof e||r){var o=t.length,i=new Array(o);i[0]=m;var l={};for(var s in n)hasOwnProperty.call(n,s)&&(l[s]=n[s]);l.originalType=e,l[c]="string"==typeof e?e:r,i[1]=l;for(var u=2;u<o;u++)i[u]=t[u];return a.createElement.apply(null,i)}return a.createElement.apply(null,t)}m.displayName="MDXCreateElement"},5162:function(e,n,t){t.d(n,{Z:function(){return i}});var a=t(7294),r=t(6010),o={tabItem:"tabItem_Ymn6"};function i(e){var n=e.children,t=e.hidden,i=e.className;return a.createElement("div",{role:"tabpanel",className:(0,r.Z)(o.tabItem,i),hidden:t},n)}},4866:function(e,n,t){t.d(n,{Z:function(){return y}});var a=t(7462),r=t(7294),o=t(6010),i=t(2466),l=t(6550),s=t(1980),u=t(7392),p=t(12);function c(e){return function(e){var n,t;return null!=(n=null==(t=r.Children.map(e,(function(e){if(!e||(0,r.isValidElement)(e)&&(n=e.props)&&"object"==typeof n&&"value"in n)return e;var n;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})))?void 0:t.filter(Boolean))?n:[]}(e).map((function(e){var n=e.props;return{value:n.value,label:n.label,attributes:n.attributes,default:n.default}}))}function f(e){var n=e.values,t=e.children;return(0,r.useMemo)((function(){var e=null!=n?n:c(t);return function(e){var n=(0,u.l)(e,(function(e,n){return e.value===n.value}));if(n.length>0)throw new Error('Docusaurus error: Duplicate values "'+n.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.')}(e),e}),[n,t])}function m(e){var n=e.value;return e.tabValues.some((function(e){return e.value===n}))}function d(e){var n=e.queryString,t=void 0!==n&&n,a=e.groupId,o=(0,l.k6)(),i=function(e){var n=e.queryString,t=void 0!==n&&n,a=e.groupId;if("string"==typeof t)return t;if(!1===t)return null;if(!0===t&&!a)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return null!=a?a:null}({queryString:t,groupId:a});return[(0,s._X)(i),(0,r.useCallback)((function(e){if(i){var n=new URLSearchParams(o.location.search);n.set(i,e),o.replace(Object.assign({},o.location,{search:n.toString()}))}}),[i,o])]}function g(e){var n,t,a,o,i=e.defaultValue,l=e.queryString,s=void 0!==l&&l,u=e.groupId,c=f(e),g=(0,r.useState)((function(){return function(e){var n,t=e.defaultValue,a=e.tabValues;if(0===a.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(t){if(!m({value:t,tabValues:a}))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+t+'" but none of its children has the corresponding value. Available values are: '+a.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");return t}var r=null!=(n=a.find((function(e){return e.default})))?n:a[0];if(!r)throw new Error("Unexpected error: 0 tabValues");return r.value}({defaultValue:i,tabValues:c})})),b=g[0],_=g[1],h=d({queryString:s,groupId:u}),k=h[0],v=h[1],y=(n=function(e){return e?"docusaurus.tab."+e:null}({groupId:u}.groupId),t=(0,p.Nk)(n),a=t[0],o=t[1],[a,(0,r.useCallback)((function(e){n&&o.set(e)}),[n,o])]),N=y[0],F=y[1],w=function(){var e=null!=k?k:N;return m({value:e,tabValues:c})?e:null}();return(0,r.useLayoutEffect)((function(){w&&_(w)}),[w]),{selectedValue:b,selectValue:(0,r.useCallback)((function(e){if(!m({value:e,tabValues:c}))throw new Error("Can't select invalid tab value="+e);_(e),v(e),F(e)}),[v,F,c]),tabValues:c}}var b=t(2389),_={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};function h(e){var n=e.className,t=e.block,l=e.selectedValue,s=e.selectValue,u=e.tabValues,p=[],c=(0,i.o5)().blockElementScrollPositionUntilNextRender,f=function(e){var n=e.currentTarget,t=p.indexOf(n),a=u[t].value;a!==l&&(c(n),s(a))},m=function(e){var n,t=null;switch(e.key){case"Enter":f(e);break;case"ArrowRight":var a,r=p.indexOf(e.currentTarget)+1;t=null!=(a=p[r])?a:p[0];break;case"ArrowLeft":var o,i=p.indexOf(e.currentTarget)-1;t=null!=(o=p[i])?o:p[p.length-1]}null==(n=t)||n.focus()};return r.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,o.Z)("tabs",{"tabs--block":t},n)},u.map((function(e){var n=e.value,t=e.label,i=e.attributes;return r.createElement("li",(0,a.Z)({role:"tab",tabIndex:l===n?0:-1,"aria-selected":l===n,key:n,ref:function(e){return p.push(e)},onKeyDown:m,onClick:f},i,{className:(0,o.Z)("tabs__item",_.tabItem,null==i?void 0:i.className,{"tabs__item--active":l===n})}),null!=t?t:n)})))}function k(e){var n=e.lazy,t=e.children,a=e.selectedValue,o=(Array.isArray(t)?t:[t]).filter(Boolean);if(n){var i=o.find((function(e){return e.props.value===a}));return i?(0,r.cloneElement)(i,{className:"margin-top--md"}):null}return r.createElement("div",{className:"margin-top--md"},o.map((function(e,n){return(0,r.cloneElement)(e,{key:n,hidden:e.props.value!==a})})))}function v(e){var n=g(e);return r.createElement("div",{className:(0,o.Z)("tabs-container",_.tabList)},r.createElement(h,(0,a.Z)({},e,n)),r.createElement(k,(0,a.Z)({},e,n)))}function y(e){var n=(0,b.Z)();return r.createElement(v,(0,a.Z)({key:String(n)},e))}},5352:function(e,n,t){t.r(n),t.d(n,{assets:function(){return f},contentTitle:function(){return p},default:function(){return b},frontMatter:function(){return u},metadata:function(){return c},toc:function(){return m}});var a=t(7462),r=t(3366),o=(t(7294),t(3905)),i=t(4866),l=t(5162),s=(t(4996),["components"]),u={id:"breedingcrosses",title:"Simulate Breeding Crosses",sidebar_label:"Breeding Crosses"},p=void 0,c={unversionedId:"simulations/breedingcrosses",id:"simulations/breedingcrosses",title:"Simulate Breeding Crosses",description:"To perfom simulations, you will need add and import the package PopGenSims.jl (available here).",source:"@site/docs/simulations/breedingcrosses.md",sourceDirName:"simulations",slug:"/simulations/breedingcrosses",permalink:"/PopGen.jl/docs/simulations/breedingcrosses",draft:!1,editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/simulations/breedingcrosses.md",tags:[],version:"current",lastUpdatedAt:1662485351,formattedLastUpdatedAt:"Sep 6, 2022",frontMatter:{id:"breedingcrosses",title:"Simulate Breeding Crosses",sidebar_label:"Breeding Crosses"},sidebar:"docs",previous:{title:"Simulating Samples",permalink:"/PopGen.jl/docs/simulations/"},next:{title:"Sibling Pairs",permalink:"/PopGen.jl/docs/simulations/sibship_simulations"}},f={},m=[{value:"Perform a cross",id:"perform-a-cross",level:2},{value:"Example",id:"example",level:4},{value:"Perform a cross/backcross",id:"perform-a-crossbackcross",level:2},{value:"Example",id:"example-1",level:4},{value:"Merge results",id:"merge-results",level:2},{value:"Example",id:"example-2",level:4}],d={toc:m},g="wrapper";function b(e){var n=e.components,t=(0,r.Z)(e,s);return(0,o.kt)(g,(0,a.Z)({},d,t,{components:n,mdxType:"MDXLayout"}),(0,o.kt)("admonition",{title:"Requires PopGenSims.jl",type:"note"},(0,o.kt)("p",{parentName:"admonition"},"To perfom simulations, you will need add and import the package ",(0,o.kt)("inlineCode",{parentName:"p"},"PopGenSims.jl")," (available ",(0,o.kt)("a",{parentName:"p",href:"https://github.com/pdimens/PopGenSims.jl"},"here"),").")),(0,o.kt)("p",null,"If you need to simulate offspring genotypes given mating between two individuals, the ",(0,o.kt)("inlineCode",{parentName:"p"},"cross()")," functions are available to simulate crosses and backcrosses."),(0,o.kt)("p",null,(0,o.kt)("strong",{parentName:"p"},"Currently, ",(0,o.kt)("inlineCode",{parentName:"strong"},"PopGenSims.jl")," can create crosses for:")),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},"haploids (ploidy = 1)"),(0,o.kt)("li",{parentName:"ul"},"diploids (ploidy = 2)"),(0,o.kt)("li",{parentName:"ul"},"tetraploids (ploidy = 4)"),(0,o.kt)("li",{parentName:"ul"},"hexaploids (ploidy = 6)"),(0,o.kt)("li",{parentName:"ul"},"octaploids (ploidy = 8)")),(0,o.kt)("h2",{id:"perform-a-cross"},"Perform a cross"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"cross(::PopData, parent1::String, parent2::String; n::Int, generation::String)\n")),(0,o.kt)("p",null,"The cross function performs a simple parental cross from individuals ",(0,o.kt)("inlineCode",{parentName:"p"},"parent1")," and ",(0,o.kt)("inlineCode",{parentName:"p"},"parent2")," in the same PopData object. The parents are strings of the names of the parents in the PopData. The keyword argument ",(0,o.kt)("inlineCode",{parentName:"p"},"n")," is the number of offspring to produce, and ",(0,o.kt)("inlineCode",{parentName:"p"},"generation")," is a keyword argument for the population name to the assign the offspring (default: ",(0,o.kt)("inlineCode",{parentName:"p"},'"F1"'),")."),(0,o.kt)("h4",{id:"example"},"Example"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'julia> cats = @nancycats;\n\njulia> f1 = cross(cats, "N111", "N107", n = 100000)\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 100000\n  Populations: 1\n')),(0,o.kt)("p",null,"Here is a look at the resulting ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")),(0,o.kt)(i.Z,{block:!0,defaultValue:"m",values:[{label:"metadata",value:"m"},{label:"genodata",value:"g"}],mdxType:"Tabs"},(0,o.kt)(l.Z,{value:"m",mdxType:"TabItem"},(0,o.kt)("p",null,"There are two things that should jump out at you:"),(0,o.kt)("ol",null,(0,o.kt)("li",{parentName:"ol"},"The ",(0,o.kt)("inlineCode",{parentName:"li"},"name")," of offspring are prepended with ",(0,o.kt)("inlineCode",{parentName:"li"},"generation")," and the ",(0,o.kt)("inlineCode",{parentName:"li"},"population")," is the ",(0,o.kt)("inlineCode",{parentName:"li"},"generation"),"."),(0,o.kt)("li",{parentName:"ol"},"There is a never-before-seen ",(0,o.kt)("inlineCode",{parentName:"li"},"parents")," column. This column exists for better record keeping of who has what parents if you are performing multiple crosses.")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre"},'julia> f1.sampleinfo\n100000\xd76 DataFrame\n\u2502 Row    \u2502 name                \u2502 ploidy \u2502 population \u2502 parents          \u2502\n\u2502        \u2502 String              \u2502 Int8   \u2502 String     \u2502 Tuple\u2026           \u2502\n\u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524\n\u2502 1      \u2502 F1_offspring_1      \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 2      \u2502 F1_offspring_2      \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 3      \u2502 F1_offspring_3      \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 4      \u2502 F1_offspring_4      \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 5      \u2502 F1_offspring_5      \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u22ee\n\u2502 99995  \u2502 F1_offspring_99995  \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 99996  \u2502 F1_offspring_99996  \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 99997  \u2502 F1_offspring_99997  \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 99998  \u2502 F1_offspring_99998  \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 99999  \u2502 F1_offspring_99999  \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n\u2502 100000 \u2502 F1_offspring_100000 \u2502 2      \u2502 F1         \u2502 ("N111", "N107") \u2502\n'))),(0,o.kt)(l.Z,{value:"g",mdxType:"TabItem"},(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre"},"julia> f1.genodata\n900000\xd74 DataFrame\n\u2502 Row    \u2502 name                \u2502 population \u2502 locus  \u2502 genotype   \u2502\n\u2502        \u2502 String              \u2502 String     \u2502 String \u2502 Tuple\u2026?    \u2502\n\u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524\n\u2502 1      \u2502 F1_offspring_1      \u2502 F1         \u2502 fca8   \u2502 (135, 135) \u2502\n\u2502 2      \u2502 F1_offspring_1      \u2502 F1         \u2502 fca23  \u2502 (146, 132) \u2502\n\u2502 3      \u2502 F1_offspring_1      \u2502 F1         \u2502 fca43  \u2502 (139, 145) \u2502\n\u2502 4      \u2502 F1_offspring_1      \u2502 F1         \u2502 fca45  \u2502 (132, 122) \u2502\n\u2502 5      \u2502 F1_offspring_1      \u2502 F1         \u2502 fca77  \u2502 (158, 150) \u2502\n\u22ee\n\u2502 899995 \u2502 F1_offspring_100000 \u2502 F1         \u2502 fca45  \u2502 (122, 128) \u2502\n\u2502 899996 \u2502 F1_offspring_100000 \u2502 F1         \u2502 fca77  \u2502 (158, 150) \u2502\n\u2502 899997 \u2502 F1_offspring_100000 \u2502 F1         \u2502 fca78  \u2502 (142, 150) \u2502\n\u2502 899998 \u2502 F1_offspring_100000 \u2502 F1         \u2502 fca90  \u2502 (201, 199) \u2502\n\u2502 899999 \u2502 F1_offspring_100000 \u2502 F1         \u2502 fca96  \u2502 (113, 103) \u2502\n\u2502 900000 \u2502 F1_offspring_100000 \u2502 F1         \u2502 fca37  \u2502 (214, 208) \u2502\n")))),(0,o.kt)("h2",{id:"perform-a-crossbackcross"},"Perform a cross/backcross"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'cross(PopData => "Parent1Name", PopData => "Parent2Name", n::Int, generation::String)\n')),(0,o.kt)("p",null,"This syntax uses the ",(0,o.kt)("inlineCode",{parentName:"p"},"Pair")," notation of ",(0,o.kt)("inlineCode",{parentName:"p"},'PopData => "Parent"')," to specify inputs. This method can be used for performing a cross like above, with the flexibility of parents allowed from two different ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," objects, which makes backcrosses possible. The keyword argument ",(0,o.kt)("inlineCode",{parentName:"p"},"n")," is the number of offspring to produce, and ",(0,o.kt)("inlineCode",{parentName:"p"},"generation")," is a keyword argument for the population name to the assign the offspring (default: ",(0,o.kt)("inlineCode",{parentName:"p"},'"F1"'),")."),(0,o.kt)("h4",{id:"example-1"},"Example"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'julia> f2_backcross = cross(cats => "N111", f1 => "F1_offspring_99", n = 100000, generation = "F2_manycats")\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 100000\n  Populations: 1\n')),(0,o.kt)("p",null,"And here you can see that ",(0,o.kt)("inlineCode",{parentName:"p"},"generation")," was again prepended to each offspring ",(0,o.kt)("inlineCode",{parentName:"p"},"name"),", along with assigned to the ",(0,o.kt)("inlineCode",{parentName:"p"},"population")," for each."),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre"},'julia> f2_backcross.sampleinfo\n100000\xd76 DataFrame\n\u2502 Row    \u2502 name                         \u2502 ploidy \u2502 population  \u2502 parents                     \u2502\n\u2502        \u2502 String                       \u2502 Int8   \u2502 String      \u2502 Tuple{String,String}        \u2502\n\u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524\n\u2502 1      \u2502 F2_manycats_offspring_1      \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 2      \u2502 F2_manycats_offspring_2      \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 3      \u2502 F2_manycats_offspring_3      \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 4      \u2502 F2_manycats_offspring_4      \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 5      \u2502 F2_manycats_offspring_5      \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u22ee\n\u2502 99995  \u2502 F2_manycats_offspring_99995  \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 99996  \u2502 F2_manycats_offspring_99996  \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 99997  \u2502 F2_manycats_offspring_99997  \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 99998  \u2502 F2_manycats_offspring_99998  \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 99999  \u2502 F2_manycats_offspring_99999  \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n\u2502 100000 \u2502 F2_manycats_offspring_100000 \u2502 2      \u2502 F2_manycats \u2502 ("N111", "F1_offspring_99") \u2502\n')),(0,o.kt)("admonition",{type:"caution"},(0,o.kt)("p",{parentName:"admonition"},"When crossing parents from different ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData"),", the parents must have the same loci. You will see error messages if they don't.")),(0,o.kt)("h2",{id:"merge-results"},"Merge results"),(0,o.kt)("p",null,"The ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," generated from breeding crosses can be combined used ",(0,o.kt)("inlineCode",{parentName:"p"},"append")," or ",(0,o.kt)("inlineCode",{parentName:"p"},"append!")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"append(::PopData, ::PopData)\nappend!(::PopData, ::PopData)\n")),(0,o.kt)("p",null,"These methods use outer joins and the ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," you are combining must have the same loci."),(0,o.kt)("h4",{id:"example-2"},"Example"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"# non mutating\ncrossed_sims = append(f1, f2_backcross)\n\n# mutating\nappend!(f1, f2_backcross)\n")))}b.isMDXComponent=!0}}]);