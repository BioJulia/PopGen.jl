"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[1409],{3905:function(e,n,a){a.d(n,{Zo:function(){return p},kt:function(){return m}});var t=a(7294);function o(e,n,a){return n in e?Object.defineProperty(e,n,{value:a,enumerable:!0,configurable:!0,writable:!0}):e[n]=a,e}function l(e,n){var a=Object.keys(e);if(Object.getOwnPropertySymbols){var t=Object.getOwnPropertySymbols(e);n&&(t=t.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),a.push.apply(a,t)}return a}function i(e){for(var n=1;n<arguments.length;n++){var a=null!=arguments[n]?arguments[n]:{};n%2?l(Object(a),!0).forEach((function(n){o(e,n,a[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(a)):l(Object(a)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(a,n))}))}return e}function r(e,n){if(null==e)return{};var a,t,o=function(e,n){if(null==e)return{};var a,t,o={},l=Object.keys(e);for(t=0;t<l.length;t++)a=l[t],n.indexOf(a)>=0||(o[a]=e[a]);return o}(e,n);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);for(t=0;t<l.length;t++)a=l[t],n.indexOf(a)>=0||Object.prototype.propertyIsEnumerable.call(e,a)&&(o[a]=e[a])}return o}var s=t.createContext({}),c=function(e){var n=t.useContext(s),a=n;return e&&(a="function"==typeof e?e(n):i(i({},n),e)),a},p=function(e){var n=c(e.components);return t.createElement(s.Provider,{value:n},e.children)},u={inlineCode:"code",wrapper:function(e){var n=e.children;return t.createElement(t.Fragment,{},n)}},d=t.forwardRef((function(e,n){var a=e.components,o=e.mdxType,l=e.originalType,s=e.parentName,p=r(e,["components","mdxType","originalType","parentName"]),d=c(a),m=o,g=d["".concat(s,".").concat(m)]||d[m]||u[m]||l;return a?t.createElement(g,i(i({ref:n},p),{},{components:a})):t.createElement(g,i({ref:n},p))}));function m(e,n){var a=arguments,o=n&&n.mdxType;if("string"==typeof e||o){var l=a.length,i=new Array(l);i[0]=d;var r={};for(var s in n)hasOwnProperty.call(n,s)&&(r[s]=n[s]);r.originalType=e,r.mdxType="string"==typeof e?e:o,i[1]=r;for(var c=2;c<l;c++)i[c]=a[c];return t.createElement.apply(null,i)}return t.createElement.apply(null,a)}d.displayName="MDXCreateElement"},8215:function(e,n,a){var t=a(7294);n.Z=function(e){var n=e.children,a=e.hidden,o=e.className;return t.createElement("div",{role:"tabpanel",hidden:a,className:o},n)}},6396:function(e,n,a){a.d(n,{Z:function(){return d}});var t=a(7462),o=a(7294),l=a(2389),i=a(9443);var r=function(){var e=(0,o.useContext)(i.Z);if(null==e)throw new Error('"useUserPreferencesContext" is used outside of "Layout" component.');return e},s=a(9521),c=a(6010),p="tabItem_1uMI";function u(e){var n,a,t,l=e.lazy,i=e.block,u=e.defaultValue,d=e.values,m=e.groupId,g=e.className,h=o.Children.map(e.children,(function(e){if((0,o.isValidElement)(e)&&void 0!==e.props.value)return e;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})),f=null!=d?d:h.map((function(e){var n=e.props;return{value:n.value,label:n.label}})),v=(0,s.lx)(f,(function(e,n){return e.value===n.value}));if(v.length>0)throw new Error('Docusaurus error: Duplicate values "'+v.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.');var k=null===u?u:null!=(n=null!=u?u:null==(a=h.find((function(e){return e.props.default})))?void 0:a.props.value)?n:null==(t=h[0])?void 0:t.props.value;if(null!==k&&!f.some((function(e){return e.value===k})))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+k+'" but none of its children has the corresponding value. Available values are: '+f.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");var y=r(),N=y.tabGroupChoices,w=y.setTabGroupChoices,b=(0,o.useState)(k),C=b[0],_=b[1],j=[],P=(0,s.o5)().blockElementScrollPositionUntilNextRender;if(null!=m){var T=N[m];null!=T&&T!==C&&f.some((function(e){return e.value===T}))&&_(T)}var D=function(e){var n=e.currentTarget,a=j.indexOf(n),t=f[a].value;t!==C&&(P(n),_(t),null!=m&&w(m,t))},S=function(e){var n,a=null;switch(e.key){case"ArrowRight":var t=j.indexOf(e.currentTarget)+1;a=j[t]||j[0];break;case"ArrowLeft":var o=j.indexOf(e.currentTarget)-1;a=j[o]||j[j.length-1]}null==(n=a)||n.focus()};return o.createElement("div",{className:"tabs-container"},o.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,c.Z)("tabs",{"tabs--block":i},g)},f.map((function(e){var n=e.value,a=e.label;return o.createElement("li",{role:"tab",tabIndex:C===n?0:-1,"aria-selected":C===n,className:(0,c.Z)("tabs__item",p,{"tabs__item--active":C===n}),key:n,ref:function(e){return j.push(e)},onKeyDown:S,onFocus:D,onClick:D},null!=a?a:n)}))),l?(0,o.cloneElement)(h.filter((function(e){return e.props.value===C}))[0],{className:"margin-vert--md"}):o.createElement("div",{className:"margin-vert--md"},h.map((function(e,n){return(0,o.cloneElement)(e,{key:n,hidden:e.props.value!==C})}))))}function d(e){var n=(0,l.Z)();return o.createElement(u,(0,t.Z)({key:String(n)},e))}},9443:function(e,n,a){var t=(0,a(7294).createContext)(void 0);n.Z=t},2699:function(e,n,a){a.r(n),a.d(n,{frontMatter:function(){return c},contentTitle:function(){return p},metadata:function(){return u},toc:function(){return d},default:function(){return g}});var t=a(7462),o=a(3366),l=(a(7294),a(3905)),i=a(6396),r=a(8215),s=["components"],c={id:"viewdata",title:"Viewing data",sidebar_label:"Viewing data"},p=void 0,u={unversionedId:"workingwithpopdata/viewdata",id:"workingwithpopdata/viewdata",isDocsHomePage:!1,title:"Viewing data",description:"A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code \"tabs\".",source:"@site/docs/workingwithpopdata/accessingelements.md",sourceDirName:"workingwithpopdata",slug:"/workingwithpopdata/viewdata",permalink:"/PopGen.jl/docs/workingwithpopdata/viewdata",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/workingwithpopdata/accessingelements.md",tags:[],version:"current",lastUpdatedAt:1639063144,formattedLastUpdatedAt:"12/9/2021",frontMatter:{id:"viewdata",title:"Viewing data",sidebar_label:"Viewing data"},sidebar:"docs",previous:{title:"Working with PopData",permalink:"/PopGen.jl/docs/workingwithpopdata/workingwithpopdata"},next:{title:"Adding sample/locus data",permalink:"/PopGen.jl/docs/workingwithpopdata/addingdata"}},d=[{value:"Loading in the data",id:"loading-in-the-data",children:[],level:2},{value:"The metadata (data about the data)",id:"the-metadata-data-about-the-data",children:[{value:"sampleinfo",id:"sampleinfo",children:[],level:3},{value:"locusinfo",id:"locusinfo",children:[],level:3}],level:2},{value:"The genotype table",id:"the-genotype-table",children:[{value:"genodata",id:"genodata",children:[],level:3}],level:2},{value:"View specific information",id:"view-specific-information",children:[{value:"sample names",id:"sample-names",children:[],level:3},{value:"locus names",id:"locus-names",children:[],level:3}],level:2},{value:"View genotypes",id:"view-genotypes",children:[{value:"all genotypes in one locus or sample",id:"all-genotypes-in-one-locus-or-sample",children:[],level:3},{value:"one sample, one locus",id:"one-sample-one-locus",children:[],level:3},{value:"many samples, one locus",id:"many-samples-one-locus",children:[],level:3},{value:"one sample, many loci",id:"one-sample-many-loci",children:[],level:3},{value:"many samples, many loci",id:"many-samples-many-loci",children:[],level:3}],level:2}],m={toc:d};function g(e){var n=e.components,a=(0,o.Z)(e,s);return(0,l.kt)("wrapper",(0,t.Z)({},m,a,{components:n,mdxType:"MDXLayout"}),(0,l.kt)("p",null,"A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code \"tabs\"."),(0,l.kt)("div",{className:"admonition admonition-danger alert alert--danger"},(0,l.kt)("div",{parentName:"div",className:"admonition-heading"},(0,l.kt)("h5",{parentName:"div"},(0,l.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,l.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"12",height:"16",viewBox:"0 0 12 16"},(0,l.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M5.05.31c.81 2.17.41 3.38-.52 4.31C3.55 5.67 1.98 6.45.9 7.98c-1.45 2.05-1.7 6.53 3.53 7.7-2.2-1.16-2.67-4.52-.3-6.61-.61 2.03.53 3.33 1.94 2.86 1.39-.47 2.3.53 2.27 1.67-.02.78-.31 1.44-1.13 1.81 3.42-.59 4.78-3.42 4.78-5.56 0-2.84-2.53-3.22-1.25-5.61-1.52.13-2.03 1.13-1.89 2.75.09 1.08-1.02 1.8-1.86 1.33-.67-.41-.66-1.19-.06-1.78C8.18 5.31 8.68 2.45 5.05.32L5.03.3l.02.01z"}))),"don't manually edit or sort")),(0,l.kt)("div",{parentName:"div",className:"admonition-content"},(0,l.kt)("p",{parentName:"div"},"There are specific relationships between the record entries in ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," objects, so ",(0,l.kt)("strong",{parentName:"p"},"do not use")," ",(0,l.kt)("inlineCode",{parentName:"p"},"sort"),", ",(0,l.kt)("inlineCode",{parentName:"p"},"sort!"),", or manually arrange/add/delete anything in PopData. There are included functions to remove samples or loci, rename things, add location data, etc. "))),(0,l.kt)("h2",{id:"loading-in-the-data"},"Loading in the data"),(0,l.kt)("p",null,"Let's keep things simple by loading in the nancycats data and calling it ",(0,l.kt)("inlineCode",{parentName:"p"},"ncats"),"."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> ncats = @nancycats\nPopData{Diploid, 9 Microsatellite loci}\n  Samples: 237\n  Populations: 17\n")),(0,l.kt)("p",null,"Now that we have nancycats loaded in, we can use standard Julia accessor conventions to view the elements within our PopData. The DataFrames uses the convention ",(0,l.kt)("inlineCode",{parentName:"p"},"dataframe.colname")," to directly access the columns we want."),(0,l.kt)("h2",{id:"the-metadata-data-about-the-data"},"The metadata (data about the data)"),(0,l.kt)("p",null,"Some critical information about the data is front-loaded into a PopData object to eliminate constantly getting these values in calculations.\nTo view this information, use ",(0,l.kt)("inlineCode",{parentName:"p"},"metadata()"),"."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"julia> metadata(ncats)\n ploidy:      2\n loci:        9\n samples:     237\n populations: 17\n biallelic:   false\n")),(0,l.kt)("p",null,"Included in ",(0,l.kt)("inlineCode",{parentName:"p"},"metadata")," are two DataFrames, one for sample information, and another for locus information."),(0,l.kt)(i.Z,{block:!0,defaultValue:"s",values:[{label:"sample information",value:"s"},{label:"locus information",value:"l"}],mdxType:"Tabs"},(0,l.kt)(r.Z,{value:"s",mdxType:"TabItem"},(0,l.kt)("h3",{id:"sampleinfo"},"sampleinfo"),(0,l.kt)("p",null,"To view the sample information, you can use ",(0,l.kt)("inlineCode",{parentName:"p"},"sampleinfo()")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> sampleinfo(ncats)\n237\xd73 DataFrame\n Row \u2502 name      population  ploidy \n     \u2502 String7\u2026  String      Int8   \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 N217      1                2\n   2 \u2502 N218      1                2\n   3 \u2502 N219      1                2\n   4 \u2502 N220      1                2\n   5 \u2502 N221      1                2\n   6 \u2502 N222      1                2\n  \u22ee  \u2502    \u22ee          \u22ee         \u22ee\n 232 \u2502 N197      14               2\n 233 \u2502 N198      14               2\n 234 \u2502 N199      14               2\n 235 \u2502 N200      14               2\n 236 \u2502 N201      14               2\n 237 \u2502 N206      14               2\n                    222 rows omitted\n\n")),(0,l.kt)("p",null,"Using the standard DataFrames ",(0,l.kt)("inlineCode",{parentName:"p"},"getindex")," methods, we can access these columns like so:"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> sinfo = sampleinfo(ncats) ;\njulia> sinfo.name\n237-element Array{String,1}:\n "N1"  \n "N2"  \n "N3"  \n "N4"  \n "N5"  \n "N6"  \n "N7"  \n "N8"  \n \u22ee     \n "N230"\n "N231"\n "N232"\n "N233"\n "N234"\n "N235"\n "N236"\n "N237"\n'))),(0,l.kt)(r.Z,{value:"l",mdxType:"TabItem"},(0,l.kt)("h3",{id:"locusinfo"},"locusinfo"),(0,l.kt)("p",null,"To view the locus information, you can use ",(0,l.kt)("inlineCode",{parentName:"p"},"locusinfo()"),". Locus information is not mandatory,\nbut present if needed for future analyses."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> locusinfo(ncats)\n9\xd74 DataFrame\n Row \u2502 chromosome  locus   cm       bp   \n     \u2502 Int8        String  Float64  Int64 \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502          0  fca8       0       0\n   2 \u2502          0  fca23      0       0\n   3 \u2502          0  fca43      0       0\n   4 \u2502          0  fca45      0       0\n   5 \u2502          0  fca77      0       0\n   6 \u2502          0  fca78      0       0\n   7 \u2502          0  fca90      0       0\n   8 \u2502          0  fca96      0       0\n   9 \u2502          0  fca37      0       0\n")))),(0,l.kt)("hr",null),(0,l.kt)("h2",{id:"the-genotype-table"},"The genotype table"),(0,l.kt)("h3",{id:"genodata"},"genodata"),(0,l.kt)("p",null,"You can view the genotype information with ",(0,l.kt)("inlineCode",{parentName:"p"},"genodata()"),"."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> genodata(ncats)\n2133\xd74 DataFrame\n  Row \u2502 name    population  locus   genotype\n      \u2502 String  String      String  Tuple\u2026?\n\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n    1 \u2502 N215    1           fca8    missing\n    2 \u2502 N216    1           fca8    missing\n    3 \u2502 N217    1           fca8    (135, 143)\n    4 \u2502 N218    1           fca8    (133, 135)\n    5 \u2502 N219    1           fca8    (133, 135)\n    6 \u2502 N220    1           fca8    (135, 143)\n  \u22ee   \u2502   \u22ee         \u22ee         \u22ee         \u22ee\n 2128 \u2502 N295    17          fca37   (208, 208)\n 2129 \u2502 N296    17          fca37   (208, 220)\n 2130 \u2502 N297    17          fca37   (208, 208)\n 2131 \u2502 N281    17          fca37   (208, 208)\n 2132 \u2502 N289    17          fca37   (208, 208)\n 2133 \u2502 N290    17          fca37   (208, 208)\n                              2121 rows omitted\n")),(0,l.kt)("p",null,'Because the genotype data is in long format (aka "tidy"), accessing genotypes in a meaningful way is fairly\nstraightforward if you have any experience with dataframe manipulation. For a deeper look into indexing ',(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),",\nread ",(0,l.kt)("a",{parentName:"p",href:"advancedindexing"},"Advanced PopData Indexing")," "),(0,l.kt)("p",null,"The functions here help you inspect your ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")," and pull information from it easily."),(0,l.kt)("h2",{id:"view-specific-information"},"View specific information"),(0,l.kt)("h3",{id:"sample-names"},"sample names"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"samplenames(data::PopData)\n")),(0,l.kt)("p",null,"View individual/sample names in a ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),". "),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> samples(sharks)\n212-element Array{String,1}:\n "cc_001" \n "cc_002" \n "cc_003" \n "cc_005" \n "cc_007" \n \u22ee        \n "seg_027"\n "seg_028"\n "seg_029"\n "seg_030"\n "seg_031"\n')),(0,l.kt)("h3",{id:"locus-names"},"locus names"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"loci(data::PopData)\n")),(0,l.kt)("p",null,"Returns a vector of strings of the loci names in a ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> loci(sharks)\n2213-element Array{String,1}:\n "contig_35208"\n "contig_23109"\n "contig_4493" \n "contig_10742"\n "contig_14898"\n \u22ee             \n "contig_43517"\n "contig_27356"\n "contig_475"  \n "contig_19384"\n "contig_22368"\n "contig_2784" \n')),(0,l.kt)("h2",{id:"view-genotypes"},"View genotypes"),(0,l.kt)("h3",{id:"all-genotypes-in-one-locus-or-sample"},"all genotypes in one locus or sample"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"genotypes(data::PopData, samplelocus::String)\n")),(0,l.kt)("p",null,"Returns a vector (view) of genotypes for a locus, or sample, depending on which the function finds in your data. Don't worry too much about\nthe wild type signature of the return vector. "),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> genotypes(sharks, "contig_2784")\n212-element view(::PooledArrays.PooledVector{Union{Missing, Tuple{Int8, Int8}}, UInt8, Vector{UInt8}}, [468097, 468098, 468099, 468100, 468101, 468102, 468103, 468104, 468105, 468106  \u2026  468299, 468300, 468301, 468302, 468303, 468304, 468305, 468306, 468307, 468308]) with eltype Union{Missing, Tuple{Int8, Int8}}:\n (1, 1)\n (1, 1)\n (1, 1)\n \u22ee\n (1, 1)\n (1, 1)\n (1, 1)\n\n\njulia> genotypes(sharks, "cc_001")\n2209-element view(::PooledArrays.PooledVector{Union{Missing, Tuple{Int8, Int8}}, UInt8, Vector{UInt8}}, [1, 213, 425, 637, 849, 1061, 1273, 1485, 1697, 1909  \u2026  466189, 466401, 466613, 466825, 467037, 467249, 467461, 467673, 467885, 468097]) with eltype Union{Missing, Tuple{Int8, Int8}}:\n (1, 2)\n (1, 1)\n (1, 2)\n \u22ee\n (2, 2)\n (1, 1)\n (1, 1)\n')),(0,l.kt)("h3",{id:"one-sample-one-locus"},"one sample, one locus"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"genotype(data::PopData, sample::String => locus::String)\n")),(0,l.kt)("p",null,"Returns the genotype of the ",(0,l.kt)("inlineCode",{parentName:"p"},"sample")," at the ",(0,l.kt)("inlineCode",{parentName:"p"},"locus"),". Uses ",(0,l.kt)("inlineCode",{parentName:"p"},"Pair")," notation."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> genotype(sharks, "cc_001" => "contig_2784")\n(1, 1)\n')),(0,l.kt)("h3",{id:"many-samples-one-locus"},"many samples, one locus"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"genotype(data::PopData, samples::Vector{String} => loci::String)\n")),(0,l.kt)("p",null,"Returns a subdataframe of the genotypes of the ",(0,l.kt)("inlineCode",{parentName:"p"},"samples")," at the ",(0,l.kt)("inlineCode",{parentName:"p"},"locus"),". Uses ",(0,l.kt)("inlineCode",{parentName:"p"},"Pair")," notation."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> genotypes(sharks, samplenames(sharks)[1:3] => "contig_2784")\n3\xd74 SubDataFrame\n Row \u2502 name     population     locus        genotype \n     \u2502 String7  String         String       Tuple\u2026?  \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 cc_001   CapeCanaveral  contig_2784  (1, 1)\n   2 \u2502 cc_002   CapeCanaveral  contig_2784  (1, 1)\n   3 \u2502 cc_003   CapeCanaveral  contig_2784  (1, 1)\n')),(0,l.kt)("h3",{id:"one-sample-many-loci"},"one sample, many loci"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"genotype(data::PopData, sample::String => loci::Vector{String})\n")),(0,l.kt)("p",null,"Returns a subdataframe of the genotypes of the ",(0,l.kt)("inlineCode",{parentName:"p"},"sample")," at the ",(0,l.kt)("inlineCode",{parentName:"p"},"loci"),". Uses ",(0,l.kt)("inlineCode",{parentName:"p"},"Pair")," notation."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> genotypes(sharks, "cc_001" => loci(sharks)[1:3])\n3\xd74 SubDataFrame\n Row \u2502 name     population     locus         genotype \n     \u2502 String7  String         String        Tuple\u2026?  \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 cc_001   CapeCanaveral  contig_35208  (1, 2)\n   2 \u2502 cc_001   CapeCanaveral  contig_23109  (1, 1)\n   3 \u2502 cc_001   CapeCanaveral  contig_4493   (1, 2)\n')),(0,l.kt)("h3",{id:"many-samples-many-loci"},"many samples, many loci"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"genotype(data::PopData, samples::Vector{String} => loci::Vector{String})\n")),(0,l.kt)("p",null,"Returns a subdataframe of the genotypes of the ",(0,l.kt)("inlineCode",{parentName:"p"},"samples")," at the ",(0,l.kt)("inlineCode",{parentName:"p"},"loci"),". Uses ",(0,l.kt)("inlineCode",{parentName:"p"},"Pair")," notation."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> genotypes(sharks, samplenames(sharks)[1:3] => loci(sharks)[1:3])\n9\xd74 SubDataFrame\n Row \u2502 name     population     locus         genotype \n     \u2502 String7  String         String        Tuple\u2026?  \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 cc_001   CapeCanaveral  contig_35208  (1, 2)\n   2 \u2502 cc_002   CapeCanaveral  contig_35208  (1, 2)\n   3 \u2502 cc_003   CapeCanaveral  contig_35208  (1, 1)\n   4 \u2502 cc_001   CapeCanaveral  contig_23109  (1, 1)\n   5 \u2502 cc_002   CapeCanaveral  contig_23109  (1, 2)\n   6 \u2502 cc_003   CapeCanaveral  contig_23109  missing  \n   7 \u2502 cc_001   CapeCanaveral  contig_4493   (1, 2)\n   8 \u2502 cc_002   CapeCanaveral  contig_4493   (1, 1)\n   9 \u2502 cc_003   CapeCanaveral  contig_4493   (1, 1)\n")))}g.isMDXComponent=!0}}]);