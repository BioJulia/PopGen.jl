"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[9503],{3905:function(e,a,t){t.d(a,{Zo:function(){return m},kt:function(){return h}});var n=t(7294);function r(e,a,t){return a in e?Object.defineProperty(e,a,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[a]=t,e}function i(e,a){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);a&&(n=n.filter((function(a){return Object.getOwnPropertyDescriptor(e,a).enumerable}))),t.push.apply(t,n)}return t}function l(e){for(var a=1;a<arguments.length;a++){var t=null!=arguments[a]?arguments[a]:{};a%2?i(Object(t),!0).forEach((function(a){r(e,a,t[a])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):i(Object(t)).forEach((function(a){Object.defineProperty(e,a,Object.getOwnPropertyDescriptor(t,a))}))}return e}function s(e,a){if(null==e)return{};var t,n,r=function(e,a){if(null==e)return{};var t,n,r={},i=Object.keys(e);for(n=0;n<i.length;n++)t=i[n],a.indexOf(t)>=0||(r[t]=e[t]);return r}(e,a);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(n=0;n<i.length;n++)t=i[n],a.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(r[t]=e[t])}return r}var o=n.createContext({}),u=function(e){var a=n.useContext(o),t=a;return e&&(t="function"==typeof e?e(a):l(l({},a),e)),t},m=function(e){var a=u(e.components);return n.createElement(o.Provider,{value:a},e.children)},p="mdxType",c={inlineCode:"code",wrapper:function(e){var a=e.children;return n.createElement(n.Fragment,{},a)}},d=n.forwardRef((function(e,a){var t=e.components,r=e.mdxType,i=e.originalType,o=e.parentName,m=s(e,["components","mdxType","originalType","parentName"]),p=u(t),d=r,h=p["".concat(o,".").concat(d)]||p[d]||c[d]||i;return t?n.createElement(h,l(l({ref:a},m),{},{components:t})):n.createElement(h,l({ref:a},m))}));function h(e,a){var t=arguments,r=a&&a.mdxType;if("string"==typeof e||r){var i=t.length,l=new Array(i);l[0]=d;var s={};for(var o in a)hasOwnProperty.call(a,o)&&(s[o]=a[o]);s.originalType=e,s[p]="string"==typeof e?e:r,l[1]=s;for(var u=2;u<i;u++)l[u]=t[u];return n.createElement.apply(null,l)}return n.createElement.apply(null,t)}d.displayName="MDXCreateElement"},5162:function(e,a,t){t.d(a,{Z:function(){return l}});var n=t(7294),r=t(6010),i={tabItem:"tabItem_Ymn6"};function l(e){var a=e.children,t=e.hidden,l=e.className;return n.createElement("div",{role:"tabpanel",className:(0,r.Z)(i.tabItem,l),hidden:t},a)}},4866:function(e,a,t){t.d(a,{Z:function(){return w}});var n=t(7462),r=t(7294),i=t(6010),l=t(2466),s=t(6550),o=t(1980),u=t(7392),m=t(12);function p(e){return function(e){var a,t;return null!=(a=null==(t=r.Children.map(e,(function(e){if(!e||(0,r.isValidElement)(e)&&(a=e.props)&&"object"==typeof a&&"value"in a)return e;var a;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})))?void 0:t.filter(Boolean))?a:[]}(e).map((function(e){var a=e.props;return{value:a.value,label:a.label,attributes:a.attributes,default:a.default}}))}function c(e){var a=e.values,t=e.children;return(0,r.useMemo)((function(){var e=null!=a?a:p(t);return function(e){var a=(0,u.l)(e,(function(e,a){return e.value===a.value}));if(a.length>0)throw new Error('Docusaurus error: Duplicate values "'+a.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.')}(e),e}),[a,t])}function d(e){var a=e.value;return e.tabValues.some((function(e){return e.value===a}))}function h(e){var a=e.queryString,t=void 0!==a&&a,n=e.groupId,i=(0,s.k6)(),l=function(e){var a=e.queryString,t=void 0!==a&&a,n=e.groupId;if("string"==typeof t)return t;if(!1===t)return null;if(!0===t&&!n)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return null!=n?n:null}({queryString:t,groupId:n});return[(0,o._X)(l),(0,r.useCallback)((function(e){if(l){var a=new URLSearchParams(i.location.search);a.set(l,e),i.replace(Object.assign({},i.location,{search:a.toString()}))}}),[l,i])]}function f(e){var a,t,n,i,l=e.defaultValue,s=e.queryString,o=void 0!==s&&s,u=e.groupId,p=c(e),f=(0,r.useState)((function(){return function(e){var a,t=e.defaultValue,n=e.tabValues;if(0===n.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(t){if(!d({value:t,tabValues:n}))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+t+'" but none of its children has the corresponding value. Available values are: '+n.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");return t}var r=null!=(a=n.find((function(e){return e.default})))?a:n[0];if(!r)throw new Error("Unexpected error: 0 tabValues");return r.value}({defaultValue:l,tabValues:p})})),k=f[0],g=f[1],b=h({queryString:o,groupId:u}),v=b[0],y=b[1],w=(a=function(e){return e?"docusaurus.tab."+e:null}({groupId:u}.groupId),t=(0,m.Nk)(a),n=t[0],i=t[1],[n,(0,r.useCallback)((function(e){a&&i.set(e)}),[a,i])]),T=w[0],N=w[1],j=function(){var e=null!=v?v:T;return d({value:e,tabValues:p})?e:null}();return(0,r.useLayoutEffect)((function(){j&&g(j)}),[j]),{selectedValue:k,selectValue:(0,r.useCallback)((function(e){if(!d({value:e,tabValues:p}))throw new Error("Can't select invalid tab value="+e);g(e),y(e),N(e)}),[y,N,p]),tabValues:p}}var k=t(2389),g={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};function b(e){var a=e.className,t=e.block,s=e.selectedValue,o=e.selectValue,u=e.tabValues,m=[],p=(0,l.o5)().blockElementScrollPositionUntilNextRender,c=function(e){var a=e.currentTarget,t=m.indexOf(a),n=u[t].value;n!==s&&(p(a),o(n))},d=function(e){var a,t=null;switch(e.key){case"Enter":c(e);break;case"ArrowRight":var n,r=m.indexOf(e.currentTarget)+1;t=null!=(n=m[r])?n:m[0];break;case"ArrowLeft":var i,l=m.indexOf(e.currentTarget)-1;t=null!=(i=m[l])?i:m[m.length-1]}null==(a=t)||a.focus()};return r.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,i.Z)("tabs",{"tabs--block":t},a)},u.map((function(e){var a=e.value,t=e.label,l=e.attributes;return r.createElement("li",(0,n.Z)({role:"tab",tabIndex:s===a?0:-1,"aria-selected":s===a,key:a,ref:function(e){return m.push(e)},onKeyDown:d,onClick:c},l,{className:(0,i.Z)("tabs__item",g.tabItem,null==l?void 0:l.className,{"tabs__item--active":s===a})}),null!=t?t:a)})))}function v(e){var a=e.lazy,t=e.children,n=e.selectedValue,i=(Array.isArray(t)?t:[t]).filter(Boolean);if(a){var l=i.find((function(e){return e.props.value===n}));return l?(0,r.cloneElement)(l,{className:"margin-top--md"}):null}return r.createElement("div",{className:"margin-top--md"},i.map((function(e,a){return(0,r.cloneElement)(e,{key:a,hidden:e.props.value!==n})})))}function y(e){var a=f(e);return r.createElement("div",{className:(0,i.Z)("tabs-container",g.tabList)},r.createElement(b,(0,n.Z)({},e,a)),r.createElement(v,(0,n.Z)({},e,a)))}function w(e){var a=(0,k.Z)();return r.createElement(y,(0,n.Z)({key:String(a)},e))}},9720:function(e,a,t){t.r(a),t.d(a,{assets:function(){return c},contentTitle:function(){return m},default:function(){return k},frontMatter:function(){return u},metadata:function(){return p},toc:function(){return d}});var n=t(7462),r=t(3366),i=(t(7294),t(3905)),l=t(4866),s=t(5162),o=["components"],u={id:"comparison",title:"Comparison",sidebar_label:"Comparison"},m=void 0,p={unversionedId:"gettingstarted/comparison",id:"gettingstarted/comparison",title:"Comparison",description:"There's a reason we started investing so many hours and so many new grey hairs into writing PopGen.jl when there was an existing ecosystem in R to perform these same tasks. Like we explain in the home page of these docs, we want a platform that is:",source:"@site/docs/gettingstarted/comparison.md",sourceDirName:"gettingstarted",slug:"/gettingstarted/comparison",permalink:"/PopGen.jl/docs/gettingstarted/comparison",draft:!1,editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/gettingstarted/comparison.md",tags:[],version:"current",lastUpdatedAt:1659126779,formattedLastUpdatedAt:"Jul 29, 2022",frontMatter:{id:"comparison",title:"Comparison",sidebar_label:"Comparison"},sidebar:"docs",previous:{title:"PopGen.jl tips",permalink:"/PopGen.jl/docs/gettingstarted/tips"},next:{title:"Provided datasets",permalink:"/PopGen.jl/docs/gettingstarted/datasets"}},c={},d=[{value:"Benchmarks",id:"benchmarks",level:2},{value:"Reading in data",id:"reading-in-data",level:3},{value:"<code>PopData</code> vs <code>genind</code> size",id:"popdata-vs-genind-size",level:3},{value:"Summary statistics",id:"summary-statistics",level:3},{value:"Chi-squared test for HWE",id:"chi-squared-test-for-hwe",level:3},{value:"Pairwise FST",id:"pairwise-fst",level:3}],h={toc:d},f="wrapper";function k(e){var a=e.components,u=(0,r.Z)(e,o);return(0,i.kt)(f,(0,n.Z)({},h,u,{components:a,mdxType:"MDXLayout"}),(0,i.kt)("p",null,"There's a reason we started investing so many hours and so many new grey hairs into writing PopGen.jl when there was an existing ecosystem in R to perform these same tasks. Like we explain in the home page of these docs, we want a platform that is:"),(0,i.kt)("ol",null,(0,i.kt)("li",{parentName:"ol"},"fast(er)"),(0,i.kt)("li",{parentName:"ol"},"written in a single language"),(0,i.kt)("li",{parentName:"ol"},"easy to use")),(0,i.kt)("p",null,"So, we'd like to prove that Julia and PopGen.jl actually achieves that by showing a few benchmarks comparing PopGen.jl to popular population genetics packages in R. It's worth mentioning that we ourselves use and have published work incorporating these packages, and are incredibly grateful for the work invested in them. We appreciate those folks and have tremendous respect and envy for the work they continue to do! Here are links to ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/thibautjombart/adegenet"},"adegenet"),", ",(0,i.kt)("a",{parentName:"p",href:"https://academic.oup.com/bioinformatics/article/26/3/419/215731/"},"pegas"),", ",(0,i.kt)("a",{parentName:"p",href:"https://cran.r-project.org/web/packages/hierfstat/index.html"},"hierfstat"),", and ",(0,i.kt)("a",{parentName:"p",href:"https://cran.r-project.org/package=ape"},"ape"),".  "),(0,i.kt)("h2",{id:"benchmarks"},"Benchmarks"),(0,i.kt)("p",null,"To make this a practical comparison, we're going to use the ",(0,i.kt)("inlineCode",{parentName:"p"},"gulfsharks")," data because it is considerably larger (212 samples x 2209 loci) than ",(0,i.kt)("inlineCode",{parentName:"p"},"nancycats"),' (237 x 9) and a bit more of a "stress test".  All benchmarks in R are performed using the ',(0,i.kt)("inlineCode",{parentName:"p"},"microbenchmark")," package, and  ",(0,i.kt)("inlineCode",{parentName:"p"},"BenchmarkTools")," are used for Julia."),(0,i.kt)(l.Z,{block:!0,defaultValue:"j",values:[{label:"load Julia packages",value:"j"},{label:"load R packages",value:"r"}],mdxType:"Tabs"},(0,i.kt)(s.Z,{value:"j",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"using BenchmarkTools, PopGen\n"))),(0,i.kt)(s.Z,{value:"r",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-r"},"library(adegenet)\nlibrary(pegas)\nlibrary(hierfstat)\nlibrary(microbenchmark)\n")))),(0,i.kt)("p",null,"As a note, the reported benchmarks are being performed on a 64-bit Manjaro Linux system on a nothing-special Lenovo Thinkbook 14S  with 8gigs of RAM and a 8th gen Intel i5 mobile processor. ",(0,i.kt)("strong",{parentName:"p"},"Note:")," all of the Julia benchmarks, unless explicitly stated, are performed single-threaded (i.e. not parallel, distributed, or GPU)."),(0,i.kt)("h3",{id:"reading-in-data"},"Reading in data"),(0,i.kt)("p",null,"Since ",(0,i.kt)("inlineCode",{parentName:"p"},"gulfsharks")," is provided in PopGenCore.jl, we will just load it with ",(0,i.kt)("inlineCode",{parentName:"p"},"genepop()"),".  If you would like to try this yourself in R, find the ",(0,i.kt)("inlineCode",{parentName:"p"},"gulfsharks.gen")," file in the package repository under ",(0,i.kt)("inlineCode",{parentName:"p"},"/data/source/gulfsharks.gen"),". Since the file importer now uses CSV.jl to read files, there is one step of the genepop parser that is multithreaded. However, the majority of the data parsing (formatting the raw data into a correct PopData structure) occurs using a single thread. This R benchmark will take a few minutes. Consider making some tea while you wait."),(0,i.kt)(l.Z,{block:!0,defaultValue:"j",values:[{label:"Julia",value:"j"},{label:"R",value:"r"}],mdxType:"Tabs"},(0,i.kt)(s.Z,{value:"j",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},'@benchmark sharks = genepop("data/source/gulfsharks.gen", silent = true)\nBenchmarkTools.Trial: 10 samples with 1 evaluation.\n Range (min \u2026 max):  472.671 ms \u2026 526.777 ms  \u250a GC (min \u2026 max): 0.00% \u2026 5.75%\n Time  (median):     507.187 ms               \u250a GC (median):    5.00%\n Time  (mean \xb1 \u03c3):   501.984 ms \xb1  17.301 ms  \u250a GC (mean \xb1 \u03c3):  3.47% \xb1 2.99%\n\n  \u2581          \u2581     \u2581   \u2581            \u2581        \u2581 \u2588   \u2581          \u2581  \n  \u2588\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2588\u2581\u2581\u2581\u2581\u2581\u2588\u2581\u2581\u2581\u2588\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2588\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2588\u2581\u2588\u2581\u2581\u2581\u2588\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2581\u2588 \u2581\n  473 ms           Histogram: frequency by time          527 ms <\n\n Memory estimate: 172.46 MiB, allocs estimate: 2246870.\n'))),(0,i.kt)(s.Z,{value:"r",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-r"},'> microbenchmark(read.genepop(file = "/home/pdimens/gulfsharks.gen", ncode = 3L, quiet = TRUE))\nUnit: seconds\n read.genepop(file = "/home/pdimens/gulfsharks.gen", ncode = 3L, quiet = FALSE)\n      min       lq     mean   median       uq      max neval\n 5.670637 6.218719 6.745065 6.387936 7.019667 9.173005   100\nmicrobenchmark)\n')))),(0,i.kt)("p",null,(0,i.kt)("img",{alt:"import plot",src:t(665).Z,width:"745",height:"372"})),(0,i.kt)("p",null,"Comparing averages, Julia (via PopGen.jl) clocks in at ",(0,i.kt)("inlineCode",{parentName:"p"},"507ms")," versus R (via adegenet) at ",(0,i.kt)("inlineCode",{parentName:"p"},"6.745s")," , so ~13.3x faster."),(0,i.kt)("h3",{id:"popdata-vs-genind-size"},(0,i.kt)("inlineCode",{parentName:"h3"},"PopData")," vs ",(0,i.kt)("inlineCode",{parentName:"h3"},"genind")," size"),(0,i.kt)("p",null,'It was pretty tricky to come up with a sensible/efficient/convenient data structure for PopGen.jl, and while the two-DataFrames design might not seem like it took a lot of effort, we ultimately decided that the column-major style and available tools, combined with careful genotype Typing was a decent "middle-ground" of ease-of-use vs performance.'),(0,i.kt)("p",null,(0,i.kt)("em",{parentName:"p"},"Anyway"),", it's important to understand how much space your data will take up in memory (your RAM) when you load it in, especially since data's only getting bigger! Keep in mind that ",(0,i.kt)("inlineCode",{parentName:"p"},"gulfsharks")," in PopGen.jl also provides lat/long data, which ",(0,i.kt)("em",{parentName:"p"},"should")," inflate the size of the object somewhat compared to the ",(0,i.kt)("inlineCode",{parentName:"p"},"genind"),", which we won't add any location data to."),(0,i.kt)(l.Z,{block:!0,defaultValue:"j",values:[{label:"Julia",value:"j"},{label:"R",value:"r"}],mdxType:"Tabs"},(0,i.kt)(s.Z,{value:"j",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"julia> Base.summarysize(sharks)\n3452870\n#bytes\n"))),(0,i.kt)(s.Z,{value:"r",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-r"},"> object.size(gen)\n5331536 bytes\n")))),(0,i.kt)("p",null,(0,i.kt)("img",{alt:"data structure plot",src:t(5268).Z,width:"759",height:"376"})),(0,i.kt)("p",null,"The original genepop file is ",(0,i.kt)("inlineCode",{parentName:"p"},"3.2mb")," (the vertical line), and our ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," object takes up ~",(0,i.kt)("inlineCode",{parentName:"p"},"3.45mb")," in memory (300kb larger than the source file) versus the ~",(0,i.kt)("inlineCode",{parentName:"p"},"5.3mb")," of a ",(0,i.kt)("inlineCode",{parentName:"p"},"genind"),", which is ~1.5x larger than the source file. That's quite a big difference!"),(0,i.kt)("h3",{id:"summary-statistics"},"Summary statistics"),(0,i.kt)("p",null,"The obvious hallmark of population genetics is heterozygosity values and F-statistics. Here we'll compare the basic summary statistics that can be produced using ",(0,i.kt)("inlineCode",{parentName:"p"},"hierfstat")," and ",(0,i.kt)("inlineCode",{parentName:"p"},"PopGen.jl"),"."),(0,i.kt)(l.Z,{block:!0,defaultValue:"j",values:[{label:"Julia",value:"j"},{label:"R",value:"r"}],mdxType:"Tabs"},(0,i.kt)(s.Z,{value:"j",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'julia> @benchmark summary(sharks, by = "global")\nBenchmarkTools.Trial:\n  memory estimate:  88.42 MiB\n  allocs estimate:  1307128\n  --------------\n  minimum time:     151.963 ms (0.00% GC)\n  median time:      171.484 ms (7.60% GC)\n  mean time:        172.456 ms (6.08% GC)\n  maximum time:     186.606 ms (7.04% GC)\n  --------------\n  samples:          29\n  evals/sample:     1\n'))),(0,i.kt)(s.Z,{value:"r",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-r"},"> microbenchmark(basic.stats(gen))\nUnit: seconds\n             expr      min       lq     mean\n basic.stats(gen) 4.276996 4.425934 4.618796\n   median       uq      max neval\n 4.609901 4.706666 5.292831   100\n")))),(0,i.kt)("p",null,(0,i.kt)("img",{alt:"summary statistics plot",src:t(2342).Z,width:"744",height:"372"})),(0,i.kt)("p",null,"Comparing averages, PopGen.jl clocks in at ~",(0,i.kt)("inlineCode",{parentName:"p"},"171ms")," versus hierfstat's ",(0,i.kt)("inlineCode",{parentName:"p"},"4.6s"),", which is ~",(0,i.kt)("strong",{parentName:"p"},"27x")," faster on these data. However, when testing on a data that was 401 samples x 5331 loci (not shown), PopGen.jl performed 36.6x faster. This gap seems to increase the larger the data is, but we have not tested the upper limits of this."),(0,i.kt)("h3",{id:"chi-squared-test-for-hwe"},"Chi-squared test for HWE"),(0,i.kt)("p",null,"This is a classic population genetics test and a relatively simple one. The R benchmark will take a while again, so if you're following along, this would be a good time to reconnect with an old friend."),(0,i.kt)(l.Z,{block:!0,defaultValue:"j",values:[{label:"Julia",value:"j"},{label:"R",value:"r"}],mdxType:"Tabs"},(0,i.kt)(s.Z,{value:"j",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"julia> @benchmark hwe_test(sharks)\nBenchmarkTools.Trial:\n  memory estimate:  46.22 MiB\n  allocs estimate:  1050525\n  --------------\n  minimum time:     145.476 ms (0.00% GC)\n  median time:      179.146 ms (4.35% GC)\n  mean time:        176.142 ms (3.56% GC)\n  maximum time:     204.813 ms (0.00% GC)\n  --------------\n  samples:          29\n  evals/sample:     1\n"))),(0,i.kt)(s.Z,{value:"r",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-r"},"> microbenchmark(hw.test(gen, B = 0))\nUnit: seconds\n                expr      min       lq     mean   median       uq      max neval\n hw.test(gen, B = 0) 5.100298 5.564807 6.265948 5.878842 6.917006 8.815179   100\n")))),(0,i.kt)("p",null,(0,i.kt)("img",{alt:"chi squared plot",src:t(3155).Z,width:"745",height:"372"})),(0,i.kt)("p",null,"Comparing averages, PopGen.jl clocks in at ~",(0,i.kt)("inlineCode",{parentName:"p"},"176ms")," versus adegenet's ",(0,i.kt)("inlineCode",{parentName:"p"},"6.3s"),", so ~",(0,i.kt)("strong",{parentName:"p"},"35.8x")," faster on these data(!)"),(0,i.kt)("h3",{id:"pairwise-fst"},"Pairwise FST"),(0,i.kt)("p",null,"You all know it, you all love it. What's population genetics without a little pairwise FST sprinkled in? This benchmark compairs the Weir & Cockerham pairwise FST calculation in ",(0,i.kt)("inlineCode",{parentName:"p"},"PopGen.jl")," against ",(0,i.kt)("inlineCode",{parentName:"p"},"hierfstat")),(0,i.kt)(l.Z,{block:!0,defaultValue:"j",values:[{label:"Julia",value:"j"},{label:"Julia (parallel)",value:"jp"},{label:"R",value:"r"}],mdxType:"Tabs"},(0,i.kt)(s.Z,{value:"j",mdxType:"TabItem"},(0,i.kt)("p",null,"We will add the extra keywords ",(0,i.kt)("inlineCode",{parentName:"p"},"samples")," and ",(0,i.kt)("inlineCode",{parentName:"p"},"seconds")," to the benchmark\nmacro so we can get a full 100 evaluations. You will need to start Julia with 1 available threads via ",(0,i.kt)("inlineCode",{parentName:"p"},"julia --threads 1")," (julia >= v1.5) or ",(0,i.kt)("inlineCode",{parentName:"p"},"JULIA_NUM_THREADS=1")," (< v1.5)."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"julia> @benchmark pairwise_fst(sharks) samples = 100 seconds = 700\nBenchmarkTools.Trial: \n  memory estimate:  869.93 MiB\n  allocs estimate:  6090633\n  --------------\n  minimum time:     557.995 ms (9.29% GC)\n  median time:      571.297 ms (11.26% GC)\n  mean time:        580.627 ms (12.59% GC)\n  maximum time:     754.451 ms (31.41% GC)\n  --------------\n  samples:          100\n  evals/sample:     1\n\n"))),(0,i.kt)(s.Z,{value:"jp",mdxType:"TabItem"},(0,i.kt)("p",null,"This is to demonstrate what the speed is like when starting Julia with 6 available threads via ",(0,i.kt)("inlineCode",{parentName:"p"},"julia --threads 6")," (julia >= v1.5) or ",(0,i.kt)("inlineCode",{parentName:"p"},"JULIA_NUM_THREADS=6")," (< v1.5)."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"julia> @benchmark pairwise_fst(sharks) samples = 100 seconds = 700\nBenchmarkTools.Trial: \n  memory estimate:  869.93 MiB\n  allocs estimate:  6090639\n  --------------\n  minimum time:     205.038 ms (0.00% GC)\n  median time:      305.189 ms (28.90% GC)\n  mean time:        299.227 ms (25.03% GC)\n  maximum time:     359.663 ms (27.05% GC)\n  --------------\n  samples:          100\n  evals/sample:     1\n"))),(0,i.kt)(s.Z,{value:"r",mdxType:"TabItem"},(0,i.kt)("p",null,"We'll need to convert ",(0,i.kt)("inlineCode",{parentName:"p"},"sharks")," into the matrix/dataframe ",(0,i.kt)("inlineCode",{parentName:"p"},"hierfstat")," needs\nto run this calculation. The conversion will be a separate step so as not\nto add unnecessary (or unfair) overhead to the benchmark. This benchmark is\ngoing to take ",(0,i.kt)("strong",{parentName:"p"},"forever")," (200s/run x 100 runs = 5.5hrs), so if you absolutely insist on\ntrying it out yourself, you may want to pop outside and enjoy some fresh\nair for a bit (I ran it overnight). Seriously, you don't want to watch this paint dry \ud83d\udd8c\ufe0f."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-r"},"> sharks_hierf <- genind2hierfstat(sharks)\n\n> microbenchmark(pairwise.WCfst(sharks_hierf))\nUnit: seconds\n                    expr      min       lq     mean   median       uq      max neval\n pairwise.WCfst(sharks2) 192.2786 192.9277 199.4861 193.5743 195.0079 301.6879   100\n")))),(0,i.kt)("p",null,(0,i.kt)("img",{alt:"pairwise fst plot",src:t(3584).Z,width:"835",height:"372"})),(0,i.kt)("p",null,"On a single thread, pairwise FST in ",(0,i.kt)("inlineCode",{parentName:"p"},"PopGen.jl")," is ",(0,i.kt)("strong",{parentName:"p"},"~340x")," faster than in ",(0,i.kt)("inlineCode",{parentName:"p"},"hierfstat"),", and a whopping ",(0,i.kt)("strong",{parentName:"p"},"665x")," faster using 6 threads with the optimized matrix-based implementation."))}k.isMDXComponent=!0},3155:function(e,a,t){a.Z=t.p+"assets/images/chisqplot-94c0d817b3b77543f63ae472bbde4d6a.png"},3584:function(e,a,t){a.Z=t.p+"assets/images/fstplot-6a672cb5538833892ee7eda34c33645a.png"},5268:function(e,a,t){a.Z=t.p+"assets/images/objectplot-ec422864760c5133e262eec866a62f7c.png"},665:function(e,a,t){a.Z=t.p+"assets/images/speedplot-a0dcb588e1923961a9aaf9e183490011.png"},2342:function(e,a,t){a.Z=t.p+"assets/images/sumstatplot-598f7c700e09791e891c6f574ce057fc.png"}}]);