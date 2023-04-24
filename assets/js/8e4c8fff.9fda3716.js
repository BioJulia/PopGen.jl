"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[2389],{3905:function(e,t,n){n.d(t,{Zo:function(){return u},kt:function(){return f}});var o=n(7294);function i(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function a(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);t&&(o=o.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,o)}return n}function r(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?a(Object(n),!0).forEach((function(t){i(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):a(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,o,i=function(e,t){if(null==e)return{};var n,o,i={},a=Object.keys(e);for(o=0;o<a.length;o++)n=a[o],t.indexOf(n)>=0||(i[n]=e[n]);return i}(e,t);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(o=0;o<a.length;o++)n=a[o],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(i[n]=e[n])}return i}var s=o.createContext({}),p=function(e){var t=o.useContext(s),n=t;return e&&(n="function"==typeof e?e(t):r(r({},t),e)),n},u=function(e){var t=p(e.components);return o.createElement(s.Provider,{value:t},e.children)},m="mdxType",c={inlineCode:"code",wrapper:function(e){var t=e.children;return o.createElement(o.Fragment,{},t)}},d=o.forwardRef((function(e,t){var n=e.components,i=e.mdxType,a=e.originalType,s=e.parentName,u=l(e,["components","mdxType","originalType","parentName"]),m=p(n),d=i,f=m["".concat(s,".").concat(d)]||m[d]||c[d]||a;return n?o.createElement(f,r(r({ref:t},u),{},{components:n})):o.createElement(f,r({ref:t},u))}));function f(e,t){var n=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var a=n.length,r=new Array(a);r[0]=d;var l={};for(var s in t)hasOwnProperty.call(t,s)&&(l[s]=t[s]);l.originalType=e,l[m]="string"==typeof e?e:i,r[1]=l;for(var p=2;p<a;p++)r[p]=n[p];return o.createElement.apply(null,r)}return o.createElement.apply(null,n)}d.displayName="MDXCreateElement"},5245:function(e,t,n){n.r(t),n.d(t,{assets:function(){return u},contentTitle:function(){return s},default:function(){return f},frontMatter:function(){return l},metadata:function(){return p},toc:function(){return m}});var o=n(7462),i=n(3366),a=(n(7294),n(3905)),r=["components"],l={id:"tips",title:"PopGen.jl tips",sidebar_label:"PopGen.jl tips"},s=void 0,p={unversionedId:"gettingstarted/tips",id:"gettingstarted/tips",title:"PopGen.jl tips",description:"Here are some useful tips to getting comfortable with PopGen.jl.",source:"@site/docs/gettingstarted/tips.md",sourceDirName:"gettingstarted",slug:"/gettingstarted/tips",permalink:"/PopGen.jl/docs/gettingstarted/tips",draft:!1,editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/gettingstarted/tips.md",tags:[],version:"current",lastUpdatedAt:1659126598,formattedLastUpdatedAt:"Jul 29, 2022",frontMatter:{id:"tips",title:"PopGen.jl tips",sidebar_label:"PopGen.jl tips"},sidebar:"docs",previous:{title:"A quick Julia primer",permalink:"/PopGen.jl/docs/gettingstarted/juliaprimer"},next:{title:"Comparison",permalink:"/PopGen.jl/docs/gettingstarted/comparison"}},u={},m=[{value:"PopGen vs PopGenCore",id:"popgen-vs-popgencore",level:3},{value:"Internal functions",id:"internal-functions",level:3},{value:"Function names",id:"function-names",level:3},{value:"Argument names",id:"argument-names",level:3},{value:"Dev tips",id:"dev-tips",level:3}],c={toc:m},d="wrapper";function f(e){var t=e.components,n=(0,i.Z)(e,r);return(0,a.kt)(d,(0,o.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,a.kt)("p",null,"Here are some useful tips to getting comfortable with PopGen.jl."),(0,a.kt)("h3",{id:"popgen-vs-popgencore"},"PopGen vs PopGenCore"),(0,a.kt)("admonition",{title:"TLDR",type:"tip"},(0,a.kt)("p",{parentName:"admonition"},"PopGen is for analyses, PopGenCore is for everything else")),(0,a.kt)("p",null,"PopGen.jl is now intended exclusively for population genetic analyses, and PopGenCore.jl is the \"core\" package\nwhere just about every utility lives relating to working with PopData that isn't an analysis. The recent split\nof PopGen and PopGenCore means PopGenCore is now a great standalone package for data viewing\nand manipulation. If you don't need higher order analyses (which is what PopGen.jl provides), then PopGenCore.jl\nshould be enough. "),(0,a.kt)("h3",{id:"internal-functions"},"Internal functions"),(0,a.kt)("admonition",{title:"TLDR",type:"tip"},(0,a.kt)("p",{parentName:"admonition"},"functions you aren't expected to use start with underscores ",(0,a.kt)("inlineCode",{parentName:"p"},"_"))),(0,a.kt)("p",null,"It's pretty common to have a series of user-facing functions in a package (the ones that get exported)\nalong with a series of unexported ones that are useful for development or are helper functions for the\nexported ones. If you see a function that starts with an underscore, like ",(0,a.kt)("inlineCode",{parentName:"p"},"_adjacency_matrix"),", then you\n(as a user) aren't expected to know or worry about it, let alone use it. You totally can, but\nnot all of them have docstrings."),(0,a.kt)("h3",{id:"function-names"},"Function names"),(0,a.kt)("admonition",{title:"TLDR",type:"tip"},(0,a.kt)("p",{parentName:"admonition"},"all user functions in PopGen.jl (not including PopGenCore.jl) are lowercase and use no underscores")),(0,a.kt)("p",null,"In an effort to be consistent with the Julia language and just consistent overall, user-facing functions\nare named with three main principals:"),(0,a.kt)("ol",null,(0,a.kt)("li",{parentName:"ol"},"as verbose of a name as is reasonable"),(0,a.kt)("li",{parentName:"ol"},"no underscores between words"),(0,a.kt)("li",{parentName:"ol"},"all lowercase, always, unless it's a DataType",(0,a.kt)("ul",{parentName:"li"},(0,a.kt)("li",{parentName:"ul"},"or method keyword option (like the ",(0,a.kt)("inlineCode",{parentName:"li"},"kinship")," or series)")))),(0,a.kt)("p",null,"These decisions are because of years-long frustration with the swirling chaos that is R function names.\nWhere it is possible and reasonable to do so, function names are verbose and descriptive. For example,\nif you wanted to perform a pairwise FST, the function is called ",(0,a.kt)("inlineCode",{parentName:"p"},"pairwisefst")," -- ",(0,a.kt)("strong",{parentName:"p"},"not")," ",(0,a.kt)("inlineCode",{parentName:"p"},"fst"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"FST"),",\n",(0,a.kt)("inlineCode",{parentName:"p"},"pairwise_fst"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"pairwise_FST"),", or ",(0,a.kt)("inlineCode",{parentName:"p"},"pairwiseFST"),". This is an ongoing effort, so if you spot something\nthat doesn't fit this mold, please submit a Pull Request!"),(0,a.kt)("h3",{id:"argument-names"},"Argument names"),(0,a.kt)("admonition",{title:"TLDR",type:"tip"},(0,a.kt)("p",{parentName:"admonition"},"all user functions in PopGen.jl have descriptive keyword argument names")),(0,a.kt)("p",null,"This too is driven by frustration with the shortened and/or truncated keyword arguments of just about\nevery R function. As if learning another language wasn't enough of a task, memorizing nonintuitive\nkeywords is a completely unneccessary stumbling block. To the best of our abilities, we try to name\nkeywords using ",(0,a.kt)("em",{parentName:"p"},"entire")," words, and hopefully words that are intuitive if you were just trying to guess\nbefore calling up the docstring. "),(0,a.kt)("h3",{id:"dev-tips"},"Dev tips"),(0,a.kt)("admonition",{title:"TLDR",type:"tip"},(0,a.kt)("ul",{parentName:"admonition"},(0,a.kt)("li",{parentName:"ul"},"Avoid allocations when possible"),(0,a.kt)("li",{parentName:"ul"},"First drafts don't need to be efficient or performant"),(0,a.kt)("li",{parentName:"ul"},"When in doubt, roll your own helper functions"))),(0,a.kt)("p",null,"If you plan on extending PopGen.jl, here are some useful tips we learned the hard way:"),(0,a.kt)("ul",null,(0,a.kt)("li",{parentName:"ul"},"Most basic (non-mathematic) or fundamental operations live in PopGenCore.jl"),(0,a.kt)("li",{parentName:"ul"},"Math and analyses live in PopGen.jl"),(0,a.kt)("li",{parentName:"ul"},"It's ok to write something functional but not efficient. Things we can always be improved later!"),(0,a.kt)("li",{parentName:"ul"},"Allocations are the enemy",(0,a.kt)("ul",{parentName:"li"},(0,a.kt)("li",{parentName:"ul"},"Ignore this recommendation if the goal is a first draft or proof of concept"),(0,a.kt)("li",{parentName:"ul"},"Creating Arrays is costly and that cost scales ",(0,a.kt)("strong",{parentName:"li"},"tremendously")," with core operations occurring thousands/millions of times"),(0,a.kt)("li",{parentName:"ul"},"If the goal is to be performant, aim for as few allocations as possible"),(0,a.kt)("li",{parentName:"ul"},"Embrace loops"))),(0,a.kt)("li",{parentName:"ul"},(0,a.kt)("inlineCode",{parentName:"li"},"if")," statements inside loops slow them down",(0,a.kt)("ul",{parentName:"li"},(0,a.kt)("li",{parentName:"ul"},"But if they are small or simple enough, can still be faster and more efficient than allocating Arrays"))),(0,a.kt)("li",{parentName:"ul"},"Sometimes ",(0,a.kt)("inlineCode",{parentName:"li"},"Base")," (or others) don't have the performance you need",(0,a.kt)("ul",{parentName:"li"},(0,a.kt)("li",{parentName:"ul"},"There's no harm in writing a small helper function if it's an improvement"))),(0,a.kt)("li",{parentName:"ul"},"Incrementing a variable inside a loop with ",(0,a.kt)("inlineCode",{parentName:"li"},"+=")," or similar (example 1) is ",(0,a.kt)("strong",{parentName:"li"},"a lot")," more performant that doing so with a list comprehension (example 2), even though it requires more lines of code.")),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-julia"},"# example 1: Incrementing inside for loop\n# performant\nn = 0\nfor i in 1:something\n    n += 1\nend\n\n# example 2: Incrementing inside list comprehension \n# not as performant\nn = 0\n[n+= 1 for i:something]\n")))}f.isMDXComponent=!0}}]);