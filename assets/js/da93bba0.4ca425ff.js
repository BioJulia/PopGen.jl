"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[3969],{3905:function(e,t,n){n.d(t,{Zo:function(){return p},kt:function(){return c}});var a=n(7294);function i(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function o(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function r(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?o(Object(n),!0).forEach((function(t){i(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):o(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,a,i=function(e,t){if(null==e)return{};var n,a,i={},o=Object.keys(e);for(a=0;a<o.length;a++)n=o[a],t.indexOf(n)>=0||(i[n]=e[n]);return i}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(a=0;a<o.length;a++)n=o[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(i[n]=e[n])}return i}var s=a.createContext({}),u=function(e){var t=a.useContext(s),n=t;return e&&(n="function"==typeof e?e(t):r(r({},t),e)),n},p=function(e){var t=u(e.components);return a.createElement(s.Provider,{value:t},e.children)},d={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},m=a.forwardRef((function(e,t){var n=e.components,i=e.mdxType,o=e.originalType,s=e.parentName,p=l(e,["components","mdxType","originalType","parentName"]),m=u(n),c=i,h=m["".concat(s,".").concat(c)]||m[c]||d[c]||o;return n?a.createElement(h,r(r({ref:t},p),{},{components:n})):a.createElement(h,r({ref:t},p))}));function c(e,t){var n=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var o=n.length,r=new Array(o);r[0]=m;var l={};for(var s in t)hasOwnProperty.call(t,s)&&(l[s]=t[s]);l.originalType=e,l.mdxType="string"==typeof e?e:i,r[1]=l;for(var u=2;u<o;u++)r[u]=n[u];return a.createElement.apply(null,r)}return a.createElement.apply(null,n)}m.displayName="MDXCreateElement"},8215:function(e,t,n){var a=n(7294);t.Z=function(e){var t=e.children,n=e.hidden,i=e.className;return a.createElement("div",{role:"tabpanel",hidden:n,className:i},t)}},6396:function(e,t,n){n.d(t,{Z:function(){return m}});var a=n(7462),i=n(7294),o=n(2389),r=n(9443);var l=function(){var e=(0,i.useContext)(r.Z);if(null==e)throw new Error('"useUserPreferencesContext" is used outside of "Layout" component.');return e},s=n(9521),u=n(6010),p="tabItem_1uMI";function d(e){var t,n,a,o=e.lazy,r=e.block,d=e.defaultValue,m=e.values,c=e.groupId,h=e.className,f=i.Children.map(e.children,(function(e){if((0,i.isValidElement)(e)&&void 0!==e.props.value)return e;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})),g=null!=m?m:f.map((function(e){var t=e.props;return{value:t.value,label:t.label}})),k=(0,s.lx)(g,(function(e,t){return e.value===t.value}));if(k.length>0)throw new Error('Docusaurus error: Duplicate values "'+k.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.');var y=null===d?d:null!=(t=null!=d?d:null==(n=f.find((function(e){return e.props.default})))?void 0:n.props.value)?t:null==(a=f[0])?void 0:a.props.value;if(null!==y&&!g.some((function(e){return e.value===y})))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+y+'" but none of its children has the corresponding value. Available values are: '+g.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");var v=l(),w=v.tabGroupChoices,b=v.setTabGroupChoices,N=(0,i.useState)(y),C=N[0],T=N[1],j=[],x=(0,s.o5)().blockElementScrollPositionUntilNextRender;if(null!=c){var I=w[c];null!=I&&I!==C&&g.some((function(e){return e.value===I}))&&T(I)}var E=function(e){var t=e.currentTarget,n=j.indexOf(t),a=g[n].value;a!==C&&(x(t),T(a),null!=c&&b(c,a))},P=function(e){var t,n=null;switch(e.key){case"ArrowRight":var a=j.indexOf(e.currentTarget)+1;n=j[a]||j[0];break;case"ArrowLeft":var i=j.indexOf(e.currentTarget)-1;n=j[i]||j[j.length-1]}null==(t=n)||t.focus()};return i.createElement("div",{className:"tabs-container"},i.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,u.Z)("tabs",{"tabs--block":r},h)},g.map((function(e){var t=e.value,n=e.label;return i.createElement("li",{role:"tab",tabIndex:C===t?0:-1,"aria-selected":C===t,className:(0,u.Z)("tabs__item",p,{"tabs__item--active":C===t}),key:t,ref:function(e){return j.push(e)},onKeyDown:P,onFocus:E,onClick:E},null!=n?n:t)}))),o?(0,i.cloneElement)(f.filter((function(e){return e.props.value===C}))[0],{className:"margin-vert--md"}):i.createElement("div",{className:"margin-vert--md"},f.map((function(e,t){return(0,i.cloneElement)(e,{key:t,hidden:e.props.value!==C})}))))}function m(e){var t=(0,o.Z)();return i.createElement(d,(0,a.Z)({key:String(t)},e))}},9443:function(e,t,n){var a=(0,n(7294).createContext)(void 0);t.Z=a},1446:function(e,t,n){n.r(t),n.d(t,{frontMatter:function(){return u},contentTitle:function(){return p},metadata:function(){return d},toc:function(){return m},default:function(){return h}});var a=n(7462),i=n(3366),o=(n(7294),n(3905)),r=n(6396),l=n(8215),s=["components"],u={id:"juliaprimer",title:"A quick Julia primer",sidebar_label:"A quick Julia primer"},p=void 0,d={unversionedId:"gettingstarted/juliaprimer",id:"gettingstarted/juliaprimer",isDocsHomePage:!1,title:"A quick Julia primer",description:"For getting the most out of this documentation",source:"@site/docs/gettingstarted/juliaprimer.md",sourceDirName:"gettingstarted",slug:"/gettingstarted/juliaprimer",permalink:"/PopGen.jl/docs/gettingstarted/juliaprimer",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/gettingstarted/juliaprimer.md",tags:[],version:"current",lastUpdatedAt:1635859759,formattedLastUpdatedAt:"11/2/2021",frontMatter:{id:"juliaprimer",title:"A quick Julia primer",sidebar_label:"A quick Julia primer"},sidebar:"docs",previous:{title:"Installation",permalink:"/PopGen.jl/docs/"},next:{title:"PopGen.jl tips",permalink:"/PopGen.jl/docs/gettingstarted/tips"}},m=[{value:"Using Julia",id:"using-julia",children:[],level:2},{value:"First-time Performance",id:"first-time-performance",children:[],level:2},{value:"Semicolons",id:"semicolons",children:[{value:"At the end of a command",id:"at-the-end-of-a-command",children:[],level:3},{value:"In between assignment commands",id:"in-between-assignment-commands",children:[],level:3}],level:2},{value:"Help mode",id:"help-mode",children:[],level:2},{value:"Type information",id:"type-information",children:[{value:"Type Unions",id:"type-unions",children:[],level:3},{value:"Subtypes",id:"subtypes",children:[{value:"where T",id:"where-t",children:[],level:4}],level:3}],level:2},{value:"Functions vs. Methods",id:"functions-vs-methods",children:[{value:"ERROR: MethodError: no method matching",id:"error-methoderror-no-method-matching",children:[],level:3}],level:2},{value:"Functions with and without keywords",id:"functions-with-and-without-keywords",children:[],level:2}],c={toc:m};function h(e){var t=e.components,n=(0,i.Z)(e,s);return(0,o.kt)("wrapper",(0,a.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,o.kt)("p",null,(0,o.kt)("em",{parentName:"p"},"For getting the most out of this documentation")),(0,o.kt)("p",null,"There is nothing inherently special about this documentation relative to other documentation, other than we really ",(0,o.kt)("em",{parentName:"p"},"really")," want you to get the most out of what's written here. This means that we need to embrace the fact that both novice and experienced Julia users may be reading these docs and using this package. So let's cover some Julia basics that will really help in navigating this package before we even get into the complicated genetics stuff. This primer is by no means \"everything you need to get started in Julia\", and is a poor substitute for actually learning the language. In general, we recommend ",(0,o.kt)("a",{parentName:"p",href:"https://benlauwens.github.io/ThinkJulia.jl/latest/book.html"},"Think Julia: How to Think Like a Data Scientist")," by Ben Lauwens to establish some solid Julia foundations. It's free online! Also, the Julia language maintains ",(0,o.kt)("a",{parentName:"p",href:"https://docs.julialang.org/en/v1/"},"its own great documentation")," that we rely on quite heavily for development."),(0,o.kt)("h2",{id:"using-julia"},"Using Julia"),(0,o.kt)("p",null,"Everyone has their own particular workflows, and if you're new to Julia, you might not have established one yet. Julia can be used rather comfortably using its built-in interpreter. For an RStudio-like experience, we recommend using the ",(0,o.kt)("a",{parentName:"p",href:"https://www.julia-vscode.org"},"VScode Julia extension"),", but you can also use ",(0,o.kt)("a",{parentName:"p",href:"https://junolab.org/"},"Atom")," add-on). If you're already a fan of Jupyter notebooks (or ",(0,o.kt)("a",{parentName:"p",href:"https://nteract.io/"},(0,o.kt)("strong",{parentName:"a"},"nteract")),"), then all you need is to install the ",(0,o.kt)("inlineCode",{parentName:"p"},"IJulia")," package in Julia and you have full Jupyter support for Julia! You can also use the new reactive notebooks provided by ",(0,o.kt)("a",{parentName:"p",href:"https://github.com/fonsp/Pluto.jl"},"Pluto.jl"),"."),(0,o.kt)("div",{className:"admonition admonition-note alert alert--secondary"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.3 5.69a.942.942 0 0 1-.28-.7c0-.28.09-.52.28-.7.19-.18.42-.28.7-.28.28 0 .52.09.7.28.18.19.28.42.28.7 0 .28-.09.52-.28.7a1 1 0 0 1-.7.3c-.28 0-.52-.11-.7-.3zM8 7.99c-.02-.25-.11-.48-.31-.69-.2-.19-.42-.3-.69-.31H6c-.27.02-.48.13-.69.31-.2.2-.3.44-.31.69h1v3c.02.27.11.5.31.69.2.2.42.31.69.31h1c.27 0 .48-.11.69-.31.2-.19.3-.42.31-.69H8V7.98v.01zM7 2.3c-3.14 0-5.7 2.54-5.7 5.68 0 3.14 2.56 5.7 5.7 5.7s5.7-2.55 5.7-5.7c0-3.15-2.56-5.69-5.7-5.69v.01zM7 .98c3.86 0 7 3.14 7 7s-3.14 7-7 7-7-3.12-7-7 3.14-7 7-7z"}))),"Trivia")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},'If you didn\'t already know,  the name "Jupyter" is actually a concatenation of ',(0,o.kt)("strong",{parentName:"p"},"Ju")," (julia) ",(0,o.kt)("strong",{parentName:"p"},"Pyt")," (python) and ",(0,o.kt)("strong",{parentName:"p"},"eR")," (R). \ud83e\udd2f"))),(0,o.kt)("h2",{id:"first-time-performance"},"First-time Performance"),(0,o.kt)("p",null,"If you're migrating to Julia from Python or R (or Matlab, etc.), you'll think Julia is slow and laggy because loading packages and running stuff has a noticeable wait time (10-40sec). It's worth mentioning that this lag is \"compilation overhead\". What this means is, Julia tries to pre-compile as much code as possible (into optimized machine code) when running something or loading a package. This lag exists ",(0,o.kt)("strong",{parentName:"p"},"only the first time")," you run something. Every subsequent run of a function, even with different parameters, will be ",(0,o.kt)("strong",{parentName:"p"},"substantially")," faster, and in most cases instant. If you want to test this yourself, try to run a line of code twice with ",(0,o.kt)("inlineCode",{parentName:"p"},"@time")," before the function and compare the results. Here's an example:"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre"},"julia> @time using PopGen\n 11.870878 seconds (27.83 M allocations: 1.389 GiB, 5.26% gc time)\n\njulia> @time using PopGen\n 0.000272 seconds (406 allocations: 21.375 KiB)\n")),(0,o.kt)("h2",{id:"semicolons"},"Semicolons"),(0,o.kt)("p",null,"Semicolons will come up a lot in Julia, probably more than you would expect if you are migrating from another language.  They mean different things depending on where they are."),(0,o.kt)("h3",{id:"at-the-end-of-a-command"},"At the end of a command"),(0,o.kt)("p",null,' When you see a semicolon after invoking a function, what that means is "don\'t show me the output".'),(0,o.kt)("p",null,"Example: "),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"julia> x = 2 + 2\n4\n\njulia> x = 10 + 2 ;\n\njulia> x\n12\n")),(0,o.kt)("p",null,"Julia will still process the command and assign ",(0,o.kt)("inlineCode",{parentName:"p"},"10 + 2")," to ",(0,o.kt)("inlineCode",{parentName:"p"},"x"),", but it won't show you the output. We sometimes include a semicolon after commands in these docs to mimic what the REPL output would look like without spitting back large volumes of text. ",(0,o.kt)("strong",{parentName:"p"},"These semicolons are optional")," "),(0,o.kt)("h3",{id:"in-between-assignment-commands"},"In between assignment commands"),(0,o.kt)("p",null,"If you see a semicolon in between two variable assignments or commands, like so:"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"julia> x = [1,2] ; y = [3,4]\n")),(0,o.kt)("p",null,"that's a Julia short-hand for making two short lines of code appear on a single line. It's the equivalent of doing:"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"julia> x = [1,2]\njulia> y = [3,4]\n")),(0,o.kt)("p",null,"We sometimes choose this writing format for very quick  and small assignments hoping to save some visual space. Use whichever method is most comfortable and sensible for you!"),(0,o.kt)("h2",{id:"help-mode"},"Help mode"),(0,o.kt)("p",null,"To enter ",(0,o.kt)("inlineCode",{parentName:"p"},"help")," mode in the REPL, simply press the question mark key ",(0,o.kt)("inlineCode",{parentName:"p"},"?")," (shift + key) and you will notice a different prompt ",(0,o.kt)("inlineCode",{parentName:"p"},"help?>")," for you to type in a function."),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"help?>population\n")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre"},"search: population populations population! populations!\n\n  population(data::PopData; listall::Bool = false)\n\nView unique population ID's and their counts in a `PopData`.\n\n\n- `listall = true` displays all samples and their `population` instead (default = `false`)\n")),(0,o.kt)("h2",{id:"type-information"},"Type information"),(0,o.kt)("p",null,"Julia encourages strong typing of variables, and the functions in ",(0,o.kt)("inlineCode",{parentName:"p"},"PopGen")," are no exception to this. However, to reduce the barrier of entry required to understand this documentation and the subsequent package, we have chosen to omit some of the ",(0,o.kt)("inlineCode",{parentName:"p"},"type")," information from functions to reduce visual clutter for newer users. As experienced users already know, if you would like to see the explicit type information, you can look at the code on github, invoke the ",(0,o.kt)("inlineCode",{parentName:"p"},"help")," system in the REPL (above), or search for a function in the Documentation pane in Juno. "),(0,o.kt)("p",null,"You'll notice types follow a specific format, which is ",(0,o.kt)("inlineCode",{parentName:"p"},"object::type"),". This format is a type declaration, so in the function ",(0,o.kt)("inlineCode",{parentName:"p"},"population"),", which looks like: "),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"population(data::PopData; counts::Bool = false)\n")),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"data")," is a positional argument of type ",(0,o.kt)("inlineCode",{parentName:"li"},"PopData")," "),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"counts")," is a keyword argument (because of the semicolon, see below) of type ",(0,o.kt)("inlineCode",{parentName:"li"},"Bool")," (Boolean) meaning it only takes ",(0,o.kt)("inlineCode",{parentName:"li"},"true")," or ",(0,o.kt)("inlineCode",{parentName:"li"},"false")," without quotes, and the default value is set to ",(0,o.kt)("inlineCode",{parentName:"li"},"false"))),(0,o.kt)("h3",{id:"type-unions"},"Type Unions"),(0,o.kt)("p",null,"You might see the type ",(0,o.kt)("inlineCode",{parentName:"p"},"Union")," appear occasionally throughout this documentation, and you can consider it a list of allowable types. For example, if something was of type ",(0,o.kt)("inlineCode",{parentName:"p"},"::Union{String,Integer}"),", that means that ",(0,o.kt)("strong",{parentName:"p"},"either")," a ",(0,o.kt)("inlineCode",{parentName:"p"},"String")," or ",(0,o.kt)("inlineCode",{parentName:"p"},"Integer")," works. "),(0,o.kt)("h3",{id:"subtypes"},"Subtypes"),(0,o.kt)("p",null,"The julia language is abound with types (and you can create your own!), and has a hierarchical system of supertypes and subtypes. As you can probably guess, a supertype can contain multiple subtypes, such as ",(0,o.kt)("inlineCode",{parentName:"p"},"Signed")," being a supertype of (among other things) ",(0,o.kt)("inlineCode",{parentName:"p"},"Int64, Int32, Int16, Int8"),". All vectors are subtypes of ",(0,o.kt)("inlineCode",{parentName:"p"},"AbstractVector"),".  If you want to try it yourself, use the ",(0,o.kt)("inlineCode",{parentName:"p"},"supertype()")," command on your favorite Type, like ",(0,o.kt)("inlineCode",{parentName:"p"},"supertype(Float32)"),". You will occasionally see ",(0,o.kt)("inlineCode",{parentName:"p"},"<:")," instead of ",(0,o.kt)("inlineCode",{parentName:"p"},"::"),', which means "is a subtype of". This is used for condtional evaluation, like ',(0,o.kt)("inlineCode",{parentName:"p"},"typeof(something) <: Signed"),", and in some function methods like ",(0,o.kt)("inlineCode",{parentName:"p"},"function(var1::T) where T <: Supertype"),", which leads us to:"),(0,o.kt)("h4",{id:"where-t"},"where T"),(0,o.kt)("p",null,"This looks weird at first,  but it's actually very simple. When we do method definitions, we can define methods with strict types, like ",(0,o.kt)("inlineCode",{parentName:"p"},"funct(data:PopData, arg1::Int8)"),", or we can generalize it with ",(0,o.kt)("inlineCode",{parentName:"p"},"where T"),", which looks like :"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"function funct1(data::PopData, thing1::T) where T\n")),(0,o.kt)("p",null,"This will auto-create a method for any possible Type for ",(0,o.kt)("inlineCode",{parentName:"p"},"thing1"),". That's really convenvient, but sometimes it's problematic, as incorrect input can lead to obscure errors (e.g. multiplying integers with strings?!). Instead, you can constrain the types for ",(0,o.kt)("inlineCode",{parentName:"p"},"T")," like this:"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"function funct2(data::PopData, thing1::T) where T <: Signed\n")),(0,o.kt)("p",null,"With the constraint above, it will generate methods for all cases where ",(0,o.kt)("inlineCode",{parentName:"p"},"thing1")," is a subtype of ",(0,o.kt)("inlineCode",{parentName:"p"},"Signed"),", which includes all the numerical Types (integers, Floats, etc.). This will make sure that the function will behave correctly for a range of input types."),(0,o.kt)("p",null,"You can also use this type of notation to clean up a method definition where multiple arguments have the same Type specification:"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"function funct3(data::PopData, thing1::T, thing2::T, thing3::T) where T<:AbstractFloat\n")),(0,o.kt)("p",null,"So, instead of writing ",(0,o.kt)("inlineCode",{parentName:"p"},"thing1::Float64, thing2::Float64, thing3::Float64"),", we just use ",(0,o.kt)("inlineCode",{parentName:"p"},"T")," as a placeholder and assign it as a subtype of something at the end.  It ends up being pretty handy!"),(0,o.kt)("h2",{id:"functions-vs-methods"},"Functions vs. Methods"),(0,o.kt)("p",null,'As part of Julia\'s type-safe paradigm and multiple dispatch (see "ERROR: MethodError: no method matching" below), type specifications in functions often reduce runtime of functions, but also establish function identity. Multiple dispatch refers to several different functions having the same name, but employing different ',(0,o.kt)("em",{parentName:"p"},"methods")," depending on the input. In Julia, it's easier to write a single function with multiple type-safe methods, rather than one mega-function that accepts any type and have a bunch of ",(0,o.kt)("inlineCode",{parentName:"p"},"if")," statements that determines what the program does depending on the input. "),(0,o.kt)("div",{className:"admonition admonition-info alert alert--info"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"}))),"Best Practice")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},"As a rule of thumb, ",(0,o.kt)("inlineCode",{parentName:"p"},"for"),"  loops with ",(0,o.kt)("inlineCode",{parentName:"p"},"if")," conditions in them slow down the compiler, so best-practice often encourages us to write type-specific methods."))),(0,o.kt)("p",null,"In practice, this looks like:"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre"},'# combine two numbers\njulia> function add(x::Integer, y::Integer)\n           x+y\n       end\nadd (generic function with 1 method)\n\n# combine two strings\njulia> function add(x::String, y::String)\n           x*y\n       end\nadd (generic function with 2 methods)\n                \njulia> add(1,2)\n3\n\njulia> add("water", "melon")\n"watermelon"\n')),(0,o.kt)("p",null,"Multiple dispatch therefor leads to a unique type of possible error: the ",(0,o.kt)("inlineCode",{parentName:"p"},"MethodError")),(0,o.kt)("h3",{id:"error-methoderror-no-method-matching"},"ERROR: MethodError: no method matching"),(0,o.kt)("p",null,"Using the function ",(0,o.kt)("inlineCode",{parentName:"p"},"add")," from the example above, let's have a look at what happens when we try to ",(0,o.kt)("inlineCode",{parentName:"p"},"add")," an ",(0,o.kt)("inlineCode",{parentName:"p"},"Integer")," with a ",(0,o.kt)("inlineCode",{parentName:"p"},"String"),":"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'julia> add(1,"melon")\nERROR: MethodError: no method matching add(::Int64, ::String)\nClosest candidates are:\n  add(::String, ::String) at none:2\n  add(::Integer, ::Integer) at none:2\nStacktrace:\n [1] top-level scope at none:0\n')),(0,o.kt)("p",null,'This error is telling us "there is no such function called ',(0,o.kt)("inlineCode",{parentName:"p"},"add"),", who's inputs are an ",(0,o.kt)("inlineCode",{parentName:"p"},"Integer")," followed by a ",(0,o.kt)("inlineCode",{parentName:"p"},"String"),'". But, it does offer us some alternatives, like the two ',(0,o.kt)("inlineCode",{parentName:"p"},"add")," functions we created earlier."),(0,o.kt)("p",null,"The functions within ",(0,o.kt)("inlineCode",{parentName:"p"},"PopGen")," are almost always explicitly typed, so if you are getting the ",(0,o.kt)("inlineCode",{parentName:"p"},"MethodError: no method matching"),' error, then you are inputting the incorrect types into the function, or perhaps your inputs for the arguments are in the wrong order (see "Functions with and without keywords" below).'),(0,o.kt)("p",null,"Sometimes you might include an argument with a keyword when there isn't one, or include an argument without a keyword when there needs to be one (honestly, we make that mistake too and we ",(0,o.kt)("em",{parentName:"p"},"wrote")," this stuff). To help minimize those mistakes, please read below about which arguments have keywords and which don't."),(0,o.kt)("div",{className:"admonition admonition-note alert alert--secondary"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.3 5.69a.942.942 0 0 1-.28-.7c0-.28.09-.52.28-.7.19-.18.42-.28.7-.28.28 0 .52.09.7.28.18.19.28.42.28.7 0 .28-.09.52-.28.7a1 1 0 0 1-.7.3c-.28 0-.52-.11-.7-.3zM8 7.99c-.02-.25-.11-.48-.31-.69-.2-.19-.42-.3-.69-.31H6c-.27.02-.48.13-.69.31-.2.2-.3.44-.31.69h1v3c.02.27.11.5.31.69.2.2.42.31.69.31h1c.27 0 .48-.11.69-.31.2-.19.3-.42.31-.69H8V7.98v.01zM7 2.3c-3.14 0-5.7 2.54-5.7 5.68 0 3.14 2.56 5.7 5.7 5.7s5.7-2.55 5.7-5.7c0-3.15-2.56-5.69-5.7-5.69v.01zM7 .98c3.86 0 7 3.14 7 7s-3.14 7-7 7-7-3.12-7-7 3.14-7 7-7z"}))),"MethodErrors")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},"MethodError's can definitely get annoying, but they are more commonly the result of incorrect inputs versus being bugs. If you double-checked your inputs and things still don't work, please submit an issue. Thanks!"))),(0,o.kt)("h2",{id:"functions-with-and-without-keywords"},"Functions with and without keywords"),(0,o.kt)("p",null,"Let's talk about semicolons some more."),(0,o.kt)("div",{className:"admonition admonition-info alert alert--info"},(0,o.kt)("div",{parentName:"div",className:"admonition-heading"},(0,o.kt)("h5",{parentName:"div"},(0,o.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,o.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,o.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"}))),"TL;DR")),(0,o.kt)("div",{parentName:"div",className:"admonition-content"},(0,o.kt)("p",{parentName:"div"},"Reading these docs, pay attention to semicolons in the function argument lists."),(0,o.kt)("ul",{parentName:"div"},(0,o.kt)("li",{parentName:"ul"},"arguments before a semicolon have no keyword and follow an explicit order"),(0,o.kt)("li",{parentName:"ul"},"arguments after a semicolon have a keyword ",(0,o.kt)("inlineCode",{parentName:"li"},"argument = value")," and their order doesn't matter"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"MethodError: no methods matching")," is often a user error and not a bug, but if it is, please open an issue!")))),(0,o.kt)("p",null,'Broadly speaking, there are two types of function declarations in Julia: ones with keywords and ones without keywords. The term "keywords" refers to an input argument that has the format ',(0,o.kt)("inlineCode",{parentName:"p"},"argument = value"),". This format is present in many of the functions in this and other packages, however there are some specifics to understand when functions use keywords and when they don't. "),(0,o.kt)(r.Z,{block:!0,defaultValue:"1",values:[{label:"1. No semicolon in argument list",value:"1"},{label:"2. Semicolon in argument list",value:"2"}],mdxType:"Tabs"},(0,o.kt)(l.Z,{value:"1",mdxType:"TabItem"},(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"function function_name(var1::type, var2::type, var3::type)\n    do stuff with vars\nend\n")),(0,o.kt)("p",null,"If a function is declared with only commas in the argument list, like shown above, then the arguments to that function ",(0,o.kt)("strong",{parentName:"p"},"must")," have no keywords and follow the exact order they appear in. If the generic example above had the typing:"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"function function_name(var1::String, var2::Float64, var3::Array{String,1})\n    do stuff with vars\nend\n")),(0,o.kt)("p",null,"then the only acceptable way to run this function without getting a ",(0,o.kt)("inlineCode",{parentName:"p"},"MethodError")," would be with arguments in the order of ",(0,o.kt)("inlineCode",{parentName:"p"},"function_name(String, Float64, Array{String,1})"),". Even if some of the arguments have a default values, like ",(0,o.kt)("inlineCode",{parentName:"p"},"var2::Float64 = 6.66"),", the order of arguments/types has to be respected as declared.")),(0,o.kt)(l.Z,{value:"2",mdxType:"TabItem"},(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"function function_name(var1::type; var2::type, var3::type)\n    do stuff with vars\nend\n")),(0,o.kt)("p",null,"In this format, everything that comes ",(0,o.kt)("strong",{parentName:"p"},"before")," the semicolon follows the strict rules from ",(0,o.kt)("strong",{parentName:"p"},"Format 1"),", and everything that comes ",(0,o.kt)("strong",{parentName:"p"},"after")," the semicolon is a keyword argument. Keyword arguments have the flexibility to not require any particular input order. However, you ",(0,o.kt)("strong",{parentName:"p"},"must")," use the keywords to declare those arguments, or you will receive another ",(0,o.kt)("inlineCode",{parentName:"p"},"MethodError: no method matching"),", which is, as we've mentioned, annoying. "))))}h.isMDXComponent=!0}}]);