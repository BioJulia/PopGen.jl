"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[2221],{3905:function(e,n,t){t.d(n,{Zo:function(){return c},kt:function(){return m}});var r=t(7294);function a(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function o(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);n&&(r=r.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,r)}return t}function i(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?o(Object(t),!0).forEach((function(n){a(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):o(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,r,a=function(e,n){if(null==e)return{};var t,r,a={},o=Object.keys(e);for(r=0;r<o.length;r++)t=o[r],n.indexOf(t)>=0||(a[t]=e[t]);return a}(e,n);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(r=0;r<o.length;r++)t=o[r],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(a[t]=e[t])}return a}var p=r.createContext({}),s=function(e){var n=r.useContext(p),t=n;return e&&(t="function"==typeof e?e(n):i(i({},n),e)),t},c=function(e){var n=s(e.components);return r.createElement(p.Provider,{value:n},e.children)},d={inlineCode:"code",wrapper:function(e){var n=e.children;return r.createElement(r.Fragment,{},n)}},u=r.forwardRef((function(e,n){var t=e.components,a=e.mdxType,o=e.originalType,p=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),u=s(t),m=a,g=u["".concat(p,".").concat(m)]||u[m]||d[m]||o;return t?r.createElement(g,i(i({ref:n},c),{},{components:t})):r.createElement(g,i({ref:n},c))}));function m(e,n){var t=arguments,a=n&&n.mdxType;if("string"==typeof e||a){var o=t.length,i=new Array(o);i[0]=u;var l={};for(var p in n)hasOwnProperty.call(n,p)&&(l[p]=n[p]);l.originalType=e,l.mdxType="string"==typeof e?e:a,i[1]=l;for(var s=2;s<o;s++)i[s]=t[s];return r.createElement.apply(null,i)}return r.createElement.apply(null,t)}u.displayName="MDXCreateElement"},4230:function(e,n,t){t.r(n),t.d(n,{frontMatter:function(){return l},contentTitle:function(){return p},metadata:function(){return s},toc:function(){return c},default:function(){return u}});var r=t(7462),a=t(3366),o=(t(7294),t(3905)),i=["components"],l={id:"popgensims_cross",title:"Cross.jl",sidebar_label:"Cross.jl"},p=void 0,s={unversionedId:"api/PopGenSims/popgensims_cross",id:"api/PopGenSims/popgensims_cross",isDocsHomePage:!1,title:"Cross.jl",description:"PopGenSims.jl/src/Cross.jl",source:"@site/docs/api/PopGenSims/Cross.md",sourceDirName:"api/PopGenSims",slug:"/api/PopGenSims/popgensims_cross",permalink:"/PopGen.jl/docs/api/PopGenSims/popgensims_cross",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGenSims/Cross.md",tags:[],version:"current",lastUpdatedAt:1635528951,formattedLastUpdatedAt:"10/29/2021",frontMatter:{id:"popgensims_cross",title:"Cross.jl",sidebar_label:"Cross.jl"},sidebar:"docs",previous:{title:"VariantCall.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/variantcall"},next:{title:"Samples.jl",permalink:"/PopGen.jl/docs/api/PopGenSims/popgensims_samples"}},c=[{value:"PopGenSims.jl/src/Cross.jl",id:"popgensimsjlsrccrossjl",children:[{value:"\u2757sample_genotype",id:"sample_genotype",children:[],level:3},{value:"\u2757haploid_cross!`",id:"haploid_cross",children:[],level:3},{value:"\u2757polyploid_cross!",id:"polyploid_cross",children:[],level:3},{value:"\u26ab cross",id:"-cross",children:[{value:"Keyword Arguments",id:"keyword-arguments",children:[],level:4},{value:"Keyword Arguments",id:"keyword-arguments-1",children:[],level:4}],level:3}],level:2}],d={toc:c};function u(e){var n=e.components,t=(0,a.Z)(e,i);return(0,o.kt)("wrapper",(0,r.Z)({},d,t,{components:n,mdxType:"MDXLayout"}),(0,o.kt)("h2",{id:"popgensimsjlsrccrossjl"},"PopGenSims.jl/src/Cross.jl"),(0,o.kt)("p",null,"\u2757 => not exported |\n\u26ab => exported by PopGenSims.jl"),(0,o.kt)("h3",{id:"sample_genotype"},"\u2757sample_genotype"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"sample_genotype(geno::T, n_alleles::Int) where T<:Genotype\n")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"sample_genotype(geno::Missing, n_alleles::Int)\n")),(0,o.kt)("hr",null),(0,o.kt)("h3",{id:"haploid_cross"},"\u2757haploid_cross!`"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"haploid_cross!(data::DataFrame, p1::T, p2::T; n::Int) where T <: GenoArray\n")),(0,o.kt)("hr",null),(0,o.kt)("h3",{id:"polyploid_cross"},"\u2757polyploid_cross!"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},"polyploid_cross!(data::DataFrame, p1::T, p2::T; n::Int, ploidy::Int) where T <: GenoArray\n")),(0,o.kt)("hr",null),(0,o.kt)("h3",{id:"-cross"},"\u26ab cross"),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'cross(data::PopData, parent1::String, parent2::String; n::Int = 100, generation::String = "F1")\n')),(0,o.kt)("p",null,"Simulate a breeding cross between individuals ",(0,o.kt)("inlineCode",{parentName:"p"},"parent1")," and ",(0,o.kt)("inlineCode",{parentName:"p"},"parent2")," from a ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," object.\nReturns PopData consisting of ",(0,o.kt)("inlineCode",{parentName:"p"},"n")," offspring resulting from the cross."),(0,o.kt)("h4",{id:"keyword-arguments"},"Keyword Arguments"),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"n")," : Integer of number of offspring to generate (default: ",(0,o.kt)("inlineCode",{parentName:"li"},"100"),")"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"generation")," : A string to assign ",(0,o.kt)("inlineCode",{parentName:"li"},"population")," identity to the offspring (default: ",(0,o.kt)("inlineCode",{parentName:"li"},'"F1"'),")")),(0,o.kt)("pre",null,(0,o.kt)("code",{parentName:"pre",className:"language-julia"},'cross(parent_1::Pair, parent_2::Pair, n::Int = 100, generation::String = "F1")\n')),(0,o.kt)("p",null,"Simulate a breeding cross between individuals ",(0,o.kt)("inlineCode",{parentName:"p"},"parent")," and ",(0,o.kt)("inlineCode",{parentName:"p"},"parent2")," from two different ",(0,o.kt)("inlineCode",{parentName:"p"},"PopData")," objects.\nReturns PopData consisting of ",(0,o.kt)("inlineCode",{parentName:"p"},"n")," offspring resulting from the cross. ",(0,o.kt)("inlineCode",{parentName:"p"},"parent_1_data")," and ",(0,o.kt)("inlineCode",{parentName:"p"},"parent_2_data"),"\nare positional arguments, therefore they must be written without keywords and in the order of parents 1, parent 2. "),(0,o.kt)("h4",{id:"keyword-arguments-1"},"Keyword Arguments"),(0,o.kt)("ul",null,(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"parent_1")," : Pair of ",(0,o.kt)("inlineCode",{parentName:"li"},'PopData => "Parent1Name"')),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"parent_2")," : Pair of ",(0,o.kt)("inlineCode",{parentName:"li"},'PopData => "Parent1Name"')),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"n")," : Integer of number of offspring to generate (default: ",(0,o.kt)("inlineCode",{parentName:"li"},"100"),")"),(0,o.kt)("li",{parentName:"ul"},(0,o.kt)("inlineCode",{parentName:"li"},"generation")," : A string to assign ",(0,o.kt)("inlineCode",{parentName:"li"},"population")," identity to the offspring (default: ",(0,o.kt)("inlineCode",{parentName:"li"},'"F1"'),")")))}u.isMDXComponent=!0}}]);