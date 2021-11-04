"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[4467],{3905:function(e,t,n){n.d(t,{Zo:function(){return m},kt:function(){return f}});var r=n(7294);function a(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function i(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);t&&(r=r.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,r)}return n}function l(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?i(Object(n),!0).forEach((function(t){a(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):i(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function o(e,t){if(null==e)return{};var n,r,a=function(e,t){if(null==e)return{};var n,r,a={},i=Object.keys(e);for(r=0;r<i.length;r++)n=i[r],t.indexOf(n)>=0||(a[n]=e[n]);return a}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(r=0;r<i.length;r++)n=i[r],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(a[n]=e[n])}return a}var p=r.createContext({}),d=function(e){var t=r.useContext(p),n=t;return e&&(n="function"==typeof e?e(t):l(l({},t),e)),n},m=function(e){var t=d(e.components);return r.createElement(p.Provider,{value:t},e.children)},c={inlineCode:"code",wrapper:function(e){var t=e.children;return r.createElement(r.Fragment,{},t)}},s=r.forwardRef((function(e,t){var n=e.components,a=e.mdxType,i=e.originalType,p=e.parentName,m=o(e,["components","mdxType","originalType","parentName"]),s=d(n),f=a,k=s["".concat(p,".").concat(f)]||s[f]||c[f]||i;return n?r.createElement(k,l(l({ref:t},m),{},{components:n})):r.createElement(k,l({ref:t},m))}));function f(e,t){var n=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var i=n.length,l=new Array(i);l[0]=s;var o={};for(var p in t)hasOwnProperty.call(t,p)&&(o[p]=t[p]);o.originalType=e,o.mdxType="string"==typeof e?e:a,l[1]=o;for(var d=2;d<i;d++)l[d]=n[d];return r.createElement.apply(null,l)}return r.createElement.apply(null,n)}s.displayName="MDXCreateElement"},603:function(e,t,n){n.r(t),n.d(t,{frontMatter:function(){return o},contentTitle:function(){return p},metadata:function(){return d},toc:function(){return m},default:function(){return s}});var r=n(7462),a=n(3366),i=(n(7294),n(3905)),l=["components"],o={id:"read",title:"ReadWrite.jl",sidebar_label:"ReadWrite.jl"},p=void 0,d={unversionedId:"api/PopGenCore/read",id:"api/PopGenCore/read",isDocsHomePage:!1,title:"ReadWrite.jl",description:"PopGenCore.jl/src/io/ReadWrite.jl",source:"@site/docs/api/PopGenCore/ReadWrite.md",sourceDirName:"api/PopGenCore",slug:"/api/PopGenCore/read",permalink:"/PopGen.jl/docs/api/PopGenCore/read",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/api/PopGenCore/ReadWrite.md",tags:[],version:"current",lastUpdatedAt:1636029729,formattedLastUpdatedAt:"11/4/2021",frontMatter:{id:"read",title:"ReadWrite.jl",sidebar_label:"ReadWrite.jl"},sidebar:"docs",previous:{title:"PopDataWrappers.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/popdatawrappers"},next:{title:"Structure.jl",permalink:"/PopGen.jl/docs/api/PopGenCore/structure"}},m=[{value:"PopGenCore.jl/src/io/ReadWrite.jl",id:"popgencorejlsrcioreadwritejl",children:[{value:"\ud83d\udce6 read",id:"-read",children:[],level:3},{value:"\ud83d\udce6 write",id:"-write",children:[],level:3}],level:2}],c={toc:m};function s(e){var t=e.components,n=(0,a.Z)(e,l);return(0,i.kt)("wrapper",(0,r.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("h2",{id:"popgencorejlsrcioreadwritejl"},"PopGenCore.jl/src/io/ReadWrite.jl"),(0,i.kt)("p",null,"\ud83d\udce6  => not exported |\n\ud83d\udfea => exported by PopGenCore.jl |\n\ud83d\udd35 => exported by PopGen.jl"),(0,i.kt)("h3",{id:"-read"},"\ud83d\udce6 read"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"PopGen.read(infile::String; kwargs...)\n")),(0,i.kt)("p",null,"Wraps the individual file importers to read a file in as a ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," object. File type is\ninferred from the file extension (case insensitive):"),(0,i.kt)("table",null,(0,i.kt)("thead",{parentName:"table"},(0,i.kt)("tr",{parentName:"thead"},(0,i.kt)("th",{parentName:"tr",align:"left"},"File Format"),(0,i.kt)("th",{parentName:"tr",align:"left"},"Extensions"),(0,i.kt)("th",{parentName:"tr",align:"left"},"Docstring"))),(0,i.kt)("tbody",{parentName:"table"},(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"delimited"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".csv"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".txt"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".tsv")),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},"?delimited"))),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"genepop"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".gen"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".genepop")),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},"?genepop"))),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"structure"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".str"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".structure")),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},"?structure"))),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"plink"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".bed"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".ped")),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},"?plink"))),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"variant call format (vcf)"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".vcf"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".vcf.gz")),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},"?vcf"))),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"variant call format (bcf)"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".bcf"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".bcf.gz")),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},"?bcf"))))),(0,i.kt)("p",null,"This function uses the same keyword arguments (and defaults) as the file importing\nfunctions it wraps; please see their respective docstrings in the Julia help console.\nfor specific usage details (e.g. ",(0,i.kt)("inlineCode",{parentName:"p"},"?genepop"),")."),(0,i.kt)("p",null,(0,i.kt)("strong",{parentName:"p"},"Examples")),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},'PopGen.read("cavernous_assfish.gen", digits = 3)\nPopGen.read("juglans_nigra.vcf")\n')),(0,i.kt)("hr",null),(0,i.kt)("h3",{id:"-write"},"\ud83d\udce6 write"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"PopGen.write(data::PopData, filename::String, kwargs...)\nPopGen.write(data::PopData; filename::String, kwargs...)\n")),(0,i.kt)("p",null,"Writes ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," to a specified file type inferred from the extension of ",(0,i.kt)("inlineCode",{parentName:"p"},"filename = ")," (case insensitive). Additional keyword\narguments ",(0,i.kt)("inlineCode",{parentName:"p"},"kwargs...")," are specific to the intended file type, and are listed in the docstrings of the specific\nfile writer with the format ",(0,i.kt)("inlineCode",{parentName:"p"},"?filetype"),". For example, to find the appropriate keywords for a conversion\nto Genepop format, call up the ",(0,i.kt)("inlineCode",{parentName:"p"},"?genepop")," docstring."),(0,i.kt)("table",null,(0,i.kt)("thead",{parentName:"table"},(0,i.kt)("tr",{parentName:"thead"},(0,i.kt)("th",{parentName:"tr",align:"left"},"File Format"),(0,i.kt)("th",{parentName:"tr",align:"left"},"Extensions"),(0,i.kt)("th",{parentName:"tr",align:"left"},"Docstring"))),(0,i.kt)("tbody",{parentName:"table"},(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"genepop"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".gen"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".genepop")),(0,i.kt)("td",{parentName:"tr",align:"left"},"?genepop")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"delimited"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".csv"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".txt"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".tsv")),(0,i.kt)("td",{parentName:"tr",align:"left"},"?delimited")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:"left"},"structure"),(0,i.kt)("td",{parentName:"tr",align:"left"},(0,i.kt)("inlineCode",{parentName:"td"},".str"),", ",(0,i.kt)("inlineCode",{parentName:"td"},".structure")),(0,i.kt)("td",{parentName:"tr",align:"left"},"?structure")))),(0,i.kt)("p",null,(0,i.kt)("strong",{parentName:"p"},"Example")),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},'cats = @nancycats;\nfewer_cats = omit(cats, name = samplenames(cats)[1:10]);\nPopGen.write(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")\n')))}s.isMDXComponent=!0}}]);