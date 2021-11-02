"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[5937],{3905:function(e,t,n){n.d(t,{Zo:function(){return c},kt:function(){return m}});var a=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function i(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?i(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):i(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,a,r=function(e,t){if(null==e)return{};var n,a,r={},i=Object.keys(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(a=0;a<i.length;a++)n=i[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var u=a.createContext({}),d=function(e){var t=a.useContext(u),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},c=function(e){var t=d(e.components);return a.createElement(u.Provider,{value:t},e.children)},s={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},p=a.forwardRef((function(e,t){var n=e.components,r=e.mdxType,i=e.originalType,u=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),p=d(n),m=r,f=p["".concat(u,".").concat(m)]||p[m]||s[m]||i;return n?a.createElement(f,o(o({ref:t},c),{},{components:n})):a.createElement(f,o({ref:t},c))}));function m(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var i=n.length,o=new Array(i);o[0]=p;var l={};for(var u in t)hasOwnProperty.call(t,u)&&(l[u]=t[u]);l.originalType=e,l.mdxType="string"==typeof e?e:r,o[1]=l;for(var d=2;d<i;d++)o[d]=n[d];return a.createElement.apply(null,o)}return a.createElement.apply(null,n)}p.displayName="MDXCreateElement"},8215:function(e,t,n){var a=n(7294);t.Z=function(e){var t=e.children,n=e.hidden,r=e.className;return a.createElement("div",{role:"tabpanel",hidden:n,className:r},t)}},5064:function(e,t,n){n.d(t,{Z:function(){return p}});var a=n(7462),r=n(7294),i=n(2389),o=n(9443);var l=function(){var e=(0,r.useContext)(o.Z);if(null==e)throw new Error('"useUserPreferencesContext" is used outside of "Layout" component.');return e},u=n(3039),d=n(6010),c="tabItem_1uMI";function s(e){var t,n,a,i=e.lazy,o=e.block,s=e.defaultValue,p=e.values,m=e.groupId,f=e.className,h=r.Children.map(e.children,(function(e){if((0,r.isValidElement)(e)&&void 0!==e.props.value)return e;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})),g=null!=p?p:h.map((function(e){var t=e.props;return{value:t.value,label:t.label}})),v=(0,u.lx)(g,(function(e,t){return e.value===t.value}));if(v.length>0)throw new Error('Docusaurus error: Duplicate values "'+v.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.');var k=null===s?s:null!=(t=null!=s?s:null==(n=h.find((function(e){return e.props.default})))?void 0:n.props.value)?t:null==(a=h[0])?void 0:a.props.value;if(null!==k&&!g.some((function(e){return e.value===k})))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+k+'" but none of its children has the corresponding value. Available values are: '+g.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");var w=l(),b=w.tabGroupChoices,y=w.setTabGroupChoices,N=(0,r.useState)(k),C=N[0],j=N[1],O=[],T=(0,u.o5)().blockElementScrollPositionUntilNextRender;if(null!=m){var P=b[m];null!=P&&P!==C&&g.some((function(e){return e.value===P}))&&j(P)}var D=function(e){var t=e.currentTarget,n=O.indexOf(t),a=g[n].value;a!==C&&(T(t),j(a),null!=m&&y(m,a))},E=function(e){var t,n=null;switch(e.key){case"ArrowRight":var a=O.indexOf(e.currentTarget)+1;n=O[a]||O[0];break;case"ArrowLeft":var r=O.indexOf(e.currentTarget)-1;n=O[r]||O[O.length-1]}null==(t=n)||t.focus()};return r.createElement("div",{className:"tabs-container"},r.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,d.Z)("tabs",{"tabs--block":o},f)},g.map((function(e){var t=e.value,n=e.label;return r.createElement("li",{role:"tab",tabIndex:C===t?0:-1,"aria-selected":C===t,className:(0,d.Z)("tabs__item",c,{"tabs__item--active":C===t}),key:t,ref:function(e){return O.push(e)},onKeyDown:E,onFocus:D,onClick:D},null!=n?n:t)}))),i?(0,r.cloneElement)(h.filter((function(e){return e.props.value===C}))[0],{className:"margin-vert--md"}):r.createElement("div",{className:"margin-vert--md"},h.map((function(e,t){return(0,r.cloneElement)(e,{key:t,hidden:e.props.value!==C})}))))}function p(e){var t=(0,i.Z)();return r.createElement(s,(0,a.Z)({key:String(t)},e))}},9443:function(e,t,n){var a=(0,n(7294).createContext)(void 0);t.Z=a},1190:function(e,t,n){n.r(t),n.d(t,{frontMatter:function(){return l},contentTitle:function(){return u},metadata:function(){return d},toc:function(){return c},default:function(){return p}});var a=n(7462),r=n(3366),i=(n(7294),n(3905)),o=(n(5064),n(8215),["components"]),l={id:"locationdata",title:"Location data",sidebar_label:"Location data"},u=void 0,d={unversionedId:"workingwithpopdata/locationdata",id:"workingwithpopdata/locationdata",isDocsHomePage:!1,title:"Location data",description:"Also known as geographic or coordinate data. The sampleinfo in the metadata can contain this information, which is not currently used in",source:"@site/docs/workingwithpopdata/locationdata.md",sourceDirName:"workingwithpopdata",slug:"/workingwithpopdata/locationdata",permalink:"/PopGen.jl/docs/workingwithpopdata/locationdata",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/docs/workingwithpopdata/locationdata.md",tags:[],version:"current",lastUpdatedAt:1635820398,formattedLastUpdatedAt:"11/2/2021",frontMatter:{id:"locationdata",title:"Location data",sidebar_label:"Location data"},sidebar:"docs",previous:{title:"Population data",permalink:"/PopGen.jl/docs/workingwithpopdata/populationdata"},next:{title:"Data exploration",permalink:"/PopGen.jl/docs/workingwithpopdata/dataexploration"}},c=[{value:"View location data",id:"view-location-data",children:[],level:2},{value:"Add geographical coordinates",id:"add-geographical-coordinates",children:[{value:"decimal minutes",id:"decimal-minutes",children:[{value:"formatting requirements",id:"formatting-requirements",children:[],level:4}],level:3},{value:"other formats",id:"other-formats",children:[{value:"formatting requirements",id:"formatting-requirements-1",children:[],level:4}],level:3}],level:2}],s={toc:c};function p(e){var t=e.components,n=(0,r.Z)(e,o);return(0,i.kt)("wrapper",(0,a.Z)({},s,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("p",null,"Also known as geographic or coordinate data. The ",(0,i.kt)("inlineCode",{parentName:"p"},"sampleinfo")," in the metadata can contain this information, which is not currently used in\nany analyses present in PopGen.jl."),(0,i.kt)("h2",{id:"view-location-data"},"View location data"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"locationdata(data::PopData)\n")),(0,i.kt)("p",null,"View location (if present) data in a ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData"),",  returning a table the longitude and latitude information in the ",(0,i.kt)("inlineCode",{parentName:"p"},"metadata"),". "),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"julia> locationdata(sharks)\n212\xd72 SubDataFrame\n Row \u2502 longitude  latitude \n     \u2502 Float64    Float64  \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502   28.3062  -80.5993\n   2 \u2502   28.3079  -80.5995\n   3 \u2502   28.3023  -80.5996\n  \u22ee  \u2502     \u22ee         \u22ee\n 210 \u2502   30.0522  -87.3662\n 211 \u2502   29.8234  -85.7143\n 212 \u2502   29.8234  -85.7143\n           206 rows omitted\n")),(0,i.kt)("h2",{id:"add-geographical-coordinates"},"Add geographical coordinates"),(0,i.kt)("h3",{id:"decimal-minutes"},"decimal minutes"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"locationdata!(data::PopData; longitude::Vector{T}, latitude::Vector{T}) where T<:AbstractFloat\n")),(0,i.kt)("p",null,"Location data can be added using one of the methods of ",(0,i.kt)("inlineCode",{parentName:"p"},"locations!"),". As indicated by the bang ",(0,i.kt)("inlineCode",{parentName:"p"},"!"),", your ",(0,i.kt)("inlineCode",{parentName:"p"},"PopData")," will be edited in place, and there will be no return output. If your data is in anything other than Decimal-Degrees format, this function will convert your long/lat into Decimal Degrees. To import those data into Julia, you'll likely want to use the wonderful ",(0,i.kt)("inlineCode",{parentName:"p"},"CSV.jl")," package first. The functions accept keywords ",(0,i.kt)("inlineCode",{parentName:"p"},"longitude")," and ",(0,i.kt)("inlineCode",{parentName:"p"},"latitude"),", or can be used without them so long as the vectors are input in that order. "),(0,i.kt)("p",null,"This method is pretty straightforward and tolerates vectors with ",(0,i.kt)("inlineCode",{parentName:"p"},"missing")," data."),(0,i.kt)("h4",{id:"formatting-requirements"},"formatting requirements"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Coordinates must be decimal-minutes as either ",(0,i.kt)("inlineCode",{parentName:"li"},"Float32")," or ",(0,i.kt)("inlineCode",{parentName:"li"},"Float64")," (e.g. ",(0,i.kt)("inlineCode",{parentName:"li"},"-21.321"),")")),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"# generate some fake location data\njulia> long = rand(212) .* 10 ; lat = rand(212) .* -10\n\njulia> locationdata!(sharks, long, lat)\n")),(0,i.kt)("h3",{id:"other-formats"},"other formats"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},"locationdata!(data::PopData; longitude::Vector{String}, latitude::Vector{String})\n")),(0,i.kt)("h4",{id:"formatting-requirements-1"},"formatting requirements"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Coordinates as ",(0,i.kt)("inlineCode",{parentName:"li"},"String")," separated by spaces (",(0,i.kt)("inlineCode",{parentName:"li"},'"11 43 41"'),") or colons (",(0,i.kt)("inlineCode",{parentName:"li"},'"11:43:41"'),")"),(0,i.kt)("li",{parentName:"ul"},"Must use negative sign (",(0,i.kt)("inlineCode",{parentName:"li"},'"-11 43.52"'),") or single-letter cardinal direction (",(0,i.kt)("inlineCode",{parentName:"li"},'"11 43.52W"'),")"),(0,i.kt)("li",{parentName:"ul"},"Missing data should be coded as the string ",(0,i.kt)("inlineCode",{parentName:"li"},'"missing"')," (can be accomplished with ",(0,i.kt)("inlineCode",{parentName:"li"},"replace!()"),")"),(0,i.kt)("li",{parentName:"ul"},"Can mix colons and spaces (although it's bad practice)")),(0,i.kt)("p",null,"If not already in decimal-minutes format, it would likely be most convenient if you imported your coordinate data as vectors of strings, which would look something like this:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-julia"},'longitude = ["-43 54 11", "22 23 11N"]\nlatitude = ["11 44 31", "-25 41 13"]\n')),(0,i.kt)("div",{className:"admonition admonition-caution alert alert--warning"},(0,i.kt)("div",{parentName:"div",className:"admonition-heading"},(0,i.kt)("h5",{parentName:"div"},(0,i.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,i.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"16",height:"16",viewBox:"0 0 16 16"},(0,i.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M8.893 1.5c-.183-.31-.52-.5-.887-.5s-.703.19-.886.5L.138 13.499a.98.98 0 0 0 0 1.001c.193.31.53.501.886.501h13.964c.367 0 .704-.19.877-.5a1.03 1.03 0 0 0 .01-1.002L8.893 1.5zm.133 11.497H6.987v-2.003h2.039v2.003zm0-3.004H6.987V5.987h2.039v4.006z"}))),"Missing values")),(0,i.kt)("div",{parentName:"div",className:"admonition-content"},(0,i.kt)("p",{parentName:"div"},"This method tolerates ",(0,i.kt)("inlineCode",{parentName:"p"},"missing")," values, but you will need to ",(0,i.kt)("inlineCode",{parentName:"p"},"replace!")," instances of ",(0,i.kt)("inlineCode",{parentName:"p"},"missing")," with the string ",(0,i.kt)("inlineCode",{parentName:"p"},'"missing"'),"."))))}p.isMDXComponent=!0}}]);