"use strict";(self.webpackChunkpop_gen_jl=self.webpackChunkpop_gen_jl||[]).push([[6524],{4137:function(e,t,a){a.d(t,{Zo:function(){return d},kt:function(){return h}});var n=a(7294);function i(e,t,a){return t in e?Object.defineProperty(e,t,{value:a,enumerable:!0,configurable:!0,writable:!0}):e[t]=a,e}function l(e,t){var a=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),a.push.apply(a,n)}return a}function r(e){for(var t=1;t<arguments.length;t++){var a=null!=arguments[t]?arguments[t]:{};t%2?l(Object(a),!0).forEach((function(t){i(e,t,a[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(a)):l(Object(a)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(a,t))}))}return e}function s(e,t){if(null==e)return{};var a,n,i=function(e,t){if(null==e)return{};var a,n,i={},l=Object.keys(e);for(n=0;n<l.length;n++)a=l[n],t.indexOf(a)>=0||(i[a]=e[a]);return i}(e,t);if(Object.getOwnPropertySymbols){var l=Object.getOwnPropertySymbols(e);for(n=0;n<l.length;n++)a=l[n],t.indexOf(a)>=0||Object.prototype.propertyIsEnumerable.call(e,a)&&(i[a]=e[a])}return i}var o=n.createContext({}),u=function(e){var t=n.useContext(o),a=t;return e&&(a="function"==typeof e?e(t):r(r({},t),e)),a},d=function(e){var t=u(e.components);return n.createElement(o.Provider,{value:t},e.children)},p="mdxType",c={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},m=n.forwardRef((function(e,t){var a=e.components,i=e.mdxType,l=e.originalType,o=e.parentName,d=s(e,["components","mdxType","originalType","parentName"]),p=u(a),m=i,h=p["".concat(o,".").concat(m)]||p[m]||c[m]||l;return a?n.createElement(h,r(r({ref:t},d),{},{components:a})):n.createElement(h,r({ref:t},d))}));function h(e,t){var a=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var l=a.length,r=new Array(l);r[0]=m;var s={};for(var o in t)hasOwnProperty.call(t,o)&&(s[o]=t[o]);s.originalType=e,s[p]="string"==typeof e?e:i,r[1]=s;for(var u=2;u<l;u++)r[u]=a[u];return n.createElement.apply(null,r)}return n.createElement.apply(null,a)}m.displayName="MDXCreateElement"},1483:function(e,t,a){a.r(t),a.d(t,{assets:function(){return c},contentTitle:function(){return d},default:function(){return b},frontMatter:function(){return u},metadata:function(){return p},toc:function(){return m}});var n=a(7462),i=a(3366),l=(a(7294),a(4137)),r=a(3992),s=a(425),o=["components"],u={slug:"relatedness",title:"Relatedness Tutorial",author:"Pavel Dimens",author_title:"Little this, little that",author_url:"https://github.com/pdimens",author_image_url:"https://avatars1.githubusercontent.com/u/19176506?s=460&u=3afad1d1ef3b09ddc4ab7108143f515be3412d5a&v=4",tags:["tutorials"]},d=void 0,p={permalink:"/PopGen.jl/blog/relatedness",editUrl:"https://github.com/BioJulia/PopGen.jl/edit/documentation/blog/2020-11-14-relatedness.md",source:"@site/blog/2020-11-14-relatedness.md",title:"Relatedness Tutorial",description:"The kinship interface has changed a bit between versions 0.7 and 0.9. This",date:"2020-11-14T00:00:00.000Z",formattedDate:"November 14, 2020",tags:[{label:"tutorials",permalink:"/PopGen.jl/blog/tags/tutorials"}],readingTime:10.265,hasTruncateMarker:!1,authors:[{name:"Pavel Dimens",title:"Little this, little that",url:"https://github.com/pdimens",imageURL:"https://avatars1.githubusercontent.com/u/19176506?s=460&u=3afad1d1ef3b09ddc4ab7108143f515be3412d5a&v=4"}],frontMatter:{slug:"relatedness",title:"Relatedness Tutorial",author:"Pavel Dimens",author_title:"Little this, little that",author_url:"https://github.com/pdimens",author_image_url:"https://avatars1.githubusercontent.com/u/19176506?s=460&u=3afad1d1ef3b09ddc4ab7108143f515be3412d5a&v=4",tags:["tutorials"]},prevItem:{title:"A bit about FST",permalink:"/PopGen.jl/blog/fst"},nextItem:{title:"Logo Graveyard",permalink:"/PopGen.jl/blog/logos"}},c={authorsImageUrls:[void 0]},m=[{value:"Getting Started",id:"getting-started",level:2},{value:"Estimators",id:"estimators",level:3},{value:"The Steps",id:"the-steps",level:2},{value:"1. Calculate pairwise relatedness",id:"1-calculate-pairwise-relatedness",level:3},{value:"2. Bootstrap to calculate CI",id:"2-bootstrap-to-calculate-ci",level:3},{value:"3. Create CI&#39;s for the sibships",id:"3-create-cis-for-the-sibships",level:3},{value:"i. simulate known sibship pairs",id:"i-simulate-known-sibship-pairs",level:4},{value:"ii. get relatedness estimates for the simulated data",id:"ii-get-relatedness-estimates-for-the-simulated-data",level:4},{value:"iii. sibship intervals",id:"iii-sibship-intervals",level:4},{value:"4. Finally, the data assessment",id:"4-finally-the-data-assessment",level:3},{value:"Closing remarks",id:"closing-remarks",level:3}],h={toc:m},f="wrapper";function b(e){var t=e.components,u=(0,i.Z)(e,o);return(0,l.kt)(f,(0,n.Z)({},h,u,{components:t,mdxType:"MDXLayout"}),(0,l.kt)("admonition",{title:"PopGen.jl <0.9.0",type:"info"},(0,l.kt)("p",{parentName:"admonition"},"The kinship interface has changed a bit between versions 0.7 and 0.9. This\npost has not yet been updated for versions 0.9.0+. To follow along, use versions 0.8.0 or lower.")),(0,l.kt)("h2",{id:"getting-started"},"Getting Started"),(0,l.kt)("p",null,"In a population genetics study, you often need to identify if there are kin in your data. This may be necessary because you are trying to remove kin from your data (because of Hardy-Weinberg assumptions), or maybe kinship is a central interest in your study. Either way, the goal of this tutorial is to provide you with a basic tutorial on using PopGen.jl to perform a relatedness analysis, which is sometimes called a ",(0,l.kt)("em",{parentName:"p"},"kinship")," analysis. To follow along, you'll need to have Julia, along with the packages PopGen.jl, PopGenSims.jl, and StatsBase.jl installed. We'll be using the ",(0,l.kt)("inlineCode",{parentName:"p"},"nancycats")," data because it's smaller than ",(0,l.kt)("inlineCode",{parentName:"p"},"gulfsharks"),", so things should be a lot quicker."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"using PopGen, PopGenSims, StatsBase\n\njulia> cats = @nancycats\nPopData{Diploid, 9 Microsatellite Loci}\n  Samples: 237\n  Populations: 17\n")),(0,l.kt)("h3",{id:"estimators"},"Estimators"),(0,l.kt)("p",null,"Like ",(0,l.kt)("inlineCode",{parentName:"p"},"Coancestry")," and the R packages that wrap it (i.e. ",(0,l.kt)("inlineCode",{parentName:"p"},"relate"),", ",(0,l.kt)("inlineCode",{parentName:"p"},"related"),"), PopGen.jl provides a whole bunch of relatedness estimators that you can choose from for your data. Unfortunately, there is no right answer and you will need to use your discretion. Some people choose an estimator based on the heterozygosity of the data, others choose one based on more liberal or conservative values, and there are yet more criteria one can consider for choosing an estimator. To keep things simple, we're going to use ",(0,l.kt)("inlineCode",{parentName:"p"},"LynchLi"),". Why? Because I'm the one writing this tutorial, and I said so \ud83d\ude01. "),(0,l.kt)("h2",{id:"the-steps"},"The Steps"),(0,l.kt)("h3",{id:"1-calculate-pairwise-relatedness"},"1. Calculate pairwise relatedness"),(0,l.kt)("p",null,"This seems pretty obvious, but it needs to be said. In order to do the analysis, you need you get the pairwise relatedness values for the individuals of interest in your data. To keep things simple, we're going to do an all-by-all comparison. But, we also want to boostrap the pairs to create confidence intervals (\"CI\") for each pair, so let's talk about that."),(0,l.kt)("h3",{id:"2-bootstrap-to-calculate-ci"},"2. Bootstrap to calculate CI"),(0,l.kt)("p",null,"It would behoove us to bootstrap the loci in a pairwise comparison ",(0,l.kt)("em",{parentName:"p"},"n")," number of times so we can create a confidence interval for the relatedness estimates for each pair. This inflates the runtime of the analysis substantially, but it's a very useful method in making sense of our estimates. If one was to be conservative (and we generally are), then we would reject an estimate for a pair whose CI includes zero. In PopGen.jl, the estimates and bootstrapping are done all at once to minimize processing time, so the command for that would be"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> rel_out  = kinship(cats, method = LynchLi, iterations = 1000)\n")),(0,l.kt)("p",null,"By default, the ",(0,l.kt)("inlineCode",{parentName:"p"},"kinship")," function uses a 95% CI (",(0,l.kt)("inlineCode",{parentName:"p"},"interval = (0.025, 0.975)"),"), but you can change that with ",(0,l.kt)("inlineCode",{parentName:"p"},"interval = (low,high)")," where ",(0,l.kt)("inlineCode",{parentName:"p"},"low")," and ",(0,l.kt)("inlineCode",{parentName:"p"},"high")," are decimals of your quantiles.\nThe ",(0,l.kt)("inlineCode",{parentName:"p"},"kinship()")," function returns a ",(0,l.kt)("inlineCode",{parentName:"p"},"NamedTuple")," of dataframes whenever you are bootstrapping, where each element shares its name with the method used, so in this case, we can access our results with ",(0,l.kt)("inlineCode",{parentName:"p"},"rel_out.LynchLi"),"."),(0,l.kt)(r.Z,{block:!0,defaultValue:"rel",values:[{label:"relatedness",value:"rel"},{label:"plotting",value:"pl"}],mdxType:"Tabs"},(0,l.kt)(s.Z,{value:"rel",mdxType:"TabItem"},(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre"},"julia> rel_out.LynchLi\n27966\xd78 DataFrame\n   Row \u2502 sample_1  sample_2  n_loci  LynchLi     LynchLi_mean  LynchLi_median  LynchLi_S \u22ef\n       \u2502 String    String    Int64   Float64?    Float64?      Float64?        Float64?  \u22ef\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n     1 \u2502 N215      N216           8   0.743535      0.747288        0.748042     0.75344 \u22ef\n     2 \u2502 N215      N217           8   0.230605      0.233593        0.240085     0.34187\n     3 \u2502 N215      N218           8   0.230605      0.230507        0.221861     0.32161\n     4 \u2502 N215      N219           8   0.230605      0.23601         0.23567      0.32782\n     5 \u2502 N215      N220           8   0.333191      0.33798         0.350492     0.39898 \u22ef\n     6 \u2502 N215      N221           8   0.589656      0.594223        0.601308     0.61945\n     7 \u2502 N215      N222           8   0.0254328     0.0347216       0.0262021    0.21408\n     8 \u2502 N215      N223           8   0.333191      0.329983        0.331411     0.38402\n     9 \u2502 N215      N224           8  -0.0258602    -0.021062       -0.0301579    0.21112 \u22ef\n    10 \u2502 N215      N7             8  -0.282325     -0.27967        -0.288337     0.33611\n    11 \u2502 N215      N141           8  -0.0771532    -0.0796867      -0.083113     0.21261\n    12 \u2502 N215      N142           8   0.0254328     0.0302549       0.0330718    0.23957\n   \u22ee   \u2502    \u22ee         \u22ee        \u22ee         \u22ee            \u22ee              \u22ee             \u22ee     \u22f1\n 27955 \u2502 N295      N289           7   0.322731      0.347021        0.34118      0.41168 \u22ef\n 27956 \u2502 N295      N290           7   0.153414      0.160102        0.164866     0.22862\n 27957 \u2502 N296      N297           7  -0.0159038    -0.0182747      -0.0187108    0.16981\n 27958 \u2502 N296      N281           7   0.0405353     0.037025        0.0422647    0.15294\n 27959 \u2502 N296      N289           7   0.322731      0.328379        0.337317     0.35578 \u22ef\n 27960 \u2502 N296      N290           7   0.153414      0.152384        0.16194      0.19131\n 27961 \u2502 N297      N281           7  -0.0159038    -0.0128349      -0.0303449    0.21030\n 27962 \u2502 N297      N289           7   0.37917       0.392517        0.392818     0.45139\n 27963 \u2502 N297      N290           7   0.435609      0.437829        0.450044     0.47027 \u22ef\n 27964 \u2502 N281      N289           8   0.20428       0.21279         0.207425     0.29611\n 27965 \u2502 N281      N290           7   0.37917       0.386583        0.390585     0.45471\n 27966 \u2502 N289      N290           7   0.209853      0.217811        0.222894     0.28649\n                                                          2 columns and 27942 rows omitted\n"))),(0,l.kt)(s.Z,{value:"pl",mdxType:"TabItem"},(0,l.kt)("p",null,"And while it's totally optional, we can plot the distribution of values for some visual data exploration. For that we'll use Plots.jl and StatsPlots.jl. "),(0,l.kt)("admonition",{title:"plotting packages",type:"note"},(0,l.kt)("p",{parentName:"admonition"},"We could have used any plotting package, but Plots.jl was chosen for simplicity. Other great\noptions are (and not limited to): Makie.jl, Gadfly.jl, VegaLite.jl, and PlotlyJS.jl.")),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'using Plots, StatsPlots\n\njulia> density(rel_out.LynchLi.LynchLi)\njulia> title!("LynchLi relatedness")\n')),(0,l.kt)("p",null,(0,l.kt)("img",{alt:"relatedness_histo",src:a(3051).Z,width:"600",height:"400"})))),(0,l.kt)("h3",{id:"3-create-cis-for-the-sibships"},"3. Create CI's for the sibships"),(0,l.kt)("admonition",{title:"There's more???",type:"note"},(0,l.kt)("p",{parentName:"admonition"},"The next set of steps seem like a lot more work, so allow me to explain. The estimators all generally give you some value between 0-1 (or 0-0.5, same idea) and you can intuit that certain values mean certain things, like that ",(0,l.kt)("inlineCode",{parentName:"p"},"0"),' is "unrelated", ',(0,l.kt)("inlineCode",{parentName:"p"},"0.25"),' is "half-sib", and ',(0,l.kt)("inlineCode",{parentName:"p"},"0.5"),' is "full-sib". However, those are fixed values, so how do we know how far we can deviate from 0.25 (for example) and still call our pair half-siblings? Instead of hand-waving, we can create confidence intervals from simulated data to act as sibship ranges for our data. If this doesn\'t make sense yet, it will below. Promise! ')),(0,l.kt)("h4",{id:"i-simulate-known-sibship-pairs"},"i. simulate known sibship pairs"),(0,l.kt)("p",null,'Next, we need to further contextualize what out estimates actually mean, and we need to devise a rigorous and defensible method to define the sibling-ship ("sibship") of a pair of samples as ',(0,l.kt)("strong",{parentName:"p"},"unrelated"),", ",(0,l.kt)("strong",{parentName:"p"},"half-sibs"),", or ",(0,l.kt)("strong",{parentName:"p"},"full-sibs"),". To do this, we're going to use PopGenSims.jl to simulate data based on the allele frequencies in our data. What ",(0,l.kt)("inlineCode",{parentName:"p"},"simulatekin")," does is simulate individuals based on the allele frequencies in your ",(0,l.kt)("inlineCode",{parentName:"p"},"PopData"),', then simulate offspring of a particular siblingship resulting from the "mating" of those individuals. For example, when you specify ',(0,l.kt)("inlineCode",{parentName:"p"},'"fullsib"'),", it generates two offspring resulting from two parents, ",(0,l.kt)("inlineCode",{parentName:"p"},"n")," number of times. We want to do this for each of the three kinds of sibships."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> kin_sims = simulatekin(cats, fullsib = 500, halfsib = 500, unrelated = 500)\nPopData{Diploid, 9 Microsatellite Loci}\n  Samples: 3000\n  Populations: 3\n")),(0,l.kt)("p",null,"We can keep all three simulated relationships together, but it will be easier to explain things (for instructional purposes) if we don't, so let's split out each into its own PopData."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'julia> fullsib = kin_sims[genodata(kin_sims).population .== "fullsib"] ;\njulia> halfsib = kin_sims[genodata(kin_sims).population .== "halfsib"] ;\njulia> unrelated = kin_sims[genodata(kin_sims).population .== "unrelated"] ;\n')),(0,l.kt)("h4",{id:"ii-get-relatedness-estimates-for-the-simulated-data"},"ii. get relatedness estimates for the simulated data"),(0,l.kt)("p",null,"Next, we want to get the relatedness estimate for each simulated pair of \"known\" sibship. We are only interested in the values for the simulated pairs and not samples across pairs. If you aren't sure why that is, think of it this way: we're trying to create a range of values where we can confidently say unknown things are full-sibs (or half-sib, etc.), so we want to know what range of values we get from a bunch of known fullsib pairs, not the unknown relationships of samples between pairs. "),(0,l.kt)("p",null,"It would a nightmare to manually specify only 2 individuals at a time into ",(0,l.kt)("inlineCode",{parentName:"p"},"kinship()"),". Instead, the function has a shortcut built into it that will recognize the ",(0,l.kt)("inlineCode",{parentName:"p"},"population")," names generated from ",(0,l.kt)("inlineCode",{parentName:"p"},"simulatekin")," and only give you relatedness estimates for those pairs. So, we just need to run it once for each of our simulations, this time ",(0,l.kt)("em",{parentName:"p"},"without")," bootstrapping because we are only interested in the estimates. Make sure to use the same estimator! The run will be a lot faster this time because it only needs to calculate estimates for 500 pairs each time (our ",(0,l.kt)("inlineCode",{parentName:"p"},"n")," from above) and without bootstrapping. When not bootstrapping, ",(0,l.kt)("inlineCode",{parentName:"p"},"kinship()")," returns a dataframe (versus a NamedTuple of dataframes)."),(0,l.kt)(r.Z,{block:!0,defaultValue:"un",values:[{label:"unrelated relatedness",value:"un"},{label:"halfsib relatedness",value:"h"},{label:"fullsib relatedness",value:"f"}],mdxType:"Tabs"},(0,l.kt)(s.Z,{value:"un",mdxType:"TabItem"},(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> un_sims_rel = kinship(unrelated_sims, method = LynchLi)\n500\xd74 DataFrame\n Row \u2502 sample_1            sample_2            n_loci  LynchLi    \n     \u2502 String              String              Int64   Float64?   \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 sim001_unrelated_1  sim001_unrelated_2       9  -0.11419\n   2 \u2502 sim002_unrelated_1  sim002_unrelated_2       9  -0.337028\n   3 \u2502 sim003_unrelated_1  sim003_unrelated_2       9  -0.0696222\n  \u22ee  \u2502         \u22ee                   \u22ee             \u22ee         \u22ee\n 498 \u2502 sim498_unrelated_1  sim498_unrelated_2       9   0.019513\n 499 \u2502 sim499_unrelated_1  sim499_unrelated_2       9   0.019513\n 500 \u2502 sim500_unrelated_1  sim500_unrelated_2       9   0.019513\n                                                  494 rows omitted\n"))),(0,l.kt)(s.Z,{value:"h",mdxType:"TabItem"},(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> half_sims_rel = kinship(halfsib_sims, method = LynchLi)\n500\xd74 DataFrame\n Row \u2502 sample_1          sample_2          n_loci  LynchLi    \n     \u2502 String            String            Int64   Float64?   \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 sim001_halfsib_1  sim001_halfsib_2       9  -0.0182636\n   2 \u2502 sim002_halfsib_1  sim002_halfsib_2       9   0.468732\n   3 \u2502 sim003_halfsib_1  sim003_halfsib_2       9   0.291643\n  \u22ee  \u2502        \u22ee                 \u22ee            \u22ee         \u22ee\n 498 \u2502 sim498_halfsib_1  sim498_halfsib_2       9   0.42446\n 499 \u2502 sim499_halfsib_1  sim499_halfsib_2       9   0.513004\n 500 \u2502 sim500_halfsib_1  sim500_halfsib_2       9   0.0702811\n                                              494 rows omitted\n"))),(0,l.kt)(s.Z,{value:"f",mdxType:"TabItem"},(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> full_sims_rel = kinship(fullsib_sims, method = LynchLi)\n500\xd74 DataFrame\n Row \u2502 sample_1          sample_2          n_loci  LynchLi  \n     \u2502 String            String            Int64   Float64? \n\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n   1 \u2502 sim001_fullsib_1  sim001_fullsib_2       9  0.732808\n   2 \u2502 sim002_fullsib_1  sim002_fullsib_2       9  0.599213\n   3 \u2502 sim003_fullsib_1  sim003_fullsib_2       9  0.376553\n  \u22ee  \u2502        \u22ee                 \u22ee            \u22ee        \u22ee\n 498 \u2502 sim498_fullsib_1  sim498_fullsib_2       9  0.376553\n 499 \u2502 sim499_fullsib_1  sim499_fullsib_2       9  0.510149\n 500 \u2502 sim500_fullsib_1  sim500_fullsib_2       9  0.109361\n                                            494 rows omitted\n")))),(0,l.kt)("p",null,"And if we wanted to plot what that looks like (optional):"),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},'using Plots, StatsPlots\n\njulia> density(rel_out.LynchLi.LynchLi, label = "real data", color = :grey, fill = (0, :grey))\njulia> density!(un_sims_rel.LynchLi, label = "unrelated")\njulia> density!(half_sims_rel.LynchLi, label = "halfsib", color = :blue)\njulia> density!(full_sims_rel.LynchLi, label = "fullsib", color = :black)\njulia> title!("relatedness estimates on simulated and real data")\n')),(0,l.kt)("p",null,(0,l.kt)("img",{alt:"relatedness_histo",src:a(1584).Z,width:"600",height:"400"})),(0,l.kt)("p",null,'Hopefully by now you are starting to contextualize why we\'re doing all of this. The distributions generated from our simulated data are giving us a better indication of what "unrelated", "halfsib", and "fullsib" estimates look like in our data.'),(0,l.kt)("h4",{id:"iii-sibship-intervals"},"iii. sibship intervals"),(0,l.kt)("p",null,"What we just did is create null distributions for each sibship relationship, so now all that's left is to get a confidence interval from each. Keep in mind that your values will be a bit different due to the randomization involved with many of these steps."),(0,l.kt)("pre",null,(0,l.kt)("code",{parentName:"pre",className:"language-julia"},"julia> unrelated_ci = quantile(un_sims_rel.LynchLi, (0.025, 0.975))\n(-0.3380513384854698, 0.33097433075726507)\n\njulia> halfsibs_ci = quantile(half_sims_rel.LynchLi, (0.025, 0.975))\n(-0.06652262584414452, 0.5556155725649398)\n\njulia> fullsibs_ci = quantile(full_sims_rel.LynchLi, (0.025, 0.975))\n(0.1989727186688347, 0.8219939374819633)\n")),(0,l.kt)("p",null,"So, given our data and the simulations we made, we can now make a reasonable assumption regarding the ranges for each sibship relationship:"),(0,l.kt)("table",null,(0,l.kt)("thead",{parentName:"table"},(0,l.kt)("tr",{parentName:"thead"},(0,l.kt)("th",{parentName:"tr",align:"left"},"Relationship"),(0,l.kt)("th",{parentName:"tr",align:"center"},"Lower Bound"),(0,l.kt)("th",{parentName:"tr",align:"center"},"Upper Bound"))),(0,l.kt)("tbody",{parentName:"table"},(0,l.kt)("tr",{parentName:"tbody"},(0,l.kt)("td",{parentName:"tr",align:"left"},"unrelated"),(0,l.kt)("td",{parentName:"tr",align:"center"},"-0.338"),(0,l.kt)("td",{parentName:"tr",align:"center"},"0.331")),(0,l.kt)("tr",{parentName:"tbody"},(0,l.kt)("td",{parentName:"tr",align:"left"},"half sibling"),(0,l.kt)("td",{parentName:"tr",align:"center"},"-0.066"),(0,l.kt)("td",{parentName:"tr",align:"center"},"0.555")),(0,l.kt)("tr",{parentName:"tbody"},(0,l.kt)("td",{parentName:"tr",align:"left"},"full sibling"),(0,l.kt)("td",{parentName:"tr",align:"center"},"0.199"),(0,l.kt)("td",{parentName:"tr",align:"center"},"0.82")))),(0,l.kt)("h3",{id:"4-finally-the-data-assessment"},"4. Finally, the data assessment"),(0,l.kt)("p",null,"Now that we have our relatedness estimates and the acceptable sibship ranges given our data, we can see where our data falls.\n",(0,l.kt)("em",{parentName:"p"},"Now"),", we can say that a particular sample pair is unrelated/halfsib/fullsib if:"),(0,l.kt)("ol",null,(0,l.kt)("li",{parentName:"ol"},"that pair's confidence interval does not include zero and "),(0,l.kt)("li",{parentName:"ol"},"that pair's estimate falls within any of the three calculate ranges")),(0,l.kt)("p",null,"Unfortunately, the ranges overlap quite a bit, which makes interpretation difficult, so it may be suitable to use a different estimator. "),(0,l.kt)("h3",{id:"closing-remarks"},"Closing remarks"),(0,l.kt)("p",null,"There is more that can be done for relatedness, like network analysis. However, this tutorial covers what we consider a reasonable way to assess kinship in population genetic studies. If removing kin is your ultimate goal, consider the effects doing that may have on the analyses you are looking to do. Additionally, consider which individual of a pair you would remove and why. If you were curious as to how many identical loci a pair of samples may share, you can check that using ",(0,l.kt)("inlineCode",{parentName:"p"},"pairwiseidentical()"),". Good luck!"))}b.isMDXComponent=!0},425:function(e,t,a){a.d(t,{Z:function(){return r}});var n=a(7294),i=a(6010),l={tabItem:"tabItem_Ymn6"};function r(e){var t=e.children,a=e.hidden,r=e.className;return n.createElement("div",{role:"tabpanel",className:(0,i.Z)(l.tabItem,r),hidden:a},t)}},3992:function(e,t,a){a.d(t,{Z:function(){return w}});var n=a(7462),i=a(7294),l=a(6010),r=a(2957),s=a(6550),o=a(5238),u=a(3609),d=a(2560);function p(e){return function(e){var t,a;return null!=(t=null==(a=i.Children.map(e,(function(e){if(!e||(0,i.isValidElement)(e)&&(t=e.props)&&"object"==typeof t&&"value"in t)return e;var t;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})))?void 0:a.filter(Boolean))?t:[]}(e).map((function(e){var t=e.props;return{value:t.value,label:t.label,attributes:t.attributes,default:t.default}}))}function c(e){var t=e.values,a=e.children;return(0,i.useMemo)((function(){var e=null!=t?t:p(a);return function(e){var t=(0,u.l)(e,(function(e,t){return e.value===t.value}));if(t.length>0)throw new Error('Docusaurus error: Duplicate values "'+t.map((function(e){return e.value})).join(", ")+'" found in <Tabs>. Every value needs to be unique.')}(e),e}),[t,a])}function m(e){var t=e.value;return e.tabValues.some((function(e){return e.value===t}))}function h(e){var t=e.queryString,a=void 0!==t&&t,n=e.groupId,l=(0,s.k6)(),r=function(e){var t=e.queryString,a=void 0!==t&&t,n=e.groupId;if("string"==typeof a)return a;if(!1===a)return null;if(!0===a&&!n)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return null!=n?n:null}({queryString:a,groupId:n});return[(0,o._X)(r),(0,i.useCallback)((function(e){if(r){var t=new URLSearchParams(l.location.search);t.set(r,e),l.replace(Object.assign({},l.location,{search:t.toString()}))}}),[r,l])]}function f(e){var t,a,n,l,r=e.defaultValue,s=e.queryString,o=void 0!==s&&s,u=e.groupId,p=c(e),f=(0,i.useState)((function(){return function(e){var t,a=e.defaultValue,n=e.tabValues;if(0===n.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(a){if(!m({value:a,tabValues:n}))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+a+'" but none of its children has the corresponding value. Available values are: '+n.map((function(e){return e.value})).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");return a}var i=null!=(t=n.find((function(e){return e.default})))?t:n[0];if(!i)throw new Error("Unexpected error: 0 tabValues");return i.value}({defaultValue:r,tabValues:p})})),b=f[0],g=f[1],k=h({queryString:o,groupId:u}),y=k[0],v=k[1],w=(t=function(e){return e?"docusaurus.tab."+e:null}({groupId:u}.groupId),a=(0,d.Nk)(t),n=a[0],l=a[1],[n,(0,i.useCallback)((function(e){t&&l.set(e)}),[t,l])]),N=w[0],_=w[1],j=function(){var e=null!=y?y:N;return m({value:e,tabValues:p})?e:null}();return(0,i.useLayoutEffect)((function(){j&&g(j)}),[j]),{selectedValue:b,selectValue:(0,i.useCallback)((function(e){if(!m({value:e,tabValues:p}))throw new Error("Can't select invalid tab value="+e);g(e),v(e),_(e)}),[v,_,p]),tabValues:p}}var b=a(1048),g={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};function k(e){var t=e.className,a=e.block,s=e.selectedValue,o=e.selectValue,u=e.tabValues,d=[],p=(0,r.o5)().blockElementScrollPositionUntilNextRender,c=function(e){var t=e.currentTarget,a=d.indexOf(t),n=u[a].value;n!==s&&(p(t),o(n))},m=function(e){var t,a=null;switch(e.key){case"Enter":c(e);break;case"ArrowRight":var n,i=d.indexOf(e.currentTarget)+1;a=null!=(n=d[i])?n:d[0];break;case"ArrowLeft":var l,r=d.indexOf(e.currentTarget)-1;a=null!=(l=d[r])?l:d[d.length-1]}null==(t=a)||t.focus()};return i.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,l.Z)("tabs",{"tabs--block":a},t)},u.map((function(e){var t=e.value,a=e.label,r=e.attributes;return i.createElement("li",(0,n.Z)({role:"tab",tabIndex:s===t?0:-1,"aria-selected":s===t,key:t,ref:function(e){return d.push(e)},onKeyDown:m,onClick:c},r,{className:(0,l.Z)("tabs__item",g.tabItem,null==r?void 0:r.className,{"tabs__item--active":s===t})}),null!=a?a:t)})))}function y(e){var t=e.lazy,a=e.children,n=e.selectedValue,l=(Array.isArray(a)?a:[a]).filter(Boolean);if(t){var r=l.find((function(e){return e.props.value===n}));return r?(0,i.cloneElement)(r,{className:"margin-top--md"}):null}return i.createElement("div",{className:"margin-top--md"},l.map((function(e,t){return(0,i.cloneElement)(e,{key:t,hidden:e.props.value!==n})})))}function v(e){var t=f(e);return i.createElement("div",{className:(0,l.Z)("tabs-container",g.tabList)},i.createElement(k,(0,n.Z)({},e,t)),i.createElement(y,(0,n.Z)({},e,t)))}function w(e){var t=(0,b.Z)();return i.createElement(v,(0,n.Z)({key:String(t)},e))}},3051:function(e,t,a){t.Z=a.p+"assets/images/nancycats_LynchLi-8b8fb132d841ad2cdff25f8b4d25c16d.png"},1584:function(e,t,a){t.Z=a.p+"assets/images/nancycats_sims-d891d2c0c1dcf2f3fb45784ca0664959.png"}}]);