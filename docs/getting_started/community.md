---
id: community
title: Community
sidebar_label: Community
---

## Join the Slack channel!

Have questions, want to talk population genetics, or have ideas?

 [![alt text](https://img.shields.io/badge/slack-join%20PopGen.jl-9d72b1?style=for-the-badge&logo=slack)](https://join.slack.com/t/popgenjl/shared_invite/zt-deam65n8-DuBs2z1oDtsbBuRplJW~Pg)

:::note Community Standards
The PopGen.jl community is committed to maintaining a welcoming, inclusive, civil, and constructive environment. We expect the following standards to be observed and upheld by all participants during all forms of communication.

#### Be respectful and inclusive
Do not attack anyone based on any aspect of personal identity, including but not limited to: gender, sexuality, politics, religion, ethnicity, race, age, or ability. Refrain from making prejudiced or sexual jokes and comments or using overtly sexual language or imagery. 

If you believe one of these standards has been violated, you can either file an issue on an appropriate repository or confidentially contact Pavel in the Slack workspace. Although many mistakes are due to ignorance rather than malice, violation of these standards may result in the removal of an individual or individuals from the workspace. Please do your part in helping the community and the package grow.

Adopted from the [Julia Community Standards](https://julialang.org/community/standards/).
:::

## Contributing
We absolutely welcome contributors to this package/project! See below for ways you can help PopGen.jl grow.

### Improving available code

By nature, most (if not all) users of PopGen.jl will be biologists, and will not come from a strong computational background. Heck, the co-creators of the package don't even have a strong computation background (we're very good at nodding and smiling). There are various [best-practices](https://docs.julialang.org/en/v1/manual/style-guide/index.html) documented in Julia that help users write code to get the best performance out of the language, and sometimes a function that gets written isn't as performant as it can be. If you come up with faster and/or more memory efficient alternatives to the functions in PopGen.jl, we wholeheartedly encourage you to open up a issue or pull request and we'll try to integrate that into the PopGen.jl ecosystem.

### Pointing out bugs 
No one is perfect (except maybe Betty White), so it's very possible that mistakes get made, especially as more users begin adopting and contributing to PopGen.jl (we hope). If something isn't functioning correctly, please open an issue **that includes a minimum working example**. We definitely want to figure out what's going on, and we need as much information as possible to try to diagnose issues as they come up.

Julia's overall performance is rooted in really clever dispatching to the LLVM compiler, and the underlying system benefits most from stricter typing, wherever we can enforce it, and we expect that this strictness may result in extra issues being raised from niche use-cases. 

#### Testing your data against the available code

During PopGen.jl's development, we generally use `nancycats` and `gulfsharks` for just about everything. Those two datasets cover the range of what we expect are common use-cases; a smaller microsatellite dataset (nancycats) and a medium sized SNP dataset (gulfsharks). However, we know better than to put all of our faith into _n_ = 2, so please use your own data and mess around with PopGen.jl and let us know if something doesn't jive.

### Writing new functions or methods for existing functions

By all means, extend PopGen.jl to include all sorts of analyses! GST, Tajima's D, AMOVA, SAMOVA, porting `BOTTLENECK`, coalescence, etc. Yes please! Part of the intent behind PopGen.jl is to have it act as a sort of sandbox to play around in, which is why most of the package's basic calculations are modularized (see the hidden API). As a design choice, it made the most sense to have written the package in this way, because many population genetic calculations are built on top of other ones (like allele frequency or heterozygosity).

### Writing or editing the docs

The documentation of PopGen.jl must be approachable and _helpful_. Helpful in the sense that if someone was tasked with "looking for siblings" in their data and didn't know exactly how, that the documentation for the `relatedness` command would provide some kind of _context_ as to what it does and how, with figures, and with helpful in-line links to the source publications or more in-depth online resources. This package is not intended for only expert-level users, which means that the documentation **needs to be accessible to entry-level users**. By no means will the documentation be the ultimate compendium of all population genetics knowledge, but it _will_ be helpful beyond simply stating a command, its arguments, and a one-liner of what it does.

If you want to contribute but don't feel comfortable with the programming side of things, then we encourage you to help grow the documentation. Clarify the language in some sections, or maybe provide a useful diagram where one doesn't exist. Or typos. There are always typoes.

### Spreading the word

The very soul of open-source projects relies on people wanting to get involved. Spread the word :smiley:. If you are into social media, give us a shoutout. We can't imagine why you'd want to, but if you insisted on using a hashtag, then #PopGenjl is probably the sensible choice. 

### Kindness and encouragement

This stuff can be hard. As the package grows, we expect that we'll be dealing with a growing number of issues/complaints. A little thumbs-up :thumbsup: or prayer-hands :pray: emoji (Pavel's personal favorite) can go a long way. Or cook up your own cactus graphic (like the one below) and send it to us. Who doesn't love adorable cactuses doing human things?

![error_cactus](/img/terminal_cactus.png)