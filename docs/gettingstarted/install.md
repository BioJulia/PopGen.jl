---
slug: /
title: Installation
sidebar_label: Installation
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

Installation is simple!

<Tabs
  block={true}
  defaultValue="repl"
  values={[
    { label: 'REPL/Juno', value: 'repl', },
    { label: 'Jupyter/nteract', value: 'jupyter', },
  ]
}>
<TabItem value="repl">

Invoke the package manager with `]` in the REPL and use

```julia
add PopGen
```

*![install](/img/install_repl.gif)*

</TabItem>
<TabItem value="jupyter">

Slightly different than the REPL, you will need to invoke `using Pkg`.

```julia
using Pkg
Pkg.add("PopGen")
```

*![install](/img/install_jupyter.gif)*

</TabItem>
</Tabs>

## Using PopGen

Like all Julia packages, to activate PopGen.jl, simply run:

```julia
julia> using PopGen
```

Feel free to play around with the test data in `/data/` or add it to your workspace with the `@nancycats` and `@gulfsharks` commands.

<details>
<summary>For ArchLinux users</summary>

If you compiled Julia from source, your PopGen.jl installation may fail due to incorrectly building `Arpack`, which is expected to be in one place, but the compilation puts in another.

**Solutions**:

- install official Julia binaries from the AUR (`julia-bin`), which includes a correctly bundled `Arpack` (recommended)
- if Julia was compiled from source: install `julia-arpack` from the AUR and make sure to delete `~/.julia/packages/Arpack` if it exists. That *should* fix things, but sometimes it still acts up.

</details>


## "Installing" the docs
The magic of Docusaurus (the engine these docs are built on) lets you "install" this documentation as a progressive web app (fancy web shortcut) if you're using a Chromium-based browser (e.g. Chrome, Chromium, Brave, Vivaldi, Edge). Doing so gives you quick-access to these docs as easy as double-clicking an application. Here's what that looks like using the Edge browser on Windows:

![edge_pwa](/img/pwa_screenshot.png)

which then shows up in your Start menu like a plain-regular application!

![pwa_app](/img/pwa_icon.png)