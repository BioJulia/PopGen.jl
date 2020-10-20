---
id: julia_primer
title: A quick Julia primer
sidebar_label: A quick Julia primer
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

*For getting the most out of this documentation*

There is nothing inherently special about this documentation relative to other documentation, other than we really *really* want you to get the most out of what's written here. This means that we need to embrace the fact that both novice and experienced Julia users may be reading these docs and using this package. So let's cover some Julia basics that will really help in navigating this package before we even get into the complicated genetics stuff. This primer is by no means "everything you need to get started in Julia", and is a poor substitute for actually learning the language. In general, we recommend [Think Julia: How to Think Like a Data Scientist](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html) by Ben Lauwens to establish some solid Julia foundations. It's free online! Also, the Julia language maintains [its own great documentation](https://docs.julialang.org/en/v1/) that we rely on quite heavily for development.


## Using Julia

Everyone has their own particular workflows, and if you're new to Julia, you might not have established one yet. Julia can be used rather comfortably using its built-in interpreter. For an RStudio-like experience, we recommend using the [VScode Julia extension](https://www.julia-vscode.org), but you can also use [Atom](https://junolab.org/) add-on). If you're already a fan of Jupyter notebooks (or [**nteract**](https://nteract.io/)), then all you need is to install the `IJulia` package in Julia and you have full Jupyter support for Julia! You can also use the new reactive notebooks provided by [Pluto.jl](https://github.com/fonsp/Pluto.jl).

:::note Trivia
If you didn't already know,  the name "Jupyter" is actually a concatenation of **Ju** (julia) **Pyt** (python) and **eR** (R). 🤯
:::

## First-time Performance
If you're migrating to Julia from Python or R (or Matlab, etc.), you'll think Julia is slow and laggy because loading packages and running stuff has a noticeable wait time (10-40sec). It's worth mentioning that this lag is "compilation overhead". What this means is, Julia tries to pre-compile as much code as possible (into optimized machine code) when running something or loading a package. This lag exists **only the first time** you run something. Every subsequent run of a function, even with different parameters, will be **substantially** faster, and in most cases instant. If you want to test this yourself, try to run a line of code twice with `@time` before the function and compare the results. Here's an example:
``` 
julia> @time using PopGen
 11.870878 seconds (27.83 M allocations: 1.389 GiB, 5.26% gc time)

julia> @time using PopGen
 0.000272 seconds (406 allocations: 21.375 KiB)
```

## Semicolons

Semicolons will come up a lot in Julia, probably more than you would expect if you are migrating from another language.  They mean different things depending on where they are.

### At the end of a command

 When you see a semicolon after invoking a function, what that means is "don't show me the output".

Example: 

```julia
julia> x = 2 + 2
4

julia> x = 10 + 2 ;

julia> x
12
```

Julia will still process the command and assign `10 + 2` to `x`, but it won't show you the output. We sometimes include a semicolon after commands in these docs to mimic what the REPL output would look like without spitting back large volumes of text. **These semicolons are optional** 



### In between assignment commands

If you see a semicolon in between two variable assignments or commands, like so:

```julia
julia> x = [1,2] ; y = [3,4]
```

that's a Julia short-hand for making two short lines of code appear on a single line. It's the equivalent of doing:

```julia
julia> x = [1,2]
julia> y = [3,4]
```

We sometimes choose this writing format for very quick  and small assignments hoping to save some visual space. Use whichever method is most comfortable and sensible for you!

## Help mode

To enter `help` mode in the REPL, simply press the question mark key `?` (shift + key) and you will notice a different prompt `help?>` for you to type in a function.

```julia
help?>population
```

```
search: population populations population! populations!

  population(data::PopData; listall::Bool = false)

View unique population ID's and their counts in a `PopData`.


- `listall = true` displays all samples and their `population` instead (default = `false`)
```

## Type information

Julia encourages strong typing of variables, and the functions in `PopGen` are no exception to this. However, to reduce the barrier of entry required to understand this documentation and the subsequent package, we have chosen to omit some of the `type` information from functions to reduce visual clutter for newer users. As experienced users already know, if you would like to see the explicit type information, you can look at the code on github, invoke the `help` system in the REPL (above), or search for a function in the Documentation pane in Juno. 

You'll notice types follow a specific format, which is `object::type`. This format is a type declaration, so in the function `population`, which looks like: 

```julia
population(data::PopData; listall::Bool = false)
```

- `data` is a variable of type `PopData` 
- `listall` is a variable of type `Bool` (Boolean) meaning it only takes `true` or `false` without quotes, and the default value is set to `false`

### Type Unions

You might see the type `Union` appear occasionally throughout this documentation, and you can consider it a list of allowable types. For example, if something was of type `::Union{String,Integer}`, that means that **either** a `String` or `Integer` works. 

### Subtypes

The julia language is abound with types (and you can create your own!), and has a hierarchical system of supertypes and subtypes. As you can probably guess, a supertype can contain multiple subtypes, such as `Signed` being a supertype of (among other things) `Int64, Int32, Int16, Int8`. All vectors are subtypes of `AbstractVector`.  If you want to try it yourself, use the `supertype()` command on your favorite Type, like `supertype(Float32)`. You will occasionally see `<:` instead of `::`, which means "is a subtype of". This is used for condtional evaluation, like `typeof(something) <: Signed`, and in some function methods like `function(var1::T) where T <: Supertype`, which leads us to:

#### where T

This looks weird at first,  but it's actually very simple. When we do method definitions, we can define methods with strict types, like `funct(data:PopData, arg1::Int8)`, or we can generalize it with `where T`, which looks like :

```julia
function funct1(data::PopData, thing1::T) where T
```

This will auto-create a method for any possible Type for `thing1`. That's really convenvient, but sometimes it's problematic, as incorrect input can lead to obscure errors (e.g. multiplying integers with strings?!). Instead, you can constrain the types for `T` like this:

```julia
function funct2(data::PopData, thing1::T) where T <: Signed
```

With the constraint above, it will generate methods for all cases where `thing1` is a subtype of `Signed`, which includes all the numerical Types (integers, Floats, etc.). This will make sure that the function will behave correctly for a range of input types.

You can also use this type of notation to clean up a method definition where multiple arguments have the same Type specification:

```julia
function funct3(data::PopData, thing1::T, thing2::T, thing3::T) where T<:AbstractFloat
```

So, instead of writing `thing1::Float64, thing2::Float64, thing3::Float64`, we just use `T` as a placeholder and assign it as a subtype of something at the end.  It ends up being pretty handy!


## Functions vs. Methods 

As part of Julia's type-safe paradigm and multiple dispatch (see "ERROR: MethodError: no method matching" below), type specifications in functions often reduce runtime of functions, but also establish function identity. Multiple dispatch refers to several different functions having the same name, but employing different *methods* depending on the input. In Julia, it's easier to write a single function with multiple type-safe methods, rather than one mega-function that accepts any type and have a bunch of `if` statements that determines what the program does depending on the input. 

:::info Best Practice
As a rule of thumb, `for`  loops with `if` conditions in them slow down the compiler, so best-practice often encourages us to write type-specific methods.
:::
In practice, this looks like:

```
# combine two numbers
julia> function add(x::Integer, y::Integer)
           x+y
       end
add (generic function with 1 method)

# combine two strings
julia> function add(x::String, y::String)
           x*y
       end
add (generic function with 2 methods)
                
julia> add(1,2)
3

julia> add("water", "melon")
"watermelon"
```

Multiple dispatch therefor leads to a unique type of possible error: the `MethodError`

### ERROR: MethodError: no method matching

Using the function `add` from the example above, let's have a look at what happens when we try to `add` an `Integer` with a `String`:

```julia
julia> add(1,"melon")
ERROR: MethodError: no method matching add(::Int64, ::String)
Closest candidates are:
  add(::String, ::String) at none:2
  add(::Integer, ::Integer) at none:2
Stacktrace:
 [1] top-level scope at none:0
```

This error is telling us "there is no such function called `add`, who's inputs are an `Integer` followed by a `String`". But, it does offer us some alternatives, like the two `add` functions we created earlier.

The functions within `PopGen` are almost always explicitly typed, so if you are getting the `MethodError: no method matching` error, then you are inputting the incorrect types into the function, or perhaps your inputs for the arguments are in the wrong order (see "Functions with and without keywords" below).

Sometimes you might include an argument with a keyword when there isn't one, or include an argument without a keyword when there needs to be one (honestly, we make that mistake too and we *wrote* this stuff). To help minimize those mistakes, please read below about which arguments have keywords and which don't.

:::note MethodErrors
MethodError's can definitely get annoying, but they are more commonly the result of incorrect inputs versus being bugs. If you double-checked your inputs and things still don't work, please submit an issue. Thanks!
:::


## Functions with and without keywords 

Let's talk about semicolons some more.

:::info TL;DR
Reading these docs, pay attention to semicolons in the function argument lists.

-  arguments before a semicolon have no keyword and follow an explicit order
-  arguments after a semicolon have a keyword `argument = value` and their order doesn't matter
-  `MethodError: no methods matching` is often a user error and not a bug, but if it is, please open an issue!

:::

Broadly speaking, there are two types of function declarations in Julia: ones with keywords and ones without keywords. The term "keywords" refers to an input argument that has the format `argument = value`. This format is present in many of the functions in this and other packages, however there are some specifics to understand when functions use keywords and when they don't. 

<Tabs
  block={true}
  defaultValue="1"
  values={[
    { label: '1. No semicolon in argument list', value: '1', },
    { label: '2. Semicolon in argument list', value: '2', },
  ]
}>
<TabItem value="1">

```julia
function function_name(var1::type, var2::type, var3::type)
    do stuff with vars
end
```
If a function is declared with only commas in the argument list, like shown above, then the arguments to that function **must** have no keywords and follow the exact order they appear in. If the generic example above had the typing:

```julia
function function_name(var1::String, var2::Float64, var3::Array{String,1})
    do stuff with vars
end
```
then the only acceptable way to run this function without getting a `MethodError` would be with arguments in the order of `function_name(String, Float64, Array{String,1})`. Even if some of the arguments have a default values, like `var2::Float64 = 6.66`, the order of arguments/types has to be respected as declared.

</TabItem>
<TabItem value="2">

```julia
function function_name(var1::type; var2::type, var3::type)
    do stuff with vars
end
```

In this format, everything that comes **before** the semicolon follows the strict rules from **Format 1**, and everything that comes **after** the semicolon is a keyword argument. Keyword arguments have the flexibility to not require any particular input order. However, you **must** use the keywords to declare those arguments, or you will receive another `MethodError: no method matching`, which is, as we've mentioned, annoying. 

</TabItem>
</Tabs>
