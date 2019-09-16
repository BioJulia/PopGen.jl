# A quick Julia primer for getting the most out of this documentation

There is nothing inherently special about this documentation relative to other documentation, other than we really *really* want you to get the most out of what's written here. This means that we need to embrace the fact that both novice and experienced Julia users may be reading these docs and using this package. So let's cover some Julia basics that will really help in navigating this package before we even get into the complicated genetics stuff. This primer is by no means "everything you need to get started in Julia", and is a poor substitute for actually learning the language. In general, we recommend [Think Julia: How to Think Like a Data Scientist](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html) by Ben Lauwens to establish some solid Julia foundations. It's free online!



## Semicolons

Semicolons will come up a lot in Julia, probably more than you would expect if you are migrating from another language.  They mean different thing depending on where they are.

### At the end of a command

 When you see a semicolon after invoking a function, what that means is "don't show me the output in the terminal window".

Example: 

```julia
julia> x = 2 + 2
4

julia> x = 10 + 2 ;

julia> x
12
```

Julia will still process the command and assign `10 + 2` to `x`, but it won't show you the output in the terminal. We sometimes include a semicolon after commands in these docs to mimic what the REPL output would look like without spitting back out an array of over 200 values. **These semicolons are optional** 



### In between assignment commands

If you see a semicolon in between two variable assignments, like so:

```julia
x = [1,2] ; y = [3,4]
```

that's a Julia short-hand for making two short lines of code appear on a single line. It's the equivalent of doing:

```julia
x = [1,2]
y = [3,4]
```

We sometimes choose this writing format for very quick  and small assignments hoping to save some visual space. Use whichever method is most comfortable and sensible for you!

## Help mode

To enter `help` mode in the REPL, simply press the question mark key `?` (shift + key) and you will notice a different prompt `help?>` for you to type in a function.

```julia
help?>popid
```

```
search: popid popid! popfirst! popdisplay precompile __precompile__ CompositeException

  popid(x::PopObj; listall::Bool = false)

  View unique population ID's in a PopObj.

  listall = true, displays ind and their popid instead (default = false).
```



## Type information

Julia encourages strong typing of variables, and the functions in `PopGen` are no exception to this. However, to reduce the barrier of entry required to understand this documentation and the subsequent package, we have chosen to omit some of the `type` information from functions to reduce visual clutter for newer users. As experienced users already know, if you would like to see the explicit type information, you can look at the code on github, invoke the `help` system in the REPL (above), or search for a function in the Documentation pane in Juno. 

You'll notice types follow a specific format, which is `object::type`. This format is a type assignment, so in the function `popid`, which looks like: `popid(x::PopObj; listall::Bool = false)`:

- `x` is a variable of type `PopObj` 
- `listall` is a variable of type `Bool` (boolean) meaning it only takes `true` or `false` without quotes

### Type Unions

You might see the type `Union` appear occasionally throughout this documentation, and you can consider it a list of types. For example, if something was of type `::Union{String,Integer}`, that means that **either** a `String` or `Integer` works. These kinds of `Union` appear in a few places, like `remove_inds!` , where the input can be either a `String` (a single individual) or an `Array{String,1}` (one-dimensional list of names). For `remove_inds!`, the second argument type appears as `::Union{String, Array{String,1}}`. Also, the order in which types appear in `Union` types don't matter.



## Functions vs. Methods 

As part of Julia's type-safe paradigm and multiple dispatch (see "ERROR: MethodError: no method matching" below), type specifications in functions often reduce runtime of functions, but also establish function identity. Multiple dispatch refers to several different functions having the same name, but employing different *methods* depending on the input. In Julia, it's easier to write a single function with multiple type-safe methods, rather than one mega-function that accepts any type and have a bunch of `if` statements that determines what the program does depending on the input. 

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

!!! Note "keep in mind"
    MethodError's can definitely get annoying, but they are usually the result of incorrect input from the user and not buggy programming by the developers. Please take that into consideration before assuming something is broken or bugged.



## Functions with and without keywords 

!!! info "TL;DR"
    Reading these docs, pay attention to semicolons in the function argument lists.
    

    -  arguments before a semicolon have no keyword and follow an explicit order
    -  arguments after a semicolon have a keyword `argument = value` and their order doesn't matter
    - `MethodError: no methods matching` is more likely an issue on your side and not on our side :smile:
        - unless we accidentally forgot to export a function! :facepalm:

Broadly speaking, there are two types of function declarations in Julia: ones with keywords and ones without keywords. The term "keywords" refers to an input argument that has the format `argument = value`. This format is present in many of the functions in this and other packages, however there are some specifics to understand when functions use keywords and when they don't. 

### Format 1: Strict argument order and no keywords - No semicolon in argument list

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

### Format 2 - semicolon in argument list

```julia
function function_name(var1::type; var2::type, var3::type)
    do stuff with vars
end
```

In this format, everything that comes **before** the semicolon follows the strict rules from **Format 1**, and everything that comes **after** the semicolon is a keyword argument. Keyword arguments have the flexibility to not require any particular input order. However, you **must** use the keywords to declare those arguments, or you will receive another `MethodError: no method matching`, which is, as we've mentioned, annoying. 

