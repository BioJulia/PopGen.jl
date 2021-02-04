fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
anyerrors = false

using PopGen
using Test

all_tests = [
    "types.jl",
    "allelefreq.jl",
    "io.jl",
    "manipulate.jl",
    "dataexploration.jl",
    "summaryinfo.jl"
]

println("Running tests:")

for a_test in all_tests
    try
        include(a_test)
        println("\t\033[1m\033[32mPASSED\033[0m: $(a_test)")
    catch e
        global anyerrors = true
        println("\t\033[1m\033[31mFAILED\033[0m: $(a_test)")
        if fatalerrors
            rethrow(e)
        elseif !quiet
            showerror(stdout, e, backtrace())
            println()
        end
    end
end

if anyerrors
    throw("Tests failed :(")
end