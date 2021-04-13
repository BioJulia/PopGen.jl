using VegaLite, DataFrames

function comparison_plot(x::Vector{String},y::Vector{Float64}, yaxis::String, main::String)
    corners = 10
    @vlplot(
        height=300,
        width=650,
        title={
          text=main,
          fontSize=20,
          fontWeight="normal"
        },
        mark={
            :bar,
            cornerRadiusTopLeft=0,
            cornerRadiusTopRight=corners,
            cornerRadiusBottomLeft=0,
            cornerRadiusBottomRight=corners
        },
        x={y,
            axis={
                title=yaxis,
                titleFontSize = 17,
                titleFontWeight = "normal",
                labelFontSize = 12,
                grid = false,
                domain = false
            }
        },
        y={x,
            axis={
                title="",
                labelAngle= 0,
                labelFontSize = 17,
                domain = false,
                ticks = false,
                labelPadding = 4
                }
            },
        color={
            x,
            scale={range=["#aa79c1","#769fd2"]},
            legend=false
        }
    )
end

pop_adeg = ["PopGen.jl", "adegenet"]
pop_hierf = ["PopGen.jl", "hierfstat"]
pop_p_hierf = ["PopGen.jl", "PopGen.jl (6 threads)", "hierfstat"]

#### Load in Data ####
import_speed = [.910, 6.745]
comparison_plot(pop_adeg, import_speed, "Seconds", "Importing a genepop file (less is better)") |> save("static/img/speedplot.png")

#### Filesize (mb) ####
obj = ["PopData\n(PopGen.jl)", "genind\n(adegenet)"]
f_size = [3.498172, 5.331536]
# I know, it's kind of obnoxiously long, but this is the syntax for getting that line on there.
DataFrame(:names => obj, :size => f_size) |>
    @vlplot() +
    @vlplot(
        height=300,
        width=650,
        title={
          text="Data structure size relative to source file",
          fontSize=22,
          fontWeight="normal"
        },
        mark={
            :bar,
            cornerRadiusTopLeft=0,
            cornerRadiusTopRight=15,
            cornerRadiusBottomLeft=0,
            cornerRadiusBottomRight=15
        },
        x={"size:q",
            axis={
                title="megabytes",
                titleFontSize = 17,
                titleFontWeight = "normal",
                labelFontSize = 12,
                grid = false,
                domain = false
            }
        },
        y={"names:n",
            axis={
                title="",
                labelAngle= 0,
                labelFontSize = 17,
                domain = false,
                ticks = false,
                labelPadding = 4
                }
            },
        color={
            "names:n",
            scale={range=["#aa79c1","#769fd2"]},
            legend=false
        }
    ) +
    @vlplot(:rule, x={datum=3.2})  |> save("static/img/objectplot.png")

#### f-stat summary ####
sumstat = [0.171, 4.6]
comparison_plot(pop_hierf, sumstat, "Seconds", "Summary Statistics (less is better)") |> save("static/img/sumstatplot.png")

#### Χ² test ####
chitest = [0.176, 6.2659]
comparison_plot(pop_adeg, chitest, "Seconds", "Hardy-Weinberg Equilibrium Χ² test (less is better)") |> save("static/img/chisqplot.png")



function comparison_plot_3x(x::Vector{String},y::Vector{Float64}, yaxis::String, main::String)
    corners = 10
    @vlplot(
        height=300,
        width=650,
        title={
          text=main,
          fontSize=20,
          fontWeight="normal"
        },
        mark={
            :bar,
            cornerRadiusTopLeft=0,
            cornerRadiusTopRight=corners,
            cornerRadiusBottomLeft=0,
            cornerRadiusBottomRight=corners
        },
        x={y,
            axis={
                title=yaxis,
                titleFontSize = 17,
                titleFontWeight = "normal",
                labelFontSize = 12,
                grid = false,
                domain = false
            }
        },
        y={x,
            axis={
                title="",
                labelAngle= 0,
                labelFontSize = 17,
                domain = false,
                ticks = false,
                labelPadding = 4
                }
            },
        color={
            x,
            scale={range=["#aa79c1","#aa79c1","#769fd2"]},
            legend=false
        }
    )
end

fst_speed = [.580, 0.299, 199.48]
tt = comparison_plot_3x(pop_p_hierf, fst_speed, "Seconds", "Weir & Cockerham Pairwise FST (less is better)")
save("static/img/fstplot.png", tt)