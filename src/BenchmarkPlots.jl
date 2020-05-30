using VegaLite

function comparison_plot(x::Vector{String},y::Vector{Float64}, yaxis::String, main::String)
    corners = 25
    @vlplot(
        height=500,
        width=450,
        title={
          text=main,
          fontSize=22,
          fontWeight="normal"
        },
        mark={
            :bar,
            cornerRadiusTopLeft=corners,
            cornerRadiusTopRight=corners,
            cornerRadiusBottomLeft=corners,
            cornerRadiusBottomRight=corners
        },
        y={y,
            axis={
                title=yaxis,
                titleFontSize = 17,
                titleFontWeight = "normal",
                labelFontSize = 12,
                grid = false,
                domain = false
            }
        },
        x={x,
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

xaxis = ["Julia", "R"]

#### Load in Data ####
import_speed = [0.910, 6.745]
speedplot = comparison_plot(xaxis,import_speed, "Seconds", "Importing a genepop file")
speedplot |> save("speedplot.png")

#### Filesize (KB) ####
obj = ["PopData (Julia)", "genind (R)"]
f_size = [3.527765, 5.331536]
objplot = comparison_plot(obj, f_size, "megabytes", "Data structure size")
objplot |> save("objectplot.png")

#### Χ² test ####
chitest = [0.176, 6.2659]
chiplot = comparison_plot(xaxis, chitest, "Seconds", "Hardy-Weinberg Equilibrium Χ² test")
chiplot |> save("chisqplot.png")
