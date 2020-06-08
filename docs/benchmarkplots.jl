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

pop_adeg = ["PopGen.jl", "adegenet"]
pop_hierf = ["PopGen.jl", "hierfstat"]

#### Load in Data ####
import_speed = [1.445, 6.745]
speedplot = comparison_plot(pop_adeg, import_speed, "Seconds", "Importing a genepop file")
speedplot |> save("speedplot.png")

#### Filesize (KB) ####
obj = ["PopData (PopGen.jl)", "genind (adegenet)"]
f_size = [3.498172, 5.331536]
objplot = comparison_plot(obj, f_size, "megabytes", "Data structure size")
objplot |> save("objectplot.png")


#### f-stat summary ####
sumstat = [0.242, 4.6]
sumstatplot = comparison_plot(pop_hierf, sumstat, "Seconds", "Summary Statistics")
sumstatplot |> save("sumstatplot.png")

#### Χ² test ####
chitest = [0.591396, 6.2659]
chiplot = comparison_plot(pop_adeg, chitest, "Seconds", "Hardy-Weinberg Equilibrium Χ² test")
chiplot |> save("chisqplot.png")


#### Makie version
using Makie
# Made in Juno, so preferring the Plot pane
popdisplay(AbstractPlotting.PlotDisplay())
AbstractPlotting.inline!(true)

# set generic X axis
xaxis = ["Julia", "R"]

# create generic plotting function
function comparison_plot(x::Vector{String},y::Vector{Float64}, yaxis::String)
    scene = barplot(
        x,
        y,
        color = ["#aa79c1","#769fd2"]
    )
    axis= scene[Axis]
    axis[:names][:axisnames] = ("", yaxis)
    axis[:grid][:linewidth] = (0, 0)
    axis[:ticks][:linewidth] = (0,0)
    axis[:frame][:frames] = ((true,false),(true,false))
    return scene
end
