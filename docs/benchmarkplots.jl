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
import_speed = [.910, 6.745]
comparison_plot(pop_adeg, import_speed, "Seconds", "Importing a genepop file") |> save("static/img/speedplot.png")

#### Filesize (mb) ####
obj = ["PopData (PopGen.jl)", "genind (adegenet)"]
f_size = [3.498172, 5.331536]
# I know, it's kind of obnoxiously long, but this is the syntax for getting that line on there.
DataFrame(:names => obj, :size => f_size) |>
    @vlplot() +
    @vlplot(
        height=500,
        width=450,
        title={
          text="Data structure size",
          fontSize=22,
          fontWeight="normal"
        },
        mark={
            :bar,
            cornerRadiusTopLeft=25,
            cornerRadiusTopRight=25,
            cornerRadiusBottomLeft=25,
            cornerRadiusBottomRight=25
        },
        y={"size:q",
            axis={
                title="megabytes",
                titleFontSize = 17,
                titleFontWeight = "normal",
                labelFontSize = 12,
                grid = false,
                domain = false
            }
        },
        x={"names:n",
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
    @vlplot(:rule, y={datum=3.2})  |> save("static/img/objectplot.png")

#### f-stat summary ####
sumstat = [0.171, 4.6]
comparison_plot(pop_hierf, sumstat, "Seconds", "Summary Statistics") |> save("static/img/sumstatplot.png")

#### Χ² test ####
chitest = [0.176, 6.2659]
comparison_plot(pop_adeg, chitest, "Seconds", "Hardy-Weinberg Equilibrium Χ² test") |> save("static/img/chisqplot.png")


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
