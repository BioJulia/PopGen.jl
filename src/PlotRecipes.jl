@userplot missing_plot

@recipe function f(h::missing_plot)
    if length(h.args) != 1 || !(typeof(h.args[1]) <: PopObj)
        error("Missingness Plot needs to be given a single PopObj")
    end
    x, y = missing(h.args[1])

    # set up the subplots
    legend := false
    framestyle := :zerolines
    xforeground_color_axis := :transparent
    yforeground_color_axis := :transparent
    xgrid := false
    ytickfontsize := 10

    group := x.population,
    color_palette := [HSL(i, .50, .50) for i in range(0, stop=300, length=unique(x.population) |> length)]
    hover := x.name .* ", " .* ["$i" for i in x.missing]
    title := "Number of missing loci per population"
    markerstrokealpha := 0.2
    markersize = 5.5

    # missing by individual
    @series begin
        seriestype := :dotplot
        group := x.population
        subplot := 1
        x.missing
    end

    # missing by locus
    @series begin
        seriestype := :histogram
        bins := length(unique(:missing))
        title := "Distribution of locus missingness"
        xguide := "# times missing"
        yguide := "# loci"
        linealpha := 0.1
        fill := "rgb(217, 217, 217)"
        xguidefontsize := 15
        xtickfontsize := 10
        ytickfontsize := 10
        yguidefontsize := 15
        subplot := 2
        y.missing
    end
end

colors = [HSL(i, .50, .50) for i in range(0, stop=300, length=unique(sh_miss[1].population) |> length)]

p1 = @df sh_miss[1] dotplot(
                :population,
                :missing,
                group = :population,
                color_palette = [HSL(i, .50, .50) for i in range(0, stop=300, length=unique(:population) |> length)],
                hover = :name .* ", " .* ["$i" for i in :missing],
                title = "Number of missing loci per population",
                xgrid = false,
                framestyle = :zerolines,
                legend = false,
                #xtickfontsize = 10,
                xforeground_color_axis = :transparent,
                ytickfontsize = 10,
                yforeground_color_axis = :transparent,
                markerstrokealpha = 0.2,
                markersize = 5.5,
)


p2 = @df sh_miss[2] histogram(
        :missing,
        bins = length(unique(:missing)),
        legend = false,
        title = "Distribution of locus missingness",
        xguide = "# times missing",
        yguide = "# loci",
        framestyle = :zerolines,
        linealpha = 0.1,
        fill = "rgb(217, 217, 217)",
        #hover = [count(i -> i ==j, :missing) for j in sort(unique(:missing))],
        xgrid = false,
        xforeground_color_axis = :transparent,
        xguidefontsize = 15,
        xtickfontsize = 10,
        yforeground_color_axis = :transparent,
        ytickfontsize = 10,
        yguidefontsize = 15,
)

plot(p1, p2, layout = (1,2))
