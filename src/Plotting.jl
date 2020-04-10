"""
    plot_missing(data::PopObj)
Return an interactive plot of the number of missing loci in individuals of a
`PopObj`, along with the number of missing individuals per locus.
"""
function plot_missing(data::PopObj)
    by_sample,by_loci = missing(data);
    #ys = Array[subdf[!, :missing] for subdf in groupby(by_sample[!, 1:3], :population)]
    #texts = Array[subdf[!, :name] for subdf in groupby(by_sample[!, 1:3], :population)]
    popnum = length(by_sample.population |> unique)
    colors = ["hsl($i, 50%, 50%)" for i in range(0, stop=300, length=popnum)]

    byind = box(
        by_sample,
        x = :missing,
        group = :population,
        marker_color=colors,
        text = :name,
        jitter = 1,
        pointpos = 0,
        marker_size = 8.5,
        boxpoints = "all",
        marker=attr(line=attr(width=0.75))
    )

    layout_ind = Layout(
        title = "Number of missing loci per population",
        hovermode = "closest",
        yaxis = attr(title = "# missing loci", zeroline = false),
        xaxis = attr(title = "Population", zeroline = false)
    )


    loci_hist = histogram(x = by_loci[!, :missing],
                          marker_color = "rgb(217, 217, 217)",
                          name = "",
                          text = "loci",
                          xbins_size = 1,
                          showlegend = false)



    layout_hist = Layout(title = "Distribution of Locus Missingness",
                         xaxis = attr(title = "# of times locus missing", zeroline = false),
                         yaxis = attr(title = "# of loci", zeroline = false)
                         )
    ind_plot =  plot(byind, layout_ind)
    loci_plot = plot(loci_hist, layout_hist)

    return [ind_plot loci_plot]
end


"""
    plot_locations(x::PopObj; region::String = "world", projection::String = "orthographic")

Returns a simple lower resolution interactive scatterplot of the individuals in a
`PopObj`. Default `region` and `projection` are "world" and "orthographic",
respectively. If a specific region is set, the plot will default to "mercator"
projection unless `projection =` is used to specify a different one.

Example:\n
`plot_locations(manatees, region = "usa", projection = "albers usa")`

[regions]\n
"usa","europe", "asia", "africa", "north america", "south america"

[projections]\n
"equirectangular", "mercator", "orthographic", "natural earth",
"kavrayskiy7", "miller", "robinson", "eckert4",
"azimuthal equal area", "azimuthal equidistant",
"conic equal area", "conic conformal" , "conic equidistant",
"gnomonic", "stereographic", "mollweide", "hammer",
"transverse mercator", "albers usa", "winkel tripel",
"aitoff", "sinusoidal"
"""
function plot_locations(x::PopObj; region::String = "world", projection::String = "orthographic")
    # test for missing?
    if projection ∉ ["equirectangular", "mercator", "orthographic", "natural earth",
                    "kavrayskiy7", "miller", "robinson", "eckert4",
                    "azimuthal equal area", "azimuthal equidistant",
                    "conic equal area", "conic conformal" , "conic equidistant",
                    "gnomonic", "stereographic", "mollweide", "hammer",
                    "transverse mercator", "albers usa", "winkel tripel",
                    "aitoff", "sinusoidal"]
        error("Projection not recognized. Please see the help doc for list of projection options")
    end
    y = PopOpt(x)
    popnum = length(y.samples.population |> unique)
    colors = ["hsl($i, 50%, 50%)" for i in range(0, stop=300, length=popnum)]
    df_split = groupby(y.samples, :population)
    map_scatter = [scattergeo(lat=df_split[i][!, :latitude],
                       lon=df_split[i][!, :longitude],
                       marker_line_color="rgb(62,90,112)", marker_line_width=1,
                       marker_color = colors[i],
                       name = df_split[i][!, :population][1],
                       text = df_split[i][!, :name]
                       ) for i in 1:popnum]
   if region == "world"
       geo = attr(scope = region,
                  projection_type = projection,
                  showcoastlines = false,
                  showcountries = true,
                  countrywidth = 0.75,
                  countrycolor = "rgb(255,255,255)",
                  subunitcolor = "rgb(255,255,255)",
                  showland = true,
                  landcolor = "rgb(217, 217, 217)",
                  )
   elseif lowercase(region) ∈ ["usa","europe", "asia", "africa", "north america", "south america"]
      if projection == "orthographic"
          proj = "mercator"
      else
          proj = projection
      end
      geo = attr(scope = lowercase(region),
                 showlakes = true,
                 lakecolor = "#a5a5a5",
                 showrivers = true,
                 rivercolor = "#a5a5a5",
                 showsubunits = true,
                 projection_type = proj,
                 showcoastlines = false,
                 showcountries = true,
                 countrywidth = 0.75,
                 countrycolor = "rgb(255,255,255)",
                 subunitcolor = "rgb(255,255,255)",
                 showland = true,
                 landcolor = "rgb(217, 217, 217)",
                 )
   else
       error("Leave \"region =\" empty or use one of: usa, europe, asia, africa, north america, south america")
   end
    layout = Layout(;title="Sample locations", showlegend=true, geo=geo)
    plot(map_scatter, layout)
end
