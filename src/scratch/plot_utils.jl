using CairoMakie
using GLMakie
using PopGen

cats = @nancycats
sharks = @gulfsharks

function plot_missing_data(data::PopData; by::String = "all")
    f = Figure()
    
    if by == "all"
      # for loci
      missing_loc = missing_data(data, by = "locus")
      Axis(f[1, 1], title = "Number of Missing Per Locus")
      hist!(f[1, 1], missing_loc.missing, color = "#769fd2")

      # for loci
      missing_samples = missing_data(data, by = "sample")
      Axis(f[2, 1], title = "Number of Missing Per Sample")
      hist!(f[2, 1], missing_samples.missing, color = "#aa79c1")
      
      # population
      missing_pop = sort(missing_data(data, by = "population"), :missing, rev = true)
      Axis(f[:, 2], title = "Missing Data Per Population", xticks = (1:length(missing_pop.population), missing_pop.population), xticklabelrotation = π/3)
      stem!(f[:, 2], missing_pop.population, missing_pop.missing, color = "#d5635c", strokewidth = 0, markersize = 12, trunkwidth = 0)
      
      Label(f[1:2, 0], "Count", rotation = pi/2)
    
    elseif by ∈ ["locus", "loci"]
      missing_dat = missing_data(data, by = "locus")
      Axis(f[1, 1], title = "Missing Data Per Locus", ylabel = "number of loci", xlabel = "number of missing genotypes")
      hist!(f[1, 1], missing_dat.missing, color ="#769fd2")
    
    elseif by ∈ ["sample", "samples"]
      missing_dat = missing_data(data, by = "sample")
      Axis(f[1, 1], title = "Missing Data Per Sample", ylabel = "number of samples", xlabel = "number of missing loci")
      hist!(f[1, 1], missing_dat.missing, color = "#aa79c1")
    
    elseif by ∈ ["population", "populations"]
      missing_pop = sort(missing_data(data, by = "population"), :missing, rev = true)
      Axis(f[1, 1], title = "Missing Genotypes Per Population", xticks = (1:length(missing_pop.population), missing_pop.population), xticklabelrotation = π/3, ylabel = "Number of missing genotypes")
      stem!(f[1, 1], missing_pop.population, missing_pop.missing, color = "#d5635c", strokewidth = 0, markersize = 12, trunkwidth = 0)
    
    else
      throw(ArgumentError("Please choose a plotting method from \"all\", \"sample\", \"locus\", or \"population\""))
    
    end
  f
end

x = summary_stats(cats)

function plot_summary(data::PopData)
  f = Figure()
  summ = summary_stats(data)
  vals = [i[1] for i in eachcol(summ)]
  stats = names(summ)
  Axis(f[1, 1], title = "Summary Statistics", xticks = (1:length(stats), stats), ylabel = "value")
  j = stem!(f[1, 1], stats, vals, color = "#aa79c1", strokewidth = 0, markersize = 12, trunkwidth = 0) 
  return j
end


function plot_sumstats(data::PopData)
  fig = Figure()
  loci_names = loci(data)
  summ = summary_stats(data, by = "locus")
  stats = names(summ)[2:end]
  menu = Menu(fig, options = loci_names)
  fig[1, 1] = vgrid!(
      Label(fig, "Locus", width = nothing),
      menu,
      tellheight = false, 
      width = 200
  )
  ax = Axis(fig[1, 2], xticks = (1:length(stats), stats))
  ylims!(ax,-0.05, 1.0)
  locs = Node{Any}(zeros(length(stats)))
  ys = @lift(collect($locs))
  stem!(ax, stats, ys, color = "#aa79c1", strokewidth = 0, markersize = 12, trunkwidth = 0) 
  
  on(menu.selection) do s
    _l = filter(:locus => x -> x == s, summ)
    locs[] = [i[1] for i in eachcol(_l)][2:end]
  end
  fig
end