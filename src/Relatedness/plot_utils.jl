using VegaLite

function relate_heat(data::DataFrame, method::Symbol)
    mth = String(method)
    data2 = select(data, :sample_2 => :sample_1, :sample_1 => :sample_2, :n_loci, :Lynch)
    append!(data2, data)
    sort!(data2, [:sample_1, :sample_2])
    transform!(data2, method => (i -> round.(i, digits = 4)) => method )
    tt_mthd = mth * ": '+datum."*mth
    data2 |> @vlplot(
        transform = [
            {calculate="''+datum.sample_1+ ' '+' '+datum.sample_2+'  '+' $tt_mthd", as="tt"}
        ],
        :rect, 
        tooltip = "tt:n",
        :sample_1, 
        :sample_2, 
        color=method
    )
end

function relate_heat2(data::DataFrame, method::Symbol)
    mth = String(method)
    data2 = select(data, :sample_2 => :sample_1, :sample_1 => :sample_2, :n_loci, :Lynch)
    append!(data2, data)
    sort!(data2, [:sample_1, :sample_2])
    transform!(data2, method => (i -> round.(i, digits = 4)) => method )
    data2 |> @vlplot(
        hconcat=[
          {
            height = 700,
            width = :container,
            selection={
              brush={
                type="interval"
              }
            },
            encoding={
              x={
                field="sample_1",
                type="nominal"
              },
              color={
                field=mth
              },
              y={
                field="sample_2",
                type="nominal"
              }
            },
            mark="rect"
          },
          {
            hconcat=[
              {
                encoding={
                  text={
                    field="sample_1",
                    type="nominal"
                  },
                  y={
                    axis=nothing,
                    field="row_number",
                    type="ordinal"
                  }
                },
                title="Sample 1",
                width=50,
                mark="text"
              },
              {
                encoding={
                  text={
                    field="sample_2",
                    type="nominal"
                  },
                  y={
                    axis=nothing,
                    field="row_number",
                    type="ordinal"
                  }
                },
                title="Sample 2",
                width=50,
                mark="text"
              },
              {
                encoding={
                  text={
                    field=mth,
                    type="nominal"
                  },
                  y={
                    axis=nothing,
                    field="row_number",
                    type="ordinal"
                  }
                },
                title="Est",
                width=50,
                mark="text"
              }
            ],
            transform=[
              {
                filter={
                  selection="brush"
                }
              },
              {
                window=[
                  {
                    as="rank",
                    op="rank"
                  }
                ]
              },
              {
                filter={
                  field="rank",
                  lt=20
                }
              }
            ]
          }
        ],
        width=300,
        resolve={
          legend={
            color="independent"
          }
        },
        transform=[
          {
            window=[
              {
                as="row_number",
                op="row_number"
              }
            ]
          }
        ]
      )
end


c = relate_heat2(a, :Lynch)

function relate_hist(data::DataFrame, method::Symbol)
    mth = String(method)
    a |> @vlplot(
        width = :container, 
        height = 350, 
        :bar, #tooltip = "true",
        x={method, bin=true, step = 5, title = mth}, 
        y={"count()", title = "Count"}
    )
end


@vlplot(
  height=300,
  hconcat=[
    {
      selection={
        brush={
          type="interval"
        }
      },
      encoding={
        x={
          field="sample_1",
          type="nominal"
        },
        color={
          field="Lynch"
        },
        y={
          field="sample_2",
          type="nominal"
        }
      },
      mark="rect"
    },
    {
      hconcat=[
        {
          encoding={
            text={
              field="sample_1",
              type="nominal"
            },
            y={
              axis=nothing,
              field="row_number",
              type="ordinal"
            }
          },
          title="Sample 1",
          width=50,
          mark="text"
        },
        {
          encoding={
            text={
              field="sample_2",
              type="nominal"
            },
            y={
              axis=nothing,
              field="row_number",
              type="ordinal"
            }
          },
          title="Sample 2",
          width=50,
          mark="text"
        },
        {
          encoding={
            text={
              field="Lynch",
              type="nominal"
            },
            y={
              axis=nothing,
              field="row_number",
              type="ordinal"
            }
          },
          title="Est",
          width=50,
          mark="text"
        }
      ],
      transform=[
        {
          filter={
            selection="brush"
          }
        },
        {
          window=[
            {
              as="rank",
              op="rank"
            }
          ]
        },
        {
          filter={
            field="rank",
            lt=20
          }
        }
      ]
    }
  ],
  width=300,
  resolve={
    legend={
      color="independent"
    }
  },
  transform=[
    {
      window=[
        {
          as="row_number",
          op="row_number"
        }
      ]
    }
  ]
)