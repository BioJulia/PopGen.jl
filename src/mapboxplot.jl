## from SGLyons. Try to implement it with gulfsharks to compare to current implementation

mapbox_access_token = "pk.eyJ1IjoiY2hlbHNlYXBsb3RseSIsImEiOiJjaXFqeXVzdDkwMHFrZnRtOGtlMGtwcGs4In0.SLidkdBMEap9POJGIe1eGw"

data = scattermapbox(
    lat=[45.5017], lon=[-73.5673], text=["Montreal"],
    mode="markers", marker_size=14, 
)

layout = Layout(
    autosize=true, hovermode=closest,
    mapbox=attr(
        accesstoken=mapbox_access_token,
        bearing=0,
        center_lat=45, center_lon=-73,
        pitch=0,
        zoom=5
    )
)

plot([data], layout)
