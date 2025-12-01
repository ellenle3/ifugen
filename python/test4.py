import plotly.graph_objects as go

# Your coordinates
x = [0, 1, 0]
y = [0.0, 0.0, 0.0]   # all y = 0
z = [0, 0, 1]

fig = go.Figure(data=[
    go.Mesh3d(
        x=x,
        y=y,
        z=z,
        i=[0],
        j=[1],
        k=[2],
        color='lightblue',
        opacity=0.7,
        flatshading=True
    )
])

fig.update_layout(
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        aspectmode='data'
    ),
    title='Single Triangle with Given Coordinates'
)

fig.show()
