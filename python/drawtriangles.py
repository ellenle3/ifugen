import numpy as np
import plotly.graph_objects as go

# ----------------------------------------------------------------------
# Load triangles from file
# Each line: x1,y1,z1,x2,y2,z2,x3,y3,z3,code
# ----------------------------------------------------------------------
def load_triangles(filename="../src/triangles.txt"):
    triangles = []
    with open(filename, "r") as f:
        for line in f:
            vals = line.strip().split(",")
            if len(vals) != 10:
                continue  # skip malformed lines

            nums = list(map(float, vals[:-1]))
            code = int(vals[-1])

            p1 = nums[0:3]
            p2 = nums[3:6]
            p3 = nums[6:9]

            triangles.append((p1, p2, p3, code))

    return triangles


# ----------------------------------------------------------------------
# Decode visibility flags
# Invisible Edges:
#   P1->P2 = 1
#   P2->P3 = 2
#   P3->P1 = 4
# ----------------------------------------------------------------------
def edge_visibility(code):
    invisible = {
        (0, 1): bool(code & 1),  # P1->P2
        (1, 2): bool(code & 2),  # P2->P3
        (2, 0): bool(code & 4),  # P3->P1
    }
    return invisible


# ----------------------------------------------------------------------
# Generate Plotly figure for triangles
# ----------------------------------------------------------------------
def plot_triangles(triangles):
    fig = go.Figure()

    for (p1, p2, p3, code) in triangles:
        pts = np.array([p1, p2, p3, p1])  # wrap back to first point for edges

        invis = edge_visibility(code)
        edges = [(0, 1), (1, 2), (2, 0)]

        # Draw edges with visibility
        for (i, j) in edges:
            is_invisible = invis[(i, j)]
            style = {
                "width": 2,
                "color": "black"
            } if not is_invisible else {
                "width": 1,
                "color": "lightgray"
            }

            fig.add_trace(go.Scatter3d(
                x=[pts[i][0], pts[j][0]],
                y=[pts[i][1], pts[j][1]],
                z=[pts[i][2], pts[j][2]],
                mode="lines",
                line=style,
                showlegend=False
            ))

        # Draw triangle face
        fig.add_trace(go.Mesh3d(
            x=pts[:3, 0],
            y=pts[:3, 1],
            z=pts[:3, 2],
            i=[0],
            j=[1],
            k=[2],
            color="lightblue",
            opacity=0.35,
            hoverinfo="skip"
        ))
    

    fig.update_layout(
        scene=dict(aspectmode='data'),
        title="Triangle Mesh Visualization",
        width=900,
        height=700
    )

    fig.show()


# ----------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("Loading triangles from triangles.txt...")
    triangles = load_triangles()

    print(f"Loaded {len(triangles)} triangles.")
    print("Generating interactive 3D plot...")
    plot_triangles(triangles)

