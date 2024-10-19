import matplotlib.pyplot as plt

def read_g2o_file(filename):
    vertices = {}
    edges = []

    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if parts[0] == 'VERTEX_SE2':
                vertex_id = int(parts[1])
                x = float(parts[2])
                y = float(parts[3])
                theta = float(parts[4])
                vertices[vertex_id] = (x, y, theta)
            elif parts[0] == 'EDGE_SE2':
                from_id = int(parts[1])
                to_id = int(parts[2])
                dx = float(parts[3])
                dy = float(parts[4])
                dtheta = float(parts[5])
                edges.append((from_id, to_id, dx, dy, dtheta))

    return vertices, edges

def plot_graph(vertices, edges):
    fig, ax = plt.subplots()

    # Plot vertices
    for vertex_id, (x, y, theta) in vertices.items():
        ax.plot(x, y, 'bo')  # Blue dot for vertices
        ax.text(x, y, str(vertex_id), fontsize=12)

    # Plot edges
    for from_id, to_id, dx, dy, dtheta in edges:
        x_from, y_from, _ = vertices[from_id]
        x_to, y_to, _ = vertices[to_id]
        ax.plot([x_from, x_to], [y_from, y_to], 'r-')  # Red line for edges

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Graph Plot from .g2o File')
    plt.show()

# Read the .g2o file
vertices, edges = read_g2o_file('build/before_optimization.g2o')

# Plot the graph
plot_graph(vertices, edges)


# Read the .g2o file
vertices, edges = read_g2o_file('build/after_optimization.g2o')

# Plot the graph
plot_graph(vertices, edges)