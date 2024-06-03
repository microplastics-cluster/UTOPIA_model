""" Sample script to visualize the massbalance diagram of Utopia. """

# Under developement

import matplotlib.pyplot as plt
import networkx as nx


def generate_mas_balance_diagram(model_lists, Results):
    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes with concentration values

    compartments = model_lists["compartmentNames_list"]
    comp_mass_concentrations = [
        sum(Results[Results["Compartment"] == c]["concentration_g_m3"])
        for c in compartments
    ]

    nodes = compartments
    concentration_values = dict(zip(compartments, comp_mass_concentrations))

    # Define the nodes (boxes) and their mass concentration values
    for node in nodes:
        G.add_node(node, concentration=concentration_values[node])

    # Add connections with flow values
    connections = [
        ("Node A", "Node B", {"flow": 10}),
        ("Node A", "Node C", {"flow": 5}),
        ("Node B", "Node C", {"flow": 7}),
        ("Node C", "Node D", {"flow": 3}),
        ("Node B", "Node D", {"flow": 8}),
    ]
    G.add_edges_from(connections)

    # Set node shapes to 'box'
    node_shapes = {node: "box" for node in nodes}

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes with mass concentration as node attribute
    for node, concentration in concentration_values.items():
        G.add_node(node, concentration=concentration)

    # Add edges with mass flux as edge attribute
    for edge in connections:
        G.add_edge(edge[0], edge[1], flux=edge[2])

    # Manually specify positions for the nodes in a grid-like layout
    pos = {"Ocean_Surface_Water": (1, 2), "Ocean_Mixed_Water": (1, 1)}

    # Draw nodes with their concentration as labels
    nx.draw_networkx_nodes(G, pos, node_size=1000, node_color="lightblue")
    nx.draw_networkx_labels(
        G,
        pos,
        font_size=12,
        font_color="black",
        font_weight="bold",
        labels={
            node: f"{node}\n{concentration}" for node, concentration in nodes.items()
        },
    )

    # Draw edges with their flux as labels
    nx.draw_networkx_edges(G, pos, arrows=True)
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels={(edge[0], edge[1]): str(edge[2]) for edge in connections},
        font_color="red",
    )

    # Set the figure size to fit the diagram nicely
    plt.figure(figsize=(8, 6))

    # Display the plot
    plt.axis("off")
    plt.show()
