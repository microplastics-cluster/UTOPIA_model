import matplotlib.pyplot as plt
import networkx as nx

# Define the nodes (boxes) and their mass concentration values
nodes = {'Ocean_Surface_Water': 5.7984E-12, 'Ocean_Mixed_Water': 5.31E-15}#, 'Ocean_Column_Water': 6.01E-15,'Coast_Surface_Water': 2.32E-13, "Coast_Column_Water":7.49E-16,"Surface_Freshwater":1.46E-07,"Bulk_Freshwater":5.99E-10,"Sediment_Freshwater":0.009963974,"Sediment_Ocean":1.58E-06,"Sediment_Coast":7.72E-08,"Urban_Soil_Surface":3.86E-26,"Urban_Soil":0,"Background_Soil_Surface":3.86E-26,"Background_Soil":0,"Agricultural_Soil_Surface":9.64E-27,"Agricultural_Soil":0,"Air":1.17E-33}

# Define the edges (connections between nodes) and their mass flux values
edges = [('Ocean_Surface_Water', 'Ocean_Mixed_Water', 0.001454389),('Ocean_Mixed_Water', 'Ocean_Surface_Water', 4.38198E-21)]#('Ocean_Surface_Water', 'Coast_Surface_Water', 0.000128146),('Ocean_Surface_Water', 'Air', 4.70741E-21)] #('Ocean_Mixed_Water', 'Ocean_Surface_Water', rising), ('Ocean_Mixed_Water', 'Ocean_Column_Water', settling)

# Create a directed graph
G = nx.DiGraph()

# Add nodes with mass concentration as node attribute
for node, concentration in nodes.items():
    G.add_node(node, concentration=concentration)

# Add edges with mass flux as edge attribute
for edge in edges:
    G.add_edge(edge[0], edge[1], flux=edge[2])

# Manually specify positions for the nodes in a grid-like layout
pos = {'Ocean_Surface_Water': (1, 2), 'Ocean_Mixed_Water': (1, 1)}

# Draw nodes with their concentration as labels
nx.draw_networkx_nodes(G, pos, node_size=1000, node_color='lightblue')
nx.draw_networkx_labels(G, pos, font_size=12, font_color='black', font_weight='bold', labels={node: f'{node}\n{concentration}' for node, concentration in nodes.items()})

# Draw edges with their flux as labels
nx.draw_networkx_edges(G, pos, arrows=True)
nx.draw_networkx_edge_labels(G, pos, edge_labels={(edge[0], edge[1]): str(edge[2]) for edge in edges}, font_color='red')

# Set the figure size to fit the diagram nicely
plt.figure(figsize=(8, 6))

# Display the plot
plt.axis('off')
plt.show()
