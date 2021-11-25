# %%
import networkx as nx
import json

file_name = "data/highways.json"

with open(file_name) as f:
    json_data = json.load(f)

G = nx.readwrite.json_graph.adjacency_graph(json_data)
nx.write_edgelist(G, 'data/highways.csv', data=False)
