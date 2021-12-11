# %%
import networkx as nx
import json

# osmium tags-filter malta-latest.osm.pbf w/highway=motorway,primary,trunk -o malta-highways.xml --overwrite

# not working:
# osmosis --read-xml file="data/croatia/croatia-highways.xml" --write-xml file="croatia-highways.xml" --networkx

# working:
# python ../libs/OsmToRoadGraph/run.py -f croatia-highways.xml --networkx

country = "croatia"

file_name = f"data/{country}/{country}-highways.json"

with open(file_name) as f:
    json_data = json.load(f)

G = nx.readwrite.json_graph.adjacency_graph(json_data)
nx.write_edgelist(G, f'data/{country}/{country}-highways.csv', data=False)
nx.info(G)
# %%
