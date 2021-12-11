import networkx as nx
import json
import sys
import subprocess
import os

if __name__ == '__main__':
    country = 'netherlands'
    if len(sys.argv) == 2:
        country = sys.argv[1]

    dir = f'./data/{country}'
    if not os.path.exists(dir):
        os.makedirs(dir)

    if not os.path.exists(f'{dir}/{country}.osm.pbf'):
        print(f'Downloading {country}...')
        subprocess.run(['python3', './download.py', country])

    filter_tags_command = f"osmium tags-filter {dir}/{country}.osm.pbf w/highway=motorway,primary,trunk -o {dir}/{country}-highways.xml --overwrite"
    subprocess.run(filter_tags_command, shell=True)

    create_networkx_json = f"python ./libs/OsmToRoadGraph/run.py -f {dir}/{country}-highways.xml --networkx"
    subprocess.run(create_networkx_json, shell=True)

    file_name = f"{dir}/{country}-highways.json"

    with open(file_name) as f:
        json_data = json.load(f)

    G = nx.readwrite.json_graph.adjacency_graph(json_data)
    nx.write_edgelist(G, f'{dir}/{country}-highways.csv', data=False)

    create_image = f"python libs/OsmToRoadGraph/examples/pycgr-to-png/pycgr-to-png.py -f {dir}/{country}-highways.pypgr -o {dir}/{country}.png"
    subprocess.run(create_image, shell=True)
