import networkx as nx
import numpy as np
import json
import sys
import subprocess
import os
import yaml


def print_if_verbose(verbose, *args):
    if verbose:
        print(*args)


def ANF(country, distance=5, r=7, k=128, verbose=False):
    dir = f'./data/{country}'
    if not os.path.exists(dir):
        os.makedirs(dir)

    if not os.path.exists(f'{dir}/{country}.osm.pbf'):
        print_if_verbose(verbose, f'Downloading {country}...')
        subprocess.run(['python3', './download.py', country], capture_output=False)
    else:
        print_if_verbose(verbose, f'{country} already downloaded.')

    print_if_verbose(verbose, f'Converting {country}...')

    if not os.path.exists(f'{dir}/{country}-highways.xml'):
        filter_tags_command = f"osmium tags-filter {dir}/{country}.osm.pbf w/highway=motorway,primary,trunk,motorway_link,trunk_link,primary_link -o {dir}/{country}-highways.xml --overwrite"
        subprocess.run(filter_tags_command, shell=True, capture_output=False)
    else:
        print_if_verbose(verbose, f'{country}-highways.xml already exists.')

    if not os.path.exists(f"{dir}/{country}-highways.json"):
        create_networkx_json = f"python ./libs/OsmToRoadGraph/run.py -f {dir}/{country}-highways.xml --networkx"
        subprocess.run(create_networkx_json, shell=True, capture_output=False)
    else:
        print_if_verbose(verbose, f'{country}-highways.json already exists.')

    file_name = f"{dir}/{country}-highways.json"

    if not os.path.exists(f'{dir}/{country}-highways.csv'):
        with open(file_name) as f:
            json_data = json.load(f)

        G = nx.readwrite.json_graph.adjacency_graph(json_data)
        nx.write_edgelist(G, f'{dir}/{country}-highways.csv', data=False)
    else:
        print_if_verbose(verbose, f'{country}-highways.csv already exists.')

    if not os.path.exists(f'{dir}/{country}.png'):
        create_image = f"python libs/OsmToRoadGraph/examples/pycgr-to-png/pycgr-to-png.py -f {dir}/{country}-highways.pypgr -o {dir}/{country}.png"
        subprocess.run(create_image, shell=True, capture_output=False)
    else:
        print_if_verbose(verbose, f'{country}.png already exists.')

    # run ANF
    anf = f"./ANF {dir}/{country}-highways.csv {distance} {r} {k}"
    print_if_verbose(verbose, f'Running ANF on {country}...')
    result = subprocess.run(anf, shell=True, stdout=subprocess.PIPE)
    x = json.loads(result.stdout.decode('utf-8'))
    print_if_verbose(verbose, f'\n\n\n---------------------{country}---------------------\n')
    print_if_verbose(verbose, yaml.dump(x))
    return x


if __name__ == '__main__':
    country = 'netherlands'
    distance = 5
    r = 7
    k = 128

    if len(sys.argv) >= 2:
        country = sys.argv[1]

    if len(sys.argv) >= 3:
        distance = sys.argv[2]

    if len(sys.argv) >= 4:
        r = sys.argv[3]

    if len(sys.argv) >= 5:
        k = sys.argv[4]

    ANF(country, distance, r, k, verbose=True)
