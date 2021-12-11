import requests
import sys
import os

if __name__ == '__main__':
    country = 'netherlands'
    if len(sys.argv) == 2:
        country = sys.argv[1]

    url = f'https://download.geofabrik.de/europe/{country}-latest.osm.pbf'
    # url = 'https://download.geofabrik.de/europe-latest.osm.pbf'

    dir = f'data/{country}'
    file_name = f'{dir}/{country}.osm.pbf'

    if not os.path.exists(dir):
        os.makedirs(dir)

    with open(file_name, "wb") as f:
        print("Downloading %s" % file_name)
        response = requests.get(url, stream=True)
        total_length = response.headers.get('content-length')

        if total_length is None:  # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50-done)))
                sys.stdout.flush()
