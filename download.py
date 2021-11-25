import requests
import sys

if __name__ == '__main__':
  url = 'https://download.geofabrik.de/europe/malta-latest.osm.pbf'
  # url = 'https://download.geofabrik.de/europe-latest.osm.pbf'

  file_name = 'data/malta.osm.pbf'

  with open(file_name, "wb") as f:
    print("Downloading %s" % file_name)
    response = requests.get(url, stream=True)
    total_length = response.headers.get('content-length')

    if total_length is None: # no content length header
        f.write(response.content)
    else:
        dl = 0
        total_length = int(total_length)
        for data in response.iter_content(chunk_size=4096):
            dl += len(data)
            f.write(data)
            done = int(50 * dl / total_length)
            sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50-done)) )
            sys.stdout.flush()
