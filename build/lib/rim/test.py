import RIM

# import other packages
import json
from shapely.geometry import shape
import csv

# define the path to the geojson file holding the data about unimelb buildings
data_path = 'rim/testdata/29_unimelb_buildings.geojson'

# read the contents of the file into a string
with open(data_path, 'r') as f:
    data_string = f.read()

# read the contents of the string into a json object
data = json.loads(data_string)
# get the features from the json object -- each feature is a building
features = data['features']

def json2wkt(feature):
    """ A helper function that takes in a feature from the features 
    object and returns its geometry as a WKT string """
    shapelyshape = shape(feature['geometry'])
    return shapelyshape.wkt

# create a dictionary from which buildings geometries (in WKT string format) can be accessed by building ids
# the keys are building ids and values are WKT geometries {bid : WKTgeom}
buildings_data = {
    feature['properties']['bid'] : json2wkt(feature) for feature in features
}

buildings_data

# read in the list of triplets of buildings for RIM calculation
filepath = 'rim/testdata/unimelb_buildings_triplets.csv'
with open(filepath, 'r') as f:
    reader = csv.DictReader(f, delimiter=';')
    triplets = [row for row in reader]
    
rimobjects = []
for row in triplets:
    a_geom = buildings_data.get(int(row['A']))
    b_geom = buildings_data.get(int(row['B']))
    o_geom = buildings_data.get(int(row['O']))
    try:
        rimobject = RIM.RIM(a_geom, b_geom, o_geom)
        rimobjects.append(rimobject.rim)
        if '[' in rimobject.rim:
            print(rimobject.rim)
    except Exception as e:
        rimobjects.append(str(e))
        print(e)

import collections

ctr = collections.Counter(rimobjects)
print(ctr)

print(rimobjects)