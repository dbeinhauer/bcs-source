# About
This directory contains program `map_reader.py` which converts the
maps from format `.osm.pbf` (which can be obtained from the server
[Openstreet Maps]()) to input format of the traffic simulator.

# Recommended Preprocessing Tool Osmium
The traffic simulator is mainly designed to simulate large scale 
traffic networks of the size comparable with the size of the
traffic network of the Czech Republic. Unfortunately the size of the 
map files from OpenStreet Maps which represents such networks are 
often very large and it is recommended to reduce them before 
converting them to the requiered input format.

In order to reduce these files it is recomended to use the tool
[Osmium]() to filter only traffic network and information about
cities from the original files.

## Example
In the example below we show the recommended reduction of
the original map file into two new reduced files with the 
representation of the road network with the roads of the choosen
types respectively with the representation of the city objects
of the choosen type.

Map reduction to the filtered road network file:

```bash
osmium tags-filter czech-republic-latest.osm.pbf \
w/highway=motorway,trunk,primary,motorway_link,trunk_link, primary_link \
-o top-roads-czech-republic.osm.pbf
```

Map reduction to the filtered cities file:

```bash
osmium tags-filter maps/czech-republic-latest.osm.pbf \
place=city,town \
-o maps/czech_cities.osm.pbf
```

# Map Convertion
In order to convert the map representation files from `.osm.pbf` 
format to the input format of the traffic simulator use the program
`map_reader.py`. It reads the files of the road network and cities 
in `osm.pbf` format and converts them to 4 files in the input format 
of the traffic simulator which represent nodes, edges, cities and 
node-city pairs for specifying the information about the context
of the individual junctions.

## Example
The following command shows the example of the program usage.

```bash
./map_reader.py \
--input_map=$(road_network_input_file_OSM_PBF_format) \
--input_cities=$(cities_input_file_OSM_PBF_format) \
--output_nodes=$(nodes_output_in_simulator_format) \
--output_edges=$(edges_output_in_simulator_format) \
--output_cities=$(cities_output_in_simulator_format) \
--output_combined=$(nodes_cities_output_in_simulator_format)
```

## Input And Output Data
Note that there are two empty subdirectories `maps/` (input) and
`prepared_graph/` (output) which should serve as the storage of the
corresponding files. These directories are empty because of the size
of the files which represent the maps.