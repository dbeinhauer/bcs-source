#!/usr/bin/env python3

import json
import pyrosm


# File headers:
CITY_FILE_HEADER = "city_id lat lon population\n"
NODE_FILE_HEADER = "new_node_id node_id lat lon\n"
EDGE_FILE_HEADER = "new_node_id_1 new_node_id_2 length_km road_type\n"


# Filters:
CITY_FILTER = {"place":True}

def parse_city_file(input_file, output_file):
    """
    Reads osm.pbf file (prefiltered recomended), filters all inhabited areas
    (specified in the `CITY_FILTER`). Then it takes only the areas with 
    specified `population` tag and stores parameters in the following format
    (each city on the separate line):
            `city_id lat lon population`

    Args:
        input_file: name of the input file
        output_file: name of the output file
    """
    
    osm = pyrosm.OSM(input_file)

    # Take only data representing cities.
    cities = osm.get_pois(custom_filter=CITY_FILTER)
    print("Data loaded!")    

    city_iterator = cities.iterfeatures()

    # Write preprocessed data into the given output file.
    with open(output_file, "w") as f:
        f.write(CITY_FILE_HEADER)

        city_id = 0
        for city in city_iterator:
            output = ""
            tags_dictionary = json.loads(city['properties']['tags'])

            # Take only areas with specified latitude, longitude and population.
            if 'population' in tags_dictionary.keys():
                if city['properties']['lat'] == None or city['properties']['lon'] == None:
                    continue

                output = " ".join([str(city_id),
                                    str(city['properties']['lat']),
                                    str(city['properties']['lon']),
                                    str(tags_dictionary['population']),
                                    "\n"])
                city_id += 1
                
            f.write(output)



def convertRoadType(road):
    """
    Checks road type and converts it into `char` representation.
    Returns:
                `m` - motorway
                `t` - trunk
                `p` - primary
                `s` - secondary
                ' ' - other 
    """
    if road == "motorway" or "motorway_link":
        return 'm'
    elif road == "trunk" or "trunk_link":
        return 't'
    elif road == "primary" or "primary_link":
        return 'p'
    elif road == "secondary":
        return 's'
    else:
        return ' '


def parse_driving_network(input_file, output_file_nodes, output_file_edges):
    """
    Reads osm.pbf file (prefiltered recomended), extracts driving network from
    the map (only the highest road types) and converts it into two files 
    (nodes and edges representation).
    Nodes file format (each node separate line):
                `new_node_id original_node_id lat lon`
                
    Edges file format (each edge separate line):
                `new_node_id_1 new_node_id_2 length road_type`
            Note that `new_node_id_1` is always smaller 
            than `new_node_id_2` (for future preprocessing).

    Args: 
        input_file: name of the input file
        output_file_nodes: name of the output file with \
            node representation
        output_file_edges: name of the output file with \
            edge representation
    """

    osm = pyrosm.OSM(input_file)

    # Take only drivable roads.
    nodes, edges = osm.get_network(network_type="driving", nodes=True)

    # Dictionary for storing key-value pair (old_id, new_id)
    nodes_new_ids = {}

    nodes_iterator = nodes.iterfeatures()
    # Write preprocessed node data into the specified file.
    with open(output_file_nodes, "w") as f:
        f.write(NODE_FILE_HEADER)

        # New id of the following node from the file.
        next_node_id = 0
        for node in nodes_iterator:
            # Take only nodes with specified latitude and longitude.
            if node['properties']['lat'] == None or node['properties']['lon'] == None:
                continue

            # Prepare each node representation.
            output = " ".join([str(next_node_id),
                                str(node['properties']['id']),
                                str(node['properties']['lat']),
                                str(node['properties']['lon']),
                                '\n'])
            f.write(output)
            
            # Add (old_id, new_id) pair into the dictionary and increase node counter.
            nodes_new_ids[node['properties']['id']] = next_node_id
            next_node_id += 1



    edges_iterator = edges.iterfeatures()
    # Write preprocessed edge data into the specified file.
    with open(output_file_edges, "w") as f:
        f.write(EDGE_FILE_HEADER)
        for edge in edges_iterator:
            # Prepare each edge

            # Find corresponding new node ids for the old ids. 
            node_u = nodes_new_ids[edge['properties']['u']]
            node_v = nodes_new_ids[edge['properties']['v']]


            # Swap nodes to have the smaller node id first in the edges file (for future preprocessing).
            if node_u > node_v:
                node_u, node_v = node_v, node_u
            
            output = " ".join([str(node_u),
                                str(node_v),
                                str(float(edge['properties']['length'] / 1000)),
                                str(convertRoadType(edge['properties']['highway'])),
                                "\n"])
            f.write(output)
