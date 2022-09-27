#!/usr/bin/env python3


import argparse

import map_parser
import node_preprocessor


parser = argparse.ArgumentParser()
parser.add_argument("--input_map", default=None, type=str, help="Input file with map representation in `osm.pbf` file format.")
parser.add_argument("--input_cities", default=None, type=str, help="Input file with cities representation in `osm.pbf` file format.")
parser.add_argument("--output_nodes", default=None, type=str, help="Output file of the nodes representation.")
parser.add_argument("--output_edges", default=None, type=str, help="Output file of the edges representation.")
parser.add_argument("--output_cities", default=None, type=str, help="Output file of the cities representation.")
parser.add_argument("--output_combined", default=None, type=str, help="Output file of the nodes-cities representation.")


# Default input filenames:
MAP_CZECH_ROADS = "maps/top-roads-czech-republic.osm.pbf"
CITY_FILE = "maps/czech_cities.osm.pbf"


# Output filenames:
NODES_OUTPUT_FILE = "prepared_graph/nodes.txt"
EDGES_OUTPUT_FILE = "prepared_graph/edges.txt"
CITIES_OUTPUT_FILE = "prepared_graph/cities.txt"
COMBINED_NODES_OUTPUT_FILE = "prepared_graph/combined_nodes.txt"

 
def main(args):
    # Initialize all filename variables.
    if args.input_map == None:
        args.input_map = MAP_CZECH_ROADS

    if args.input_cities == None:
        args.input_cities = CITY_FILE
    
    if args.output_nodes == None:
        args.output_nodes = NODES_OUTPUT_FILE
    
    if args.output_edges == None:
        args.output_edges = EDGES_OUTPUT_FILE

    if args.output_cities == None:
        args.output_cities = CITIES_OUTPUT_FILE

    if args.output_combined == None:
        args.output_combined = COMBINED_NODES_OUTPUT_FILE

    # Parse only necessary info about the network and cities 
    # -> store it into the corresponding file.
    map_parser.parse_driving_network(args.input_map, args.output_nodes, args.output_edges)
    map_parser.parse_city_file(args.input_cities, args.output_cities)

    print("Nodes and edges loaded!s")

    # Find the nearest city for each node and store the combination into corresponding file.
    node_proc = node_preprocessor.NodePreprocessor()
    node_proc.load_cities(args.output_cities)
    node_proc.combine_nodes_cities(args.output_nodes, args.output_combined)



if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
