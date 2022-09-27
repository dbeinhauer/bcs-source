#!/usr/bin/env python3

import numpy as np

# File informations:
CITY_COLUMNS = 4
NODE_COLUMNS = 4

CITY_ID_COLUMN = 0
CITY_LAT_COLUMN = 1
CITY_LON_COLUMN = 2
CITY_POPULATION_COLUMN = 3

NEW_NODE_ID_COLUMN = 0
OLD_NODE_ID_COLUMN = 1
NODE_LAT_COLUMN = 2
NODE_LON_COLUMN = 3


NODE_CITY_FILE_HEADER = "node_id latitude longitude city_id distance_km\n"


class NodePreprocessor:
    """
    Class for combining node a city files and computing the closest city of each node.
    """
    def __init__(self):
        # Array of the cities coordinates (each row is `(lat, lon)`).
        self.cities = []

    def load_cities(self, filename):
        """
        Loads coordinates of all cities and stores them into numpy array
        with shape (num_cities, 2).
        Args:   
            filename: name of the file where the city \
                                representation is stored
        """
        with open(filename, "r") as cities_file:
            # skip header
            line = cities_file.readline()
            
            while True:
                # Read until end of the file.
                line = cities_file.readline()

                # End of the file
                if not line:
                    break

                splited_line = line.split()
                if len(splited_line) == CITY_COLUMNS:
                    # City is represented by correct number of column.
                    self.cities.append((float(splited_line[CITY_LAT_COLUMN]), float(splited_line[CITY_LON_COLUMN])))

        self.cities = np.array(self.cities)

    
    def spherical_dist(self, pos1, pos2, r=6371):
        """
        Computes spherical distance of the two numpy matrices of the coordinates 
        (2 columns (lat, lon)).
        Args:  
            pos1: numpy matrices of coordinates \
                (computes distance between each 2 points from the matrices)
            pos2: numpy matrices of coordinates \
                (computes distance between each 2 points from the matrices)
            r: radius of the sphere in kilometers (by default earth radius)

        Returns:
        Returs matrix of the distances of between each 2 points 
        from `pos1` and `pos2` in kilometers. 

        Inspired by:
        https://stackoverflow.com/questions/19413259/efficient-way-to-calculate-distance-matrix-given-latitude-and-longitude-data-in
        """

        pos1 = pos1 * np.pi / 180
        pos2 = pos2 * np.pi / 180
        cos_lat1 = np.cos(pos1[..., 0])
        cos_lat2 = np.cos(pos2[..., 0])
        cos_lat_d = np.cos(pos1[..., 0] - pos2[..., 0])
        cos_lon_d = np.cos(pos1[..., 1] - pos2[..., 1])
        return r * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))

    
    def find_closest_city(self, coordinates):
        """
        Returns:
        Returns id and distance of the closest city from the given coordinates. 
        """
        distance_matrix = self.spherical_dist(self.cities, coordinates)

        # Find the closest city and its distance from the point.
        best_id = np.argmin(distance_matrix)
        min_distance = np.min(distance_matrix)

        return best_id, min_distance


    
    def combine_nodes_cities(self, input_filename_nodes, output_filename):
        """
        Finds for each node its closest city and saves the info into the 
        corresponding file, where each row is in format: 
        `node_id latitude longitude city_id city_distance_km`.
        Args: 
            input_filename_nodes: name of the file with \ 
                nodes representation
            output_filename: name of the output file with node + city representation
        """
        with open(output_filename, 'w') as output_file:
            # Write header.
            output_file.write(NODE_CITY_FILE_HEADER)
            
            # Read node file line by line and find corresponding city.
            with open(input_filename_nodes, 'r') as input_file:
                # Skip header.
                line = input_file.readline()
                counter = 0

                while True:
                    # Read until end of the file.

                    # Just to check progres of the program.
                    if counter % 1000 == 0:
                        print(f"Current number: {counter}")


                    line = input_file.readline()

                    if not line:
                        # End of the file -> finish reading
                        break


                    splited = line.split()
                    if len(splited) == NODE_COLUMNS:
                        # Correct format of the node
                        node_id, lat, lon = splited[NEW_NODE_ID_COLUMN], \
                                            float(splited[NODE_LAT_COLUMN]), \
                                            float(splited[NODE_LON_COLUMN])

                        best_city_id, distance = self.find_closest_city(np.array([(lat, lon)]))

                        output_text = " ".join([str(node_id),
                                                str(lat),
                                                str(lon),
                                                str(best_city_id),
                                                str(distance),
                                                "\n"])
                                                
                        output_file.write(output_text)
                        counter += 1