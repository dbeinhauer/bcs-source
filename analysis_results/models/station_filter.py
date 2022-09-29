#!/usr/bin/env python3
"""
Filters all used charging stations and to write them to the stdin.
"""

import argparse
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--inputStations",
    default="genetic.txt", 
    type=str,
    help="Filename of all charging stations.")
parser.add_argument("--inputUsed", 
    default="genetic_usage.txt", 
    type=str, 
    help="Filename of stations usage.")


def main(args):
    stations = []
    usage = []
    header = ""

    with open(args.inputStations, 'r') as f:
        header = f.readline()
        for line in f:
            stations.append(line)

    with open(args.inputUsed, 'r') as f:
        for line in f:
            if line.startswith("Number of waitings:"):
                usage.append(int(line.split(":")[1].split(" ")[1]))

    stations = np.array(stations)
    usage = np.array(usage)
    indices = np.where(usage != 0)

    # All at least once used charging stations.
    stations = stations[indices]

    # Print all used charging stations.
    print(header, end="")
    for line in stations:
        print(line, end="")

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)