#!/usr/bin/env python3

from inspect import GEN_CLOSED
from turtle import color
from unittest import result
import matplotlib
import matplotlib.pyplot as plt

import argparse
import numpy as np

import copy
# from sklearn.cluster import k_means




parser = argparse.ArgumentParser()

parser.add_argument("--input", default="result_overview/parameters/run2.txt", 
    type=str, help="Input file with the results.")
parser.add_argument("--first", default="numClosestStations", 
    type=str, help="Input file with the results.")
parser.add_argument("--second", default="runDownBat", 
    type=str, help="Second parameter to plot.")
parser.add_argument("--plot_type", default="runDown", 
    type=str, help="")
parser.add_argument("--output", default=None, 
    type=str, help="Input file with the results.")

# red, blue, orange, light blue, yellow
# COLORS = ['#d7191c','#2c7bb6','#fdae61','#abd9e9', '#ffffbf']
COLORS = ['#d7191c', '#2c7bb6', '#fdae61', '#abd9e9', '#ffffbf', '#abd9e9']

RANDOM_FILE = "results/random.txt"
GREEDY_FILE = "results/greedy.txt"
GENETIC_FILE = "results/genetic.txt"
KMEANS_FILE = "results/kmeans.txt"

# RANDOM_FILE = "../analysis_results/random_results/results.txt"
# GREEDY_FILE = "../analysis_results/greedy_results/results.txt"
# GENETIC_FILE = "../analysis_results/genetic_results/results.txt"
# KMEANS_FILE = "../analysis_results/kmeans_results/results.txt"


def ReadData(filename, first, second, specify_par=None):

    results = {}

    first_col, second_col = 0, 0
    specify_col = 0

    with open(filename, 'r') as fp:
        header = fp.readline().split()
        for i, key in enumerate(header):
            if key == first:
                first_col = i           
            elif key == second:
                second_col = i
            elif key == specify_par:
                specify_col = i
        
        fp.readline()

        cur_data = ([], [])
        for line in fp:
            splited = line.split()

            if specify_par == None:
                first_data = cur_data[0]
                second_data = cur_data[1]
                first_data.append(splited[first_col])
                second_data.append(splited[second_col])
                cur_data = (first_data, second_data)

            else:
                cur_data = ([], [])
                if splited[specify_col] in results:
                    cur_data = results[splited[specify_col]]

                first_data = cur_data[0]
                second_data = cur_data[1]
                first_data.append(splited[first_col])
                second_data.append(splited[second_col])
                results[splited[specify_col]] = (first_data, second_data)

    if specify_par == None:
        results['all'] = cur_data

    return results



def PlotRunDownStations(args, num_stations=None):

    # random_file = "results/random.txt"
    # greedy_file = "results/greedy.txt"
    # genetic_file = "results/genetic.txt"
    # kmeans_file = "results/kmeans.txt"

    first_column = "numStations"
    second_column = "runDownBat"

    first_column = "numUsedStat"
    # second_column = "waitTime"

    results_greedy = ReadData(GREEDY_FILE, first_column, second_column)
    results_genetic = ReadData(GENETIC_FILE, first_column, second_column)
    results_kmeans = ReadData(KMEANS_FILE, first_column, second_column)
    results_random = ReadData(RANDOM_FILE, first_column, second_column)

    all_results = [results_greedy, results_genetic, results_kmeans, results_random]
    all_results = [results_greedy, results_genetic, results_kmeans]

    labels = ["Hladový přístup - průměr", "Genetický algoritmus - průměr", "K-Means - průměr", "Náhodný přístup - průměr"]
    scatter_labels = ["Hladový přístup - experimenty", "Genetický algoritmus - experimenty", "K-Means - experimenty", "Náhodný přístup - - experimenty"]

    font = {'size'   : 13}

    matplotlib.rc('font', **font)

    fig, ax = plt.subplots()
    fig.set_figheight(7.1)
    fig.set_figwidth(9.15)
    # ax.set_xlabel('Počet umístěných nabíjecích stanic', size=22)
    ax.set_xlabel('Počet použitých nabíjecích stanic', size=22)
    ax.set_ylabel('Počet vybitých vozidel', size=22)
    

    # ax.set_title('Porovnání optimalizačních metod s ohledem na\npočet rozmístěných stanic (původní experimenty)', fontweight="bold", size=23)
    # ax.set_title('Porovnání optimalizačních metod s ohledem na\npočet rozmístěných stanic', fontweight="bold", size=23)
    ax.set_title('Porovnání optimalizačních metod s ohledem\nna počet použitých stanic', fontweight="bold", size=23)
    

    # Show optimization results.
    for i, results in enumerate(all_results):
        for _, item in enumerate(results.items()):
            key = item[0]
            value = item[1]
            x = np.array([float(m) for m in value[0]])
            y = np.array([float(m) for m in value[1]])

            indices = np.argsort(x)
            x = x[indices]
            y = y[indices]

            unique_x = np.unique(x)
            y_means = copy.deepcopy(y)
            for un in unique_x:
                indices = np.where(x == un)
                mean_y = np.mean(y[indices])  
                y_means[indices] = mean_y

            ax.plot(x, y_means, color=COLORS[i], label=labels[i])
            # ax.scatter(x, y_means, color=COLORS[i])#, label=scatter_labels[i])
            
            # ax.plot(x, y, color=COLORS[i], label=labels[i])
            # ax.scatter(x, y, color=COLORS[i], label=scatter_labels[i])
            # curr_plot = ax.plot(x, y, color=COLORS[i], label=key)


    # Show random results
    for _, item in enumerate(results_random.items()):
            key = item[0]
            value = item[1]
            x = np.array([float(m) for m in value[0]])
            y = np.array([float(m) for m in value[1]])
            curr_plot = ax.plot(x, y, color=COLORS[3], linewidth=2, linestyle='dashdot', label=labels[3])
            curr_plot = ax.scatter(x, y, color=COLORS[3])#, label=scatter_labels[3])

    ax.legend(fontsize=16.2)

    if args.output != None:
        plt.savefig(args.output)
    else:
        plt.show()



def main(args):
    # PlotRunDown(args)
    # PlotRunDownStations(args, num_stations=500)
    PlotRunDownStations(args)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)