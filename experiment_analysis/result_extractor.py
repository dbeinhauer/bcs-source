#!/usr/bin/env python3
"""
Prepares the plots of the experiment results.
"""

import matplotlib
import matplotlib.pyplot as plt

import argparse
import numpy as np

import copy


parser = argparse.ArgumentParser()
parser.add_argument("--input", 
    default="result_overview/parameters/run2.txt", 
    type=str,
    help="Input file with the results.")
parser.add_argument("--first", 
    default="numClosestStations", 
    type=str,
    help="Input file with the results.")
parser.add_argument("--second", 
    default="runDownBat", 
    type=str, 
    help="Second parameter to plot.")
parser.add_argument("--plot_type", 
    default="runDown", 
    type=str, 
    help="Select plot to be displayed(")
parser.add_argument("--output", default=None, 
    type=str, 
    help="Input file with the results.")


# Used colors in the plots.
COLORS = ['#d7191c', '#2c7bb6', '#fdae61', '#abd9e9', '#ffffbf', '#abd9e9']

# Input filenames.
RANDOM_FILE = "../analysis_results/results/random_results/results0.txt"
GREEDY_FILE = "../analysis_results/results/greedy_results/results0.txt"
GENETIC_FILE = "../analysis_results/results/genetic_results/results0.txt"
KMEANS_FILE = "../analysis_results/results/kmeans_results/results0.txt"


def ReadData(filename, first, second, specify_par=None):
    """
    Reads experiment results from the input file and filters choosen columns to
    prepare data for the plots.
    Args:
        filename: name of the input file (in format of result overview  
            from `metacentrum_grid_search.py` output)
        first: key of the first data column
        second: key of the second data column
        specify_par: parameter to filter only choosen values (currently not working)

    Returns prepared 1 item dictionary with 2D array of filtered results.
    """

    results = {}

    first_col, second_col = 0, 0
    specify_col = 0

    with open(filename, 'r') as fp:
        # Find indices of the choosen columns (first line is header).
        header = fp.readline().split()
        for i, key in enumerate(header):
            if key == first:
                first_col = i           
            elif key == second:
                second_col = i
            elif key == specify_par:
                specify_col = i
        
        # Skip header separator.
        fp.readline()

        # Take the choosen columns and create data for the plot.
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
    """
    Creates plot of the run down vehicles and number of stations (multiple variants). 
    """


    # First variant:
    first_column = "numStations"
    second_column = "runDownBat"


    # Second variant (uncomment to use it):
    # first_column = "numUsedStat"
    # second_column = "waitTime"


    results_greedy = ReadData(GREEDY_FILE, first_column, second_column)
    results_genetic = ReadData(GENETIC_FILE, first_column, second_column)
    results_kmeans = ReadData(KMEANS_FILE, first_column, second_column)
    results_random = ReadData(RANDOM_FILE, first_column, second_column)

    all_results = [results_greedy, results_genetic, results_kmeans]

    # Labels of the legend.
    labels = ["Hladový přístup - průměr", 
        "Genetický algoritmus - průměr", 
        "K-Means - průměr", 
        "Náhodný přístup - průměr"]
    scatter_labels = ["Hladový přístup - experimenty",
        "Genetický algoritmus - experimenty", 
        "K-Means - experimenty", 
        "Náhodný přístup - experimenty"]


    # Default font of the plot.
    font = {'size' : 13}
    matplotlib.rc('font', **font)

    
    # Plot size.
    fig, ax = plt.subplots()
    fig.set_figheight(7.1)
    fig.set_figwidth(9.15)
    

    # Plot title ((uncomment satisfactory).
    ax.set_title('Porovnání optimalizačních metod s ohledem na\npočet rozmístěných stanic (původní experimenty)',
        fontweight="bold",
        size=23)
    # ax.set_title('Porovnání optimalizačních metod s ohledem na\npočet rozmístěných stanic', 
    #     fontweight="bold", 
    #     size=23)
    # ax.set_title('Porovnání optimalizačních metod s ohledem\nna počet použitých stanic', 
    #     fontweight="bold", 
    #     size=23)


    # Plot axes labels (uncomment satisfactory).

    # x axis label
    ax.set_xlabel('Počet umístěných nabíjecích stanic', size=22)
    # ax.set_xlabel('Počet použitých nabíjecích stanic', size=22)

    # y axis label
    ax.set_ylabel('Počet vybitých vozidel', size=22)

    
    # Show optimization results.
    for i, results in enumerate(all_results):
        for _, item in enumerate(results.items()):
            # Prepare data for the plot.
            value = item[1]
            x = np.array([float(m) for m in value[0]])
            y = np.array([float(m) for m in value[1]])

            indices = np.argsort(x)
            x = x[indices]
            y = y[indices]

            # Compute mean values for the plot.
            unique_x = np.unique(x)
            y_means = copy.deepcopy(y)
            for un in unique_x:
                indices = np.where(x == un)
                mean_y = np.mean(y[indices])  
                y_means[indices] = mean_y

            
            # Display plots (uncomment satisfacory).
            
            # Plots of the mean.
            ax.plot(x, y_means, color=COLORS[i], label=labels[i])
            # ax.scatter(x, y_means, color=COLORS[i])#, label=scatter_labels[i])
            # Plots of the real values.
            
            # ax.plot(x, y, color=COLORS[i], label=labels[i])
            # ax.scatter(x, y, color=COLORS[i], label=scatter_labels[i])
            

    # Show random results
    for _, item in enumerate(results_random.items()):
            value = item[1]
            x = np.array([float(m) for m in value[0]])
            y = np.array([float(m) for m in value[1]])
            ax.plot(x, y, color=COLORS[3], linewidth=2, linestyle='dashdot', label=labels[3])
            ax.scatter(x, y, color=COLORS[3])#, label=scatter_labels[3])


    # Show the label legend.
    ax.legend(fontsize=16.2)

    # Decide whether save the plot or display it.
    if args.output != None:
        plt.savefig(args.output)
    else:
        plt.show()



def main(args):
    PlotRunDownStations(args)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)