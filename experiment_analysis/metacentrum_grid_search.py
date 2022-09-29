#!/usr/bin/env python3

import sys
import itertools
import argparse

import json
import numpy as np

from tabulate import tabulate


parser = argparse.ArgumentParser()

parser.add_argument("--simulationParamFile", default="simulation_parameters.json", 
    type=str, help="Filename of the simulation parameters")
parser.add_argument("--lossParamFile", default="loss_parameters.json", 
    type=str, help="Filename of the loss parameters")
parser.add_argument("--optimizerParamFile", default="optimizer_parameters.json", 
    type=str, help="Filename of the optimizer parameters")

parser.add_argument("--jobPrefix", default="trafficSimulator_", 
    type=str, help="Prefix of the job name.")
parser.add_argument("--outputFilePrefix", default="trafficSimulator_output_", 
    type=str, help="Prefix of the job name.")
parser.add_argument("--outputDir", default="prepared_commands/", 
    type=str, help="Output directory of the commands.")
parser.add_argument("--resourceSelection", default="select=1:ncpus=4:mem=8gb:scratch_local=30gb", 
    type=str, help="What resources ask for the job.")
parser.add_argument("--walltime", default="walltime=12:00:00",
    type=str, help="How long should the job last.")
parser.add_argument("--dataDir", default="/storage/praha1/home/$(whoami)", 
    type=str, help="What is the data directory.")
parser.add_argument("--copyCommand", default="-r $DATADIR/Traffic_Simulator/*",
    type=str, help="What files copy to scratch directory.")
parser.add_argument("--optimizer", default="randomModels", 
    type=str, help="Which optimizer choose (`randomModels`, `greedy`, `genetic`, `kMeans`)")
parser.add_argument("--programOutputPath", default="$DATADIR/Simulator_outputs", 
    type=str, help="Where the output of the program will be stored.")
# parser.add_argument("--saveBest", default=False, 
#     type=bool, help="Whether save the best model to file.")

parser.add_argument("--repeatCommand", default=False, 
    type=bool, help="Whether to create same command multiple times (for statistical research).")
parser.add_argument("--repeatNum", default=10,
    type=int, help="Number of repetitions of the same command creation.")

parser.add_argument("--writeResults", default=False, 
    type=bool, help="Whether write grid search results.")
parser.add_argument("--inputResults", default="Simulator_outputs/correct_outputs.txt", 
    type=str, help="Input file with the results.")
parser.add_argument("--rigitParameters", nargs='+', type=str, 
    help="Specify some rigit parameters which we are interested in\
         (we expect all values are also specified in `rigitValues`).")
parser.add_argument("--rigitValues", nargs='+', type=str, 
    help="Values of the specified rigit parameters.(in the same order as `rigitParameters`).")
parser.add_argument("--sortBy", default="bestLoss", type=str,
        help="By which parameter sort the results.")


# String constants:
FILE_POSTFIX = ".sh"
PROGRAM_COMMAND ="./Traffic_Simulator"
EDGE_FILE = "prepared_graph/prepared_edges.txt"
NODE_FILE = "prepared_graph/prepared_combined_nodes.txt" 
CITY_FILE = "prepared_graph/cities.txt"

# Parameters dictionaries:
SIMULATION_PARAMS = {}
LOSS_PARAMS = {}
OPTIMIZER_PARAMS = {}


def loadDictionary(filename):
    """
    Loads dictionary from JSON file.
    """
    with open(filename, 'r') as handle:
        json_data = json.load(handle)
    
    return json_data


def getPermutations(params, skip_prompt):
    """
    Creates permutations of all specified parameter values.
    """
    keys, values = zip(*params.items())
    permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]

    print(f"There will be {str(len(permutations_dicts))} instances.")
    if not skip_prompt:
        print("Do you want to continue? [y/n]")

        answer = None
        while answer != 'y' and answer != 'n':
            if answer != None:
                print("Wrong answer. Asking again: Do you want to continue? [y/n]")
            answer = input()

        if answer == 'n':
            sys.exit(0)

    return permutations_dicts



def PrepareParameters(args, params, commandPrefix, jobPrefix, outputPrefix):
    """
    Prepares parameters for all permutations and writes all scripts
    into its proper files.
    Args:
        args: command line arguments
        params: all specified parameters of the simulator
        commandPrefix: prefix of the command to start the simulator.
        jobPrefix: prefix of the metacentrum job name
        outputPrefix: prefix of the output file
    """

    permutations_dicts = {}

    if not args.repeatCommand:
        permutations_dicts = getPermutations(params, False)
    else:
        permutations_dicts = getPermutations(params, True)
    

    for _, param_dict in enumerate(permutations_dicts):
        # For each combination of parameters create command a job.

        outputFilename = outputPrefix
        jobName = jobPrefix
        command = commandPrefix

        filename = args.outputDir
        for key, val in param_dict.items():
            # Prepare command and proper names.

            command += " --" + str(key) + "=" + str(val)
            if not args.repeatCommand:
                # Filename specifies which parameters were choosen.
                nameFormat = "_" + str(key[0:2]) + str(val)
                outputFilename += nameFormat
                jobName +=nameFormat

        outputFilename += ".txt"
        command += " >" + outputFilename
        filename += jobName + FILE_POSTFIX

        # Write prepared metacentrum job into proper file.
        with open(filename, 'w') as f:
            f.write(WriteHeader(args, jobName))
            f.write(WriteDataPreparation(args))
            f.write(command + "\n")
            f.write(WriteTail(args, outputFilename))



def WriteHeader(args, name):
    """
    Returns string representation of the script header.
    Args:
        name: job name
    """

    return "\n".join(
        [
            "#!/bin/bash",
            f"#PBS -N {str(name)}",
            f"#PBS -l {args.resourceSelection}",
            f"#PBS -l {args.walltime}",
            "#PBS -m ae",
            "\n"
        ])


def WriteDataPreparation(args):
    """
    Returns string representation of the script data preparation.
    """

    return "\n".join(
        [
            f"DATADIR={args.dataDir}",
            "\n",
            "echo \"$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR\" >> $DATADIR/jobs_info.txt",
            "\n",
            "test -n \"$SCRATCHDIR\" || { echo >&2 \"Variable SCRATCHDIR is not set!\"; exit 1; }",
            "\n",
            "cp " + args.copyCommand + " $SCRATCHDIR || { echo >&2 \"Error while copying input file(s)!\"; exit 2; }",
            "\n",
            "cd $SCRATCHDIR",
            "\n",
            "g++ -Ofast *.cpp -o Traffic_Simulator -lboost_program_options -lstdc++",
            "\n"
        ])


def PrepareCommand(optimizer):
    """
    Returns string representation of prefix of the main command for all scripts.
    """

    # best_model_part = ""

    return " ".join(["\n"
        f"{PROGRAM_COMMAND}",
        f" --edgesFile={EDGE_FILE}",
        f"--nodeCityFile={NODE_FILE}",
        f"--citiesFile={CITY_FILE}",
        f"--{optimizer}",
        "--saveBestModel",
        f"--modelFile={optimizer}.txt"
        ])


def WriteTail(args, outputFile):
    """
    Returns string representation of the tail of the script.
    """
    return "\n".join(
        [
            "cp " + outputFile + " "  + args.programOutputPath + " || { echo >&2 \"Result file(s) copying failed (with a code $?) !!\"; exit 4; }",
            "cp " + f"{args.optimizer}.txt" + " "  + args.programOutputPath + " || { echo >&2 \"Result file(s) copying failed (with a code $?) !!\"; exit 4; }",
            "\n",
            "clean_scratch"
        ])



def ReadResults(args):
    """
    Loads results of the optimizations. 
    Input file is in format of output of command:
            `tail -n10 dir_of_all_correct_simulation_results\*`
    We expect header for each simulation is in format:
    `==> filename`
    Which is followed by 10 last lines of the Traffic_Simulator output.

    Returns:
        Returns optimization results.
    """

    # Results of the optimizations.
    filenames = []
    numStations = []
    runDownBatteries = []
    travelDurations = []
    batteryDifferences = []
    waitingTimes = []
    numFinished = []
    numNotFinished = []
    numUsedStations = []
    totalNotFinished = []

    bestLosses = []

    with open(args.inputResults, 'r') as fp:
        # Parse input results in format:
        for line in fp:
            if line.startswith("==>"):
                filenames.append(line.split(" ")[1].split("__")[2].removesuffix(".txt"))
            
            elif line.startswith("Number of stations:"):
                numStations.append(int(line.split(":")[1]))

            elif line.startswith("Number of run down batteries:"):
                runDownBatteries.append(int(line.split(":")[1]))
            
            elif line.startswith("Average traveling duration:"):
                travelDurations.append(float(line.split(":")[1].split(" ")[1]))
            
            elif line.startswith("Average battery difference between start and end:"):
                batteryDifferences.append(float(line.split(":")[1]))

            elif line.startswith("Average waiting times in charging station:"):
                waitingTimes.append(float(line.split(":")[1].split(" ")[1]))

            elif line.startswith("Number of vehicles which reached the end:"):
                numFinished.append(int(line.split(":")[1]))

            elif line.startswith("Number of vehicles which didn't finish before end of simulation:"):
                numNotFinished.append(int(line.split(":")[1]))

            elif line.startswith("Total number of used charging stations:"):
                numUsedStations.append(int(line.split(":")[1]))


            elif line.startswith("Total number of vehicles which didn't reach the end:"):
                totalNotFinished.append(int(line.split(":")[1]))

            elif line.startswith("Model loss:"):
                bestLosses.append(float(line.split(":")[1]))




    return filenames, numStations, runDownBatteries, travelDurations, batteryDifferences, waitingTimes, numFinished, \
        numNotFinished, numUsedStations, totalNotFinished, bestLosses


def WriteResults(args, params):
    """
    Writes results of the simulations in the specified format from the `args`.
    """
    filenames, numStations, runDownBatteries, travelDurations, batteryDifferences, waitingTimes, numFinished, \
        numNotFinished, numUsedStations, totalNotFinished, bestLosses = ReadResults(args)

    allResults = np.array(list(zip(numStations, runDownBatteries, travelDurations, batteryDifferences, waitingTimes, numFinished, \
        numNotFinished, numUsedStations, totalNotFinished, bestLosses)))

    # Get rigit parameters:
    rigitParams = {}
    if args.rigitParameters != None:
        rigitParams = dict(zip(args.rigitParameters, args.rigitValues)) 


    def testNonAlphabetChar(c):
        return not str.isalpha(c)

    indicesToChoose = []
    choosenParameterValues = []
    choosenParameterNames = []
    for i, filename in enumerate(filenames):
        # Search all models

        parametersPart = filename.split('_')
        parameterValues = [''.join(filter(testNonAlphabetChar, x)) for x in parametersPart]
        parameterValuesToPrint = []
        parametersCorrect = True

        for j, key in enumerate(params.keys()):
            # Search all parameters of the commands

            if key in rigitParams:
                if parameterValues[j] != rigitParams[key]:
                    # Rigit parameter with wrong value -> skip results
                    parametersCorrect = False

                # Rigit parameter with correct value -> continue with search
                continue

            elif len(params[key]) == 1:
                # Parameter has only 1 possible value -> don't write it to the table
                continue

            elif i == 0:
                # Save parameter name.
                choosenParameterNames.append(key)
            
            parameterValuesToPrint.append(parameterValues[j])

        if parametersCorrect:
            # Parameters passed the filter -> add it to the table list
            indicesToChoose.append(i)
            choosenParameterValues.append(parameterValuesToPrint)


    choosenParameterNames = np.array(choosenParameterNames)
    choosenParameterValues = np.array(choosenParameterValues)

    # Take only filtered values.
    allResults= np.take(allResults, indicesToChoose, axis=0)

    # Concatenate parameters with loss info.
    choosenResults = np.concatenate((allResults, choosenParameterValues), axis=1).astype(float)


    # Select by which parameter should we sort:

    offset = 0
    lossOffset = 9

    # By some loss info:
    if args.sortBy == "numStat":
        choosenResults = choosenResults[choosenResults[:, offset].argsort()]

    elif args.sortBy == "runDownBat":
        choosenResults = choosenResults[choosenResults[:, offset + 1].argsort()]

    elif args.sortBy == "travDur":
        choosenResults = choosenResults[choosenResults[:, offset + 2].argsort()]

    elif args.sortBy == "battDiff":
        choosenResults = choosenResults[choosenResults[:, offset + 3].argsort()]

    elif args.sortBy == "waitTime":
        choosenResults = choosenResults[choosenResults[:, offset + 4].argsort()]

    elif args.sortBy == "numFinish":
        choosenResults = choosenResults[choosenResults[:, offset + 5].argsort()]

    elif args.sortBy == "numNotFinish":
        choosenResults = choosenResults[choosenResults[:, offset + 6].argsort()]

    elif args.sortBy == "numUsedStat":
        choosenResults = choosenResults[choosenResults[:, offset + 7].argsort()]
    
    elif args.sortBy == "notFinis":
        choosenResults = choosenResults[choosenResults[:, offset + 8].argsort()]

    elif args.sortBy == "bestLoss":
        choosenResults = choosenResults[choosenResults[:, offset + 9].argsort()]

    else:
        # By some simulator parameter:
        index = np.where(choosenParameterNames == args.sortBy)
        choosenResults = choosenResults[choosenResults[:, index[0][0] + offset + lossOffset + 1].argsort()]

    # Create table header:
    header = np.array(["numStat", "runDownBat", "travDur", "battDiff", "waitTime", "numFinish",
        "numNotFinish", "numUsedStat", "notFinis", "bestLoss"])
    header = np.concatenate((header, choosenParameterNames), axis=0)

    # Print table of results
    print(tabulate(choosenResults, headers=header))



def main(args):
    # Load parameters:
    SIMULATION_PARAMS = loadDictionary(args.simulationParamFile)
    LOSS_PARAMS = loadDictionary(args.lossParamFile)

    # Prepare variables
    params = SIMULATION_PARAMS | LOSS_PARAMS

    jobPrefix = args.jobPrefix
    commandPrefix = ""
    filenamePrefix = args.outputDir
    outputPrefix = args.outputFilePrefix

    # Prepare correct optimizer parameters:
    if args.optimizer == "greedy":
        OPTIMIZER_PARAMS = loadDictionary(args.optimizerParamFile)
        params = params | OPTIMIZER_PARAMS
        jobPrefix += "_greedy_"
        commandPrefix = PrepareCommand(args.optimizer)
        outputPrefix += "_greedy_"
        filenamePrefix += "_greedy_"

    elif args.optimizer == "genetic":
        OPTIMIZER_PARAMS = loadDictionary(args.optimizerParamFile)
        params = params | OPTIMIZER_PARAMS
        jobPrefix += "_genetic_"
        commandPrefix = PrepareCommand(args.optimizer)
        outputPrefix += "_genetic_"
        filenamePrefix += "_genetic_"

    elif args.optimizer == "kMeans":
        OPTIMIZER_PARAMS = loadDictionary(args.optimizerParamFile)
        params = params | OPTIMIZER_PARAMS
        jobPrefix += "_kMeans_"
        commandPrefix = PrepareCommand(args.optimizer)
        outputPrefix += "_kMeans_"
        filenamePrefix += "_kMeans_"
 
    else:
        # Random Values
        jobPrefix += "_simulation_"
        commandPrefix = PrepareCommand("randomModels")
        outputPrefix += "_simulation_"
        filenamePrefix += "_simulation_"


    if args.writeResults:
        # Write results of the grid search
        WriteResults(args, params)

    else:
        # Prepare scripts for the grid search
        if args.repeatCommand:
            for i in range(args.repeatNum):
                jobPrefix_new = jobPrefix + f"_{i}_"
                outputPrefix_new = outputPrefix + f"_{i}_"
                PrepareParameters(args, params, commandPrefix, jobPrefix_new, outputPrefix_new)

        else:
            PrepareParameters(args, params, commandPrefix, jobPrefix, outputPrefix)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)