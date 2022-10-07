# About
This directory contains program `metacentrum_grid_search.py` which prepares bash scripts that are 
used to run jobs with experiments on [Metacentrum](https://metavo.metacentrum.cz/). The program also
provides functionality to display overview of the experiment results.

Program `result_extractor.py` is used to prepare plots of the experiment results.


# Experiment Preparation
In order to prepare bash scripts which executes experiments as jobs on 
[Metacentrum](https://metavo.metacentrum.cz/) we use the program 
`metacentrum_grid_search.py` with flag `--writeResults` set to `False`.

## Simulator Parameters
Program reads parameters for the simulator for the choosen experiments from 
the given JSON files.
The keys of the files are simulator parameters, the values are 
arrays of all possible values of the parameter.

Simulator parameters should be splitted into 3 files by its purpose, these files are specified
by corresponding switch of the program. The switches are the following:

```bash
--simulationParamFile SIMULATIONPARAMFILE
                      Filename of the simulation parameters
--lossParamFile LOSSPARAMFILE
                      Filename of the loss parameters
--optimizerParamFile OPTIMIZERPARAMFILE
                      Filename of the optimizer parameters
```

Though it is not neccesary to split the parameter we strongly recommend 
doing so. To see how to correctly write parameters files see 
arbitrary `JSON` file from `analysis_results/results`
(see [here](https://github.com/dbeinhauer/bcs-source/tree/main/analysis_results/results)) 
which can serve as an example of the correct format.

If more than 1 parameter value is specified then all combinations of the 
parameters values are used for the experiments. 

## Metacentrum Jobs Properties
In order to specify metacentrum jobs properties there can be used parameters 
from below.

```
--jobPrefix JOBPREFIX
                      Prefix of the job name.
--outputFilePrefix OUTPUTFILEPREFIX
                      Prefix of the job name.
--outputDir OUTPUTDIR
                      Output directory of the commands.
--resourceSelection RESOURCESELECTION
                      What resources ask for the job.
--walltime WALLTIME   How long should the job last.
--dataDir DATADIR     What is the data directory.
--copyCommand COPYCOMMAND
                      What files copy to scratch directory.
--optimizer OPTIMIZER
                      Which optimizer choose (`randomModels`, `greedy`,
                      `genetic`, `kMeans`)
--programOutputPath PROGRAMOUTPUTPATH
                      Where the output of the program will be stored.
```

### Multiple Same Experiments
The program also offers possibility to create jobs that executes multiple
experiments with the same parameters setup. In order to use this functionality
it is neccesary to specify exactly 1 value for each of the simulator parameters.
in the parameters files mentioned above. The switches to adjust properties of 
this functionality are:

```
--repeatCommand REPEATCOMMAND
                      Whether to create same command multiple times (for
                      statistical research).
--repeatNum REPEATNUM
                      Number of repetitions of the same command creation.
```

### Output Filenames and Result Overview
The prepared jobs store the experiment outputs in the files with names in 
the specific format which allows the program write experiment results in number
of different forms. If you want to use this functionality, never rename the 
output files of the experiments.


## Example
The following command gives an example of the program usage to generate
metacentrum jobs.

```
./metacentrum_grid_search.py \
--simulationParamFile=simulation_parameters.json \
--lossParamFile=loss_parameters.json\
--optimizerParamFile=optimizer_parameters.json \
--optimizer=genetic \
--outputDir=prepared_commands/prepared_genetic_search/ \
--programOutputPath='$DATADIR/Simulator_outputs/genetic_search' \
--resourceSelection=select=1:ncpus=8:mem=4gb:scratch_local=30gb \
--walltime=walltime=48:00:00
```

# Results Analysis
The program `metacentrum_grid_search.py` also offers functionality
to write a structured overview of the specific experiment results.

## Required Experiment Properties
The result overview is only possible when the following is fullfiled:
1. All of the examined simulator results are from the same optimizer.
2. All result files are located in the same directory with no other redundant
files.
3. None of the file with the results was renamed.
4. The possibility of the multiple same experimets was not used.

## Results Preprocessing
The program takes as an input one file in the following format of the text:

```
==> <output_filename> <==
Number of stations: <NUM>
Number of run down batteries: <NUM>
Average traveling duration: <NUM> minutes
Average battery difference between start and end: <NUM>
Average waiting times in charging station: <NUM> minutes
Number of vehicles which reached the end: <NUM>
Number of vehicles which didn't finish before end of simulation: <NUM>
Total number of used charging stations: <NUM>
Total number of charging stations usage: <NUM>
Run down batteries during going to station: <NUM>
Total number of vehicles which didn't reach the end: <NUM>
------------------------------------------------------
Model loss: <NUM>

```

This text represents the results of the one experiment. The input file can be
structured of an arbitrary number of text sequences structured like the example above.

The structure of the input is in the output format of the bash command 

```bash
tail -n14
```

The recommended way to prepare the input file is the following:
1. Store all the output files with the results to be examined in the 
empty directory.
2. In this directory run command:

```bash
tail -n14 * >results_input.txt
```

3. The prepared input file is stored in the file `results_input.txt`.

## Result analysis usage
In order to write experiments overview it is necessary to set the flag
`--writeResults` to `True`. Then it is necessary to specify files with the used 
simulator parameters in the same format as when creating the experiments
(for more info see the section Experiment Preparation above).

It is also neccesary specify input file with the preprocessed results mentioned
in the section Results Preprocessing above). 
In addition it is possible to specify rigit value of some parameters to
reduce the number of results and to analyse only choosen properties of the
experiments. It is also possible to specify by which parameter should be the 
result table sorted. The parameters mentioned in this paragraph are the 
following:

```
--inputResults INPUTRESULTS
                      Input file with the results.
--rigitParameters RIGITPARAMETERS [RIGITPARAMETERS ...]
                      Specify some rigit parameters which we are interested
                      in (we expect all values are also specified in
                      `rigitValues`).
--rigitValues RIGITVALUES [RIGITVALUES ...]
                      Values of the specified rigit parameters.(in the same
                      order as `rigitParameters`).
--sortBy SORTBY       By which parameter sort the results.
```

The output of the command is the table of the choosen experiment results
with the columns which specify given experiment parameters (if not rigit)
or values of the examined property and rows which specify exact experiments.

## Example
The following command gives an example of how to use the program to write 
experiment results:

```bash
./metacentrum_grid_search.py \
--simulationParamFile=simulation_parameters.json \
--lossParamFile=loss_parameters.json\
--optimizerParamFile=optimizer_parameters.json \
--optimizer=genetic \
--writeResults=True \
--inputResults=Simulator_outputs/tails_genetic_search.txt \
--rigitParameters simulationTime \
--rigitValues 1000 \
--sortBy=numStations
```

The following could be the output of the experiment results overview:

```bash
numStations    runDownBat    travelDur    batteryDiff    waitingTime     bestLoss    numStations
-------------  ------------  -----------  -------------  -------------  -----------  -------------
        4000        586865      93.7943      -0.312037        163.055  5.90868e+07           4000
        2000        590083      93.1925      -0.316216        161.792  5.92086e+07           2000
        1000        631222      95.7939      -0.314245        166.954  6.32225e+07           1000
        500        687167      95.9103      -0.31716         174.007  6.8767e+07             500
        300        753241      80.9886      -0.317146        182.513  7.53544e+07            300
```

# Plot Preparation
In order to plot the results of the experiments it is possible to
use the program `results_extractor.py`. Note that the program is currently in 
the provisory state and different plots can be choosen only by commenting
and uncommenting the appropriate parts of the code. For more information
about the used of the program use the command:

```bash
./result_extractor.py --help
```

# Příprava úloh metacentra a analýza výsledků simulací

Adresář obsahuje pomocný program pro přípravu úloh metacentra pro
vybrané parametry zapsané ve formátu JSON. 
Program zároveň umožňuje vypisovat výsledky experimentů
ve formě tabulky.

## Požadavky
Program vyžaduje knihovny Numpy a tabulate.

## Příklad použití:

Vytváření úloh pro metacentrum:
```
./metacentrum_grid_search.py \
--simulationParamFile=simulation_parameters.json \
--lossParamFile=loss_parameters.json\
--optimizerParamFile=optimizer_parameters.json \
--optimizer=genetic \
--outputDir=prepared_commands/prepared_genetic_search/ \
--programOutputPath='$DATADIR/Simulator_outputs/genetic_search' \
--resourceSelection=select=1:ncpus=8:mem=4gb:scratch_local=30gb \
--walltime=walltime=48:00:00
```

Výpis výsledků experimentů:
```
./metacentrum_grid_search.py \
--simulationParamFile=simulation_parameters.json \
--lossParamFile=loss_parameters.json\
--optimizerParamFile=optimizer_parameters.json \
--optimizer=genetic \
--writeResults=True \
--inputResults=Simulator_outputs/tails_genetic_search.txt \
--rigitParameters simulationTime \
--rigitValues 1000 \
--sortBy=numStations
```

## Příklad výpisu výsledků
```
  numStations    runDownBat    travelDur    batteryDiff    waitingTime     bestLoss    numStations
-------------  ------------  -----------  -------------  -------------  -----------  -------------
         4000        586865      93.7943      -0.312037        163.055  5.90868e+07           4000
         2000        590083      93.1925      -0.316216        161.792  5.92086e+07           2000
         1000        631222      95.7939      -0.314245        166.954  6.32225e+07           1000
          500        687167      95.9103      -0.31716         174.007  6.8767e+07             500
          300        753241      80.9886      -0.317146        182.513  7.53544e+07            300
```