# About
This directory contains program `metacentrum_grid_search.py` which prepares bash scripts that are 
used to run jobs with experiments on [Metacentrum](https://metavo.metacentrum.cz/). The program also
provides functionality to display overview of the experiment results.

Program `result_extractor.py` is used to prepare plots of the experiment results.


# Metacentrum Jobs Preparation
In order to prepare bash scripts which executes experiments as jobs on 
[Metacentrum](https://metavo.metacentrum.cz/) we use the program 
`metacentrum_grid_search.py` with flag `--writeResults` set to `False`.

## Simulator Parameters
Program reads parameters of the simulator for the experiments from the JSON files.
Keys are simulator parameters, values are arrays of all possible values of the parameter.

Simulator parameters should be splitted into 3 files by its purpose, these files are specified
by corresponding switch of the program (use `--help` switch for more info). Though it is not
neccesary to split the parameter we strongly recommend doing so. To see how to correctly
write parameters files see arbitrary `JSON` file from `analysis_results/results`
(see [here](https://github.com/dbeinhauer/bcs-source/tree/main/analysis_results/results)).

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