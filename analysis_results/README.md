# About
This directory is composed of `results/` subdirectory which contains organized system of `JSON` files which specify 
parameter choice for the experiments and `.txt` files with results of the 
experiments. Subdirectory `plots/` contains plots of the 
choosen results (note that the plots are in Czech) and subdirectory 
`models/` contains `.txt` files with choosen model 
descriptions and simple Python program for filtering the used stations.

# Organization of Results Subdirectory
Subdirectory `analysis_results/results/` contains for each used optimizer 
another subdirectory with coresponding results and parameter descriptions.
Each optimizer subdirectory contains `.txt` files with results of the 
experiments and subdirectories with prefix `parameters` containing 
`.json` files which specify choosen parameters for the corresponding 
experiment.

Structure of the content of the `analysis_results/results/` subdirectory:
```bash
├── genetic_results
│   ├── parameters0
│   │   ├── loss_parameters.json
│   │   ├── optimizer_parameters.json
│   │   └── simulation_parameters.json
│   ├── parameters1
│   │   ├── loss_parameters.json
│   │   ├── optimizer_parameters.json
│   │   └── simulation_parameters.json
│   ├── results0.txt
│   └── results1.txt
├── greedy_results
│   ├── parameters0
│   │   ├── loss_parameters.json
│   │   ├── optimizer_parameters.json
│   │   └── simulation_parameters.json
│   ├── parameters1
│   │   ├── loss_parameters.json
│   │   ├── optimizer_parameters.json
│   │   └── simulation_parameters.json
│   ├── results0.txt
│   └── results1.txt
├── kmeans_results
│   ├── parameters0
│   │   ├── loss_parameters.json
│   │   ├── optimizer_parameters.json
│   │   └── simulation_parameters.json
│   ├── parameters1
│   │   ├── loss_parameters.json
│   │   ├── optimizer_parameters.json
│   │   └── simulation_parameters.json
│   ├── results0.txt
│   └── results1.txt
├── random_results
│   ├── parameters0
│   │   ├── loss_parameters.json
│   │   └── simulation_parameters.json
│   ├── parameters1
│   │   ├── loss_parameters.json
│   │   └── simulation_parameters.json
│   ├── results0.txt
│   └── results1.txt
└── simulation_parameters_results
    ├── parameters0
    │   ├── loss_parameters.json
    │   └── simulation_parameters.json
    ├── parameters1
    │   ├── loss_parameters.json
    │   └── simulation_parameters.json
    ├── parameters2
    │   ├── loss_parameters.json
    │   └── simulation_parameters.json
    ├── results0.txt
    ├── results1.txt
    └── results2.txt
```

# Example
As an example we choose parameters and results of the experiments on 
genetic algorith located in the subdirectory `results/genetic_results/`.

The parameters of the run 0 of the experiments are located in the 
subdirectory `results/genetic_results/parameters0`. This subdirectory 
contains `.json` files with values of the parameters of the simulator.
Parameters are distributed by its purpose.

Results of the choosen experiment are stored in the file 
`results/genetic_results/results0.txt`.