# Experiment analysis results
This directory contains organized system of `JSON` files which specify parameter choice for the experiments and `.txt` files with results of the experiments. 

# Content Organization
Directory contains subdirectory for each used optimizer. Each optimizer 
subdirectory contains `.txt` files with results of the experiments and 
subdirectories with prefix `parameters` containing `.json` files which specify 
choosen parameters for the corresponding experiment.

Structure of the content:
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
For example parameters and results of the experiments on genetic algorith are 
located in the subdirectory `genetic_results/`.

The parameters of the run 0 of the experiments are located in the subdirectory 
`genetic_results/parameters0`. This subdirectory contains `.json` files with
values of the parameters of the simulator. Parameters are distributed by its 
purpose.

Results of this experiment are store in the file `genetic_results/results0.txt`.