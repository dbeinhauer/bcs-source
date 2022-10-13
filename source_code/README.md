# Compilation
The requirements are:
- CMake 3.11 or better; 3.14+ highly recommended.
- A C++14 compatible compiler
- The Boost libraries

To configure:
```bash
cmake -S . -B build
```
Add `-GNinja` if you have Ninja.

To build:
```bash
cmake --build build
```

The final executable will be found in ``build/apps/Traffic_simulator``

# Run Simulator
To run the traffic simulator without the additional parameter 
specification (in default setup) copy the compiled program to 
the root directory of the repository with the following command:

```bash
cp build/apps/Traffic_simulator ..
```

Then just move to the root directory and run the simulator with the
specified optimizer swich or `--help` switch to get more info about 
the program parameters.

## Example
Example of the program execution in the root directory of the 
repository with choosen `--genetic` optimizer:

```bash
./Traffic_Simulator --genetic
```

# Documentation
Documentation of the program can be found either in the file 
[Traffic_Simulator.pdf]() or in the text of the thesis in the 
repository [bsc-thesis]() (note that this text is only in Czech).