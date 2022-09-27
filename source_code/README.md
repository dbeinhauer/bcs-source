The requirements are:
- CMake 3.11 or better; 3.14+ highly recommended.
- A C++11 compatible compiler
- The Boost libraries
- Git

To configure:
```bash
cmake -S . -B build
```
Add `-GNinja` if you have Ninja.

To build:
```bash
cmake --build build
```

The final executable will be found in ``build/apps/traffic_simulator``