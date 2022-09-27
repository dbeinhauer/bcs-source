cmake_minimum_required(VERSION 3.1...3.23)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# Setup the project
project(TrafficSimulator
    VERSION 0.1
    DESCRIPTION "Simulate traffic and optimize charging station placement"
    LANGUAGES CXX)

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add the boost library
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package( Boost REQUIRED COMPONENTS program_options graph random )

# Compiled library code
add_subdirectory(src)

# Executable code
add_subdirectory(apps)

target_compile_features(Traffic_Simulator PUBLIC cxx_std_17)