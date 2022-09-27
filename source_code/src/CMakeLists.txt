
set(HEADER_LIST 
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/boost_header.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/Car.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/CarIterator_header.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/ChargingStation.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/Edge.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/global_parameters.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/GraphAdjuster.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/Map.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/MapPosition.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/MapReader.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/ModelRepresentation.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/Node.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/Optimizer.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/OptimizerParameters.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/SimulationParameters.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/StationParameters.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/TableEvent.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/TimeTable.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/TrafficSimulator.hpp
    ${TrafficSimulator_SOURCE_DIR}/include/optimizer/Vehicle.hpp
)

set(SOURCE_LIST
    ${TrafficSimulator_SOURCE_DIR}/src/Car.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/ChargingStation.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/City.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/Edge.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/GraphAdjuster.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/Map.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/MapPosition.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/MapReader.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/ModelRepresentation.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/Node.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/Optimizer.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/OptimizerParameters.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/SimulationParameters.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/StationParameters.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/TableEvent.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/TimeTable.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/TrafficSimulator.cpp
    ${TrafficSimulator_SOURCE_DIR}/src/Vehicle.cpp
)


# Make an automatic library - will be static or dynamic based on user setting
add_library(optimizer ${SOURCE_LIST} ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(optimizer PUBLIC ../include)

# This depends on (header only) boost
target_link_libraries(optimizer PRIVATE Boost::graph Boost::random)

# All users of this library will need at least C++11
target_compile_features(optimizer PUBLIC cxx_std_11)
target_compile_options(optimizer PUBLIC -Ofast)



# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST}
)
