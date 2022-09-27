#ifndef MAPREADER_H_
#define MAPREADER_H_


#include <iostream>
#include <fstream>
#include <sstream>

#include "Map.hpp"
#include "GraphAdjuster.hpp"
#include "SimulationParameters.hpp"
#include "ModelRepresentation.hpp"


class MapReader
{
public:
	MapReader();
	
	
	// Map Loading:
	bool LoadMap(Map& map,
		SimulationParameters& simulationParameters);

	bool LoadNodesCities(Map& map, std::istream& nodeCityStream);
	bool LoadEdges(Map& map, std::istream& edgeStream);
	bool LoadCities(Map& map, std::istream& cityStream);
	bool LoadChargingStations(Map& map, std::istream& chargingStationsStream);


	// Map writing:
	bool WriteAddedEdges(std::ostream& addedEdgesStream);
	bool WritePreparedCombinedNodes(Map& map, std::ostream& preparedNodesStream);
	bool WritePreparedEdges(Map& map, std::ostream& preparedEdgesStream);
	bool WriteStations(ModelRepresentation& model, Map& map, std::ostream& stationsStream);


	// Get information used during loading.
	GraphAdjuster& GetGraphAdjuster();

private:
	// Error flag.
	bool error_ = false;

	// Object to preprocess the graph.
	GraphAdjuster graphAdjuster_;
	// Usage to degree 2 vetices elimination.
	int32_t lastFirstNode_;


	// Help loading functions.
	bool loadNodeCityPair(Map& map, const std::string& line);
	bool loadEdge(Map& map, const std::string& line);
	bool loadCity(Map& map, const std::string& line);
	bool loadChargingStation(Map& map, const std::string& line);
};

#endif