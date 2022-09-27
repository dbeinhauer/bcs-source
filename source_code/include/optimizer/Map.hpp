#ifndef MAP_H_
#define MAP_H_

#include <iostream>
#include <vector>

#include <math.h>
#include <random>

#include "global_parameters.hpp"
#include "boost_header.hpp"

#include "ChargingStation.hpp"
#include "MapPosition.hpp"


//---------------------------------------------------------------------------
// Part of the boost A-Star algorithm, which is copied from the boost library 
// documentation with minor changes to adapt to our example.

// Euclidean distance heuristic - for A-Star algorithm.
template <class m_graph, class CostType, class LocMap>
class distance_heuristic : public boost::astar_heuristic<Graph, CostType>
{
public:
	distance_heuristic(LocMap l, vertex_t goal)
		: m_location(l), m_goal(goal) {}

	CostType operator()(vertex_t u)
	{
		CostType dx = m_location[m_goal].Longitude - m_location[u].Longitude;
		CostType dy = m_location[m_goal].Latitude - m_location[u].Latitude;
		return ::sqrt(dx * dx + dy * dy);
	}
private:
	LocMap m_location;
	vertex_t m_goal;
};


// Exception for termination the A-star algorothm computing.
struct found_goal {};

// Visitor that terminates when we find the goal.
template <class Vertex>
class astar_goal_visitor : public boost::default_astar_visitor
{
public:
	astar_goal_visitor(Vertex goal) : m_goal(goal) {}
	template <class Graph>
	void examine_vertex(Vertex u, Graph& g) {
		if (u == m_goal)
			throw found_goal();
	}
private:
	Vertex m_goal;
};

//-----------------------------------------------------------------------------


/// <summary>
/// Class to store and do basic operation with map objects.
/// </summary>
class Map
{
public:
	Map(int32_t numClosestStations, double segmentLength);

	void InitHeuristicLocations();

	void ResetSimulation();


	// Map setup:
	void AddNode(int32_t newID, double latitude, double longitude,
		int32_t cityID, double distance, int32_t oldID = -1);
	void AddEdge(int32_t firstID, int32_t secondID, char type, double length);
	void AddEdge(vertex_t firstID, vertex_t secondID, PermitedRoadTypes type, double length);
	void AddCity(int32_t cityID, int32_t population);
	void AddChargingStation(ChargingStation chargingStation);
	void SetCitiesDistances(Double2DMatrix citiesDistances);


	// Get information from the map:
	const Graph& GetConstGraph();
	Graph& GetGraph();
	const Node& GetNode(int32_t nodeID);
	std::vector<std::unique_ptr<ChargingStation>>& GetChargingStations();
	std::unique_ptr<ChargingStation>& GetChargingStation(int32_t stationID);
	double GetCityPairDistance(int32_t firstCityID, int32_t secondCityID);
	const std::vector<int32_t>& GetCitiesPopulation();
	const std::vector<vertex_t>& GetCityNodes(int32_t cityID);
	std::unique_ptr<ChargingStation>& GetLastChargingStation();
	std::map<edge_t, std::vector<BateryPair>&> GetAllSegments();

	bool ExistNode(int32_t nodeID);


	// Map operations:
	double ComputeSphericalDistance(double firstLatitude, double firstLongitude,
		double secondLatitude, double secondLongitude);

	void FindCityNodes();
	void FindCitiesNearestChargingStations();
	void FindVerticesNearestChargingStations();

	MapPosition FindClosestPositionToCoordinates(double latitude, double longitude);

	LengthPathPair FindVehiclePath(vertex_t start, vertex_t end, bool goCharging);
	LengthPathPair FindShortestPath(const vertex_t& start, const vertex_t& goal);
	std::vector<double> DijkstraFindDistances(vertex_t start);

	void TestComponents(std::tuple<int32_t, std::vector<vertex_t>>& numComponentsTuple);
	void SimulateVehiclePassedThroughEdge(edge_t edge, 
		std::pair<double, double> batteryPair, int32_t startSegment,
		bool segmentIncrease, int32_t endSegment);

private:
	// Random generator variables:
	std::random_device rd;
	std::mt19937 gen;

	// Number of closest charging station to check.
	int32_t numClosestStations_;
	// Length of the edge segment.
	double segmentLength_;
	// Graph representing the road network.
	Graph graph_;
	
	// Vector of populations of each city (index is city ID).
	std::vector<int32_t> cityPopulations_;
	// Vector of vector of nodes which belongs to the corresponding city (index is city ID).
	std::vector<std::vector<vertex_t>> cityNodes_;

	// Distances between each pair of the cities.
	Double2DMatrix citiesDistances_;

	// Vector of all charging stations objects placed on the map.
	std::vector<std::unique_ptr<ChargingStation>> allChargingStations_;
	
	// Vector of the given number of closest charging stations IDs for each City.
	std::vector<std::vector<int32_t>> closestStations_;
	std::vector<int32_t> vertexClosestStations_;

	// ID of the charging station where leads the last computed path 
	// (to eventualy set to the proper vehicle).
	int32_t lastClosestStation_;

	// Locations of each vertex in the map (for A-Star heuristic).
	std::vector<location> heuristicLocations_;


	std::vector<int32_t> findNearestChargingStations(int32_t cityID);
	int32_t findNearestStationForVertex(vertex_t vertexID);
	std::pair<vertex_t, double> findClosestVertex(
		double latitude, double longitude, vertex_t from, bool wholeMap);
};

#endif