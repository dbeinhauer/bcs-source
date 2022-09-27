#ifndef GRAPHADJUSTER_H_
#define GRAPHADJUSTER_H_

#include <vector>
#include <cmath>

#include "global_parameters.hpp"

#include "Map.hpp"


// Pair of distance and vertex 
// (typicaly for storing distance of the vertex from the given city).
typedef std::pair<double, vertex_t> City_vertex_distance;
// Map of vertices and its distances of all cities.
typedef std::unordered_map<int32_t, City_vertex_distance> Component_city_vertex;
// Vector of all distances of all vertices from all cities.
typedef std::vector<Component_city_vertex> Graph_components;


/// <summary>
/// Struct to store city coordinates.
/// </summary>
struct City
{
public:
	// City center latitude and longitude.
	double Lat;
	double Lon;

	City(double lat, double lon)
		: Lat(lat), Lon(lon) {}
};


/// <summary>
/// Struct to store necessary info about the pair of vertices.
/// </summary>
struct ComponentsDistance
{
public:
	double Distance;
	vertex_t FirstVertexID;
	vertex_t SecondVertexID;

	ComponentsDistance()
	{
		this->Distance = 0;
		this->FirstVertexID = 0;
		this->SecondVertexID = 0;
	}

	ComponentsDistance(double distance, vertex_t firstID, vertex_t secondID)
		: Distance(distance), FirstVertexID(firstID), SecondVertexID(secondID) {}
};


/// <summary>
/// Struct to store all necessary info about edge.
/// </summary>
struct AddedEdge
{
public:
	vertex_t FirstVertex;
	vertex_t SecondVertex;
	char Type;
	double Length;

	AddedEdge(vertex_t first, vertex_t second, char type, double length)
		: FirstVertex(first), SecondVertex(second), Length(length), Type(type)
	{}
};


/// <summary>
/// Class to do preprocessing of the graph.
/// </summary>
class GraphAdjuster
{
public:
	GraphAdjuster();

	void AddCityCoordinate(double lat, double lon);

	std::vector<City>& GetCitiesCoordinates();
	std::vector<AddedEdge> GetAddedEdges();
	std::vector< std::vector<double>> GetCitiesDistances();

	void ComputeCitiesDistances();

	bool MergeTwoComponents(Map& map, int32_t numComponents, std::vector<vertex_t>&& verticesComponents);
	bool MergeDegreeTwoVertices(Map& map);
	bool MergeDegreeTwoVertex(vertex_t vertex_to_merge, Map& map);


private:
	std::vector<AddedEdge> addedEdges_;

	std::vector<City> citiesCoordinates_;
	std::vector<std::vector<double>> citiesDistances_;
	//Graph_components graphComponents_;

	//void computeCitiesDistances();
	double computeSphericalDistance(const City& firstCity, const City& secondCity);
	Graph_components computeCityVertexPairs(const Graph& graph, int numComponents, std::vector<vertex_t>& verticesComponents);
	double computeCitiesComponentDistance(const std::pair<int, int>& citiesID, double firstVertexDistance, double secondVertexDistance);

	Component_city_vertex findComponentCityVertexPairs(const Graph& graph, const std::vector<vertex_t>& component);
	ComponentsDistance findClosestCitiesPair(int firstComponentID, int secondComponentID, const Graph_components& graph_components);

};

#endif