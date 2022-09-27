#include "optimizer/GraphAdjuster.hpp"


GraphAdjuster::GraphAdjuster() {}


/// <summary>
/// Adds City coordinates to the proper container.
/// </summary>
/// <param name="lat">Latitude of the city.</param>
/// <param name="lon">Longitude of the city.</param>
void GraphAdjuster::AddCityCoordinate(double lat, double lon)
{
	this->citiesCoordinates_.push_back(City(lat, lon));
}


/// <summary>
/// </summary>
/// <returns>Returns reference to coordinates of all cities.</returns>
std::vector<City>& GraphAdjuster::GetCitiesCoordinates()
{
	return citiesCoordinates_;
}


/// <summary>
/// Computes distances between all pairs of the cities.
/// </summary>
void GraphAdjuster::ComputeCitiesDistances()
{
	int32_t numCities = this->citiesCoordinates_.size();
	this->citiesDistances_ = std::vector<std::vector<double>>(numCities, std::vector<double>(numCities));

	for (size_t i = 0; i < numCities; i++)
	{
		for (size_t j = i; j < numCities; j++)
		{
			if (i == j)
			// Same city -> distance is 0
			{
				this->citiesDistances_[i][j] = 0;
				continue;
			}

			double distance = this->computeSphericalDistance(this->citiesCoordinates_[i], this->citiesCoordinates_[j]);

			this->citiesDistances_[i][j] = distance;
			this->citiesDistances_[j][i] = distance;
		}
	}
}


/// <summary>
/// Finds closest component from the component with index 0 and merges them. 
/// Distance is computed with the closest cities of each component and 
/// distances of the closest vertices from the cities center.
/// Merges components with the new edge of type `OTHER` between the closest
/// vertices of the components with distance:
///			`distance = cities_distance + vertex1_city_distance + vertex2_city_distance`
/// </summary>
/// <param name="graph">Graph to work with.</param>
/// <param name="numComponents">Number of components of the graph.</param>
/// <param name="verticesComponents">Vector of IDs of each vertex from the component.</param>
/// <returns>Returns `true` if merging was sucessfull, else `false`.</returns>
bool GraphAdjuster::MergeTwoComponents(Map& map, int32_t numComponents, 
	std::vector<vertex_t>&& verticesComponents)
{
	// Graph has already merged all components.
	if (numComponents <= 1)
		return false;

	Graph_components graphComponents = 
		this->computeCityVertexPairs(map.GetGraph(), numComponents, verticesComponents);

	// ID of the closest component from component with id 0.
	int32_t closestComponentID = 0;
	ComponentsDistance closestPair;


	for (size_t i = 1; i < numComponents; i++)
	// Find nearest component from the first component in the container and its nearest vertices.
	{
		ComponentsDistance citiesPair = this->findClosestCitiesPair(0, i, graphComponents);

		// Error during finding closest pair -> end computing.
		if (citiesPair.FirstVertexID == -1)
			return false;

		if (citiesPair.Distance < closestPair.Distance || closestComponentID == 0)
		// First pair of components tested or current components are closer then currently closest
		// - > update closest components info
		{
			closestPair = citiesPair;
			closestComponentID = i;
		}
	}

	// Add best edge to the graph (with type `OTHER` and distance computed above).
	map.AddEdge(closestPair.FirstVertexID, closestPair.SecondVertexID, 'o', closestPair.Distance);
	this->addedEdges_.push_back(AddedEdge(closestPair.FirstVertexID, closestPair.SecondVertexID, 'o', closestPair.Distance));

	return true;
}


/// <summary>
/// For every vertex of degree 2 conected with edges of the same road type, 
/// merges its edges (to avoid current vertex) and finally deletes the vertex.
/// </summary>
/// <param name="map">Map object where to do the graph changes.</param>
bool GraphAdjuster::MergeDegreeTwoVertices(Map& map)
{
	bool wasChange = true;
	bool wasMerged = false;
	while (wasChange)
	// While there exist vertex of degree 2 conected with edges of the same road type.
	{
		wasChange = false;

		Graph& graph = map.GetGraph();
		auto allVertices = boost::vertices(graph);

		// Vector of the vertices to be deleted (sorted from the smallest to the greatest).
		std::vector<vertex_t> deletedVertices;

		for (auto it = allVertices.first; it != allVertices.second; it++)
		// Iterate through all the vertices.
		{
			if (this->MergeDegreeTwoVertex(*it, map))
			// Vertex merge was successfull.
			{
				wasChange = true;
				deletedVertices.push_back(*it);
			}

			//auto vertex_edges = graph.m_vertices[*it].m_out_edges;
			//if (vertex_edges.size() == 2)
			//// Vertex of the degree 2.
			//{
			//	vertex_t v1 = vertex_edges[0].m_target;
			//	vertex_t v2 = vertex_edges[1].m_target;
			//	edge_t e1 = boost::edge(*it, v1, graph).first;
			//	edge_t e2 = boost::edge(*it, v2, graph).first;
			//	Edge& edge1 = graph[e1];
			//	Edge& edge2 = graph[e2];

			//	if (e1 == e2)
			//		std::cout << "Rovnost hran" << std::endl;

			//	if (edge1.GetRoadType() == edge2.GetRoadType())
			//	// Both edges of the same type -> vertex to be deleted.
			//	{
			//		// Set change flag.
			//		wasChange = true;

			//		// Get new edge parameters.
			//		double newLength = edge1.GetLength() + edge2.GetLength();
			//		PermitedRoadTypes newRoadType = edge1.GetRoadType();

			//		// Remove the old edges.
			//		boost::remove_edge(e1, graph);
			//		boost::remove_edge(e2, graph);

			//		deletedVertices.push_back(*it);


			//		// Check whether edge to add already exists -> if so take the better one
			//		std::pair<edge_t, bool> edgeToAdd = boost::edge(v1, v2, graph);
			//		if (edgeToAdd.second)
			//		// Edge to add already exists.
			//		{
			//			Edge& duplicateEdge = graph[edgeToAdd.first];

			//			// New edge has lower road type -> keep the original edge.
			//			if (duplicateEdge.GetRoadType() > newRoadType)
			//				continue;

			//			// Edges has same road type, but original edge is shorter -> keep the original edge.
			//			if (duplicateEdge.GetRoadType() == newRoadType && duplicateEdge.GetLength() < newLength)
			//				continue;

			//			// Remove the original edge (end replace it with the new one).
			//			boost::remove_edge(edgeToAdd.first, graph);
			//		}

			//		// Add merged edge.
			//		map.AddEdge(v1, v2, newRoadType, newLength);
			//	}
			//}
		}


		int counter = 0;
		for (auto rit = deletedVertices.rbegin(); rit != deletedVertices.rend(); rit++)
		// Remove all found inconvenient vertices from the graph.
		// Going in reverse order because during deleting vertex renaming -> vertices with higher ID than deleted 
		// has ID diminished by 1 (if going in reverse, then only vertices that we don't want to delete).
		{
			if (counter % 10000 == 0)
				std::cout << "Iteration: " << counter << " Vertices of degree 2 merged and deleted." << std::endl;
			boost::remove_vertex(*rit, graph);
			counter++;
			wasMerged = true;
		}

		std::cout << "Iteration DONE" << std::endl;
		return wasMerged;
	}
}


/// <summary>
/// Tries to merge given vertex (if degree two and the same road type).
/// </summary>
/// <param name="vertex_to_merge">Vertex id to be merged.</param>
/// <param name="map">Map object where to do the graph changes.</param>
/// <returns>Returns `true` if vertex should be deleted (network was reconected), 
/// `false` do not delete the node.</returns>
bool GraphAdjuster::MergeDegreeTwoVertex(vertex_t vertex_to_merge, Map& map)
{
	Graph& graph = map.GetGraph();

	auto vertex_edges = graph.m_vertices[vertex_to_merge].m_out_edges;
	if (vertex_edges.size() == 2)
		// Vertex of the degree 2.
	{
		vertex_t v1 = vertex_edges[0].m_target;
		vertex_t v2 = vertex_edges[1].m_target;
		edge_t e1 = boost::edge(vertex_to_merge, v1, graph).first;
		edge_t e2 = boost::edge(vertex_to_merge, v2, graph).first;
		Edge& edge1 = graph[e1];
		Edge& edge2 = graph[e2];

		if (e1 == e2)
		{
			std::cout << "Rovnost hran" << std::endl;
			
		}

		if (edge1.GetRoadType() == edge2.GetRoadType())
			// Both edges of the same type -> vertex to be deleted.
		{
			// Get new edge parameters.
			double newLength = edge1.GetLength() + edge2.GetLength();
			PermitedRoadTypes newRoadType = edge1.GetRoadType();

			// Remove the old edges.
			boost::remove_edge(e1, graph);
			boost::remove_edge(e2, graph);

			// Check whether edge to add already exists -> if so take the better one
			std::pair<edge_t, bool> edgeToAdd = boost::edge(v1, v2, graph);
			if (edgeToAdd.second)
				// Edge to add already exists.
			{
				Edge& duplicateEdge = graph[edgeToAdd.first];

				// New edge has lower road type -> keep the original edge.
				if (duplicateEdge.GetRoadType() > newRoadType)
					// Delete the vertex (merge was successfull).
					return true;

				// Edges has same road type, but original edge is shorter -> keep the original edge.
				if (duplicateEdge.GetRoadType() == newRoadType && duplicateEdge.GetLength() < newLength)
					// Delete the vertex (merge was successfull).
					return true;

				// Remove the original edge (end replace it with the new one).
				boost::remove_edge(edgeToAdd.first, graph);
			}

			// Add merged edge.
			map.AddEdge(v1, v2, newRoadType, newLength);

			// Delete the vertex (merge was successfull).
			return true;
		}
	}

	// Do not delete the vertex.
	return false;
}


/// <summary>
/// Moves added edges container to the user. From that moment it is 
/// not permited to use any function from `GraphAdjuster` which uses
/// container `this->addedEdges_`. Typical usage is right before 
/// destroying `this` object to retrieve information about added edges.
/// </summary>
/// <returns>Return container of all added edges to the graph in order
/// to merge all graph components.</returns>
std::vector<AddedEdge> GraphAdjuster::GetAddedEdges()
{
	return std::move(this->addedEdges_);
}


/// <summary>
/// Moves container of the distances between each city pair to the user.
/// After calling this function it is not permited to use any function
/// from `GraphAdjuster` which uses container `this->citiesDistances_`.
/// Typical usage is right before destroying `this` object to retrieve
/// information about cities distances.
/// </summary>
/// <returns>Returns container of the distanes between each pair of cities.</returns>
std::vector< std::vector<double>> GraphAdjuster::GetCitiesDistances()
{
	return std::move(this->citiesDistances_);
}


/// <summary>
/// Computes distance between 2 cities based on longitude and latitude.
/// </summary>
/// <param name="firstCity">First city coordinates.</param>
/// <param name="secondCity">Second city coordinates.</param>
/// <returns>Returns computed distance between 2 cities in meters.</returns>
double GraphAdjuster::computeSphericalDistance(const City& firstCity, const City& secondCity)
{
	double dist;
	dist = std::sin(firstCity.Lat) * std::sin(secondCity.Lat)
		+ std::cos(firstCity.Lat) * std::cos(secondCity.Lat)
		* std::cos(firstCity.Lon - secondCity.Lon);

	dist = std::acos(dist);
	dist = (EARTH_RADIUS * PI * dist) / 180;

	return dist;
}

/// <summary>
/// Finds vector of cities-(distace_from_city, vertex) maps for each graph component.
/// </summary>
/// <param name="graph">Graph to get vertices from.</param>
/// <param name="numComponents">Number of components of the graph.</param>
/// <param name="verticesComponents">Vector of componentID for each vertex in graph.</param>
/// <returns>Returns vector of maps of cities-(distace_from_city, vertex) 
/// pairs for each component of the graph.</returns>
Graph_components GraphAdjuster::computeCityVertexPairs(const Graph& graph, int numComponents, std::vector<vertex_t>& verticesComponents)
{

	std::vector<std::vector<vertex_t>> components(numComponents);

	auto vert = boost::vertices(graph);

	// Gather vertices based on the components.
	for (auto it = vert.first; it != vert.second; ++it)
		components[verticesComponents[*it]].push_back(*it);


	// Find cities-(distace_from_city, vertex) pairs for each component, 
	// each city inside the component and its nearest vertex.
	Graph_components graphComponents(numComponents);
	for (size_t i = 0; i < numComponents; i++)
	{
		graphComponents[i] = this->findComponentCityVertexPairs(graph, components[i]);
	}

	return std::move(graphComponents);
}


/// <summary>
/// Finds city-(distance, vertex) pairs for each city and its closest vertex 
/// (and distance) for given component of the graph.
/// </summary>
/// <param name="graph">Graph to find vertices from</param>
/// <param name="component">Vector of vertices inside the given component.</param>
/// <returns>Returns computes Map of city and pair (distance_from_city, vertex) of 
/// closest vertices from all cities of the given component.</returns>
Component_city_vertex GraphAdjuster::findComponentCityVertexPairs(const Graph& graph, const std::vector<vertex_t>& component)
{
	Component_city_vertex citiesVertices;

	for (auto it = component.begin(); it != component.end(); it++)
		// For each component find all pairs of cities and its closest vetrtex (with its distance).
	{
		int cityID = graph[*it].CityID_;
		double cityDistance = graph[*it].CityDistance_;
		auto city_it = citiesVertices.find(cityID);

		if (city_it == citiesVertices.end())
			// City is not in the map of cities of the current component -> add city-actual_vertex pair
		{
			citiesVertices.insert(std::make_pair(cityID, std::make_pair(cityDistance, *it)));
			continue;
		}

		if ((*city_it).second.first > cityDistance)
			// If current vertex is closer to the city, than the current closest -> update best city-vertex pair
		{
			(*city_it).second.first = cityDistance;
			(*city_it).second.second = *it;
		}
	}

	return std::move(citiesVertices);
}



/// <summary>
/// Finds closest vertex pairs for the given pair of components.
/// </summary>
/// <param name="firstComponentID">ID of the first component to compute distance from.</param>
/// <param name="secondComponentID"></param>
/// <param name="graph_components"></param>
/// <returns></returns>
ComponentsDistance GraphAdjuster::findClosestCitiesPair(int firstComponentID, int secondComponentID, const Graph_components& graph_components)
{
	// Bad component IDs -> returns error Component pair.
	if (firstComponentID >= graph_components.size() || secondComponentID >= graph_components.size())
		return ComponentsDistance(-1, -1, -1);

	auto firstComponentMap = graph_components[firstComponentID];
	auto secondComponentMap = graph_components[secondComponentID];

	ComponentsDistance closestPair(0, -1, -1);
	ComponentsDistance currentPair;

	for (auto it1 = firstComponentMap.begin(); it1 != firstComponentMap.end(); it1++)
	{
		for (auto it2 = secondComponentMap.begin(); it2 != secondComponentMap.end(); it2++)
		{
			// Compute current city pair distance.
			currentPair.FirstVertexID = (*it1).second.second;
			currentPair.SecondVertexID = (*it2).second.second;
			currentPair.Distance =
				this->computeCitiesComponentDistance(
					std::make_pair((*it1).first, (*it2).first), (*it1).second.first, (*it2).second.first);

			if (closestPair.FirstVertexID == -1 || currentPair.Distance < closestPair.Distance);
			// If first pair or closer cities found -> update closest cities.
			closestPair = currentPair;
		}
	}

	return std::move(closestPair);
}


/// <summary>
/// Computes distance between the components (from cities distance and distance between its two closest vertices).
/// </summary>
/// <param name="citiesID">Pair of the cities to compute distance .</param>
/// <param name="firstVertexDistance">Distance of the closest vertex from the first city.</param>
/// <param name="secondVertexDistance">Distance of the closest vertex from the second city.</param>
/// <returns>Returns computed distance between two closest vertices from given cities 
/// from two different components of the graph.</returns>
double GraphAdjuster::computeCitiesComponentDistance(const std::pair<int, int>& citiesID, double firstVertexDistance, double secondVertexDistance)
{
	return this->citiesDistances_[citiesID.first][citiesID.second] + firstVertexDistance + secondVertexDistance;
}