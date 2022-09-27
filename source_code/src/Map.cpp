#include "optimizer/Map.hpp"


/// <summary>
/// Initializes map object.
/// </summary>
/// <param name="numClosestStations">Number of closest charging station to 
/// conside while planning the vehicle route.</param>
/// <param name="segmentLength">Length of the segment.</param>
Map::Map(int32_t numClosestStations, double segmentLength)
	: numClosestStations_(numClosestStations), segmentLength_(std::move(segmentLength))
{
	this->gen = std::mt19937(rd());

	this->lastClosestStation_ = 0;
}


/// <summary>
/// Initializes vector of locations for each vertex (for a-star heuristic function).
/// </summary>
void Map::InitHeuristicLocations()
{
	auto vert = boost::vertices(this->graph_);
	for (auto it = vert.first; it != vert.second; it++)
	{
		this->heuristicLocations_.push_back(this->graph_[*it].Location_);
	}
}


/// <summary>
/// Resets simulation information in the object (deletes charging stations and 
/// info about them and resets edges information about its battery levels).
/// </summary>
void Map::ResetSimulation()
{
	this->allChargingStations_.clear();
	this->closestStations_.clear();
	
	this->vertexClosestStations_.clear();

	
	auto allEdges = boost::edges(this->graph_);
	for (auto it = allEdges.first; it != allEdges.second; it++)
	{
		this->graph_[*it].ResetSimulation();
	}
}


/// <summary>
/// Adds node on the graph and sets its properties.
/// </summary>
/// <param name="newID">Current ID of the node.</param>
/// <param name="latitude">Latitude of the node in the real world.</param>
/// <param name="longitude">longitude of the node in the real world.</param>
/// <param name="cityID">ID of the nearest city of the node.</param>
/// <param name="distance">Distance between the node and its nearest city.</param>
/// <param name="oldID">Original ID of the node (from the source map), 
/// if not specified `-1` by default.</param>
void Map::AddNode(int32_t newID, double latitude, double longitude,
	int32_t cityID, double distance, int32_t oldID)
{
	vertex_t vertex = boost::add_vertex(this->graph_);
	this->graph_[vertex].NewID_ = newID;
	this->graph_[vertex].Location_.Latitude = latitude;
	this->graph_[vertex].Location_.Longitude = longitude;
	this->graph_[vertex].CityID_ = cityID;
	this->graph_[vertex].CityDistance_ = distance;
	this->graph_[vertex].OldID_ = oldID;
}


/// <summary>
/// Adds edge on the map and sets its properties.
/// </summary>
/// <param name="firstID">ID of the first node of the edge.</param>
/// <param name="secondID">ID of the second node of the edge.</param>
/// <param name="type">Type of the edge 
/// (`m` - motorway, `t` - trunk, `p` - primary, `o` - other).</param>
/// <param name="length">Length of the edge (in kilometers).</param>
void Map::AddEdge(int32_t firstID, int32_t secondID, char type, double length)
{
	edge_t edge; bool b;
	auto num_vertices = this->graph_.m_vertices.size();

	// Check whether the vertex with given id exists (if not, ignore the edge).
	if (num_vertices <= firstID || num_vertices <= secondID)
		return;

	// Set edge properties.
	boost::tie(edge, b) = boost::add_edge(firstID, secondID, this->graph_);

	this->graph_[edge].AddType(std::move(type));
	this->graph_[edge].SetLength(std::move(length), this->segmentLength_);
}


/// <summary>
/// Adds edge on the map and sets its properties.
/// </summary>
/// <param name="firstID">ID of the first node of the edge.</param>
/// <param name="secondID">ID of the second node of the edge.</param>
/// <param name="type">Road type of the edge.</param>
/// <param name="length">Length of the edge (in kilometers).</param>
void Map::AddEdge(vertex_t firstID, vertex_t secondID, PermitedRoadTypes type, double length)
{
	edge_t edge; bool b;
	auto num_vertices = this->graph_.m_vertices.size();

	// Check whether the vertex with given id exists (if not, ignore the edge).
	if (num_vertices <= firstID || num_vertices <= secondID)
		return;

	// Set edge properties.
	boost::tie(edge, b) = boost::add_edge(firstID, secondID, this->graph_);

	this->graph_[edge].AddType(std::move(type));
	this->graph_[edge].SetLength(std::move(length), this->segmentLength_);
}


/// <summary>
/// Adds population of the given city into the desired storage. 
/// If city already initialized, then rewrites the info. If cities with lower 
/// cityID not defined, defines it with zero population.
/// </summary>
/// <param name="cityID">ID of the city to add.</param>
/// <param name="population">Population of the city.</param>
void Map::AddCity(int32_t cityID, int32_t population)
{
	auto numCities = this->cityPopulations_.size();
	if (numCities > cityID)
	// If city already in the set (just actualize the city population).
	{
		this->cityPopulations_[cityID] = population;
		return;
	}

	if (numCities < cityID)
	// If some cities missing (fill vector with them (with population 0)).
	// To ensure city is stored in vector on index (`cityID`).
	{
		for (size_t i = 0; i < cityID - numCities - 1; i++)
			// Add missing cities with ID less than `cityID`.
		{
			this->cityPopulations_.push_back(0);
		}
	}

	// Add given city.
	this->cityPopulations_.push_back(std::move(population));
}


/// <summary>
/// Adds given charging station object to the container of all charging stations of the map.
/// </summary>
/// <param name="chargingStation">Object representing charging station to be added.</param>
void Map::AddChargingStation(ChargingStation chargingStation)
{
	this->allChargingStations_.push_back(std::make_unique<ChargingStation>(chargingStation));
}


/// <summary>
/// Stores cities distances information into proper variable.
/// </summary>
/// <param name="citiesDistances">Matrix of distances between each 
/// pair of the city (to be stored in the proper variable).</param>
void Map::SetCitiesDistances(Double2DMatrix citiesDistances)
{
	this->citiesDistances_ = std::move(citiesDistances);
}


/// <summary>
/// </summary>
/// <returns>Returns constant reference to graph.</returns>
const Graph& Map::GetConstGraph()
{
	return this->graph_;
}


/// <summary>
/// </summary>
/// <returns>Returns reference to graph of the map.</returns>
Graph& Map::GetGraph()
{
	return this->graph_;
}


/// <summary>
/// </summary>
/// <param name="nodeID">ID of the node to get reference.</param>
/// <returns>Returns const reference to the node with given ID.</returns>
const Node& Map::GetNode(int32_t nodeID)
{
	return this->graph_[nodeID];
}


/// <summary>
/// </summary>
/// <returns>Returns reference to the container of all charging stations.</returns>
std::vector<std::unique_ptr<ChargingStation>>& Map::GetChargingStations()
{
	return this->allChargingStations_;
}


/// <summary>
/// </summary>
/// <param name="stationID">ID of the wanted charging station.</param>
/// <returns>Returns reference to the charging station with given ID.</returns>
std::unique_ptr<ChargingStation>& Map::GetChargingStation(int32_t stationID)
{
	return this->allChargingStations_[stationID];
}


/// <summary>
/// </summary>
/// <param name="firstCityID">ID of the first city.</param>
/// <param name="secondCityID">ID of the second city.</param>
/// <returns>Returns distance between two cities.</returns>
double Map::GetCityPairDistance(int32_t firstCityID, int32_t secondCityID)
{
	return this->citiesDistances_[firstCityID][secondCityID];
}


/// <summary>
/// </summary>
/// <returns>Returns const reference to each city population container.</returns>
const std::vector<int32_t>& Map::GetCitiesPopulation()
{
	return this->cityPopulations_;
}


/// <summary>
/// </summary>
/// <param name="cityID">City ID to get vertices from.</param>
/// <returns>Returns constant reference to the container of all city 
/// nodes for the given city.</returns>
const std::vector<vertex_t>& Map::GetCityNodes(int32_t cityID)
{
	return this->cityNodes_[cityID];
}


/// <summary>
/// </summary>
/// <returns>Returns reference to charging station which was used 
/// in the last route planning through charging station (typicaly used to get
/// information about position of the station for vehicle object).</returns>
std::unique_ptr<ChargingStation>& Map::GetLastChargingStation()
{
	return this->allChargingStations_[this->lastClosestStation_];
}


/// <summary>
/// </summary>
/// <returns>Returns map of edge identifier and its corresponding battery levels data.</returns>
std::map<edge_t, std::vector<BateryPair>&> Map::GetAllSegments()
{
	// Create the wanted map.
	std::map<edge_t, std::vector<BateryPair>&> allSegments;
	auto edges = boost::edges(this->graph_);
	for (auto it = edges.first; it != edges.second; it++)
	{
		allSegments.insert(std::pair<edge_t, std::vector<BateryPair>&>
			(*it, this->graph_[*it].GetBatteryLevels()));
	}

	return allSegments;
}


/// <summary>
/// Checks whether node with given ID exists.
/// </summary>
/// <param name="nodeID">Node ID to check.</param>
/// <returns>Returns `true` if node exists, else `false`.</returns>
bool Map::ExistNode(int32_t nodeID)
{
	// Node with given ID exists (indexing from 0 sequentaly).
	if (this->graph_.m_vertices.size() > nodeID)
		return true;

	// Node does not exist.
	return false;
}

/// <summary>
/// Computes distance between 2 cities based on longitude and latitude 
/// (using the 'Haversine' formula).
/// </summary>
/// <param name="firstCity">First city coordinates.</param>
/// <param name="secondCity">Second city coordinates.</param>
/// <returns>Returns distance between 2 cities in kilometers.</returns>
double Map::ComputeSphericalDistance(double firstLatitude, double firstLongitude,
	double secondLatitude, double secondLongitude)
{
	double dist;

	dist = std::sin(firstLatitude) * std::sin(secondLatitude)
		+ std::cos(firstLatitude) * std::cos(secondLatitude)
		* std::cos(firstLongitude - secondLongitude);

	dist = std::acos(dist);
	dist = (EARTH_RADIUS * PI * dist) / 180;

	return dist;
}


/// <summary>
/// Finds all nodes which belongs to the city and stores them into proper container.
/// </summary>
void Map::FindCityNodes()
{
	this->cityNodes_ = std::vector < std::vector<vertex_t>>(
		this->cityPopulations_.size());
	auto nodes = boost::vertices(this->graph_);
	for (auto it = nodes.first; it < nodes.second; it++)
	{
		this->cityNodes_[this->graph_[*it].CityID_].push_back(*it);
	}
}


/// <summary>
/// Finds `this->numClosestStations_` nearest stations for each city in the map and
/// stores them in ascending order (by distance from city) in the proper container.
/// For path searching through charging stations.
/// </summary>
void Map::FindCitiesNearestChargingStations()
{
	// --------------------------------UNCOMMENT TO GET PREVIOUS VARIANT:

	// for (size_t i = 0; i < this->cityPopulations_.size(); i++)
	// {
	// 	this->closestStations_.push_back(this->findNearestChargingStations(i));
	// }
	
	// --------------------------------END OF UNCOMMENT TO GET PREVIOUS VARIANT:

	// --------------------------------COMMENT TO GET PREVIOUS VARIANT:

	this->FindVerticesNearestChargingStations();

	// --------------------------------END OF COMMENT TO GET PREVIOUS VARIANT:
}


void Map::FindVerticesNearestChargingStations()
{
	for (size_t i =0; i < this->graph_.m_vertices.size(); i++)
	{
		this->vertexClosestStations_.push_back(this->findNearestStationForVertex(i));
	}
}

/// <summary>
/// Finds the closest vertex from the given coordinates, then its neighbor closest
/// to the given coordinates and on the edge connecting these two vertices randomly generates 
/// segment closer to the closest vertex. It (partly randomly) 
/// approximates the position of the given coordinates 
/// (for K-Means algorithm and finding centroids there).
/// </summary>
/// <param name="latitude">Latitude of the searched position.</param>
/// <param name="longitude">Longitude of the searched position.</param>
/// <returns>Returns the approximately nearest position on the map from the given position.</returns>
MapPosition Map::FindClosestPositionToCoordinates(double latitude, double longitude)
{
	// Find the closest vertex from the given coordinates.
	std::pair<vertex_t, double> closestVertexValues = 
		this->findClosestVertex(latitude, longitude, -1, true);

	// Find the closest neighbor of the closest vertex from the given coordinates.
	std::pair<vertex_t, double> closestNeighborValues = 
		this->findClosestVertex(latitude, longitude, closestVertexValues.first, false);

	// Get edge where should be the closest point from the given coordinates.
	Edge& edge = this->graph_[
		boost::edge(closestVertexValues.first, closestNeighborValues.first, this->graph_).first];

	int32_t numSegments = edge.GetNumSegments();

	// Uniformly randomly generate segment on the edge candidate 
	// (the furthest segment from the vertex can be in the middle of the edge).
	std::uniform_int_distribution<> segmentDistribution(0, numSegments / 2);
	int32_t segmentID = segmentDistribution(this->gen);

	// Closer vertex has higher ID -> segmentID should be in the second half of the segment IDs 
	// (closer vertex should be closer to the given segment).
	if (closestNeighborValues.first < closestVertexValues.first)
		segmentID += numSegments / 2;

	// Segment ID overflow fix -> set segment ID to highest possible.
	if (segmentID >= numSegments)
		segmentID = numSegments - 1;

	return MapPosition(closestVertexValues.first, closestNeighborValues.first, segmentID, this->graph_);
}



/// <summary>
/// Finds estimated shortest path between two given vertices. If we want to go charging, 
/// then finds estimated shorted path through some charging station (returns only path
/// to the charging station, rest of the path should be computed after the vehicle is charged).
/// </summary>
/// <param name="start">Start vertex of the path.</param>
/// <param name="end">Final vertex of the path.</param>
/// <param name="goCharging">Flag if car should go charging then `true`, else `false`.</param>
/// <returns>Returns the pair of estimated relative path duration and the 
/// container of the reversed path.</returns>
LengthPathPair Map::FindVehiclePath(vertex_t start, vertex_t end, bool goCharging)
{
	if (!goCharging)
	// Don't go charging -> just find the shortest path to the end.
	{
		return this->FindShortestPath(start, end);

	}
	else
	// Go charging, rearrange the shortest path to go through charging station.
	// Choose the approximately best variant from the closest stations for current city.
	{
		// --------------------------UNCOMENT TO GET PREVIOUS VARIANT:

		// LengthPathPair bestPath, currentPath;
		// double bestLength, helpLength;

		// for (size_t i = 0; i < this->numClosestStations_; i++)
		// // Check only `this->numClosestStations_` closest stations of the current city.
		// {
		// 	std::unique_ptr<ChargingStation>& currentStation =
		// 		this->allChargingStations_[this->closestStations_[
		// 			this->graph_[start].CityID_][i]];


		// 	currentPath = this->FindShortestPath(start,
		// 		(*currentStation).Position_.CloserVertexID);

		// 	helpLength = currentPath.first;

		// 	// Compute length of the path from the station to the end 
		// 	// (reckon that car will pass the edge (not returning on the same segments 
		// 	// (approximation of the path length)).
		// 	helpLength +=
		// 		this->FindShortestPath((*currentStation).Position_.CloserVertexID, end).first;
		// 	helpLength +=
		// 		this->graph_[(*currentStation).Position_.EdgeID].GetLength();

		// 	if (i == 0 || helpLength < bestLength)
		// 	// If first tested path or current length is smaller than 
		// 	// currently best -> update best path values
		// 	{
		// 		bestLength = helpLength;
		// 		bestPath = currentPath;

		// 		this->lastClosestStation_ = (*currentStation).StationID_;
		// 	}
		// }

		// --------------------------------END OF UNCOMMENT TO GET PREVIOUS VARIANT:

		// --------------------------------COMMENT TO GET PREVIOUS VARIANT:

		LengthPathPair bestPath = this->FindShortestPath(start, 
			this->allChargingStations_[this->vertexClosestStations_[start]]->Position_.CloserVertexID);
		this->lastClosestStation_ = this->vertexClosestStations_[start];
		
		// -------------------------------END OF COMMENT TO GET PREVIOUS VARIANT:


		// Returns best path only to the station 
		// (path from the charging station will be computed later).
		return bestPath;
	}
}


/// <summary>
/// Finds shortest path between two choosen vertices using A-Star algorithm.
/// </summary>
/// <param name="start">Start vertex.</param>
/// <param name="goal">Final vertex.</param>
/// <returns>Returns pair of estimated path length and reversed path of the vertices
///  from `start` to `goal`.</returns>
LengthPathPair Map::FindShortestPath(const vertex_t& start, const vertex_t& goal)
{
	std::vector<vertex_t> predecessors(boost::num_vertices(this->graph_));

	std::vector<vertex_t> path;
	double pathLength = 0;

	try 
	// Use A-Star algorithm to compute distances (until it finds the `goal` vertex).
	{
		boost::astar_search
		(this->graph_, start,
			distance_heuristic<Graph, float, std::vector<location>> 
			(this->heuristicLocations_, goal),
			boost::weight_map(boost::get(&Edge::TransitTime, this->graph_))
			.predecessor_map(boost::make_iterator_property_map(predecessors.begin(),
				boost::get(boost::vertex_index, this->graph_)))
			.visitor(astar_goal_visitor<vertex_t>(goal)));
	}
	catch (found_goal fg) 
	// Found a path to the goal -> stop searching other paths.
	{
		for (vertex_t v = goal;; v = predecessors[v])
		{
			if (predecessors[v] == v && v != start)
			// Test path don't exist.
			{
				std::cout << "Bad path in edge " << v << std::endl;
				return std::make_pair(-1, path);
			}

			if (v != start)
			// Still not in the starting vertex -> update path.
			{
				pathLength += this->graph_
					[boost::edge(v, predecessors[v], this->graph_).first].TransitTime;
			}

			path.push_back(v);

			// Whole path found.
			if (predecessors[v] == v)
				break;
		}
	}

	return std::make_pair(pathLength, std::move(path));
}


/// <summary>
/// Performs Dijkstra algorithm to find distance of each vertex from the start vertex.
/// </summary>
/// <param name="start">Start vertex to compute distance from.</param>
/// <returns>Returns container of distances of each vertex from the `start` vertex.</returns>
std::vector<double> Map::DijkstraFindDistances(vertex_t start)
{
	std::vector<double> distances(boost::num_vertices(this->graph_));

	// Perform dijkstra algorithm and store distances.
	boost::dijkstra_shortest_paths(this->graph_, start,
		boost::weight_map(boost::get(&Edge::TransitTime, this->graph_))
		.distance_map(boost::make_iterator_property_map(distances.begin(),
			boost::get(boost::vertex_index, this->graph_))));

	return std::move(distances);
}


/// <summary>
/// Finds components of the graph.
/// </summary>
/// <param name="numComponentsTuple">Tuple where to store number of components (first parameter)
/// and vector of component IDs for each vertex.</param>
void Map::TestComponents(std::tuple<int32_t, std::vector<vertex_t>>& numComponentsTuple)
{
	std::get<1>(numComponentsTuple) = std::vector<vertex_t>(boost::num_vertices(this->graph_));
	std::get<0>(numComponentsTuple)
		= boost::connected_components(this->graph_, &std::get<1>(numComponentsTuple)[0]);

	std::cout << "Total number of components: " << std::get<0>(numComponentsTuple) << std::endl;
}


/// <summary>
/// Simulates vehicle passage through the given edge (updates battery usage info).
/// </summary>
/// <param name="edge">Identifier of the edge to be updated.</param>
/// <param name="batteryPair">Pair of in and out battery level.</param>
/// <param name="startSegment">Starting segment of the vehicle.</param>
/// <param name="segmentIncrease">Flag whether vehicle increasing in segments
/// (`true` if so, else `false`).</param>
/// <param name="endSegment">End Segment of the path 
/// (if `-1` then automaticaly finds appropriate end of the edge).</param>
void Map::SimulateVehiclePassedThroughEdge(edge_t edge,
	std::pair<double, double> batteryPair, int32_t startSegment,
	bool segmentIncrease, int32_t endSegment)
{
	this->graph_[edge].UpdateSegmentData(batteryPair, startSegment, segmentIncrease, endSegment);
}


/// <summary>
/// For given city finds its `this->numClosestStations_` nearest station IDs.
/// </summary>
/// <param name="CityID">ID of the city to find nearest stations from.</param>
/// <returns>Returns array of the `this->numClosestStations_` 
/// nearest stations IDs from the given city.</returns>
std::vector<int32_t> Map::findNearestChargingStations(int32_t cityID)
{
	std::vector<int32_t> nearestStations(this->numClosestStations_, -1);

	// Help array to remember distance of each station to the city center.
	std::vector<double> stationsDistances(this->numClosestStations_);
	for (auto&& station : this->allChargingStations_)
	// Check all stations to find the nearest ones.
	{
		Node stationVertex = this->graph_[(*station).Position_.CloserVertexID];
		double stationDistance = stationVertex.CityDistance_;

		if (stationVertex.CityID_ != cityID)
		// If station in the other city -> add cities distances to total distance.
		{
			stationDistance += this->citiesDistances_[stationVertex.CityID_][cityID];
		}

		for (size_t i = 0; i < this->numClosestStations_; i++)
		// Check the position in the currently best charging stations (and update info).
		{
			if (nearestStations[i] == -1)
			// Not enough charging stations found -> put current station into the container.
			{
				nearestStations[i] = (*station).StationID_;
				stationsDistances[i] = stationDistance;
				break;
			}

			if (stationsDistances[i] > stationDistance)
			// Current station is closer than the currently i-th nearest 
			// -> move all further stations one position behind.
			{
				// Define help variables to rearrange array.
				int helpID1, helpID2;
				helpID1 = nearestStations[i];
				double helpDistance1, helpDistance2;
				helpDistance1 = stationsDistances[i];

				// Set current station to i-th nearest position.
				nearestStations[i] = (*station).StationID_;
				stationsDistances[i] = stationDistance;

				for (size_t j = i + 1; j < this->numClosestStations_; j++)
				// Move all further stations 1 position higher (more distant).
				{
					// Save currently j-th station info.
					helpID2 = nearestStations[j];
					helpDistance2 = stationsDistances[j];

					// Move (j-1)-th station to j-th position. 
					nearestStations[j] = helpID1;
					stationsDistances[j] = helpDistance1;

					// Rearrange help variables.
					helpID1 = helpID2;
					helpDistance1 = helpDistance2;
				}

				// Don't check rest of the further stations.
				break;
			}
		}
	}

	return std::move(nearestStations);
}



int32_t Map::findNearestStationForVertex(vertex_t vertexID)
{
	auto distances = this->DijkstraFindDistances(vertexID);
	int32_t bestStationID = -1;
	double bestDistance = 0;
	for (auto&& station : this->allChargingStations_)
	{
		if (bestStationID == -1 || bestDistance > distances[station->Position_.CloserVertexID])
		{
			bestStationID = station->StationID_;
			bestDistance = distances[station->Position_.CloserVertexID];
		}
	}

	return bestStationID;
}




/// <summary>
/// Finds closest vertex from the given coordinates and from the specified neighborhood.
/// </summary>
/// <param name="latitude">Latitude of the point to find the closest vertex from.</param>
/// <param name="longitude">Longitude of the point to find the closest vertex from.</param>
/// <param name="from">Vertex which neighbors to check 
/// (when checking only neighbors of the vertex).</param>
/// <param name="wholeMap">Flag signifying that all vertices from the map should be 
/// considered as a closest vertex candidates (`false` - check only neighbors of the vertex `from`).</param>
/// <returns>Returns pair of the closest vertex and its distance from the given point.</returns>
std::pair<vertex_t, double> Map::findClosestVertex(
	double latitude, double longitude, vertex_t from, bool wholeMap)
{
	double smallestDistance = -1;
	vertex_t closestVertex;

	
	if (!wholeMap)
	// Check only neighbors of the vertex `from`.
	{
		auto iteratorRange = 
			boost::make_iterator_range(boost::adjacent_vertices(from, this->graph_));

		for (auto it : iteratorRange)
		// Go through all vertices and find the closest one.
		{
			auto currentCoordinates = this->graph_[it].Location_;
			double currentDistance = this->ComputeSphericalDistance(latitude, longitude,
				currentCoordinates.Latitude, currentCoordinates.Longitude);

			if (smallestDistance < 0 || currentDistance < smallestDistance)
			// First vertex to check or closer than currently closest vertex 
			// -> update the best.
			{
				closestVertex = it;
				smallestDistance = currentDistance;
			}
		}
	}
	else
	// Check the whole map.
	{
		auto vert = boost::vertices(this->graph_);
		auto iteratorRange = boost::make_iterator_range(vert);

		for (auto it : iteratorRange)
		// Go through all vertices and find the closest one.
		{
			auto currentCoordinates = this->graph_[it].Location_;
			double currentDistance = this->ComputeSphericalDistance(latitude, longitude,
				currentCoordinates.Latitude, currentCoordinates.Longitude);

			if (smallestDistance < 0 || currentDistance < smallestDistance)
			// First vertex to check or closer than currently closest vertex 
			// -> update best.
			{
				closestVertex = it;
				smallestDistance = currentDistance;
			}
		}
	}

	return std::make_pair(closestVertex, smallestDistance);
}