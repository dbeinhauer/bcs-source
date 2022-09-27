#include "optimizer/MapReader.hpp"


/// <summary>
/// Initializes the reader.
/// </summary>
MapReader::MapReader() 
{
	this->lastFirstNode_ = 0;
}


/// <summary>
/// Loads map from the file streams from given parameters in the appropriate format.
/// </summary>
/// <param name="map">Object of the map to store the loaded information.</param>
/// <param name="simulationParameters">Object which clusters all simulation parameters.</param>
/// <returns>Return `true` if loading was successful, else `false` (error happened).</returns>
bool MapReader::LoadMap(Map& map, SimulationParameters& simulationParameters)
{
	// Initialize the file streams to read map represantation from.
	std::ifstream nodeCityStream(simulationParameters.NodeCityFile);
	std::ifstream edgeStream(simulationParameters.EdgesFile);
	std::ifstream cityStream(simulationParameters.CitiesFile);
	std::ofstream addedEdgesStream(simulationParameters.AddedEdgesFile);


	// Load nodes.
	if (!this->LoadNodesCities(map, nodeCityStream))
	// Error happened -> end function
	{
		nodeCityStream.close();
		return false;
	}
	nodeCityStream.close();


	// Load edges.
	if (!(this->LoadEdges(map, edgeStream)))
	// Error happened -> end function
	{
		edgeStream.close();
		return false;
	}
	edgeStream.close();



	if (this->graphAdjuster_.MergeDegreeTwoVertices(map) && simulationParameters.SavePrepared)
	// Some vertices of degree 2 was merged -> save the preprocessed map (because we want to).
	{
		// Prepare output file streams.
		std::ofstream preparedNodesStream(simulationParameters.PreparedNodesFile);
		std::ofstream preparedEdgesStream(simulationParameters.PreparedEdgesFile);

		// Save prepared nodes.
		if (!this->WritePreparedCombinedNodes(map, preparedNodesStream))
			// Error happened -> end function
		{
			preparedNodesStream.close();
			return false;
		}
		preparedNodesStream.close();

		// Save prepared edges.
		if (!this->WritePreparedEdges(map, preparedEdgesStream))
		// Error happened -> end function
		{
			preparedEdgesStream.close();
			return false;
		}
		preparedEdgesStream.close();
	}


	// Load cities.
	if (!this->LoadCities(map, cityStream))
	// Error happened -> end function
	{
		cityStream.close();
		return false;
	}
	cityStream.close();


	// Make the graph connected.
	std::tuple<int32_t, std::vector<vertex_t>> components;
	map.TestComponents(components);

	if (std::get<0>(components) != 1)
	// If there are multiple components of the graph -> merge them
	{
		this->graphAdjuster_.ComputeCitiesDistances();
		while (true)
		// Step by step merge all components to final connected graph.
		{
			std::tuple<int32_t, std::vector<vertex_t>> components;

			// Find components and get its properties.
			map.TestComponents(components);

			int32_t numComponents = std::get<0>(components);
			std::vector<vertex_t> vertexComponentPairs = std::get<1>(components);


			// Error happened -> end function.
			if (numComponents < 1)
				return false;

			// Graph is connected -> end merging
			if (numComponents == 1)
				break;

			// Merge some components.
			this->graphAdjuster_.MergeTwoComponents(map, numComponents, std::move(vertexComponentPairs));
		}


		if (simulationParameters.SavePrepared)
		// If we want to write added edges to make the graph connected to the proper file.
		{
			if (!this->WriteAddedEdges(addedEdgesStream))
			// Error happened -> end function.
			{
				addedEdgesStream.close();
				return false;
			}
			addedEdgesStream.close();
		}
	}

	// Find all nodes which belongs to the corresponding city.
	map.FindCityNodes();
	map.InitHeuristicLocations();

	return true;
}


/// <summary>
/// Loads node-city pairs from the input.
/// </summary>
/// <param name="map">Object to store info about node-city pair.</param>
/// <param name="nodeCityStream">Input stream.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool MapReader::LoadNodesCities(Map& map, std::istream& nodeCityStream)
{
	std::string line;

	// Skip header -> if not possible -> Error
	if (!std::getline(nodeCityStream, line))
		return false;

	// Read whole file line by line.
	while (std::getline(nodeCityStream, line))
	{
		// Error during the loading.
		if (!(this->loadNodeCityPair(map, line)))
			return false;
	}

	// Loading was successful.
	return true;
}


/// <summary>
/// Loads edges from the input.
/// </summary>
/// <param name="map">Object to store info about the edges.</param>
/// <param name="edgeStream">Input stream.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool MapReader::LoadEdges(Map& map, std::istream& edgeStream)
{
	std::string line;

	// Skip header -> if not possible -> Error
	if (!std::getline(edgeStream, line))
		return false;

	// Read whole file line by line.
	while (std::getline(edgeStream, line))
	{
		// Error during the loading.
		if (!this->loadEdge(map, line))
			return false;
	}

	// Loading was successful.
	return true;
}


/// <summary>
/// Loads cities from the input.
/// </summary>
/// <param name="map">Object to store info about the cities.</param>
/// <param name="cityStream">Input stream.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool MapReader::LoadCities(Map& map, std::istream& cityStream)
{
	std::string line;

	// Skip header -> if not possible -> Error
	if (!std::getline(cityStream, line))
		return false;

	// Read whole file line by line.
	while (std::getline(cityStream, line))
	{
		// Error during the loading.
		if (!this->loadCity(map, line))
			return false;
	}

	// Loading was successful.
	return true;
}


/// <summary>
/// Loads charging stations.
/// </summary>
/// <param name="map">Object to store info about the stations.</param>
/// <param name="chargingStationsStream">Input stream.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool MapReader::LoadChargingStations(Map& map, std::istream& chargingStationsStream)
{
	std::string line;

	// Skip header (2 lines) -> if not possible -> Error
	if (!std::getline(chargingStationsStream, line))
		return false;
	if (!std::getline(chargingStationsStream, line))
		return false;


	// Read whole file line by line.
	while (std::getline(chargingStationsStream, line))
	{
		// Error during the loading.
		if (!this->loadChargingStation(map, line))
			return false;
	}

	// Loading was successful.
	return true;
}


/// <summary>
/// Writes added edges used for components merging.
/// In format:
///		`first_vertex second_vertex length type`
/// </summary>
/// <param name="addedEdgesStream">Stream to write info about added adges.</param>
/// <returns>Returns `true` if writing was successful, else `false`.</returns>
bool MapReader::WriteAddedEdges(std::ostream& addedEdgesStream)
{
	auto addedEdges = this->graphAdjuster_.GetAddedEdges();

	for (auto&& edge : addedEdges)
	{
		addedEdgesStream << edge.FirstVertex << " " << edge.SecondVertex << " "
			<< edge.Length << " " << edge.Type << std::endl;
	}

	return true;
}


/// <summary>
/// Writes preprocessed combined nodes (after removing degree 2 vertices from the original graph).
/// In format:
///		`node_ID latitude longitude nearest_city_ID city_distance`
/// </summary>
/// <param name="map">Map object to get graph info.</param>
/// <param name="preparedNodesStream">Stream to write node-city info.</param>
/// <returns>Returns `true` if writing was successful, else `false`.</returns>
bool MapReader::WritePreparedCombinedNodes(Map& map, std::ostream& preparedNodesStream)
{
	// Write header.
	preparedNodesStream << "node_id latitude longitude city_id distance_km" << std::endl;
	Graph& graph = map.GetGraph();

	// Write info about each vertex.
	for (size_t i = 0; i < boost::num_vertices(graph); i++)
	{
		preparedNodesStream << i << " " //<< graph[i].OldID_ << " " 
			<< graph[i].Location_.Latitude << " " << graph[i].Location_.Longitude << " "
			<< graph[i].CityID_ << " " << graph[i].CityDistance_ << std::endl;
	}

	return true;
}


/// <summary>
/// Writes preprocessed edges (after removing degree 2 vertices from the original graph).
/// In format:
///		`node_id_1 node_id_2 length road_type` 
/// </summary>
/// <param name="map">Map object to get graph info.</param>
/// <param name="preparedEdgesStream">Stream to write edges info.</param>
/// <returns>Returns `true` if writing was successful, else `false`.</returns>
bool MapReader::WritePreparedEdges(Map& map, std::ostream& preparedEdgesStream) 
{
	// Write header.
	preparedEdgesStream << "new_node_id_1 new_node_id_2 length road_type" << std::endl;

	Graph& graph = map.GetGraph();
	auto edges_iter = boost::edges(graph);

	// Write info about each edge.
	for (auto it = edges_iter.first; it != edges_iter.second; it++)
	{
		preparedEdgesStream << (*it).m_source << " " << (*it).m_target << " " << 
			graph[*it].GetLength() << " " << graph[*it].GetRoadTypeChar() << std::endl;
	}

	return true;
}


/// <summary>
///	Writes all charging stations representations into the given stream.
/// In format:
///		`station_ID closer_vertex further_vertex segment capacity city_ID waiting_time estimated_charging`
/// </summary>
/// <param name="map">Map object to load charging stations from.</param>
/// <param name="stationsStream">Stream to write charging stations configurations.</param>
/// <returns>Returns `true` if writing was successfull, else `false`.</returns>
bool MapReader::WriteStations(ModelRepresentation& model, Map& map, std::ostream& stationsStream)
{
	// auto& chargingStations = map.GetChargingStations();
	auto& chargingStations = model.AllChargingStations_;

	// Write header.
	// stationsStream << "Total number of charging stations in the map: " << chargingStations.size() << std::endl;

	stationsStream << "Station_ID Latitude Longitude Closer_Vertex_ID Further_Vertex_ID Segment_ID NumTotalSegments Capacity" << std::endl;

	// Write info about each station.
	for (auto&& station : chargingStations)
	{
		auto edge_t = boost::edge(station.Position_.CloserVertexID, station.Position_.FurtherVertexID, map.GetGraph());
		Edge& edge = map.GetGraph()[edge_t.first];

		double coordinatesOffset = (double)(station).Position_.EdgeSegmentID / (double)edge.GetNumSegments();

		vertex_t lowerIDVertex = station.Position_.CloserVertexID;
		vertex_t higherIDVertex = station.Position_.FurtherVertexID;
		if (lowerIDVertex > higherIDVertex)
			std::swap(lowerIDVertex, higherIDVertex);

		double latitudeDifference =
			map.GetGraph()[higherIDVertex].Location_.Latitude - 
			map.GetGraph()[lowerIDVertex].Location_.Latitude;

		double longitudeDifference = 
			map.GetGraph()[higherIDVertex].Location_.Longitude - 
			map.GetGraph()[lowerIDVertex].Location_.Longitude;
		
		double latitude = map.GetGraph()[lowerIDVertex].Location_.Latitude + latitudeDifference * coordinatesOffset;
		double longitude = map.GetGraph()[higherIDVertex].Location_.Longitude + longitudeDifference * coordinatesOffset;
		
		stationsStream << station.StationID_ << " " << latitude << " " << longitude << " " << station.Position_.CloserVertexID << " "
			<< station.Position_.FurtherVertexID << " " << (station).Position_.EdgeSegmentID << " "
			<< edge.GetNumSegments() << " " << (station).Capacity_ << std::endl;
	}

	return true;
}


/// <summary>
/// </summary>
/// <returns>Returns reference to currently used `GraphAdjuster` object.</returns>
GraphAdjuster& MapReader::GetGraphAdjuster()
{
	return this->graphAdjuster_;
}



/// <summary>
/// Loads one pair of the node-city pair and stores info in the map.
/// The input is in format:
///			`nodeID latitude longitude cityID distance`
/// </summary>
/// <param name="map">Object to store info about the pair.</param>
/// <param name="line">String with pair representation.</param>
/// <returns>Returns `true` if loading is successful, else `false`.</returns>
bool MapReader::loadNodeCityPair(Map& map, const std::string& line)
{
	std::stringstream stream(line);

	int32_t nodeID, cityID;
	double latitude, longitude, distance;

	// Skip empty line.
	if (std::all_of(line.begin(), line.end(), isspace))
		return true;

	// Bad input format -> Error
	if (!(stream >> nodeID >> latitude >> longitude >> cityID >> distance))
		return false;

	// Node to add already exists -> bad format of the input -> Error
	if (map.ExistNode(nodeID))
		return false;

	// Some of the previous nodes are missing in the input -> bad format -> Error
	if (nodeID > 0 && !map.ExistNode(nodeID - 1))
		return false;

	// Correct format -> add node
	map.AddNode(std::move(nodeID), std::move(latitude), std::move(longitude),
		std::move(cityID), std::move(distance));

	return true;
}


/// <summary>
/// Loads one edge from the input representation in format:
///		`firstNodeID secondNodeID length roadType`
/// </summary>
/// <param name="map">Object to store info about the node.</param>
/// <param name="line">Input string with edge representation.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool MapReader::loadEdge(Map& map, const std::string& line)
{
	std::stringstream stream(line);

	int32_t firstNodeID, secondNodeID;
	double length;
	char roadType;

	// Skip empty line.
	if (std::all_of(line.begin(), line.end(), isspace))
		return true;

	// Bad input format -> Error
	if (!(stream >> firstNodeID >> secondNodeID >> length >> roadType))
		return false;

	// Given nodes are not in the graph representation -> Error
	if (!map.ExistNode(firstNodeID) || !map.ExistNode(secondNodeID))
		return false;

	if (boost::edge(firstNodeID, secondNodeID, map.GetGraph()).second)
		return true;


	// Correct format -> add edge
	map.AddEdge(std::move(firstNodeID), std::move(secondNodeID), std::move(roadType), std::move(length));

	return true;
}


/// <summary>
/// Loads one city from the input representation in format:
///		`cityID latitude longitude population`
/// </summary>
/// <param name="map">Object to store info about the node.</param>
/// <param name="line">Input string with city representation.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool MapReader::loadCity(Map& map, const std::string& line)
{
	std::stringstream stream(line);

	int cityID, population;
	double lat, lon;

	// Skip empty line.
	if (std::all_of(line.begin(), line.end(), isspace))
		return true;

	// Bad input format -> Error
	if (!(stream >> cityID >> lat >> lon >> population))
		return false;

	// Correct format -> add city
	// It is ok to rewrite existing city or skip some city (lower index not defined).
	map.AddCity(cityID, population);
	this->graphAdjuster_.AddCityCoordinate(lat, lon);

	return true;
}


/// <summary>
/// Loads one charging station from the input representation in format:
///	`stationID closer_vertexID furtherVertexID segmentID station_capacity 
/// closest_cityID charging_waiting_time estimated_charging_level`
/// </summary>
/// <param name="map">Object to store info about the station.</param>
/// <param name="line">Input string with station representation in format described above.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool MapReader::loadChargingStation(Map& map, const std::string& line)
{
	std::stringstream stream(line);
	
	int stationID, closerVertex, furtherVertex, segment, capacity, closestCityID;
	double chargingWaitingTime, estimatedChargingLevel;

	// Bad input format -> Error
	if (!(stream >> stationID >> closerVertex >> furtherVertex >> segment >> capacity >>
		closestCityID >> chargingWaitingTime >> estimatedChargingLevel))
		return false;

	// Correct format -> add station
	MapPosition stationPosition(closerVertex, furtherVertex, segment, map.GetGraph());
	map.AddChargingStation(ChargingStation(stationID, stationPosition, capacity, 
		closestCityID, chargingWaitingTime, estimatedChargingLevel));

	return true;
}