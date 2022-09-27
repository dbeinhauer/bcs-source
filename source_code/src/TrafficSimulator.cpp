#include "optimizer/TrafficSimulator.hpp"


/// <summary>
/// Intitilizes the object.
/// </summary>
/// <param name="numClosestStations">Number of closest charging stations to 
/// consider during optimal path finding.</param>
/// <param name="segmentLength">Length of the edge segment.</param>
TrafficSimulator::TrafficSimulator(
	int32_t numClosestStations, double segmentLength)
	:map_(numClosestStations, segmentLength)
{
	this->totalPopulation_ = 0;
	this->highestCarID_ = 0;
}

/// <summary>
/// Load map and stores all other necessary precomputed data 
/// (precomputed distances between each pair of cities).
/// </summary>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <returns>Returns `true` if loading was successful, else `false`.</returns>
bool TrafficSimulator::LoadMap(SimulationParameters& simulationParameters)
{
	MapReader mapReader;

	// Load map
	bool readingSuccesfull = 
		mapReader.LoadMap(this->map_, simulationParameters);

	// Move information about cities distances into the graph
	this->map_.SetCitiesDistances(
		mapReader.GetGraphAdjuster().GetCitiesDistances());

	// Compute total population (for probabilistic generator).
	this->totalPopulation_ = this->precomputeTotalPopulation();

	// Compute probability of each city (for car generation).
	this->computeCityProbabilities();
	// Compute probability of all start and end city pair of cities 
	// (for path generating).
	this->computeAllPairsStartEndCityProbabilities(
		simulationParameters.EndCityRatio,
		simulationParameters.ExponentialLambdaCities);
	// Initialize random generators.
	this->initRandomGenerators(simulationParameters.ExponentialLambdaDepartures,
		simulationParameters.CarBatteryMean, 
		simulationParameters.CarBatteryDeviation,
		simulationParameters.BatteryTresholdLambda);
	
	// Return reading status.
	return readingSuccesfull;
}


/// <summary>
/// Reset simulation parameters.
/// </summary>
void TrafficSimulator::ResetSimulation()
{
	this->allCars_.clear();
	this->highestCarID_ = 0;

	this->map_.ResetSimulation();
}


/// <summary>
/// Randomly generates charging stations (based on city probability 
/// and its distance from the center of the city) and stores it.
/// </summary>
/// <param name="numStations">Number of stations to generate.</param>
/// <param name="simulationParameters">Parameters of the simulation.</param>
void TrafficSimulator::GenerateChargingStations(
	int numStations, SimulationParameters& simulationParameters)
{
	// Get container of all charging stations from the map and generate them.
	auto& allChargingStations = this->map_.GetChargingStations();

	// Offset of the station ID (number of already inserted charging stations).
	int32_t offsetStationID = allChargingStations.size();

	for (int32_t i = 0; i < numStations; i++)
	// Randomly generate given number of stations and store them.
	{
		// Compute station ID.
		int32_t stationID = i + offsetStationID;
		
		// Generate charging station.
		allChargingStations.push_back(std::make_unique<ChargingStation>(
				this->GenerateChargingStation(
					stationID, simulationParameters)));
	}

	// Find nearest charging stations for each city.
	this->map_.FindCitiesNearestChargingStations();
}


/// <summary>
/// Randomly generates charging station.
/// </summary>
/// <param name="stationID">ID of the station.</param>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <returns></returns>
ChargingStation TrafficSimulator::GenerateChargingStation(
	int32_t stationID, SimulationParameters simulationParameters)
{
	// Generate city ID.
	int32_t cityID = this->generateRandomCityID();

	return ChargingStation(stationID, 
		this->generateRandomCityMapPosition(cityID), 
		simulationParameters.StationCapacity,
		cityID,
		simulationParameters.ChargingWaitingTime,
		simulationParameters.MeanChargingLevel);
}


/// <summary>
/// Generates car and adds it to the vector of all currently used cars.
/// </summary>
/// <returns>Returns ID of the new car.</returns>
int64_t TrafficSimulator::GenerateCar(
	double startTime, SimulationParameters& simulationParameters)
{
	// Generate car object
	Car car = this->generateCarObject(startTime,
		simulationParameters.CarConsumption,
		simulationParameters.CarStartBatteryBottomLimit,
		simulationParameters.CarVelocity);


	// Car ID is the lowest non-used car ID.
	int64_t carID = this->highestCarID_;
	
	// Update next car ID.
	this->highestCarID_++;


	// Computes currently fastest path.
	car.UpdateVehiclePath(car.GetPosition().CloserVertexID,
		car.GetEndPosition().CloserVertexID,
		this->map_,
		this->GenerateBatteryTreshold(),
		false,
		simulationParameters);


	// Add car to the vector of all current cars.
	auto carIterator = 
		this->allCars_.emplace(std::make_pair(carID, std::move(car)));

	if (!carIterator.second)
	// Car couldn't be added to the  vector -> Error
	{
		std::cout << "Bad Car ID: " << std::endl;
		throw "Bad id exception";
	}

	return carID;
}

/// <summary>
/// Deletes the car from the container of all currently used cars.
/// </summary>
/// <param name="carID">ID of the car to delete.</param>
void TrafficSimulator::DeleteCar(int64_t carID)
{
	if (this->allCars_.find(carID) == this->allCars_.end())
	// Car ID is not in the container of all currently used cars 
	// -> probable error
	{
		std::cout << "Error: CarID " << carID << " doesn't exist" << std::endl;
		return;
	}

	this->allCars_.erase(this->allCars_.find(carID));
}


int64_t TrafficSimulator::GetNumCars()
{
	return this->allCars_.size();
}



/// <summary>
/// </summary>
/// <returns>Returns reference to map.</returns>
Map& TrafficSimulator::GetMap()
{
	return this->map_;
}

/// <summary>
/// Finds queried car and returns its reference as `CarIterator` value.
/// </summary>
/// <param name="carID">ID of the car to find if `-1`, then return end of 
/// the container (we want to propagate information that this object 
/// doesn't exist). </param>
/// <returns>Returns iterator to the given car into the container of all 
/// currently used cars or `this->allCars_.end()` if this object doesn't exist.</returns>
CarIterator TrafficSimulator::GetCarIterator(int64_t carID)
{
	if (carID == -1)
	// Searched car doesn't exist (we know it already when we call the function)
	// -> return the `this->allCars_.end()` iterator
	{
		return this->allCars_.end();
	}

	// Find iterator of the car.
	CarIterator carIterator = this->allCars_.find(carID);
	if (carIterator == this->allCars_.end())
	// Searched car doesn't exist -> probable error.
	{
		std::cout << "Car does not exist" << std::endl;
		throw "Car does not exist error";
	}


	return carIterator;
}


/// <summary>
/// Initializes random generators used by the object.
/// </summary>
/// <param name="lambdaDepartures">Parameter of exponential distribution of 
/// car departure time.</param>
/// <param name="carBatteryMean">Mean of the car battery level at the 
/// start of the route.</param>
/// <param name="carBatteryDeviation">Deviation of the car battery 
/// level at the start of the route.</param>
/// <param name="lambdaBatteryTreshold">Paremeter of exponential distribution of
/// battery treshold to consider battery charging.</param>
void TrafficSimulator::initRandomGenerators(double lambdaDepartures, 
	double carBatteryMean, double carBatteryDeviation, 
	double lambdaBatteryTreshold)
{
	this->gen = std::mt19937(rd());

	// Init start city distribution (weighted by the population).
	this->startCityDistribution_ = 
		DiscreteDistribution(this->citiesProbabilities_);

	// Compute end city probability distribution
	// (counts with population and distance from the start).
	for (size_t i = 0; i < this->citiesProbabilities_.size(); i++)
	{
		this->endCityDistributions_.push_back(
			DiscreteDistribution(this->startEndCitiesProbabilities_[i]));
	}

	// Init distribution of the departure and waiting times.
	this->timeDistribution_ =
		std::exponential_distribution<>(lambdaDepartures);

	// Init distribution of the start battery level.
	this->carBatteryDistribution_ =
		std::normal_distribution<double>(carBatteryMean, carBatteryDeviation);

	// Init distribution of the battery treshold (to consider charging).
	this->batteryTresholdDistribution_ = 
		std::exponential_distribution<>(lambdaBatteryTreshold);

}


/// <summary>
/// Computes total population of map (for further probabilistic car deployment).
/// </summary>
/// <returns>Returns total population of the map (sum of all cites).</returns>
int32_t TrafficSimulator::precomputeTotalPopulation()
{
	auto cities = this->map_.GetCitiesPopulation();
	int32_t totalPopulation = 0;
	for (auto&& city : cities)
	{
		totalPopulation += city;
	}

	return totalPopulation;
}



/// <summary>
/// Computes probability of city for each city (infered from population).
/// </summary>
void TrafficSimulator::computeCityProbabilities()
{
	auto citiesPopulations = this->map_.GetCitiesPopulation();
	for (size_t i = 0; i < citiesPopulations.size(); i++)
	{
		this->citiesProbabilities_.push_back(std::log(citiesPopulations[i]));
	}
}


/// <summary>
/// Computes probability distribution of each city based on its population and 
/// its distance from the start city (exponential distribution) based on formula:
///		`probability = city_probability * (exponentialImportaneParameter * distance_probability)
/// </summary>
/// <param name="startCityID">ID of the start city.</param>
/// <param name="exponentiaImportance">Importance parameter of the distance of 
/// the city versus the city population.
/// ( `< 1` - less important, `== 1` - same importance, `> 1` - more important).</param>
/// <param name="lambdaCities">Parameter of the exponential distribution 
/// of the city distance.</param>
/// <returns>Returns probability distribution of the end cities based on 
/// the start city.</returns>
std::vector<double> TrafficSimulator::computeEndCitiesProbabilites(
	int32_t startCityID, double distanceImportance, double lambdaCities)
{
	std::vector<double> citiesProbabilities(this->citiesProbabilities_.size());

	// Sum of the pseudo probabilites, to normalize them to get real probabilities.
	double pseudoProbabilitySum = 0;


	// Compute pseudo probabilities of final city as:
	//		`city_probability + (parameter * city_distance_probability)`
	// Not normalized.
	for (size_t i = 0; i < this->citiesProbabilities_.size(); i++)
	{
		double pseudoProbability = this->citiesProbabilities_[i] +
			(distanceImportance *
				this->getExponentialProbability(
					this->map_.GetCityPairDistance(startCityID, i),
					lambdaCities));

		citiesProbabilities[i] = pseudoProbability;
		pseudoProbabilitySum += pseudoProbability;
	}

	// Normalize pseudo probabilities to get actual probabilites.
	for (size_t i = 0; i < this->citiesProbabilities_.size(); i++)
	{
		citiesProbabilities[i] /= pseudoProbabilitySum;
	}

	return std::move(citiesProbabilities);
}


/// <summary>
/// Computes probability that the city is end city if known start city 
/// (for each start-end city pair). Be aware that this probability is not 
/// symetric (it is conditional probability)!
/// </summary>
/// <param name="exponentiaImportance">Importance parameter of the distance of 
/// the city versus the city population.
/// ( `< 1` - less important, `== 1` - same importance, `> 1` - more important).</param>
/// <param name="lambdaCities">Parameter of the exponential distribution 
/// of the city distance.</param>
void TrafficSimulator::computeAllPairsStartEndCityProbabilities(
	double distanceImportance, double lambdaCities)
{
	this->startEndCitiesProbabilities_ = Double2DMatrix(this->citiesProbabilities_.size());

	for (size_t i = 0; i < this->citiesProbabilities_.size(); i++)
	{
		this->startEndCitiesProbabilities_[i] =
			this->computeEndCitiesProbabilites(i, distanceImportance, lambdaCities);
	}
}



/// <summary>
/// Randomly generates city using uniform distribution 
/// with weights proportional to its population.
/// </summary>
/// <returns>Returns ID of the randomly generated city.</returns>
int32_t TrafficSimulator::generateRandomCityID()
{
	return this->startCityDistribution_(this->gen_);
}


/// <summary>
/// Randomly generates end city with distribution that counts with city 
/// population and its distance from the start city.
/// </summary>
/// <param name="startCityID">ID of the start city.</param>
/// <returns>Returns randomly generated ID of the end city.</returns>
int32_t TrafficSimulator::generateRandomFinalCityID(int32_t startCityID)
{
	return this->endCityDistributions_[startCityID](this->gen_);
}


/// <summary>
/// Uniformly randomly generates segment ID of the given edge.
/// </summary>
/// <param name="numSegments">Number of segments of the edge.</param>
/// <returns>Returns randomly generated segment index.</returns>
int32_t TrafficSimulator::generateRandomSegment(int32_t numSegments)
{
	return std::rand() % numSegments;
}

/// <summary>
/// Uniformly randomly generates node index from the given city.
/// </summary>
/// <param name="cityID">City ID to generate node from.</param>
/// <returns>Returns randomly generated node ID.</returns>
vertex_t TrafficSimulator::generateRandomVertexInCity(int32_t cityID)
{
	auto cityNodes = this->map_.GetCityNodes(cityID);
	
	while (cityNodes.size() == 0)
	// To make sure that city has at least 1 node 
	// (low probability in real example - probably never happens).
	// When there are no available vertex -> generate vertex from other city
	{
		cityID = this->generateRandomCityID();
		cityNodes = this->map_.GetCityNodes(cityID);
	}

	int32_t randomIndex = std::rand() % (cityNodes.size());
	return cityNodes[randomIndex];
}


/// <summary>
/// Uniformly randomly generates neighbor of the vertex.
/// </summary>
/// <param name="vertex">Vertex identificator to get its neighbor.</param>
/// <returns>Return ID of the randomly generated neighbor or 
/// `-1` if error occured.</returns>
vertex_t TrafficSimulator::generateRandomNeighbor(vertex_t vertex)
{
	const Graph& graph = this->map_.GetConstGraph();
	// All neighbors of the vertex.
	auto neighbors = boost::adjacent_vertices(vertex, graph);
	int32_t numNeighbors = std::distance(neighbors.first, neighbors.second);
	int32_t randomIndex = std::rand() % numNeighbors;

	// Iterate through the adjacency list to get the generated neighbor.
	int32_t actIndex = 0;
	for (auto it = neighbors.first; it < neighbors.second; it++)
	{
		if (actIndex == randomIndex)
			return *it;

		actIndex++;
	}

	// Error ocurred.
	return -1;
}


/// <summary>
/// Randomly generates position in map belonging to given city.
/// </summary>
/// <param name="cityID">ID of the city to generate position.</param>
/// <returns>Returns randomly generated `MapPosition` object.</returns>
MapPosition TrafficSimulator::generateRandomCityMapPosition(int32_t cityID)
{
	// Generate random edge.
	vertex_t vertex = this->generateRandomVertexInCity(cityID);
	vertex_t neighborVertex = this->generateRandomNeighbor(vertex);
	edge_t edge =
		boost::edge(vertex, neighborVertex, this->map_.GetGraph()).first;

	if (this->map_.GetGraph()[edge].GetNumSegments() == 0)
	// Number of segments of the edge is zero -> error happened
	{
		std::cout << "Number of segments is 0 in edge: " << edge << std::endl;
		throw "Number of segments is 0 in edge: ";
	}

	// Generate random segment.
	int32_t segmentID = this->generateRandomSegment(
		this->map_.GetGraph()[edge].GetNumSegments());

	return MapPosition(vertex, neighborVertex, segmentID, this->map_.GetGraph());

}


/// <summary>
/// Randomly generates battery level from normal distribution.
/// </summary>
/// <param name="batteryBottomLimit">Bottom limit of the battery level.</param>
/// <returns>Randomly generated battery level.</returns>
double TrafficSimulator::generateRandomBatteryLevel(double batteryBottomLimit)
{
	double batteryLevel = this->carBatteryDistribution_(gen);

	// If generated level higher than 100 % -> return full battery level.
	if (batteryLevel > 1)
		return 1;

	// If battery level lower than bottom limit -> return bottom limit
	if (batteryLevel < batteryBottomLimit)
		return batteryBottomLimit;

	// Return generated battery level.
	return batteryLevel;
}


/// <summary>
/// Generates random additional part to battery treshold 
/// (for random part of the charging decision).
/// </summary>
/// <returns>Returns randomly additional battery level above the bottom treshold.</returns>
double TrafficSimulator::GenerateBatteryTreshold()
{
	return this->batteryTresholdDistribution_(gen);
}


/// <summary>
/// Returns probability of the travel distance computed 
/// from the pmf of the exponential distribution.
/// </summary>
/// <param name="value">Distance to compute probability from.</param>
/// <param name="lambda">Parameter of the exponential distribution.</param>
/// <returns>Returns probability of the `value` in choosen 
/// exponentional distribution.</returns>
double TrafficSimulator::getExponentialProbability(double value, double lambda)
{
	auto distr = boost::math::exponential_distribution<>{ lambda };
	return boost::math::pdf(distr, value);
}


/// <summary>
/// Creates `Car` object and randomly generates its start and end position.
/// </summary>
/// <param name="startTime">Time of the start of the car.</param>
/// <param name="carConsumption">Battery consumption of the car.</param>
/// <param name="batteryBottomLimit">Battery bottom limit of the car.</param>
/// <param name="relativeSpeed">Speed of the car (in km/minute).</param>
/// <returns>Returns randomly generated `Car` object.</returns>
Car TrafficSimulator::generateCarObject(double startTime, 
	double carConsumption, double batteryBottomLimit, double relativeSpeed)
{
	int startCityID, endCityID;
	MapPosition startPosition, endPosition;

	// Generate start and end city.
	startCityID = this->generateRandomCityID();
	endCityID = this->generateRandomFinalCityID(startCityID);

	// std::cout << startCityID << "    " << endCityID << std::endl; 

	// Generate start and end position.
	startPosition = this->generateRandomCityMapPosition(startCityID);
	endPosition = this->generateRandomCityMapPosition(endCityID);

	// std::cout << 
	// 	this->map_.FindShortestPath(startPosition.CloserVertexID,
	// 	 endPosition.CloserVertexID).first << std::endl;

	// Generate waiting time and battery level.
	double waitingTime = this->GenerateNextDepartureTime();
	double batteryLevel = this->generateRandomBatteryLevel(batteryBottomLimit);

	// std::cout << batteryLevel << std::endl;


	return Car(std::move(startTime), std::move(waitingTime), 
		std::move(startPosition), std::move(endPosition), carConsumption,
		batteryLevel, relativeSpeed);
}


/// <summary>
/// Generates next departure or waiting time.
/// </summary>
/// <returns>Returns exponentialy randomly generated time.</returns>
double TrafficSimulator::GenerateNextDepartureTime()
{
	return timeDistribution_(gen);
}
