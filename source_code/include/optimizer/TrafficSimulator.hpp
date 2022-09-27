#ifndef TRAFFICSIMULATOR_H_
#define TRAFFICSIMULATOR_H_


#include <random>
#include <algorithm>
#include <boost/math/distributions/exponential.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>


#include "CarIterator_header.hpp"
#include "SimulationParameters.hpp"


#include "MapReader.hpp"
#include "Map.hpp"
#include "ChargingStation.hpp"


typedef boost::random::discrete_distribution<> DiscreteDistribution;


/// <summary>
/// Class for managing simulator operations.
/// </summary>
class TrafficSimulator
{
public:
	TrafficSimulator(int32_t numClosestStations, double segmentLength);

	bool LoadMap(SimulationParameters& simulationParameters);

	void ResetSimulation();

	// Random generations:
	void GenerateChargingStations(
		int32_t numStations, SimulationParameters& simulationParameters);
	ChargingStation GenerateChargingStation(
		int32_t stationID, SimulationParameters simulationParameters);
	int64_t GenerateCar(
		double startTime, SimulationParameters& simulationParameters);
	double GenerateNextDepartureTime();
	double GenerateBatteryTreshold();


	void DeleteCar(int64_t CarID);

	int64_t GetNumCars();

	// Getting info:
	Map& GetMap();
	CarIterator GetCarIterator(int64_t carID);

private:
	// All map information.
	Map map_;
	// Map of all currently active `Cars`.
	CarMap allCars_;
	// Highest ID of the used car.
	int64_t highestCarID_;

	// Total population of all map cities.
	int32_t totalPopulation_;

	// Probability of each city (based on the population).
	std::vector<double> citiesProbabilities_;
	// Probability of each start-end city pair
	// (first item is start city, second is end city).
	Double2DMatrix startEndCitiesProbabilities_;

	// Random generators:
	std::random_device rd;
	std::mt19937 gen;
	boost::mt19937 gen_;

	// Distributions:
	// Distribution of the start city (based on the city population).
	DiscreteDistribution startCityDistribution_;
	// Distribution of the end city (based on the city population and 
	// distance from the start).
	std::vector<DiscreteDistribution> endCityDistributions_;
	// Distribution to generate next car generation time or waiting time of 
	// the car in the end (before returning back to start).
	std::exponential_distribution<> timeDistribution_;
	// Distribution to compute battery treshold in not certain case 
	// (not certain what to do with the vehicle).
	std::exponential_distribution<> batteryTresholdDistribution_;
	// Distribution of car battery level 
	std::normal_distribution<double> carBatteryDistribution_;


	void initRandomGenerators(double lambdaDepartures,
		double carBatteryMean,
		double carBatteryDeviation,
		double lambdaBatteryTreshold);


	// Probability computations:
	int32_t precomputeTotalPopulation();
	void computeCityProbabilities();
	std::vector<double> computeEndCitiesProbabilites(int32_t startCityID,
		double exponentiaImportance, double lambdaCities);
	void computeAllPairsStartEndCityProbabilities(
		double distanceImportance, double lambdaCities);

	// Random generators:
	int32_t generateRandomCityID();
	int32_t generateRandomFinalCityID(int32_t startCityID);
	vertex_t generateRandomVertexInCity(int32_t cityID);
	vertex_t generateRandomNeighbor(vertex_t vertex);
	int32_t generateRandomSegment(int32_t numSegments);
	MapPosition generateRandomCityMapPosition(int32_t cityID);
	double generateRandomBatteryLevel(double batteryBottomLimit);
	double getExponentialProbability(double value, double lambda);
	Car generateCarObject(double startTime, double carConsumption,
		double batteryBottomLimit, double relativeSpeed);

};

#endif