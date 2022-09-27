#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <tuple>


#include "TimeTable.hpp"
#include "StationParameters.hpp"
#include "OptimizerParameters.hpp"
#include "ModelRepresentation.hpp"


// To compute aberage of the vector objects.
template<typename T>
double getAverage(std::vector<T> const& v) {
	if (v.empty()) {
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}


/// <summary>
/// Class to optimize number of charging stations and its locations.
/// </summary>
class Optimizer
{
public:
	Optimizer(SimulationParameters simulationParameters, 
		OptimizerParameters optimizerParamters);

	void SaveBestModel(std::ostream& stationsStream);

	double ModelLoss(double stationNumberParameter, double runDownParameter, 
		double durationParameter, double batteryDifferenceParameter, 
		double waitingTimesParameter);

	// Run optimizations:
	void RunMultipleSimulations(int16_t numIterations, bool logs);
	void GreedyAlgorithm(bool logs);
	void GeneticAlgorithm(bool logs);
	void KMeansAlgorithm(bool logs);

private:
	// Object to manage simulations.
	TimeTable timeTable_;
	// Parameters of the simulator.
	SimulationParameters simulationParameters_;
	// Parameters of the optimizer.
	OptimizerParameters optimizerParameters_;

	// Current number of charging stations (for some optimization algorithms).
	int32_t numStations_;
	// Loss of the last model.
	double lastLoss_;

	// Random generator variables.
	std::random_device rd;
	std::mt19937 gen;
	boost::mt19937 gen_;


	// Best loss parameters:
	double bestLoss_;
	int32_t bestNumStations_;
	int32_t bestCarsBatteryRunDown_;
	double bestAverageTravelDuration_;
	double bestAverageBatteryDifference_;
	double bestAverageChargingWaitingTimes_;
	int32_t bestNumReturned_;
	int32_t bestUsedVehicles_;
	int32_t bestRunDownCharging_;
	int32_t bestNumGoingCharing_;
	int32_t bestNotFinished_;
	int32_t bestNonNullStations_;
	int32_t bestTotalStationCustomers_;
	ModelRepresentation bestModel;


	// Loss paramters sum for averaging.
	double lossSum_;
	int64_t batteryRunDownSum_;
	double averageTravelDurationSum_;
	double averageBatteryDifferemceSum_;
	double averageChargingWaitingTimesSum_;


	void initGenerators();

	void printBestLossInfo();

	// General operations on the model:
	void addChoosenStations(
		std::vector<StationParameters>& stationsToMaintain, bool updateSimulator);
	bool randomlyChooseNumberOfStationsChange(double changeProbability);

	// Greedy optimization functions:
	void greedyFillRestStations(bool updateSimulator);
	void greedyChooseStations(
		std::vector<StationParameters>& stationsToMaintain, int32_t numThrowAway);
	void greedyTakeUsedChargingStations(std::vector<StationParameters>& stationsToMaintain);
	void greedyGenerateStationsBatteryRunDown(
		std::vector<StationParameters>& stationsToMaintain, int32_t numThrowAway);

	// Genetic optimization functions:
	std::vector<ModelRepresentation> geneticGetNewPopulation(
		std::vector<ModelRepresentation>& population, int32_t numBest);
	void geneticRandomlyInitPopulation(std::vector<ModelRepresentation>& population);
	int32_t geneticTournamentSelection(
		std::vector<ModelRepresentation>& population, double betterTreshold);
	std::pair<ModelRepresentation, ModelRepresentation> geneticPointCrossover(
		std::vector<ModelRepresentation>& population, 
		int32_t firstParentID, int32_t secondParentID);
	int32_t geneticRandomlyChooseNumberOfStationsChange(int32_t mean, double variance);
	void geneticMutation(ModelRepresentation& model, double mutationProbability);
	void geneticMutationRandomlyChangeStations(
		ModelRepresentation& model, double mutationProbability);
	void geneticMutationChangeNumberOfStations(ModelRepresentation& model, 
		int32_t numberStationsMean, double numberStationsVariance);
	
	// K-Means optimization functions:
	void kMeansFindClusters(std::vector<std::vector<vertex_t>>& allClusters,
		std::vector<MapPosition>& allCentroidsPositions);
	void kMeansComputeVerticesWeights(std::map<vertex_t, double>& verticesWeights);
	void kMeansFindNewCentroids(std::vector<std::vector<vertex_t>>& allClusters,
		std::vector<MapPosition>& allCentroidsPositions, 
		std::map<vertex_t, double>& verticesWeights);
	MapPosition kMeansFindNewCentroid(
		std::vector<vertex_t>& cluster, std::map<vertex_t, double>& verticesWeights);
};


#endif
