#include "optimizer/OptimizerParameters.hpp"


/// <summary>
/// Initializes the optimizer parameters.
/// </summary>
OptimizerParameters::OptimizerParameters()
{
	this->ScoreDifferenceTreshold_ = 10;

	// Loss parameters:
	this->StationNumberParameter_ = 1000;
	this->RunDownParameter_ = 100;
	this->DurationParameter_ = 0.01;
	this->BatteryDifferenceParameter_ = 0.001;
	this->WaitingTimesParameter_ = 0.1;
}


/// <summary>
/// Initializes the optimizer parameters.
/// </summary>
/// <param name="greedyParameters">Parameters for greedy optimization.</param>
/// <param name="geneticParameters">Parameters for genetic optimization.</param>
/// <param name="kMeansParameters">Parameters for K-Means optimization.</param>
/// <param name="scoreDifference">Treshold of score difference to continue in 
/// optimization (in some algorithms).</param>
/// <param name="stationNumberParameter">Loss parameter for number of stations.</param>
/// <param name="runDownParameter">Loss parameter for run down vehicles.</param>
/// <param name="durationParameter">Loss parameter for travel duration.</param>
/// <param name="batteryDifferenceParameter">Loss parameter for battery difference.</param>
/// <param name="waitingTimesParameter"><Loss parameter for waiting times in stations./param>
OptimizerParameters::OptimizerParameters(GreedyParameters greedyParameters, 
	GeneticParameters geneticParameters, KMeansParameters kMeansParameters, 
	double scoreDifference, double stationNumberParameter, 
	double runDownParameter, double durationParameter,
	double batteryDifferenceParameter, double waitingTimesParameter)
{
	this->GreedyParameters_ = greedyParameters;
	this->GeneticParameters_ = geneticParameters;
	this->KMeansParameters_ = kMeansParameters;

	this->ScoreDifferenceTreshold_ = scoreDifference;

	// Loss parameters:
	this->StationNumberParameter_ = stationNumberParameter;
	this->RunDownParameter_ = runDownParameter;
	this->DurationParameter_ = durationParameter;
	this->BatteryDifferenceParameter_ = batteryDifferenceParameter;
	this->WaitingTimesParameter_ = waitingTimesParameter;
}