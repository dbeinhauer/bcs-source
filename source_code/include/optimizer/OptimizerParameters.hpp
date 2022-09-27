#ifndef OPTIMIZERPARAMETERS_H_
#define OPTIMIZERPARAMETERS_H_

#include <iostream>

/// <summary>
/// Parameters for the greedy algorithm.
/// </summary>
struct GreedyParameters
{
public:
	// Number of iterations of the greedy algorithm.
	int32_t MaxIterations_;
	// Number of stations to throw away in each iteration.
	int32_t NumThrowAway_;

	GreedyParameters() 
	{
		this->MaxIterations_ = 0;
		this->NumThrowAway_ = 0;
	}

	GreedyParameters(int32_t maxIterations, int32_t numThrowAway)
	{
		this->MaxIterations_ = maxIterations;
		this->NumThrowAway_ = numThrowAway;
	}
};


/// <summary>
/// Parameters for the genetic algorithm.
/// </summary>
struct GeneticParameters
{
public:
	// Size of the population.
	int32_t PopulationSize_;
	// Number of generations of the algorithm.
	int32_t NumGenerations_;
	// Number of best members to automaticaly choose to next population.
	int32_t NumBestSelection_;
	// Probability of choosing the better member from tournament selection.
	double TournamentSelectionTreshold_;
	// Probability of the mutation of one station.
	double MutationTreshold_;
	// Variance of the member size (for mutation).
	double MemberSizeVariance_;

	GeneticParameters() 
	{
		this->PopulationSize_ = 0;
		this->NumGenerations_ = 0;
		this->NumBestSelection_ = 0;
		this->TournamentSelectionTreshold_ = 0;
		this->MutationTreshold_ = 0;
		this->MemberSizeVariance_ = 0;
	}

	GeneticParameters(int32_t populationSize, int32_t numGenerations, 
		int32_t numBestSelection, double tournamentSelectionTreshold,
		double mutationTreshold, double memberSizeVariance)
	{
		this->PopulationSize_ = populationSize;
		this->NumGenerations_ = numGenerations;
		this->NumBestSelection_ = numBestSelection;
		this->TournamentSelectionTreshold_ = tournamentSelectionTreshold;
		this->MutationTreshold_ = mutationTreshold;
		this->MemberSizeVariance_ = memberSizeVariance;
	}

};


/// <summary>
/// Parameters for the K-Means optimalization.
/// </summary>
struct KMeansParameters
{
public:
	// Number of iterations of the K-Means algorithm.
	int32_t NumIterationsOneRun_;
	// Number of generations of the stations (optimized by K-Means).
	int32_t NumGenerations_;

	KMeansParameters() 
	{
		this->NumIterationsOneRun_ = 0;
		this->NumGenerations_ = 0;
	}
	
	KMeansParameters(int32_t numIterationsOneRun, int32_t numGenerations)
	{
		this->NumIterationsOneRun_ = numIterationsOneRun;
		this->NumGenerations_ = numGenerations;
	}
};



class OptimizerParameters
{
public:
	// Parameters for greedy optimization.
	GreedyParameters GreedyParameters_;
	// Parameters for genetic optimization.
	GeneticParameters GeneticParameters_;
	// Parameters for K-Means optimization.
	KMeansParameters KMeansParameters_;

	// Treshold of score difference to continue in optimization 
	// (in some algorithms).
	double ScoreDifferenceTreshold_;

	// Loss parameters:
	double StationNumberParameter_;
	double RunDownParameter_;
	double DurationParameter_;
	double BatteryDifferenceParameter_;
	double WaitingTimesParameter_;

	OptimizerParameters();
	OptimizerParameters(GreedyParameters greedyParameters, 
		GeneticParameters geneticParameters, KMeansParameters kMeansParameters,
		double scoreDifference, double stationNumberParameter, 
		double runDownParameter, double durationParameter,
		double batteryDifferenceParameter, double waitingTimesParameter);
};


#endif // !OPTIMIZERPARAMETERS_H_

