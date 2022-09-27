#include "optimizer/Optimizer.hpp"


/// <summary>
/// Initilizes the optimizer object.
/// </summary>
/// <param name="simulationParameters">Parameters of the traffic simulator.</param>
/// <param name="optimizerParameters">Parameters of the optimizer.</param>
Optimizer::Optimizer(SimulationParameters simulationParameters,
	OptimizerParameters optimizerParameters)
	:timeTable_(simulationParameters), 
	numStations_(simulationParameters.NumStations),
	simulationParameters_(std::move(simulationParameters)),
	optimizerParameters_(std::move(optimizerParameters))
{
	this->lastLoss_ = 0;
	this->timeTable_.LoadMap(this->simulationParameters_);

	// Best loss information init:
	this->bestLoss_ = -1;
	this->bestNumStations_ = 0;
	this->bestCarsBatteryRunDown_ = 0;
	this->bestAverageTravelDuration_ = 0;
	this->bestAverageBatteryDifference_ = 0;
	this->bestAverageChargingWaitingTimes_ = 0;
	this->bestNumReturned_ = 0;
	this->bestUsedVehicles_ = 0;
	this->bestRunDownCharging_ = 0;
	this->bestNumGoingCharing_ = 0;
	this->bestNotFinished_ = 0;
	this->bestNonNullStations_ = 0;
	this->bestTotalStationCustomers_ = 0;

	this->bestModel = ModelRepresentation(-1, std::vector<StationParameters>());

	// Loss info sums init:
	this->lossSum_ = 0;
	this->batteryRunDownSum_ = 0;
	this->averageTravelDurationSum_ = 0;
	this->averageBatteryDifferemceSum_ = 0;
	this->averageChargingWaitingTimesSum_ = 0;
}


void Optimizer::SaveBestModel(std::ostream& stationsStream)
{
	MapReader mapReader;
	mapReader.WriteStations(this->bestModel, this->timeTable_.GetTrafficSimulator().GetMap(), stationsStream);
}



/// <summary>
/// Computes loss of the current station model.
/// Score is just linear combination of: 
/// 
/// (number_of_charging_stations,
/// number_of_battery_run_downs, 
/// average_travel_duration,
/// average_battery_difference (difference between start and end battery level 
/// (negative - battery increase, positive - battery decrease)),
/// average_waiting_time)
/// 
/// Each part of the linear combination has its apropriate multiplication contants.
/// </summary>
/// <param name="stationNumberParameter">Multiplication constant of 
/// the number of charging stations.</param>
/// <param name="runDownParameter">Multiplication constant of 
/// the number of  run down vehicles.</param>
/// <param name="durationParameter">Multiplication constant of 
/// the average duration (in minutes).</param>
/// <param name="batteryDifferenceParameter">Multiplication constant of 
/// the average battery level difference (values from -1 to 1 (negative values 
/// means battery level increases, positive - decreases)).</param>
/// <param name="waitingTimesParameter">Multiplication constant for average 
/// waiting time in the charging station.</param>
/// <returns>Returns loss value of the current charging stations model.</returns>
double Optimizer::ModelLoss(double stationNumberParameter, double runDownParameter,
	double durationParameter, double batteryDifferenceParameter, double waitingTimesParameter)
{
	// Compute parameters needed to compute loss.

	auto& durations = this->timeTable_.GetTravelDurations();

	double averageTravelDuration = getAverage(this->timeTable_.GetTravelDurations());
	double averageBatteryDifference = getAverage(this->timeTable_.GetBatteryDifferences());
	double averageChargingWaitingTimes = 0;

	// Comute average waiting times on the station.
	auto& waitingTimes = this->timeTable_.GetWaitingTimes();
	int32_t numCustomers = 0;
	double waitingTimesSum = 0;

	int32_t nonNullStations = 0;
	for (auto&& stationPair : waitingTimes)
	{
		if (stationPair.first > 0)
			nonNullStations++;

		numCustomers += stationPair.first;
		waitingTimesSum += stationPair.second;
		std::cout << "Number of waitings: " << stationPair.first << "    Avg. waiting time: " << stationPair.second / stationPair.first << std::endl;
	}

	// Dividing by zero check.
	if (numCustomers != 0)
		averageChargingWaitingTimes = waitingTimesSum / numCustomers;


	// Loss computation.
	double loss = this->numStations_ * stationNumberParameter
		+ this->timeTable_.CarsBatteryRunDown_ * runDownParameter
		+ averageTravelDuration * durationParameter
		+ averageBatteryDifference * batteryDifferenceParameter
		+ averageChargingWaitingTimes * waitingTimesParameter;


	// Loss logs:
	std::cout << "Number of stations: " << this->numStations_ << std::endl;
	std::cout << "Number of run down batteries: " << this->timeTable_.CarsBatteryRunDown_ << std::endl;
	std::cout << "Average traveling duration: " << averageTravelDuration <<" minutes" << std::endl;
	std::cout << "Average battery difference between start and end: " << averageBatteryDifference << std::endl;
	std::cout << "Average waiting times in charging station: " << averageChargingWaitingTimes << " minutes" << std::endl;
	std::cout << "Number of vehicles which reached the end: " << this->timeTable_.GetNumReturnedVehicles() << std::endl;
	std::cout << "Number of vehicles which didn't finish before end of simulation: " << this->timeTable_.NumNotFinishedCars << std::endl;
	std::cout << "Total number of used charging stations: " << nonNullStations << std::endl;
	// std::cout << "Total number of charging stations usage: " << numCustomers << std::endl;
	std::cout << "Run down batteries during going to station: " << this->timeTable_.RunDownCharging << std::endl;
	std::cout << "Total number of vehicles which didn't reach the end: " << this->timeTable_.CarsBatteryRunDown_ + this->timeTable_.NumNotFinishedCars << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "Model loss: " << loss << std::endl;
	std::cout << std::endl;


	// Update loss sum:
	this->lossSum_ += loss;
	this->batteryRunDownSum_ += this->timeTable_.CarsBatteryRunDown_;
	this->averageTravelDurationSum_ += averageTravelDuration;
	this->averageBatteryDifferemceSum_ += averageBatteryDifference;
	this->averageChargingWaitingTimesSum_ += averageChargingWaitingTimes;


	if (this->bestLoss_ < 0 || this->bestLoss_ > loss)
	// Update total best loss info.
	{
		this->bestLoss_ = loss;
		this->bestNumStations_ = this->numStations_;
		this->bestCarsBatteryRunDown_ = this->timeTable_.CarsBatteryRunDown_;
		this->bestAverageTravelDuration_ = averageTravelDuration;
		this->bestAverageBatteryDifference_ = averageBatteryDifference;
		this->bestAverageChargingWaitingTimes_ = averageChargingWaitingTimes;
		this->bestNumReturned_ = this->timeTable_.GetNumReturnedVehicles();
		this->bestUsedVehicles_ = nonNullStations;
		this->bestRunDownCharging_ = this->timeTable_.RunDownCharging;
		this->bestNumGoingCharing_ = this->timeTable_.NumCharging_;
		this->bestNotFinished_ = this->timeTable_.NumNotFinishedCars;
		this->bestNonNullStations_ = nonNullStations;
		this->bestTotalStationCustomers_ = numCustomers;
	}


	return loss;
}


/// <summary>
/// Runs traffic simulation multiple times with randomly placed charing stations.
/// </summary>
/// <param name="numIterations">Number of simulations to run.</param>
/// <param name="logs">Whether to write simulator logs.</param>
void Optimizer::RunMultipleSimulations(int16_t numIterations, bool logs)
{
	for (int32_t i = 0; i < numIterations; i++)
	// For each iteration generate new charging stations and run the simulation.
	{
		this->timeTable_.ResetSimulation();

		this->timeTable_.RandomlyGenerateChargingStations(
			this->numStations_,
			this->simulationParameters_);

		this->timeTable_.RunSimulation(
			this->numStations_, this->simulationParameters_, logs);


		// Compute current loss.
		double currentLoss =
			this->ModelLoss(this->optimizerParameters_.StationNumberParameter_,
				this->optimizerParameters_.RunDownParameter_,
				this->optimizerParameters_.DurationParameter_,
				this->optimizerParameters_.BatteryDifferenceParameter_,
				this->optimizerParameters_.WaitingTimesParameter_);

		if (currentLoss < this->bestModel.Loss_ || i == 0)
		// Current loss is recently the best or first iteration of the algorithm 
		// -> update best loss information
		{
			auto& map = this->timeTable_.GetTrafficSimulator().GetMap();
			this->bestModel.SaveModelRepresentation(map.GetChargingStations());
			// std::cout << this->bestModel.AllChargingStations_.size() << std::endl;
			this->bestModel.Loss_ = currentLoss;
		}
	}

	// Write average loss info, not the best one (best doesn't make sense).
	this->bestLoss_ = this->lossSum_ / numIterations;
	this->bestCarsBatteryRunDown_ = this->batteryRunDownSum_ / numIterations;
	this->bestAverageTravelDuration_ = 
		this->averageTravelDurationSum_ / numIterations;
	this->bestAverageBatteryDifference_ = 
		this->averageBatteryDifferemceSum_ / numIterations;
	this->bestAverageChargingWaitingTimes_ = 
		this->averageChargingWaitingTimesSum_ / numIterations;

	this->printBestLossInfo();
}


/// <summary>
/// Optimizes charging station positions and its number based on the station 
/// usage and positions of the battery run down vehicles.
/// Maintains all at least once time used stations and rest randomly 
/// generates on the positions of battery run down and also some stations 
/// just doesn't use (to optimize number of stations too).
/// </summary>
/// <param name="logs">Whether to write simulator logs.</param>
void Optimizer::GreedyAlgorithm(bool logs)
{
	// ModelRepresentation bestModel;

	for (int32_t i = 0; i < this->optimizerParameters_.GreedyParameters_.MaxIterations_; i++)
	// Run given number of iteations of the algorithm.
	{
		auto& map = this->timeTable_.GetTrafficSimulator().GetMap();

		if (i == 0)
		// First iteration -> just randomly place the stations
		{
			this->timeTable_.RandomlyGenerateChargingStations(
				this->numStations_,
				this->simulationParameters_);
		}
		else
		// Model tuning (maintains used stations and rest either randomly 
		// generates or doesn't add to the map at all).
		{
			std::vector<StationParameters> stationsToMaintain;
			this->greedyChooseStations(stationsToMaintain, 
				this->optimizerParameters_.GreedyParameters_.NumThrowAway_);

			this->timeTable_.ResetSimulation();

			// Add choosen stations and fill rest with the randomly generated ones.
			this->addChoosenStations(stationsToMaintain, false);
			this->greedyFillRestStations(false);

			// Update info about the number of stations.
			this->numStations_ = map.GetChargingStations().size();

			// Find nearest charging stations for each city.
			map.FindCitiesNearestChargingStations();
		}

		// Run simulation with new charging stations.
		this->timeTable_.RunSimulation(
			this->numStations_, this->simulationParameters_, logs);

		// Compute current loss.
		double currentLoss =
			this->ModelLoss(this->optimizerParameters_.StationNumberParameter_,
				this->optimizerParameters_.RunDownParameter_,
				this->optimizerParameters_.DurationParameter_,
				this->optimizerParameters_.BatteryDifferenceParameter_,
				this->optimizerParameters_.WaitingTimesParameter_);


		if (std::abs(currentLoss -  this->lastLoss_) 
			< this->optimizerParameters_.ScoreDifferenceTreshold_)
		// Loss difference is neglieble -> end optimalization
		{
			this->lastLoss_ = currentLoss;
			std::cout << "Loss didn't change much -> End of the optimalization" 
				<< std::endl;
			break;
		}

		// Update last loss value.
		this->lastLoss_ = currentLoss;

		if (currentLoss < this->bestModel.Loss_ || i == 0)
		// Current loss is recently the best or first iteration of the algorithm 
		// -> update best loss information
		{
			this->bestModel.SaveModelRepresentation(map.GetChargingStations());
			this->bestModel.Loss_ = currentLoss;

			// std::cout << this->bestModel.AllChargingStations_.size() << std::endl;

		}
	}

	std::cout << "Best loss: " << this->bestModel.Loss_ << std::endl;

	this->printBestLossInfo();
}



/// <summary>
/// Optimizes charging stations position and number using the genetic algorithm approach.
/// One model is considered to be an individual. All operations works only with
/// the loss (fitness) and all charging stations of the model. For selection we
/// take given number of best models and rest of the population we choose using
/// the tournament selection. We then use one point crossover (changes stations 
/// from the given index in the vector of all charging stations). And in mutation
/// we randomly delete stations and replace them with randomly generated one or
/// we randomly delete or add some stations (mutation in position 
/// and also in number of stations).
/// </summary>
/// <param name="logs">Whether to write simulator logs.</param>
void Optimizer::GeneticAlgorithm(bool logs)
{
	// Randomly generate first generation of the population.
	std::vector<ModelRepresentation> population;
	this->geneticRandomlyInitPopulation(population);


	// In each vector is a vector representing one generation of the algorithm
	// and losses of its members (for final loss printing).
	std::vector<std::vector<double>> lossesForGenerations;
	// ModelRepresentation bestModel(-1, std::vector<StationParameters>());

	for (int32_t i = 0; i < this->optimizerParameters_.GeneticParameters_.NumGenerations_; i++)
	// Each generation of the genetic algorithm.
	{
		std::cout << "Iteration number: " << i << std::endl;
		std::cout << "------------------------------------" << std::endl;

		double generationBestLoss = -1;
		int32_t generationBestModelID = 0;

		// Help flag to know whether it is the first member of the population.
		bool firstIteration = true;
		for (auto&& model : population)
		// Perform simulation on each model of the current population (to evaluate loss).
		{
			this->addChoosenStations(model.AllChargingStations_, true);

			this->timeTable_.RunSimulation(
				model.AllChargingStations_.size(), this->simulationParameters_, logs);

			// Compute current loss.
			double currentLoss =
				this->ModelLoss(this->optimizerParameters_.StationNumberParameter_,
					this->optimizerParameters_.RunDownParameter_,
					this->optimizerParameters_.DurationParameter_,
					this->optimizerParameters_.BatteryDifferenceParameter_,
					this->optimizerParameters_.WaitingTimesParameter_);

			model.Loss_ = currentLoss;
			if (firstIteration)
			// First member of the population.
			{
				firstIteration = false;
				lossesForGenerations.push_back(std::vector<double>());
				generationBestLoss = currentLoss;
			}
			lossesForGenerations[i].push_back(currentLoss);

			// Current loss is lower than currently the best.
			if (currentLoss < generationBestLoss)
			{
				generationBestLoss = currentLoss;
				generationBestModelID = i;
			}

			this->timeTable_.ResetSimulation();
		}

		if (generationBestLoss < this->bestModel.Loss_ || this->bestModel.Loss_ < 0)
		// Update best loss of all the models (either current model is 
		// currently the best model or best model is not initialized.
		{
			this->bestModel.Loss_ = generationBestLoss;
			this->bestModel.AllChargingStations_.clear();
			for (auto&& station : population[generationBestModelID].AllChargingStations_)
			{
				this->bestModel.AllChargingStations_.push_back(
					StationParameters(station.StationID_,
						station.Position_,
						station.Capacity_,
						station.ClosestCityID_,
						station.ChargingWaitingTime_,
						station.EstimatedChargingLevel_
					));
			}

			// std::cout << this->bestModel.AllChargingStations_.size() << std::endl;

		}

		// Generate new population.
		population = this->geneticGetNewPopulation(
			population, this->optimizerParameters_.GeneticParameters_.NumBestSelection_);

		std::cout << "---------------------------------------" << std::endl;
		std::cout << "Best loss is: " << generationBestLoss << std::endl;
		std::cout << std::endl;
	}

	std::cout << "Algorithm finished" << std::endl;


	// Print all losses during the simulation.
	for (int32_t i = 0; i < this->optimizerParameters_.GeneticParameters_.NumGenerations_; i++)
	{
		std::cout << "Losses for the generation: " << i << std::endl;

		for (auto&& loss : lossesForGenerations[i])
		{
			std::cout << "Loss: " << loss << std::endl;
		}
	}

	std::cout << std::endl;
	std::cout << "Total best loss is: " << this->bestModel.Loss_ << std::endl;

	this->printBestLossInfo();
	
}


/// <summary>
/// Optimizes the position of the charging stations using the algorithm 
/// inspired by K-Means algorithm. Randomly generates centroids finds its clusters
/// (approximately using Dijsktra algorithm) and then approximately computes new
/// clusters (using geographical coordinates of the nodes).
/// </summary>
/// <param name="logs">Whether to write simulator logs.</param>
void Optimizer::KMeansAlgorithm(bool logs)
{
	// All current clusters.
	std::vector<std::vector<vertex_t>> allClusters;
	// Positions of all centroids.
	std::vector<MapPosition> allCentroidsPositions;

	// Currently the best model.
	// ModelRepresentation bestModel(-1, std::vector<StationParameters>());

	auto& map = this->timeTable_.GetTrafficSimulator().GetMap();
	for (int32_t i = 0; i < this->optimizerParameters_.KMeansParameters_.NumGenerations_; i++)
	// For each generation first simulate and then perform K-Means algorithm 
	// (new station positions is centroids from previous K-Means results).
	{
		if (i == 0)
		// First iteration -> randomly generate charging stations (centroids).
		{
			this->timeTable_.RandomlyGenerateChargingStations(
				this->simulationParameters_.NumStations, this->simulationParameters_);
			
			// Initialize centroids.
			auto& allChargingStations = map.GetChargingStations();
			for (auto&& station : allChargingStations)
			{
				allCentroidsPositions.push_back(MapPosition((*station).Position_));
			}
		}
		else
		// Positions of the charging stations (and begining centroids) 
		// are result centroids of the previous K-Means algorithm.
		{
			for (size_t j = 0; j < allCentroidsPositions.size(); j++)
			// For all centroids place a charging station.
			{
				// std::cout << allCentroidsPositions[j].CloserVertexID << std::endl;
				map.AddChargingStation(
					ChargingStation(j, allCentroidsPositions[j], 
						this->simulationParameters_.StationCapacity,
						map.GetGraph()[allCentroidsPositions[j].CloserVertexID].CityID_,
						this->simulationParameters_.ChargingWaitingTime, 
						this->simulationParameters_.MeanChargingLevel
					));
			}

			// Find nearest charging stations for each city.
			map.FindCitiesNearestChargingStations();
		}

		// Run simulation and compute the model loss.
		this->timeTable_.RunSimulation(
			this->simulationParameters_.NumStations, this->simulationParameters_, logs);

		double currentLoss = this->ModelLoss(
			this->optimizerParameters_.StationNumberParameter_,
			this->optimizerParameters_.RunDownParameter_,
			this->optimizerParameters_.DurationParameter_,
			this->optimizerParameters_.BatteryDifferenceParameter_,
			this->optimizerParameters_.WaitingTimesParameter_);
		

		if (currentLoss < this->bestModel.Loss_ || this->bestModel.Loss_ < 0)
		// Currently the best loss or first model -> update best model
		{
			this->bestModel.Loss_ = currentLoss;
			this->bestModel.AllChargingStations_.clear();
			for (auto&& station : map.GetChargingStations())
			{
				this->bestModel.AllChargingStations_.push_back(
					StationParameters((*station).StationID_,
						(*station).Position_,
						(*station).Capacity_,
						(*station).CityID_,
						(*station).GetChargingTime(0),
						(*station).GetExpectedWaitingTime()
					));
			}

			// std::cout << this->bestModel.AllChargingStations_.size() << std::endl;

		}

		std::cout << "Current loss is: " << currentLoss << std::endl;

		// Perform K-Means algorithm to find new positions of the stations.

		// Compute weights of the vertices based on the average battery 
		// levels after the simulation run.
		std::map<vertex_t, double> verticesWeights;
		this->kMeansComputeVerticesWeights(verticesWeights);

		this->timeTable_.ResetSimulation();

		for (int32_t j = 0; j < this->optimizerParameters_.KMeansParameters_.NumIterationsOneRun_; j++)
		// Run given number of iteration of the K-Means algorithm (find clusters and new centroids).
		{
			allClusters.clear();
			this->kMeansFindClusters(allClusters, allCentroidsPositions);
			allCentroidsPositions.clear();
			this->kMeansFindNewCentroids(
				allClusters, allCentroidsPositions, verticesWeights);
		}
	}

	std::cout << "-------------------------------------------------" << std::endl;
	std::cout << "Best loss is: " << this->bestModel.Loss_ << std::endl;

	this->printBestLossInfo();
}



/// <summary>
/// Initilizes the random generators.
/// </summary>
void Optimizer::initGenerators()
{
	this->gen = std::mt19937(rd());
}


/// <summary>
/// Prints all information about the best model.
/// </summary>
void Optimizer::printBestLossInfo()
{
	// Loss logs:
	std::cout << std::endl;
	std::cout << "Best model info:" << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "Number of stations: " << this->bestNumStations_ << std::endl;
	std::cout << "Number of run down batteries: " << this->bestCarsBatteryRunDown_ << std::endl;
	std::cout << "Average traveling duration: " << this->bestAverageTravelDuration_ << std::endl;
	std::cout << "Average battery difference between start and end: " << this->bestAverageBatteryDifference_ << std::endl;
	std::cout << "Average waiting times in charging station: " << this->bestAverageChargingWaitingTimes_ << std::endl;
	std::cout << "Number of vehicles which reached the end: " << this->bestNumReturned_ << std::endl;
	std::cout << "Number of used stations: " << this->bestUsedVehicles_ << std::endl;
	std::cout << "Number of vehicles which didn't finish before end of simulation: " << this->bestNotFinished_ << std::endl;
	std::cout << "Total number of used charging stations: " << this->bestNonNullStations_ << std::endl;
	// std::cout << "Total number of charging stations usage: " << this->bestTotalStationCustomers_ << std::endl;
	std::cout << "Total number of vehicles which didn't reach the end: " << this->bestCarsBatteryRunDown_ + this->bestNotFinished_ << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "Model loss: " << this->bestLoss_ << std::endl;
	std::cout << std::endl;

}


/// <summary>
/// Adds given stations to the map.
/// </summary>
/// <param name="stationsToMaintain">Contatiner of stations 
/// parameters to be added to the map.</param>
/// <param name="updateSimulator>Whether to update whole simulator
/// (if `true` then prepare simulator for immediate usage).</param>
void Optimizer::addChoosenStations(
	std::vector<StationParameters>& stationsToMaintain, bool updateSimulator)
{
	auto& map = this->timeTable_.GetTrafficSimulator().GetMap();
	int32_t stationID = map.GetChargingStations().size();

	for (auto&& station : stationsToMaintain)
	// Create charging stations from the given representation.
	{
		map.AddChargingStation(ChargingStation(stationID,
			station.Position_,
			station.Capacity_,
			station.ClosestCityID_,
			station.ChargingWaitingTime_,
			station.EstimatedChargingLevel_
		));

		stationID++;
	}

	if (updateSimulator)
	// Update simulator information based on the new stations 
	// (new number of stations and find closest stations to each city).
	{
		this->numStations_ = map.GetChargingStations().size();

		// Find nearest charging stations for each city.
		map.FindCitiesNearestChargingStations();
	}
}


/// <summary>
/// Does one decision from bernouli distribution.
/// </summary>
/// <param name="changeProbability">Probability of the choosen variable 
/// (if `true` then take it).</param>
/// <returns>Returns `true` if random variable with given probability was choosen,
///  else `false`. </returns>
bool Optimizer::randomlyChooseNumberOfStationsChange(double changeProbability)
{
	double choice = (double)rand() / RAND_MAX;

	// Choose variable with `changeProbability`.
	if (choice < changeProbability)
		return true;

	// Choose the other variable.
	return false;
}


/// <summary>
/// Randomly adds remaining charging stations when using the greedy optimization 
/// (if already not enough charging stations).
/// </summary>
/// <param name="updateSimulator>Whether to update whole simulator
/// (if `true` then prepare simulator for immediate usage).</param>
void Optimizer::greedyFillRestStations(bool updateSimulator)
{
	auto& map = this->timeTable_.GetTrafficSimulator().GetMap();
	int32_t numStations = map.GetChargingStations().size();

	// Number of stations to just not generate.
	int32_t numThrowAway = this->optimizerParameters_.GreedyParameters_.NumThrowAway_;

	// Number of the stations to throw away is greather or 
	// equal to total number of stations -> don't throw away any other station
	if (numThrowAway >= this->numStations_)
		numThrowAway = 0;


	// Rest of the stations to be generated.
	int32_t numStationsToGenerate = this->numStations_ - numStations - numThrowAway;
	if (numStationsToGenerate > 0)
	// Not enough stations generated using the greedy algorithm 
	//	-> generate randomly the rest.
	{
		this->timeTable_.RandomlyGenerateChargingStations(
			numStationsToGenerate, this->simulationParameters_);
	}


	if (updateSimulator)
	// Update simulator information based on the new stations 
	// (new number of stations and find closest stations to each city).
	{
		this->numStations_ = map.GetChargingStations().size();

		// Find nearest charging stations for each city.
		map.FindCitiesNearestChargingStations();
	}
}


/// <summary>
/// Takes used stations from the last simulation and randomly generates 
/// rest of the stations on randomly choosen positions where car battery ran down.
/// </summary>
/// <param name="stationsToMaintain">Container of the new stations representations.</param>
/// <param name="numThrowAway">Number of stations to not place on the new map 
/// (if not sufficient number of them used).</param>
void Optimizer::greedyChooseStations(
	std::vector<StationParameters>& stationsToMaintain, int32_t numThrowAway)
{
	this->greedyTakeUsedChargingStations(stationsToMaintain);
	this->greedyGenerateStationsBatteryRunDown(stationsToMaintain, numThrowAway);
}


/// <summary>
/// Finds all charging stations that have been used at least once and 
/// stores its parameters (for future usage).
/// </summary>
/// <param name="stationsToMaintain">Container to store info 
/// about used charging stations.</param>
void Optimizer::greedyTakeUsedChargingStations(
	std::vector<StationParameters>& stationsToMaintain)
{
	auto& chargingLevels = this->timeTable_.GetChargingLevels();
	TrafficSimulator& trafficSimulator = this->timeTable_.GetTrafficSimulator();
	auto& allChargingStations = trafficSimulator.GetMap().GetChargingStations();

	for (size_t i = 0; i < chargingLevels.size(); i++)
	// Saves all charging stations which have been used 
	//	-> rest throw away (and replace them with other station).
	{
		if (chargingLevels[i].size() > 0)
		// Used charging station -> save its properties
		{
			stationsToMaintain.push_back(
				StationParameters((*allChargingStations[i]).StationID_,
					(*allChargingStations[i]).Position_,
					(*allChargingStations[i]).Capacity_,
					(*allChargingStations[i]).CityID_,
					(*allChargingStations[i]).GetChargingTime(0),
					(*allChargingStations[i]).GetExpectedWaitingTime()
				));
		}
	}
}


/// <summary>
/// Generates new charging stations on some position of the 
/// previous battery run down.
/// </summary>
/// <param name="stationsToMaintain">Container to store info 
/// about used charging stations.</param>
/// <param name="numThrowAway">Number of stations to throw away 
/// (to optimize number of used stations).</param>
void Optimizer::greedyGenerateStationsBatteryRunDown(
	std::vector<StationParameters>& stationsToMaintain, int32_t numThrowAway)
{
	// Prepare cities probabilities (weighted by the number of 
	// run down positions in it).
	std::vector<double> citiesProbabilities;
	auto& runDownPositions = this->timeTable_.GetRunDownPositions();
	for (auto&& cityPositionsPair : runDownPositions)
	{
		citiesProbabilities.push_back(cityPositionsPair.second.size());
	}

	// Init distribution of the cities to generate charging station.
	DiscreteDistribution citiesDistribution(citiesProbabilities);

	// Number of stations to randomly generate
	int32_t newStationsTreshold = 
		this->numStations_ - stationsToMaintain.size() - numThrowAway;

	for (int32_t i = 0; i < newStationsTreshold; i++)
	// Generate new stations on some positions of the previous battery run down.
	{
		// Check if there are some run down positions to take.
		bool positionExists = false;
		for (auto&& runDownPair : runDownPositions)
		{
			if (runDownPair.second.size() > 0)
			// There exists some position to choose.
			{
				positionExists = true;
				break;
			}
		}

		// No candidate position for the next charging station 
		// -> don't add any other station.
		if (!positionExists)
			break;


		// Generate station's city:
		int32_t cityIndex = citiesDistribution(this->gen);
		auto it = runDownPositions.begin();
		std::advance(it, cityIndex);
		while ((*it).second.size() == 0)
		// Repeat city choice until choosen city with some charging station 
		// candidates (still not empty list of candidates). 
		{
			cityIndex = citiesDistribution(this->gen);
			it = runDownPositions.begin();
			std::advance(it, cityIndex);
		}

		// Uniformly randomly generate position of the station:
		std::uniform_int_distribution<> 
			positionDistribution(0, (*it).second.size() - 1);
		int32_t positionIndex = positionDistribution(this->gen);


		// Add new station to list of the new charging stations.
		stationsToMaintain.push_back(
			StationParameters(-1,
				(*it).second[positionIndex],
				this->simulationParameters_.StationCapacity,
				(*it).first,
				this->simulationParameters_.ChargingWaitingTime,
				this->simulationParameters_.MeanChargingLevel
			));

		// Remove choosen position from the list of candidates.
		auto removeIterator = (*it).second.begin() + positionIndex;
		(*it).second.erase(removeIterator);
	}
}


/// <summary>
/// Creates new population of the models using the genetic algorithm approach.
/// Individual is encodes as vector of charging stations.
/// Part of the selection is just to take given number of best models, rest is
/// choosen using the tournament selection. We use one poin crossover (index of station).
/// Mutation is either replacing station with new randomly generated 
/// (position mutation) or randomly add or delete station (number of stations mutation).
/// </summary>
/// <param name="population">Current population of parents.</param>
/// <param name="numBest">Number of best models to be immediately copy to the
/// next generation of models.</param>
/// <returns>Returns new generation of children (for next iteration).</returns>
std::vector<ModelRepresentation> Optimizer::geneticGetNewPopulation(
	std::vector<ModelRepresentation>& population, int32_t numBest)
{
	std::vector<ModelRepresentation> newPopulation;

	// Sort models by its loss.
	std::map<double, int32_t> lossIndexPair;
	for (size_t i = 0; i < population.size(); i++)
	{
		lossIndexPair.insert(std::make_pair(population[i].Loss_, i));
	}


	// Take the `numBest` best models based on the loss and 
	// copy them into the new population.
	int32_t numBestToTake = numBest;
	if (numBest > this->optimizerParameters_.GeneticParameters_.PopulationSize_)
	// Check if `numBest` is not higher than population size 
	// -> if so, change it to the population size.
	{
		numBestToTake = 
			this->optimizerParameters_.GeneticParameters_.PopulationSize_;
	}


	for (int32_t i = 0; i < numBestToTake; i++)
	// Copy `numBestToTake` best models to next generation.
	{
		auto it = lossIndexPair.begin();
		std::advance(it, i);
		newPopulation.push_back(population[(*it).second]);
	}


	// Compute number of remaining models to generate using genetic algorithm.
	int32_t numToTake =
		this->optimizerParameters_.GeneticParameters_.PopulationSize_ - numBestToTake;

	for (int32_t i = 0; i < numToTake / 2; i++)
	// In each step generate 2 children using genetic algoritm techniques.
	{
		// Selection:
		int32_t firstParentID =
			this->geneticTournamentSelection(population, 
				this->optimizerParameters_.GeneticParameters_.TournamentSelectionTreshold_);
		int32_t secondParentID =
			this->geneticTournamentSelection(population, 
				this->optimizerParameters_.GeneticParameters_.TournamentSelectionTreshold_);
	

		// Crossover:
		auto newPairChildren = 
			this->geneticPointCrossover(population, firstParentID, secondParentID);
		ModelRepresentation firstChild = std::move(newPairChildren.first);
		ModelRepresentation secondChild = std::move(newPairChildren.second);
		

		// Mutation:
		this->geneticMutation(firstChild, 
			this->optimizerParameters_.GeneticParameters_.MutationTreshold_);
		this->geneticMutation(secondChild, 
			this->optimizerParameters_.GeneticParameters_.MutationTreshold_);


		// Adding children to new population:
		newPopulation.push_back(std::move(firstChild));
		newPopulation.push_back(std::move(secondChild));
	}

	if (this->optimizerParameters_.GeneticParameters_.PopulationSize_ % 2 == 1)
	// If population has odd number of members 
	// -> fill last member with the first member from the last population
	{
		newPopulation.push_back(std::move(population[0]));
	}

	return std::move(newPopulation);
}


/// <summary>
/// Randomly generate future number of stations (for number of stations mutation).
/// Using normal distribution to generate the count.
/// </summary>
/// <param name="mean">Mean of the normal distribution to generate from.</param>
/// <param name="variance">Variance of the normal distribution to generate from.</param>
/// <returns>Retunrs randomly generated new number of stations.</returns>
int32_t Optimizer::geneticRandomlyChooseNumberOfStationsChange(
	int32_t mean, double variance)
{
	// Initialize the distribution:
	std::normal_distribution<double> distribution(mean, variance);
	double number = distribution(this->gen);
	

	// Generate the count:
	int32_t choice = (int32_t)std::round(number);
	
	// If choice is non-positive number, choose 1.
	if (choice < 1)
		choice = 1;

	return choice;
}


/// <summary>
/// Randomly initializes first population of the genetic algorithm.
/// </summary>
/// <param name="population">Container to store generated 
/// members of the population.</param>
void Optimizer::geneticRandomlyInitPopulation(
	std::vector<ModelRepresentation>& population)
{
	for (int32_t i = 0; i < this->optimizerParameters_.GeneticParameters_.PopulationSize_; i++)
	// Generate all members of the population.
	{
		std::vector<StationParameters> allStations;

		// Get size of the member.
		int32_t numberOfStations = this->geneticRandomlyChooseNumberOfStationsChange(
			this->simulationParameters_.NumStations, 
			this->optimizerParameters_.GeneticParameters_.MemberSizeVariance_);
		
		for (int32_t j = 0; j < numberOfStations; j++)
		// Generate all charging stations of the model.
		{
			ChargingStation station = 
				this->timeTable_.GetTrafficSimulator()
				.GenerateChargingStation(j, this->simulationParameters_);

			allStations.push_back(
				StationParameters(station.StationID_,
					station.Position_,
					station.Capacity_,
					station.CityID_,
					station.GetChargingTime(0),
					station.GetExpectedWaitingTime()
				));
		}
		population.push_back(ModelRepresentation(0, std::move(allStations)));
	}
}


/// <summary>
/// Select one parent using tournament selection (compare only 2 uniformly randomly 
/// choosen members of the population and randomly choose parent from them).
/// Higher probability has the one with lower loss. 
/// </summary>
/// <param name="population">Container of all members of the current population.</param>
/// <param name="betterTreshold">Probability of choosing the better
/// member (from 0 to 1).</param>
/// <returns>Returns the member ID using the tournament selection.</returns>
int32_t Optimizer::geneticTournamentSelection(
	std::vector<ModelRepresentation>& population, double betterTreshold)
{
	if (population.size() == 1)
	// Just one member of the population 
	// (doesn't make much sense -> probably error).
	{
		return 0;
	}

	// Generate randomly 2 members of the population.
	std::uniform_int_distribution<int32_t> distr(0, population.size() - 1);
	int32_t first = distr(this->gen);
	int32_t second = distr(this->gen);
	while (second == first)
	// Repeat until second parent is different from the first one.
	{
		second = distr(this->gen);
	}


	if (population[second].Loss_ < population[first].Loss_)
	// Make sure that in `first` is a member with lower loss, than in `second`.
	{
		std::swap(first, second);
	}


	// Choose one member (higher probability has the one with the lower loss).
	double tournamentChoice = (double)rand() / RAND_MAX;

	if (tournamentChoice < 
		this->optimizerParameters_.GeneticParameters_.TournamentSelectionTreshold_)
	// Choose the more probable member (with lower loss).
	{
		return first;
	}

	// Choose the less probable member (with higher loss).
	return second;
}


/// <summary>
/// Performs one point crossover on the parents 
/// (after crossover point children switches parameters (charging stations)).
/// </summary>
/// <param name="population">Container of all member from the current population.</param>
/// <param name="firstParentID">ID of the first parent for crossover.</param>
/// <param name="secondParentID">ID of the second parent for crossover.</param>
/// <returns>Return pair of two new children for the future population (before mutation).</returns>
std::pair<ModelRepresentation, ModelRepresentation> Optimizer::geneticPointCrossover(
	std::vector<ModelRepresentation>& population, 
	int32_t firstParentID, int32_t secondParentID)
{
	// Find parents.
	ModelRepresentation& firstParent = population[firstParentID];
	ModelRepresentation& secondParent = population[secondParentID];

	if (firstParent.AllChargingStations_.size() > 
		secondParent.AllChargingStations_.size())
	// Make sure first parent is smaller 
	// (has at most same number of charging stations than second.
	{
		std::swap(firstParent, secondParent);
	}

	// Randomly generate croossover point.
	std::uniform_int_distribution<int32_t> 
		distr(0, firstParent.AllChargingStations_.size() - 1);
	int32_t crossoverPoint = distr(this->gen);


	// Do a one point crossover.
	ModelRepresentation firstChild, secondChild;
	for (size_t i = 0; i < secondParent.AllChargingStations_.size(); i++)
	// Before crossover point copy stations, 
	// after it switch parents and do the same.
	{
		if (i < crossoverPoint)
		// Before crossover point -> let parameters same.
		{
			firstChild.AllChargingStations_.push_back(
				std::move(firstParent.AllChargingStations_[i]));
			secondChild.AllChargingStations_.push_back(
				std::move(secondParent.AllChargingStations_[i]));
		}
		else
		// After crossover point 
		// -> change parameters (first has from second, second has from first).
		{
			// We know that second parent is always greather or equal than first 
			// -> we can just copy station (without check)
			firstChild.AllChargingStations_.push_back(
				std::move(secondParent.AllChargingStations_[i]));
			
			if (i < firstParent.AllChargingStations_.size())
			// If there already exists stations of the smaller parent
			// -> add them to the child.
			{
				secondChild.AllChargingStations_.push_back(
					std::move(firstParent.AllChargingStations_[i]));
			}
		}
	}

	return std::make_pair(std::move(firstChild), std::move(secondChild));
}


/// <summary>
/// Performs mutation on one model. 
/// Randomly changes some stations and increases or decreases number of stations.
/// </summary>
/// <param name="model">Model to be changed.</param>
/// <param name="mutationProbability">Probability of the mutation for each station 
/// (probability that the station will be replaced by the random one).</param>
void Optimizer::geneticMutation(ModelRepresentation& model, double mutationProbability)
{	
	// Perform mutation on each station.
	this->geneticMutationRandomlyChangeStations(model, mutationProbability);

	// Perform mutation of number of stations.
	this->geneticRandomlyChooseNumberOfStationsChange(
		model.AllChargingStations_.size(), 
		this->optimizerParameters_.GeneticParameters_.MemberSizeVariance_);
}


/// <summary>
/// Randomly replaces stations from the given model with the 
/// new randomly generated station.
/// </summary>
/// <param name="model">Model to be changed.</param>
/// <param name="mutationProbability">Probability of the mutation for each station
/// (probability that the station will be replaced by the random one).</param>
void Optimizer::geneticMutationRandomlyChangeStations(
	ModelRepresentation& model, double mutationProbability)
{
	for (size_t i = 0; i < model.AllChargingStations_.size(); i++)
	// Try mutation on each station.
	{
		double mutate = (double)rand() / RAND_MAX;

		if (mutate < mutationProbability)
		// Mutation of the station is active -> replace it with the random one.
		{
			ChargingStation newStation =
				this->timeTable_.GetTrafficSimulator()
				.GenerateChargingStation(-1, this->simulationParameters_);

			model.AllChargingStations_[i] = StationParameters(-1,
				newStation.Position_,
				newStation.Capacity_,
				newStation.CityID_,
				newStation.GetChargingTime(0),
				newStation.GetExpectedWaitingTime());
		}
	}
}

/// <summary>
/// Randomly adds or deletes stations from the given model using normal 
/// distribution to determine new number of stations.
/// </summary>
/// <param name="model">Model to be changed.</param>
/// <param name="numberStationsMean">Mean of the new number of stations of 
/// the given model.</param>
/// <param name="numberStationsVariance">Variance of the new number of 
/// stations of the given model.</param>
void Optimizer::geneticMutationChangeNumberOfStations(ModelRepresentation& model,
	int32_t numberStationsMean, double numberStationsVariance)
{
	// Generated new number of stations.
	uint32_t newNumberOfStations =
		(uint32_t)this->geneticRandomlyChooseNumberOfStationsChange(
			numberStationsMean, numberStationsVariance);

	if (newNumberOfStations > model.AllChargingStations_.size())
	// New number of stations is higher -> randomly generate new stations
	{
		for (size_t i = 0; i < newNumberOfStations - model.AllChargingStations_.size(); i++)
		// Radnomly generate the overflowing stations.
		{
			ChargingStation newStation =
				this->timeTable_.GetTrafficSimulator()
				.GenerateChargingStation(-1, this->simulationParameters_);

			model.AllChargingStations_.push_back(
				StationParameters(-1,
					newStation.Position_,
					newStation.Capacity_,
					newStation.CityID_,
					newStation.GetChargingTime(0),
					newStation.GetExpectedWaitingTime()));
		}
	}

	while (newNumberOfStations < model.AllChargingStations_.size())
	// New number of station is lower -> randomly delete stations
	{
		// Randomly delete station.
		std::uniform_int_distribution<int32_t> 
			distr(0, model.AllChargingStations_.size() - 1);
		int32_t deletionPoint = distr(this->gen);

		model.AllChargingStations_.erase(
			model.AllChargingStations_.begin() + deletionPoint);
	}
}



/// <summary>
/// Finds cluster of vertices for the current centroids. 
/// The distances from the clusters are computed only approximately using 
/// distances of other vertices from the closest vertex from the centroid 
/// and its distance from the centroid.
/// Because of that there could be some vertices situated in the wrong cluster 
/// (some centroid is closer). However we assume the differences of the distances 
/// are neglieble and doesn't happen very often.
/// </summary>
/// <param name="allClusters">Container to store all clusters (each item is a 
/// container representing one cluster and containg IDs of the vertices 
/// which belong there).</param>
/// <param name="allCentroidsPositions">Container of all current centroids positions.</param>
void Optimizer::kMeansFindClusters(std::vector<std::vector<vertex_t>>& allClusters,
	std::vector<MapPosition>& allCentroidsPositions)
{
	// For each vertex we remember closest centroid (`int32_t` 
	// signifies centroid ID) and its distance from the closest centroid.
	std::unordered_map<vertex_t, std::pair<int32_t, double>> closestCentroidsInfo;

	for (size_t centroidID = 0; centroidID < allCentroidsPositions.size(); centroidID++)
	// For all centroids find distance to all the vertices.
	{
		// Find distance of each vertex from the cetroid using Dijkstra algorithm.
		std::vector<double> distances = 
			this->timeTable_.GetTrafficSimulator().GetMap()
				.DijkstraFindDistances(
					allCentroidsPositions[centroidID].CloserVertexID);

		for (size_t vertexID = 0; vertexID < distances.size(); vertexID++)
		// For each vertex check whether the distance from the centroid 
		// is currently the smallest, if so update distances.
		{
			// Compute current distance as distance between current vertex and 
			// the closest vertex from the centroid plus the distance to the 
			// vertex from the centroid.
			double currentDistance = distances[vertexID] + 
				allCentroidsPositions[centroidID].GetDistanceFromCloserVertex(
					this->timeTable_.GetTrafficSimulator().GetMap().GetGraph(),
					this->simulationParameters_.SegmentLength);

			if (centroidID == 0)
			// First centroid -> initialize all vertex and centroid distances 
			// (currently best are from this centroid).
			{
				closestCentroidsInfo.insert(std::make_pair(vertexID, 
					std::make_pair(centroidID, currentDistance)));

				continue;
			}

			if (currentDistance < closestCentroidsInfo[vertexID].second)
			// Current centroid is closer to the vertex -> update distances
			{
				closestCentroidsInfo[vertexID].first = centroidID;
				closestCentroidsInfo[vertexID].second = currentDistance;
			}
		}
	}

	// Create the clusters:
	allClusters = std::vector<std::vector<vertex_t>>
		(allCentroidsPositions.size(), std::vector<vertex_t>());
	for (auto&& vertexPair : closestCentroidsInfo)
	// Put each vertex in the correct cluster.
	{
		allClusters[vertexPair.second.first].push_back(vertexPair.first);
	}
}


/// <summary>
/// Computes weights for each vertex for centroids searching (based on the 
/// average battery levels on the edges connected with the vertex). 
/// </summary>
/// <param name="verticesWeights">Container of weights for each vertex.</param>
void Optimizer::kMeansComputeVerticesWeights(
	std::map<vertex_t, double>& verticesWeights)
{
	auto& graph = this->timeTable_.GetTrafficSimulator().GetMap().GetGraph();
	verticesWeights.clear();

	auto verticesIterators = boost::vertices(graph);
	for (auto vert : boost::make_iterator_range(verticesIterators))
	// For each vertex compute weight.
	{
		verticesWeights.insert(std::make_pair(vert, 1));

		auto neighborsRange = boost::adjacent_vertices(vert, graph);
		for (auto neighbor : boost::make_iterator_range(neighborsRange))
		// For each edge appropriately increase vertex weight based on the average 
		// battery level on each segment and its distance from the vertex.
		{
			Edge& edge = graph[boost::edge(vert, neighbor, graph).first];

			// Part of the vertex weight belonging to current edge.
			double edgeWeight = 0;

			// How much the significance of the segment value decreases 
			// when moving one segment away of the vertex.
			double discountParameter = 1 / edge.GetNumSegments();

			auto& segmentsBatteryLevels = edge.GetBatteryLevels();
			
			for (int32_t i = 0; i < edge.GetNumSegments(); i++)
			// Go through all segments and appropriately add its part of the weight.
			{
				// Average battery level of the segment.
				double segmentAverageLevel = 
					segmentsBatteryLevels[i].second / segmentsBatteryLevels[i].first;

				if (segmentAverageLevel == 0)
					segmentAverageLevel = 1;


				// Segment significance, if further from the vertex, 
				// then its significance decreases.
				double segmentWeight = (1 - i * discountParameter);
				if (neighbor < vert)
				// Order of the segments from the `vert` is in decreasing order 
				// (significance with rising `i` rises).
				{ 
					segmentWeight = i * discountParameter;
				}

				// segmentWeight *= 100000;
				// std::cout << segmentWeight << "    " << (1 - segmentAverageLevel) * segmentWeight << std::endl;
				// If lower average battery level, then greather part of the 
				// weight should be added to the edge weight.
				edgeWeight += (1 - segmentAverageLevel);// * segmentWeight;
			}
			// std::cout << edgeWeight << std::endl;

			// Increase vertex weight by the weight of the edge.
			verticesWeights[vert] += edgeWeight;
		}
	}
}


/// <summary>
/// Finds new centroids based on the clusters.
/// </summary>
/// <param name="allClusters">Vector of all current clusters.</param>
/// <param name="allCentroidsPositions">Positions of all centroids.</param>
/// <param name="verticesWeights">Weight of each vertex.</param>
void Optimizer::kMeansFindNewCentroids(std::vector<std::vector<vertex_t>>& allClusters, 
	std::vector<MapPosition>& allCentroidsPositions, std::map<vertex_t, double>& verticesWeights)
{
	allCentroidsPositions.clear();
	for (auto&& cluster : allClusters)
	// For each cluster find the new centroid.
	{
		allCentroidsPositions.push_back(
			this->kMeansFindNewCentroid(cluster, verticesWeights));

		// std::cout << cluster.size() << std::endl;
	}
}


/// <summary>
/// Finds the centroid of the cluster. Just computes weighted average of the
/// coordinates of the vertices from the cluster. Weights are derived from the
/// average battery levels on the edges connected with the vertex. For the new
/// position (average) approximately finds the closest position on the map.
/// </summary>
/// <param name="cluster">Cluster to find centroid for.</param>
/// <param name="verticesWeights">Weights of the vertices.</param>
/// <returns>Returns position of the new centroid of the given cluster.</returns>
MapPosition Optimizer::kMeansFindNewCentroid(std::vector<vertex_t>& cluster, 
	std::map<vertex_t, double>& verticesWeights)
{
	double latitudeSum = 0;
	double longitudeSum = 0;

	auto& graph = this->timeTable_.GetTrafficSimulator().GetMap().GetGraph();

	// Sum of all weights (to compute weighted average).
	double normalizationParameter = 0;

	for (auto&& vertex : cluster)
	// Compute weighted average of coordinates for each vertex from the cluster.
	{
		double vertexWeight = verticesWeights[vertex];
		// std::cout << vertexWeight << std::endl;
		latitudeSum += graph[vertex].Location_.Latitude * vertexWeight;
		longitudeSum += graph[vertex].Location_.Longitude * vertexWeight;
		normalizationParameter += vertexWeight;
	}

	// Get centroid coordinates.
	double centroidLatitude = latitudeSum / normalizationParameter;
	double centroidLongitude = longitudeSum / normalizationParameter;

	// std::cout << centroidLatitude << " " << centroidLongitude << std::endl;

	// Find approximately the nearest position of the map from the centroid.
	return this->timeTable_.GetTrafficSimulator().GetMap()
		.FindClosestPositionToCoordinates(centroidLatitude, centroidLongitude);
}
