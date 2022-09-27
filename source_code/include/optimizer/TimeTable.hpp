#ifndef TIMETABLE_H_
#define TIMETABLE_H_


#include <queue>

#include "TableEvent.hpp"
#include "TrafficSimulator.hpp"
#include "SimulationParameters.hpp"


// Structure for priority queue (to compare event times).
struct  MyGreater 
{
	bool operator()(const TableEvent& lhs, const TableEvent& rhs)
	{
		return lhs.ActionTime_ > rhs.ActionTime_;
	}
};


/// <summary>
/// Class to manage all high level simulation logic.
/// </summary>
class TimeTable
{
public:
	// Current simulation time.
	double ActualTime_;
	// Nummber of cars which ran down with the battery level.
	int32_t CarsBatteryRunDown_ = 0;

	int32_t RunDownCharging = 0;
	int32_t NumCharging_ = 0;
	int32_t NumNotFinishedCars = 0;


	TimeTable(SimulationParameters& simulationParameters);

	bool LoadMap(SimulationParameters& simulationParameters);
	
	void RandomlyGenerateChargingStations(
		int32_t numStations, SimulationParameters& simulationParameters);

	// Managing the simulation:
	void ResetSimulation();

	void RunSimulation(int32_t numStations,
		SimulationParameters& simulationParameters, bool logs);

	void AddEvent(TableEvent tableEvent);
	const TableEvent& GetNextEvent();
	void RunRestSimulation(SimulationParameters& simulationParameters, bool logs);


	// Obtaining simulation information:
	std::vector<double>& GetTravelDurations();
	std::vector<double>& GetBatteryDifferences();
	std::vector<std::vector<double>>& GetChargingLevels();
	std::vector<std::pair<int32_t, double>>& GetWaitingTimes();
	std::map<edge_t, std::vector<BateryPair>&> GetAllMapSegmentsInfo();
	TrafficSimulator& GetTrafficSimulator();
	std::map<int32_t, std::vector<MapPosition>>& GetRunDownPositions();
	int32_t GetNumReturnedVehicles();

private:
	// Object to manage simulation operations.
	TrafficSimulator trafficSimulator_;
	// Queue of all scheduled events (for discrete simulation).
	std::priority_queue<TableEvent, std::vector<TableEvent>, MyGreater> allEvents_;


	// CARS info:

	// Container of the traveling durations of all paths 	
	// (separately for first part of the path and returning).
	std::vector<double> travelDurations_;

	// Contatiner of all the differences between start and end battery levels 
	// (separately for first part of the path and returning).
	// Negative values means that vehicle reaches end with higher battery level
	// than when it started, positive values means battery level decreased.
	std::vector<double> batteryDifferences_;


	// CHARGING STATIONS info:

	// Container of container of charged battery levels (for each charging
	// station 1 container) (its ID is index in the container).
	std::vector<std::vector<double>> chargingLevels_;
	
	// Container of pairs of number of cars and sum of waiting times in the 
	// line of all cars in the given charging station (each charging station
	// has 1 pair). ID of the station is index in the container.
	std::vector<std::pair<int32_t, double>> waitingTimes_;

	// Battery run down positions of each city.
	std::map<int32_t, std::vector<MapPosition>> runDownPositions_;
	int32_t numReturnedVehicles_ = 0;


	// Discrete simulator actions:
	void reaction_START(SimulationParameters& simulationParameters);
	void reaction_RUNNIG(SimulationParameters& simulationParameters, bool logs);
	void reaction_CHARGING();
	void reaction_END_CHARGING(SimulationParameters& simulationParameters);
	void reaction_WAITING_LINE();
	void reaction_WAITING_RETURN(bool logs);

	void addRunDownCar(Car& car);
};

#endif