#ifndef CHARGINGSTATION_H_
#define CHARGINGSTATION_H_

#include <queue>

#include "MapPosition.hpp"

/// <summary>
/// Class representing one charging station.
/// </summary>
class ChargingStation
{
public:
	// ID is the same as index in vector in `TrafficSimulator` class.
	int32_t StationID_;
	// Position of the station.
	MapPosition Position_;
	// ID of the nearest city.
	int32_t CityID_;
	// Capacity of the station.
	int32_t Capacity_;
	// Current number of customers at the station.
	int32_t NumCustomers_;

	ChargingStation();
	ChargingStation(int32_t stationID, MapPosition position, int32_t capacity,
		int32_t closestCityID, double chargingWaitingTime, double estimatedChargingLevel);

	void AddCustomer(int32_t carID);
	void RemoveCustomer();

	int NextCustomer();
	
	double GetChargingTime(double batteryLevel);
	double GetExpectedWaitingTime();

private:
	// Queue of the waiting cars.
	std::queue<int> waitingCars_;
	// Sum of all waiting times.
	double chargingWaitingTime_;
	// Estimated level of battery to be charged.
	double estimatedChargingLevel_;
	// Total number of customer who visited the station.
	int64_t totalCustomers_;
};

#endif // !CHARGINGSTATION_H_