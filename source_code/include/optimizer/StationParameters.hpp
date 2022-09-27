#ifndef	STATIONPARAMETERS_H_
#define STATIONPARAMETERS_H_


#include "MapPosition.hpp"


/// <summary>
/// Struct to store information about the charging station.
/// </summary>
class StationParameters
{
public:
	// ID of the station.
	int32_t StationID_;
	// Position of the station.
	MapPosition Position_;
	// Capacity of the station.
	int32_t Capacity_;
	// ID of the nearest city of the station.
	int32_t ClosestCityID_;
	// Waiting time for full charge.
	double ChargingWaitingTime_;
	// Expected level to be charged on the station.
	double EstimatedChargingLevel_;


	StationParameters();
	StationParameters(int32_t stationID, MapPosition position, 
		int32_t capacity, int32_t closestCityID, 
		double chargingWaitingTime, double estimatedChargingLevel);
};


#endif