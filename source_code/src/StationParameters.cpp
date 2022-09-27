#include "optimizer/StationParameters.hpp"


/// <summary>
/// Initializes station parameters.
/// </summary>
StationParameters::StationParameters()
{
	this->StationID_ = -1;
	this->Position_ = MapPosition();
	this->Capacity_ = -1;
	this->ClosestCityID_ = -1;
	this->ChargingWaitingTime_ = 0;
	this->EstimatedChargingLevel_ = 0;
}


/// <summary>
/// Initializes station parameters.
/// </summary>
/// <param name="stationID">ID of the station.</param>
/// <param name="position">Position of the station.</param>
/// <param name="capacity">Capacity of the station.</param>
/// <param name="closestCityID">ID of the nearest city of the station.</param>
/// <param name="chargingWaitingTime">Waiting time for full charge.</param>
/// <param name="estimatedChargingLevel">Expected level to be 
/// charged on the station.</param>
StationParameters::StationParameters(int32_t stationID, MapPosition position,
	int32_t capacity, int32_t closestCityID, double chargingWaitingTime, 
	double estimatedChargingLevel)
{
	this->StationID_ = stationID;
	this->Position_ = position;
	this->Capacity_ = capacity;
	this->ClosestCityID_ = closestCityID;
	this->ChargingWaitingTime_ = chargingWaitingTime;
	this->EstimatedChargingLevel_ = estimatedChargingLevel;
}