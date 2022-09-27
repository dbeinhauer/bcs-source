#include "optimizer/SimulationParameters.hpp"


/// <summary>
/// Initializes simulation parameters.
/// </summary>
SimulationParameters::SimulationParameters()
{
	this->SavePrepared = false;

	this->SimulationTime = 0;
	this->NumClosestStations = 0;
	this->NumStations = 0;
	this->StationCapacity = 0;

	this->SegmentLength = 0;
	this->CarConsumption = 0;
	this->ExponentialLambdaCities = 0;
	this->ExponentialLambdaDepartures = 0;
	this->EndCityRatio = 0;
	this->BatteryTresholdLambda = 0;
	this->CarBatteryMean = 0;
	this->CarBatteryDeviation = 0;
	this->CarStartBatteryBottomLimit = 0;
	this->ChargingTreshold = 0;
	this->NotChargingTreshold = 0;
	this->BatteryTolerance = 0;
	this->CarVelocity = 0;
	this->ChargingWaitingTime = 0;
	this->MeanChargingLevel = 0;
}