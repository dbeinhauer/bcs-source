#include "optimizer/Car.hpp"


/// <summary>
/// Initializes car object.
/// </summary>
/// <param name="startTime">Time of the car departure.</param>
/// <param name="waitingTime">Waiting time for start returning (after reaching end position).</param>
/// <param name="start">Start position.</param>
/// <param name="end">End position.</param>
/// <param name="consumption">Car consumption (battery percentages per minute).</param>
/// <param name="batteryLevel">Starting battery level.</param>
/// <param name="relativeSpeed">Speed of the car (in kilometers per minute (1/60 * km/hr)).</param>
Car::Car(double startTime, double waitingTime, MapPosition start,
	MapPosition end, double consumption, double batteryLevel, double relativeSpeed)
	: Vehicle(std::move(startTime), std::move(start), std::move(end), 
		std::move(consumption), std::move(batteryLevel))
{
	this->relativeVelocity_ = relativeSpeed;
	this->waitingTime_ = std::move(waitingTime);
	this->returning_ = false;
}


/// <summary>
/// Appropriately changes car setup to returning to start.
/// Changes start and end position and sets flags of the car to returning.
/// </summary>
void Car::StartReturning(double actualTime)
{
	// Appropriately sets flags.
	this->startTime_ = actualTime;
	this->startBatteryLevel_ = this->batteryLevel_;
	this->returning_ = true;
	this->insideEdge_ = true;
	this->chargingTested_ = false;

	// Swaps start and end positions.
	MapPosition helpPosition = this->startPosition_;
	this->startPosition_ = this->endPosition_;
	this->endPosition_ = helpPosition;
	this->currentPosition_ = this->startPosition_;

	this->batteryLevel_ = 1;
}


/// <summary>
/// </summary>
/// <returns>Returns `true` if car is returning (already has been in the original finish).</returns>
bool Car::IsReturning()
{
	return this->returning_;
}


/// <summary>
/// </summary>
/// <returns>Returns waiting time of the car in the end to start returnning.</returns>
double Car::GetWaitingTime()
{
	return this->waitingTime_;
}
