#include "optimizer/ChargingStation.hpp"

/// <summary>
/// Initializes charging station.
/// </summary>
ChargingStation::ChargingStation()
{
	this->StationID_ = -1;
	this->Position_ = MapPosition();
	this->CityID_ = -1;
	this->Capacity_ = 0;
	this->NumCustomers_ = 0;

	this->chargingWaitingTime_ = 0;
	this->estimatedChargingLevel_ = 0;
	this->totalCustomers_ = 0;
}


/// <summary>
/// Initializes charging station.
/// </summary>
/// <param name="stationID">ID of the station.</param>
/// <param name="position">Station position.</param>
/// <param name="capacity">Station capacity.</param>
/// <param name="closestCityID">Closest city center from the station.</param>
/// <param name="chargingWaitingTime">Waiting time of the complete charge of the battery.</param>
/// <param name="estimatedChargingLevel">Estimated charing level of the customers.</param>
ChargingStation::ChargingStation(int32_t stationID, MapPosition position,
	int32_t capacity, int32_t closestCityID, double chargingWaitingTime,
	double estimatedChargingLevel)
	: StationID_(std::move(stationID)), Position_(std::move(position)),
	Capacity_(std::move(capacity)), chargingWaitingTime_(std::move(chargingWaitingTime)),
	estimatedChargingLevel_(std::move(estimatedChargingLevel))
{
	this->CityID_ = closestCityID;
	this->NumCustomers_ = 0;
	this->totalCustomers_ = 0;
}


/// <summary>
/// Adds customer into the waiting line.
/// </summary>
/// <param name="carID">ID of the customers car.</param>
void ChargingStation::AddCustomer(int32_t carID)
{
	this->NumCustomers_++;
	this->waitingCars_.push(carID);
	this->totalCustomers_++;
}


/// <summary>
/// Updates information about customers when some customer left the charging station.
/// </summary>
void ChargingStation::RemoveCustomer()
{
	this->NumCustomers_--;

	// Validity check (if no customer in the station -> no change).
	if (this->NumCustomers_ < 0)
		this->NumCustomers_ = 0;
}



/// <summary>
/// Returns ID of the next customer waiting in line for charging and 
/// removes it from the waiting queue.
/// </summary>
/// <returns>Returns ID of the next waiting customer or 
/// `-1` if no next customer available.</returns>
int ChargingStation::NextCustomer()
{
	// Check whether no waiting cars -> return `-1`
	if (this->waitingCars_.empty())
		return -1;

	// Get next customers ID and removes it from the queue.
	int nextCustomerID = this->waitingCars_.front();
	this->waitingCars_.pop();

	return nextCustomerID;
}


/// <summary>
/// Computes waiting time needed for the vehicle to be fully charged.
/// </summary>
/// <param name="vehicle">Vehicle object to compute waiting time 
/// (for possible different types of the vehicle extension).</param>
/// <returns>Returns computed waiting time.</returns>
double ChargingStation::GetChargingTime(double batteryLevel)
{
	return this->chargingWaitingTime_ * (1 - batteryLevel);
}


/// <summary>
/// Estimates waiting time of the vehicle which just aproached charging 
/// station until it starts charging.
/// </summary>
/// <returns>Return estimated waiting time in the queue for the charging.</returns>
double ChargingStation::GetExpectedWaitingTime()
{
	if (this->NumCustomers_ > this->Capacity_)
	// Number of customers is above the capacity -> add to estimated waiting 
	// time till start charging.
	{
		return this->chargingWaitingTime_ *
			(this->NumCustomers_ - this->Capacity_) * this->estimatedChargingLevel_;
	}

	// Vehicle comes straight to the charging station.
	return 0;
}