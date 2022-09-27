#include "optimizer/ModelRepresentation.hpp"


/// <summary>
/// Initializes the object.
/// </summary>
ModelRepresentation::ModelRepresentation() 
{
	this->Loss_ = 0;
}


/// <summary>
/// Initilizes the object.
/// </summary>
/// <param name="loss">Loss of the model.</param>
/// <param name="allStations">Vector of all charging stations of the model.</param>
ModelRepresentation::ModelRepresentation(double loss,
	std::vector<StationParameters> allStations)
	: Loss_(std::move(loss)), AllChargingStations_(std::move(allStations))
{ }


/// <summary>
/// Stores representation of the given model to `this` object.
/// </summary>
/// <param name="allStations">Container of all charging stations of the model.</param>
void ModelRepresentation::SaveModelRepresentation(
	std::vector<std::unique_ptr<ChargingStation>>& allStations)
{
	this->AllChargingStations_.clear();
	for (auto&& station : allStations)
	// Store all stations representations.
	{
		this->AllChargingStations_.push_back(
			StationParameters((*station).StationID_,
				(*station).Position_,
				(*station).Capacity_,
				(*station).CityID_,
				(*station).GetChargingTime(0),
				(*station).GetExpectedWaitingTime()));
	}
}