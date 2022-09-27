#ifndef MODELREPRESENTATION_H_
#define MODELREPRESENTATION_H_

#include "StationParameters.hpp"
#include "ChargingStation.hpp"


/// <summary>
/// Class to represent one model.
/// </summary>
class ModelRepresentation
{
public:
	// Loss of the model.
	double Loss_;
	// Vector of all charging stations of the model.
	std::vector<StationParameters> AllChargingStations_;

	ModelRepresentation();
	ModelRepresentation(double loss, std::vector<StationParameters> allStations);

	void SaveModelRepresentation(std::vector<std::unique_ptr<ChargingStation>>& allStations);
};


#endif // !MODELREPRESENTATION_H_
