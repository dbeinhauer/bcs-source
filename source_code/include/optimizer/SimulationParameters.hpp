#ifndef SIMULATIONPARAMETERS_H_
#define SIMULATIONPARAMETERS_H_

#include <iostream>

/// <summary>
/// Class to manage program arguments.
/// </summary>
class SimulationParameters
{
public:
    bool SavePrepared;
    int32_t SimulationTime, NumClosestStations, NumStations, StationCapacity;
    double SegmentLength, CarConsumption, ExponentialLambdaCities, 
        ExponentialLambdaDepartures, EndCityRatio, BatteryTresholdLambda, 
        CarBatteryMean, CarBatteryDeviation, CarStartBatteryBottomLimit,
        ChargingTreshold, NotChargingTreshold, BatteryTolerance, CarVelocity,
        ChargingWaitingTime, MeanChargingLevel;
    std::string EdgesFile, NodeCityFile, CitiesFile, AddedEdgesFile, 
        PreparedNodesFile, PreparedEdgesFile;

    SimulationParameters();

};


#endif // !SIMULATIONPARAMETERS_H_

