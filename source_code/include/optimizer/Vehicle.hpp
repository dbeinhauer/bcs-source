#ifndef VEHICLE_H_
#define VEHICLE_H_


#include <random>

#include "SimulationParameters.hpp"
#include "Map.hpp"
#include "MapPosition.hpp"


/// <summary>
/// Parent class of all vehicle objects for the simulation.
/// </summary>
class Vehicle
{
public:
	Vehicle(double startTime, MapPosition start, MapPosition end, 
		double consumption, double batteryLevel);

	~Vehicle();

	// Get information about the vehicle:
	double GetBatteryLevel();
	const MapPosition& GetPosition();
	const MapPosition& GetStationPosition();
	const MapPosition& GetEndPosition();

	double GetVelocity();

	double GetNextNodeTransitTime(
		Map& map, int32_t segment, bool increasing, int32_t secondSegment);
	double GetExpectedPathTransitTime(double pathLength);
	double GetExpectedEdgeTransitBatteryLevel(
		Map& map, int32_t segment, bool increasing, int32_t secondSegment);
	double GetExpectedPathBatteryLevel(double estimatedTime);

	int32_t GetChargingStationID();

	int32_t GetStartingSegment(bool isDecreasing, Edge& edge);
	bool IsIncreasing(vertex_t start, vertex_t end);

	bool ChangedCity(Map& map);
	bool GoingToStation();

	double GetStartTime();
	double GetStartBatteryLevel();
	double GetStartLineWaitingTime();

	int32_t GetPathSize();
	bool TriedCharging();

	bool IsInsideEdge();
	bool IsFinishing();


	// Set properties of vehicle:
	void StartGoingCharging();
	void SetStartLineWaitingTime(double startTime);
	void SetChargingTry(bool value);


	// Vehicle get decisions:
	bool GoCharging(vertex_t start, vertex_t end, Map& map, 
		double batteryTreshold, SimulationParameters& simulationParameters);
	bool GoCharging(
		double batteryTreshold, SimulationParameters& simulationParameters);
	bool BatteryLevelToGoCharging(double bottomTreshold,
		double upperTreshold, double additionalTreshold);


	// Vehicle actions:
	void UpdateVehiclePath(vertex_t start, vertex_t end, Map& map, 
		double batteryTreshold, bool removeFirst, 
		SimulationParameters& simulationParameters);

	void ChargeBattery(double chargeLevel);

	bool MoveFirstSegment(Map& map, bool logs);
	bool MoveToNextNode(Map& map, bool logs);
	bool MoveToFinalSegment(Map& map, bool logs);

	void PopFirstPathVertex();

protected:
	// Edges and its segments of the start, end and current position.
	// Start or end vertex for closest path algorithms is the closser 
	// one (based on the segments).
	MapPosition startPosition_;
	MapPosition stationPosition_;
	MapPosition endPosition_;
	MapPosition currentPosition_;

	// Departure time.
	double startTime_;
	// Time of the arrival to the station (and starts waiting in the line).
	double startLineWaitingTime_;

	// ID of the station where the vehicle is heading.
	int32_t chargingStationID_;

	// Vertex from where is the vehicle currently heading.
	vertex_t previousVertex_;

	// Flag indicating if car is currently inside the segment 
	// (starting, finishing or going from the charging station).
	bool insideEdge_;

	// Battery info:
	double startBatteryLevel_;
	double batteryConsumption_;
	double batteryLevel_;

	// Speed of the vehicle.
	double relativeVelocity_;

	// Flag if vehicle is already going to the charging station.
	bool headingToChargingStation_;

	// Flag if vehicle already testing charging variant 
	// (won't try to go charging anymore).
	bool chargingTested_;

	// Precomputed reversed path of the vehicle (last member is next edge).
	std::vector<vertex_t> path_;

	// Aproximate path travel time (considering traffic during path finding).
	double pathLength_;
};

#endif