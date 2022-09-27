#ifndef EDGE_H_
#define EDGE_H_


#include <iostream>
#include <vector>
#include <tuple>
#include <float.h>



/// <summary>
/// Supported road types and its constants for road capacity computation 
/// (higher values means higher capacity).
/// </summary>
enum PermitedRoadTypes
{
	MOTORWAY = 10,
	TRUNK = 4,
	PRIMARY = 2,
	OTHER = 1
};


// For storing number of cars and sum of battery levels.
typedef std::pair<uint64_t, double> BateryPair;


/// <summary>
/// Class for edge representation (road between the junctions).
/// </summary>
class Edge
{
public:
	// Constant for computing transition time when traffic is
	// higher than capacity (simulation of the traffic jams).
	static double CapacityOverflowConstant;

	// Current time needed by the vehicle to pass through the edge
	// (time is in minutes, we consider vehicle velocity is 1 km/minute,
	// its possible to appropriately adjust transition time for different
	// types and speeds of the vehicle).
	double TransitTime;


	Edge();

	// Reset simulation info.
	void ResetSimulation();

	// Edge properties setup.
	void AddType(char type);
	void AddType(PermitedRoadTypes type);
	void SetLength(double length, double segmentLength);

	// Get information about the edge.
	double GetLength();
	int32_t GetNumSegments();
	PermitedRoadTypes GetRoadType();
	char GetRoadTypeChar();
	std::vector<BateryPair>& GetBatteryLevels();

	// Update edge information.
	void UpdateTransitTime(bool increase);
	bool UpdateSegmentData(std::pair<double, double> batteryPair,
		int32_t startSegment, bool segmentIncrease, int32_t endSegment);

private:
	// Length of the edge in kilometers.
	double length_;
	// Number of segments the edge is divided into.
	int32_t numSegments_;
	// Current number of cars on the edge.
	uint32_t currentTraffic_;
	// Capacity of the whole edge 
	// (num of cars on the edge until traffic jam starts to create).
	uint32_t capacity_;
	// Edge road type.
	PermitedRoadTypes type_;

	// Vector of the battery levels for each segment (indexing from the vertex with lower ID).
	// (to compute average battery level on the position).
	std::vector<BateryPair> bateryLevels_;


	void setCapacity();

	PermitedRoadTypes getRoadType(char type);

	double computeTransitTime();
};

#endif