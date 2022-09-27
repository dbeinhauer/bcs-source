#include "optimizer/Vehicle.hpp"


/// <summary>
/// Initializes object.
/// </summary>
/// <param name="startTime">Time of the departure.</param>
/// <param name="start">Starting position.</param>
/// <param name="end">Ending position.</param>
/// <param name="consumption">Battery consumption.</param>
/// <param name="batteryLevel">Start battery level.</param>
Vehicle::Vehicle(double startTime, MapPosition start, MapPosition end, 
	double consumption, double batteryLevel)
	: startTime_(std::move(startTime)), startPosition_(std::move(start)),
	endPosition_(std::move(end)), batteryConsumption_(std::move(consumption)),
	batteryLevel_(std::move(batteryLevel))
{
	// Current position is in the start position.
	this->currentPosition_ = this->startPosition_;
	// Starting in some segment.
	this->insideEdge_ = true;
	// By default we don't want to charge immediately.
	this->headingToChargingStation_ = false;
	// We haven't tried charging.
	this->chargingTested_ = false;

	this->chargingStationID_ = -1;
	this->previousVertex_ = -1;
	this->relativeVelocity_ = 0;
	this->startBatteryLevel_ = batteryLevel;
	this->startLineWaitingTime_ = 0;
	this->pathLength_ = 0;
}

Vehicle::~Vehicle() { }


/// <summary>
/// </summary>
/// <returns>Returns current battery level.</returns>
double Vehicle::GetBatteryLevel()
{
	return this->batteryLevel_;
}


/// <summary>
/// </summary>
/// <returns>Returns current vehicle position.</returns>
const MapPosition& Vehicle::GetPosition()
{
	return this->currentPosition_;
}


/// <summary>
/// </summary>
/// <returns>Returns charging station position.</returns>
const MapPosition& Vehicle::GetStationPosition()
{
	return this->stationPosition_;
}


/// <summary>
/// </summary>
/// <returns>Returns end position of the vehicle.</returns>
const MapPosition& Vehicle::GetEndPosition()
{
	return this->endPosition_;
}


double Vehicle::GetVelocity()
{
	return this->relativeVelocity_;
}



/// <summary>
/// Computes expected transition time for current traffic.
/// </summary>
/// <param name="map">Map object to compute transit time from.</param>
/// <param name="segment">ID of the segment from which we are begining.</param>
/// <param name="increasing">Flag if increasing in segment ID 
/// (if we want to get from end of the edge to the segment, then opposite way).</param>
/// <param name="secondSegment">ID of the final segment in the edge 
/// (if movent to the end of the edge, then set attribute to `-1`)</param>
/// <returns>Returns transition time the vehicle needs to get from current 
/// node to the next node.</returns>
double Vehicle::GetNextNodeTransitTime(
	Map& map, int32_t segment, bool increasing, int32_t secondSegment)
{
	Graph& graph = map.GetGraph();
	double currentRelativeTime =
		graph[this->currentPosition_.EdgeID].TransitTime;

	int32_t numSegments = graph[this->currentPosition_.EdgeID].GetNumSegments();
	if (segment != 0 && segment != numSegments - 1)
	// If vehicle not in the first or last segment -> transit time is lower.
	{
		// If decreasing in segment -> pass exactly `segment` segments.
		int32_t segmentsToPass = segment;
		if (increasing)
		// Increasing in segment ID 
		// (minus 1, because we are already in one of the segments).
		{
			segmentsToPass = numSegments - segment - 1;
		}

		if (secondSegment != -1)
		// Movement inside the edge.
		{
			segmentsToPass = std::abs(segment - secondSegment);
		}

		// Proportionaly diminish the transition time.
		currentRelativeTime *= ((double)segmentsToPass / (double)numSegments);
	}

	// Return travel time based on the Vehicle speed.
	return currentRelativeTime / this->relativeVelocity_;
}

/// <summary>
/// </summary>
/// <param name="pathLength">Lenght of the path to pass.</param>
/// <returns>Returns expected transition time based on the path length.</returns>
double Vehicle::GetExpectedPathTransitTime(double pathLength)
{
	return pathLength / this->relativeVelocity_;
}


/// <summary>
/// Computes expected battery level, if crossing the edge from given segment.
/// Works also if we want to get to given segment, but in `incresing` should be 
/// flag for the direction from the `segment` to end of the edge. 
/// </summary>
/// <param name="map">Map to find corresponding edge info.</param>
/// <param name="segment">ID of the segment from which we are begining.</param>
/// <param name="increasing">Flag if increasing in segment ID 
/// (if we want to get from end of the edge to the segment, then opposite way).</param>
/// <param name="secondSegment">ID of the final segment in the edge 
/// (if movent to the end of the edge, then set attribute to `-1`)</param>
/// <returns>Returns expected battery level after crossing the edge.</returns>
double Vehicle::GetExpectedEdgeTransitBatteryLevel(
	Map& map, int32_t segment, bool increasing, int32_t secondSegment)
{
	return this->batteryLevel_ -
		this->batteryConsumption_ * 
		this->GetNextNodeTransitTime(map, segment, increasing, secondSegment);
}


/// <summary>
/// Computes estimated battery level after reaching the path end in estimated time.
/// </summary>
/// <param name="estimatedTime">Estimated time to reach the end of the road.</param>
/// <returns>Returns estimated end battery level.</returns>
double Vehicle::GetExpectedPathBatteryLevel(double estimatedTime)
{
	return this->batteryLevel_ - this->batteryConsumption_ * estimatedTime;
}


/// <summary>
/// </summary>
/// <returns>Returns ID of the charging station where the vehicle 
/// is heading or where currently is.</returns>
int32_t Vehicle::GetChargingStationID()
{
	return this->chargingStationID_;
}

/// <summary>
/// Computes starting segment ID (when not starting inside the edge).
/// Decides whether start from the begin or the end side of the edge.
/// </summary>
/// <param name="isDecreasing">Flag indicating whether the vehicle is moving in 
/// decreasing segment IDs.</param>
/// <param name="edge">The current edge object (to get number of segments).</param>
/// <returns>Returns starting segment ID.</returns>
int32_t Vehicle::GetStartingSegment(bool isDecreasing, Edge& edge)
{
	// Vehicle moving from the lower vertex to higher 
	// -> first index of the segment is 0.
	if (isDecreasing)
		return 0;

	// Vehicle moving from higher to lower vertex 
	// -> first index of the segment is `num_segments - 1`
	return edge.GetNumSegments() - 1;
}


/// <summary>
/// </summary>
/// <param name="start">Start vertex identifier.</param>
/// <param name="end">End vertex identifier.</param>
/// <returns>Returns `true` if vehicle is moving in increasing order 
/// with the segments, else if decreasing `false`</returns>
bool Vehicle::IsIncreasing(vertex_t start, vertex_t end)
{
	// The vehicle is moving from the vertex with smaller ID 
	// to higher one.
	if (start < end)
		return true;

	// Vehicle is moving towards vertex with lower ID.
	return false;
}


/// <summary>
/// Finds out whether car changed the city.
/// </summary>
/// <param name="previous">Previous vertex.</param>
/// <param name="current">Current vertex.</param>
/// <param name="map">Map object to get info from.</param>
/// <returns>Returns `true` if vehicle crossed the border of the city,
///  else `false`.</returns>
bool Vehicle::ChangedCity(Map& map)
{
	Graph& graph = map.GetGraph();

	// Vehicle crossed the city border.
	if (graph[this->previousVertex_].CityID_ !=
		graph[currentPosition_.FurtherVertexID].CityID_)
		return true;

	// Vehicle didn't cross the city border.
	return false;
}


/// <summary>
/// </summary>
/// <returns>Returns `true` if vehicle is heading to the charging station,
///  else `false`.</returns>
bool Vehicle::GoingToStation()
{
	return this->headingToChargingStation_;
}


/// <summary>
/// </summary>
/// <returns>Returns start time of the vehicle 
/// (either of the start or when it starts returning).</returns>
double Vehicle::GetStartTime()
{
	return this->startTime_;
}


/// <summary>
/// </summary>
/// <returns>Returns battery level at the start of the path.</returns>
double Vehicle::GetStartBatteryLevel()
{
	return this->startBatteryLevel_;
}


/// <summary>
/// </summary>
/// <returns>Returns time of the start waiting to go charging.</returns>
double Vehicle::GetStartLineWaitingTime()
{
	return this->startLineWaitingTime_;
}


/// <summary>
/// </summary>
/// <returns>Returns number of vertices stored in the precomputed path.</returns>
int32_t Vehicle::GetPathSize()
{
	return this->path_.size();
}


/// <summary>
/// </summary>
/// <returns>Returns `true` if car already tested if it should go charging,
///  else `false`.</returns>
bool Vehicle::TriedCharging()
{
	return this->chargingTested_;
}

/// <summary>
/// </summary>
/// <returns>Returns `true` if vehicle is inside the edge 
/// (segment not last or first), else `false`.</returns>
bool Vehicle::IsInsideEdge()
{
	return this->insideEdge_;
}



/// <summary>
/// </summary>
/// <returns>Returns `true` if vehicle is in the last step of the road 
/// or in the penultimate step, which is the further vertex 
/// (just get to the correct segment), (finish or charging station).</returns>
bool Vehicle::IsFinishing()
{
	// Get correct finishing position.
	MapPosition finishingPosition = this->endPosition_;
	if (this->GoingToStation())
		// Vehicle heading to the charging station.
	{
		finishingPosition = this->stationPosition_;
	}

	if (this->currentPosition_.EdgeID == finishingPosition.EdgeID)
		// If current edge same as end -> just go to the correct segment.
	{
		// Vehicle is finishing.
		return true;
	}

	// If no vertex in the path -> we get to the end -> get to correct segment
	return this->path_.empty();
}


/// <summary>
/// Sets proper flags to start going to the charging station.
/// </summary>
void Vehicle::StartGoingCharging()
{
	this->headingToChargingStation_ = true;
}


/// <summary>
/// Sets start time of the waiting in the charging station queue.
/// </summary>
/// <param name="startTime">Start time of the waiting.</param>
void Vehicle::SetStartLineWaitingTime(double startTime)
{
	this->startLineWaitingTime_ = startTime;
}

/// <summary>
/// Sets whether vehicle tried to go to the charging station.
/// </summary>
/// <param name="value">Flag if tried `true`, else `false`.</param>
void Vehicle::SetChargingTry(bool value)
{
	this->chargingTested_ = value;
}


/// <summary>
/// Decides whether go charging or not (based on the battery level, using 
/// randomness and estimated battery level in the end (if enough -> don't go charging). 
/// </summary>
/// <param name="start">Start vertex of the path.</param>
/// <param name="end">End vertex of the path.</param>
/// <param name="map">Map object to find the path from.</param>
/// <param name="batteryTreshold">Randomly generated additional battery 
/// level for decision if go charging or not 
/// (if level above current level then go).</param>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <returns>Returns `true` if go to the charging station, else `false`.</returns>
bool Vehicle::GoCharging(vertex_t start, vertex_t end, Map& map, 
	double batteryTreshold, SimulationParameters& simulationParameters)
{
	// Finds shortest path to the end.
	auto lengthPathPair = map.FindShortestPath(start, end);
	// Estimates battery level in the end of the road.
	double expectedBatteryLevel = this->GetExpectedPathBatteryLevel(
		this->GetExpectedPathTransitTime(lengthPathPair.first));

	if (expectedBatteryLevel > simulationParameters.BatteryTolerance)
	// Expected battery level in the end of the road is above the tolerance
	// -> go straight to the finish.
	{
		return false;
	}

	// Go charging.
	return true;
}


/// <summary>
/// Decides whether go charging or not (based on the battery level, using 
/// randomness and estimated battery level in the end 
/// (if enough -> don't go charging). Works without need to compute the path.
/// </summary>
/// <param name="batteryTreshold">Randomly generated additional battery 
/// level for decision if go charging or not 
/// (if level above current level then go).</param>
/// <param name="simulationParameters">Parameters of the simulation.</param>
/// <returns>Returns whether go to the charging station `true` or not `false`.</returns>
bool Vehicle::GoCharging(
	double batteryTreshold, SimulationParameters& simulationParameters)
{
	if (this->GetExpectedPathBatteryLevel(
		this->GetExpectedPathTransitTime(this->pathLength_)) > 
		simulationParameters.BatteryTolerance)
	// Expected battery level in the end of the road is above the tolerance
	// -> go straight to the finish.
	{
		return false;
	}

	// Go charging.
	return true;
}


/// <summary>
/// Decides whether it is time to charge the battery or not 
/// (starts moving to the charging station) only based on the battery level 
/// (doesn't count with the expected battery level in the end).
/// </summary>
/// <param name="batteryTreshold">Randomly generated additional battery level 
/// for decision if go charging or not (if level above then go).</param>
/// <returns>Returns `true` if vehicle should go charging, else `false`.</returns>
bool Vehicle::BatteryLevelToGoCharging(
	double bottomTreshold, double upperTreshold, double additionalTreshold)
{
	// Battery level under treshold to always charge.
	if (this->batteryLevel_ < bottomTreshold + additionalTreshold)
		return true;

	// Go to the end (don't charge).
	return false;
}


/// <summary>
/// Updates vehicle path and decides whether it is neccesary to go 
/// to the charging station. If so then finds the path to the best 
/// estimated charging station.
/// </summary>
/// <param name="start">Start vertex of the path.</param>
/// <param name="end">Final vertex of the path.</param>
/// <param name="batteryTreshold">Randomly generated additional battery level
/// for decision if go charging or not (if level above then go).</param>
/// <param name="map">Map object to compute path from.</param>
/// <param name="removeFirst">Flag if remove first (current) vertex from the 
/// path (useful when not inside edge). If remove then `true` 
/// (if path has length 1, then don't remove the vertex (we want to know that
/// it is end of the road), else `false`.</param>
void Vehicle::UpdateVehiclePath(vertex_t start, vertex_t end, Map& map,
	double batteryTreshold, bool removeFirst, 
	SimulationParameters& simulationParameters)
{
	// Decide whether go charging or not.
	bool goCharging = 
		this->GoCharging(start, end, map, batteryTreshold, simulationParameters);

	// Find new path.
	auto lengthPathPair = map.FindVehiclePath(start, end, goCharging);
	// Position of the charging station where we are heading.
	this->stationPosition_ = (*map.GetLastChargingStation()).Position_;
	// ID of the charging station where we are heading.
	this->chargingStationID_ = (*map.GetLastChargingStation()).StationID_;
	// Save the path information.
	this->pathLength_ = lengthPathPair.first;
	this->path_ = std::move(lengthPathPair.second);


	if (goCharging)
	// Start heading to the charging station.
	{
		this->headingToChargingStation_ = true;
	}

	if (removeFirst && this->path_.size() > 1)
	// Remove redundant start vertex and change next vehicle move
	// (based on the path update).
	{
		this->path_.pop_back();
		// Change next edge based on the update.
		this->currentPosition_.ChangeParameters(this->path_.back(),
			this->GetPosition().FurtherVertexID,
			-1,
			map.GetGraph());
	}
}


/// <summary>
///	Increases the battery level (maximal capacity is 1) and stops charging.
/// </summary>
/// <param name="chargeLevel">Battery percentage to be charged 
/// (needs to be non-negative value).</param>
void Vehicle::ChargeBattery(double chargeLevel)
{
	// Articulary decrease battery level (when waiting on the station).
	if (chargeLevel < 0)
	{
		this->batteryLevel_ += chargeLevel;
		return;
	}

	// Charge the battery.
	this->batteryLevel_ += chargeLevel;


	if (this->batteryLevel_ > 1)
	// Battery level overflow check 
	// -> if overflow then set battery level to 1 (100%). 
	{
		this->batteryLevel_ = 1;
	}

	// Stops charging.
	this->headingToChargingStation_ = false;
}



/// <summary>
/// Moves the vehicle from the starting position to the closest node 
/// (either starting or leaving the charging station).
/// </summary>
/// <param name="map">Map to get information about the route.</param>
/// <param name="logs">Whether to write simulation logs (for debugging).</param>
/// <returns>Returns `true` if moving was successful, 
/// else `false (battery run down).</returns>
bool Vehicle::MoveFirstSegment(Map& map, bool logs)
{
	// Shouldn't be there (path is empty) -> possible error
	if (this->path_.empty())
		return false;

	// Check whether vehicle is in the start of the road -> else error
	if (this->path_.back() != this->startPosition_.CloserVertexID)
		return false;


	// If path has at least 2 nodes (not path inside the edge) 
	// -> pop starting vertex
	if (this->path_.size() > 1)
		this->path_.pop_back();


	// Expected battery level after crossing the edge.
	double expectedBattery = this->GetExpectedEdgeTransitBatteryLevel(
		map,
		this->currentPosition_.EdgeSegmentID,
		// Opposite direction than the reality 
		// (because moving to the middle of the edge).
		this->IsIncreasing(this->currentPosition_.FurtherVertexID,
			this->currentPosition_.CloserVertexID),
		-1);


	// Move to the end of the edge (closest vertex).
	this->insideEdge_ = false;


	// Update info about the battery levels.
	map.SimulateVehiclePassedThroughEdge(this->currentPosition_.EdgeID,
		std::make_pair(this->batteryLevel_, expectedBattery),
		this->currentPosition_.EdgeSegmentID,
		// Opposite direction than the reality 
		// (because moving to the middle of the edge).
		this->IsIncreasing(this->currentPosition_.FurtherVertexID,
			this->currentPosition_.CloserVertexID),
		-1);


	// Write info.
	if (logs)
		std::cout << "Battery level: " << expectedBattery << std::endl;


	if (expectedBattery < 0)
	// Vehicle run down with battery, but we still want to propagate information 
	// about the battery levels (if it is negative -> more severe problem and with
	// higher probability ideal place of the station).
	{
		// Write info.
		if (logs)
			std::cout << "BATTERY RUN DOWN!" << std::endl;
		
		// Battery run down (crossing failed).
		return false;
	}

	// Update previous vertex.
	this->previousVertex_ = this->currentPosition_.CloserVertexID;


	// Update current position.

	// Decide whether the final position is the end of the road or station.
	MapPosition currentLastPosition = this->endPosition_;
	if (this->GoingToStation())
	// Vehicle going to the charging station.
	{
		currentLastPosition = this->stationPosition_;
	}


	if (this->path_.back() == this->currentPosition_.CloserVertexID)
	// Last part of the path.
	{
		this->currentPosition_.ChangeParameters(
			currentLastPosition.FurtherVertexID,
			currentLastPosition.CloserVertexID,
			-1,
			map.GetGraph());
	}
	else
	// Not last part of the path.
	{
		this->currentPosition_.ChangeParameters(
			this->path_.back(),
			this->currentPosition_.CloserVertexID,
			-1,
			map.GetGraph());
	}

	// Update battery level of the car.
	this->batteryLevel_ = expectedBattery;

	// Car successfully passed the edge.
	return true;
}


/// <summary>
/// Moves the vehicle to the next node (without the last and first move).
/// </summary>
/// <param name="map">Map to get information about the route.</param>
/// <param name="logs">Whether to write simulation logs (for debugging).</param>
/// <returns>Returns `true` if moving was successful, 
/// else `false (battery run down).</returns>
bool Vehicle::MoveToNextNode(Map& map, bool logs)
{
	// In the end of the road -> vehicle shouldn't be in this function 
	// -> probably error
	if (this->path_.empty())
		return false;


	// Expected battery level after crossing the edge.
	double expectedBattery = this->GetExpectedEdgeTransitBatteryLevel(
		map, 
		this->currentPosition_.EdgeSegmentID,
		this->IsIncreasing(this->currentPosition_.FurtherVertexID,
			this->currentPosition_.CloserVertexID), 
		-1);


	// Always moves to the end of the edge 
	// (final part of the road is checked by other function).
	this->insideEdge_ = false;

	// Update info about the battery levels.
	map.SimulateVehiclePassedThroughEdge(this->currentPosition_.EdgeID,
		std::make_pair(this->batteryLevel_, expectedBattery),
		this->currentPosition_.EdgeSegmentID,
		this->IsIncreasing(this->currentPosition_.FurtherVertexID,
			this->currentPosition_.CloserVertexID),
		-1);

	// Write info.
	if (logs)
		std::cout << "Battery level: " << expectedBattery << std::endl;


	if (expectedBattery < 0)
	// Vehicle run down with battery, but we still want to propagate information 
	// about the battery levels (if it is negative -> more severe problem and with
	// higher probability ideal place of the station).
	{
		// Write info.
		if (logs)
			std::cout << "BATTERY RUN DOWN!" << std::endl;
		
		// Battery run down (crossing failed).
		return false;
	}


	// Update previous vertex.
	this->previousVertex_ = this->currentPosition_.CloserVertexID;
	
	// Update path length.
	this->pathLength_ -= this->GetNextNodeTransitTime(map,
		this->currentPosition_.EdgeSegmentID,
		this->IsIncreasing(this->currentPosition_.FurtherVertexID,
			this->currentPosition_.CloserVertexID),
		-1);


	// Decide whether the final position is the end of the road or station.
	MapPosition currentEndPosition = this->endPosition_;
	if (this->GoingToStation())
	// Vehicle going to the charging station.
	{
		currentEndPosition = this->stationPosition_;
	}


	// If path has at least 2 nodes (not path inside the edge) 
	// -> pop current vertex
	if (this->path_.size() > 1)
		this->path_.pop_back();




	if (this->path_.back() == this->currentPosition_.CloserVertexID)
	// Last part of the path.
	{
		// Update current position.
		this->currentPosition_.ChangeParameters(
			currentEndPosition.FurtherVertexID,
			currentEndPosition.CloserVertexID,
			-1,
			map.GetGraph());
	}
	else
	// Part inside the path.
	{
		// Update current position.
		this->currentPosition_.ChangeParameters(
			this->path_.back(),
			this->currentPosition_.CloserVertexID,
			-1,
			map.GetGraph());
	}

	// Car successfully passed the edge.
	return true;
}



/// <summary>
/// Moves the vehicle on the final segment of the path 
/// (moving to the middle the edge).
/// </summary>
/// <param name="map">Map to get information about the route.</param>
/// <param name="logs">Whether to write simulation logs (for debugging).</param>
/// <returns>Returns `true` if moving was successful, 
/// else `false (battery run down).</returns>
bool Vehicle::MoveToFinalSegment(Map& map, bool logs)
{
	// Decide whether the final position is the end of the road or station.
	MapPosition currentEndPosition = this->endPosition_;
	if (this->GoingToStation())
		// Vehicle going to the charging station.
	{
		currentEndPosition = this->stationPosition_;
	}

	// Check movement inside the edge.
	int secondSegment = -1;
	int numSegments = 
		map.GetGraph()[this->currentPosition_.EdgeID].GetNumSegments();
	if (this->currentPosition_.EdgeSegmentID != numSegments - 1
		&& this->currentPosition_.EdgeSegmentID != 0)
	// Moving inside of the edge -> set also `secondSegment`
	{
		secondSegment = this->currentPosition_.EdgeSegmentID;
	}


	// Expected battery level after crossing the edge.
	double expectedBattery = this->GetExpectedEdgeTransitBatteryLevel(
		map,
		currentEndPosition.EdgeSegmentID,
		// Opposite direction than the reality
		// (because moving to the middle of the edge).
		this->IsIncreasing(currentEndPosition.FurtherVertexID,
			currentEndPosition.CloserVertexID),
		secondSegment);


	// Move to the middle of the edge (charging station or the end).
	this->insideEdge_ = true;


	// Update info about the battery levels.
	map.SimulateVehiclePassedThroughEdge(currentEndPosition.EdgeID,
		std::make_pair(this->batteryLevel_, expectedBattery),
		currentEndPosition.EdgeSegmentID,
		// Opposite direction than the reality 
		// (because moving to the middle of the edge).
		this->IsIncreasing(currentEndPosition.FurtherVertexID,
			currentEndPosition.CloserVertexID),
		secondSegment);


	// Write info.
	if (logs)
		std::cout << "Battery level: " << expectedBattery << std::endl;

	if (expectedBattery < 0)
	// Vehicle run down with battery, but we still want to propagate information 
	// about the battery levels (if it is negative -> more severe problem and with
	// higher probability ideal place of the station).
	{
		// Write info.
		if (logs)
			std::cout << "BATTERY RUN DOWN!" << std::endl;

		// Battery run down (crossing failed).
		return false;
	}


	// Update previous vertex.
	this->previousVertex_ = currentEndPosition.CloserVertexID;
	// Update path length -> end of the precomputed road.
	this->pathLength_ = 0;
	// Update current position.
	this->currentPosition_ = std::move(currentEndPosition);
	// Update battery level of the car.
	this->batteryLevel_ = expectedBattery;

	// Car successfully passed the edge.
	return true;
}


/// <summary>
/// Removes first vertex of the path 
/// (only use when we want to remove current vertex of the path start).
/// </summary>
void Vehicle::PopFirstPathVertex()
{
	this->path_.pop_back();
}