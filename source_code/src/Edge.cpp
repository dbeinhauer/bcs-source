#include "optimizer/Edge.hpp"


double Edge::CapacityOverflowConstant = 1;


/// <summary>
/// Initializes the edge.
/// </summary>
Edge::Edge()
{
	// Just init the variables.
	this->TransitTime = 0;
	this->length_ = 0;
	this->numSegments_ = 0;
	this->currentTraffic_ = 0;
	this->type_ = PermitedRoadTypes::OTHER;

	// By default is number of segments -1 to check whether the edge is initialized.
	this->numSegments_ = -1;
}


/// <summary>
/// Resets all simulation information.
/// </summary>
void Edge::ResetSimulation()
{
	this->currentTraffic_ = 0;
	std::fill(this->bateryLevels_.begin(), this->bateryLevels_.end(), BateryPair());
	this->computeTransitTime();
}


/// <summary>
/// Sets proper edge type based on the input.
/// </summary>
/// <param name="type">Type of the road
/// (supported types: `m` - motorway, `t` - trunk, `p` - primary, `o` - other).</param>
void Edge::AddType(char type)
{

	this->type_ = this->getRoadType(std::move(type));
}


/// <summary>
/// Sets proper edge type based on the input.
/// </summary>
/// <param name="type">Type of the road.</param>
void Edge::AddType(PermitedRoadTypes type)
{
	this->type_ = type;
}


/// <summary>
/// Sets length of the edge in kilometers, prepares segmentation of the edge 
/// (for battery levels and car generations) and sets capacity of the road.
/// </summary>
/// <param name="length">Length of the edge in kilometers.</param>
void Edge::SetLength(double length, double segmentLength)
{
	this->length_ = length;

	this->numSegments_ = this->length_ / segmentLength;

	// Make sure there is at least 1 segment (if edge is too short).
	if (this->numSegments_ == 0)
		this->numSegments_ = 1;

	this->bateryLevels_ = std::vector<BateryPair>(
		this->numSegments_, std::make_pair(0, 0.0));
	this->setCapacity();

}


/// <summary>
/// </summary>
/// <returns>Returns length of the edge.</returns>
double Edge::GetLength()
{
	return this->length_;
}


/// <summary>
/// </summary>
/// <returns>Returns number of segments of the edge.</returns>
int32_t Edge::GetNumSegments()
{
	return this->numSegments_;
}


/// <summary>
/// </summary>
/// <returns>Returns type of the road.</returns>
PermitedRoadTypes Edge::GetRoadType()
{
	return this->type_;
}


/// <summary>
/// </summary>
/// <returns>Returns type of the road in char representation 
/// (`m` - motorway, `t` - trunk, `p` - primary, `o` - other).</returns>
char Edge::GetRoadTypeChar()
{
	switch (this->type_)
	{
		case MOTORWAY:
			return 'm';
		case TRUNK:
			return 't';
		case PRIMARY:
			return 'p';
		default:
			return 'o';
	}
}


/// <summary>
/// </summary>
/// <returns>Returns reference to contatiner of battery levels data of the edge.</returns>
std::vector<BateryPair>& Edge::GetBatteryLevels()
{
	return this->bateryLevels_;
}


/// <summary>
/// Updates transition time based on the change of the traffic 
/// (vehicle either left or enter the edge).
/// </summary>
/// <param name="increase">If vehicle entered the edge `true`, 
/// else if vehicle left `false`.</param>
void Edge::UpdateTransitTime(bool increase)
{
	if (increase)
	// Vehicle entered the edge.
	{
		this->currentTraffic_++;
	}

	else
	// Vehicle left the edge.
	{
		// Check underflow error (negative number of the cars).
		if (this->currentTraffic_ > 0)
			this->currentTraffic_--;
	}

	// Update transition time.
	this->TransitTime = this->computeTransitTime();
}


/// <summary>
/// Updates segments information about the battery levels.
/// </summary>
/// <param name="batteryPair">Pair of start and end battery levels.</param>
/// <param name="startSegment">Index of the start segment.</param>
/// <param name="segmentIncrease">Flag if vehicle moving in increasing 
/// order of the segment indices `true` or decreasing order `false`.</param>
/// <param name="endSegment">Index of the end segment 
/// (if vehicle movement ends in the middle of the edge),
/// (if -1, then automaticaly choosen the appropriate end of the edge).</param>
/// <returns>Returns `true` if update was successful, else `false`.</returns>
bool Edge::UpdateSegmentData(std::pair<double, double> batteryPair,
	int32_t startSegment, bool segmentIncrease, int32_t endSegment)
{
	double batteryChange = batteryPair.first - batteryPair.second;

	// Compute number of segments to go through.
	int32_t numPassingSegments;
	if (segmentIncrease)
		// Going to higher segment indices,
		numPassingSegments = this->numSegments_ - startSegment;
	else
		// Going to lower segment indices.
		numPassingSegments = startSegment + 1;

	if (endSegment != -1)
	// Segment update should end in the middle of the edge.
	{
		// Bad end segment identifier -> error
		if (segmentIncrease && startSegment > endSegment)
			return false;

		// Bad end segment identifier -> error
		if (!segmentIncrease && startSegment < endSegment)
			return false;

		// Compute number of segments to go throught based on the start and end segment ID 
		// (we add 1 because we want also include segment with the ID `endSegment`).
		numPassingSegments = std::abs(startSegment - endSegment) + 1;
	}

	// Decrease of the battery level for each segment.
	double batterySegmentChange = batteryChange / numPassingSegments;


	int32_t segmentID = startSegment;
	double currentBattery = batteryPair.first;

	// Compute battery level for each segment and adds it with
	// the sum of battery levels of all vehicles that crossed the edge.
	for (size_t i = 0; i < numPassingSegments; i++)
	{
		// Error battery level (can't be negative) -> Error
		if (currentBattery < 0)
			return false;

		// Update the levels.
		this->bateryLevels_[segmentID].first++;
		this->bateryLevels_[segmentID].second += currentBattery;

		// Move to next segment.
		if (segmentIncrease)
			segmentID++;
		else
			segmentID--;

		// Update battery level for the next segment.
		currentBattery -= batterySegmentChange;
	}

	// Update was successful.
	return true;
}



/// <summary>
/// Computes and sets the capacity of the road.
/// Computes capacity using equation:
///		capacity = road_type_costant * length
/// </summary>
void Edge::setCapacity()
{
	this->capacity_ = this->type_ * this->length_;
	
	// Too small capacity (set it to 1 vehicle cacity).
	if (this->capacity_ < 1)
		this->capacity_ = 1;

	this->TransitTime = this->computeTransitTime();
}




/// <summary>
/// Based on given character returns apropriate road type.
/// </summary>
/// <param name="type">Type of the road
/// (supported types: `m` - motorway, `t` - trunk, `p` - primary).</param>
/// <returns>Returns apropriate road type.</returns>
PermitedRoadTypes Edge::getRoadType(char type)
{
	switch (type)
	{
		case 'm':
			return MOTORWAY;
		case 't':
			return TRUNK;
		case 'p':
			return PRIMARY;
		default:
			return OTHER;
	}
}


/// <summary>
/// Computes current transition time.
/// Computed:
///		transition_time = length + capacity_overflow_constant * ((current_traffic - capacity) / capacity)
/// </summary>
/// <returns>Returns current transition time in minutes of the 
/// vehicle with speed 1 kilimeter per minute.</returns>
double Edge::computeTransitTime()
{
	// Incorrect edge -> transiton length is infinite (ignore the road).
	if (this->capacity_ <= 0)
		return DBL_MAX;

	//// If not above capacity -> just return length (speed is in km/minutes).
	if (this->currentTraffic_ <= this->capacity_)
		return this->length_;

	// Computed how many vehicles there are above the capacity of the road.
	uint32_t aboveCapacity = this->currentTraffic_ - this->capacity_;

	// Above capacity -> compute transition time as follows.
	return this->length_ + this->CapacityOverflowConstant * (aboveCapacity / this->capacity_);

}