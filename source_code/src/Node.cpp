#include "optimizer/Node.hpp"


/// <summary>
/// Initializes the object. 
/// </summary>
Node::Node() 
{
	this->NewID_ = -1;
	this->OldID_ = -1;
	this->CityID_ = -1;
	this->CityDistance_ = 0;
	this->Location_ = location();
}


/// <summary>
/// Initializes the object.
/// </summary>
/// <param name="newID">Current node ID.</param>
/// <param name="oldID">Old node ID (in the original map). 
/// To potentially project the node to the real map).</param>
/// <param name="cityID">ID of the nearest city.</param>
/// <param name="cityDistance">Distance between the node and its nearest city center.</param>
/// <param name="latitude">Node latitude.</param>
/// <param name="longitude">Node longitude.</param>
Node::Node(int32_t newID, int32_t oldID, int32_t cityID, double cityDistance,
	double latitude, double longitude)
{
	this->NewID_ = newID;
	this->OldID_ = oldID;
	this->CityID_ = cityID;
	this->CityDistance_ = cityDistance;
	location loc;
	loc.Latitude = latitude;
	loc.Longitude = longitude;
	this->Location_ = loc;
}