#ifndef NODE_H_
#define NODE_H_

#include <iostream>

// Struct representing the coordinates.
struct location
{
	double Latitude, Longitude;
};

/// <summary>
/// Class representing road junction.
/// </summary>
class Node
{
public:
	// Current node ID.
	int32_t NewID_;
	// Old node ID (in the original map). 
	// To potentially project the node to the real map).
	int32_t OldID_;
	// ID of the nearest city.
	int32_t CityID_;
	// Distance between the node and its nearest city center.
	double CityDistance_;
	// Geographic coordinates of the node.
	location Location_;

	Node();
	Node(int32_t newID, int32_t oldID, int32_t cityID, 
		double cityDistance, double latitude, double longitude);
};

#endif