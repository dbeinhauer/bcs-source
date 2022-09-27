#ifndef MAPPOSITION_H_
#define MAPPOSITION_H_


#include "boost_header.hpp"


/// <summary>
/// Class to store info about the position in the map (edge and its segment).
/// </summary>
class MapPosition
{
public:
	// ID of the edge.
	edge_t EdgeID;
	// Always vertex to get (the vehicle is moving towards the `closerVertexID`).
	vertex_t CloserVertexID;
	// Either the start vertex of the edge or the part which 
	// will not be visited (in case of start inside the segment).
	vertex_t FurtherVertexID;
	// ID of the segment.
	int32_t EdgeSegmentID;

	MapPosition();
	MapPosition(vertex_t firstVertex, vertex_t secondVertex, int32_t segment, Graph& graph);

	void ChangeParameters(vertex_t closer, vertex_t further, int32_t segment, Graph& graph);

	double GetDistanceFromCloserVertex(Graph& graph, double segmentLength);
};

#endif