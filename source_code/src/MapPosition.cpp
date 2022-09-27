#include "optimizer/MapPosition.hpp"


/// <summary>
/// Intitilizes object.
/// </summary>
MapPosition::MapPosition()
	: CloserVertexID(0), FurtherVertexID(0), EdgeSegmentID(0)
{ }


/// <summary>
/// Intitilizes object.
/// </summary>
/// <param name="firstVertex">ID of the first vertex of the edge.</param>
/// <param name="secondVertex">ID of the second vertex of the edge.</param>
/// <param name="segment">ID of the segment on the edge.</param>
/// <param name="graph">Graph object where the position is (to get edge ID).</param>
MapPosition::MapPosition(vertex_t firstVertex, vertex_t secondVertex, 
	int32_t segment, Graph& graph)
	: EdgeSegmentID(std::move(segment))
{
	this->EdgeSegmentID = segment;
	this->EdgeID = boost::edge(firstVertex, secondVertex, graph).first;

	// Make sure `firstVertex` is smaller than `secondVertex`.
	if (secondVertex < firstVertex)
		std::swap(firstVertex, secondVertex);

	if (segment < graph[this->EdgeID].GetNumSegments() / 2)
	// Segment is closer to the vertex with smaller index.
	{
		this->CloserVertexID = firstVertex;
		this->FurtherVertexID = secondVertex;
	}
	else
	// Segment is closer to the vertex with higher index
	{
		this->CloserVertexID = secondVertex;
		this->FurtherVertexID = firstVertex;
	}
}

/// <summary>
/// Changes parameters of the object and computes other (if specified).
/// </summary>
/// <param name="closer">New closer vertex ID.</param>
/// <param name="further">New further vertex ID.</param>
/// <param name="segment">New segment ID (if `-1` then let it 
/// precompute by the function (the last or the first segment)).</param>
/// <param name="graph">Graph object to get all necessary info.</param>
void MapPosition::ChangeParameters(vertex_t closer, vertex_t further,
	int32_t segment, Graph& graph)
{
	this->CloserVertexID = closer;
	this->FurtherVertexID = further;
	this->EdgeID = boost::edge(closer, further, graph).first;

	this->EdgeSegmentID = segment;

	if (segment == -1)
	// Segment to be computed by the function (last segment in the moving direction).
	{
		// If increasing in order, start ID is 0.
		this->EdgeSegmentID = 0;

		// Decreasing in the segment IDs (from higher to lower).
		if (closer < further)
			this->EdgeSegmentID = graph[this->EdgeID].GetNumSegments() - 1;
	}
}


/// <summary>
/// Computes distance from the closer vertex.
/// Distance is only approximate (there could be difference of 1 segment length).
/// We assume this function will be used only as a rough estimate of the 
/// distance from closer vertex.
/// </summary>
/// <param name="graph">Graph object to get the distance from.</param>
/// <param name="segmentLength">Length of one segment of the edge.</param>
/// <returns>Returns estimated distance of the position from the closer vertex.</returns>
double MapPosition::GetDistanceFromCloserVertex(Graph& graph, double segmentLength)
{
	int totalNumSegments = graph[this->EdgeID].GetNumSegments();
	int numSegmentsFromCloser = this->EdgeSegmentID;
	
	if (this->CloserVertexID > this->FurtherVertexID)
	// Closer vertex has higher ID than further -> appropriately compute distance
	// (edge segment ID is in rising order from the vertex with smaller ID).
	{
		numSegmentsFromCloser = totalNumSegments - this->EdgeSegmentID;
	}

	return (double)numSegmentsFromCloser * segmentLength;
}