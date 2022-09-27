#ifndef BOOST_HEADER_H_
#define BOOST_HEADER_H_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/graph/astar_search.hpp>

#include <boost/graph/graphviz.hpp>

#include "Edge.hpp"
#include "Node.hpp"


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Node, Edge > Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

typedef std::vector<std::vector<double>> Double2DMatrix;

typedef std::pair<double, std::vector<vertex_t>> LengthPathPair;


#endif // !BOOST_HEADER_H_


#pragma once
