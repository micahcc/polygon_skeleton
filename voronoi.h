#pragma once

#include <set>
#include <cassert>
#include <algorithm>
#include <tuple>
#include <vector>
#include <memory>
#include <iostream>

#include "geometry.h"

using std::sqrt;
using std::tuple;
using std::get;

struct Voronoi
{
public:
    typedef std::shared_ptr<Voronoi> Ptr;
    typedef std::shared_ptr<const Voronoi> ConstPtr;

    struct Node;

    struct Edge
    {
        typedef std::shared_ptr<Edge> Ptr;
        static constexpr size_t n_parents = 2;

        uint8_t n_neighbors;

        size_t parents[n_parents]; // original points that this edge separates
        std::shared_ptr<Node> nodes[2];      // endpoints for the edge
        std::shared_ptr<Edge> neighbors[6];  // other edges adjacent to this one
    };

    struct Node
    {
        typedef std::shared_ptr<Node> Ptr;

        uint8_t n_parents;   // # original points that this node separates (2 or 3)
        uint8_t n_edges;     // # edges attached to this node (2 or 3)
        uint8_t n_neighbors; // # other nodes connected to this by an edge 2 or 3

        size_t parents[3];  // original points that this node separates (2 or 3)
        std::shared_ptr<Edge> edges[3];     // edges attached to this node
        std::shared_ptr<Node> neighbors[3]; // other nodes attached to this one
                                            // by an edge
        float x, y;         // position
    };

    Voronoi(const std::vector<Point>& points);

    const std::vector<Edge::Ptr> getEdges() const
    {
        return m_edges;
    }

    const std::vector<Node::Ptr> getNodes() const
    {
        return m_nodes;
    }

private:

    class Implementation;

    std::vector<Edge::Ptr> m_edges;
    std::vector<Node::Ptr> m_nodes;

};

//Voronoi computeVoronoi(const std::vector<Point>& points);
Voronoi::Ptr computeVoronoi(const std::vector<Point>& points);

