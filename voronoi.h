#pragma once

#include <set>
#include <cassert>
#include <algorithm>
#include <tuple>
#include <vector>
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


    Voronoi(const std::vector<Point>& points);

    std::vector<tuple<Point, Point>> getLines()
    {
        std::vector<tuple<Point, Point>> out;
        for(const auto& segment : m_edges) {
            tuple<Point, Point> tmp{
                m_points[std::get<0>(segment)],
                m_points[std::get<1>(segment)],
            };
            out.push_back(tmp);
        }

        return out;
    }

    const std::vector<Edge> getEdges() const
    {
        return m_edges;
    }

    const std::vector<Node> getNodes() const
    {
        return m_nodes;
    }

    struct Edge
    {
        typedef std::shared_ptr<Edge> Ptr;

        constexpr size_t n_parents = 2;
        uint8_t n_neighbors;

        size_t parents[2];   // original points that this edge separates
        Node::Ptr nodes[2];      // endpoints for the edge
        Edge::Ptr neighbors[6];  // other edges adjacent to this one
    };

    struct Node
    {
        typedef std::shared_ptr<Node> Ptr;

        uint8_t n_parents;   // # original points that this node separates (2 or 3)
        uint8_t n_edges;     // # edges attached to this node (2 or 3)
        uint8_t n_neighbors; // # other nodes connected to this by an edge 2 or 3

        size_t parents[3];  // original points that this node separates (2 or 3)
        Edge::Ptr edges[3];     // edges attached to this node
        Node::Ptr neighbors[3]; // other nodes attached to this one by an edge
        float x, y;         // position
    };

private:

    class Implementation;

    std::vector<Edge> m_edges;
    std::vector<Nodes> m_nodes;

};

//Voronoi computeVoronoi(const std::vector<Point>& points);
Voronoi::Ptr computeVoronoi(const std::vector<Point>& points);

