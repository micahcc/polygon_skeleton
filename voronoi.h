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
private:

    class Implementation;

public:
    typedef std::shared_ptr<Voronoi> Ptr;
    typedef std::shared_ptr<const Voronoi> ConstPtr;

    struct Node;

    class Edge
    {
    public:
        typedef std::shared_ptr<Edge> Ptr;

        // original points that this edge separates
        std::set<size_t> parents;

        // endpoints for the edge
        std::shared_ptr<Node> nodes[2];

        // other edges adjacent to this one
        std::set<std::shared_ptr<Edge>> neighbors;
    };

    class Node
    {
    public:
        typedef std::shared_ptr<Node> Ptr;

        // position
        float x, y;

        // original points that this node separates (2 or 3)
        std::set<size_t> parents;

        // edges attached to this node
        std::set<std::shared_ptr<Edge>> edges;

        // other nodes attached to this one by an edge
        std::set<std::shared_ptr<Node>> neighbors;

    private:
        friend Voronoi::Implementation;
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

    std::vector<Edge::Ptr> m_edges;
    std::vector<Node::Ptr> m_nodes;

};

//Voronoi computeVoronoi(const std::vector<Point>& points);
Voronoi::Ptr computeVoronoi(const std::vector<Point>& points);

