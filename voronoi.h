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

class VoronoiComputer;

struct Voronoi
{
public:

    std::vector<tuple<Point, Point>> getLines()
    {
        std::vector<tuple<Point, Point>> out;
        for(const auto& segment : m_segments) {
            tuple<Point, Point> tmp{
                m_points[std::get<0>(segment)],
                m_points[std::get<1>(segment)],
            };
            out.push_back(tmp);
        }

        return out;
    }

    const std::vector<Point> getPoints() const { return m_points; };


    std::vector<Point> m_points;
    std::vector<tuple<size_t, size_t>> m_segments;

private:
    friend VoronoiComputer;
};

//Voronoi computeVoronoi(const std::vector<Point>& points);
std::vector<Line> computeVoronoi(const std::vector<Point>& points);

