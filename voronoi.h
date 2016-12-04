#pragma once

#include <tuple>
#include <vector>

#include "geometry.h"

struct Voronoi
{
    std::vector<Point> m_points;
    std::vector<std::tuple<size_t, size_t>> m_segments;

    std::vector<std::tuple<Point, Point>> getLines()
    {
        std::vector<std::tuple<Point, Point>> out;
        for(const auto& segment : m_segments) {
            std::tuple<Point, Point> tmp{
                m_points[std::get<0>(segment)],
                m_points[std::get<1>(segment)],
            };
            out.push_back(tmp);
        }

        return out;
    }
};


