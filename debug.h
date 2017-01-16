#pragma once

#include <random>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "simple_svg.hpp"

template <typename IntersectionContainer, typename EventContainer>
void draw_state(const IntersectionContainer& intersections,
        const EventContainer& events, double sweep_y)
{
    svg::Dimensions dimensions(1200, 1200);

    static int count = 0;
    std::ostringstream oss;
    oss << "state_" << std::setfill('0') << std::setw(5) << count++ << ".svg";
    svg::Document doc(oss.str(), svg::Layout(dimensions, svg::Layout::BottomLeft));

    if(intersections.empty())
        return;

    auto it = intersections.end();
    it--;

    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    for(const auto& curr_int : intersections) {
        if(!curr_int.pt_right)
            break;

        min_x = std::min<double>(min_x, curr_int.pt_right->x) - 100;
        max_x = std::max<double>(max_x, curr_int.pt_right->x) + 100;
    }

    for(const auto& curr_int : intersections) {
        if(!curr_int.pt_right)
            break;

       draw_parabola(doc, min_x, max_x, *curr_int.pt_right, sweep_y);
    }

    for(auto it = events.cbegin(); it != events.cend(); ++it) {
        doc << svg::Circle(svg::Point(it->circle.center.x, it->circle.center.y),
                5, svg::Fill(svg::Color::Red));

        doc << svg::Circle(svg::Point(it->circle.center.x, it->circle.center.y),
                2*it->circle.radius, svg::Fill(svg::Color::Transparent),
                svg::Stroke(1, svg::Color::Red));
    }

    doc.save();
}


template <typename IntersectionContainer, typename EventContainer>
void draw_state(const IntersectionContainer& intersections,
        const EventContainer& events, double start_sweep_y, double stop_sweep_y)
{
    if(std::isnan(start_sweep_y) || std::isnan(stop_sweep_y))
        return;

    for(double y = start_sweep_y; y > stop_sweep_y; y--)
        draw_state(intersections, events, y);
}

inline
void draw_parabola(svg::Document& doc, double min_x, double max_x,
    const Point& pt, double sweep_y)
{
    size_t h1 = std::hash<double>{}(pt.x);
    size_t h2 = std::hash<double>{}(pt.y);
    size_t r = ((h1 >> 1 ^ h2) & 0x0000FF);
    size_t g = ((h1 >> 1 ^ h2) & 0x00FF00) >> 8;
    size_t b = ((h1 >> 1 ^ h2) & 0xFF0000) >> 16;

    std::ostringstream tmp;
    svg::Color color(r, g, b);

    doc << svg::Circle(svg::Point(pt.x, pt.y), 5, svg::Fill(color));

    tmp << &pt;
    doc << svg::Text(svg::Point(pt.x, pt.y), tmp.str(), svg::Fill(color));

    // line defined by
    // (x - p.x)^2 + (y - p.y)^2 == (y - sweep_y)
    svg::Polyline parabola(svg::Fill(svg::Color::Transparent), svg::Stroke(1, color));
    if(sweep_y == pt.y) {
        for(double y = sweep_y; y < sweep_y + 1000; y += 1) {
            parabola << svg::Point(pt.x, y);
        }
    } else {
        for(double x = min_x; x < max_x; x += 1) {
            double y = 0.5*(pt.x*pt.x + pt.y*pt.y - 2*pt.x*x + x*x - sweep_y*sweep_y)/
                (pt.y - sweep_y);
            parabola << svg::Point(x, y);
        }
    }

    doc << parabola;
    doc.save();
}


