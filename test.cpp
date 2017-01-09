#include <vector>
#include <tuple>
#include "geometry.h"
#include "simple_svg.hpp"
#include "voronoi.h"

int main()
{
    const size_t POINT_RADIUS = 5;
    std::vector<Point> points{
        {600, 600},
        {900, 500},
        {500, 100},
        {100, 300},
        {1000, 100},
    };

    Voronoi graph(points);
    //auto vpoints = result.getPoints();

    //result.getLines();
    svg::Dimensions dimensions(1000, 1000);
    svg::Document doc("my_svg.svg", svg::Layout(dimensions, svg::Layout::BottomLeft));

    for(const auto& pt: points) {
        std::cerr << "(" << pt.x << ", " << pt.y << ")" << std::endl;
        doc << svg::Circle(svg::Point(pt.x, pt.y), 5, svg::Fill(svg::Color::Black));
    }

    // Condensed notation, parenthesis isolate temporaries that are inserted into parents.
    //doc << svg::LineChart(svg::Dimensions(1000, 1000));
    for(const auto& edge: graph.getEdges()) {
        const auto& pt0 = *edge->nodes[0];
        const auto& pt1 = *edge->nodes[1];
        std::cerr << "(" << pt0.x << ", " << pt0.y << ")"
            << " -- (" << pt1.x << ", " << pt1.y << ")" << std::endl;
        doc << svg::Line(
                svg::Point(pt0.x, pt0.y),
                svg::Point(pt1.x, pt1.y),
                svg::Stroke(2, svg::Color::Black));
    }

    doc.save();
}
