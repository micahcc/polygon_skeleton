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
    };

    auto lines = computeVoronoi(points);
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
    for(const auto& line: lines) {
        std::cerr << "(" << line.pt0.x << ", " << line.pt0.y << ")"
            << " -- (" << line.pt1.x << ", " << line.pt1.y << ")" << std::endl;
        doc << svg::Line(
                svg::Point(line.pt0.x, line.pt0.y),
                svg::Point(line.pt1.x, line.pt1.y),
                svg::Stroke(.1, svg::Color::Black));
    }

    doc.save();
}
