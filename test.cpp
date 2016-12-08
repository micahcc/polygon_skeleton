#include <vector>
#include "geometry.h"
#include "simple_svg.hpp"
#include "voronoi.h"

int main()
{
    const size_t POINT_RADIUS = 5;
    std::vector<Point> points{
        {10, 20},
        {40, 80},
        {22, 70},
        {22, 33},
    };

    Voronoi result = computeVoronoi(points);
    auto vpoints = result.getPoints();

    //result.getLines();
    svg::Dimensions dimensions(100, 100);
    svg::Document doc("my_svg.svg", svg::Layout(dimensions, svg::Layout::BottomLeft));

    for(const auto& pt : points) {
        doc << svg::Circle(svg::Point(pt.x, pt.y), POINT_RADIUS,
                svg::Fill(svg::Color::Blue), svg::Stroke(svg::Color::Black));
    }

    for(const auto& pt : vpoints) {
        doc << svg::Circle(svg::Point(pt.x, pt.y), POINT_RADIUS,
                svg::Fill(svg::Color::Red), svg::Stroke(svg::Color::Black));
    }

//    // Long notation.  Local variable is created, children are added to varaible.
//    LineChart chart(5.0);
//    Polyline polyline_a(Stroke(.5, Color::Blue));
//    Polyline polyline_b(Stroke(.5, Color::Aqua));
//    Polyline polyline_c(Stroke(.5, Color::Fuchsia));
//    polyline_a << Point(0, 0) << Point(10, 30)
//        << Point(20, 40) << Point(30, 45) << Point(40, 44);
//    polyline_b << Point(0, 10) << Point(10, 22)
//        << Point(20, 30) << Point(30, 32) << Point(40, 30);
//    polyline_c << Point(0, 12) << Point(10, 15)
//        << Point(20, 14) << Point(30, 10) << Point(40, 2);
//    chart << polyline_a << polyline_b << polyline_c;
//    doc << chart;
//
//    // Condensed notation, parenthesis isolate temporaries that are inserted into parents.
//    doc << (LineChart(Dimensions(65, 5))
//        << (Polyline(Stroke(.5, Color::Blue)) << Point(0, 0) << Point(10, 8) << Point(20, 13))
//        << (Polyline(Stroke(.5, Color::Orange)) << Point(0, 10) << Point(10, 16) << Point(20, 20))
//        << (Polyline(Stroke(.5, Color::Cyan)) << Point(0, 5) << Point(10, 13) << Point(20, 16)));
//
//    doc << Circle(Point(80, 80), 20, Fill(Color(100, 200, 120)), Stroke(1, Color(200, 250, 150)));
//
//    doc << Text(Point(5, 77), "Simple SVG", Color::Silver, Font(10, "Verdana"));
//
//    doc << (Polygon(Color(200, 160, 220), Stroke(.5, Color(150, 160, 200))) << Point(20, 70)
//        << Point(25, 72) << Point(33, 70) << Point(35, 60) << Point(25, 55) << Point(18, 63));
//
//    doc << Rectangle(Point(70, 55), 20, 15, Color::Yellow);
//
    doc.save();
}
