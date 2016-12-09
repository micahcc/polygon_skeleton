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
};

Voronoi computeVoronoi(const std::vector<Point>& points);

// Typedefs
struct Intersection;
struct BeachCompare;
struct Circle;
typedef std::set<Intersection, BeachCompare> BeachLineT;

// Helper Functions
Point getIntersection(float sweep_y, const Point& p, const Point& r, float sign);
Circle solveCircle(const Point& p, const Point& q, const Point& r);

// Helper Structures
struct Circle
{
    Point center;
    float radius;
};

struct Intersection
{
    Intersection(const Point* pt0, const Point* pt1, float sign) :
        pt0(pt0), pt1(pt1), sign(sign) {};
    Intersection() : pt0(nullptr), pt1(nullptr), sign(1) {};

    const Point* pt0;
    const Point* pt1;
    float sign;  // add or subtract the radical?
};

struct Line
{
    Point pt0;
    Point pt1;
};

struct BeachCompare
{
    float sweep_y;

    bool operator()(const Intersection& lhs, const Intersection& rhs) const
    {
        // Compares the x index of thwo intersections. In the case where a point
        // is missing (nullptr), there is no intersection (or intersection is at
        // positive or negative infinity (depending on sign)
        std::cerr << "Comparing: ("<< lhs.pt0 << ", " << lhs.pt1 <<")"
            << " to (" << rhs.pt0 << ", " << rhs.pt1 << std::endl;
        std::cerr << "Using sweep = " << sweep_y << std::endl;
        if(lhs.pt0 == nullptr && rhs.pt0 == nullptr) {
            // shouldn't really happen, maybe on the first point
            std::cerr << "Comparing two nulls!" << std::endl;
            return lhs.sign < rhs.sign;
        } else if(lhs.pt1 == nullptr) {
            // left point is the lowest parabola (goes forever)
            return true;
        } else if(rhs.pt1 == nullptr) {
            // right point is largest parabola (goes forever)
            return false;
        } else {
            // get intersection of left two parabolas, and compare x with
            // intersection of right two
            Point right = getIntersection(sweep_y, *rhs.pt0, *rhs.pt1, rhs.sign);
            Point left  = getIntersection(sweep_y, *lhs.pt0, *lhs.pt1, lhs.sign);
            return left.x < right.x;
        }
    }
};

struct CircleEvent
{
    Intersection int0;
    Intersection int1;
    Circle circle;

    bool operator<(const CircleEvent& rhs) const
    {
        return circle.center.y - circle.radius <
            rhs.circle.center.y - rhs.circle.radius;
    }
};


class CircleQueue
{
public:
    bool empty() const
    {
        return m_queue.empty();
    }

    const CircleEvent& front() const
    {
        auto it = m_queue.begin();
        return *it;
    };

    const CircleEvent& back() const
    {
        auto it = m_queue.end();
        it--;
        return *it;
    };

    void pop_back()
    {
        auto it = m_queue.end();
        it--;
        m_queue.erase(it);
    };

    void pop_front()
    {
        auto it = m_queue.begin();
        m_queue.erase(it);
    };

    void insert(const Intersection& int0, const Intersection& int1)
    {
        assert(int0.pt1 == int1.pt0);
        auto pt0 = int0.pt0;
        auto pt1 = int0.pt1;
        auto pt2 = int1.pt1;
        assert(pt0 != pt1 && pt1 != pt2);

        CircleEvent evt;
        evt.circle = solveCircle(*pt0, *pt1, *pt2);
        evt.int0 = int0;
        evt.int1 = int1;
        m_queue.insert(evt);
    }

    void erase(const Intersection& int0, const Intersection& int1)
    {
        assert(int0.pt1 == int1.pt0);
        Circle circle = solveCircle(*int0.pt0, *int0.pt1, *int1.pt1);
        float end = circle.center.y - circle.radius;
        CircleEvent tmp{int0, int1, circle};
        m_queue.erase(tmp);
    }


private:

    std::set<CircleEvent> m_queue;
};

