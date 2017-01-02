#include "voronoi.h"

#include <cmath>

#include "geometry.h"

// Types
struct Intersection;
struct BeachCompare;
struct Circle;
typedef std::set<Intersection, BeachCompare> BeachLineT;

// Helper Functions
Circle solveCircle(const Point& p, const Point& q, const Point& r);
Point getIntersection(float sweep_y, const Point& p, const Point& r, float sign);
Point getIntersection(float sweep_y, const Point& p, double x);
double sqr(double v);

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


struct Boundary
{
    // a boundary connects the midline between two points with a circle point of
    // the original two points and a third point.
    const Point* start_pt0;
    const Point* start_pt1;

    const Point* circle_pt;

    bool operator<(const Boundary& rhs)
    {
        return std::min(start_pt0->y, start_pt1->y) <
            std::min(rhs.start_pt0->y, rhs.start_pt1->y);
    }

    bool operator==(const Boundary& rhs)
    {
        return ((start_pt0 == rhs.start_pt0 &&
                 start_pt1 == rhs.start_pt1) ||
                (start_pt0 == rhs.start_pt1 &&
                 start_pt1 == rhs.start_pt0));
    }
};

struct BeachCompare
{
    BeachCompare(float* sweep_y) : sweep_y(sweep_y) {} ;
    float* sweep_y;

    bool operator()(const Intersection& lhs, const Intersection& rhs) const
    {
        // Compares the x index of thwo intersections. In the case where a point
        // is missing (nullptr), there is no intersection (or intersection is at
        // positive or negative infinity (depending on sign)
        std::cerr << "<<<Comparing: ("
            << lhs.pt0 << ", " << lhs.pt1 << ", " << lhs.sign << ") to ("
            << rhs.pt0 << ", " << rhs.pt1 << ", " << rhs.sign << ")" << std::endl;
        if(lhs.pt0)
            std::cerr << "<<<Left Point 0: " << *lhs.pt0 << std::endl;
        if(lhs.pt1)
            std::cerr << "<<<Left Point 1: " << *lhs.pt1 << std::endl;
        if(rhs.pt0)
            std::cerr << "<<<Right Point 0: " << *rhs.pt0 << std::endl;
        if(rhs.pt1)
            std::cerr << "<<<Right Point 1: " << *rhs.pt1 << std::endl;
        std::cerr << "<<<Using sweep = " << *sweep_y << std::endl;
        bool left_infinite = (lhs.pt0 == nullptr || lhs.pt1 == nullptr);
        bool right_infinite = (rhs.pt0 == nullptr || rhs.pt1 == nullptr);
        if(left_infinite && right_infinite)  {
            // shouldn't really happen, maybe on the first point
            std::cerr << "<<<Comparing two nulls!" << std::endl;
            std::cerr << "<<<" << (lhs.sign < rhs.sign) << std::endl;
            return lhs.sign < rhs.sign;
        } else if(left_infinite) {
            // left point is either the lowest parabola (goes forever) or
            // highest parabola depending on the sign. so if this is the
            // left end (sign = -1) everything is greater (return false) and if
            // it is the right end (sign = +1) everything is less (return true)
            if(lhs.sign < 0) {
                // -infinity < rhs => true
                std::cerr << "<<<" << true << std::endl;
                return true;
            } else {
                // +infinity < rhs => false
                std::cerr << "<<<" << false << std::endl;
                return false;
            }
        } else if(right_infinite) {
            // right point is iether lowest or highest parabola (depending on
            // sign). End parabolas can never be overtaken, so if this is the
            // left end (sign = -1) everything is greater (return false) and if
            // it is the right end (sign = +1) everything is less (return true)
            if(rhs.sign < 0) {
                // lhs < -infinity => false
                std::cerr << "<<<" << false << std::endl;
                return false;
            } else {
                // lhs < +infinity => true
                std::cerr << "<<<" << true << std::endl;
                return true;
            }
        } else if((lhs.pt0 == rhs.pt0 && lhs.pt1 == rhs.pt1) ||
                (lhs.pt1 == rhs.pt0 && lhs.pt0 == rhs.pt1)) {
            // Special case, two intersections of the same two parabolas
            return lhs.sign < rhs.sign;
        } else if(lhs.pt0 == lhs.pt1) {
            // Special case, intersection of two identical points is assumed to
            // be just the x value of the double-point intersection
            assert(rhs.pt0 != rhs.pt1);
            Point right = getIntersection(*sweep_y, *rhs.pt0, *rhs.pt1, rhs.sign);
            return lhs.pt0->x < right.x;
        } else if(rhs.pt0 == rhs.pt1) {
            // Special case, intersection of two identical points is assumed to
            // be just the x value of the double-point intersection
            assert(lhs.pt0 != lhs.pt1);
            Point left  = getIntersection(*sweep_y, *lhs.pt0, *lhs.pt1, lhs.sign);
            return left.x < rhs.pt0->x;
        } else {
            // get intersection of left two parabolas, and compare x with
            // intersection of right two
            std::cerr << "<<<Computing intersections" << std::endl;
            Point right = getIntersection(*sweep_y, *rhs.pt0, *rhs.pt1, rhs.sign);
            Point left  = getIntersection(*sweep_y, *lhs.pt0, *lhs.pt1, lhs.sign);
            std::cerr << "<<<" << (left.x < right.x) << std::endl;
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
        std::cerr << "<<<Inserting Event: ("
            << int0.pt0 << ", " << int0.pt1 << ", " << int0.sign << ") and ("
            << int1.pt0 << ", " << int1.pt1 << ", " << int1.sign << ")"
            << std::endl;
        assert(int0.pt0);
        assert(int0.pt1);
        assert(int1.pt0);
        assert(int1.pt1);
        std::cerr << "<<<Left Point 0: " << *int0.pt0 << std::endl;
        std::cerr << "<<<Left Point 1: " << *int0.pt1 << std::endl;
        std::cerr << "<<<Right Point 0: " << *int1.pt0 << std::endl;
        std::cerr << "<<<Right Point 1: " << *int1.pt1 << std::endl;

        //assert(int0.pt1 == int1.pt0);
        const Point *ptA, *ptB, *ptC;
        if(int0.pt0 != int1.pt0 && int0.pt1 != int1.pt0) {
            ptA = int0.pt0;
            ptB = int0.pt1;
            ptC = int1.pt0;
        } else if(int0.pt0 != int1.pt1 && int0.pt1 != int1.pt1) {
            ptA = int0.pt0;
            ptB = int0.pt1;
            ptC = int1.pt1;
        } else {
            throw -8;
        }

        CircleEvent evt;
        evt.circle = solveCircle(*ptA, *ptB, *ptC);
        evt.int0 = int0;
        evt.int1 = int1;
        m_queue.insert(evt);
    }

    void erase(const Intersection& left_int, const Intersection& right_int)
    {
        if(left_int.pt0 == nullptr)
            return;
        if(right_int.pt1 == nullptr)
            return;

        std::cerr << "<<<Erasing Event: ("
            << left_int.pt0 << ", " << left_int.pt1 << ", " << left_int.sign << ") to ("
            << right_int.pt0 << ", " << right_int.pt1 << ", " << right_int.sign << ")"
            << std::endl;

        assert(left_int.pt0);
        assert(left_int.pt1);
        assert(right_int.pt0);
        assert(right_int.pt1);

        std::cerr << "<<<Left Point 0: " << *left_int.pt0 << std::endl;
        std::cerr << "<<<Left Point 1: " << *left_int.pt1 << std::endl;
        std::cerr << "<<<Right Point 0: " << *right_int.pt0 << std::endl;
        std::cerr << "<<<Right Point 1: " << *right_int.pt1 << std::endl;

        assert(left_int.pt1 == right_int.pt0);
        Circle circle = solveCircle(*left_int.pt0, *left_int.pt1, *right_int.pt1);
        float end = circle.center.y - circle.radius;
        CircleEvent tmp{left_int, right_int, circle};
        m_queue.erase(tmp);
    }


private:

    std::set<CircleEvent> m_queue;
};


// State Holder (beach line and output voronoi diagram)
class VoronoiComputer
{
    public:
    VoronoiComputer() : m_beach_compare(&sweep_y), m_beach(m_beach_compare),
        m_min_x(std::numeric_limits<double>::infinity()),
        m_max_x(-std::numeric_limits<double>::infinity()),
        m_min_y(std::numeric_limits<double>::infinity()),
        m_max_y(-std::numeric_limits<double>::infinity())
    {
    }

    void compute(const std::vector<Point>& points);
    std::vector<Line> getLines();

    private:
    void processPoint(const Point& pt);
    void processEvent(const CircleEvent& event);

    float sweep_y;
    BeachCompare m_beach_compare;
    BeachLineT m_beach;
    CircleQueue m_events;

    double m_min_x, m_max_x, m_min_y, m_max_y;

    std::vector<Boundary> m_bounds;
};

/**
 * Functions
 */

//Voronoi computeVoronoi(const std::vector<Point>& points)
std::vector<Line> computeVoronoi(const std::vector<Point>& points)
{
    VoronoiComputer computer;
    computer.compute(points);
    return computer.getLines();
}

Point getIntersection(float sweep_y, const Point& p, const Point& r, float sign)
{
    // Given any point q, its distance from sweep is:
    //
    // (q.y - y_s)
    //
    // Its distance from an event p is:
    //
    // sqrt((q.x - p.x)^2 + (q.y - p.y)^2)
    //
    // the parabola equation is then:
    //
    // (q.y - y_s)^2 = (q.x - p.x)^2 + (q.y - p.y)^2
    //
    // For two event points, r and p, we can solve q.x -- the intersection
    // location must satisfy:
    //
    // sqrt( sqr(q_x - p_x) + sqr(q_y - p_y) ) ==
    // sqrt( sqr(q_x - r_x) + sqr(q_y - r_y) )
    //
    // and
    //
    // sqrt( sqr(q_x - r_x) + sqr(q_y - r_y) ) == q_y - sweep_y
    std::cerr << "<<<<Intersection of:\n<<<<" << p << "\n<<<<" << r << "\n<<<<" <<
        sweep_y << std::endl;

    // Solve for x first
    float y_s = sweep_y;
    Point q;
    if(std::abs(p.y - sweep_y) < 0.00001) {
        // parabola around p has no width, just select point on parabola r at
        // p.x
        std::cerr << "<<<<p_y == sweep_y" << std::endl;
        q.x = p.x;
        q.y = 0.5*( sqr(q.x) - 2*q.x*r.x + sqr(r.x) + sqr(r.y) - sqr(y_s))/(r.y - y_s);
    } else if(std::abs(r.y - sweep_y) < 0.00001) {
        // parabola around r has no width, just select point on parabola q at
        // r.x
        std::cerr << "<<<<r_y == sweep_y" << std::endl;
        q.x = r.x;
        q.y = 0.5*(p.x*p.x + p.y*p.y - 2*p.x*q.x + q.x*q.x - y_s*y_s)/(p.y - y_s);
    } else {
        float term1 = (p.y*r.x - p.x*r.y + (p.x - r.x)*y_s) / (p.y - r.y);
        float rad =
            sqrt(p.x*p.x + p.y*p.y - 2*p.x*r.x + r.x*r.x - 2*p.y*r.y + r.y*r.y)*
            sqrt(p.y - y_s)*sqrt(r.y - y_s)/(p.y - r.y);

        std::cerr << "<<<<"
            << sqrt(p.x*p.x + p.y*p.y - 2*p.x*r.x + r.x*r.x - 2*p.y*r.y + r.y*r.y)
            << ", " << sqrt( p.y - y_s) << ", " << sqrt(r.y - y_s)/(p.y - r.y) << std::endl;;
        std::cerr << "<<<<" << term1 << " + " << sign << " * " << rad << std::endl;
        // choose +- radical to be between

        assert(!std::isnan(term1));
        assert(!std::isnan(rad));
        q.x = term1 + sign*rad;

        // use one of the parabolas to find y
        q.y = 0.5*(p.x*p.x + p.y*p.y - 2*p.x*q.x + q.x*q.x - y_s*y_s)/(p.y - y_s);
    }

    std::cerr << "Solution: " << q << std::endl;
    std::cerr << "<<<<Sweep line distance0: " << (q.y - sweep_y) << "\n"
        << "<<<<Solution distance0: "
        << std::sqrt( sqr(p.x - q.x ) + sqr( p.y - q.y )) << std::endl;
    std::cerr << "<<<<Sweep line distance1: " << (q.y - sweep_y) << "\n"
        << "<<<<Solution distance1: "
        << std::sqrt( sqr(r.x - q.x ) + sqr( r.y - q.y )) << std::endl;
    return q;
}

double sqr(double v) {
    return v*v;
}

Circle solveCircle(const Point& p, const Point& q, const Point& r)
{
    // find minimal point
    Circle circle;
    circle.center.x =
        0.5*(p.y*sqr(q.x) + p.y*sqr(q.y) - (p.y - q.y)*sqr(r.x) -
                (p.y - q.y)*sqr(r.y) - (sqr(p.x) + sqr(p.y))*q.y +
                (sqr(p.x) + sqr(p.y) - sqr(q.x) - sqr(q.y))*r.y) /
        (p.y*q.x - p.x*q.y - (p.y - q.y)*r.x + (p.x - q.x)*r.y);

    circle.center.y =
        -0.5*(p.x*sqr(q.x) + p.x*sqr(q.y) - (p.x - q.x)*sqr(r.x) -
                (p.x - q.x)*sqr(r.y) - (sqr(p.x) + sqr(p.y))*q.x +
                (sqr(p.x) + sqr(p.y) - sqr(q.x) - sqr(q.y))*r.x) /
        (p.y*q.x - p.x*q.y - (p.y - q.y)*r.x + (p.x - q.x)*r.y);

    circle.radius = sqrt(sqr(p.x - circle.center.x) + sqr(p.y - circle.center.y));
    return circle;
};

// VoronoiComputer Implementation
void VoronoiComputer::processEvent(const CircleEvent& event)
{
    std::cerr << "Processing Event for: [" << *event.int0.pt0 << " -- "
        << *event.int0.pt1 << "], " << *event.int1.pt0 << " -- "
        << *event.int1.pt1 << "]\n";

    // This essentially locks in the results of a single point (the middle part
    // of the two intersections, that means we must remove all events related to
    // this point

    // find intersections to the left and right on the beach line, so we can
    // create a new event for when they meet
    auto it = m_beach.find(event.int0);
    std::cerr << "Left Int: [" << *(*it).pt0 << " -- " << *(*it).pt1 << std::endl;
    it--;
    auto new_left_int = *it;
    it++;
    it++;
    std::cerr << "Right Int: [" << *(*it).pt0 << " -- " << *(*it).pt1 << std::endl;
    it++;
    auto new_right_int = *it;

    // We need to remove events involing either of the original
    // intersections, since this intersection is destroyed when it meets another
    // intersection (at a circle point). The events affected are the meeting of
    // the two intersections and their left and right neighbors. If the
    // intersections are dummy intersections (null point at one side, then don't
    // erase)
    m_events.erase(new_left_int, event.int0);
    m_events.erase(event.int1, new_right_int);

    std::cerr << "New Left/Right: " << new_left_int.pt0 << ", "
        << new_left_int.pt1 << " -- " << new_right_int.pt0 << ", "
        << new_right_int.pt1 << std::endl;
    if(new_left_int.pt0 != nullptr && new_right_int.pt1 != nullptr) {
        std::cerr << "Inserting event for: " << *new_left_int.pt0
            << ", " << *new_left_int.pt1 << ", -- " << *new_right_int.pt0
            << ", " << *new_right_int.pt1 << std::endl;
        m_events.insert(new_left_int, new_right_int);
    }

    // delete arc (i.e. erase both intersections related to it
    std::cerr << "Erasing from beach" << std::endl;
    m_beach.erase(event.int0);
    m_beach.erase(event.int1);

//    Line line0{*event.int0.pt0, *event.int0.pt1};
//    Line line1{*event.int1.pt0, *event.int1.pt1};
//    lines.push_back(line0);
//    lines.push_back(line1);
//
//    // finish off boundaries related to the middle point
//    auto itb = m_bounds.find(new_left_int.pt0, new_left_int.pt1);
//    assert(itb != m_bounds.end());
//    itb->circle_pt = new_right_int.pt1;
//
//    itb = m_bounds.find(new_right_int.pt0, new_right_int.pt1);
//    assert(itb != m_bounds.end());
//    itb->circle_pt = new_left_int.pt0;

    const Point *ptA, *ptB, *ptC;
    if(event.int0.pt0 != event.int1.pt0 && event.int0.pt1 != event.int1.pt0) {
        ptA = event.int0.pt0;
        ptB = event.int0.pt1;
        ptC = event.int1.pt0;
    } else if(event.int0.pt0 != event.int1.pt1 && event.int0.pt1 != event.int1.pt1) {
        ptA = event.int0.pt0;
        ptB = event.int0.pt1;
        ptC = event.int1.pt1;
    } else {
        throw -8;
    }

    // The new center point connects to bisectors of each of the individual
    // pairs of points, these are rays from the center of the event circle to
    // each of the bisectors
    m_bounds.emplace_back(Boundary{ptA, ptB, ptC});
    m_bounds.emplace_back(Boundary{ptB, ptA, ptC});
    m_bounds.emplace_back(Boundary{ptC, ptA, ptC});
}

void VoronoiComputer::processPoint(const Point& pt)
{
    std::cerr << "<----------------------" << std::endl;
    std::cerr << "<Processing point: " << pt << std::endl;

    // Update sweep location in beach line so that insertion takes place at the
    // right location

    *m_beach_compare.sweep_y = pt.y;

    // insert two new intersections in between existing intersections
    Intersection dummy{&pt, &pt, -1};
    bool success;
    BeachLineT::iterator it1, it2, it_new;
    const Point* ptA = nullptr;
    const Point* ptB = nullptr;
    const Point* ptC = nullptr;
    const Point* ptD = nullptr;
    if(m_beach.empty()) {
        std::cerr << "<<<Beach empty, inserting special" << std::endl;
        // add null intersection
        // no intersections to erase
        m_beach.emplace(nullptr, &pt, -1);
        m_beach.emplace(&pt, nullptr, +1);
    } else {
        // In between two previous intersections, on the parabolar for the
        // shared point
        //
        //  iterator:      it1       it2
        //  struct:        pt0 pt1   pt0 pt1
        //  points:         A   B     B   C
        //  new inter:          B  D  B
        // intersection >= so take the first point
        std::cerr << "<<Finding beach location" << std::endl;
        it1 = m_beach.lower_bound(dummy);
        std::cerr << "<<Lower bound: (" << it1->pt0 << " -- " << it1->pt1 <<
            ", " << it1->sign << ")" << std::endl;
        if(it1->pt0) {
            std::cerr << "<<Pt0: " << *it1->pt0 << std::endl;
        }
        if(it1->pt1) {
            std::cerr << "<<Pt1: " << *it1->pt1 << std::endl;
        }
        std::cerr << "<<Done" << std::endl;
        it2 = it1; it1--;
        ptB = it1->pt1;
        ptD = &pt;

        std::cerr << "B: " << ptB << std::endl
            << "D: " << ptD << std::endl;
        std::cerr << "<<<Inserting intersection of " << *ptB
            << " -- " << *ptD << "into beach" << std::endl;
        // Insert new intersection into beach, then create an event for the old
        // left and the new intersection point
        std::tie(it_new, success) = m_beach.emplace(ptB, ptD, -1);
        assert(success);
        if(it1->pt0 != nullptr)
            m_events.insert(*it1, *it_new);

        // Insert new intersection int beach, then create a new event for the
        // old upper intersection and the new one
        std::tie(it_new, success) = m_beach.emplace(ptD, ptB, +1);
        assert(success);
        if(it2->pt1 != nullptr)
            m_events.insert(*it_new, *it2);

        //// Erase the event that involved the meeting of our previous left and
        //// right intersections (since we got in the middle)
        //if(it1->pt0 != nullptr && it2->pt1 != nullptr)
        //    m_events.erase(*it1, *it2);
    }


    std::cerr << "<......................" << std::endl;
}

void VoronoiComputer::compute(const std::vector<Point>& points)
{
    for(const auto& pt : points) {
        m_min_x = std::min<double>(pt.x, m_min_x);
        m_max_x = std::max<double>(pt.x, m_max_x);
        m_min_y = std::min<double>(pt.y, m_min_y);
        m_max_y = std::max<double>(pt.y, m_max_y);
    }

    std::cerr << "Sorting points" << std::endl;
    // Sort by decreasing y
    std::vector<size_t> ordered(points.size());
    for(size_t ii = 0; ii < points.size(); ii++) ordered[ii] = ii;
    std::sort(ordered.begin(), ordered.end(),
            [&](size_t ii, size_t jj) { return points[ii].y > points[jj].y; });

    std::cerr << "Ordered points: " << std::endl;
    for(size_t ii : ordered) {
        std::cerr << points[ii] << std::endl;
    }
    std::cerr << std::endl;

    // Travel downward so at each step take
    size_t ii = 0;
    while(!m_events.empty() || ii < ordered.size()) {
        std::cerr << "Next Point" << points[ordered[ii]] << std::endl;

        if(m_events.empty()) {
            std::cerr << "Events Empty, processing next point" << std::endl;
            processPoint(points[ordered[ii]]);
            ii++;
        } else if(ii == ordered.size()) {
            std::cerr << "Points Done, processing next event" << std::endl;
            auto evt = m_events.front();
            m_events.pop_front();
            processEvent(evt);
        } else {
            auto evt = m_events.front();
            std::cerr << "Next point: " << points[ordered[ii]].y
                << ", Next Event: " << evt.circle.center.y - evt.circle.radius
                << std::endl;
            if(points[ordered[ii]].y > evt.circle.center.y - evt.circle.radius) {
                processPoint(points[ordered[ii]]);
                ii++;
            } else {
                m_events.pop_front();
                processEvent(evt);
            }
        }

        std::cerr << "Final Beach: " << std::endl;
        for(const auto& inter: m_beach) {
            std::cerr << "(" << inter.pt0 << ", " << inter.pt1 << ", "
                << inter.sign << ")" << std::endl;
            if(inter.pt0) std::cerr << "Point 0: " << *inter.pt0 << std::endl;
            if(inter.pt1) std::cerr << "Point 1: " << *inter.pt1 << std::endl;
        }
    }

    //return voronoi;
}

std::vector<Line> VoronoiComputer::getLines()
{
    std::vector<Line> out;
    for(auto it = m_bounds.begin(); it != m_bounds.end(); ++it) {
        // start point is a bisector
        Point pt0{0.5 * (it->start_pt0->x + it->start_pt1->x),
                  0.5 * (it->start_pt0->y + it->start_pt1->y)};

        // was closed
        auto c = solveCircle(*it->start_pt0, *it->start_pt1, *it->circle_pt);
        Point pt1 = c.center;
        out.emplace_back(Line{pt0, pt1});
    }
    return out;
}
