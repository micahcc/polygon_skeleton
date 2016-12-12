#include "voronoi.h"
#include "geometry.h"

// Types
struct Intersection;
struct BeachCompare;
struct Circle;
typedef std::set<Intersection, BeachCompare> BeachLineT;

// Helper Functions
Circle solveCircle(const Point& p, const Point& q, const Point& r);
Point getIntersection(float sweep_y, const Point& p, const Point& r, float sign);
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
        if(int0.pt0)
            std::cerr << "<<<Left Point 0: " << *int0.pt0 << std::endl;
        if(int0.pt1)
            std::cerr << "<<<Left Point 1: " << *int0.pt1 << std::endl;
        if(int1.pt0)
            std::cerr << "<<<Right Point 0: " << *int1.pt0 << std::endl;
        if(int1.pt1)
            std::cerr << "<<<Right Point 1: " << *int1.pt1 << std::endl;

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
        std::cerr << "<<<Erasing Event: ("
            << int0.pt0 << ", " << int0.pt1 << ", " << int0.sign << ") to ("
            << int1.pt0 << ", " << int1.pt1 << ", " << int1.sign << ")"
            << std::endl;
        if(int0.pt0)
            std::cerr << "<<<Left Point 0: " << *int0.pt0 << std::endl;
        if(int0.pt1)
            std::cerr << "<<<Left Point 1: " << *int0.pt1 << std::endl;
        if(int1.pt0)
            std::cerr << "<<<Right Point 0: " << *int1.pt0 << std::endl;
        if(int1.pt1)
            std::cerr << "<<<Right Point 1: " << *int1.pt1 << std::endl;

        assert(int0.pt1 == int1.pt0);
        Circle circle = solveCircle(*int0.pt0, *int0.pt1, *int1.pt1);
        float end = circle.center.y - circle.radius;
        CircleEvent tmp{int0, int1, circle};
        m_queue.erase(tmp);
    }


private:

    std::set<CircleEvent> m_queue;
};


// State Holder (beach line and output voronoi diagram)
class VoronoiComputer
{
    public:
    VoronoiComputer() : m_beach_compare(&sweep_y), m_beach(m_beach_compare)
    {
    }

    std::vector<Line> compute(const std::vector<Point>& points);

    // debug
    std::vector<Line> lines;

    private:
    void processPoint(const Point& pt);
    void processEvent(const CircleEvent& event);

    float sweep_y;
    BeachCompare m_beach_compare;
    BeachLineT m_beach;
    CircleQueue m_events;
};

/**
 * Functions
 */

//Voronoi computeVoronoi(const std::vector<Point>& points)
std::vector<Line> computeVoronoi(const std::vector<Point>& points)
{
    VoronoiComputer computer;
    return computer.compute(points);
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
    // For two event points, r and p, we can solve q.x -- the interction
    // location
    std::cerr << "<<<<Intersection of:\n<<<<" << p << "\n<<<<" << r << "\n<<<<" <<
        sweep_y << std::endl;

    Point q;
    if(std::abs(p.y - sweep_y) < 0.00001) {
        // parabola around p has no width, just select point on parabola r at
        // p.x
        q.x = p.x;
        q.y = r.y - std::sqrt( sqr(r.y - sweep_y) - sqr(r.x - p.x) );
    } else if(std::abs(r.y - sweep_y) < 0.00001) {
        // parabola around r has no width, just select point on parabola q at
        // r.x
        q.x = r.x;
        q.y = p.y - std::sqrt( sqr(p.y - sweep_y) - sqr(p.x - r.x) );
    } else {
        float y_s = sweep_y;
        float term1 = (p.y*r.x - p.x*r.y + (p.x - r.x)*y_s) / (p.y - r.y);
        float rad =
            sqrt(p.x*p.x + p.y*p.y - 2*p.x*r.x + r.x*r.x - 2*p.y*r.y + r.y*r.y)*
            sqrt(-p.y + y_s)*sqrt(-r.y + y_s)/(p.y - r.y);

        std::cerr << "<<<<"
            << sqrt(p.x*p.x + p.y*p.y - 2*p.x*r.x + r.x*r.x - 2*p.y*r.y + r.y*r.y)
            << ", " << sqrt(-p.y + y_s) << ", " << sqrt(-r.y + y_s)/(p.y - r.y) << std::endl;;
        std::cerr << "<<<<" << term1 << " + " << sign << " * " << rad << std::endl;
        // choose +- radical to be between
        q.x = term1 + sign*rad;

        // use one of the parabolas to find y
        q.y = 0.5*(p.x*p.x + p.y*p.y - 2*p.x*q.x + q.x*q.x - y_s*y_s)/(p.y - y_s);
    }

    std::cerr << "Solution: " << q << std::endl;
    std::cerr << "<<<<Sweep line distance0: " << (p.y - sweep_y) << "\n"
        << "<<<<Solution distance0: "
        << std::sqrt( sqr(p.x - q.x ) + sqr( p.y - q.y )) << std::endl;
    std::cerr << "<<<<Sweep line distance1: " << (r.y - sweep_y) << "\n"
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
        1/2*(p.y*sqr(q.x) + p.y*sqr(q.y) - (p.y - q.y)*sqr(r.x) -
                (p.y - q.y)*sqr(r.y) - (sqr(p.x) + sqr(p.y))*q.y +
                (sqr(p.x) + sqr(p.y) - sqr(q.x) - sqr(q.y))*r.y) /
        (p.y*q.x - p.x*q.y - (p.y - q.y)*r.x + (p.x - q.x)*r.y);

    circle.center.y =
        -1/2*(p.x*sqr(q.x) + p.x*sqr(q.y) - (p.x - q.x)*sqr(r.x) -
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
    it--;
    auto left_int = *it;

    it = m_beach.find(event.int1);
    it++;
    auto right_int = *it;

    // there are still possible events involving the
    m_events.erase(left_int, event.int0);
    m_events.erase(event.int1, right_int);
    m_events.insert(left_int, right_int);

    // delete arc (i.e. erase both intersections related to it
    m_beach.erase(event.int0);
    m_beach.erase(event.int1);

    Line line0{*event.int0.pt0, *event.int0.pt1};
    Line line1{*event.int1.pt0, *event.int1.pt1};
    lines.push_back(line0);
    lines.push_back(line1);
}

void VoronoiComputer::processPoint(const Point& pt)
{
    std::cerr << "<----------------------" << std::endl;
    std::cerr << "<Processing point: " << pt << std::endl;

    // Update sweep location in beach line so that insertion takes place at the
    // right location

    *m_beach_compare.sweep_y = pt.y;

    // insert two new intersections
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
        //  iterator:      it        it2
        //  struct:        pt0 pt1   pt0 pt1
        //  points:         A   B     B   C
        //  new inter:          B  D  B
        // intersection >= so take the first point
        std::cerr << "<<Finding beach location" << std::endl;
        it1 = m_beach.lower_bound(dummy);
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
        std::tie(it_new, success) = m_beach.emplace(ptB, ptD, +1);
        assert(success);
        if(it2->pt1 != nullptr)
            m_events.insert(*it_new, *it2);

        // Erase the event that involved the meeting of our previous left and
        // right intersections (since we got in the middle)
        if(it1->pt0 != nullptr && it2->pt1 != nullptr)
            m_events.erase(*it1, *it2);

        //// TODO add dangling too
        //voronoi.m_points.emplace_back();
        //voronoi.m_points.back().x = 0.5*(ptB->x + pt.x);
        //voronoi.m_points.back().y = 0.5*(ptB->y + pt.y);
    }

    std::cerr << "<......................" << std::endl;
}

std::vector<Line> VoronoiComputer::compute(const std::vector<Point>& points)
{
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
        } else {
            const auto& evt = m_events.front();
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

    return lines;
    //return voronoi;
}

