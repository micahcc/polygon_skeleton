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
Point getIntersection(float sweep_y, const Point& p, const Point& r);
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
    Intersection(const Point* pt_left, const Point* pt_right) :
        pt_left(pt_left), pt_right(pt_right) {} ;
    Intersection() : pt_left(nullptr), pt_right(nullptr) {};

    const Point* pt_left;
    const Point* pt_right;
    float sign;  // add or subtract the radical?
};


struct Boundary
{
    // a boundary connects the midline between two points with a circle point of
    // the original two points and a third point.
    const Point* start_pt_left;
    const Point* start_pt_right;

    const Point* circle_pt;

    bool operator<(const Boundary& rhs)
    {
        return std::min(start_pt_left->y, start_pt_right->y) <
            std::min(rhs.start_pt_left->y, rhs.start_pt_right->y);
    }

    bool operator==(const Boundary& rhs)
    {
        return ((start_pt_left == rhs.start_pt_left &&
                 start_pt_right == rhs.start_pt_right) ||
                (start_pt_left == rhs.start_pt_right &&
                 start_pt_right == rhs.start_pt_left));
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
            << lhs.pt_left << ", " << lhs.pt_right << ", " << lhs.sign << ") to ("
            << rhs.pt_left << ", " << rhs.pt_right << ", " << rhs.sign << ")" << std::endl;
        if(lhs.pt_left)
            std::cerr << "<<<Left Point 0: " << *lhs.pt_left << std::endl;
        if(lhs.pt_right)
            std::cerr << "<<<Left Point 1: " << *lhs.pt_right << std::endl;
        if(rhs.pt_left)
            std::cerr << "<<<Right Point 0: " << *rhs.pt_left << std::endl;
        if(rhs.pt_right)
            std::cerr << "<<<Right Point 1: " << *rhs.pt_right << std::endl;
        std::cerr << "<<<Using sweep = " << *sweep_y << std::endl;
        bool left_infinite = (lhs.pt_left == nullptr || lhs.pt_right == nullptr);
        bool right_infinite = (rhs.pt_left == nullptr || rhs.pt_right == nullptr);
        assert(!left_infinite || !right_infinite);
        if(left_infinite) {
            // -infinity < rhs => true
            std::cerr << "<<<" << true << std::endl;
            return true;
        } else if(right_infinite) {
            // lhs < -infinity => false
            std::cerr << "<<<" << false << std::endl;
            return false;
        } else if((lhs.pt_left == rhs.pt_left && lhs.pt_right == rhs.pt_right) ||
                (lhs.pt_right == rhs.pt_left && lhs.pt_left == rhs.pt_right)) {
            // Special case, two intersections of the same two parabolas
            std::cerr << "parabolas are identical!" << std::endl;
            throw -1;
        } else if(lhs.pt_left == lhs.pt_right) {
            // Special case, intersection of two identical points is assumed to
            // be just the x value of the double-point intersection
            assert(rhs.pt_left != rhs.pt_right);
            Point right = getIntersection(*sweep_y, *rhs.pt_left, *rhs.pt_right);
            return lhs.pt_left->x < right.x;
        } else if(rhs.pt_left == rhs.pt_right) {
            // Special case, intersection of two identical points is assumed to
            // be just the x value of the double-point intersection
            assert(lhs.pt_left != lhs.pt_right);
            Point left  = getIntersection(*sweep_y, *lhs.pt_left, *lhs.pt_right);
            return left.x < rhs.pt_left->x;
        } else {
            // get intersection of left two parabolas, and compare x with
            // intersection of right two
            std::cerr << "<<<Computing intersections" << std::endl;
            Point right = getIntersection(*sweep_y, *rhs.pt_left, *rhs.pt_right);
            Point left  = getIntersection(*sweep_y, *lhs.pt_left, *lhs.pt_right);
            std::cerr << "<<<" << (left.x < right.x) << std::endl;
            return left.x < right.x;
        }
    }
};


struct CircleEvent
{
    Intersection left_int;
    Intersection right_int;
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

    void insert(const Intersection& left_int, const Intersection& right_int)
    {
        std::cerr << "<<<Inserting Event: ("
            << left_int.pt_left << ", " << left_int.pt_right << ", " << left_int.sign << ") and ("
            << right_int.pt_left << ", " << right_int.pt_right << ", " << right_int.sign << ")"
            << std::endl;
        assert(left_int.pt_left);
        assert(left_int.pt_right);
        assert(right_int.pt_left);
        assert(right_int.pt_right);
        std::cerr << "<<<Left Point 0: " << *left_int.pt_left << std::endl;
        std::cerr << "<<<Left Point 1: " << *left_int.pt_right << std::endl;
        std::cerr << "<<<Right Point 0: " << *right_int.pt_left << std::endl;
        std::cerr << "<<<Right Point 1: " << *right_int.pt_right << std::endl;

        //assert(left_int.pt_right == right_int.pt_left);
        const Point *ptA, *ptB, *ptC;
        if(left_int.pt_left != right_int.pt_left && left_int.pt_right != right_int.pt_left) {
            ptA = left_int.pt_left;
            ptB = left_int.pt_right;
            ptC = right_int.pt_left;
        } else if(left_int.pt_left != right_int.pt_right && left_int.pt_right != right_int.pt_right) {
            ptA = left_int.pt_left;
            ptB = left_int.pt_right;
            ptC = right_int.pt_right;
        } else {
            throw -8;
        }

        CircleEvent evt;
        evt.circle = solveCircle(*ptA, *ptB, *ptC);
        evt.left_int = left_int;
        evt.right_int = right_int;
        m_queue.insert(evt);
    }

    void erase(const Intersection& left_int, const Intersection& right_int)
    {
        if(left_int.pt_left == nullptr) return;
        if(right_int.pt_right == nullptr) return;

        std::cerr << "<<<Erasing Event: ("
            << left_int.pt_left << ", " << left_int.pt_right << ", " << left_int.sign << ") to ("
            << right_int.pt_left << ", " << right_int.pt_right << ", " << right_int.sign << ")"
            << std::endl;

        assert(left_int.pt_left);
        assert(left_int.pt_right);
        assert(right_int.pt_left);
        assert(right_int.pt_right);

        std::cerr << "<<<Left Point 0: " << *left_int.pt_left << std::endl;
        std::cerr << "<<<Left Point 1: " << *left_int.pt_right << std::endl;
        std::cerr << "<<<Right Point 0: " << *right_int.pt_left << std::endl;
        std::cerr << "<<<Right Point 1: " << *right_int.pt_right << std::endl;

        assert(left_int.pt_right == right_int.pt_left);
        Circle circle = solveCircle(*left_int.pt_left, *left_int.pt_right, *right_int.pt_right);
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

/**
 *  Find the intersection of two parabolas created from points and a sweep line
 *  (directrix)
 *
 *  Note: order is important
 *
 *  When a parabola is initially inserted, it is essentially a line:
 *  \  A  /
 *   \   /
 *   |\ /
 *   | v
 *   |
 * --B---------
 *
 *  \  A  /
 *   \   /
 *    \ /|
 *     v |
 *       |
 * ------B------
 *
 * The beach line then has intersections: A:B, B:A. where the sign of
 * the radical is -1 for A:B and +1 for B:A. In other words, the sign of the
 * radical determines the side of B that we are on
 *
 * @param sweep_y Position of sweep line
 * @param left parab
 * @param right_parab
 */
Point getIntersection(float sweep_y, const Point& left_parab,
        const Point& right_parab)
{
    if(left_parab.y <= right_parab.y) {
        // A = right_parab
        // B = left_parab
        // Transitioning form B to A (B:A)
        return getIntersection(sweep_y, left_parab, right_parab, +1);
    } else {
        // A = left_parab
        // B = right_parab
        // Transitioning from A to B (A:B)
        return getIntersection(sweep_y, left_parab, right_parab, -1);
    }
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
    std::cerr << "Processing Event for: [" << *event.left_int.pt_left << " -- "
        << *event.left_int.pt_right << "], " << *event.right_int.pt_left << " -- "
        << *event.right_int.pt_right << "]\n";

    // This essentially locks in the results of a single point (the middle part
    // of the two intersections, that means we must remove all events related to
    // this point


    // Suppose that intersections A and B meet (when the middle peak is beat out
    // by the side peaks):
    // l |   |     | r
    //   \  A/ |   |
    //    \ / v\B  /
    //     v    \ /
    //           v

    // |     |
    // |     |      |
    // \    C/      |
    //  \   / \     /
    //   \ /   \   /
    //    v     \ /
    //           v
    //
    // First we will remove the intersections A and B from the beach line
    //
    // Now if there are other intersections to the left of A, there would be an
    // event associated with A and that (l), same for an event to the right of
    // B (r). Those events can no longer occur because the intersections A (of
    // the left and center peak) and B (of the center and right peaks) will no
    // longer exist when the center peak has been hidden. Therefore we must
    // remove those events.
    //
    // We then have a new intersection on the beach line (C) which can have
    // events with l and r (when they exist)

    // find intersections to the left and right on the beach line, so we can
    // create a new event for when they meet

    bool success;
    BeachLineT::iterator it_new;
    auto it = m_beach.find(event.left_int);
    std::cerr << "Left Int: [" << *(*it).pt_left << " -- " << *(*it).pt_right << std::endl;
    it--;
    auto new_left_int = *it;
    it++;
    it++;
    std::cerr << "Right Int: [" << *(*it).pt_left << " -- " << *(*it).pt_right << std::endl;
    it++;
    auto new_right_int = *it;

    // We need to remove events involing either of the original
    // intersections, since this intersection is destroyed when it meets another
    // intersection (at a circle point). The events affected are the meeting of
    // the two intersections and their left and right neighbors. If the
    // intersections are dummy intersections (null point at one side, then don't
    // erase)
    m_events.erase(new_left_int, event.left_int);
    m_events.erase(event.right_int, new_right_int);

    // delete arc (i.e. erase both intersections related to it
    std::cerr << "Erasing from beach" << std::endl;
    m_beach.erase(event.left_int);
    m_beach.erase(event.right_int);

    // TODO add intersection of previous peaks to beach line
    throw -1; // the issue is that the radical sign is not very clear
    std::tie(it_new, success) = m_beach.emplace(event.left_int.pt_left,
            event.right_int.pt_right);

    std::cerr << "New Left/Right: " << new_left_int.pt_left << ", "
        << new_left_int.pt_right << " -- " << new_right_int.pt_left << ", "
        << new_right_int.pt_right << std::endl;
    if(new_left_int.pt_left != nullptr) {
        m_events.insert(new_left_int, *it_new);
    }
    if(new_right_int.pt_right != nullptr) {
        m_events.insert(*it_new, new_right_int);
    }

//    Line line0{*event.left_int.pt_left, *event.left_int.pt_right};
//    Line line1{*event.right_int.pt_left, *event.right_int.pt_right};
//    lines.push_back(line0);
//    lines.push_back(line1);
//
//    // finish off boundaries related to the middle point
//    auto itb = m_bounds.find(new_left_int.pt_left, new_left_int.pt_right);
//    assert(itb != m_bounds.end());
//    itb->circle_pt = new_right_int.pt_right;
//
//    itb = m_bounds.find(new_right_int.pt_left, new_right_int.pt_right);
//    assert(itb != m_bounds.end());
//    itb->circle_pt = new_left_int.pt_left;

    const Point *ptA, *ptB, *ptC;
    if(event.left_int.pt_left != event.right_int.pt_left &&
            event.left_int.pt_right != event.right_int.pt_left) {
        ptA = event.left_int.pt_left;
        ptB = event.left_int.pt_right;
        ptC = event.right_int.pt_left;
    } else if(event.left_int.pt_left != event.right_int.pt_right &&
            event.left_int.pt_right != event.right_int.pt_right) {
        ptA = event.left_int.pt_left;
        ptB = event.left_int.pt_right;
        ptC = event.right_int.pt_right;
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
    Intersection dummy{&pt, &pt};
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
        m_beach.emplace(nullptr, &pt);
        m_beach.emplace(&pt, nullptr);
    } else {
        // In between two previous intersections, on the parabolar for the
        // shared point
        //
        //  iterator:      it1       it2
        //  struct:        pt_left pt_right   pt_left pt_right
        //  points:         A   B     B   C
        //  new inter:          B  D  B
        // intersection >= so take the first point
        std::cerr << "<<Finding beach location" << std::endl;
        it1 = m_beach.lower_bound(dummy);
        std::cerr << "<<Lower bound: (" << it1->pt_left << " -- " << it1->pt_right <<
            ", " << it1->sign << ")" << std::endl;
        if(it1->pt_left) {
            std::cerr << "<<pt_left: " << *it1->pt_left << std::endl;
        }
        if(it1->pt_right) {
            std::cerr << "<<pt_right: " << *it1->pt_right << std::endl;
        }
        std::cerr << "<<Done" << std::endl;
        it2 = it1; it1--;
        ptB = it1->pt_right;
        ptD = &pt;

        std::cerr << "B: " << ptB << std::endl
            << "D: " << ptD << std::endl;
        std::cerr << "<<<Inserting intersection of " << *ptB
            << " -- " << *ptD << "into beach" << std::endl;
        // Insert new intersection into beach, then create an event for the old
        // left and the new intersection point
        std::tie(it_new, success) = m_beach.emplace(ptB, ptD);
        assert(success);
        if(it1->pt_left != nullptr)
            m_events.insert(*it1, *it_new);

        // Insert new intersection int beach, then create a new event for the
        // old upper intersection and the new one
        std::tie(it_new, success) = m_beach.emplace(ptD, ptB);
        assert(success);
        if(it2->pt_right != nullptr)
            m_events.insert(*it_new, *it2);

        //// Erase the event that involved the meeting of our previous left and
        //// right intersections (since we got in the middle)
        //if(it1->pt_left != nullptr && it2->pt_right != nullptr)
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
            std::cerr << "(" << inter.pt_left << ", " << inter.pt_right << ", "
                << inter.sign << ")" << std::endl;
            if(inter.pt_left) std::cerr << "Point 0: " << *inter.pt_left << std::endl;
            if(inter.pt_right) std::cerr << "Point 1: " << *inter.pt_right << std::endl;
        }
    }

    //return voronoi;
}

std::vector<Line> VoronoiComputer::getLines()
{
    std::vector<Line> out;
    for(auto it = m_bounds.begin(); it != m_bounds.end(); ++it) {
        // start point is a bisector
        Point pt_left{0.5 * (it->start_pt_left->x + it->start_pt_right->x),
                  0.5 * (it->start_pt_left->y + it->start_pt_right->y)};

        // was closed
        auto c = solveCircle(*it->start_pt_left, *it->start_pt_right, *it->circle_pt);
        Point pt_right = c.center;
        out.emplace_back(Line{pt_left, pt_right});
    }
    return out;
}
