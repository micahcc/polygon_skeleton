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


// State Holder (beach line and output voronoi diagram)
class VoronoiComputer
{
    public:
    VoronoiComputer() : m_beach(m_beach_compare)
    {
    }

    std::vector<Line> compute(const std::vector<Point>& points);

    // debug
    std::vector<Line> lines;

    private:
    void processPoint(const Point& pt);
    void processEvent(const CircleEvent& event);

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
    float y_s = sweep_y;
    float term1 = (p.y*r.x - p.x*r.y + (p.x - r.x)*y_s) / (p.y - r.y);
    float rad = sqrt(p.x*p.x + p.y*p.y - 2*p.x*r.x + r.x*r.x - 2*p.y*r.y + r.y*r.y)*
            sqrt(-p.y + y_s)*sqrt(-r.y + y_s)/(p.y - r.y);

    // choose +- radical to be between
    Point q;
    q.x = term1 + sign*rad;

    // use one of the parabolas to find y
    q.y = 0.5*(p.x*p.x + p.y*p.y - 2*p.x*q.x + q.x*q.x - y_s*y_s)/(p.y - y_s);

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
    std::cerr << "Processing point: " << pt << std::endl;

    // Update sweep location in beach line so that insertion takes place at the
    // right location
    m_beach_compare.sweep_y = pt.y;

    // insert two new intersections
    Intersection dummy{&pt, nullptr, -1};
    bool success;
    BeachLineT::iterator it1, it2;
    it1 = m_beach.lower_bound(dummy);
    const Point* ptA = nullptr;
    const Point* ptB = nullptr;
    const Point* ptC = nullptr;
    const Point* ptD = nullptr;
    if(m_beach.empty()) {
        // add null intersection
        // no intersections to erase
        m_beach.emplace(&pt, nullptr, -1);
        m_beach.emplace(&pt, nullptr, +1);
    } else if(it1 == m_beach.end()) {
        //  Since two intersections coming together forms an event and on the
        //  end there is no second intersection to merge, nothing to erase
        //
        //  iterator:      it--      it
        //  struct:        pt0 pt1   NOTHING
        //  points:         A   B
        //  intersection:       B    C (new)

        // shouldn't happen
        throw -1;
        it1--;
        ptA = it1->pt0;
        ptB = it1->pt1;
        ptC = &pt;

        std::tie(it2, success) = m_beach.emplace(ptB, ptC, -1);
        it1 = it2; it1--;
        m_events.insert(*it1, *it2);
        m_beach.emplace(ptB, ptC, +1);
    } else if(it1 == m_beach.begin()) {
        //  Since two intersections coming together forms an event and on the
        //  end there is no second intersection to merge, nothing to erase
        //
        //  iterator:                it
        //  struct:        NOTHING   pt0 pt1
        //  points:                   B   C
        //  intersection:  (new) A    B
        // intersection >= so take the first point
        throw -1;

        // shouldn't happen
        ptA = &pt;
        ptB = it1->pt0;
        ptC = it1->pt1;
        std::tie(it1, success) = m_beach.emplace(ptA, ptB, +1);
        it2 = it1; it2++;
        m_events.insert(*it1, *it2);
        m_beach.emplace(ptA, ptB, -1);
    } else {
        // In between two previous intersections, on the parabolar for the
        // shared point
        //
        //  iterator:      it--      it
        //  struct:        pt0 pt1   pt0 pt1
        //  points:         A   B     B   C
        //  new inter:          B  D  B
        // intersection >= so take the first point
        it2 = it1; it2++;
        ptA = it1->pt0;
        ptB = it1->pt1;
        ptC = it2->pt1;
        ptD = &pt;

        // Insert new intersection into beach, then create an event for the old
        // left and the new intersection point
        std::tie(it2, success) = m_beach.emplace(ptB, ptD, -1);
        it1 = it2; it1--;
        m_events.insert(*it1, *it2);
        auto left_event = *it1;

        // Insert new intersection int beach, then create a new event for the
        // old upper intersection and the new one
        std::tie(it1, success) = m_beach.emplace(ptB, ptD, +1);
        it2 = it1; it2++;
        m_events.insert(*it1, *it2);
        auto right_event = *it2;

        // Erase the event that involved the meeting of our previous left and
        // right intersections (since we got in the middle)
        m_events.erase(left_event, right_event);

        //// TODO add dangling too
        //voronoi.m_points.emplace_back();
        //voronoi.m_points.back().x = 0.5*(ptB->x + pt.x);
        //voronoi.m_points.back().y = 0.5*(ptB->y + pt.y);
    }
}

std::vector<Line> VoronoiComputer::compute(const std::vector<Point>& points)
{
    // Sort by decreasing y
    std::vector<size_t> ordered(points.size());
    for(size_t ii = 0; ii < points.size(); ii++) ordered[ii] = ii;
    std::sort(ordered.begin(), ordered.end(),
            [&](size_t ii, size_t jj) { return points[ii].y > points[jj].y; });

    // Travel downward so at each step take
    size_t ii = 0;
    while(!m_events.empty() && ii < ordered.size()) {

        if(m_events.empty()) {
            processPoint(points[ordered[ii]]);
            ii++;
        } else {
            const auto& evt = m_events.front();
            if(points[ordered[ii]].y > evt.circle.center.y - evt.circle.radius) {
                processPoint(points[ordered[ii]]);
                ii++;
            } else {
                m_events.pop_front();
                processEvent(evt);
            }
        }
    }

    return lines;
    //return voronoi;
}

