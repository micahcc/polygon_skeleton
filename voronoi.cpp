#include "voronoi.h"
#include "geometry.h"

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



void processEvent(const CircleEvent& event, CircleQueue& events,
        BeachLineT& beach, std::vector<Line>& lines)
{
    // This essentially locks in the results of a single point (the middle part
    // of the two intersections, that means we must remove all events related to
    // this point

    // find intersections to the left and right on the beach line, so we can
    // create a new event for when they meet
    auto it = beach.find(event.int0);
    it--;
    auto left_int = *it;

    it = beach.find(event.int1);
    it++;
    auto right_int = *it;

    // there are still possible events involving the
    events.erase(left_int, event.int0);
    events.erase(event.int1, right_int);
    events.insert(left_int, right_int);

    // delete arc (i.e. erase both intersections related to it
    beach.erase(event.int0);
    beach.erase(event.int1);

    Line line0{*event.int0.pt0, *event.int0.pt1};
    Line line1{*event.int1.pt0, *event.int1.pt1};
    lines.push_back(line0);
    lines.push_back(line1);
}

void processPoint(BeachLineT& beach, BeachCompare& beach_compare,
        CircleQueue& events,
        const Point& pt, Voronoi& voronoi)
{
    // Update sweep location in beach line

    beach_compare.sweep_y = pt.y;

    // insert two new intersections
    Intersection dummy{&pt, nullptr, -1};
    bool success;
    BeachLineT::iterator it1, it2;
    it1 = beach.lower_bound(dummy);
    const Point* ptA = nullptr;
    const Point* ptB = nullptr;
    const Point* ptC = nullptr;
    const Point* ptD = nullptr;
    if(beach.empty()) {
        // add null intersection
        // no intersections to erase
        beach.emplace(&pt, nullptr, -1);
        beach.emplace(&pt, nullptr, +1);
    } else if(it1 == beach.end()) {
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

        std::tie(it2, success) = beach.emplace(ptB, ptC, -1);
        it1 = it2; it1--;
        events.insert(*it1, *it2);
        beach.emplace(ptB, ptC, +1);
    } else if(it1 == beach.begin()) {
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
        std::tie(it1, success) = beach.emplace(ptA, ptB, +1);
        it2 = it1; it2++;
        events.insert(*it1, *it2);
        beach.emplace(ptA, ptB, -1);
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
        std::tie(it2, success) = beach.emplace(ptB, ptD, -1);
        it1 = it2; it1--;
        events.insert(*it1, *it2);
        auto left_event = *it1;

        // Insert new intersection int beach, then create a new event for the
        // old upper intersection and the new one
        std::tie(it1, success) = beach.emplace(ptB, ptD, +1);
        it2 = it1; it2++;
        events.insert(*it1, *it2);
        auto right_event = *it2;

        // Erase the event that involved the meeting of our previous left and
        // right intersections (since we got in the middle)
        events.erase(left_event, right_event);

        // TODO add dangling too
        voronoi.m_points.emplace_back();
        voronoi.m_points.back().x = 0.5*(ptB->x + pt.x);
        voronoi.m_points.back().y = 0.5*(ptB->y + pt.y);
    }
}

Voronoi computeVoronoi(const std::vector<Point>& points)
{
    // Sort by decreasing y
    std::vector<size_t> ordered(points.size());
    for(size_t ii = 0; ii < points.size(); ii++) ordered[ii] = ii;
    std::sort(ordered.begin(), ordered.end(),
            [&](size_t ii, size_t jj) { return points[ii].y > points[jj].y; });

    std::vector<Line> vlines;
    Voronoi voronoi;
    BeachCompare beach_compare;
    BeachLineT beach(beach_compare);
    CircleQueue events;

    // Travel downward so at each step take
    size_t ii = 0;
    while(!events.empty() && ii < ordered.size()) {

        if(events.empty()) {
            processPoint(beach, beach_compare, events, points[ordered[ii]],
                    voronoi);
            ii++;
        } else {
            const auto& evt = events.front();
            if(points[ordered[ii]].y > evt.circle.center.y - evt.circle.radius) {
                processPoint(beach, beach_compare, events, points[ordered[ii]],
                        voronoi);
                ii++;
            } else {
                events.pop_front();
                processEvent(evt, events, beach, vlines);
            }
        }
    }

    return voronoi;
}

