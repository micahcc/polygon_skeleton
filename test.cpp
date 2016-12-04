
int main()
{
    std::vector<Point> points{
        {1, 2},
        {4, 8},
        {2, 7},
        {2.2, 3},
    };

    Voronoi result = voronoi(points);
    result.getLines();
}
