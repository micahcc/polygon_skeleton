
test: test.o voronoi.o
	clang++ $^ -o $@ -std=c++11

%.o: %.cpp
	clang++ $< -c -o $@ -std=c++11

clean:
	rm -f test.o voronoi.o test
